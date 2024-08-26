// Shasta.
#include "mode3-AnchorGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "MarkerGraph.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>



// Create the AnchorGraph and its vertices given a vector of AnchorIds.
AnchorGraph::AnchorGraph(const vector<AnchorId>& anchorIds) :
    anchorIds(anchorIds)
{

    // Check that the AnchorIds are sorted and distinct.
    for(uint64_t i=1; i<anchorIds.size(); i++) {
        SHASTA_ASSERT(anchorIds[i-1] < anchorIds[i]);
    }

    // Create the vertices.
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        AnchorGraphVertex vertex;
        vertex.localAnchorId = localAnchorId;
        const vertex_descriptor v = add_vertex(vertex, *this);
        vertexDescriptors.push_back(v);
    }
}



void AnchorGraph::addEdgeFromLocalAnchorIds(
    uint64_t localAnchorId0,
    uint64_t localAnchorId1,
    const MarkerGraphEdgePairInfo& info,
    uint64_t coverage)
{
    boost::add_edge(
        vertexDescriptors[localAnchorId0],
        vertexDescriptors[localAnchorId1],
        AnchorGraphEdge(info, coverage), *this);
}



// Write a AnchorGraph in graphviz format.
void AnchorGraph::writeGraphviz(
    const string& name,
    const AnchorGraphDisplayOptions& options,
    const MarkerGraph& markerGraph) const
{
    ofstream out(name + ".dot");

    const AnchorGraph& graph = *this;
    out << "digraph " << name << " {\n";

    BGL_FORALL_VERTICES(v, graph, AnchorGraph) {
        out << getAnchorId(v);

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "[";
        }

        if(options.labels) {
            out << "label=\"";
            out << getAnchorId(v) << "\\n" << markerGraph.edgeCoverage(getAnchorId(v));
            out << "\" ";
        }

        if(options.tooltips) {
            out << "tooltip=\"";
            out << getAnchorId(v);
            out << "\" ";
        }

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "]";
        }
        out << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];
        if(not options.showNonTransitiveReductionEdges and edge.isNonTransitiveReductionEdge) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        out <<
            getAnchorId(v0) << "->" <<
            getAnchorId(v1);

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << " [";
        }

        if(edge.isNonTransitiveReductionEdge) {
            out << "style=dashed ";
        }

        if(options.tooltips) {
            out <<
                "tooltip=\"" <<
                getAnchorId(v0) << "->" <<
                getAnchorId(v1) << " ";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << " " <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << " " <<
                edge.info.offsetInBases << "\" ";
        }

        if(options.labels) {
            out <<
                "label=\"";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << "\\n" <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << "\\n" <<
                edge.info.offsetInBases << "\" ";

        }

        // Color.
        if(options.colorEdges) {
            const double correctedJaccard = edge.info.correctedJaccard();
            if(correctedJaccard <= options.redJ) {
                out << " color=red ";
            } else if(correctedJaccard >= options.greenJ) {
                out << " color=green ";
            } else {
                const double hue = (correctedJaccard - options.redJ) / (3. * (options.greenJ - options.redJ));
                out << " color=\"" << hue << ",1,1\" ";
            }
        }

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << "]";
        }
        out << ";\n";
    }

    out << "}\n";
}



void AnchorGraph::writeEdgeCoverageHistogram(const string& fileName) const
{
    const AnchorGraph& primaryGraph = *this;

    // Create a histogram indexed by histogram[coverage][commonCount].
    vector< vector<uint64_t> > histogram;

    // Loop over all edges.
    BGL_FORALL_EDGES(e, primaryGraph, AnchorGraph) {
        const AnchorGraphEdge& edge = primaryGraph[e];
        const uint64_t coverage = edge.coverage;
        const uint64_t commonCount = edge.info.common;
        SHASTA_ASSERT(coverage <= commonCount);

        // Increment the histogram, making space as necessary.
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1);
        }
        vector<uint64_t>& h = histogram[coverage];
        if(commonCount >= h.size()) {
            h.resize(commonCount + 1, 0);
        }
        ++h[commonCount];
    }

    // Write out the histogram.
    ofstream csv(fileName);
    csv << "Coverage,Common count,Loss,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        const vector<uint64_t>& h = histogram[coverage];
        for(uint64_t commonCount=0; commonCount<h.size(); commonCount++) {
            const uint64_t frequency = h[commonCount];

            if(frequency > 0) {
                const uint64_t loss = commonCount - coverage;
                csv << coverage << ",";
                csv << commonCount << ",";
                csv << loss << ",";
                csv << frequency << "\n";
            }
        }
    }
}



// Remove cross-edges.
// This removes an edge v0->v1 if the following are all true:
// - It is not marked as removed by transitive reduction.
// - Its coverage is at most lowCoverageThreshold.
// - Its estimated offset is at least minOffset.
// - v0 has at least one out-edge with coverage at least highCoverageThreshold
//   (ignoring edges marked as removed by transitive reduction).
// - v1 has at least one in-edge with coverage at least highCoverageThreshold.
//   (ignoring edges marked as removed by transitive reduction).
void AnchorGraph::removeCrossEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    uint64_t minOffset,
    bool debug)
{
    AnchorGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];

        // If it is marked as removed by transitive reduction, skip it.
        if(edge.isNonTransitiveReductionEdge) {
            continue;
        }

        // Check coverage.
        if(edge.coverage > lowCoverageThreshold) {
            continue;
        }

        // Check estimated offset.
        if(edge.info.offsetInBases < int64_t(minOffset)) {
            continue;
        }

        // Check out-edges of v0.
        const vertex_descriptor v0 = source(e, graph);
        bool v0HasStrongOutEdge = false;
        BGL_FORALL_OUTEDGES(v0, e0, graph, AnchorGraph) {
            // If it is marked as removed by transitive reduction, ignore it.
            if(graph[e0].isNonTransitiveReductionEdge) {
                continue;
            }
            if(graph[e0].coverage >= highCoverageThreshold) {
                v0HasStrongOutEdge = true;
                break;
            }
        }
        if(not v0HasStrongOutEdge) {
            continue;
        }

        // Check in-edges of v1.
        const vertex_descriptor v1 = target(e, graph);
        bool v1HasStrongOutEdge = false;
        BGL_FORALL_INEDGES(v1, e1, graph, AnchorGraph) {
            // If it is marked as removed by transitive reduction, ignore it.
            if(graph[e1].isNonTransitiveReductionEdge) {
                continue;
            }
            if(graph[e1].coverage >= highCoverageThreshold) {
                v1HasStrongOutEdge = true;
                break;
            }
        }
        if(not v1HasStrongOutEdge) {
            continue;
        }

        // If all above checks passed, this edge will be removed.
        edgesToBeRemoved.push_back(e);
        if(debug) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            cout << "Removing cross edge " <<
                getAnchorId(v0) << "->" <<
                getAnchorId(v1) << endl;
        }
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// Remove edges for which loss = (commonCount - coverage) / commonCount > maxLoss
void AnchorGraph::removeWeakEdges(double maxLoss, bool debug)
{
    AnchorGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];
        const double loss = double(edge.info.common - edge.coverage) / double(edge.info.common);
        if(loss > maxLoss) {
            edgesToBeRemoved.push_back(e);

            if(debug) {
                const vertex_descriptor v0 = source(e, graph);
                const vertex_descriptor v1 = target(e, graph);
                cout << "Removing weak edge " <<
                    getAnchorId(v0) << "->" <<
                    getAnchorId(v1) << ", loss " << loss << endl;
            }
        }
    }



    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



void AnchorGraph::separateStrands(
    const Anchors& anchors,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
{
    AnchorGraph& anchorGraph = *this;

    findReverseComplementAnchors(anchors, markers);
    findReverseComplementVertices();
    findReverseComplementEdges();

    // Gather pairs of reverse complemented edges, by coverage.
    vector< vector< pair<edge_descriptor, edge_descriptor> > > edgePairs;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        const edge_descriptor eRc = anchorGraph[e].eRc;
        if(e < eRc) {
            const uint64_t coverage = anchorGraph[e].coverage;
            SHASTA_ASSERT(coverage == anchorGraph[eRc].coverage);
            if(coverage >= edgePairs.size()) {
                edgePairs.resize(coverage + 1);
            }
            edgePairs[coverage].push_back({e, eRc});
        }
    }

    /*
    cout << "Number of edge pairs by coverage:" << endl;
    uint64_t edgePairCount = 0;
    for(uint64_t coverage=0; coverage<edgePairs.size(); coverage++) {
        const uint64_t n = edgePairs[coverage].size();
        if(n > 0) {
            edgePairCount += n;
            cout << coverage << "," << n << endl;
        }
    }
    SHASTA_ASSERT(2 * edgePairCount == num_edges(anchorGraph));
    */


    // The main strand separation code begins here.
    // We create a disjoint sets data structure and add edges
    // in order of decreasing coverage.
    // Edges that would cause two reverse complemented vertices to end up
    // in the same connected component are skipped.
    // In the end, we are left with two connected components, one per strand.

    SHASTA_ASSERT(0);
}



void AnchorGraph::findReverseComplementAnchors(
    const Anchors& anchors,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
{

    // Index the AnchorIds by their first MarkerInterval.
    class AnchorInfo {
    public:
        uint64_t localAnchorId;
        MarkerInterval firstMarkerInterval;
        AnchorInfo(
            uint64_t localAnchorId,
            MarkerInterval firstMarkerInterval) :
            localAnchorId(localAnchorId), firstMarkerInterval(firstMarkerInterval) {}
        bool operator<(const AnchorInfo& that) const
        {
            return firstMarkerInterval < that.firstMarkerInterval;
        }
    };
    vector<AnchorInfo> anchorInfos;
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        const AnchorId anchorId = anchorIds[localAnchorId];
        const Anchor anchor = anchors[anchorId];
        anchorInfos.push_back(AnchorInfo(localAnchorId, anchor[0]));
    }
    sort(anchorInfos.begin(), anchorInfos.end());



    // Find the local anchor id corresponding to the reverse complement
    // of a given local anchor id.
    reverseComplementAnchor.resize(anchorIds.size());
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        const AnchorId anchorId = anchorIds[localAnchorId];
        const Anchor anchor = anchors[anchorId];
        const MarkerInterval& firstMarkerInterval = anchor[0];

        // Get the reverse complemented OrientedReadId.
        const OrientedReadId orientedReadId = firstMarkerInterval.orientedReadId;
        OrientedReadId orientedReadIdRc = orientedReadId;
        orientedReadIdRc.flipStrand();


        // Get the reverse complemented MarkerInterval.
        const uint32_t markerCount = uint32_t(markers[orientedReadId.getValue()].size());
        MarkerInterval firstMarkerIntervalRc;
        firstMarkerIntervalRc.orientedReadId = orientedReadIdRc;
        firstMarkerIntervalRc.ordinals[0] = markerCount - 1 - firstMarkerInterval.ordinals[1];
        firstMarkerIntervalRc.ordinals[1] = markerCount - 1 - firstMarkerInterval.ordinals[0];

        // Look it up in our AnchorInfo vector.
        const AnchorInfo anchorInfoRc(invalid<uint64_t>, firstMarkerIntervalRc);
        const auto it = lower_bound(anchorInfos.begin(), anchorInfos.end(), anchorInfoRc);
        SHASTA_ASSERT(it != anchorInfos.end());
        SHASTA_ASSERT(it->firstMarkerInterval == firstMarkerIntervalRc);

        // Store the local anchor id of the reverse complement anchor.
        reverseComplementAnchor[localAnchorId] = it->localAnchorId;
    }

    // Sanity check.
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        const uint64_t localAnchorIdRc = reverseComplementAnchor[localAnchorId];
        SHASTA_ASSERT(localAnchorIdRc != localAnchorId);
        SHASTA_ASSERT(reverseComplementAnchor[localAnchorIdRc] == localAnchorId);
    }
}



void AnchorGraph::findReverseComplementVertices()
{
    AnchorGraph& anchorGraph = *this;

    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        const vertex_descriptor v = vertexDescriptors[localAnchorId];
        const uint64_t localAnchorIdRc = reverseComplementAnchor[localAnchorId];
        const vertex_descriptor vRc = vertexDescriptors[localAnchorIdRc];
        anchorGraph[v].vRc = vRc;
    }

    // Sanity check.
    BGL_FORALL_VERTICES(v, anchorGraph, AnchorGraph) {
        const AnchorGraphVertex& vertex = anchorGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA_ASSERT(vRc != v);
        const AnchorGraphVertex& vertexRc = anchorGraph[vRc];
        SHASTA_ASSERT(vertexRc.vRc == v);
    }

}



void AnchorGraph::findReverseComplementEdges()
{
    AnchorGraph& anchorGraph = *this;

    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        const vertex_descriptor v0 = source(e, anchorGraph);
        const vertex_descriptor v1 = target(e, anchorGraph);
        const vertex_descriptor v0Rc = anchorGraph[v0].vRc;
        const vertex_descriptor v1Rc = anchorGraph[v1].vRc;

        // The reverse complement edge is v1Rc->v0Rc.
        edge_descriptor eRc;
        bool edgeWasFound = false;
        tie(eRc, edgeWasFound) = boost::edge(v1Rc, v0Rc, anchorGraph);
        SHASTA_ASSERT(edgeWasFound);

        // Store it.
        anchorGraph[e].eRc = eRc;
    }

    // Sanity check.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        const AnchorGraphEdge& edge = anchorGraph[e];
        const edge_descriptor eRc = edge.eRc;
        SHASTA_ASSERT(eRc != e);
        const AnchorGraphEdge& edgeRc = anchorGraph[eRc];
        SHASTA_ASSERT(edgeRc.eRc == e);
    }
}



#if 0
// Given sets of two primary in-edges and two primary out-edges,
// find primary mid-edges in-between that can be used for detangling.
void GlobalPathGraph::searchForDetangling(
    const array<MarkerGraphEdgeId, 2>& in,
    const array<MarkerGraphEdgeId, 2>& out,
    uint64_t highCommonCountThreshold,
    uint64_t lowCommonCountThreshold,
    const Assembler& assembler,
    array<array<vector<MarkerGraphEdgeId>, 2>, 2>& mid)
{
    // Loop over the primary journeys of oriented reads in the "in" primary edges.
    // Only use the journey portion following the "in" primary edges.
    array<vector<MarkerGraphEdgeId>, 2> inFollowers;
    array<vector<uint64_t>, 2> inFollowersCommonCount;
    for(uint64_t i=0; i<2; i++) {
        assembler.markerGraph.followPrimaryJourneysForward(in[i], inFollowers[i], inFollowersCommonCount[i]);
    }



    // Find inFollowers that have high common count with in[0]
    // and low common count with in[1], or vice versa.
    array<vector<MarkerGraphEdgeId>, 2> inCandidates;
    {
        uint64_t i0 = 0;
        uint64_t i1 = 0;
        uint64_t end0 = inFollowers[0].size();
        uint64_t end1 = inFollowers[1].size();
        while(i0<end0 and i1<end1) {
            const MarkerGraphEdgeId edgeId0 = inFollowers[0][i0];
            const MarkerGraphEdgeId edgeId1 = inFollowers[1][i1];

            if(edgeId0 < edgeId1) {
                // edgeId0 is in inFollowers[0] but not in inFollowers[1].
                if(inFollowersCommonCount[0][i0] >= highCommonCountThreshold) {
                    inCandidates[0].push_back(edgeId0);
                }
                ++i0;
            }

            else if(edgeId1 < edgeId0) {
                // edgeId1 is in inFollowers[1] but not in inFollowers[0].
                if(inFollowersCommonCount[1][i1] >= highCommonCountThreshold) {
                    inCandidates[1].push_back(edgeId1);
                }
                ++i1;
            }

            else {
                // edgeId0 is in inFollowers[0] and in inFollowers[1].
                const uint64_t common0 = inFollowersCommonCount[0][i0];
                const uint64_t common1 = inFollowersCommonCount[1][i1];
                if(common0 >= highCommonCountThreshold and common1 <= lowCommonCountThreshold) {
                    inCandidates[0].push_back(edgeId0);
                }
                else if(common1 >= highCommonCountThreshold and common0 <= lowCommonCountThreshold) {
                    inCandidates[1].push_back(edgeId1);
                }
                ++i0;
                ++i1;
            }
        }
    }



    // Loop over the primary journeys of oriented reads in the "out" primary edges.
    // Only use the journey portion preceding the "out" primary edges.
    array<vector<MarkerGraphEdgeId>, 2> outPreceders;
    array<vector<uint64_t>, 2> outPrecedersCommonCount;
    for(uint64_t i=0; i<2; i++) {
        assembler.markerGraph.followPrimaryJourneysBackward(out[i], outPreceders[i], outPrecedersCommonCount[i]);
    }



    // Find outPreceders that have high common count with out[0]
    // and low common count with out[1], or vice versa.
    array<vector<MarkerGraphEdgeId>, 2> outCandidates;
    {
        uint64_t i0 = 0;
        uint64_t i1 = 0;
        uint64_t end0 = outPreceders[0].size();
        uint64_t end1 = outPreceders[1].size();
        while(i0<end0 and i1<end1) {
            const MarkerGraphEdgeId edgeId0 = outPreceders[0][i0];
            const MarkerGraphEdgeId edgeId1 = outPreceders[1][i1];

            if(edgeId0 < edgeId1) {
                // edgeId0 is in outPreceders[0] but not in outPreceders[1].
                if(outPrecedersCommonCount[0][i0] >= highCommonCountThreshold) {
                    outCandidates[0].push_back(edgeId0);
                }
                ++i0;
            }

            else if(edgeId1 < edgeId0) {
                // edgeId1 is in outPreceders[1] but not in outPreceders[0].
                if(outPrecedersCommonCount[1][i1] >= highCommonCountThreshold) {
                    outCandidates[1].push_back(edgeId1);
                }
                ++i1;
            }

            else {
                // edgeId0 is in outPreceders[0] and in outPreceders[1].
                const uint64_t common0 = outPrecedersCommonCount[0][i0];
                const uint64_t common1 = outPrecedersCommonCount[1][i1];
                if(common0 >= highCommonCountThreshold and common1 <= lowCommonCountThreshold) {
                    outCandidates[0].push_back(edgeId0);
                }
                else if(common1 >= highCommonCountThreshold and common0 <= lowCommonCountThreshold) {
                    outCandidates[1].push_back(edgeId1);
                }
                ++i0;
                ++i1;
            }
        }
    }



    // Find MarkerGraphEdgeIds that are both inCandidates and outCandidates.
    for(uint64_t i0=0; i0<2; i0++) {
        for(uint64_t i1=0; i1<2; i1++) {
            mid[i0][i1].clear();
            std::set_intersection(
                inCandidates[i0].begin(), inCandidates[i0].end(),
                outCandidates[i1].begin(), outCandidates[i1].end(),
                back_inserter(mid[i0][i1]));
        }
    }
}
#endif



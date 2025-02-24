// Shasta.
#include "mode3-AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "Marker.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "weightedShuffle.hpp"
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



// Create the AnchorGraph and its vertices and edges given a vector of AnchorIds.
AnchorGraph::AnchorGraph(const Anchors& anchors, span<const AnchorId> anchorIds, uint64_t minEdgeCoverage) :
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

    // Create the edges.
    vector<AnchorId> children;
    vector<uint64_t> counts;
    for(uint64_t localAnchorId0=0; localAnchorId0<anchorIds.size(); localAnchorId0++) {
        const AnchorId anchorId0 = anchorIds[localAnchorId0];
        anchors.findChildren(anchorId0, children, counts, minEdgeCoverage);
        const uint64_t n = children.size();
        SHASTA_ASSERT(n == counts.size());
        for(uint64_t i=0; i<n; i++) {
            const AnchorId anchorId1 = children[i];
            const uint64_t coverage = counts[i];
            const AnchorId localAnchorId1 = anchors.getLocalAnchorIdInComponent(anchorId1);

            AnchorPairInfo info;
            anchors.analyzeAnchorPair(anchorId0, anchorId1, info);

            addEdgeFromLocalAnchorIds(localAnchorId0, localAnchorId1, info, coverage);
        }
    }
}



void AnchorGraph::addEdgeFromLocalAnchorIds(
    uint64_t localAnchorId0,
    uint64_t localAnchorId1,
    const AnchorPairInfo& info,
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
    const Anchors& anchors) const
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
            out << getAnchorId(v) << "\\n" << anchors[getAnchorId(v)].coverage();
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



void AnchorGraph::writeEdgeDetails(
    const string& fileName,
    const Anchors& anchors) const
{
    const AnchorGraph& anchorGraph = *this;

    ofstream csv(fileName);
    csv << "AnchorId0,AnchorId1,Coverage,Common,Common with positive offset,Offset,\n";


    // Loop over all edges.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        const AnchorGraphEdge& edge = anchorGraph[e];
        const uint64_t coverage = edge.coverage;
        const uint64_t commonCount = edge.info.common;
        SHASTA_ASSERT(coverage <= commonCount);

        const vertex_descriptor v0 = source(e, anchorGraph);
        const vertex_descriptor v1 = target(e, anchorGraph);
        const AnchorId anchorId0 = getAnchorId(v0);
        const AnchorId anchorId1 = getAnchorId(v1);

        const uint64_t commonCountPositiveOffset = anchors.countCommon(anchorId0, anchorId1, true);
        SHASTA_ASSERT(commonCountPositiveOffset <= commonCount);

        csv << anchorIdToString(anchorId0) << ",";
        csv << anchorIdToString(anchorId1) << ",";
        csv << coverage << ",";
        csv << commonCount << ",";
        csv << commonCountPositiveOffset << ",";
        csv << edge.info.offsetInBases << ",";
        csv << "\n";

    }

}



void AnchorGraph::removeNegativeOffsetEdges()
{
    AnchorGraph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];
        SHASTA_ASSERT(edge.info.common > 0);
        if(edge.info.offsetInBases < 0) {
            edgesToBeRemoved.push_back(e);
        }
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
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

// Shasta.
#include "mode3-PrimaryGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "mode3-AssemblyPath.hpp"
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



PrimaryGraph::vertex_descriptor PrimaryGraph::addVertex(MarkerGraphEdgeId edgeId)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    const vertex_descriptor v = add_vertex({edgeId}, *this);
    vertexMap.insert({edgeId, v});
    return v;
}



void PrimaryGraph::addEdgeFromVertexDescriptors(
    vertex_descriptor v0,
    vertex_descriptor v1,
    const MarkerGraphEdgePairInfo& info,
    uint64_t coverage)
{
    add_edge(v0, v1, {info, coverage}, *this);
}



void PrimaryGraph::addEdge(
    MarkerGraphEdgeId edgeId0,
    MarkerGraphEdgeId edgeId1,
    const MarkerGraphEdgePairInfo& info,
    uint64_t coverage)
{
    auto it0 = vertexMap.find(edgeId0);
    auto it1 = vertexMap.find(edgeId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    addEdgeFromVertexDescriptors(v0, v1, info, coverage);
}



// Write a PrimaryGraph in graphviz format.
void PrimaryGraph::writeGraphviz(
    const string& name,
    const PrimaryGraphDisplayOptions& options,
    const MarkerGraph& markerGraph) const
{
    ofstream out(name + ".dot");

    const PrimaryGraph& graph = *this;
    out << "digraph " << name << " {\n";

    BGL_FORALL_VERTICES(v, graph, PrimaryGraph) {
        const PrimaryGraphVertex& vertex = graph[v];
        out << vertex.edgeId;

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "[";
        }

        if(options.labels) {
            out << "label=\"";
            out << vertex.edgeId << "\\n" << markerGraph.edgeCoverage(vertex.edgeId);
            out << "\" ";
        }

        if(options.tooltips) {
            out << "tooltip=\"";
            out << vertex.edgeId;
            out << "\" ";
        }

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "]";
        }
        out << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, PrimaryGraph) {
        const PrimaryGraphEdge& edge = graph[e];
        if(not options.showNonTransitiveReductionEdges and edge.isNonTransitiveReductionEdge) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        out <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId;

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << " [";
        }

        if(edge.isNonTransitiveReductionEdge) {
            out << "style=dashed ";
        }

        if(options.tooltips) {
            out <<
                "tooltip=\"" <<
                graph[v0].edgeId << "->" <<
                graph[v1].edgeId << " ";
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



void PrimaryGraph::writeEdgeCoverageHistogram(const string& fileName) const
{
    const PrimaryGraph& primaryGraph = *this;

    // Create a histogram indexed by histogram[coverage][commonCount].
    vector< vector<uint64_t> > histogram;

    // Loop over all edges.
    BGL_FORALL_EDGES(e, primaryGraph, PrimaryGraph) {
        const PrimaryGraphEdge& edge = primaryGraph[e];
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



// Create the connected components of this PrimaryGraph,
// without changing the PrimaryGraph itself.
vector< shared_ptr<PrimaryGraph> > PrimaryGraph::createConnectedComponents(
    uint64_t minComponentSize) const
{
    const PrimaryGraph& graph = *this;

    // Compute connected components.
    // We can't use boost::connected_components because it only works
    // for undirected graphs.
    const uint64_t n = num_vertices(graph);
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    BGL_FORALL_EDGES(e, graph, PrimaryGraph) {
        const PrimaryGraph::vertex_descriptor v0 = source(e, graph);
        const PrimaryGraph::vertex_descriptor v1 = target(e, graph);
        disjointSets.union_set(v0, v1);
    }


    // Gather the vertices in each connected component.
    vector< shared_ptr<PrimaryGraph> > allComponentPointers(num_vertices(graph));
    BGL_FORALL_VERTICES(v, graph, PrimaryGraph) {
        const PrimaryGraphVertex& vertex = graph[v];
        const uint64_t componentId = disjointSets.find_set(v);
        shared_ptr<PrimaryGraph>& componentPointer = allComponentPointers[componentId];
        if(not componentPointer) {
            componentPointer = make_shared<PrimaryGraph>();
        }
        PrimaryGraph& component = *componentPointer;
        component.addVertex(vertex.edgeId);
    }


    // Gather the edges in each connected component.
    BGL_FORALL_EDGES(e, graph, PrimaryGraph) {
        const PrimaryGraph::vertex_descriptor v0 = source(e, graph);
        const PrimaryGraph::vertex_descriptor v1 = target(e, graph);
        const uint64_t edgeId0 = graph[v0].edgeId;
        const uint64_t edgeId1 = graph[v1].edgeId;
        const uint64_t componentId = disjointSets.find_set(v0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(v1));
        shared_ptr<PrimaryGraph>& componentPointer = allComponentPointers[componentId];
        SHASTA_ASSERT(componentPointer);
        PrimaryGraph& component = *componentPointer;
        component.addEdge(
            edgeId0,
            edgeId1,
            graph[e].info,
            graph[e].coverage);
    }



    // Keep only the components with at least minComponentSize vertices
    // and sort them by size.
    vector< pair<shared_ptr<PrimaryGraph>, uint64_t> > componentPointersWithSizes;
    for(const shared_ptr<PrimaryGraph>& p: allComponentPointers) {
        if(p) {
            const uint64_t componentSize = num_vertices(*p);
            if(componentSize >= minComponentSize) {
                componentPointersWithSizes.push_back({p, componentSize});
            }
        }
    }
    sort(componentPointersWithSizes.begin(), componentPointersWithSizes.end(),
        OrderPairsBySecondOnlyGreater<shared_ptr<PrimaryGraph>, uint64_t>());


    // For now return all components, including the empty ones.
    // But we want to remove the small ones and sort them by size.
    vector< shared_ptr<PrimaryGraph> > componentPointers;
    for(const auto& p: componentPointersWithSizes) {
        componentPointers.push_back(p.first);
    }
    return componentPointers;
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
void PrimaryGraph::removeCrossEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    uint64_t minOffset)
{
    PrimaryGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PrimaryGraph) {
        const PrimaryGraphEdge& edge = graph[e];

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
        BGL_FORALL_OUTEDGES(v0, e0, graph, PrimaryGraph) {
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
        BGL_FORALL_INEDGES(v1, e1, graph, PrimaryGraph) {
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
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// Remove edges for which loss = (commonCount - coverage) / commonCount > maxLoss
void PrimaryGraph::removeWeakEdges(double maxLoss)
{
    PrimaryGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PrimaryGraph) {
        const PrimaryGraphEdge& edge = graph[e];
        const double loss = double(edge.info.common - edge.coverage) / double(edge.info.common);
        if(loss > maxLoss) {
            edgesToBeRemoved.push_back(e);
        }
    }



    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
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


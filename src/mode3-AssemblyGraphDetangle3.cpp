// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detangle3.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>



void AssemblyGraph::run3(
    uint64_t /* threadCount */,
    bool /* assembleSequence */,
    bool /* debug */)
{
    cout << "AssemblyGraph::run3 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    write("A");

    detangle3();
}



void AssemblyGraph::detangle3()
{
    AssemblyGraph& assemblyGraph = *this;

    Detangle3Graph detangle3Graph(assemblyGraph);
    cout << "The Detangle3Graph has " << num_vertices(detangle3Graph) <<
        " vertices and " << num_edges(detangle3Graph) << " edges." << endl;
}



Detangle3Graph::Detangle3Graph(AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    Detangle3Graph& detangle3Graph = *this;

    createVertices();
    createEdges();
    writeGraphviz("Detangle3Graph-Complete.dot");

    // Remove edges that have either or both hasMaximumCommon flags set to false.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphEdge& edge = detangle3Graph[e];
        if(not (edge.hasMaximumCommon[0] and edge.hasMaximumCommon[1])) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle3Graph);
    }

    // Remove isolated vertices.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        if(
            (in_degree(v, detangle3Graph) == 0) and
            (out_degree(v, detangle3Graph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, detangle3Graph);
    }


    writeGraphviz("Detangle3Graph.dot");

}



// Each assembly graph edge, which must consist of a single chain, can generate
// a vertex of the Detangle3Graph.
void Detangle3Graph::createVertices()
{
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        createVertex(e);
    }
}



void Detangle3Graph::createVertex(AssemblyGraph::edge_descriptor e)
{

    Detangle3Graph& detangle3Graph = *this;

    const BubbleChain& bubbleChain = assemblyGraph[e];
    SHASTA_ASSERT(bubbleChain.isSimpleChain());
    const Chain& chain = bubbleChain.getOnlyChain();

    // A Chain must have at least one internal anchor to generate
    // a vertex of the Detangle3Graph.
    if(chain.size() <= 2) {
        return;
    }

    // Compute average coverage of the internal anchors.
    double sum = 0.;
    for(uint64_t i=1; i<chain.size()-1; i++) {
        const AnchorId anchorId = chain[i];
        sum += double(assemblyGraph.anchors[anchorId].size());
    }
    sum /= double(chain.size() - 2);
    const uint64_t coverage = uint64_t(std::round(sum));

    // Compute the total offset between the first and last internal anchors.
    uint64_t offset = 0;
    for(uint64_t i=1; i<chain.size()-2; i++) {
        const AnchorId anchorId0 = chain[i];
        const AnchorId anchorId1 = chain[i + 1];
        AnchorPairInfo info;
        assemblyGraph.anchors.analyzeAnchorPair(anchorId0, anchorId1, info);
        offset += info.offsetInBases;
    }

    // Add the vertex.
    const AnchorId firstInternalAncorId = chain[1];
    const AnchorId lastInternalAncorId = chain[chain.size() - 2];
    const vertex_descriptor v = add_vertex(
        Detangle3GraphVertex(e, firstInternalAncorId, lastInternalAncorId, coverage, offset),
        detangle3Graph);

    // Store it in the vertexMap.
    vertexMap.insert(make_pair(e, v));
}



void Detangle3Graph::createEdges()
{
    Detangle3Graph& detangle3Graph = *this;

    BGL_FORALL_VERTICES(v0, detangle3Graph, Detangle3Graph) {
        createEdges(v0, 0);
        createEdges(v0, 1);
    }



    // Set the hasMinimumOffset and hasMaximumCommon flags.
    array<vector<edge_descriptor>, 2> edges;
    vector< pair<edge_descriptor, uint64_t> > edgesWithOffsets;
    vector< pair<edge_descriptor, uint64_t> > edgesWithCommonCount;
    BGL_FORALL_VERTICES(v0, detangle3Graph, Detangle3Graph) {

        // Gather the out-edges (direction=0) and in-edges (direction=1).
        edges[0].clear();
        BGL_FORALL_OUTEDGES(v0, e, detangle3Graph, Detangle3Graph) {
            edges[0].push_back(e);
        }
        edges[1].clear();
        BGL_FORALL_INEDGES(v0, e, detangle3Graph, Detangle3Graph) {
            edges[1].push_back(e);
        }

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            if(edges[direction].empty()) {
                continue;
            }
            edgesWithOffsets.clear();
            edgesWithCommonCount.clear();
            for(const edge_descriptor e: edges[direction]) {
                const Detangle3GraphEdge& edge = detangle3Graph[e];
                const uint64_t offset = edge.info.offsetInBases;
                const uint64_t common = edge.info.common;
                edgesWithOffsets.push_back(make_pair(e, offset));
                edgesWithCommonCount.push_back(make_pair(e, common));
            }

            sort(edgesWithOffsets.begin(), edgesWithOffsets.end(),
                OrderPairsBySecondOnly<edge_descriptor, uint64_t>());
            const uint64_t minOffset = edgesWithOffsets.front().second;
            for(const auto& p: edgesWithOffsets) {
                if(p.second > minOffset) {
                    break;
                }
                detangle3Graph[p.first].hasMinimumOffset[direction] = true;
            }

            sort(edgesWithCommonCount.begin(), edgesWithCommonCount.end(),
                OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
            const uint64_t maxCommon = edgesWithCommonCount.front().second;
            for(const auto& p: edgesWithCommonCount) {
                if(p.second < maxCommon) {
                    break;
                }
                detangle3Graph[p.first].hasMaximumCommon[direction] = true;

            }
        }
    }

#if 0
        edgesWithOffsets.clear();
        edgesWithCommonCount.clear();
        BGL_FORALL_OUTEDGES(v0, e, detangle3Graph, Detangle3Graph) {
            const Detangle3GraphEdge& edge = detangle3Graph[e];
            const uint64_t offset = edge.info.offsetInBases;
            const uint64_t common = edge.info.common;
            edgesWithOffsets.push_back(make_pair(e, offset));
            edgesWithCommonCount.push_back(make_pair(e, common));
        }
        if(not edgesWithOffsets.empty()) {
            sort(edgesWithOffsets.begin(), edgesWithOffsets.end(),
                OrderPairsBySecondOnly<edge_descriptor, uint64_t>());
            const uint64_t minOffset = edgesWithOffsets.front().second;
            for(const auto& p: edgesWithOffsets) {
                if(p.second > minOffset) {
                    break;
                }
                detangle3Graph[p.first].hasMinimumOffset[0] = true;
            }

            sort(edgesWithCommonCount.begin(), edgesWithCommonCount.end(),
                OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
            const edge_descriptor eMaxCommon = edgesWithCommonCount.front().first;
            detangle3Graph[eMaxCommon].hasMaximumCommon[0] = true;
        }
#endif

}



// Create the edges starting at v0:
// If direction is 0, move forward in the AssemblyGraph.
// If direction is 1, move backward in the AssemblyGraph.
// The code is similar to Mode3Assembler::exploreReadFollowingAssemblyGraph.
void Detangle3Graph::createEdges(vertex_descriptor v0, uint64_t direction)
{
    // EXPOSE WHEN CODE STABILIZES.
    uint64_t minCommon = 4;
    double minJaccard = 0.;
    double minCorrectedJaccard = 0.8;

    Detangle3Graph& detangle3Graph = *this;
    const Detangle3GraphVertex& vertex0 = detangle3Graph[v0];
    const AssemblyGraph::edge_descriptor e0 = vertex0.e;
    const AnchorId anchorId0 = (direction == 0) ? vertex0.lastInternalAncorId : vertex0.firstInternalAncorId;

    // Do a BFS in the AssemblyGraph, starting at e0.
    // It is a forward BFS if direction is 0 and a backward BFS if direction is 1.
    std::queue<AssemblyGraph::edge_descriptor> q;
    q.push(e0);
    std::set<AssemblyGraph::edge_descriptor> s;
    s.insert(e0);



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a Chain.
        const AssemblyGraph::edge_descriptor e1 = q.front();
        q.pop();
        const BubbleChain& bubbleChain1 = assemblyGraph[e1];
        const Chain& chain1 = bubbleChain1.getOnlyChain();
        bool hasInternalAnchors = (chain1.size() > 2);


        // If the Chain has internal anchors, see if we can
        // generate an edge v0->v1 (if direction is 0) or v1->v0 (if direction is 1).
        if((e1 != e0) and hasInternalAnchors) {
            const auto it1 = vertexMap.find(e1);
            SHASTA_ASSERT(it1 != vertexMap.end());
            const vertex_descriptor v1 = it1->second;
            const Detangle3GraphVertex& vertex1 = detangle3Graph[v1];

            const AnchorId anchorId1 = (direction == 0) ? vertex1.firstInternalAncorId : vertex1.lastInternalAncorId;

            // Analyze this pair of anchors.
            AnchorPairInfo info;
            if(direction == 0) {
                assemblyGraph.anchors.analyzeAnchorPair(anchorId0, anchorId1, info);
            } else {
                assemblyGraph.anchors.analyzeAnchorPair(anchorId1, anchorId0, info);
            }



            // If good enough, generate an edge if we don't already have it.
            if(
                (info.common >= minCommon) and
                (info.jaccard() >= minJaccard) and
                (info.correctedJaccard() >= minCorrectedJaccard)
                ) {

                if(direction == 0) {
                    // Generate an edge v0->v1 if we don't already have it.
                    edge_descriptor e;
                    bool edgeExists = false;
                    tie(e, edgeExists) = boost::edge(v0, v1, detangle3Graph);
                    if(not edgeExists) {
                        add_edge(v0, v1, Detangle3GraphEdge(info), detangle3Graph);
                    }
                } else {
                    // Generate an edge v1->v0 if we don't already have it.
                    edge_descriptor e;
                    bool edgeExists = false;
                    tie(e, edgeExists) = boost::edge(v1, v0, detangle3Graph);
                    if(not edgeExists) {
                        add_edge(v1, v0, Detangle3GraphEdge(info), detangle3Graph);
                    }
                }

            }

            // If there are no common reads with anchorId0, don't continue the BFS past this point
            // (but AssemblyGraph edges still in the queue will continue to be processed.
            if(info.common == 0) {
                continue;
            }

        }


        // Enqueue the children (if direction=0) or parents (if direction=1).
        if(direction == 0) {
            const AssemblyGraph::vertex_descriptor v2 = target(e1, assemblyGraph);
            BGL_FORALL_OUTEDGES(v2, e2, assemblyGraph, AssemblyGraph) {
                if(not s.contains(e2)) {
                    q.push(e2);
                    s.insert(e2);
                }
            }
        } else {
            const AssemblyGraph::vertex_descriptor v2 = source(e1, assemblyGraph);
            BGL_FORALL_INEDGES(v2, e2, assemblyGraph, AssemblyGraph) {
                if(not s.contains(e2)) {
                    q.push(e2);
                    s.insert(e2);
                }
            }

        }
    }

}



void Detangle3Graph::writeGraphviz(const string& fileName) const
{
    const Detangle3Graph& detangle3Graph = *this;

    ofstream dot(fileName);
    dot << "digraph detangle3Graph {\n";



    // Write the vertices.
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphVertex& vertex = detangle3Graph[v];
        const AssemblyGraph::edge_descriptor e = vertex.e;
        const Chain& chain = assemblyGraph[e].getOnlyChain();

        dot << "\"" << assemblyGraph.bubbleChainStringId(e) << "\"";

        // Begin attributes.
        dot << "[";

        // Label
        dot << "label=\"" <<
            assemblyGraph.bubbleChainStringId(e) << "\\n" <<
            "n = " << chain.size() - 2 << "\\n" <<
            "c = " << vertex.coverage << "\\n" <<
            "o = " << vertex.offset <<
            "\"";

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot <<";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphEdge& edge = detangle3Graph[e];
        const vertex_descriptor v0 = source(e, detangle3Graph);
        const vertex_descriptor v1 = target(e, detangle3Graph);
        const AssemblyGraph::edge_descriptor e0 = detangle3Graph[v0].e;
        const AssemblyGraph::edge_descriptor e1 = detangle3Graph[v1].e;

        // Write the source and target vertices.
        dot << "\"" << assemblyGraph.bubbleChainStringId(e0) <<
            "\"->\"" << assemblyGraph.bubbleChainStringId(e1) << "\"";

        // Begin attributes.
        dot << "[";

        // Label.
        dot << "label=\"" <<
            edge.info.common << "\\n" <<
            std::fixed << std::setprecision(2)  << edge.info.correctedJaccard() << "\\n" <<
            edge.info.offsetInBases <<
            "\"";

        // Style.
        if(edge.hasMaximumCommon[0]) {
            if(edge.hasMaximumCommon[1]) {
                // Leave if black.
            } else {
                dot << " color=green";
            }
        } else {
            if(edge.hasMaximumCommon[1]) {
                dot << " color=blue";
            } else {
                dot << " color=red";
            }
        }

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";
    }

    dot << "}\n";
}

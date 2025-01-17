// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detangle3.hpp"
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
    createVertices();
    createEdges();
    writeGraphviz("Detangle3Graph.dot");
}



void Detangle3Graph::createVertices()
{

    Detangle3Graph& detangle3Graph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = assemblyGraph[e];
        SHASTA_ASSERT(bubbleChain.isSimpleChain());
        const Chain& chain = bubbleChain.getOnlyChain();
        if(chain.size() > 2) {
            const AnchorId firstInternalAncorId = chain[1];
            const AnchorId lastInternalAncorId = chain[chain.size() - 2];
            const vertex_descriptor v = add_vertex(
                Detangle3GraphVertex(e, firstInternalAncorId, lastInternalAncorId), detangle3Graph);
            vertexMap.insert(make_pair(e, v));
        }
    }
}



void Detangle3Graph::createEdges()
{
    Detangle3Graph& detangle3Graph = *this;

    BGL_FORALL_VERTICES(v0, detangle3Graph, Detangle3Graph) {
        createEdges(v0, 0);
        createEdges(v0, 1);
    }

}



// Create the edges starting at v0:
// If direction is 0, move forward in the AssemblyGraph.
// If direction is 1, move backward in the AssemblyGraph.
// The code is similar to Mode3Assembler::exploreReadFollowingAssemblyGraph.
void Detangle3Graph::createEdges(vertex_descriptor v0, uint64_t direction)
{
    // EXPOSE WHEN CODE STABILIZES.
    uint64_t minCommon = 6;
    double minJaccard = 0.;
    double minCorrectedJaccard = 0.9;

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

    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        const AssemblyGraph::edge_descriptor e = detangle3Graph[v].e;
        dot << "\"" << assemblyGraph.bubbleChainStringId(e) << "\";\n";
    }

    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        const vertex_descriptor v0 = source(e, detangle3Graph);
        const vertex_descriptor v1 = target(e, detangle3Graph);
        const AssemblyGraph::edge_descriptor e0 = detangle3Graph[v0].e;
        const AssemblyGraph::edge_descriptor e1 = detangle3Graph[v1].e;
        dot << "\"" << assemblyGraph.bubbleChainStringId(e0) <<
            "\"->\"" << assemblyGraph.bubbleChainStringId(e1) << "\";\n";
    }

    dot << "}\n";
}

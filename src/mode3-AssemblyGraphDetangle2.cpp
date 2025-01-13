// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detangle2.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"

// Standard library.
#include <queue>
#include <set>



// Detangling with path following.
// When this is called, all the BubbleChains must consist of a single Chain.
// This can be achieved by calling expand first.
void AssemblyGraph::detangle2()
{
    AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t chainLengthThreshold = 100000;

    cout << "AssemblyGraph::detangle2 called." << endl;
    cout << "Component " << componentId << endl;
    cout << orientedReadIds.size() << " oriented reads." << endl;
    cout << anchorIds.size() << " anchors." << endl;
    cout << num_vertices(assemblyGraph) << " vertices." << endl;
    cout << num_edges(assemblyGraph) << " edges." << endl;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());
    }

    // Make sure we have anchor annotations.
    annotateAnchors();

    // Create the Detangle2Graph.
    Detangle2Graph detangle2Graph(assemblyGraph, chainLengthThreshold);
    detangle2Graph.addEdges();
    detangle2Graph.writeGraphviz();
}



// Construct the vertices from the long Chains of the AssemblyGraph.
Detangle2Graph::Detangle2Graph(
    const AssemblyGraph& assemblyGraph,
    uint64_t chainLengthThreshold) :
    assemblyGraph(assemblyGraph)
{
    Detangle2Graph& detangle2Graph = *this;

    // Loop over all AssemblyGraph edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {

        // Get the BubbleChain corresponding to thsi AssemblyGraph edge.
        const BubbleChain& bubbleChain = assemblyGraph[e];

        // It must consist of a single Chain.
        SHASTA_ASSERT(bubbleChain.isSimpleChain());
        const Chain& chain = assemblyGraph[e].getOnlyChain();

        // If long enough, generate a vertex.
        const uint64_t offset = assemblyGraph.chainOffset(chain);
        if(offset >= chainLengthThreshold) {
            const vertex_descriptor v = add_vertex(Detangle2GraphVertex(e, offset), detangle2Graph);
            vertexMap.insert(make_pair(e, v));
        }
    }
}



void Detangle2Graph::addEdges()
{
    Detangle2Graph& detangle2Graph = *this;

    BGL_FORALL_VERTICES(v, detangle2Graph, Detangle2Graph) {
        cout << "Working on " << assemblyGraph.bubbleChainStringId(detangle2Graph[v].e) << endl;
        findForwardPath(v);
        findBackwardPath(v);
    }
}



void Detangle2Graph::findForwardPath(vertex_descriptor v0)
{
    findPath(v0, 0);
}



void Detangle2Graph::findBackwardPath(vertex_descriptor v0)
{
    findPath(v0, 1);
}



// Find a path starting at v0 and create edges if appropriate.
// Direction is 0 for a forward path and 1 for a backward path.
void Detangle2Graph::findPath(vertex_descriptor v0, uint64_t direction)
{
    Detangle2Graph& detangle2Graph = *this;

    const AssemblyGraph::edge_descriptor e0 = detangle2Graph[v0].e;
    const Chain& chain = assemblyGraph[e0].getOnlyChain();

    std::set<edge_descriptor> longChainsFound;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t seedAnchorCount = 10;
    uint64_t minCommonCount = 4;
    double minJaccard = 0;
    double minCorrectedJaccard = 0.8;

    // We do a recursive search starting from the last/first seedAnchorCount internal anchors of e0.
    std::vector<AnchorInfo> h;
    std::set<AnchorId> anchorIdsEncountered;
    if(direction == 0) {
        for(uint64_t i=0; i<seedAnchorCount; i++) {
            const uint64_t positionInChain = chain.size() - 2 - i;
            if(positionInChain == 0) {
                break;
            }
            const AnchorId anchorId = chain[positionInChain];
            h.push_back(AnchorInfo(anchorId, 0));
            std::push_heap(h.begin(), h.end());
            anchorIdsEncountered.insert(anchorId);
        }
    } else {
        for(uint64_t i=0; i<seedAnchorCount; i++) {
            const uint64_t positionInChain = i + 1;
            if(positionInChain == chain.size() - 1) {
                break;
            }
            const AnchorId anchorId = chain[positionInChain];
            h.push_back(AnchorInfo(anchorId, 0));
            std::push_heap(h.begin(), h.end());
            anchorIdsEncountered.insert(anchorId);
        }
    }



    // Main recursive loop.
    vector< pair<AnchorId, AnchorPairInfo> > anchorIds1;
    while(not h.empty()) {
        std::pop_heap(h.begin(), h.end());
        const AnchorInfo& anchorInfo0 = h.back();
        const AnchorId anchorId0 = anchorInfo0.anchorId;
        const uint64_t offset0 = anchorInfo0.offset;
        h.pop_back();

        // Path following starting at anchorId0;
        assemblyGraph.anchors.followOrientedReads(anchorId0, direction,
            minCommonCount, minJaccard, minCorrectedJaccard, anchorIds1);

        bool foundLongChain = false;
        for(const auto& p: anchorIds1) {
            const AnchorId anchorId1 = p.first;
            const AnchorPairInfo& info1 = p.second;
            if(not anchorIdsEncountered.contains(anchorId1)) {
                anchorIdsEncountered.insert(anchorId1);

                // Check the annotation for this anchor.
                const uint64_t localAnchorId1 = assemblyGraph.anchors.getLocalAnchorIdInComponent(anchorId1);
                const auto& internalChainInfo = assemblyGraph.anchorAnnotations[localAnchorId1].internalChainInfo;
                if(internalChainInfo.size() == 1)  {
                    const pair<ChainIdentifier, uint64_t>& p = internalChainInfo.front();
                    const ChainIdentifier& chainIdentifier = p.first;
                    SHASTA_ASSERT(chainIdentifier.positionInBubbleChain == 0);
                    SHASTA_ASSERT(chainIdentifier.indexInBubble == 0);
                    const edge_descriptor e1 = chainIdentifier.e;
                    if((e1 != e0) and vertexMap.contains(e1)) {
                        longChainsFound.insert(e1);
                        if(longChainsFound.size() >= 1) {
                            const auto it1 = vertexMap.find(e1);
                            SHASTA_ASSERT(it1 != vertexMap.end());
                            const vertex_descriptor v1 = it1->second;
                            if(direction == 0) {
                                add_edge(v0, v1, detangle2Graph);
                            } else {
                                add_edge(v1, v0, detangle2Graph);
                            }
                            foundLongChain = true;
                            break;
                        }
                    }
                }

                h.push_back(AnchorInfo(anchorId1, offset0 + info1.offsetInBases));
                std::push_heap(h.begin(), h.end());

            }
            if(foundLongChain) {
                break;
            }
        }
        if(foundLongChain) {
            break;
        }
    }


}



void Detangle2Graph::writeGraphviz() const
{
    const Detangle2Graph& detangle2Graph = *this;

    ofstream dot("Detangle2.dot");
    dot << "digraph Detangle2 {\n";

    BGL_FORALL_VERTICES(v, detangle2Graph, Detangle2Graph) {
        const AssemblyGraph::edge_descriptor e = detangle2Graph[v].e;
        dot << "\"" << assemblyGraph.bubbleChainStringId(e) << "\";\n";
    }

    BGL_FORALL_EDGES(e, detangle2Graph, Detangle2Graph) {
        const vertex_descriptor v0 = source(e, detangle2Graph);
        const vertex_descriptor v1 = target(e, detangle2Graph);
        const AssemblyGraph::edge_descriptor ae0 = detangle2Graph[v0].e;
        const AssemblyGraph::edge_descriptor ae1 = detangle2Graph[v1].e;
        dot << "\"" << assemblyGraph.bubbleChainStringId(ae0) << "\""
            "->\"" << assemblyGraph.bubbleChainStringId(ae1) << "\";\n";
    }

    dot << "}\n";

}

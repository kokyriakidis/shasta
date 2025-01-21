// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detangle2.hpp"
#include "findLinearChains.hpp"
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
    const uint64_t chainLengthThreshold = 50000;

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
    detangle2Graph.writeGraphviz("Detangle2Graph-Initial.dot");
    detangle2Graph.removeWeakEdges();
    detangle2Graph.writeGraphviz("Detangle2Graph-Final.dot");
    detangle2Graph.findLinearChains();
    detangle2Graph.createChains();

    assemblyGraph.write("Detangle2");
}



// Construct the vertices from the long Chains of the AssemblyGraph.
Detangle2Graph::Detangle2Graph(
    AssemblyGraph& assemblyGraph,
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

    cout << "The Detangle2Graph has " << num_vertices(detangle2Graph) << " vertices." << endl;
}



void Detangle2Graph::addEdges()
{
    Detangle2Graph& detangle2Graph = *this;

    BGL_FORALL_VERTICES(v, detangle2Graph, Detangle2Graph) {
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

    cout << "Looking for a " <<
        (direction == 0 ? "forward" : "backward") <<
        " path starting at " << assemblyGraph.bubbleChainStringId(detangle2Graph[v0].e) << endl;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t seedAnchorCount = 10;
    uint64_t minCommonCount = 4;
    double minJaccard = 0;
    double minCorrectedJaccard = 0.7;

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
    std::map<AnchorId, AnchorId> predecessorMap;
    while(not h.empty()) {

        // Get from h the AnchorId with the lowest total offset.
        std::pop_heap(h.begin(), h.end());
        const AnchorInfo& anchorInfo0 = h.back();
        const AnchorId anchorId0 = anchorInfo0.anchorId;
        const uint64_t totalOffset0 = anchorInfo0.totalOffset;
        h.pop_back();

        // Path following starting at anchorId0;
        assemblyGraph.anchors.followOrientedReads(anchorId0, direction,
            minCommonCount, minJaccard, minCorrectedJaccard, anchorIds1);

        // Check all the AnchorIds we reached.
        for(const auto& p: anchorIds1) {
            const AnchorId anchorId1 = p.first;
            const AnchorPairInfo& info1 = p.second;

            if(not anchorIdsEncountered.contains(anchorId1)) {
                anchorIdsEncountered.insert(anchorId1);
                predecessorMap.insert(make_pair(anchorId1, anchorId0));

                // Get the annotation for this anchor.
                const uint64_t localAnchorId1 = assemblyGraph.anchors.getLocalAnchorIdInComponent(anchorId1);
                const auto& internalChainInfo = assemblyGraph.anchorAnnotations[localAnchorId1].internalChainInfo;

                // If this AnchorId is internal to a single Chain,
                // we have found a good path that joins the long Chains
                // corresponding to two vertices of the Detangle2Graph.
                if(internalChainInfo.size() == 1)  {
                    const pair<ChainIdentifier, uint64_t>& p = internalChainInfo.front();
                    const ChainIdentifier& chainIdentifier = p.first;
                    SHASTA_ASSERT(chainIdentifier.positionInBubbleChain == 0);
                    SHASTA_ASSERT(chainIdentifier.indexInBubble == 0);
                    const AssemblyGraph::edge_descriptor e1 = chainIdentifier.e;
                    if((e1 != e0) and vertexMap.contains(e1)) {

                        cout << "Reached " << assemblyGraph.bubbleChainStringId(e1) << endl;
                        const auto it1 = vertexMap.find(e1);
                        SHASTA_ASSERT(it1 != vertexMap.end());
                        const vertex_descriptor v1 = it1->second;

                        vertex_descriptor u0 = v0;
                        vertex_descriptor u1 = v1;
                        if(direction == 1) {
                            std::swap(u0, u1);
                        }

                        edge_descriptor e;
                        bool edgeWasFound = false;
                        tie(e, edgeWasFound) = boost::edge(u0, u1, detangle2Graph);
                        if(not edgeWasFound) {
                            tie(e, ignore) = add_edge(u0, u1, detangle2Graph);
                        }
                        Detangle2GraphEdge& edge = detangle2Graph[e];
                        edge.found[direction] = true;

                        // Walk back the predecessor map to find the path that led us here.
                        vector<AnchorId> path;
                        AnchorId pathAnchorId = anchorId1;
                        while(true) {
                            path.push_back(pathAnchorId);
                            auto it = predecessorMap.find(pathAnchorId);
                            if(it == predecessorMap.end()) {
                                break;
                            } else
                            {
                                pathAnchorId = it->second;
                            }
                        }
                        if(direction == 0) {
                            reverse(path.begin(), path.end());
                        }
                        cout << "Found a path with " << path.size() << " anchors starting at " <<
                            anchorIdToString(path.front()) << " and ending at " <<
                            anchorIdToString(path.back()) << endl;

                        edge.paths[direction] = path;

                        return;
                    }
                }

                h.push_back(AnchorInfo(anchorId1, totalOffset0 + info1.offsetInBases));
                std::push_heap(h.begin(), h.end());

            }
        }
    }

    cout << "No path found." << endl;
}



void Detangle2Graph::writeGraphviz(const string& fileName) const
{
    const Detangle2Graph& detangle2Graph = *this;

    ofstream dot(fileName);
    dot << "digraph Detangle2 {\n";

    BGL_FORALL_VERTICES(v, detangle2Graph, Detangle2Graph) {
        const AssemblyGraph::edge_descriptor e = detangle2Graph[v].e;
        dot << "\"" << assemblyGraph.bubbleChainStringId(e) << "\";\n";
    }

    BGL_FORALL_EDGES(e, detangle2Graph, Detangle2Graph) {
        const Detangle2GraphEdge& edge = detangle2Graph[e];

        const vertex_descriptor v0 = source(e, detangle2Graph);
        const vertex_descriptor v1 = target(e, detangle2Graph);
        const AssemblyGraph::edge_descriptor ae0 = detangle2Graph[v0].e;
        const AssemblyGraph::edge_descriptor ae1 = detangle2Graph[v1].e;
        dot << "\"" << assemblyGraph.bubbleChainStringId(ae0) << "\""
            "->\"" << assemblyGraph.bubbleChainStringId(ae1) << "\"";

        if(edge.found[0] and (not edge.found[1])) {
            dot << " [style=dashed color=green]";
        }
        if(edge.found[1] and (not edge.found[0])) {
            dot << " [style=dashed color=red]";
        }

        dot << ";\n";
    }

    dot << "}\n";

}



// Remove edges found in one direction only.
void Detangle2Graph::removeWeakEdges()
{
    Detangle2Graph& detangle2Graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, detangle2Graph, Detangle2Graph) {
        const Detangle2GraphEdge& edge = detangle2Graph[e];
        if(not (edge.found[0] and edge.found[1])) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle2Graph);
    }
}



void Detangle2Graph::findLinearChains()
{
    const Detangle2Graph& detangle2Graph = *this;

    shasta::findLinearChains(detangle2Graph, 1, linearChains);

    cout << "Found " << linearChains.size() << " linear chains:" << endl;

    for(const auto& chain: linearChains) {
        SHASTA_ASSERT(not chain.empty());

        const edge_descriptor eFirst = chain.front();
        const vertex_descriptor vFirst = source(eFirst, detangle2Graph);
        cout << assemblyGraph.bubbleChainStringId(detangle2Graph[vFirst].e);

        for(const edge_descriptor e: chain) {
            const vertex_descriptor v = target(e, detangle2Graph);
            cout << " " << assemblyGraph.bubbleChainStringId(detangle2Graph[v].e);
        }
        cout << endl;
    }


    // By construction, each vertex should only appear in a single linear chain.
    // Check that this is the case.
    vector<vertex_descriptor> chainVertices;
    for(const auto& chain: linearChains) {
        SHASTA_ASSERT(not chain.empty());
        const edge_descriptor eFirst = chain.front();
        const vertex_descriptor vFirst = source(eFirst, detangle2Graph);
        chainVertices.push_back(vFirst);

        for(const edge_descriptor e: chain) {
            const vertex_descriptor v = target(e, detangle2Graph);
            chainVertices.push_back(v);
        }
    }
    sort(chainVertices.begin(), chainVertices.end());
    SHASTA_ASSERT(std::adjacent_find(chainVertices.begin(), chainVertices.end()) == chainVertices.end());

}



// Use each linear chain to create a new Chain in the assembly graph.
void Detangle2Graph::createChains()
{
    for(const vector<edge_descriptor>& linearChain: linearChains) {
        createChain(linearChain);
    }
}



void Detangle2Graph::createChain(const vector<edge_descriptor>& linearChain)
{
    Detangle2Graph& detangle2Graph = *this;

    // Gather the vertices of the linear chain.
    vector<vertex_descriptor> linearChainVertices;
    linearChainVertices.push_back(source(linearChain.front(), detangle2Graph));
    for(const edge_descriptor e: linearChain) {
        linearChainVertices.push_back(target(e, detangle2Graph));
    }

    Chain newChain;
    for(uint64_t i=0; /* Check later */; i++) {
        const bool isFirstVertex = (i == 0);
        const bool isLastVertex = (i == linearChainVertices.size() - 1);
        const vertex_descriptor v = linearChainVertices[i];
        const Detangle2GraphVertex& vertex = detangle2Graph[v];
        const Chain& oldChain = assemblyGraph[vertex.e].getOnlyChain();



        // Add the AnchorIds of the AssemblyGraph Chain corresponding to this vertex,
        // but excluding the AnchorIds in the preceding and following edges, if they exist.
        if(isFirstVertex) {

            // This is the first vertex.
            // We have to exclude the AnchorIds of the next edge.
            const edge_descriptor eNext = linearChain[i];
            const vector<AnchorId>& nextPath = detangle2Graph[eNext].shortestPath();
            const AnchorId nextAnchorId = nextPath.front();
            const uint64_t nextLocalAnchorId =
                assemblyGraph.anchors.getLocalAnchorIdInComponent(nextAnchorId);
            const auto& nextInternalChainInfo = assemblyGraph.anchorAnnotations[nextLocalAnchorId].internalChainInfo;
            SHASTA_ASSERT(nextInternalChainInfo.size() == 1);
            const pair<ChainIdentifier, uint64_t>& pNext = nextInternalChainInfo.front();
            SHASTA_ASSERT(pNext.first.e == vertex.e);
            const uint64_t nextPosition = pNext.second;

            copy(oldChain.begin(), oldChain.begin() + nextPosition, back_inserter(newChain));

        } else if(isLastVertex) {

            // This is the last vertex.
            // We have to exclude the AnchorIds of the preceding edge.
            const edge_descriptor ePrevious = linearChain[i-1];
            const vector<AnchorId>& previousPath = detangle2Graph[ePrevious].shortestPath();
            const AnchorId previousAnchorId = previousPath.back();
            const uint64_t previousLocalAnchorId =
                assemblyGraph.anchors.getLocalAnchorIdInComponent(previousAnchorId);
            const auto& previousInternalChainInfo = assemblyGraph.anchorAnnotations[previousLocalAnchorId].internalChainInfo;
            SHASTA_ASSERT(previousInternalChainInfo.size() == 1);
            const pair<ChainIdentifier, uint64_t>& pPrevious = previousInternalChainInfo.front();
            SHASTA_ASSERT(pPrevious.first.e == vertex.e);
            const uint64_t previousPosition = pPrevious.second;

            copy(oldChain.begin() + previousPosition + 1, oldChain.end(), back_inserter(newChain));


        } else {

            // This is a vertex internal to the linear chain.
            // We have to exclude AnchorIds in the preceding and following edge.
            const edge_descriptor ePrevious = linearChain[i-1];
            const edge_descriptor eNext = linearChain[i];

            const vector<AnchorId>& previousPath = detangle2Graph[ePrevious].shortestPath();
            const vector<AnchorId>& nextPath = detangle2Graph[eNext].shortestPath();

            const AnchorId previousAnchorId = previousPath.back();
            const AnchorId nextAnchorId = nextPath.front();

            const uint64_t previousLocalAnchorId =
                assemblyGraph.anchors.getLocalAnchorIdInComponent(previousAnchorId);
            const uint64_t nextLocalAnchorId =
                assemblyGraph.anchors.getLocalAnchorIdInComponent(nextAnchorId);

            const auto& previousInternalChainInfo = assemblyGraph.anchorAnnotations[previousLocalAnchorId].internalChainInfo;
            const auto& nextInternalChainInfo = assemblyGraph.anchorAnnotations[nextLocalAnchorId].internalChainInfo;

            SHASTA_ASSERT(previousInternalChainInfo.size() == 1);
            SHASTA_ASSERT(nextInternalChainInfo.size() == 1);

            const pair<ChainIdentifier, uint64_t>& pPrevious = previousInternalChainInfo.front();
            const pair<ChainIdentifier, uint64_t>& pNext = nextInternalChainInfo.front();

            SHASTA_ASSERT(pPrevious.first.e == vertex.e);
            SHASTA_ASSERT(pNext.first.e == vertex.e);

            const uint64_t previousPosition = pPrevious.second;
            const uint64_t nextPosition = pNext.second;
            SHASTA_ASSERT(nextPosition >= previousPosition);

            copy(oldChain.begin() + previousPosition + 1, oldChain.begin() + nextPosition, back_inserter(newChain));

        }



        // If this is the last vertex, there is no following edge and we are done.
        if(isLastVertex) {
            break;
        }

        // Add the AnchorIds of one of the path found in this DetangleGraphEdge.
        const edge_descriptor e = linearChain[i];
        const Detangle2GraphEdge& edge = detangle2Graph[e];
        const vector<AnchorId>& path = edge.shortestPath();
        copy(path.begin(), path.end(), back_inserter(newChain));
    }



    // Add this chain to the assembly graph.
    {
        // Find the source and target AssemblyGraph vertices.
        const vertex_descriptor v0 = linearChainVertices.front();
        const vertex_descriptor v1 = linearChainVertices.back();
        const AssemblyGraph::edge_descriptor ae0 = detangle2Graph[v0].e;
        const AssemblyGraph::edge_descriptor ae1 = detangle2Graph[v1].e;
        const AssemblyGraph::vertex_descriptor av0 = source(ae0, assemblyGraph);
        const AssemblyGraph::vertex_descriptor av1 = target(ae1, assemblyGraph);

        // Add the edge to the AssemblyGraph.
        AssemblyGraph::edge_descriptor aeNew;
        tie(aeNew, ignore) = boost::add_edge(av0, av1, assemblyGraph);
        AssemblyGraphEdge& newEdge = assemblyGraph[aeNew];
        newEdge.id = assemblyGraph.nextEdgeId++;

        // Store the Chain in the edge.
        BubbleChain& bubbleChain = newEdge;
        bubbleChain.resize(1);
        Bubble& bubble = bubbleChain.front();
        bubble.resize(1);
        bubble.front().swap(newChain);
    }

    // Now we can remove the old Chains.
    for(const vertex_descriptor v: linearChainVertices) {
        const AssemblyGraph::edge_descriptor e = detangle2Graph[v].e;
        boost::remove_edge(e, assemblyGraph);
    }
}

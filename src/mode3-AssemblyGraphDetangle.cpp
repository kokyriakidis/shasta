// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Tangle.hpp"
#include "AssemblerOptions.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



// Detangle with read following.
// This requires all bubble chains to be trivial
// (that is, to consist of just one haploid bubble).
void AssemblyGraph::detangleSuperbubblesWithReadFollowing(
    bool debug,
    uint64_t maxOffset,
    double maxLoss,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    AssemblyGraph& assemblyGraph = *this;

    // Check that all bubble chains are trivial.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());
    }

    if(debug) {
        cout << "Superbubble detangling with read following begins." << endl;
        cout << "Maximum offset to define superbubbles is " << maxOffset << "." << endl;
    }

    // Find the superbubbles.
    Superbubbles superbubbles(assemblyGraph, maxOffset);
    if(debug) {
        cout << "Found " << superbubbles.size() << " superbubbles." << endl;
    }

    // Loop over the superbubbles.
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        detangleSuperbubbleWithReadFollowing(debug, superbubbles, superbubbleId, maxOffset, maxLoss,
            lowCoverageThreshold, highCoverageThreshold);
        writeGraphviz("Z-" + to_string(superbubbleId), false);
        writeGfa("Z-" + to_string(superbubbleId));
    }

    if(debug) {
        cout << "Superbubble detangling with read following ends." << endl;
    }

}



void AssemblyGraph::detangleSuperbubbleWithReadFollowing(
    bool debug,
    const Superbubbles& superbubbles,
    uint64_t superbubbleId,
    uint64_t maxOffset,
    double maxLoss,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    if(debug) {
        cout << "detangleSuperbubbleWithReadFollowing begins for superbubble " << superbubbleId << endl;
    }
    AssemblyGraph& assemblyGraph = *this;
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

    // Use Tangle and TangleGraph to compute detangled chains.
    shared_ptr<Tangle> tangle = make_shared<Tangle>(debug, superbubbleId, *this, maxOffset, superbubble);

    // Gather AssemblyGraph edges that are both an entrance and an exit.
    vector<edge_descriptor> entranceExits;
    tangle->findEntranceExits(entranceExits);

    // If there are any entrance/exits, we have to split them and recreate the Tangle.
    if(not entranceExits.empty()) {
        // For now just give up.
        if(debug) {
            cout << "Found the following entrance/exits:";
            for(const AssemblyGraph::edge_descriptor e: entranceExits) {
                cout << " " << assemblyGraph.bubbleChainStringId(e);
            }
            cout << endl;
        }
        return;
    }

    vector< vector<AnchorId> > anchorChains;
    tangle->detangle(debug, superbubbleId, *this, maxLoss,
        lowCoverageThreshold, highCoverageThreshold,
        anchorChains);

    if(not tangle->success) {
        if(debug) {
            cout << "Could not detangle superbubble " << superbubbleId << endl;
        }
        return;
    }

    // Create a local AssemblyGraph from the detangled anchorChains.
    // This is a detangled representation of this superbubble.
    AssemblyGraph localAssemblyGraph(
        anchors,
        componentId,
        k,
        orientedReadIds,
        anchorIds,
        anchorChains,
        1,  // threadCount
        options,
        debug);
    if(debug) {
        cout << "The initial local assembly graph for this superbubble has " <<
            localAssemblyGraph.totalChainCount() << " chains." << endl;
    }

    // Cleanup steps similar to what happens for the global assembly graph.
    localAssemblyGraph.compress();
    for(uint64_t iteration=0; ; iteration ++) {
        const uint64_t oldEdgeCount = num_edges(localAssemblyGraph);
        localAssemblyGraph.cleanupBubbles(
            false,
            options.assemblyGraphOptions.bubbleCleanupMaxOffset,
            options.assemblyGraphOptions.chainTerminalCommonThreshold,
            1   // threadCount
            );
        localAssemblyGraph.compressBubbleChains();
        localAssemblyGraph.compress();
#if 0
        localAssemblyGraph.cleanupSuperbubbles(false,
            options.assemblyGraphOptions.superbubbleLengthThreshold1,
            options.assemblyGraphOptions.chainTerminalCommonThreshold);
        localAssemblyGraph.compress();
        localAssemblyGraph.removeShortSuperbubbles(false,
            options.assemblyGraphOptions.superbubbleLengthThreshold2,
            options.assemblyGraphOptions.superbubbleLengthThreshold3);
        localAssemblyGraph.compress();
        localAssemblyGraph.compressBubbleChains();
#endif

        if(num_edges(localAssemblyGraph) == oldEdgeCount) {
            break;
        }
    }
    if(debug) {
        cout << "After bubble/superbubble removal, the local assembly graph for this superbubble has " <<
            localAssemblyGraph.totalChainCount() << " chains." << endl;
    }

    // Also do a pass of vertex detangling.
    localAssemblyGraph.expand();
    localAssemblyGraph.detangleVertices(false,
        options.assemblyGraphOptions.detangleToleranceLow,
        options.assemblyGraphOptions.detangleToleranceHigh,
        true, // useBayesianModel
        options.assemblyGraphOptions.epsilon,
        options.assemblyGraphOptions.minLogP);
    localAssemblyGraph.compress();
    localAssemblyGraph.compressBubbleChains();

    if(debug) {
        cout << "The final local assembly graph for this superbubble has " <<
            localAssemblyGraph.totalChainCount() << " chains." << endl;
        localAssemblyGraph.writeGfa("Tangle-" + to_string(componentId) + "-" + to_string(superbubbleId));
        localAssemblyGraph.writeGfaExpanded("Tangle-" + to_string(componentId) + "-" + to_string(superbubbleId),
            false, false);
    }



    // Now we have to add these detangled chains to the global assembly graph,
    // in place of the superbubble. The detangled chains can be connected
    // to the entrance/exit chains after clipping their first/last AnchorId,
    // or to new vertices that we create as needed and store in this vertex map.
    std::map<AnchorId, vertex_descriptor> vertexMap;

    // Clip the entrances and get the corresponding vertices.
    for(const auto& entrance: tangle->entrances) {
        const vertex_descriptor v = cloneAndTruncateAtEnd(debug, entrance.e);
        vertexMap.insert(make_pair(entrance.anchorId, v));
    }

    // Clip the exits and get the corresponding vertices.
    for(const auto& exit: tangle->exits) {
        const vertex_descriptor v = cloneAndTruncateAtBeginning(debug, exit.e);
        vertexMap.insert(make_pair(exit.anchorId, v));
    }



    // Add a Chain for each Chain in the local AssemblyGraph.
    BGL_FORALL_EDGES(eLocal, localAssemblyGraph, AssemblyGraph) {
        const BubbleChain& localBubbleChain = localAssemblyGraph[eLocal];
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<localBubbleChain.size(); positionInBubbleChain++) {
            const Bubble& localBubble = localBubbleChain[positionInBubbleChain];
            for(uint64_t indexInBubble=0; indexInBubble<localBubble.size(); indexInBubble++) {
                const Chain& localChain = localBubble[indexInBubble];

                const AnchorId anchorId0 = localChain.front();
                const AnchorId anchorId1 = localChain.back();
                const vertex_descriptor v0 = getVertex(anchorId0, vertexMap);
                const vertex_descriptor v1 = getVertex(anchorId1, vertexMap);

                // Add the edge.
                edge_descriptor eGlobal;
                bool edgeWasAdded = false;
                tie(eGlobal, edgeWasAdded) = add_edge(v0, v1, assemblyGraph);
                SHASTA_ASSERT(edgeWasAdded);
                AssemblyGraphEdge& globalEdge = assemblyGraph[eGlobal];
                globalEdge.id = nextEdgeId++;

                // Make it a trivial BubbleChain consisting of a single Chain.
                BubbleChain& globalBubbleChain = globalEdge;
                globalBubbleChain.resize(1);
                Bubble& globalBubble = globalBubbleChain.front();
                globalBubble.resize(1);
                Chain& globalChain = globalBubble.front();

                // Build the chain.
                for(const AnchorId anchorId: localChain) {
                    globalChain.push_back(anchorId);
                }

                if(debug) {
                    cout << "Added detangled chain " << bubbleChainStringId(eGlobal) <<
                        " with " << globalChain.size() <<
                        " anchors to the global assembly graph." << endl;
                }
            }
        }
    }

    // Now we can remove all the vertices and in the superbubble
    // and all of their edges.
    for(const vertex_descriptor v: superbubble) {
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }

}

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
uint64_t AssemblyGraph::detangleSuperbubblesWithReadFollowing(
    bool debug,
    SuperbubbleCreationMethod superbubbleCreationMethod,
    uint64_t maxOffset,
    double maxLoss,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    AssemblyGraph& assemblyGraph = *this;

    // Check that all bubble chains are trivial.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(not assemblyGraph[e].isSimpleChain()) {
            cout << "Assertion failure for " << bubbleChainStringId(e) << endl;
            write("Failure");
        }
        SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());
    }

    if(debug) {
        cout << "Superbubble detangling with read following begins." << endl;
        cout << "Maximum offset to define superbubbles is " << maxOffset << "." << endl;
    }

    // Find the superbubbles.
    shared_ptr<Superbubbles> superbubbles;
    if(superbubbleCreationMethod == SuperbubbleCreationMethod::ByLength) {
        superbubbles = make_shared<Superbubbles>(assemblyGraph, maxOffset);
        // cout << "Found " << superbubbles->size() << " superbubbles for detangling." << endl;
    } else if(superbubbleCreationMethod == SuperbubbleCreationMethod::SingleEdges) {
        superbubbles = make_shared<Superbubbles>(assemblyGraph, Superbubbles::FromTangledEdges());
        // cout << "Found " << superbubbles->size() <<
        //     " superbubbles for detangling, each consisting of a single edge." << endl;
    } else {
        SHASTA_ASSERT(0);
    }
    if(debug) {
        cout << "Found " << superbubbles->size() << " superbubbles." << endl;
    }

    // Loop over the superbubbles.
    uint64_t successCount = 0;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles->size(); superbubbleId++) {
        if(detangleSuperbubbleWithReadFollowing(debug, *superbubbles, superbubbleId, maxOffset, maxLoss,
            lowCoverageThreshold, highCoverageThreshold)) {
            ++successCount;
        }
        if(debug) {
            writeGraphviz("Z-" + to_string(superbubbleId), false);
            writeGfa("Z-" + to_string(superbubbleId));
        }
    }

    if(debug) {
        cout << "Superbubble detangling with read following ends." << endl;
    }
    return successCount;
}



bool AssemblyGraph::detangleSuperbubbleWithReadFollowing(
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
        if(debug) {
            cout << "Found the following entrance/exits:";
            for(const AssemblyGraph::edge_descriptor e: entranceExits) {
                cout << " " << assemblyGraph.bubbleChainStringId(e);
            }
            cout << endl;
        }

        // If all of them are at least 4 anchors long, we can split them,
        // then recreate the tangle. Otherwise we have to give up.
        bool canDo = true;
        for(const edge_descriptor e: entranceExits) {
            const Chain& chain = assemblyGraph[e].getOnlyChain();
            if(chain.size() < 4) {
                canDo = false;
                if(debug) {
                    cout << "Chain for " << bubbleChainStringId(e) <<
                        " is only " << chain.size() << " anchors long. Cannot split." << endl;
                }
                break;
            }
        }
        if(not canDo) {
            if(debug) {
                cout << "Cannot split all entrances which are also exits. "
                    " Skipping detangling for superbubble " << superbubbleId << endl;
            }
            return false;
        }

        // Split in two all the entrances that are also exits.
        for(const edge_descriptor e: entranceExits) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const Chain& chain = assemblyGraph[e].getOnlyChain();
            auto splitPosition = chain.begin() + chain.size() / 2;
            // The AnchorId at splitPosition will be included in both chainA and chainB.
            Chain chainA;
            copy(chain.begin(), splitPosition + 1, back_inserter(chainA));
            Chain chainB;
            copy(splitPosition, chain.end(), back_inserter(chainB));
            const vertex_descriptor vMiddle = add_vertex(AssemblyGraphVertex(*splitPosition), assemblyGraph);

            // Add an edge for chainA
            edge_descriptor eA;
            tie(eA, ignore) = add_edge(v0, vMiddle, assemblyGraph);
            AssemblyGraphEdge& edgeA = assemblyGraph[eA];
            edgeA.id = nextEdgeId++;
            BubbleChain& bubbleChainA = edgeA;
            bubbleChainA.resize(1);
            Bubble& bubbleA = bubbleChainA.front();
            bubbleA.resize(1);
            bubbleA.front() = chainA;

            // Add an edge for chainB
            edge_descriptor eB;
            tie(eB, ignore) = add_edge(vMiddle, v1, assemblyGraph);
            AssemblyGraphEdge& edgeB = assemblyGraph[eB];
            edgeB.id = nextEdgeId++;
            BubbleChain& bubbleChainB = edgeB;
            bubbleChainB.resize(1);
            Bubble& bubbleB = bubbleChainB.front();
            bubbleB.resize(1);
            bubbleB.front() = chainB;

            // Remove the edge we split.
            boost::remove_edge(e, assemblyGraph);
        }

        // Now we can recreate the Tangle.
        tangle = make_shared<Tangle>(debug, superbubbleId, *this, maxOffset, superbubble);
        tangle->findEntranceExits(entranceExits);
        SHASTA_ASSERT(entranceExits.empty());
    }



    vector< vector<AnchorId> > anchorChains;
    tangle->detangle(debug, superbubbleId, *this, maxLoss,
        lowCoverageThreshold, highCoverageThreshold,
        anchorChains);

    if(not tangle->success) {
        if(debug) {
            cout << "Could not detangle superbubble " << superbubbleId << endl;
        }
        return false;
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

    if(debug) {
        cout << "Successfully complete detangling for superbubble " << superbubbleId << endl;
    }
    return true;
}

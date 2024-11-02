// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Anchor.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

const uint64_t haploidCoverageThreshold = 18;


// This function removes short hanging bubble chains
void AssemblyGraph::prune(
    bool debug,
    uint64_t pruneLength)
{

    AssemblyGraph& assemblyGraph = *this;

    // Variable to keep track of the edges to be removed. 
    // These edges are haploid bubble chains with no coverage.
    vector<edge_descriptor> edgesToRemove;

    // Loop over edges of the AssemblyGraph. Each edge corresponds to a
    // BubbleChain (not just a single Chain).
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        // e is the edge descriptor for this edge (boost graph library).

        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const uint64_t inDegree0 = in_degree(v0, assemblyGraph);
        const uint64_t outDegree1 = out_degree(v1, assemblyGraph);

        if ((inDegree0 == 0 and outDegree1 == 1) or ((inDegree0 == 1 and outDegree1 == 0))) {
            uint64_t averageOffset;
            uint64_t minOffset;
            uint64_t maxOffset;
            assemblyGraph.bubbleChainOffset(assemblyGraph[e], averageOffset, minOffset, maxOffset);
            
            if (averageOffset <= pruneLength) {
                if (debug) {
                    cout << "Edge " << bubbleChainStringId(e) << " is a short hanging segment with total length " << averageOffset << ". Removing it." << endl;
                }
                edgesToRemove.push_back(e);
            }
        }
    }

    // Remove the edges that we marked for removal
    for (const edge_descriptor& edgeToRemove : edgesToRemove) {
        boost::remove_edge(edgeToRemove, assemblyGraph);
    }

}





// This function removes the chains in the bubbles that have no coverage
void AssemblyGraph::removeChainsInBubblesWithNoInternalAnchors(bool debug)
{

    AssemblyGraph& assemblyGraph = *this;

    // Loop over edges of the AssemblyGraph. Each edge corresponds to a
    // BubbleChain (not just a single Chain).
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        // e is the edge descriptor for this edge (boost graph library).

        if (debug) {
            cout << "Working on BubbleChain " << bubbleChainStringId(e) << endl;
        }

        // Access the edge.
        AssemblyGraphEdge& edge = assemblyGraph[e];

        // Access the bubble chain.
        BubbleChain& bubbleChain = edge;  // Because BubbleChain is derived from AssemblyGraphEdge

        // If the bubble chain contains only one haploid bubble, skip it.
        if(bubbleChain.isSimpleChain()) {
            continue;
        }

        // Loop over all bubbles in the bubble chain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {

            // Access the bubble at this position.
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Get its ploidy.
            const uint64_t ploidy = bubble.size();

            if(debug) {
                cout << "Working on Bubble " << bubbleStringId(e, positionInBubbleChain) <<
                    " with ploidy " << ploidy << endl;
            }
            
            // If the bubble is haploid, skip it.
            if(bubble.isHaploid()) {
                if(debug) {
                    cout << "Skipped it because it is a haploid bubble..." << endl;
                }
                continue;
            }


            // Find the chains of the bubble that have coverage and keep them
            vector<u_int64_t> indexesToKeep;
            for(uint64_t indexInBubble = 0; indexInBubble < bubble.size(); indexInBubble++) {
                const Chain& chain = bubble[indexInBubble];
                if (chain.size() > 2) {
                    double coverage = primaryCoverage(chain);
                    indexesToKeep.push_back(indexInBubble);
                    if (debug) {
                        cout << "  Chain at index " << indexInBubble << " has coverage " << coverage << " . We keep it." << endl;
                    }
                }
                else {
                    if (debug) {
                        cout << "  Chain at index " << indexInBubble << " has only 2 anchors. We remove it." << endl;
                    }
                }
            }

            // If there are no chains with coverage, we skip this bubble
            if (indexesToKeep.empty()) {
                if (debug) {
                    cout << "No chains with coverage found. We skip this bubble." << endl;
                }
                continue;
            }

            // Create a new bubble with the chains that have coverage
            Bubble newBubble;
            for (const u_int64_t index : indexesToKeep) {
                newBubble.push_back(bubble[index]);
            }

            // Replace the original bubble in the bubble chain with the new bubble
            bubbleChain[positionInBubbleChain] = newBubble;
            
        }
    }
}




// This function fixes the bubbles that are haploid but appear as polyploid with low coverage.
void AssemblyGraph::haplotizeWronglyPolyploidBubbles(bool debug)
{

    AssemblyGraph& assemblyGraph = *this;

    // Loop over edges of the AssemblyGraph. Each edge corresponds to a
    // BubbleChain (not just a single Chain).
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        // e is the edge descriptor for this edge (boost graph library).

        if (debug) {
            cout << "Working on BubbleChain " << bubbleChainStringId(e) << endl;
        }

        // Access the edge.
        AssemblyGraphEdge& edge = assemblyGraph[e];

        // Access the bubble chain.
        BubbleChain& bubbleChain = edge;  // Because BubbleChain is derived from AssemblyGraphEdge

        // If the bubble chain contains only one haploid bubble, skip it.
        if(bubbleChain.isSimpleChain()) {
            continue;
        }

        // Loop over all bubbles in the bubble chain.
        for (uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {

            // Access the bubble at this position.
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Get its ploidy.
            const uint64_t ploidy = bubble.size();

            if(debug) {
                cout << "Working on Bubble " << bubbleStringId(e, positionInBubbleChain) <<
                    " with ploidy " << ploidy << endl;
            }
            
            // If the bubble is haploid, skip it.
            if(bubble.isHaploid()) {
                continue;
            }

            // Check the first bubble of the bubble chain            
            if (positionInBubbleChain == 0) {
                // Check coverage of the next bubble only
                bool nextBubbleIsHaploidAndLowCoverage = false;
                if (positionInBubbleChain < bubbleChain.size() - 1) {
                    const Bubble& nextBubble = bubbleChain[positionInBubbleChain + 1];
                    if(nextBubble.isHaploid()) {
                        uint64_t indexInNextHaploidBubble = 0;
                        const Chain& nextChain = nextBubble[indexInNextHaploidBubble];
                        if (nextChain.size() > 2) {
                            double nextPrimaryCoverage = primaryCoverage(nextChain);
                            if (nextPrimaryCoverage <= haploidCoverageThreshold) {
                                nextBubbleIsHaploidAndLowCoverage = true;
                            }
                        }
                    }
                }

                if(not nextBubbleIsHaploidAndLowCoverage) {
                    continue;
                }
                
                //
                // Remove the bubble by turning it into a trivial bubble chain with just two anchors
                //

                // Access the first Chain of the bubble (it does not matter which one we choose).
                const Chain& firstChainInBubble = bubble.front();

                // The Chain must have at least two "anchors", one at its beginning
                // and one at its end.
                SHASTA_ASSERT(firstChainInBubble.size() >= 2);

                // Get the first and last anchor of the bubble.
                const AnchorId firstAnchor = firstChainInBubble.front();
                const AnchorId lastAnchor = firstChainInBubble.back();

                // Check if there are common oriented reads between the first and last anchors.
                const uint64_t commonCount = anchors.countCommon(firstAnchor, lastAnchor);

                if (debug) {
                    cout << "Anchors " << firstAnchor << " " << lastAnchor <<
                            ", common count " << commonCount << endl;
                    }

                // If there are no supporting reads, we skip
                if (commonCount == 0) {
                    if (debug) {
                        cout << "No supporting reads found between the anchors." << endl;
                    }
                    continue;
                }

                // Create a new bubble with a single chain
                Bubble newBubble;
                Chain newChain;

                // Add only the first and last anchors to the new chain
                newChain.push_back(firstAnchor);
                newChain.push_back(lastAnchor);

                // Add the new chain to the new bubble
                newBubble.push_back(newChain);

                // Replace the original bubble in the bubble chain with the new bubble
                bubbleChain[positionInBubbleChain] = newBubble;
                
            } else if (positionInBubbleChain == bubbleChain.size() - 1) {
                // Check coverage of the previous bubble only
                bool previousBubbleIsHaploidAndLowCoverage = false;
                const Bubble& previousBubble = bubbleChain[positionInBubbleChain - 1];
                if(previousBubble.isHaploid()) {
                    uint64_t indexInPreviousHaploidBubble = 0;
                    const Chain& previousChain = previousBubble[indexInPreviousHaploidBubble];
                    if (previousChain.size() > 2) {
                        double previousPrimaryCoverage = primaryCoverage(previousChain);
                        if (previousPrimaryCoverage <= haploidCoverageThreshold) {
                            previousBubbleIsHaploidAndLowCoverage = true;
                        }
                    }
                }
                

                if(not previousBubbleIsHaploidAndLowCoverage) {
                    continue;
                }
                
                //
                // Remove the bubble by turning it into a trivial bubble chain with just two anchors
                //

                // Access the first Chain of the bubble (it does not matter which one we choose).
                const Chain& firstChainInBubble = bubble.front();

                // The Chain must have at least two "anchors", one at its beginning
                // and one at its end.
                SHASTA_ASSERT(firstChainInBubble.size() >= 2);

                // Get the first and last anchor of the bubble.
                const AnchorId firstAnchor = firstChainInBubble.front();
                const AnchorId lastAnchor = firstChainInBubble.back();

                // Check if there are common oriented reads between the first and last anchors.
                const uint64_t commonCount = anchors.countCommon(firstAnchor, lastAnchor);

                if (debug) {
                    cout << "Anchors " << firstAnchor << " " << lastAnchor <<
                            ", common count " << commonCount << endl;
                    }

                // If there are no supporting reads, we skip
                if (commonCount == 0) {
                    if (debug) {
                        cout << "No supporting reads found between the anchors." << endl;
                    }
                    continue;
                }

                // Create a new bubble with a single chain
                Bubble newBubble;
                Chain newChain;

                // Add only the first and last anchors to the new chain
                newChain.push_back(firstAnchor);
                newChain.push_back(lastAnchor);

                // Add the new chain to the new bubble
                newBubble.push_back(newChain);

                // Replace the original bubble in the bubble chain with the new bubble
                bubbleChain[positionInBubbleChain] = newBubble;

            } else {
                //
                // Check coverage in both previous and next bubble
                //

                // Check coverage of the next bubble
                bool nextBubbleIsHaploidAndLowCoverage = false;
                if (positionInBubbleChain < bubbleChain.size() - 1) {
                    const Bubble& nextBubble = bubbleChain[positionInBubbleChain + 1];
                    if(nextBubble.isHaploid()) {
                        uint64_t indexInNextHaploidBubble = 0;
                        const Chain& nextChain = nextBubble[indexInNextHaploidBubble];
                        if (nextChain.size() > 2) {
                            double nextPrimaryCoverage = primaryCoverage(nextChain);
                            if (nextPrimaryCoverage <= haploidCoverageThreshold) {
                                nextBubbleIsHaploidAndLowCoverage = true;
                            }
                        }
                    }
                }

                // Check coverage of the previous bubble
                bool previousBubbleIsHaploidAndLowCoverage = false;
                const Bubble& previousBubble = bubbleChain[positionInBubbleChain - 1];
                if(previousBubble.isHaploid()) {
                    uint64_t indexInPreviousHaploidBubble = 0;
                    const Chain& previousChain = previousBubble[indexInPreviousHaploidBubble];
                    if (previousChain.size() > 2) {
                        double previousPrimaryCoverage = primaryCoverage(previousChain);
                        if (previousPrimaryCoverage <= haploidCoverageThreshold) {
                            previousBubbleIsHaploidAndLowCoverage = true;
                        }
                    }
                }
                

                if(not previousBubbleIsHaploidAndLowCoverage and not nextBubbleIsHaploidAndLowCoverage) {
                    continue;
                }

                //
                // Remove the bubble by turning it into a trivial bubble chain with just two anchors
                //

                // Access the first Chain of the bubble (it does not matter which one we choose).
                const Chain& firstChainInBubble = bubble.front();

                // The Chain must have at least two "anchors", one at its beginning
                // and one at its end.
                SHASTA_ASSERT(firstChainInBubble.size() >= 2);

                // Get the first and last anchor of the bubble.
                const AnchorId firstAnchor = firstChainInBubble.front();
                const AnchorId lastAnchor = firstChainInBubble.back();

                // Check if there are common oriented reads between the first and last anchors.
                const uint64_t commonCount = anchors.countCommon(firstAnchor, lastAnchor);

                if (debug) {
                    cout << "Anchors " << firstAnchor << " " << lastAnchor <<
                            ", common count " << commonCount << endl;
                    }

                // If there are no supporting reads, we skip
                if (commonCount == 0) {
                    if (debug) {
                        cout << "No supporting reads found between the anchors." << endl;
                    }
                    continue;
                }

                // Create a new bubble with a single chain
                Bubble newBubble;
                Chain newChain;

                // Add only the first and last anchors to the new chain
                newChain.push_back(firstAnchor);
                newChain.push_back(lastAnchor);

                // Add the new chain to the new bubble
                newBubble.push_back(newChain);

                // Replace the original bubble in the bubble chain with the new bubble
                bubbleChain[positionInBubbleChain] = newBubble;

            }
            
        }
    }
}

// void AssemblyGraph::removeSimpleBubbleChainsWithNoInternalAnchors(bool debug)
// {

//     AssemblyGraph& assemblyGraph = *this;

//     // Change the type to store edge descriptors instead of AssemblyGraphEdge
//     vector<edge_descriptor> edgesToRemove;

//     // Loop over edges of the AssemblyGraph. Each edge corresponds to a
//     // BubbleChain (not just a single Chain).
//     BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
//         // e is the edge descriptor for this edge (boost graph library).

//         if (debug) {
//             cout << "Working on BubbleChain " << bubbleChainStringId(e) << endl;
//         }

//         // Access the edge.
//         AssemblyGraphEdge& edge = assemblyGraph[e];

//         // Access the bubble chain.
//         BubbleChain& bubbleChain = edge;  // Because BubbleChain is derived from AssemblyGraphEdge

//         // Process the bubble chain if it contains only one haploid bubble.
//         if(bubbleChain.isSimpleChain()) {
//             // Access the bubble at this position.
//             const Bubble& bubble = bubbleChain[0];

//             if(debug) {
//                 cout << "Working on Bubble " << bubbleStringId(e, 0) <<
//                     " with ploidy " << 1 << endl;
//             }
            
//             const Chain& chain = bubble[0];
//             if (chain.size() <= 2) {
//                 if (debug) {
//                     cout << "  Chain at index " << 0 << " has only 2 anchors" << endl;
//                     cout << "  Removing haploid bubble chain " << bubbleChainStringId(e) << " due to no coverage" << endl;
//                 }
//                 // Store the edge descriptor instead of the edge itself
//                 edgesToRemove.push_back(e);
//             }
//         }
//     }

//     // Remove the edges using the stored edge descriptors
//     for (const edge_descriptor& edgeToRemove : edgesToRemove) {
//         boost::remove_edge(edgeToRemove, assemblyGraph);
//     }
// }


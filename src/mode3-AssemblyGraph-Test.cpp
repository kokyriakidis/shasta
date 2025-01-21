// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Anchor.hpp"
#include "timestamp.hpp"
// Standard library.
#include <stack>

using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

const uint64_t haploidCoverageThreshold = 10;
const uint64_t haploidLengthThreshold = 5000000;



/*******************************************************************************
 * Remove short hanging bubble chains from the assembly graph.
 * A hanging bubble chain is one that is connected to the rest of the graph only 
 * at one end (either source or target vertex). A bubble chain is considered short
 * if its average offset is less than or equal to the pruneLength threshold.
 * 
 * The function iteratively removes short hanging chains and checks for new hanging
 * chains that may be created after each removal. This process continues until no
 * more short hanging chains are found.
 *
 * @param debug If true, enables debug output
 * @param pruneLength The length threshold below which hanging chains are removed
 */
void AssemblyGraph::prune(
    bool /* debug */,
    uint64_t pruneLength)
{
    AssemblyGraph& assemblyGraph = *this;
    std::stack<edge_descriptor> edgesToRemove;

    // Find short leaf edges
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const bool isLeaf = in_degree(v0, assemblyGraph) == 0 || out_degree(v1, assemblyGraph) == 0;
        if(!isLeaf) {
            continue;
        }

        uint64_t averageOffset, minOffset, maxOffset;
        assemblyGraph.bubbleChainOffset(assemblyGraph[e], averageOffset, minOffset, maxOffset);
        
        if (averageOffset <= pruneLength) {
            edgesToRemove.push(e);
        }
    }

    // Remove short leaves and check for new ones
    uint64_t pruneCount = 0;
    while(!edgesToRemove.empty()) {
        const edge_descriptor e = edgesToRemove.top();
        edgesToRemove.pop();

        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        // Check for new potential leaves after removing this edge
        if(in_degree(v0, assemblyGraph) == 0 && out_degree(v1, assemblyGraph) > 0) {
            if(in_degree(v1, assemblyGraph) == 1) {
                // Check children edges
                BGL_FORALL_OUTEDGES(v1, child, assemblyGraph, AssemblyGraph) {
                    if(child == e) {
                        continue;
                    }
                    
                    uint64_t avgOffset, minOffset, maxOffset;
                    assemblyGraph.bubbleChainOffset(assemblyGraph[child], avgOffset, minOffset, maxOffset);
                    
                    if(avgOffset < pruneLength && out_degree(target(child, assemblyGraph), assemblyGraph) != 0) {
                        edgesToRemove.push(child);
                    }
                }
            }
        } else if(in_degree(v0, assemblyGraph) > 0 && out_degree(v1, assemblyGraph) == 0) {
            if(out_degree(v0, assemblyGraph) == 1) {
                // Check parent edges  
                BGL_FORALL_INEDGES(v0, parent, assemblyGraph, AssemblyGraph) {
                    if(parent == e) {
                        continue;
                    }
                    
                    uint64_t avgOffset, minOffset, maxOffset; 
                    assemblyGraph.bubbleChainOffset(assemblyGraph[parent], avgOffset, minOffset, maxOffset);
                    
                    if(avgOffset < pruneLength && in_degree(source(parent, assemblyGraph), assemblyGraph) != 0) {
                        edgesToRemove.push(parent);
                    }
                }
            }
        }

        boost::remove_edge(e, assemblyGraph);
        ++pruneCount;
    }

    if(pruneCount > 0) {
        cout << "Pruned " << pruneCount << " edges." << endl;
        cout << timestamp << "mode3-AssemblyGraph::prune ends" << endl;
    }
}















// /*******************************************************************************
//  * Check if a bubble at the given position in a bubble chain is haploid and has 
//  * low coverage.
//  *
//  * A bubble is considered to have low coverage if its primary coverage is below
//  * haploidCoverageThreshold. The function checks several conditions:
//  * 1. Position must be valid within the bubble chain
//  * 2. Bubble must be haploid (contain only one chain)
//  * 3. Chain must have more than 2 anchors
//  * 4. Chain's primary coverage must be below the threshold
//  *
//  * @param assemblyGraph The assembly graph containing the bubble chain
//  * @param bubbleChain The bubble chain to check
//  * @param position Position of the bubble within the chain to check
//  * @return true if the bubble is haploid and has low coverage, false otherwise
//  */
// bool hasLowCoverageHaploidBubble(AssemblyGraph& assemblyGraph, const BubbleChain& bubbleChain, uint64_t position)
// {
//     // Check if position is valid
//     if (position >= bubbleChain.size()) {
//         return false;
//     }

//     const Bubble& bubble = bubbleChain[position];

//     // Check if bubble is haploid
//     if (!bubble.isHaploid()) {
//         return false;
//     }

//     // Get the only chain
//     const Chain& chain = bubble.front();

//     // Skip chains with only 2 anchors
//     if (chain.size() <= 2) {
//         return false;
//     }

//     // Check if coverage is low
//     const double coverage = assemblyGraph.primaryCoverage(chain);
//     return coverage < haploidCoverageThreshold;
// }



bool hasLowCoverageOrHighLengthHaploidBubble(AssemblyGraph& assemblyGraph, const BubbleChain& bubbleChain, uint64_t position)
{
    // Check if position is valid
    if (position >= bubbleChain.size()) {
        return false;
    }

    const Bubble& bubble = bubbleChain[position];

    // Check if bubble is haploid
    if (!bubble.isHaploid()) {
        return false;
    }

    // Get the only chain
    const Chain& chain = bubble.front();

    // Skip chains with only 2 anchors
    if (chain.size() <= 2) {
        return false;
    }

    uint64_t avgOffset, minOffset, maxOffset; 
    assemblyGraph.bubbleChainOffset(bubbleChain, avgOffset, minOffset, maxOffset);
    
    // Check if length is high
    const bool hasHighlength = (avgOffset>=haploidLengthThreshold);
    
    // Check if coverage is low
    const double coverage = assemblyGraph.primaryCoverage(chain);
    const bool hasLowCoverage = (coverage <= haploidCoverageThreshold);

    return (hasHighlength or hasLowCoverage);
}


/*******************************************************************************
 * Fix bubbles that appear polyploid but are likely haploid with low coverage.
 * 
 * This function identifies polyploid bubbles that are likely actually haploid
 * based on the coverage patterns of their neighboring bubbles. If a polyploid
 * bubble has a neighboring haploid bubble with low coverage (based on haploidCoverageThreshold),
 * it is considered a candidate for simplification.
 * 
 * For each candidate bubble:
 * 1. Takes the first chain's start and end anchors
 * 2. Verifies there are common reads between these anchors
 * 3. Creates a new simplified haploid bubble with just those two anchors
 * 4. Replaces the original polyploid bubble with this simplified version
 * 
 * @param debug If true, outputs debug information during processing
 */
void AssemblyGraph::haplotizeWronglyPolyploidBubbles(bool debug)
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over edges of the AssemblyGraph.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if (debug) {
            cout << "Working on BubbleChain " << bubbleChainStringId(e) << endl;
        }

        // Access the bubble chain.
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Skip if simple chain
        if(bubbleChain.isSimpleChain()) {
            continue;
        }

        // Loop over all bubbles in the chain
        for (uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            if(debug) {
                cout << "Working on Bubble " << bubbleStringId(e, positionInBubbleChain) <<
                    " with ploidy " << bubble.size() << endl;
            }
            
            // Skip haploid bubbles
            if(bubble.isHaploid()) {
                continue;
            }

            // Check if neighboring bubbles have low coverage
            bool shouldSimplify = false;
            
            if (positionInBubbleChain == 0) {
                // Check next bubble only
                shouldSimplify = hasLowCoverageOrHighLengthHaploidBubble(assemblyGraph, bubbleChain, positionInBubbleChain + 1);
            }
            else if (positionInBubbleChain == bubbleChain.size() - 1) {
                // Check previous bubble only  
                shouldSimplify = hasLowCoverageOrHighLengthHaploidBubble(assemblyGraph, bubbleChain, positionInBubbleChain - 1);
            }
            else {
                // Check both neighbors
                shouldSimplify = hasLowCoverageOrHighLengthHaploidBubble(assemblyGraph, bubbleChain, positionInBubbleChain - 1) ||
                                hasLowCoverageOrHighLengthHaploidBubble(assemblyGraph, bubbleChain, positionInBubbleChain + 1);
            }

            if (!shouldSimplify) {
                continue;
            }

            // Get first and last anchors
            const Chain& firstChain = bubble.front();
            SHASTA_ASSERT(firstChain.size() >= 2);
            const AnchorId firstAnchor = firstChain.front();
            const AnchorId lastAnchor = firstChain.back();

            // Check for common reads
            const uint64_t commonCount = anchors.countCommon(firstAnchor, lastAnchor);
            if (debug) {
                cout << "Anchors " << firstAnchor << " " << lastAnchor <<
                        ", common count " << commonCount << endl;
            }

            if (commonCount == 0) {
                if (debug) {
                    cout << "No supporting reads found between the anchors." << endl;
                }
                continue;
            }

            // Create simplified bubble with single chain
            Bubble newBubble;
            Chain newChain;
            newChain.push_back(firstAnchor);
            newChain.push_back(lastAnchor);
            newBubble.push_back(newChain);

            // Replace original bubble
            bubbleChain[positionInBubbleChain] = newBubble;
        }
    }
}












































/*******************************************************************************
 * Remove chains in bubbles that have no internal anchors (length <= 2).
 * 
 * This function processes each bubble chain in the assembly graph and filters out
 * chains that have 2 or fewer anchors, since these have no internal anchors.
 * For each non-haploid bubble in a chain:
 * 1. Creates a new filtered bubble containing only chains with >2 anchors
 * 2. Replaces the original bubble with the filtered version if any chains remain
 * 3. Skips haploid bubbles and empty filtered bubbles
 *
 * @param debug If true, prints debug information about skipped/kept chains
 */
void AssemblyGraph::removeChainsInBubblesWithNoInternalAnchors(bool debug) 
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over edges (BubbleChains) in the assembly graph
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Skip simple chains (containing only one haploid bubble)
        if(bubbleChain.isSimpleChain()) {
            continue;
        }

        // Process each bubble in the chain
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Skip haploid bubbles
            if(bubble.isHaploid()) {
                if(debug) {
                    cout << "Skipping haploid bubble " << bubbleStringId(e, positionInBubbleChain) << endl;
                }
                continue;
            }

            // Keep only chains with >2 anchors
            Bubble newBubble;
            for(uint64_t i=0; i<bubble.size(); i++) {
                const Chain& chain = bubble[i];
                if(chain.size() > 2) {
                    newBubble.push_back(chain);
                    if(debug) {
                        cout << "Keeping chain " << i << " with " << chain.size() << " anchors" << endl;
                    }
                }
            }

            // Skip if no chains remain
            if(newBubble.empty()) {
                if(debug) {
                    cout << "No chains with >2 anchors found in bubble " << bubbleStringId(e, positionInBubbleChain) << endl;
                }
                continue;
            }

            // Replace original bubble with filtered version
            bubble = newBubble;
        }
    }
}

































































/*******************************************************************************
 * Check if a vertex has any outgoing edges containing chains with internal anchors.
 * A chain has internal anchors if it contains more than 2 anchors.
 * Returns true if any outgoing edge has a chain with internal anchors in its first bubble.
 *
 * @param v The vertex descriptor to check
 * @param assemblyGraph The assembly graph containing the vertex
 * @return true if the vertex has an outgoing chain with internal anchors, false otherwise
 */
bool vertexHasOutgoingChainWithInternalAnchors(AssemblyGraph::vertex_descriptor v, AssemblyGraph& assemblyGraph) {
    BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
        // Access the edge.
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Access the bubble chain.
        const BubbleChain& bubbleChain = edge;

        // Access the first bubble of the bubble chain.
        const Bubble& bubble = bubbleChain[0];

        for(uint64_t indexInBubble = 0; indexInBubble < bubble.size(); indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            if (chain.size() > 2) {
                return true;
            }
        }
    }
    return false;
}


/*******************************************************************************
 * Check if a vertex has any incoming edges containing chains with internal anchors.
 * A chain has internal anchors if it contains more than 2 anchors.
 * Returns true if any incoming edge has a chain with internal anchors in its last bubble.
 *
 * @param v The vertex descriptor to check
 * @param assemblyGraph The assembly graph containing the vertex
 * @return true if the vertex has an incoming chain with internal anchors, false otherwise
 */
bool vertexHasIncomingChainWithInternalAnchors(AssemblyGraph::vertex_descriptor v, AssemblyGraph& assemblyGraph) {
    BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
        // Access the edge.
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Access the bubble chain.
        const BubbleChain& bubbleChain = edge;
        
        // Access the last bubble of the bubble chain.
        const Bubble& bubble = bubbleChain[bubbleChain.size() - 1];

        for(uint64_t indexInBubble = 0; indexInBubble < bubble.size(); indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            if (chain.size() > 2) {
                return true;
            }
        }
    }
    return false;
}


/*******************************************************************************
 * Check if a bubble chain is simple (contains only one haploid bubble) and has no internal anchors.
 * A chain has no internal anchors if it contains exactly 2 anchors.
 *
 * @param e The edge descriptor for the bubble chain to check
 * @param assemblyGraph The assembly graph containing the bubble chain
 * @return true if the bubble chain is simple and has no internal anchors, false otherwise
 */
bool isSimpleBubbleChainWithNoInternalAnchors(AssemblyGraph::edge_descriptor e, AssemblyGraph& assemblyGraph) {
    // Access the edge.
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    // Access the bubble chain.
    const BubbleChain& bubbleChain = edge;

    // If the bubble chain contains only one haploid bubble, skip it.
    if(not bubbleChain.isSimpleChain()) {
        return false;
    }
    
    // Access the bubble of the simple bubble chain.
    const Bubble& bubble = bubbleChain[0];

    // Access the chain of the simple bubble chain.
    const Chain& chain = bubble[0];

    if (chain.size() == 2) {
        return true;
    }

    return false;
}



/*******************************************************************************
 * Remove cross-edges in the AssemblyGraph.
 * 
 * This function identifies and removes edges that appear to be "cross-edges" - 
 * edges that likely represent incorrect connections in the assembly graph.
 * An edge Z:v0->v1 is considered a cross-edge and removed if all of the following
 * conditions are met:
 * 
 * 1. Z has no internal anchors (contains exactly 2 anchors)
 * 2. v0 (source vertex) has at least one outgoing chain with internal anchors
 * 3. v1 (target vertex) has at least one incoming chain with internal anchors
 * 
 * The rationale is that edges with no internal anchors connecting vertices that
 * have other well-supported connections (chains with internal anchors) are likely
 * to be spurious cross-connections that should be removed to improve the graph.
 * 
 * @param debug If true, outputs debug information during processing
 */
void AssemblyGraph::removeCrossEdgesInAssemblyGraph(
    bool debug)
{
    AssemblyGraph& assemblyGraph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;

    // Loop over edges of the AssemblyGraph. Each edge corresponds to a
    // BubbleChain (not just a single Chain).
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        // e is the edge descriptor for this edge (boost graph library).

        if (debug) {
            cout << "Working on BubbleChain " << bubbleChainStringId(e) << endl;
        }

        bool isPotentialCrossEdge = isSimpleBubbleChainWithNoInternalAnchors(e, assemblyGraph);
        if(not isPotentialCrossEdge) {
            continue;
        }

        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        bool hasReliableOutEdge = vertexHasOutgoingChainWithInternalAnchors(v0, assemblyGraph);
        bool hasReliableInEdge = vertexHasIncomingChainWithInternalAnchors(v1, assemblyGraph);

        if(not hasReliableOutEdge and not hasReliableInEdge) {
            continue;
        }

        edgesToBeRemoved.push_back(e);
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }
}





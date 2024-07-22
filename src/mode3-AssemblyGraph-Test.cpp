// Shasta.
#include "mode3-AssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.



// This test code uses const references for now, because it does not modify anything.
// Some of the "const" keywords will have to go away when the code starts making
// changes to the graph.
void AssemblyGraph::test()
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over edges of the AssemblyGraph. Each edge corresponds to a
    // BubbleChain (not just a single Chain).
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        // e is the edge descriptor for this edge (boost graph library).

        cout << "Working on BubbleChain " << bubbleChainStringId(e) << endl;

        // Access the edge.
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Access the bubble chain.
        const BubbleChain& bubbleChain = edge;  // Because BubbleChain is derived from AssemblyGraphEdge

        // Loop over all bubbles in the bubble chain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {

            // Access the bubble at this position.
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Get its ploidy.
            const uint64_t ploidy = bubble.size();

            cout << "Working on Bubble " << bubbleStringId(e, positionInBubbleChain) <<
                " with ploidy " << ploidy << endl;

            // Loop over all chains in this bubble.
            for(uint64_t indexInBubble=0; indexInBubble<ploidy; indexInBubble++) {

                // Access this Chain.
                const Chain& chain = bubble[indexInBubble];

                // The Chain must have at least two "anchors", one at its beginning
                // and one at its end.
                SHASTA_ASSERT(chain.size() >= 2);

                // All Chains of the bubble must have the same first and last "anchor".
                const Chain& firstChainInBubble = bubble.front();
                SHASTA_ASSERT(chain.front() == firstChainInBubble.front());
                SHASTA_ASSERT(chain.back() == firstChainInBubble.back());

                // Write its name.
                cout << "Found a chain in a non-haploid bubble: " <<
                    chainStringId(e, positionInBubbleChain, indexInBubble);

                // If it has any internal "anchors", also write its primary coverage.
                if(chain.size() > 2) {
                    cout << " with primary coverage " << primaryCoverage(chain) << endl;
                }
                cout << endl;
            }
        }
    }
}

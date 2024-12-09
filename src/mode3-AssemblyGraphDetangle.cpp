// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Tangle.hpp"
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
    }

    if(debug) {
        cout << "Superbubble detangling with read following begins." << endl;
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
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

    vector< vector<AnchorId> > anchorChains;
    Tangle tangle(debug, superbubbleId, *this, maxOffset, maxLoss,
        lowCoverageThreshold, highCoverageThreshold,
        superbubble, anchorChains);
}

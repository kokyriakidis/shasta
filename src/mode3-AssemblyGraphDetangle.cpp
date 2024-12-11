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

    // Cleanup steps similar to what happens for the global assembly graph.
    localAssemblyGraph.compress();
    for(uint64_t iteration=0; ; iteration ++) {
        const uint64_t oldEdgeCount = num_edges(localAssemblyGraph);
        localAssemblyGraph.cleanupBubbles(
            debug,
            options.assemblyGraphOptions.bubbleCleanupMaxOffset,
            options.assemblyGraphOptions.chainTerminalCommonThreshold,
            1   // threadCount
            );
        localAssemblyGraph.compressBubbleChains();
        localAssemblyGraph.compress();
        localAssemblyGraph.cleanupSuperbubbles(false,
            options.assemblyGraphOptions.superbubbleLengthThreshold1,
            options.assemblyGraphOptions.chainTerminalCommonThreshold);
        localAssemblyGraph.compress();
        localAssemblyGraph.removeShortSuperbubbles(false,
            options.assemblyGraphOptions.superbubbleLengthThreshold2,
            options.assemblyGraphOptions.superbubbleLengthThreshold3);
        localAssemblyGraph.compress();
        localAssemblyGraph.compressBubbleChains();

        if(num_edges(localAssemblyGraph) == oldEdgeCount) {
            break;
        }
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

    localAssemblyGraph.writeGfa("Tangle-" + to_string(componentId) + "-" + to_string(superbubbleId));
}

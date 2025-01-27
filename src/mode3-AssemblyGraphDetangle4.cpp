// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "AssemblerOptions.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;



void AssemblyGraph::run4(
    uint64_t threadCount,
    bool /* assembleSequence */,
    bool debug)
{
    cout << "AssemblyGraph::run4 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    write("A");

    compress();
    for(uint64_t iteration=0; ; iteration ++) {
        performanceLog << timestamp << "Iteration " << iteration <<
            " of bubble cleanup begins." << endl;
        const uint64_t cleanedUpBubbleCount = cleanupBubbles(
            debug,
            options.assemblyGraphOptions.bubbleCleanupMaxOffset,
            options.assemblyGraphOptions.chainTerminalCommonThreshold,
            threadCount);
        if(cleanedUpBubbleCount == 0) {
            break;
        }
        cout << "Cleaned up " << cleanedUpBubbleCount << " bubbles." << endl;
        compressBubbleChains();
        compress();
    }

    expand();
    // Using the Bayesian model can sometimes generate adjacent anchor pairs without common reads.
    const bool useBayesianModel = true;
    write("B");
    detangleEdges(false,
        options.assemblyGraphOptions.detangleToleranceLow,      // Use 1
        options.assemblyGraphOptions.detangleToleranceHigh,     // Use 3
        useBayesianModel,
        options.assemblyGraphOptions.epsilon,
        options.assemblyGraphOptions.minLogP);
    write("C");
    while(compressSequentialEdges());
    compressBubbleChains();

    // Expand, so we can do read following on the AssemblyGraph in the http server.
    expand();
    write("D");

    // Assemble sequence.
    assembleAllChainsMultithreaded(
        options.assemblyGraphOptions.chainTerminalCommonThreshold,
        threadCount);
    writeAssemblyDetails();
    write("Final", true);
}

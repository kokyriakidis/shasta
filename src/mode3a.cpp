// Shasta.
#include "mode3a.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "mode3a-BubbleCleaner.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "Reads.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

// Standard library.
#include "fstream.hpp"
#include <map>



Assembler::Assembler(
    uint64_t threadCount,
    uint64_t k, // Marker length
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MarkerGraph& markerGraph) :
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    markers(markers),
    markerGraph(markerGraph)
{
    // EXPOSE WHEN CODE STABILIZES.
    // These are used for detangling.
#if 0
    const uint64_t minLinkCoverage = 4;
    const uint64_t minTangleCoverage = 4;
    const uint64_t detangleIterationCount = 3;
#endif
#if 1
    // These are used to compute partial paths.
    const uint64_t segmentCoverageThreshold1ForPaths = 3;
    const uint64_t segmentCoverageThreshold2ForPaths = 6;
    const uint64_t minLinkCoverageForPaths = 3;
#endif

    // This requires the marker length k to be even.
    SHASTA_ASSERT((k % 2) == 0);

    // This does not work with RLE.
    SHASTA_ASSERT(reads.representation == 0);

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Clear all the superbubble flags in marker graph edges -
    // just in case we already ran this before.
    for( MarkerGraph::Edge& edge: markerGraph.edges) {
        edge.isSuperBubbleEdge = 0;
    }

    // Create the initial PackedMarkerGraph.
    packedMarkerGraph = make_shared<PackedMarkerGraph>(
        MappedMemoryOwner(*this), "Mode3a-PackedMarkerGraph-Initial", k, reads, markers, markerGraph, false);
    cout << "The initial PackedMarkerGraph has " <<
        packedMarkerGraph->segments.size() << " segments, " <<
        packedMarkerGraph->links.size() << " links, and " <<
        packedMarkerGraph->totalSegmentLength() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph->writeGfa();

    // Clean up the bubbles causes by errors.
    // This keeps one branch of each bubble.
    // The marker graph edges of the remaining branches are flagged as removed.
    BubbleCleaner cleaner(*packedMarkerGraph);
    cleaner.cleanup(markerGraph);
    packedMarkerGraph->remove();
    packedMarkerGraph = 0;

    // Create the final PackedMarkerGraph.
    packedMarkerGraph = make_shared<PackedMarkerGraph>(
        MappedMemoryOwner(*this), "Mode3a-PackedMarkerGraph", k, reads, markers, markerGraph, false);
    cout << "After bubble cleanup, the PackedMarkerGraph has " <<
        packedMarkerGraph->segments.size() << " segments, " <<
        packedMarkerGraph->links.size() << " links, and " <<
        packedMarkerGraph->totalSegmentLength() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph->writeSegments();
    packedMarkerGraph->writeGfa();

    // For the final PackedMarkerGraph we also need to compute the oriented reads journeys.
    packedMarkerGraph->computeJourneys(threadCount);
    packedMarkerGraph->writeJourneys();

    // Create the AssemblyGraph.
    AssemblyGraph assemblyGraph(*packedMarkerGraph);
     cout << "The initial AssemblyGraph has " <<
        num_vertices(assemblyGraph) << " segments and " <<
        num_edges(assemblyGraph) << " links." << endl;

    // Create a snapshot.
    AssemblyGraphSnapshot snapshot0(assemblyGraph, "Mode3a-AssemblyGraphSnapshot-0", *this);
    snapshot0.write();

#if 0
    // Simple detangling.
    for(uint64_t iteration=0; iteration<detangleIterationCount; iteration++) {
        assemblyGraph.simpleDetangle(minLinkCoverage, minTangleCoverage);
        cout << "After simple detangling, the AssemblyGraph has " <<
           num_vertices(assemblyGraph) << " segments and " <<
           num_edges(assemblyGraph) << " links." << endl;

        // Create a snapshot.
        AssemblyGraphSnapshot snapshot1(assemblyGraph, "Mode3a-AssemblyGraphSnapshot-" + to_string(iteration+1), *this);
        snapshot1.write();
    }
#endif

    // Follow reads to compute partial paths.
    assemblyGraph.computePartialPaths(threadCount,
        segmentCoverageThreshold1ForPaths, segmentCoverageThreshold2ForPaths, minLinkCoverageForPaths);
    assemblyGraph.writePartialPaths();
    assemblyGraph.analyzePartialPaths(threadCount);

    // Find TangledAssemblyPaths.
    assemblyGraph.computeTangledAssemblyPaths(threadCount);
}



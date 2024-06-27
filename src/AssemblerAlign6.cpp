// Shasta.
#include "Assembler.hpp"
#include "Align6.hpp"
#include "Align6Marker.hpp"
#include "KmerCounter.hpp"
#include "longestPath.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"


// Version that uses banded alignments.
void Assembler::alignOrientedReads6(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    Align6& align6)
{

    SHASTA_ASSERT(kmerCounter and kmerCounter->isAvailable());
    const array<OrientedReadId, 2> orientedReadIds = {orientedReadId0, orientedReadId1};

    // Align6 needs spans of Align6Markers for the two oriented reads.
    array<vector<Align6Marker>, 2> orientedReadAlign6MarkersVectors;
    array<span<Align6Marker>, 2> orientedReadAlign6MarkersSpans;

    if(align6Markers.isOpen()) {

        // We precomputed the Align6Markers for all oriented reads.
        // Make the spans point to them.
        for(uint64_t i=0; i<2; i++) {
            orientedReadAlign6MarkersSpans[i] = align6Markers[orientedReadIds[i].getValue()];
        }

    } else {

        // We don't have precomputed Align6Markers.
        // Conpute them here, storing in the two vectors,
        // and making the spans point to the vectors.
        for(uint64_t i=0; i<2; i++) {
            const uint64_t markerCount = markers[orientedReadIds[i].getValue()].size();
            orientedReadAlign6MarkersVectors[i].resize(markerCount);
            orientedReadAlign6MarkersSpans[i]  = span<Align6Marker>(orientedReadAlign6MarkersVectors[i].begin(), markerCount);
            getOrientedReadAlign6Markers(orientedReadIds[i], orientedReadAlign6MarkersSpans[i]);
        }
    }

    align6.align(orientedReadAlign6MarkersSpans, alignment, alignmentInfo);

}



void Assembler::computeAlign6Markers(uint64_t threadCount)
{
    // Check that we have what we need.
    checkMarkersAreOpen();
    const uint64_t orientedReadCount = markers.size();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Make space for the Align6Markers for each oriented read.
    align6Markers.createNew(largeDataName("tmp-Align6Markers"), largeDataPageSize);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        align6Markers.appendVector(markers[i].size());
    }

    // Compute them.
    const uint64_t batchSize = 10;
    setupLoadBalancing(orientedReadCount, batchSize);
    runThreads(&Assembler::computeAlign6MarkersThreadFunction, threadCount);
}



void Assembler::accessAlign6Markers()
{
    align6Markers.accessExistingReadOnly(largeDataName("tmp-Align6Markers"));
}



void Assembler::computeAlign6MarkersThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));

            // Compute the Align6Markers for this oriented read.
            getOrientedReadAlign6Markers(orientedReadId, align6Markers[i]);
        }
    }
}

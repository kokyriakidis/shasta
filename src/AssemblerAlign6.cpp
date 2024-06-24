// Shasta.
#include "Assembler.hpp"
#include "Align6.hpp"
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
    uint64_t maxSkip,
    uint64_t maxDrift,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    ostream& html)
{
#if 0
    // EXPOSE WHEN CODE STABILIZES
    const uint64_t maxLocalFrequency = 1000000000;
    const uint64_t minGlobalFrequency = 10;
    const uint64_t maxGlobalFrequency = 50;
    const double driftRateTolerance = 0.05;

    // Get the length of marker k-mers.
    const uint64_t k = assemblerInfo->k;
#endif

    SHASTA_ASSERT(kmerCounter and kmerCounter->isAvailable());

    // Get the marker KmerIds for the two oriented reads.
    array<span<KmerId>, 2> orientedReadKmerIds;
    array<vector<KmerId>, 2> orientedReadKmerIdsVectors;
    if(markerKmerIds.isOpen()) {
        if(html) {
            html << "<br>Using stored marker KmerIds." << endl;
        }

        orientedReadKmerIds[0] = markerKmerIds[orientedReadId0.getValue()];
        orientedReadKmerIds[1] = markerKmerIds[orientedReadId1.getValue()];

    } else {
        if(html) {
            html << "<br>Stored sorted marker KmerIds are not available - computing them." << endl;
        }

        // This is slower and will happen if markerKmerIds is not available.
        // Resize the vectors and make the spans point to the vectors.
        // Then call getOrientedReadMarkerKmerIds to fill them in.
        orientedReadKmerIdsVectors[0].resize(markers.size(orientedReadId0.getValue()));
        orientedReadKmerIdsVectors[1].resize(markers.size(orientedReadId1.getValue()));
        orientedReadKmerIds[0] = span<KmerId>(orientedReadKmerIdsVectors[0]);
        orientedReadKmerIds[1] = span<KmerId>(orientedReadKmerIdsVectors[1]);
        getOrientedReadMarkerKmerIds(orientedReadId0, orientedReadKmerIds[0]);
        getOrientedReadMarkerKmerIds(orientedReadId1, orientedReadKmerIds[1]);
    }



    // Align6 needs markers sorted by KmerId.
    // Use the ones from sortedMarkers if available, or else compute them.
    array<span< pair<KmerId, uint32_t> >, 2> orientedReadSortedMarkersSpans;
    array<vector< pair<KmerId, uint32_t> >, 2> orientedReadSortedMarkers;
    if(sortedMarkers.isOpen()) {
        if(html) {
            html << "<br>Using stored sorted markers." << endl;
        }

        // Make the spans point to the stored sorted markers.
        orientedReadSortedMarkersSpans[0] = sortedMarkers[orientedReadId0.getValue()];
        orientedReadSortedMarkersSpans[1] = sortedMarkers[orientedReadId1.getValue()];

    } else {

        // We don't have the sorted markers. we have to compute them here.
        if(html) {
            html << "<br>Stored sorted markers are not available - computing them." << endl;
        }

        for(uint64_t i=0; i<2; i++) {

            // Unsorted markers for this oriented read.
            const span<const KmerId>& km = orientedReadKmerIds[i];

            // Sorted markers for this oriented read.
            vector<pair<KmerId, uint32_t> >& sm = orientedReadSortedMarkers[i];

            // Copy the unsorted markers.
            const uint64_t n = km.size();
            sm.resize(n);
            for(uint64_t ordinal=0; ordinal<n; ordinal++) {
                sm[ordinal] = make_pair(km[ordinal], uint32_t(ordinal));
            }

            // Sort them.
            sort(sm.begin(), sm.end(), OrderPairsByFirstOnly<KmerId, uint32_t>());

            // Make the span point to the data in the vector.
            orientedReadSortedMarkersSpans[i] = sm;
        }
    }



    Align6(
        orientedReadSortedMarkersSpans,
        assemblerInfo->k, *kmerCounter,
        maxSkip, maxDrift,
        alignment, alignmentInfo,
        html);

}

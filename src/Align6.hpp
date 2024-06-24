#pragma once

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Align6;

    class Alignment;
    class AlignmentInfo;
    class KmerCounter;
}



class shasta::Align6 {
public:

    Align6(
        const array<span< pair<KmerId, uint32_t> >, 2>& orientedReadSortedMarkersSpans,
        uint64_t k,
        const KmerCounter&,
        uint64_t maxSkip,
        uint64_t maxDrift,
        Alignment& alignment,
        AlignmentInfo& alignmentInfo,
        ostream& html);

private:

    // References/copies for constructor arguments.
    const array<span< pair<KmerId, uint32_t> >, 2>& orientedReadSortedMarkersSpans;
    uint64_t k;
    const KmerCounter& kmerCounter;
    uint64_t maxSkip;
    uint64_t maxDrift;
    Alignment& alignment;
    AlignmentInfo& alignmentInfo;
    ostream& html;

    // CONSTANTS TO BE EXPOSED WHEN CODE STABILIZES.
    const uint64_t maxLocalFrequency = 6;
    const uint64_t minGlobalFrequency = 10;
    const uint64_t maxGlobalFrequency = 50;
    const uint64_t minLowFrequencyCount = 4;
    const double driftRateTolerance = 0.05;

    // Class used to store information about a pair of markers in the
    // two oriented reads that have the same KmerId.
    class MarkerPairInfo {
    public:
        uint32_t ordinal0;
        uint32_t ordinal1;
        KmerId kmerId;
        uint64_t globalFrequency;
        uint64_t localFrequency0;
        uint64_t localFrequency1;
        uint64_t ordinalSum() const
        {
            return ordinal0 + ordinal1;
        }
        int64_t ordinalOffset() const
        {
            return int64_t(ordinal0) - int64_t(ordinal1);
        }

        // Sort by ordinalSum.
        bool operator<(const MarkerPairInfo& that) const
        {
            return ordinalSum() < that.ordinalSum();
        }
    };



    // The offsets of the low frequency marker pairs.
    // The low frequency marker pairs consists of pairs (ordinal0, ordinal1)
    // such that KmerId(orientedReadId0, ordinal0) == KmerId(orientedReadId1, ordinal1),
    // and that satisfy the requirements on local and global frequency.
    vector<int64_t> lowFrequencyMarkerPairOffsets;
    void computeLowFrequencyMarkerPairOffsets();

    // This is called when we give up.
    // Stores an empty alignment and the corresponding AlignmentInfo.
    void storeEmptyAlignment();
};

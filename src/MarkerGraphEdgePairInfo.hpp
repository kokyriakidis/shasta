#ifndef SHASTA_MARKER_GRAPH_EDGE_PAIR_INFO_HPP
#define SHASTA_MARKER_GRAPH_EDGE_PAIR_INFO_HPP

// Shasta.
#include "invalid.hpp"

// Standard library.
#include "algorithm.hpp"
#include "cstdint.hpp"

namespace shasta {
    class MarkerGraphEdgePairInfo;
}



// Information abut the read similarity composition of two marker graph edges.
class shasta::MarkerGraphEdgePairInfo {
public:

    // The total number of OrientedReadIds in each of the edges A and B.
    uint64_t totalA = 0;
    uint64_t totalB = 0;

    // The number of common oriented reads.
    uint64_t common = 0;

    // The number of oriented reads present in A but not in B.
    uint64_t onlyA = 0;

    // The number of oriented reads present in B but not in A.
    uint64_t onlyB = 0;

    // The rest of the statistics are only valid if the number
    // of common oriented reads is not 0.

    // The estimated offset between the two edges.
    // The estimate is done using the common oriented reads.
    int64_t offsetInMarkers = invalid<int64_t>;
    int64_t offsetInBases = invalid<int64_t>;

    // The number of onlyA reads which are too short to be on edge B,
    // based on the above estimated offset.
    uint64_t onlyAShort = invalid<uint64_t>;

    // The number of onlyB reads which are too short to be on edge A,
    // based on the above estimated offset.
    uint64_t onlyBShort = invalid<uint64_t>;

    uint64_t intersectionCount() const
    {
        return common;
    }
    uint64_t unionCount() const {
        return totalA + totalB - common;
    }
    uint64_t correctedUnionCount() const
    {
        return unionCount() - onlyAShort - onlyBShort;
    }
    double jaccard() const
    {
        return double(intersectionCount()) / double(unionCount());
    }
    double correctedJaccard() const
    {
        return double(intersectionCount()) / double(correctedUnionCount());
    }

    // Order them by number of common oriented reads.
    bool operator<(const MarkerGraphEdgePairInfo& that) const
    {
        return correctedJaccard() < that.correctedJaccard();
    }
    bool operator>(const MarkerGraphEdgePairInfo& that) const
    {
        return correctedJaccard() > that.correctedJaccard();
    }

    void reverse()
    {
        swap(totalA, totalB);
        swap(onlyA, onlyB);
        swap(onlyAShort, onlyBShort);
        offsetInMarkers = - offsetInMarkers;
        offsetInBases = - offsetInBases;
    }

};

#endif

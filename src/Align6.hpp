#pragma once

// Shasta.
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Align6;

    class Align6Marker;
    class Align6Options;
    class Alignment;
    class AlignmentInfo;
    class KmerDistributionInfo;
}



class shasta::Align6 {
public:

    // All work areas get allocated here so they can be reused.
    // This reduces memory allocation activity.
    Align6(
        uint64_t k,
        uint64_t maxSkip,
        uint64_t maxDrift,
        const Align6Options& align6Options,
        const KmerDistributionInfo&,
        ostream& html);


    // The orientedReadMarkers for the two reads must be sorted by KmerId.
    // This is not checked.
    void align(
        const array<span<Align6Marker>, 2>& orientedReadMarkers,
        Alignment& alignment,
        AlignmentInfo& alignmentInfo);

    void writeGlobalFrequencyCriteria(ostream&) const;

private:

    uint64_t k;
    uint64_t maxSkip;
    uint64_t maxDrift;
    const Align6Options& align6Options;
    ostream& html;

    // The minGlobalFrequency and maxLocalFrequency actually used.
    // If align6Options.minGlobalFrequency and align6Options.maxGlobalFrequency
    // are both 0 (the default), these are obtained from the
    // KmerDistributionInfo. Otherwise, they are taken from the Align6Options.
    uint64_t minGlobalFrequency;
    uint64_t maxGlobalFrequency;

    uint64_t maxOffsetSumDelta;

    // Class used to store information about a pair of markers in the
    // two oriented reads that have the same KmerId.
    // This corresponds to an element of the alignment matrix
    // that is set.
    class MarkerPair {
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
        bool operator<(const MarkerPair& that) const
        {
            return ordinalSum() < that.ordinalSum();
        }
    };

    void clear();
    void free();

    // The number of times align was called.
    uint64_t alignCount = 0;
    void freeOrClear();

    // The ordinal offsets of the low frequency marker pairs.
    // The low frequency marker pairs consists of pairs (ordinal0, ordinal1)
    // such that KmerId(orientedReadId0, ordinal0) == KmerId(orientedReadId1, ordinal1),
    // and that satisfy the requirements on local and global frequency.
    vector<int64_t> lowFrequencyMarkerPairOffsets;
    void computeLowFrequencyMarkerPairOffsets(
        const array<span<Align6Marker>, 2>& orientedReadMarkers);

    // Histogram of ordinal offsets for the low frequency marker pairs.
    // Contains pairs(offset, frequency).
    vector< pair<int64_t, uint64_t> > offsetHistogram;
    int64_t minHistogramOffset() const
    {
        return offsetHistogram.front().first;
    }
    int64_t maxHistogramOffset() const
    {
        return offsetHistogram.back().first;
    }
    void computeOffsetHistogram();


    // The alignment matrix band that contains the most low frequency marker pairs.
    int64_t bandLow;
    int64_t bandHigh;
    int64_t bandCenter;
    void computeBand(const array<span<Align6Marker>, 2>& orientedReadMarkers);
    vector< pair<int64_t, int64_t> > derivativeChanges;

    // The marker pairs in contained in this band.
    vector<MarkerPair> inBandMarkerPairs;
    void gatherMarkerPairsInBand(const array<span<Align6Marker>, 2>& orientedReadMarkers);

    // Find out if two MarkerPairs can be connected compatibly with maxSkip and maxDrift.
    bool canBeConnected(const MarkerPair&, const MarkerPair&) const;

    // Compute connected components of the marker pairs in the band.
    // Two marker pairs belong to the same component if
    // their ordinals are compatible with maxSkip and maxDrift.
    // This vector has one entry for each entry in the
    // inBandMarkerPairs vector.
    vector<uint64_t> rank;
    vector<uint64_t> parent;
    vector<uint64_t> component;
    void computeComponents();

    // Write the marker pairs in the band and the component they belong to.
    void writeMarkerPairsInBand();

    // The best component is the one with the most low frequency markers.
    // The call to findBestComponent() returns the number of low frequency
    // markers in the best component.
    uint64_t bestComponent;
    uint64_t findBestComponent();
    vector<uint64_t> count;

    // The active marker pairs are the ones  in the best component.
    vector<MarkerPair> activeMarkerPairs;
    void gatherActiveMarkerPairs();

    // Use the active marker pairs to compute the alignment.
    vector<uint64_t> longestPath;
    void computeAlignment(
        const array<span<Align6Marker>, 2>& orientedReadMarkers,
        Alignment& alignment,
        AlignmentInfo& alignmentInfo);
    using Graph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::bidirectionalS,
        boost::no_property, boost::no_property, boost::no_property,
        boost::vecS>;
    Graph graph;
    void computeAlignment1(
        const array<span<Align6Marker>, 2>& orientedReadMarkers,
        Alignment& alignment,
        AlignmentInfo& alignmentInfo);
    vector<uint64_t> length;
    vector<uint64_t> predecessor;

    // This is called when we give up.
    // Stores an empty alignment and the corresponding AlignmentInfo.
    void storeEmptyAlignment(
        const array<span<Align6Marker>, 2>& orientedReadMarkers,
        Alignment& alignment,
        AlignmentInfo& alignmentInfo);
};

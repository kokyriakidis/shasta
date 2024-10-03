#pragma once

#include "invalid.hpp"
#include "MarkerInterval.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"

#include "cstdint.hpp"
#include "memory.hpp"
#include "span.hpp"

namespace shasta {

    class Base;
    class CompressedMarker;
    class MarkerGraph;
    class MarkerInterval;
    class Reads;


    // The main input to mode 3 assembly is a set of anchors.
    // Each anchor consists of a span of MarkerIntervals, with the following requirements:
    // - All MarkerIntervals correspond to exactly the same sequence in the corresponding oriented reads, and:
    //      * Those portions of the oriented reads are believed to be aligned.
    //      * They apear in a low number of copies in the genome being sequenced.
    // - There are no duplicate oriented reads in an anchor.
    // - The anchor coverage (number of oriented reads) is in [minPrimaryCoverage, maxPrimaryCoverage].

    namespace mode3 {

        using AnchorId = uint64_t;
        class Anchor;
        class AnchorMarkerInterval;
        class Anchors;
        class AnchorPairInfo;

    }
}



// For now the AnchorMarkerInterval is the same as MarkerInterval,
// but we want a separate class so we can add to it and remove from it.
class shasta::mode3::AnchorMarkerInterval {
public:
    OrientedReadId orientedReadId;
    array<uint32_t, 2> ordinals;

    AnchorMarkerInterval() {}

    AnchorMarkerInterval(const MarkerInterval& markerInterval) :
        orientedReadId(markerInterval.orientedReadId),
        ordinals(markerInterval.ordinals)
    {}

    bool operator==(const AnchorMarkerInterval& that) const
    {
        return
            tie(orientedReadId, ordinals[0], ordinals[1])
            ==
            tie(that.orientedReadId, that.ordinals[0], that.ordinals[1]);
    }
    bool operator<(const AnchorMarkerInterval& that) const
    {
        return
            tie(orientedReadId, ordinals[0], ordinals[1])
            <
            tie(that.orientedReadId, that.ordinals[0], that.ordinals[1]);
    }
};



// An Anchor defines a set of MarkerIntervals accessible via this public interface.
// Internals are kept private to facilitate restructuring.
class shasta::mode3::Anchor : private span<const AnchorMarkerInterval> {
private:
    using BaseClass = span<const AnchorMarkerInterval>;
public:

    Anchor(const span<const AnchorMarkerInterval>& s) : span<const AnchorMarkerInterval>(s) {}

    const AnchorMarkerInterval& operator[](uint64_t i) const
    {
        return BaseClass::operator[](i);
    }

    void check() const;

    using BaseClass::size;
    using BaseClass::begin;
    using BaseClass::end;
    using BaseClass::front;
    using BaseClass::empty;

    uint64_t coverage() const
    {
        return size();
    }

    // Return the number of common oriented reads with another Anchor.
    uint64_t countCommon(const Anchor& that) const;
};



class shasta::mode3::Anchors :
    public MultithreadedObject<Anchors>,
    public MappedMemoryOwner {
public:

    // This constructor creates the Anchors from marker graph edges.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage,
        uint64_t threadCount);

    // This constructor access existing Anchors.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    Anchor operator[](AnchorId) const;
    uint64_t size() const;

    span<const Base> anchorSequence(AnchorId) const;

    // Return the number of common oriented reads between two Anchors.
    uint64_t countCommon(AnchorId, AnchorId) const;

    // Analyze the oriented read composition of two anchors.
    void analyzeAnchorPair(AnchorId, AnchorId, AnchorPairInfo&) const;

    // Return true if the second Anchor is adjacent to the first one.
    // For precise definition see the code.
    bool areAdjacentAnchors(AnchorId, AnchorId) const;

private:
    MemoryMapped::VectorOfVectors<AnchorMarkerInterval, uint64_t> anchorMarkerIntervals;

public:
    const Reads& reads;
    uint64_t k;
    uint64_t kHalf;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

    // The sequences of the anchors.
    // We assume that marker length k is even, and anchor sequence
    // starts at the midpoint of the first marker of the anchor
    // and ends at the midpoint of the second marker of the anchor.
    // This means:
    // - If the Anchor ordinal are consecutive (as it happens
    //   when getting Anchors from marker graph edges),
    //   its sequence is guaranteed to have at least one base.
    // - If the Anchor ordinal are identical (as it may happen
    //   in a future alignment free formulation), its sequence
    //   is empty.
    MemoryMapped::VectorOfVectors<Base, uint64_t> anchorSequences;

private:
    void check() const;



    // Data and functions used when constructing the Anchors from the MarkerGraph.
    class ConstructFromMarkerGraphData {
    public:
        uint64_t minPrimaryCoverage;
        uint64_t maxPrimaryCoverage;

        const MarkerGraph* markerGraphPointer;

        // The marker intervals of the edges found by each thread.
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadMarkerIntervals;

        // The corresponding sequences
        vector< shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > > threadSequences;
    };
    ConstructFromMarkerGraphData constructFromMarkerGraphData;
    void constructFromMarkerGraphThreadFunction(uint64_t threadId);
};



// Information about the read composition similarity of two anchors A and B.
class shasta::mode3::AnchorPairInfo {
public:

    // The total number of OrientedReadIds in each of the anchors A and B.
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

    // The estimated offset between the two Anchors.
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

    void reverse()
    {
        swap(totalA, totalB);
        swap(onlyA, onlyB);
        swap(onlyAShort, onlyBShort);
        offsetInMarkers = - offsetInMarkers;
        offsetInBases = - offsetInBases;
    }

};

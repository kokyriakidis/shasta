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
    // Each anchor consists of a span of AnchorMarkerInterval, with the following requirements:
    // - All AnchorMarkerInterval correspond to exactly the same sequence in the corresponding oriented reads, and:
    //      * Those portions of the oriented reads are believed to be aligned.
    //      * They appear in a low number of copies in the genome being sequenced.
    // - There are no duplicate oriented reads in an anchor.
    // - The anchor coverage (number of oriented reads) is in [minPrimaryCoverage, maxPrimaryCoverage].

    namespace mode3 {

        using AnchorId = uint64_t;
        class Anchor;
        class AnchorMarkerInterval;
        class Anchors;
        class AnchorInfo;
        class AnchorPairInfo;

        using AnchorBaseClass = span<const AnchorMarkerInterval>;

        string anchorIdToString(AnchorId);
        AnchorId anchorIdFromString(const string&);
    }
}



// The second ordinal of the marker interval is not stored.
// It can be obtained by adding to ordinal0 the value returned
// by Anchors::ordinalOffset(AnchorId).
// This value is the same for all marker intervals of an Anchor, by construction.
// Currently this is also the same value for all Anchors and equal to 1,
// but this could change.
class shasta::mode3::AnchorMarkerInterval {
public:
    OrientedReadId orientedReadId;
    uint32_t ordinal0;
    uint32_t positionInJourney = invalid<uint32_t>;

    AnchorMarkerInterval() {}

    AnchorMarkerInterval(
        OrientedReadId orientedReadId,
        uint32_t ordinal0) :
        orientedReadId(orientedReadId),
        ordinal0(ordinal0)
    {}
};



class shasta::mode3::AnchorInfo {
public:
    uint64_t componentId = invalid<uint64_t>;
    uint64_t localAnchorIdInComponent = invalid<uint64_t>;
};



// An Anchor is a set of AnchorMarkerIntervals.
class shasta::mode3::Anchor : public AnchorBaseClass {
public:

    Anchor(const AnchorBaseClass& s) : AnchorBaseClass(s) {}

    void check() const;

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
    void writeHtml(AnchorId, AnchorId, AnchorPairInfo&, ostream&) const;

    // Return true if the second Anchor is adjacent to the first one.
    // For precise definition see the code.
    bool areAdjacentAnchors(AnchorId, AnchorId) const;

    // The offset to be added to ordinal0 of an Anchor to obtaine ordinal1.
    // Currently this is the same for all Anchors and always equal to 1,
    // but this could change.
    uint32_t ordinalOffset(AnchorId) const
    {
        return 1;
    }

private:
    MemoryMapped::VectorOfVectors<AnchorMarkerInterval, uint64_t> anchorMarkerIntervals;

public:

    // Get the first ordinal for the AnchorMarkerInterval corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInterval
    // for the requested OrientedReadId.
    uint32_t getFirstOrdinal(AnchorId, OrientedReadId) const;

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


    // The journey of each oriented read is the sequence of AnchorIds
    // encountered by the oriented read.
    MemoryMapped::VectorOfVectors<AnchorId, uint64_t> journeys;
    void computeJourneys(uint64_t threadCount);
    void writeJourneys() const;
private:
    void computeJourneysThreadFunction1(uint64_t threadId);
    void computeJourneysThreadFunction2(uint64_t threadId);
    void computeJourneysThreadFunction12(uint64_t pass);
    void computeJourneysThreadFunction3(uint64_t threadId);
    void computeJourneysThreadFunction4(uint64_t threadId);

    // Temporary storage of journeys with ordinals.
    MemoryMapped::VectorOfVectors<pair<uint64_t, uint32_t>, uint64_t> journeysWithOrdinals;

    void check() const;

public:

    // For a given AnchorId, follow the read journeys forward/backward by one step.
    // Return a vector of the AnchorIds reached in this way.
    // The count vector is the number of oriented reads each of the AnchorIds.
    void findChildren(
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count) const;
    void findParents(
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count) const;


    // In addition to the marker intervals, we also store an AnchorInfo for each Anchor.
    MemoryMapped::Vector<AnchorInfo> anchorInfos;
public:
        void storeAnchorInfo(
            AnchorId anchorId,
            uint64_t componentId,
            uint64_t localAnchorIdInComponent)
        {
            AnchorInfo& anchorInfo = anchorInfos[anchorId];
            anchorInfo.componentId =  componentId;
            anchorInfo.localAnchorIdInComponent =  localAnchorIdInComponent;
        }
        uint64_t getLocalAnchorIdInComponent(AnchorId anchorId) const
        {
            return anchorInfos[anchorId].localAnchorIdInComponent;
        }
private:



    // Data and functions used when constructing the Anchors from the MarkerGraph.
    class ConstructFromMarkerGraphData {
    public:
        uint64_t minPrimaryCoverage;
        uint64_t maxPrimaryCoverage;

        const MarkerGraph* markerGraphPointer;

        // The marker intervals of the anchors found by each thread.
        class ThreadMarkerInterval {
        public:
            OrientedReadId orientedReadId;
            uint32_t ordinal0;
        };
        vector< shared_ptr< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> > > threadMarkerIntervals;

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

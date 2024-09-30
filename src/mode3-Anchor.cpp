#include "mode3-Anchor.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;

#include <cmath>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Anchors>;



// This constructor creates the Anchor MarkerIntervals from marker graph edges.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    markers(markers)
{

    // For now copy the marker intervals from the marker graph.
    anchorMarkerIntervals.createNew(largeDataName("AnchorMarkerIntervals"), largeDataPageSize);
    for(uint64_t anchorId=0; anchorId<markerGraph.edgeMarkerIntervals.size(); anchorId++) {
        const auto v = markerGraph.edgeMarkerIntervals[anchorId];
        anchorMarkerIntervals.appendVector(v.begin(), v.end());
    }

    // Also copy the Anchor sequences from the marker graph.
    anchorSequences.createNew(largeDataName("AnchorSequences"), largeDataPageSize);
    for(uint64_t anchorId=0; anchorId<markerGraph.edgeSequence.size(); anchorId++) {
        const auto v = markerGraph.edgeSequence[anchorId];
        anchorSequences.appendVector(v.begin(), v.end());
    }

    SHASTA_ASSERT(anchorSequences.size() == anchorMarkerIntervals.size());

    check();
}



// This constructor access existing Anchors.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    markers(markers)
{
    anchorMarkerIntervals.accessExistingReadOnly(largeDataName("AnchorMarkerIntervals"));
    anchorSequences.accessExistingReadOnly(largeDataName("AnchorSequences"));
}



Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerIntervals[anchorId];
}



span<const Base> Anchors::anchorSequence(AnchorId anchorId) const
{
    return anchorSequences[anchorId];
}



uint64_t Anchors::size() const
{
    return anchorMarkerIntervals.size();
}



void Anchors::check() const
{
    const Anchors& anchors = *this;

    for(AnchorId anchorId=0; anchorId<size(); anchorId++) {
        const Anchor& anchor = anchors[anchorId];
        anchor.check();
    }
}



void Anchor::check() const
{
    const Anchor& anchor = *this;

    for(uint64_t i=1; i<size(); i++) {
        SHASTA_ASSERT(anchor[i-1].orientedReadId.getReadId() < anchor[i].orientedReadId.getReadId());
    }
}



// Return the number of common oriented reads between two Anchors.
uint64_t Anchors::countCommon(AnchorId anchorId0, AnchorId anchorId1) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    return anchor0.countCommon(anchor1);
}



// Return the number of common oriented reads with another Anchor.
// Oriented reads in each Anchor are sorted and not duplicated.
uint64_t Anchor::countCommon(const Anchor& that) const
{
    const Anchor& anchor0 = *this;
    const Anchor& anchor1 = that;

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();

    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            ++count;
            ++it0;
            ++it1;
        }
    }

    return count;
}



void Anchors::analyzeAnchorPair(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info
    ) const
{
    const Anchors& anchors = *this;

    // Prepare for the joint loop over OrientedReadIds of the two Anchors.
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];
    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    // Store the total number of OrientedReadIds on the two edges.
    info.totalA = endA - beginA;
    info.totalB = endB - beginB;


    // Joint loop over the MarkerIntervals of the two Anchors,
    // to count the common oreiented reads and compute average offsets.
    info.common = 0;
    int64_t sumMarkerOffsets = 0;
    int64_t sumTwiceBaseOffsets = 0;
    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        ++info.common;
        const OrientedReadId orientedReadId = itA->orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Compute the offset in markers.
        // This assumes an ordinal offset of 1.
        // If this changes, this code will need some changes.
        SHASTA_ASSERT(itA->ordinals[1] == itA->ordinals[0] + 1);
        SHASTA_ASSERT(itB->ordinals[1] == itB->ordinals[0] + 1);
        const uint32_t ordinalA = itA->ordinals[0];
        const uint32_t ordinalB = itB->ordinals[0];
        const int64_t markerOffset = int64_t(ordinalB) - int64_t(ordinalA);
        sumMarkerOffsets += markerOffset;

        // Compute the offset in bases.
        const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA].position);
        const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA+1].position);
        const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB].position);
        const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB+1].position);
        sumTwiceBaseOffsets -= positionA0;
        sumTwiceBaseOffsets -= positionA1;
        sumTwiceBaseOffsets += positionB0;
        sumTwiceBaseOffsets += positionB1;

        // Continue the joint loop.
        ++itA;
        ++itB;

    }
    info.onlyA = info.totalA - info.common;
    info.onlyB = info.totalB - info.common;

    // If there are no common reads, this is all we can do.
    if(info.common == 0) {
        info.offsetInMarkers = invalid<int64_t>;
        info.offsetInBases = invalid<int64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        return;
    }

    // Compute the estimated offsets.
    info.offsetInMarkers = int64_t(std::round(double(sumMarkerOffsets) / double(info.common)));
    info.offsetInBases = int64_t(0.5 * std::round(double(sumTwiceBaseOffsets) / double(info.common)));



    // Now do the joint loop again, and count the onlyA and onlyB oriented reads
    // that are too short to appear in the other edge.
    itA = beginA;
    itB = beginB;
    uint64_t onlyACheck = 0;
    uint64_t onlyBCheck = 0;
    info.onlyAShort = 0;
    info.onlyBShort = 0;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            ++onlyACheck;
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinals[0];
            const uint32_t ordinalA1 = itA->ordinals[1];
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Find the hypothetical positions of edge B, assuming the estimated base offset.
            const int64_t positionB0 = positionA0 + info.offsetInBases;
            const int64_t positionB1 = positionA1 + info.offsetInBases;

            // If this ends up outside the read, this counts as onlyAShort.
            if(positionB0 < 0 or positionB1 >= lengthInBases) {
                ++info.onlyAShort;
            }

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            ++onlyBCheck;
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinals[0];
            const uint32_t ordinalB1 = itB->ordinals[1];
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA0 = positionB0 - info.offsetInBases;
            const int64_t positionA1 = positionB1 - info.offsetInBases;

            // If this ends up outside the read, this counts as onlyBShort.
            if(positionA0 < 0 or positionA1 >= lengthInBases) {
                ++info.onlyBShort;
            }

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both edges. In this loop, we
            // don't need to do anything.
            ++itA;
            ++itB;
        }
    }
    SHASTA_ASSERT(onlyACheck == info.onlyA);
    SHASTA_ASSERT(onlyBCheck == info.onlyB);
}



// This is only needed for checking during the transition.
void AnchorPairInfo::checkIdentical(const MarkerGraphEdgePairInfo& info) const
{
    // The two have identical memory layouts so we can just check the bytes.
    const uint64_t n = sizeof(AnchorPairInfo);
    static_assert(sizeof(MarkerGraphEdgePairInfo) == n);
    SHASTA_ASSERT(std::memcmp(this, &info, n) == 0);

}



// Return true if the second Anchor is adjacent to the first one,
// as seen by at least one of the common oriented reads.
// THIS CRITERION MAY BE TOO LOSE.
bool Anchors::areAdjacentAnchors(AnchorId anchorId0, AnchorId anchorId1) const
{
    const auto markerIntervals0 = anchorMarkerIntervals[anchorId0];
    const auto markerIntervals1 = anchorMarkerIntervals[anchorId1];

    // Joint loop over the marker intervals.
    auto it0 = markerIntervals0.begin();
    auto it1 = markerIntervals1.begin();
    const auto end0 = markerIntervals0.end();
    const auto end1 = markerIntervals1.end();
    while(it0 != end0 and it1 != end1) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;

        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            if(it0->ordinals[1] == it1->ordinals[0]) {
                return true;
            }
            ++it0;
            ++it1;
        }
    }

    return false;

}

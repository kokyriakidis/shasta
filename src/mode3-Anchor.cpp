#include "mode3-Anchor.hpp"
#include "deduplicate.hpp"
#include "Marker.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

#include <cmath>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Anchors>;



// This constructor accesses existing Anchors.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    SHASTA_ASSERT((k %2) == 0);
    kHalf = k / 2;

    anchorMarkerIntervals.accessExistingReadOnly(largeDataName("AnchorMarkerIntervals"));
    anchorSequences.accessExistingReadOnly(largeDataName("AnchorSequences"));
    anchorInfos.accessExistingReadWrite(largeDataName("AnchorInfos"));
    journeys.accessExistingReadOnly(largeDataName("Journeys"));
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
        SHASTA_ASSERT(ordinalOffset(anchorIdA) == 1);
        SHASTA_ASSERT(ordinalOffset(anchorIdB) == 1);
        const uint32_t ordinalA = itA->ordinal0;
        const uint32_t ordinalB = itB->ordinal0;
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
            const uint32_t ordinalA0 = itA->ordinal0;
            const uint32_t ordinalA1 = ordinalA0 + ordinalOffset(anchorIdA);
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
            const uint32_t ordinalB0 = itB->ordinal0;
            const uint32_t ordinalB1 = ordinalB0 + ordinalOffset(anchorIdB);
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
            if(it0->ordinal0 + ordinalOffset(anchorId0) == it1->ordinal0) {
                return true;
            }
            ++it0;
            ++it1;
        }
    }

    return false;

}


// Anchors are numbered such that each pair of reverse complemented
// AnchorIds are numbered (n, n+1), where n is even, n = 2*m.
// We represent an AnchorId as a string as follows:
// - AnchorId n is represented as m+
// - AnchorId n+1 is represented as m-
// For example, the reverse complemented pair (150, 151) is represented as (75+, 75-).

string shasta::mode3::anchorIdToString(AnchorId n)
{
    std::ostringstream s;

    const AnchorId m = (n >> 1);
    s << m;

    if(n & 1) {
        s << "-";
    } else {
        s << "+";
    }

    return s.str();
}



AnchorId shasta::mode3::anchorIdFromString(const string& s)
{
    SHASTA_ASSERT(s.size() >= 2);

    const char cLast = s.back();

    uint64_t lastBit;
    if(cLast == '+') {
        lastBit = 0;
    } else if(cLast == '-') {
        lastBit = 1;
    } else {
        SHASTA_ASSERT(0);
    }

    const uint64_t m = std::stoul(s.substr(0, s.size()-1));

    return 2* m + lastBit;
}



void Anchors::computeJourneys(uint64_t threadCount)
{
    performanceLog << timestamp << "Anchors::computeJourneys begins." << endl;

    const uint64_t orientedReadCount = 2 * reads.readCount();

    // Pass1: make space for the journeysWithOrdinals.
    journeysWithOrdinals.createNew(largeDataName("tmp-JourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(orientedReadCount);
    const uint64_t anchorBatchCount = 1000;
    setupLoadBalancing(size(), anchorBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithOrdinals.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(size(), anchorBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3:sort the journeysWithOrdinals and make space for the journeys
    journeys.createNew(largeDataName("Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    const uint64_t orientedReadBatchCount = 1000;
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction3, threadCount);

    // Pass 4: copy the sorted journeysWithOrdinals to the journeys.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction4, threadCount);
    journeys.endPass2(false, true);

    journeysWithOrdinals.remove();

    performanceLog << timestamp << "Anchors::computeJourneys ends." << endl;

    writeJourneys();
}



void Anchors::computeJourneysThreadFunction1(uint64_t /* threadId */)
{
    computeJourneysThreadFunction12(1);
}


void Anchors::writeJourneys() const
{
    const ReadId readCount = reads.readCount();

    ofstream csv("Journeys.csv");
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            csv << orientedReadId << ",";
            const auto journey = journeys[orientedReadId.getValue()];
            for(const AnchorId anchorId: journey) {
                csv << anchorIdToString(anchorId) << ",";
            }
            csv << "\n";
        }
    }
}



void Anchors::computeJourneysThreadFunction2(uint64_t /* threadId */)
{
    computeJourneysThreadFunction12(2);
}



void Anchors::computeJourneysThreadFunction12(uint64_t pass)
{
    Anchors& anchors = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all AnchorIds in this batch.
        for(AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            Anchor anchor = anchors[anchorId];

            // Loop over the marker intervals of this Anchor.
            for(const auto& anchorMarkerInterval: anchor) {
                const auto orientedReadIdValue = anchorMarkerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithOrdinals.storeMultithreaded(
                        orientedReadIdValue, {anchorId, anchorMarkerInterval.ordinal0});
                }

            }
        }

    }
}




void Anchors::computeJourneysThreadFunction3(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            auto v = journeysWithOrdinals[orientedReadValue];
            sort(v.begin(), v.end(), OrderPairsBySecondOnly<uint64_t, uint32_t>());
            journeys.incrementCountMultithreaded(orientedReadValue, v.size());
        }
    }
}



void Anchors::computeJourneysThreadFunction4(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orientedReadValue));

            // Copy the journeysWithOrdinals to the journeys.
            const auto v = journeysWithOrdinals[orientedReadValue];
            const auto journey = journeys[orientedReadValue];
            SHASTA_ASSERT(journey.size() == v.size());
            for(uint64_t i=0; i<v.size(); i++) {
                journey[i] = v[i].first;
            }

            // Store journey information for this oriented read in the marker interval.
            for(uint64_t position=0; position<journey.size(); position++) {
                const AnchorId anchorId = journey[position];
                span<AnchorMarkerInterval> markerIntervals = anchorMarkerIntervals[anchorId];
                bool found = false;
                for(AnchorMarkerInterval& markerInterval: markerIntervals) {
                    if(markerInterval.orientedReadId == orientedReadId) {
                        markerInterval.positionInJourney = uint32_t(position);
                        found = true;
                        break;
                    }
                }
                SHASTA_ASSERT(found);
            }
        }
    }
}



// For a given AnchorId, follow the read journeys forward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findChildren(
    AnchorId anchorId,
    vector<AnchorId>& children,
    vector<uint64_t>& count) const
{
    children.clear();
    for(const auto& markerInterval: anchorMarkerIntervals[anchorId]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto journey = journeys[orientedReadId.getValue()];
        const uint64_t position = markerInterval.positionInJourney;
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size()) {
            const AnchorId nextAnchorId = journey[nextPosition];
            children.push_back(nextAnchorId);
        }
    }

    deduplicateAndCount(children, count);
}

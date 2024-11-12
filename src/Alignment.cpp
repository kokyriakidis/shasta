#include "Alignment.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
using namespace shasta;

uint32_t Alignment::maxSkip() const
{
    uint32_t returnValue = 0;
    for(uint64_t i=1; i<ordinals.size(); i++) {
        const uint32_t skip = max(
            abs(int32_t(ordinals[i][0]) - int32_t(ordinals[i-1][0])),
            abs(int32_t(ordinals[i][1]) - int32_t(ordinals[i-1][1]))
            );
        returnValue = max(returnValue, skip);
    }
    return returnValue;
}



uint32_t Alignment::maxDrift() const
{
    uint32_t returnValue = 0;
    for(uint64_t i=1; i<ordinals.size(); i++) {
        const int32_t offset =  int32_t(ordinals[i][0]) - int32_t(ordinals[i][1]);
        const int32_t previousOffset =  int32_t(ordinals[i-1][0]) - int32_t(ordinals[i-1][1]);
        const uint32_t drift = abs(offset - previousOffset);
        returnValue = max(returnValue, drift);
    }
    return returnValue;
}



void Alignment::swap()
{
    for(auto& o: ordinals) {
        std::swap(o[0], o[1]);
    }
}



void Alignment::reverseComplement(
    uint32_t markerCount0,
    uint32_t markerCount1)
{
    for(auto& o: ordinals) {
        o[0] = markerCount0 - 1 - o[0];
        o[1] = markerCount1 - 1 - o[1];
    }
    reverse(ordinals.begin(), ordinals.end());
}



void Alignment::checkStrictlyIncreasing() const
{
    for(uint64_t i=1; i<ordinals.size(); i++) {
        const auto& previous = ordinals[i-1];
        const auto& current = ordinals[i];
        SHASTA_ASSERT(previous[0] < current[0]);
        SHASTA_ASSERT(previous[1] < current[1]);
    }
}



void AlignmentInfo::create(
    const Alignment& alignment,
    const array<uint32_t, 2>& markerCounts)
{
    // Store the number of markers in the alignment.
    markerCount = uint32_t(alignment.ordinals.size());

    // Store alignment information for each of the two oriented reads.
    for(size_t i=0; i<2; i++) {
        data[i] = Data(
            markerCounts[i],
            (markerCount == 0) ? 0 : alignment.ordinals.front()[i],
            (markerCount == 0) ? 0 : alignment.ordinals.back()[i]);
        data[i].check();
    }



    // Compute minimum, maximum, and average ordinal offset.
    // Also compute maxSkip and maxDrift.
    minOrdinalOffset = std::numeric_limits<int32_t>::max();
    maxOrdinalOffset = std::numeric_limits<int32_t>::min();
    maxSkip = 0;
    maxDrift = 0;
    double sum = 0.;
    for(uint64_t i=0; i<alignment.ordinals.size(); i++) {
        auto& ordinals = alignment.ordinals[i];
        const int32_t offset =
            int32_t(ordinals[0]) - int32_t(ordinals[1]);
        minOrdinalOffset = min(minOrdinalOffset, offset);
        maxOrdinalOffset = max(maxOrdinalOffset, offset);
        sum += double(offset);

        if(i != 0) {
            auto& previousOrdinals = alignment.ordinals[i-1];
            maxSkip = max(maxSkip, uint32_t(abs(
                int32_t(ordinals[0]) - int32_t(previousOrdinals[0]))));
            maxSkip = max(maxSkip, uint32_t(abs(
                int32_t(ordinals[1]) - int32_t(previousOrdinals[1]))));
            maxDrift = max(maxDrift, uint32_t(abs(
                (int32_t(ordinals[0]) - int32_t(ordinals[1])) -
                (int32_t(previousOrdinals[0]) - int32_t(previousOrdinals[1]))
                )));
        }
    }
    averageOrdinalOffset = int32_t(std::round(sum / double(markerCount)));
}



void AlignmentInfo::create(
    const Alignment& alignment,
    uint32_t markerCount0,
    uint32_t markerCount1)
{
    create(alignment, array<uint32_t, 2>({markerCount0, markerCount1}));
}



// Given an AlignmentData, return its AlignmentInfo,
// after swapping and/or reverse complementing it
// to make sure it refers to the given OrientedReadId's,
// in that order.
AlignmentInfo AlignmentData::orient(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1) const
{
    OrientedReadId alignmentOrientedReadId0(readIds[0], 0);
    OrientedReadId alignmentOrientedReadId1(readIds[1], isSameStrand ? 0 : 1);
    AlignmentInfo alignmentInfo = info;

    // Do a swap, if needed.
    if(alignmentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
        alignmentInfo.swap();
        swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
    }
    SHASTA_ASSERT(alignmentOrientedReadId0.getReadId() == orientedReadId0.getReadId());

    // Reverse complement, if needed.
    if(alignmentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
        alignmentOrientedReadId0.flipStrand();
        alignmentOrientedReadId1.flipStrand();
        alignmentInfo.reverseComplement();
    }
    SHASTA_ASSERT(alignmentOrientedReadId0 == orientedReadId0);
    SHASTA_ASSERT(alignmentOrientedReadId1 == orientedReadId1);

    return alignmentInfo;

}



// Return the number of bases in the range covered by the alignment.
uint32_t AlignmentInfo::Data::baseRange(
    uint64_t k,
    OrientedReadId orientedReadId,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) const
{
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];
    const uint32_t position0 = orientedReadMarkers[firstOrdinal].position;
    const uint32_t position1 = orientedReadMarkers[lastOrdinal].position + uint32_t(k);
    return position1 - position0;
}


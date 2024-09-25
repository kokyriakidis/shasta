#include "mode3-Anchor.hpp"
#include "MarkerGraph.hpp"
using namespace shasta;
using namespace mode3;


// This constructor creates the Anchor MarkerIntervals from marker graph edges.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const MarkerGraph& markerGraph) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    anchorMarkerIntervals.createNew(largeDataName("AnchorMarkerIntervals"), largeDataPageSize);

    // For now copy the marker intervals from the marker graph graph.
    for(uint64_t anchorId=0; anchorId<markerGraph.edgeMarkerIntervals.size(); anchorId++) {
        const auto v = markerGraph.edgeMarkerIntervals[anchorId];
        anchorMarkerIntervals.appendVector(v.begin(), v.end());
    }
}



// This constructor access existing Anchors.
Anchors::Anchors(const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    anchorMarkerIntervals.accessExistingReadOnly(largeDataName("AnchorMarkerIntervals"));
}



Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerIntervals[anchorId];
}



uint64_t Anchors::size() const
{
    return anchorMarkerIntervals.size();
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

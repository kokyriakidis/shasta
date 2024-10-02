/********************************************************************************

Creation of Anchors from the MarkerGraph.

Each Anchor corresponds to a "primary marker graph edge"
defined as follows:

- All contributing oriented reads have exactly the same sequence.
  If more than one distinct sequence is present, the edge (anchor)
  is split into two edges (anchors).
- Edge coverage is >= minPrimaryCoverage and <= maxPrimaryCoverage.
- Both vertices have no duplicate ReadIds, and as a result the
  resulting anchor has no duplicate ReadIds.

This uses as input the MarkerGraph vertices. MarkerGraph edges
are created implicitly and never stored. Instead, the corresponding
information is stored directly in the Anchors.

********************************************************************************/

#include "mode3-Anchor.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;


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




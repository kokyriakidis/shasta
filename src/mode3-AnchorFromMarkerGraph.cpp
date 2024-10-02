// Creation of Anchors from the MarkerGraph.

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




// Shasta.
#include "mode3b-PathFinder.hpp"
#include "invalid.hpp"
#include "MarkerGraph.hpp"
using namespace shasta;
using namespace mode3b;

// Standard library.
#include "iostream.hpp"
#include <set>
#include "stdexcept.hpp"



PathFinder::PathFinder(
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction              // 0=forward, 1=backward
    ) :
    markers(markers),
    markerGraph(markerGraph)
{
    // Move in the specified direction until we find an edge with similar
    // read composition which is suitable to becoe the next primary vertex.
    const MarkerGraphEdgeId nextEdgeId = findNextPrimaryEdge(startEdgeId, direction);
    cout << "Next primary edge is " << nextEdgeId << endl;

    throw runtime_error("Not implemented.");
}



// Move in the specified direction until we find an edge with similar
// read composition which is suitable to becoe the next primary vertex.
MarkerGraphEdgeId PathFinder::findNextPrimaryEdge(
    MarkerGraphEdgeId edgeId0,
    uint64_t direction)
{
    cout << "PathFinder::findNextPrimaryEdge begins for edge " << edgeId0 << endl;
    const auto markerIntervals0 = markerGraph.edgeMarkerIntervals[edgeId0];

    // The edges we already checked.
    std::set<MarkerGraphEdgeId> edgeIds;

    // Try increasing marker offsets in the given direction.
    for(uint32_t ordinalOffset=1; /* Check for termination later */ ; ++ordinalOffset) {

        // Try this offset on all of our marker intervals.
        for(const MarkerInterval& markerInterval0: markerIntervals0) {

            // Find the edge that contains the offset marker interval.
            const MarkerGraphEdgeId edgeId1 = markerGraph.locateMarkerIntervalWithOffset(
                markers, markerInterval0, ordinalOffset, direction);

            // If this goes outside the ordinal range for this oriented read, skip it.
            if(edgeId1 == invalid<MarkerGraphEdgeId>) {
                continue;
            }

            // If we already checked it, skip it.
            if(edgeIds.contains(edgeId1)) {
                continue;
            }
            edgeIds.insert(edgeId1);

            cout << "Checking edge " << edgeId1 << endl;

        }

        // Temporary.
        if(ordinalOffset > 3) {
            SHASTA_ASSERT(0);
        }
    }


    throw runtime_error("Not implemented.");
}

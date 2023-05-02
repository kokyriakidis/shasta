// Shasta.
#include "mode3b-PathFinder.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
using namespace shasta;
using namespace mode3b;

// Standard library.
#include "iostream.hpp"
#include <set>
#include "stdexcept.hpp"



PathFinder::PathFinder(
    uint64_t k,
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
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0)) {
        return invalid<MarkerGraphEdgeId>;
    }

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

            // If this edge is visited by an oriented read more than once, skip it.
            if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeId1)) {
                continue;
            }

            edgeIds.insert(edgeId1);

            const bool isGood = haveSimilarCompositions(edgeId0, edgeId1);
        }

        // Temporary.
        if(ordinalOffset > 200) {
            SHASTA_ASSERT(0);
        }
    }


    throw runtime_error("Not implemented.");
}



// Given two marker graph edges, figure out if they have similar read compositions.
// This assumes that markerGraph.edgeHasDuplicateOrientedReadIds is false
// for both edges.
bool PathFinder::haveSimilarCompositions(
    MarkerGraphEdgeId edgeId0,
    MarkerGraphEdgeId edgeId1) const
{
    cout << "Checking read compositions of edges " << edgeId0 << " and " << edgeId1 << endl;

    const auto markerIntervals0 = markerGraph.edgeMarkerIntervals[edgeId0];
    const auto markerIntervals1 = markerGraph.edgeMarkerIntervals[edgeId1];
    auto begin0 = markerIntervals0.begin();
    auto begin1 = markerIntervals1.begin();
    auto end0 = markerIntervals0.end();
    auto end1 = markerIntervals1.end();

    cout << "Number of marker intervals " << markerIntervals0.size() << " " << markerIntervals1.size() << endl;

    // Joint loop over the marker intervals of the two reads.
    // This assumes that there are no duplicate OrientedReadIds.
    uint64_t commonCount = 0;
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 and it1!=end1) {
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
            continue;
        }
        if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
            continue;
        }

        // If getting here, we found two marker intervals for the same OrientedReadId.
        ++commonCount;
        const OrientedReadId orientedReadId = it0->orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Compute the gap in bases between the two MarkerIntervals.
        const uint64_t position0 = orientedReadMarkers[it0->ordinals[1]].position + k;
        const uint64_t position1 = orientedReadMarkers[it1->ordinals[0]].position + k;
        const uint64_t gap = position1 - position0;

        cout << orientedReadId << " gap " << gap << endl;

        // Update for the next iteration of the joint loop.
        ++it0;
        ++it1;
    }
    cout << "Common count " << commonCount << endl;

    return false;
}


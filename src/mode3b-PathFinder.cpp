// Shasta.
#include "mode3b-PathFinder.hpp"
#include "Assembler.hpp"
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
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction              // 0=forward, 1=backward
    ) :
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxMarkerOffset = 10000;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.75;

    // Move in the specified direction until we find an edge with similar
    // read composition which is suitable to becoe the next primary vertex.
    MarkerGraphEdgeId edgeId = startEdgeId;
    while(true) {
        const MarkerGraphEdgeId nextEdgeId = findNextPrimaryEdge(
            edgeId, direction, maxMarkerOffset, minCommonCount, minCorrectedJaccard);
        if(nextEdgeId == invalid<MarkerGraphEdgeId>) {
            break;
        }
        // cout << "Next primary edge is " << nextEdgeId << endl;
        edgeId = nextEdgeId;
    }
}



// Move in the specified direction until we find an edge with similar
// read composition which is suitable to becoe the next primary vertex.
MarkerGraphEdgeId PathFinder::findNextPrimaryEdge(
    MarkerGraphEdgeId edgeId0,
    uint64_t direction,
    uint64_t maxMarkerOffset,
    uint64_t minCommonCount,
    double minCorrectedJaccard) const
{

    if(assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0)) {
        return invalid<MarkerGraphEdgeId>;
    }

    const auto markerIntervals0 = assembler.markerGraph.edgeMarkerIntervals[edgeId0];
    /*
    cout << "PathFinder::findNextPrimaryEdge begins for edge " << edgeId0 <<
        " with coverage " << markerIntervals0.size() << endl;
    */

    // The edges we already checked.
    std::set<MarkerGraphEdgeId> edgeIds;

    // Try increasing marker offsets in the given direction.
    for(uint32_t ordinalOffset=1; ordinalOffset<=maxMarkerOffset ; ++ordinalOffset) {

        // Try this offset on all of our marker intervals.
        for(const MarkerInterval& markerInterval0: markerIntervals0) {

            // Find the edge that contains the offset marker interval.
            const MarkerGraphEdgeId edgeId1 = assembler.markerGraph.locateMarkerIntervalWithOffset(
                assembler.markers, markerInterval0, ordinalOffset, direction);

            // If this goes outside the ordinal range for this oriented read, skip it.
            if(edgeId1 == invalid<MarkerGraphEdgeId>) {
                continue;
            }

            // If we already checked it, skip it.
            if(edgeIds.contains(edgeId1)) {
                continue;
            }

            // If this edge is visited by an oriented read more than once, skip it.
            if(assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId1)) {
                continue;
            }

            edgeIds.insert(edgeId1);

            Assembler::MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            const double correctedJaccard = info.correctedJaccard();
            if(info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
                cout << edgeId0 << " " << edgeId1 << ": base offset " <<
                    info.offsetInBases << ", common " << info.common <<
                    ", corrected jaccard " <<
                    " " << info.correctedJaccard() << endl;
                return edgeId1;
            }
        }

    }


    return invalid<MarkerGraphEdgeId>;
}


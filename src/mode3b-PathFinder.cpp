// Shasta.
#include "mode3b-PathFinder.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3b;

// Standard library.
#include "fstream.hpp"
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



void PathFinder::findNextPrimaryEdges(
    MarkerGraphEdgeId edgeId0,
    uint64_t direction,
    uint64_t minCoverage,
    uint64_t maxCoverage,
    uint64_t maxEdgeCount,
    uint64_t maxMarkerOffset,
    uint64_t minCommonCount,
    double minCorrectedJaccard,
    vector<MarkerGraphEdgeId>& nextPrimaryEdges) const
{
    nextPrimaryEdges.clear();

    if(assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0)) {
        return;
    }

    const auto markerIntervals0 = assembler.markerGraph.edgeMarkerIntervals[edgeId0];

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

            // If coverage is not in the required range, skip it.
            const uint64_t coverage = assembler.markerGraph.edgeCoverage(edgeId1);
            if(coverage<minCoverage or coverage>maxCoverage) {
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
                nextPrimaryEdges.push_back(edgeId1);
                if(nextPrimaryEdges.size() >= maxEdgeCount) {
                    return;
                }
            }
        }

    }


}



// This version finds edge pairs and use them to create a graph.
PathFinder::PathFinder(
    const Assembler& assembler) :
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 6;
    const uint64_t maxCoverage = 18;
    const uint64_t maxEdgeCount = 6;
    const uint64_t maxMarkerOffset = 10000;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.8;

    // Loop over all marker graph edges.
    // This is slow and should run in multithreaded code.
    uint64_t primaryEdgeCount = 0;
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgeId> > edgePairs;
    vector<MarkerGraphEdgeId> nextPrimaryEdges;
    for(MarkerGraphEdgeId edgeId0=0; edgeId0<assembler.markerGraph.edges.size(); edgeId0++) {
        if((edgeId0 % 100000) == 0) {
            cout << timestamp << edgeId0 << "/" << assembler.markerGraph.edges.size() << endl;
        }

        // If coverage is not in the required range, skip it.
        const uint64_t coverage = assembler.markerGraph.edgeCoverage(edgeId0);
        if(coverage<minCoverage or coverage>maxCoverage) {
            continue;
        }

        // If the edge has duplicate oriented reads, skip it.
        if(assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0)) {
            continue;
        }
        ++primaryEdgeCount;

        // Find nearby edges with similar read composition, in both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            findNextPrimaryEdges(
                edgeId0,
                direction,
                minCoverage,
                maxCoverage,
                maxEdgeCount,
                maxMarkerOffset,
                minCommonCount,
                minCorrectedJaccard,
                nextPrimaryEdges);
            for(const MarkerGraphEdgeId edgeId1: nextPrimaryEdges) {
                if(direction == 0) {
                    edgePairs.push_back({edgeId0, edgeId1});
                } else {
                    edgePairs.push_back({edgeId1, edgeId0});
                }
            }
        }
    }
    cout << "Total number of marker graph edges is " << assembler.markerGraph.edges.size() << endl;
    cout << "Total number of primary edges is " << primaryEdgeCount << endl;
    cout << "Number of primary edge pairs before deduplication is " << edgePairs.size() << endl;
    deduplicate(edgePairs);
    cout << "Number of primary edge pairs after deduplication is " << edgePairs.size() << endl;


    // Write out the pairs in Graphviz format.
    {
        ofstream out("PathGraph.dot");
        out << "digraph PathGraph {\n";
        for(const auto& p: edgePairs) {
            out << p.first << "->" << p.second << ";\n";
        }
        out << "}\n";
    }


}

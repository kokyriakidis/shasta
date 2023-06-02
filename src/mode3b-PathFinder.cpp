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



#if 0
// Old version without backtracking.
PathFinder::PathFinder(
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction,            // 0=forward, 1=backward
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >& primaryEdges
    ) :
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxMarkerOffset = 10000;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.75;

    // To create the primary edges, move in the specified direction
    // until we find an edge with similar read composition
    // that is suitable to become the next primary vertex.
    primaryEdges.clear();
    MarkerGraphEdgeId edgeId = startEdgeId;
    while(true) {
        const pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> p = findNextPrimaryEdge(
            edgeId, direction, maxMarkerOffset, minCommonCount, minCorrectedJaccard);
        const MarkerGraphEdgeId nextEdgeId = p.first;
        if(nextEdgeId == invalid<MarkerGraphEdgeId>) {
            break;
        }
        edgeId = nextEdgeId;
        primaryEdges.push_back(p);
    }
    cout << "Found " << primaryEdges.size() + 1 << " primary edges including the starting edge." << endl;
}
#endif



// New version with backtracking.
PathFinder::PathFinder(
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction,            // 0=forward, 1=backward
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >& primaryEdges
    ) :
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxMarkerOffset = 10000;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.75;
    const uint64_t maxBacktrackStreakLength = 6;

    std::set<MarkerGraphEdgeId> forbiddenEdgeIds;

    const bool debug = true;



    // To create the primary edges, move in the specified direction
    // until we find an edge with similar read composition
    // that is suitable to become the next primary vertex.
    // If we get stuck, backtrack.
    uint64_t backTrackingStreakLength = 0;
    primaryEdges.clear();
    while(true) {
        if(backTrackingStreakLength > maxBacktrackStreakLength) {
            break;
        }

        // The last edge we have.
        const MarkerGraphEdgeId edgeId =
            (primaryEdges.empty() ? startEdgeId : primaryEdges.back().first);
        if(debug) {
            cout << "Looking for next edge after " << edgeId << endl;
        }

        // Look for the next edge.
        const pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> p = findNextPrimaryEdge(
            edgeId, direction, maxMarkerOffset, minCommonCount, minCorrectedJaccard,
            forbiddenEdgeIds);
        const MarkerGraphEdgeId nextEdgeId = p.first;

        if(nextEdgeId == invalid<MarkerGraphEdgeId>) {

            // We did not find a next edge. Backtrack if possible.
            if(primaryEdges.empty()) {
                break;
            }
            primaryEdges.pop_back();
            forbiddenEdgeIds.insert(edgeId);
            ++backTrackingStreakLength;
            if(debug) {
                cout << "Did not find a next edge. Backtracking. Current length is " <<
                    primaryEdges.size() << endl;
            }

        } else {

            // We found a next edge. Store it.
            backTrackingStreakLength = 0;
            primaryEdges.push_back(p);
            if(debug) {
                cout << "Added " << nextEdgeId <<
                    ", offset " << p.second.offsetInBases <<
                    ", current length is " << primaryEdges.size()  << endl;
            }
        }
    }
    cout << "Found " << primaryEdges.size() + 1 << " primary edges including the starting edge." << endl;
}



// Move in the specified direction until we find an edge with similar
// read composition which is suitable to becoe the next primary vertex.
pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> PathFinder::findNextPrimaryEdge(
    MarkerGraphEdgeId edgeId0,
    uint64_t direction,
    uint64_t maxMarkerOffset,
    uint64_t minCommonCount,
    double minCorrectedJaccard) const
{
    const MarkerGraph& markerGraph = assembler.markerGraph;
    const auto& markers = assembler.markers;

    // Check for duplicate oriented reads on edgeId0 or its vertices.
    const MarkerGraph::Edge& edge0 = markerGraph.edges[edgeId0];
    if(
        markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0) or
        markerGraph.vertexHasDuplicateOrientedReadIds(edge0.source, markers) or
        markerGraph.vertexHasDuplicateOrientedReadIds(edge0.target, markers)) {
        return make_pair(invalid<MarkerGraphEdgeId>, MarkerGraphEdgePairInfo());
    }

    const auto markerIntervals0 = markerGraph.edgeMarkerIntervals[edgeId0];
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

            // Check for duplicate oriented reads on edgeId1 or its vertices.
            const MarkerGraph::Edge& edge1 = markerGraph.edges[edgeId1];
            if(
                markerGraph.edgeHasDuplicateOrientedReadIds(edgeId1) or
                markerGraph.vertexHasDuplicateOrientedReadIds(edge1.source, markers) or
                markerGraph.vertexHasDuplicateOrientedReadIds(edge1.target, markers)) {
                continue;
            }

            edgeIds.insert(edgeId1);

            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            const double correctedJaccard = info.correctedJaccard();
            if(info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
                cout << edgeId0 << " " << edgeId1 << ": base offset " <<
                    info.offsetInBases << ", common " << info.common <<
                    ", corrected jaccard " <<
                    " " << info.correctedJaccard() << endl;
                return make_pair(edgeId1, info);
            }
        }

    }


    return make_pair(invalid<MarkerGraphEdgeId>, MarkerGraphEdgePairInfo());
}



// Same, but with a forbidden list that can be used for backtracking.
pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> PathFinder::findNextPrimaryEdge(
    MarkerGraphEdgeId edgeId0,
    uint64_t direction,
    uint64_t maxMarkerOffset,
    uint64_t minCommonCount,
    double minCorrectedJaccard,
    const std::set<MarkerGraphEdgeId>& forbiddedEdgeIds) const
{
    const MarkerGraph& markerGraph = assembler.markerGraph;
    const auto& markers = assembler.markers;

    // Check for duplicate oriented reads on edgeId0 or its vertices.
    const MarkerGraph::Edge& edge0 = markerGraph.edges[edgeId0];
    if(
        markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0) or
        markerGraph.vertexHasDuplicateOrientedReadIds(edge0.source, markers) or
        markerGraph.vertexHasDuplicateOrientedReadIds(edge0.target, markers)) {
        return make_pair(invalid<MarkerGraphEdgeId>, MarkerGraphEdgePairInfo());
    }

    const auto markerIntervals0 = markerGraph.edgeMarkerIntervals[edgeId0];
    /*
    cout << "PathFinder::findNextPrimaryEdge begins for edge " << edgeId0 <<
        " with coverage " << markerIntervals0.size() << endl;
    */

    // The edges we already checked.
    std::set<MarkerGraphEdgeId> alreadyEncounteredEdgeIds;

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
            if(alreadyEncounteredEdgeIds.contains(edgeId1)) {
                continue;
            }

            // If it is in the forbidden list, skip it.
            if(forbiddedEdgeIds.contains(edgeId1)) {
                continue;
            }

            // Check for duplicate oriented reads on edgeId1 or its vertices.
            const MarkerGraph::Edge& edge1 = markerGraph.edges[edgeId1];
            if(
                markerGraph.edgeHasDuplicateOrientedReadIds(edgeId1) or
                markerGraph.vertexHasDuplicateOrientedReadIds(edge1.source, markers) or
                markerGraph.vertexHasDuplicateOrientedReadIds(edge1.target, markers)) {
                continue;
            }

            alreadyEncounteredEdgeIds.insert(edgeId1);

            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            const double correctedJaccard = info.correctedJaccard();
            if(info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
                cout << edgeId0 << " " << edgeId1 << ": base offset " <<
                    info.offsetInBases << ", common " << info.common <<
                    ", corrected jaccard " <<
                    " " << info.correctedJaccard() << endl;
                return make_pair(edgeId1, info);
            }
        }

    }


    return make_pair(invalid<MarkerGraphEdgeId>, MarkerGraphEdgePairInfo());

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
    vector< pair<MarkerGraphEdgeId, int64_t> >& nextPrimaryEdges // With offset in bases
    ) const
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

            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            const double correctedJaccard = info.correctedJaccard();
            if(info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
                nextPrimaryEdges.push_back({edgeId1, info.offsetInBases});
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
    const uint64_t minCoverage = 8;
    const uint64_t maxCoverage = 15;
    const uint64_t maxEdgeCount = 8;
    const uint64_t maxMarkerOffset = 30000;
    const uint64_t minCommonCount = 5;
    const double minCorrectedJaccard = 0.9;

    class EdgePair {
    public:
        MarkerGraphEdgeId edgeId0;
        MarkerGraphEdgeId edgeId1;
        uint64_t offsetInBases;
        bool operator==(const EdgePair& that) const
        {
            return tie(edgeId0, edgeId1) == tie(that.edgeId0, that.edgeId1);
        }
        bool operator<(const EdgePair& that) const
        {
            return tie(edgeId0, edgeId1) < tie(that.edgeId0, that.edgeId1);
        }
    };
    vector<EdgePair> edgePairs;


    // Loop over all marker graph edges.
    // This is slow and should run in multithreaded code.
    uint64_t primaryEdgeCount = 0;
    vector< pair<MarkerGraphEdgeId, int64_t> > nextPrimaryEdges; // With offset in bases
    for(MarkerGraphEdgeId edgeId0=0; edgeId0<assembler.markerGraph.edges.size(); edgeId0++) {
        if((edgeId0 % 100000) == 0) {
            cout << timestamp << edgeId0 << "/" << assembler.markerGraph.edges.size() << endl;
        }

        // If coverage is not in the required range, skip it.
        const uint64_t coverage = assembler.markerGraph.edgeCoverage(edgeId0);
        if(coverage<minCoverage or coverage>maxCoverage) {
            continue;
        }

        // Check for duplicate oriented reads on edgeId0 or its vertices.
        const MarkerGraph::Edge& edge0 = assembler.markerGraph.edges[edgeId0];
        if(
            assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId0) or
            assembler.markerGraph.vertexHasDuplicateOrientedReadIds(edge0.source, assembler.markers) or
            assembler.markerGraph.vertexHasDuplicateOrientedReadIds(edge0.target, assembler.markers)) {
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
            for(const auto& p: nextPrimaryEdges) {
                const MarkerGraphEdgeId edgeId1 = p.first;
                const uint64_t offsetInBases = p.second;
                if(direction == 0) {
                    edgePairs.push_back({edgeId0, edgeId1, offsetInBases});
                } else {
                    edgePairs.push_back({edgeId1, edgeId0, -offsetInBases});
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
        for(const EdgePair& edgePair: edgePairs) {
            out << edgePair.edgeId0 << "->" << edgePair.edgeId1 << ";\n";
        }
        out << "}\n";
    }


}

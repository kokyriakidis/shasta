// Shasta.
#include "mode3b-PathFinder.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3b;

// Standard library.
#include <cstdlib>
#include "fstream.hpp"
#include "iostream.hpp"
#include <set>
#include "stdexcept.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<PathFinder>;



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



// Similar to the above, but with backtracking.
PathFinder::PathFinder(
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction,            // 0=forward, 1=backward
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >& primaryEdges
    ) :
    MultithreadedObject<PathFinder>(*this),
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxMarkerOffset = 30000;
    const uint64_t minCoverage = 8;
    const uint64_t maxCoverage = 35;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.8;
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
            edgeId, direction,
            minCoverage, maxCoverage,
            maxMarkerOffset, minCommonCount, minCorrectedJaccard,
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



#if 0
// This version finds a few possible choices for the next primary edge at each step,
// then chooses the best.
PathFinder::PathFinder(
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction,            // 0=forward, 1=backward
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >& primaryEdges
    ) :
    assembler(assembler)
{

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 8;
    const uint64_t maxCoverage = 30;
    const uint64_t maxEdgeCount = 6;
    const uint64_t maxMarkerOffset = 30000;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.6;


    const bool debug = true;



    primaryEdges.clear();
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > nextPrimaryEdges;
    while(true) {

        // The last edge we have.
        const MarkerGraphEdgeId edgeId =
            (primaryEdges.empty() ? startEdgeId : primaryEdges.back().first);
        if(debug) {
            cout << "Looking for next edge after " << edgeId << endl;
        }

        // Find possible choices for the next primary edge.
        findNextPrimaryEdges(
            edgeId, direction,
            minCoverage,
            maxCoverage,
            maxEdgeCount,
            maxMarkerOffset,
            minCommonCount,
            minCorrectedJaccard,
            nextPrimaryEdges);

        // If none found, we are done. The path ends here.
        if(nextPrimaryEdges.empty()) {
            break;
        }

        // Sort them by decreasing number of common oriented reads.
        sort(nextPrimaryEdges.begin(), nextPrimaryEdges.end(),
            OrderPairsBySecondOnlyGreater<MarkerGraphEdgeId, MarkerGraphEdgePairInfo>());

        if(debug) {
            cout << "Found " << nextPrimaryEdges.size() << " possible choices:" << endl;
            for(const auto& p: nextPrimaryEdges) {
                cout << p.first << " " <<
                    p.second.common << " " <<
                    p.second.offsetInBases << " " <<
                    p.second.correctedJaccard() << endl;
            }
        }

        // Add the best choice to our primary edges.
        primaryEdges.push_back(nextPrimaryEdges.front());
        if(debug) {
            cout << "Added " << primaryEdges.back().first <<
                ", offset " << primaryEdges.back().second.offsetInBases <<
                ", current length is " << primaryEdges.size()  << endl;
        }
    }
    cout << "Found " << primaryEdges.size() + 1 << " primary edges including the starting edge." << endl;
}
#endif



// Move in the specified direction until we find an edge with similar
// read composition which is suitable to become the next primary vertex.
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
    uint64_t minCoverage,
    uint64_t maxCoverage,
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

            const uint64_t coverage = assembler.markerGraph.edgeCoverage(edgeId1);
            if(coverage < minCoverage or coverage > maxCoverage) {
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



// This returns a few possible choices for the next primary edge.
// It starts on edgeId0 and moves in the specified direction
// (0 = forward, 1 = backward).
void PathFinder::findNextPrimaryEdges(
    MarkerGraphEdgeId edgeId0,
    uint64_t direction,
    uint64_t minCoverage,
    uint64_t maxCoverage,
    uint64_t maxEdgeCount,
    uint64_t maxMarkerOffset,
    uint64_t minCommonCount,
    double minCorrectedJaccard,
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >& nextPrimaryEdges
    ) const
{
    const bool debug = false;
    if(debug) {
        cout << "Looking for next primary edge after " << edgeId0 <<
            ", coverage " << assembler.markerGraph.edgeCoverage(edgeId0) <<
            ", direction " << direction << endl;
    }

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
                if(false) {
                    cout << edgeId1 << " discarded due to coverage " << coverage << endl;
                }
                continue;
            }

            // If we already checked it, skip it.
            if(edgeIds.contains(edgeId1)) {
                continue;
            }

            // If this edge or one of its vertices is visited by an oriented read more than once, skip it.
            const MarkerGraph::Edge& edge1 = assembler.markerGraph.edges[edgeId1];
            if(
                assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId1) or
                assembler.markerGraph.vertexHasDuplicateOrientedReadIds(edge1.source, assembler.markers) or
                assembler.markerGraph.vertexHasDuplicateOrientedReadIds(edge1.target, assembler.markers)) {
                continue;
            }

            edgeIds.insert(edgeId1);

            // Analyze its read composition.
            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            const double correctedJaccard = info.correctedJaccard();
            if(info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
                // Store it.
                nextPrimaryEdges.push_back({edgeId1, info});
                if(nextPrimaryEdges.size() >= maxEdgeCount) {
                    return;
                }
            }
        }

    }


}



// This version finds edge pairs and use them to create a graph.
PathFinder::PathFinder(
    const Assembler& assembler,
    uint64_t threadCount) :
    MultithreadedObject<PathFinder>(*this),
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxMarkerOffset = 10000; // 30000, reduced for speed
    // const uint64_t minCoverage = 15; // Change back to 8 when done debugging.
    // const uint64_t maxCoverage = 20; // Change back to 35 when done debugging.
    const uint64_t minCoverage = 8;
    const uint64_t maxCoverage = 35;
    const uint64_t minCommonCount = 3;      // 6, reduced for speed
    const double minCorrectedJaccard = 0.8;
    const uint64_t maxEdgeCount = 3;    // 6, reduced for speed

    // Store the parameters so the threads can see them.
    threadFunction1Data.maxMarkerOffset = maxMarkerOffset;
    threadFunction1Data.minCoverage = minCoverage;
    threadFunction1Data.maxCoverage = maxCoverage;
    threadFunction1Data.minCommonCount =  minCommonCount;
    threadFunction1Data.minCorrectedJaccard = minCorrectedJaccard;
    threadFunction1Data.maxEdgeCount = maxEdgeCount;

    // Make space for the EdgePairs found by each thread.
    threadFunction1Data.threadEdgePairs.clear();
    threadFunction1Data.threadEdgePairs.resize(threadCount);

    // Find the EdgePairs.
    const uint64_t batchSize = 100;
    setupLoadBalancing(assembler.markerGraph.edges.size(), batchSize);
    runThreads(&PathFinder::threadFunction1, threadCount);

    // Merge the EdgePairs found by all threads.
    edgePairs.clear();
    for(const vector<EdgePair>& thisThreadEdgePairs: threadFunction1Data.threadEdgePairs) {
        copy(thisThreadEdgePairs.begin(), thisThreadEdgePairs.end(), back_inserter(edgePairs));
    }
    threadFunction1Data.threadEdgePairs.clear();

    cout << "Total number of marker graph edges is " << assembler.markerGraph.edges.size() << endl;
    cout << "Number of primary edge pairs before deduplication is " << edgePairs.size() << endl;
    deduplicate(edgePairs);
    cout << "Number of primary edge pairs after deduplication is " << edgePairs.size() << endl;


    // Write out the pairs in Graphviz format.
    // The len attribute is only honored by neato and fdp.
    {
        ofstream out("PathGraph.dot");
        out << "digraph PathGraph {\n";
        for(const EdgePair& edgePair: edgePairs) {
            out << edgePair.edgeId0 << "->" << edgePair.edgeId1 <<
                " [len=" << max(1UL, uint64_t(0.001 * double(edgePair.offsetInBases))) << "]"
                ";\n";
        }
        out << "}\n";
    }
}



void PathFinder::threadFunction1(uint64_t threadId)
{
    // Get the parameters.
    const uint64_t maxMarkerOffset = threadFunction1Data.maxMarkerOffset;
    const uint64_t minCoverage = threadFunction1Data.minCoverage;
    const uint64_t maxCoverage = threadFunction1Data.maxCoverage;
    const uint64_t minCommonCount = threadFunction1Data.minCommonCount;
    const double minCorrectedJaccard = threadFunction1Data.minCorrectedJaccard;
    const uint64_t maxEdgeCount = threadFunction1Data.maxEdgeCount;

    // The vector where we will store the EdgePairs found by this thread.
    vector<EdgePair>& thisThreadEdgePairs = threadFunction1Data.threadEdgePairs[threadId];

    // Work vector used in the main loop but defined here to reduce memory allocation activity.
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > nextPrimaryEdges;

    // Look over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all marker graph edges assigned to this batch.
        for(MarkerGraphEdgeId edgeId0=begin; edgeId0<end; edgeId0++) {
            if((edgeId0 % 100000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
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
                    const uint64_t offsetInBases = p.second.offsetInBases;
                    if(direction == 0) {
                        thisThreadEdgePairs.push_back({edgeId0, edgeId1, offsetInBases});
                    } else {
                        thisThreadEdgePairs.push_back({edgeId1, edgeId0, -offsetInBases});
                    }
                }
            }
        }
    }
}


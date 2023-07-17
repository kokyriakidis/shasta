// Shasta.
#include "mode3b-PathFinder.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "longestPath.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

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
    MappedMemoryOwner(assembler),
    MultithreadedObject<PathFinder>(*this),
    assembler(assembler)
{
    SHASTA_ASSERT(direction < 2);

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
            cout << "Looking for next edge after " << edgeId << ", direction " << direction << endl;
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
            const bool hasConsistentOffset =
                (direction==0 and info.offsetInBases >= 0) or
                (direction==1 and info.offsetInBases <= 0);
            if(hasConsistentOffset and info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
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
            const bool hasConsistentOffset =
                (direction==0 and info.offsetInBases >= 0) or
                (direction==1 and info.offsetInBases <= 0);
            if(hasConsistentOffset and info.common >= minCommonCount and correctedJaccard >= minCorrectedJaccard) {
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

            // If it satisfies our criteria, store it.
            const bool hasConsistentOffset =
                (direction==0 and info.offsetInBases >= 0) or
                (direction==1 and info.offsetInBases <= 0);
            if(
                hasConsistentOffset and
                info.common >= minCommonCount and
                correctedJaccard >= minCorrectedJaccard) {
                nextPrimaryEdges.push_back({edgeId1, info});
                if(nextPrimaryEdges.size() >= maxEdgeCount) {
                    return;
                }
            }
        }

    }


}



// This returns a few possible choices for the next primary edge.
// It starts on edgeId0 and moves in the specified direction
// (0 = forward, 1 = backward).
// Fast version that used the edge table.
void PathFinder::findNextPrimaryEdgesFast(
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
            const OrientedReadId orientedReadId = markerInterval0.orientedReadId;

            // Compute the offset ordinal.
            // If we end up outside the read, do nothing.
            uint32_t ordinal0WithOffset;
            if(direction == 0) {
                ordinal0WithOffset = markerInterval0.ordinals[0] + ordinalOffset;
                if(ordinal0WithOffset >= assembler.markers.size(orientedReadId.getValue())) {
                    continue;
                }
            } else {
                if(ordinalOffset > markerInterval0.ordinals[0]) {
                    continue;
                }
                ordinal0WithOffset = markerInterval0.ordinals[0] - ordinalOffset;
            }

            // Find the edgeId that begins at this ordinal.
            const MarkerGraphEdgeId edgeId1 = markerGraphEdgeTable[orientedReadId.getValue()][ordinal0WithOffset];
            if(edgeId1 == invalid<MarkerGraphEdgeId>) {
                continue;
            }

            // If we already checked it, skip it.
            if(edgeIds.contains(edgeId1)) {
                continue;
            }

            edgeIds.insert(edgeId1);

            // Analyze its read composition.
            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            const double correctedJaccard = info.correctedJaccard();

            // If it satisfies our criteria, store it.
            const bool hasConsistentOffset =
                (direction==0 and info.offsetInBases >= 0) or
                (direction==1 and info.offsetInBases <= 0);
            if(
                hasConsistentOffset and
                info.common >= minCommonCount and
                correctedJaccard >= minCorrectedJaccard) {
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
    MappedMemoryOwner(assembler),
    MultithreadedObject<PathFinder>(*this),
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxMarkerOffset = 1000;
    const uint64_t minCoverage = 15;
    const uint64_t maxCoverage = 35;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.8;
    const uint64_t maxEdgeCount = 2;

    // Find edge pairs with similar read compositions.
    findEdgePairs(
        threadCount,
        maxMarkerOffset,
        minCoverage,
        maxCoverage,
        minCommonCount,
        minCorrectedJaccard,
        maxEdgeCount);
    writeEdgePairsGraphviz();

    // Find connected components of marker graph edges.
    findComponents();


    // Create an AssemblyPath for the largest few components.
    const uint64_t componentCount = 10;
    for(uint64_t componentRank=0;
        componentRank < min(componentCount, componentIndex.size());
        componentRank++) {

        // Access this component.
        const uint64_t componentId = componentIndex[componentRank].first;
        const Graph& component = components[componentId];

        // Compute a longest path.
        vector<MarkerGraphEdgeId> primaryEdges;
        vector<MarkerGraphEdgePairInfo> infos;
        component.getLongestPath(primaryEdges, infos);

        // Create an assembly path.
        AssemblyPath assemblyPath(assembler, primaryEdges, infos);
        vector<Base> sequence;
        assemblyPath.getSequence(sequence);

        ofstream fasta("AssemblyPath-" + to_string(componentRank) + ".fasta");
        assemblyPath.writeFasta(fasta);
        ofstream csv("AssemblyPath.csv");
        assemblyPath.writeCsv(csv);
        cout << "Component " << componentRank << " has " << num_vertices(component) <<
            " vertices (primary marker graph edges)." << endl;
        cout << "\tIts longest path has " << primaryEdges.size() <<
            " vertices (primary marker graph edges)." << endl;
        cout << "\tLinearity ratio is " <<
            double(primaryEdges.size()) / double(num_vertices(component)) << endl;
        cout << "\tAssembled sequence length " << sequence.size() << endl;
    }

}



// Write out the edge pairs in Graphviz format.
// The len attribute is only honored by neato and fdp.
void PathFinder::writeEdgePairsGraphviz() const
{

    ofstream out("PathGraph.dot");
    out << "digraph PathGraph {\n";
    for(const EdgePair& edgePair: edgePairs) {
        out << edgePair.edgeId0 << "->" << edgePair.edgeId1 <<
            " [len=" << max(1UL, uint64_t(0.001 * double(edgePair.info.offsetInBases))) << "]"
            ";\n";
    }
    out << "}\n";
}



void PathFinder::findEdgePairs(
    uint64_t threadCount,
    uint64_t maxMarkerOffset,
    uint64_t minCoverage,
    uint64_t maxCoverage,
    uint64_t minCommonCount,
    double minCorrectedJaccard,
    uint64_t maxEdgeCount)
{
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
    createMarkerGraphEdgeTable(threadCount, minCoverage, maxCoverage);
    setupLoadBalancing(assembler.markerGraph.edges.size(), batchSize);
    runThreads(&PathFinder::threadFunction1, threadCount);
    markerGraphEdgeTable.remove();

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
            if((edgeId0 % 1000000) == 0) {
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
                findNextPrimaryEdgesFast(
                    edgeId0,
                    direction,
                    minCoverage,
                    maxCoverage,
                    maxEdgeCount,
                    maxMarkerOffset,
                    minCommonCount,
                    minCorrectedJaccard,
                    nextPrimaryEdges);
                for(auto& p: nextPrimaryEdges) {
                    const MarkerGraphEdgeId edgeId1 = p.first;
                    MarkerGraphEdgePairInfo& info = p.second;
                    if(direction == 0) {
                        thisThreadEdgePairs.push_back({edgeId0, edgeId1, info});
                    } else {
                        info.reverse();
                        thisThreadEdgePairs.push_back({edgeId1, edgeId0, info});
                    }
                }
            }
        }
    }
}



PathFinder::PathFinder(
    const Assembler& assembler) :
    MappedMemoryOwner(assembler),
    MultithreadedObject<PathFinder>(*this),
    assembler(assembler)
{
}



void PathFinder::createMarkerGraphEdgeTable(
    uint64_t threadCount,
    uint64_t minCoverage,
    uint64_t maxCoverage
    )
{
    createNew(markerGraphEdgeTable, "MarkerGraphEdgeTable");
    for(uint64_t i=0; i<assembler.markers.size(); i++) {
        markerGraphEdgeTable.appendVector(assembler.markers.size(i));
    }
    fill(markerGraphEdgeTable.begin(), markerGraphEdgeTable.end(), invalid<MarkerGraphEdgeId>);

    createMarkerGraphEdgeTableData.minCoverage = minCoverage;
    createMarkerGraphEdgeTableData.maxCoverage = maxCoverage;
    setupLoadBalancing(assembler.markerGraph.edges.size(), 1000);
    runThreads(&PathFinder::createMarkerGraphEdgeTableThreadFunction, threadCount);
}


void PathFinder::createMarkerGraphEdgeTableThreadFunction(uint64_t threadId)
{
    const uint64_t minCoverage = createMarkerGraphEdgeTableData.minCoverage;
    const uint64_t maxCoverage = createMarkerGraphEdgeTableData.maxCoverage;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerGraphEdgeId edgeId=begin; edgeId!=end; ++edgeId) {

            // Check coverage.
            const uint64_t coverage = assembler.markerGraph.edgeCoverage(edgeId);
            if(coverage<minCoverage or coverage>maxCoverage) {
                continue;
            }

            // Check for duplicate oriented reads.
            const MarkerGraph::Edge& edge = assembler.markerGraph.edges[edgeId];
            if(
                assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeId) or
                assembler.markerGraph.vertexHasDuplicateOrientedReadIds(edge.source, assembler.markers) or
                assembler.markerGraph.vertexHasDuplicateOrientedReadIds(edge.target, assembler.markers)) {
                continue;
            }

            const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                const OrientedReadId orientedReadId = markerInterval.orientedReadId;
                const uint32_t ordinal0 = markerInterval.ordinals[0];
                markerGraphEdgeTable[orientedReadId.getValue()][ordinal0] = edgeId;
            }
        }
    }
}



void PathFinder::findComponents()
{

    // Gather all the MarkerGraphEdgeIds that appear in edge pairs.
    vector<MarkerGraphEdgeId> primaryEdges;
    for(const EdgePair& edgePair: edgePairs) {
        primaryEdges.push_back(edgePair.edgeId0);
        primaryEdges.push_back(edgePair.edgeId1);
    }
    deduplicate(primaryEdges);
    cout << "Found " << primaryEdges.size() <<
        " primary edges that appear in at least one edge pair." << endl;

    // Compute connected components.
    const uint64_t n = primaryEdges.size();
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
    for(const EdgePair& edgePair: edgePairs) {
        const auto it0 = lower_bound(primaryEdges.begin(), primaryEdges.end(), edgePair.edgeId0);
        const auto it1 = lower_bound(primaryEdges.begin(), primaryEdges.end(), edgePair.edgeId1);
        SHASTA_ASSERT(it0 != primaryEdges.end());
        SHASTA_ASSERT(it1 != primaryEdges.end());
        const uint64_t i0 = it0 - primaryEdges.begin();
        const uint64_t i1 = it1 - primaryEdges.begin();
        disjointSets.union_set(i0, i1);
    }

    // Generate vertices of each connected component.
    components.clear();
    components.resize(n);
    for(uint64_t i=0; i<n; i++) {
        const uint64_t componentId = disjointSets.find_set(i);
        components[componentId].addVertex(primaryEdges[i]);
    }

    // Compute a histogram of component sizes.
    vector<uint64_t> histogram;
    for(const Graph& component: components) {
        const uint64_t componentSize = num_vertices(component);
        if(componentSize >= histogram.size()) {
            histogram.resize(componentSize + 1, 0);
        }
        ++histogram[componentSize];
    }
    {
        ofstream csv("ComponentSizeHistogram.csv");
        csv << "Size,Frequency,TotalSize\n";
        for(uint64_t componentSize=1; componentSize<histogram.size(); componentSize++) {
            const uint64_t frequency = histogram[componentSize];
            if(frequency) {
                csv << componentSize << ",";
                csv << frequency << ",";
                csv << frequency * componentSize << "\n";
            }
        }
    }

    // Generate edges.
    for(const EdgePair& edgePair: edgePairs) {
        const auto it0 = lower_bound(primaryEdges.begin(), primaryEdges.end(), edgePair.edgeId0);
        const auto it1 = lower_bound(primaryEdges.begin(), primaryEdges.end(), edgePair.edgeId1);
        SHASTA_ASSERT(it0 != primaryEdges.end());
        SHASTA_ASSERT(it1 != primaryEdges.end());
        const uint64_t i0 = it0 - primaryEdges.begin();
        const uint64_t i1 = it1 - primaryEdges.begin();
        const uint64_t componentId = disjointSets.find_set(i0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(i1));
        components[componentId].addEdge(edgePair.edgeId0, edgePair.edgeId1, edgePair.info);
    }


    // Create the componentIndex.
    componentIndex.clear();
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const auto& component = components[componentId];
        const uint64_t componentSize = num_vertices(component);
        if(componentSize > 0) {
            componentIndex.push_back({componentId, componentSize});
        }
    }
    sort(componentIndex.begin(), componentIndex.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());
}



void PathFinder::Graph::addVertex(MarkerGraphEdgeId edgeId)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    vertexMap.insert({edgeId, add_vertex(Vertex({edgeId}), *this)});
}



void PathFinder::Graph::addEdge(
    MarkerGraphEdgeId edgeId0,
    MarkerGraphEdgeId edgeId1,
    const MarkerGraphEdgePairInfo& info)
{
    auto it0 = vertexMap.find(edgeId0);
    auto it1 = vertexMap.find(edgeId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    add_edge(it0->second, it1->second, Edge({info}), *this);
}



void PathFinder::Graph::getLongestPath(
    vector<MarkerGraphEdgeId>& primaryEdges,
    vector<MarkerGraphEdgePairInfo>& infos) const
{
    const Graph& graph = *this;
    SHASTA_ASSERT(num_vertices(graph));

    // Compute the longest path.
    vector<vertex_descriptor> pathVertices;
    longestPath(graph, pathVertices);

    // Fill in the primary edges.
    primaryEdges.clear();
    for(const vertex_descriptor v: pathVertices) {
        primaryEdges.push_back(graph[v].edgeId);
    }

    // Fill in the infos.
    infos.clear();
    for(uint64_t i=1; i<pathVertices.size(); i++) {
        const vertex_descriptor v0 = pathVertices[i-1];
        const vertex_descriptor v1 = pathVertices[i];
        edge_descriptor e;
        bool edgeExists = false;
        tie(e, edgeExists) = edge(v0, v1, graph);
        SHASTA_ASSERT(edgeExists);
        infos.push_back(graph[e].info);
    }
}


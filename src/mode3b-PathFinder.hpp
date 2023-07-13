#ifndef SHASTA_MODE3B_PATH_FINDER_HPP
#define SHASTA_MODE3B_PATH_FINDER_HPP

#include "MappedMemoryOwner.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"

#include <set>
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3b {
        class PathFinder;
    }

    class Assembler;
};


// Find and assemble a path in the complete marker graph.
class shasta::mode3b::PathFinder :
    public MappedMemoryOwner,
    public MultithreadedObject<PathFinder> {
public:

    PathFinder(
        const Assembler&,
        MarkerGraphEdgeId startEdgeId,  // The path starts here.
        uint64_t direction,             // 0=forward, 1=backward
        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >& primaryEdges
        );

    PathFinder(const Assembler&, uint64_t threadCount);
    PathFinder(const Assembler&);
private:

    // Things we get from the constructor.
    const Assembler& assembler;


    // Move in the specified direction until we find an edge with similar
    // read composition which is suitable to become the next primary vertex.
    pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> findNextPrimaryEdge(
        MarkerGraphEdgeId,
        uint64_t direction,
        uint64_t maxMarkerOffset,
        uint64_t minCommonCount,
        double minCorrectedJaccard) const;

    // Same, but with a forbidden list that can be used for backtracking.
    pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> findNextPrimaryEdge(
        MarkerGraphEdgeId,
        uint64_t direction,
        uint64_t minCoverage,
        uint64_t maxCoverage,
        uint64_t maxMarkerOffset,
        uint64_t minCommonCount,
        double minCorrectedJaccard,
        const std::set<MarkerGraphEdgeId>& forbiddedEdgeIds) const;

public:
    void findNextPrimaryEdges(
        MarkerGraphEdgeId,
        uint64_t direction,
        uint64_t minCoverage,
        uint64_t maxCoverage,
        uint64_t maxEdgeCount,
        uint64_t maxMarkerOffset,
        uint64_t minCommonCount,
        double minCorrectedJaccard,
        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >&
        ) const;
    void findNextPrimaryEdgesFast(
        MarkerGraphEdgeId,
        uint64_t direction,
        uint64_t minCoverage,
        uint64_t maxCoverage,
        uint64_t maxEdgeCount,
        uint64_t maxMarkerOffset,
        uint64_t minCommonCount,
        double minCorrectedJaccard,
        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> >&
        ) const;
private:


    // Multithreaded code used to create a global graph.
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
    void threadFunction1(uint64_t threadId);
    class ThreadFunction1Data {
    public:
        uint64_t maxMarkerOffset;
        uint64_t minCoverage;
        uint64_t maxCoverage;
        uint64_t minCommonCount;
        double minCorrectedJaccard;
        uint64_t maxEdgeCount;
        // The EdgePairs found by each thread.
        vector< vector<EdgePair> > threadEdgePairs;
    };
    ThreadFunction1Data threadFunction1Data;
    // All the EdgePairs found by all threads.
    vector<EdgePair> edgePairs;



    // A table that gives the MarkerGraphEdgeId of the edge
    // whose source vertex contains the marker at a specified
    // oriented read and ordinal:
    // const MarkerGraphEdgeId edgeId = markerGraphEdgeTable[orientedReadId.getValue()][ordinal0];
    // If such a marker graph edge does not exist or is outside the coverage range
    // specified when calling createMarkerGraphEdgeTable,
    // the table stores invalid<MarkerGraphEdgeId>.
    // The table also stores invalid<MarkerGraphEdgeId>
    // If the edge or one of its vertices has duplicate oriented read ids.
    // This is used to speed up findNextPrimaryEdges.
    MemoryMapped::VectorOfVectors<MarkerGraphEdgeId, uint64_t> markerGraphEdgeTable;
    void createMarkerGraphEdgeTable(
        uint64_t threadCount,
        uint64_t minCoverage,
        uint64_t maxCoverage
        );
    void createMarkerGraphEdgeTableThreadFunction(uint64_t threadId);
    class CreateMarkerGraphEdgeTableData {
public:
        uint64_t minCoverage;
        uint64_t maxCoverage;
    };
    CreateMarkerGraphEdgeTableData createMarkerGraphEdgeTableData;

};

#endif

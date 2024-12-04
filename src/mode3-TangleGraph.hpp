#pragma once

// Class TangleGraph is used for read following in the AnchorGraph
// with a given set of entrances and exits.

// Shasta.
#include "mode3-Anchor.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <cstdint.hpp>
#include <vector.hpp>

namespace shasta {
    namespace mode3 {

        class TangleGraphVertex;
        class TangleGraphEdge;
        class TangleGraph;
        using TangleGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            TangleGraphVertex,
            TangleGraphEdge>;
    }
}



class shasta::mode3::TangleGraphVertex {
public:
    AnchorId anchorId;
    uint64_t coverage = 0;
    TangleGraphVertex(AnchorId anchorId) : anchorId(anchorId) {}

    // Each vertex can appear in zero or one Entrances
    // and zero or one Exits.
    uint64_t entranceIndex = invalid<uint64_t>;
    uint64_t exitIndex = invalid<uint64_t>;
};



class shasta::mode3::TangleGraphEdge {
public:
    uint64_t coverage = 0;
};



class shasta::mode3::TangleGraph : public TangleGraphBaseClass {
public:
    TangleGraph(
        bool debug,
        const Anchors&,
        const vector<AnchorId>& entranceAnchors,
        const vector<AnchorId>& exitAnchors,
        bool bidirectional
    );

private:
    bool debug;
    const Anchors& anchors;
    bool bidirectional;

    // A base class to describe an Entrance or Exit to the TangleGraph.
    class EntranceOrExit {
    public:
        AnchorId anchorId;

        // The AnchorMarkerIntervals on that AnchorId.
        // These are initially copied from class Anchors.
        // But later, for entrances we remove AnchorMarkerIntervals
        // for which the same OrientedReadId appears in another entrance;
        // and for exits we remove AnchorMarkerIntervals
        // for which the same OrientedReadId appears in another exit.
        vector<AnchorMarkerInterval> anchorMarkerIntervals;

        // The AnchorIds encountered during read following.
        // For an entrance, read following moves forward, starting at the entrance.
        // For an exit, read following moves backward, starting at the exit.
        // However, if bidirectional is true read following moves in both directions
        // for both entrances and exits.
        vector<AnchorId> journeyAnchorIds;

        // The AnchorIds encountered during read following starting from this Entrance/Exit
        // and no other Entrance/Exit.
        vector<AnchorId> uniqueJourneyAnchorIds;

        EntranceOrExit(AnchorId, const Anchor&);
    };

    class Entrance : public EntranceOrExit {
    public:
        using EntranceOrExit::EntranceOrExit;
        void readFollowing(bool debug, const Anchors&, bool bidirectional);
    };

    class Exit : public EntranceOrExit {
    public:
        using EntranceOrExit::EntranceOrExit;
        void readFollowing(bool debug, const Anchors&, bool bidirectional);
    };

    vector<Entrance> entrances;
    vector<Exit> exits;
    void constructEntrances(const vector<AnchorId>& entranceAnchors);
    void constructExits(const vector<AnchorId>& entranceAnchors);



    // The oriented reads used in this TangleGraph..
    // Each oriented read can appear in at most one entrance and one exit.
    // They are stored sorted by OrientedReadId.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        OrientedReadInfo(OrientedReadId orientedReadId) : orientedReadId(orientedReadId) {}
        bool operator<(const OrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }

        // These are set if this oriented read appears in one entrance.
        uint64_t entranceIndex = invalid<uint64_t>;
        uint32_t entrancePositionInJourney = invalid<uint32_t>;

        // These are set if this oriented read appears in one exit.
        uint64_t exitIndex = invalid<uint64_t>;
        uint32_t exitPositionInJourney = invalid<uint32_t>;

        // The journey of this oriented read in the tangle graph.
        // This is obtained from the journey of the oriented read in the Anchors,
        // possibly clipped at the beginning/end if bidirectional==false.
        // Only AnchorIds that correspond to a vertex are entered.
        vector<vertex_descriptor> tangleJourney;
    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    OrientedReadInfo& getOrientedReadInfo(OrientedReadId);



    // Create TangleGraph vertices.
    // There is a vertex for each AnchorId that is unique to one Entrance and/or one Exit.
    void createVertices();

    // Create edges.
    // This uses the OrientedReadInfo::tangleJourney.
    void createEdges();

    // The tangle matrix is the number of common unique anchors between each
    // entrance/exit pair.
    // Indexed by [entranceIndex][exitIndex]
    vector< vector<uint64_t> > tangleMatrix;

    // The vertexTable contains pairs of AnchorIds with the corresponding
    // file descriptors. Sorted by AnchorId so findVertex can use std::lowerBound.
    vector< pair<AnchorId, vertex_descriptor> > vertexTable;
    vertex_descriptor getVertex(AnchorId) const;

    void writeGraphviz() const;

};

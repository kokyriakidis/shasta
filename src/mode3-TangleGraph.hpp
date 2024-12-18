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
    TangleGraphVertex(AnchorId anchorId) : anchorId(anchorId) {}

    // The OrientedReadIds that visit this vertex.
    // This only includes OrientedReadIds used in the TangleGraph
    // and therefore is generally a subset of the OrientedReadIds
    // for the same AnchorId in the global AnchorGraph.
    // The OrientedReadIds are stored sorted.
    vector<OrientedReadId> orientedReadIds;
    uint64_t coverage() const
    {
        return orientedReadIds.size();
    }

    bool wasSeenByForwardBfs = false;
    bool wasSeenByBackwardBfs = false;
};



class shasta::mode3::TangleGraphEdge {
public:
    // The OrientedReadIds that contribute to this edge.
    // This only includes OrientedReadIds used in the TangleGraph
    // and therefore is generally a subset of the OrientedReadIds
    // for the same edge of the global AnchorGraph.
    // The OrientedReadIds are stored sorted.
    vector<OrientedReadId> orientedReadIds;
    uint64_t coverage() const
    {
        return orientedReadIds.size();
    }
};



class shasta::mode3::TangleGraph : public TangleGraphBaseClass {
public:
    TangleGraph(
        bool debug,
        uint64_t tangleId,
        const Anchors&,
        const vector<AnchorId>& entranceAnchors,
        const vector<AnchorId>& exitAnchors,
        bool bidirectional,
        double maxLoss,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold);

    // Return true if successful, that is, all Entrances are
    // connected to at least one Exit, and all Exits are
    // connected to at least one Entrance, and the failure flag is not set.
    bool failure = false;
    bool isSuccessful() const;

    void getChains(vector< vector<AnchorId> >&) const;

private:
    bool debug;
    uint64_t tangleId;
    const Anchors& anchors;
    bool bidirectional;

    // A base class to describe an Entrance or Exit to the TangleGraph.
    class EntranceOrExit {
    public:
        AnchorId anchorId;
#if 0
        // The AnchorIds encountered during read following.
        // For an entrance, read following moves forward, starting at the entrance.
        // For an exit, read following moves backward, starting at the exit.
        // However, if bidirectional is true read following moves in both directions
        // for both entrances and exits.
        vector<AnchorId> journeyAnchorIds;

        // The AnchorIds encountered during read following starting from this Entrance/Exit
        // and no other Entrance/Exit.
        vector<AnchorId> uniqueJourneyAnchorIds;
#endif
        EntranceOrExit(AnchorId);
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

    // Find out if a given AnchorId is an entrance or exit.
    bool isEntrance(AnchorId) const;
    bool isExit(AnchorId) const;



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

        // The portion of this oriented read journey used in the TangleGraph.
        uint64_t journeyBegin;
        uint64_t journeyEnd;

        // The journey of this oriented read in the TangleGraph.
        // This is the sequence of vertices it encounters.
        // On ther initial TangleGraph, it is also a path.
        vector<vertex_descriptor> tangleJourney;
    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    OrientedReadInfo* getOrientedReadInfo(OrientedReadId);

    bool createVertices(uint64_t anchorCoverageThreshold, uint64_t minVertexCoverage);
    void createEdges();

    // The vertexTable contains pairs of AnchorIds with the corresponding
    // file descriptors. Sorted by AnchorId so findVertex can use std::lowerBound.
    vector< pair<AnchorId, vertex_descriptor> > vertexTable;
    vertex_descriptor getVertex(AnchorId) const;

    double edgeLoss(edge_descriptor) const;
    void removeWeakEdges(double maxLoss);
    void removeCrossEdges(
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold);

    // This removes vertices that are not reachable by at least one entrance
    // and at least one exit. If any entrance or exit would be removed in this way,
    // this returns false ant sets failure to true.
    bool removeUnreachable();

    void writeGraphviz(const string& name) const;
};

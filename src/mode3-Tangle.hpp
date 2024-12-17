#pragma once

// A Tangle is essentially the same as an AssemblyGraph::Superbubble,
// but it includes code for detangling by read following.
// It requires all BubbleChains in the assembly graph to consist of a single Chain.

#include "mode3-AssemblyGraph.hpp"

namespace shasta {
    namespace mode3 {
        class Tangle;
    }
}


class shasta::mode3::Tangle {
public:

    Tangle(
        bool debug,
        uint64_t tangleId,
        AssemblyGraph&,
        uint64_t maxOffset,
        const vector<AssemblyGraph::vertex_descriptor>& tangleVertices);

    void detangle(
        bool debug,
        uint64_t tangleId,
        AssemblyGraph&,
        double maxLoss,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold,
        vector< vector<AnchorId> >& anchorChains);

    bool success = false;

    bool debug;
    AssemblyGraph& assemblyGraph;

    // The AssemblyGraph vertices that define this tangle.
    // They are a connected component of the AssemblyGraph subgraph
    // constructed using only edges with offset up to maxOffset.
    // They are stored sorted so we can quickly check
    // if a vertex belongs to the tangle.
    vector<AssemblyGraph::vertex_descriptor> tangleVertices;
    void writeTangleVertices() const;
    bool isTangleVertex(AssemblyGraph::vertex_descriptor) const;

    // The AssemblyGraph edges in the tangle.
    // These are edges in which both the source and target vertex
    // are in the tangle, and the edge offset is no more than maxOffset.
    // They are stored sorted so we can quickly check
    // if a vertex belongs to the tangle.
    vector<AssemblyGraph::edge_descriptor> tangleEdges;
    void findTangleEdges(uint64_t maxOffset);
    void writeTangleEdges() const;
    bool isTangleEdge(AssemblyGraph::edge_descriptor) const;



    // The entrances are AssemblyGraph edges that are not in the Tangle
    // but whose target vertex is in the Tangle.
    // The exits are AssemblyGraph edges that are not in the Tangle
    // but whose source vertex is in the Tangle.
    // Note that an AssemblyGraph edge can be both an entrance and an exit at the
    // same time. This can happen if the corresponding Chain has an offset
    // larger than maxOffset, because in that case the edge is not part of the Tangle.
    // The AnchorId of an entrance is the second to last AnchorId in the corresponding Chain.
    // The AnchorId of an exit is the second AnchorId in the corresponding Chain.
    class EntranceOrExit {
    public:
        AssemblyGraph::edge_descriptor e;
        AnchorId anchorId;

#if 0
        // The AnchorMarkerIntervals on that AnchorId.
        // These are initially copies from class Anchors.
        // But later, for entrances we remove AnchorMarkerIntervals
        // for which the same OrientedReadId appears in another entrance;
        // and for exits we remove AnchorMarkerIntervals
        // for which the same OrientedReadId appears in another exit.
        vector<AnchorMarkerInterval> anchorMarkerIntervals;

        // The AnchorIds encountered during read following.
        // For an entrance, read following moves forward, starting at the entrance.
        // For an exit, read following moves backward, starting at the exit.
        vector<AnchorId> journeyAnchorIds;

        // The AnchorIds encountered during read following starting from this Entrance
        // and no other entrance.
        vector<AnchorId> uniqueJourneyAnchorIds;
#endif

        EntranceOrExit(
            AssemblyGraph::edge_descriptor,
            AnchorId);
    };
    class Entrance : public EntranceOrExit {using EntranceOrExit::EntranceOrExit;};
    class Exit : public EntranceOrExit {using EntranceOrExit::EntranceOrExit;};

    vector<Entrance> entrances;
    vector<Exit> exits;

    // Find Assembly graph edges that are both an entrance and an exit.
    void findEntranceExits(vector<AssemblyGraph::edge_descriptor>&) const;

    private:

    void findEntrances();
    void findExits();
    bool isEntrance(AssemblyGraph::edge_descriptor) const;
    bool isExit(AssemblyGraph::edge_descriptor) const;
    bool isEntrance(AnchorId) const;
    bool isExit(AnchorId) const;
    void writeEntrances() const;
    void writeExits() const;

    // Use read following to fill in the journeyAnchorIds of each entrance/exit.
    void readFollowingFromEntrances();
    void readFollowingFromEntrance(Entrance&);
    void readFollowingFromExits();
    void readFollowingFromExit(Exit&);

};

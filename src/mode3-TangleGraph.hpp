#pragma once

// Class TangleGraph is used for read following in the AnchorGraph
// with a given set of entrances and exits.

// Shasta.
#include "mode3-Anchor.hpp"

// Standard library.
#include <cstdint.hpp>
#include <vector.hpp>

namespace shasta {
    namespace mode3 {
        class TangleGraph;
    }
}



class shasta::mode3::TangleGraph {
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

    // The tangle matrix is the number of common unique anchors between each
    // entrance/exit pair.
    // Indexed by [entranceIndex][exitIndex]
    vector< vector<uint64_t> > tangleMatrix;
    void computeTangleMatrix();

};

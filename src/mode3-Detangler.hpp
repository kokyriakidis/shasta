#pragma once

/******************************************************************

mode3::Detangler is an abstract class representing an object
that knows how to detangle (or at least try to) a superbubble
of an AssemblyGraph.

From the point of view of the detangler, a superbubble is simply
a set of vertices that is disjoint from any other superbubble.

If successful, operator() returns true and is only allowed
to make AssemblyGraph changes that affect superbubble vertice and edge,
plus its incoming/outgoing edges.

There are two derived abstract classes:
- ChainDetangler works on the incoming/outgoing chains of a superbubble.
- VertexDetangler works on the incoming/outgoing vertices of a superbubble.

******************************************************************/

// Shasta.
#include "mode3-Anchor.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class Detangler;
        class ChainDetangler;
        class VertexDetangler;

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            AssemblyGraphVertex,
            AssemblyGraphEdge>;
    }
}



class shasta::mode3::Detangler {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;
    using edge_descriptor = AssemblyGraphBaseClass::edge_descriptor;

    Detangler(bool debug, AssemblyGraph&);
    virtual bool operator()(const vector<vertex_descriptor>& superbubble) = 0;

protected:
    bool debug;
    AssemblyGraph& assemblyGraph;

    void writeInitialMessage(const vector<vertex_descriptor>& superbubble) const;

    void removeAllSuperbubbleVertices(const vector<vertex_descriptor>& superbubble) const;
};



// A detangler that works on the entrance/exit Chains of a superbubble.
class shasta::mode3::ChainDetangler: public Detangler {
public:

    ChainDetangler(bool debug, AssemblyGraph&);

    virtual bool operator()(const vector<vertex_descriptor>& superbubble) = 0;

protected:

    // This fills in entrances and exits and constructs the tangle matrix.
    void prepare(const vector<vertex_descriptor>& superbubble);

    class Entrance {
    public:

        // The AssemblyGraph edge for this Chain.
        // This must be the only Chain of this BubbleChain.
        edge_descriptor e;

        // The second to last AnchorId of that Chain.
        AnchorId anchorId = invalid<AnchorId>;

        // Total common coverage for this Entrance.
        // This is the sum of tangle matrix elements
        // pertaining to this EntranceChain.
        uint64_t commonCoverage = 0;

        Entrance(edge_descriptor, const AssemblyGraph&);
        Entrance() {}
    };

    class Exit {
    public:

        // The AssemblyGraph edge for this Chain.
        // This must be the only Chain of this BubbleChain.
        edge_descriptor e;

        // The second AnchorId of that Chain.
        AnchorId anchorId = invalid<AnchorId>;

        // Total common coverage for this ExitChain.
        // This is the sum of tangle matrix elements
        // pertaining to this EntranceChain.
        uint64_t commonCoverage = 0;

        Exit(edge_descriptor, const AssemblyGraph&);
        Exit() {}
    };

    vector<Entrance> entrances;
    vector<Exit> exits;
    void writeEntrancesAndExits() const;

    // The tangle matrix is computed using the second to last AnchorId of
    // each incoming Chain and the second AnchorId of each outgoing Chain.
    // Indexed by [iEntrance][iExit].
    vector< vector<uint64_t> > tangleMatrix;
    uint64_t totalCommonCoverage = 0;
    void computeTangleMatrix();
    void writeTangleMatrix() const;

    // Return true if there is one or more Entrance/Exit pair
    // with the same edge_descriptor.
    bool commonChainsBetweenEntrancesAndExitsExists() const;

    // Return true if there is one or more Entrance/Exit pair
    // with the same AnchorId.
    bool commonAnchorsBetweenEntrancesAndExitsExists() const;

    // Connect the Chains of an Entrance and Exit and
    // remove the AssemblyGraph edges for the entrance and exit.
    void connect(const Entrance&, const Exit&) const;

};



// A detangler that works on the entrance/exit vertices of a superbubble.
class shasta::mode3::VertexDetangler: public Detangler {
public:

};


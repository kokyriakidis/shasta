#pragma once

/******************************************************************

mode3::Detangler is an abstract class representing an object
that knows how to detangle (or at least try to) a superbubble
of an AssemblyGraph.

From the point of view of the detangler, a superbubble is simply
a set of vertices that is disjoint from any other superbubble.

If successful, operator() returns true and is only allowed
to make the AssemblyGraph that affect superbubble vertice and edge,
plus its incoming/outgoing edges.

******************************************************************/

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class Detangler;

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

    virtual bool operator()(
        AssemblyGraph&,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& superbubble
        ) = 0;
};

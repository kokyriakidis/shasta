#pragma once

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "vector.hpp"



namespace shasta {
    namespace mode3 {
        class Superbubble;
        class Superbubbles;

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;

        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            AssemblyGraphVertex,
            AssemblyGraphEdge>;

        class AssemblyGraphEdgePredicate;
    }
}



class shasta::mode3::Superbubble : public vector<AssemblyGraphBaseClass::vertex_descriptor> {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;
    vector<vertex_descriptor> entrances;
    vector<vertex_descriptor> exits;

    // Fill in the superbubble given a single entrance and exit.
    void fillInFromEntranceAndExit(const AssemblyGraph&);
};



class shasta::mode3::Superbubbles {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;

    // This computes connected components using only edges with length up to maxOffset1.
    Superbubbles(
        AssemblyGraph&,
        uint64_t maxOffset1
        );

    // This uses dominator trees.
    // It only finds superbubbles with one entrance and one exit.
    Superbubbles(AssemblyGraph&);

    // This constructs superbubbles consisting of a single tangled vertex.
    // An vertex v0 is tangle if:
    // inDegree(v0)>1, outDegree(v0)>1.
    class FromTangledVertices{};
    Superbubbles(AssemblyGraph&, const FromTangledVertices&);

    // This constructs superbubbles consisting of a single tangled edge.
    // An edge v0->v1 is tangled if
    // inDegree(v0)>1, outDegree(v0)==1, inDegree(v1)==1, outDegree(v1)>1.
    class FromTangledEdges{};
    Superbubbles(AssemblyGraph&, const FromTangledEdges&);

    // This computes connected components using the set of edges
    // for which the AssemblyGraphEdgePredicate returns true.
    // Each connected component with more than one vertex becomes a Superbubble.
    // This does not compute entrances and exits of each Superbubble.
    Superbubbles(AssemblyGraph&, const AssemblyGraphEdgePredicate&);

    // This creates an empty Superbubbles object.
    class Empty{};
    Superbubbles(AssemblyGraph&, const Empty&);

    ~Superbubbles();

    // Return the number of superbubbbles.
    uint64_t size() const
    {
        return superbubbles.size();
    }

    // Return the vertices in the specified superbubble.
    Superbubble& getSuperbubble(uint64_t superBubbleId)
    {
        return superbubbles[superBubbleId];
    }
    const Superbubble& getSuperbubble(uint64_t superBubbleId) const
    {
        return superbubbles[superBubbleId];
    }

    // Figure out if a vertex is in the specified superbubble.
    bool isInSuperbubble(uint64_t superbubbleId, vertex_descriptor cv) const;

    AssemblyGraph& assemblyGraph;

    // The superbubbles are the connected components with size at least 2,
    // computed using only the edges with offset up to maxOffset1.
    vector<Superbubble> superbubbles;
};




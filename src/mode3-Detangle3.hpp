#pragma once

/******************************************************************************

Class mode3::Detangle3Graph and related classes used by mode3::AssemblyGraph::detangle3.
This code requires all BubbleChains of the AssemblyGraph to consist of a single Chain.
It is intended to run immediately after construction of the AssemblGraph,
as in AssemblyGraph::run3.

In the Detangle3Graph, each vertex corresponds to a Segment=Chain=BubbleChain=Edge
of the Assembly graph. Only Chains with at least one internal AnchorId
generate a Detangle3Graph vertex.

Given two vertices v0 and v1 corresponding to Chains chain0 and chain1,
we create a directed edge v0->v1 if there is sufficient similarity
between the read compositions of the last AnchorId of Chain0
and the first AnchorId of Chain1.

*******************************************************************************/

// Shasta.
#include "mode3-AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>

namespace shasta {
    namespace mode3 {

        class Detangle3Graph;
        class Detangle3GraphVertex;
        class Detangle3GraphEdge;

        using Detangle3GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Detangle3GraphVertex,
            Detangle3GraphEdge>;

    }
}



class shasta::mode3::Detangle3GraphVertex {
public:
    AssemblyGraph::edge_descriptor e;

    // The first and last internal AnchorIds in the Chain.
    AnchorId firstInternalAncorId;
    AnchorId lastInternalAncorId;

    Detangle3GraphVertex(
        AssemblyGraph::edge_descriptor e,
        AnchorId firstInternalAncorId,
        AnchorId lastInternalAncorId
        ) :
        e(e),
        firstInternalAncorId(firstInternalAncorId),
        lastInternalAncorId(lastInternalAncorId)
        {}
};



class shasta::mode3::Detangle3GraphEdge {
public:
    AnchorPairInfo info;
    Detangle3GraphEdge(const AnchorPairInfo& info) : info(info) {}
};



class shasta::mode3::Detangle3Graph : public Detangle3GraphBaseClass {
public:

    Detangle3Graph(AssemblyGraph&);

private:
    AssemblyGraph& assemblyGraph;

    void createVertices();
    std::map<AssemblyGraph::edge_descriptor, vertex_descriptor> vertexMap;

    // Create all edges.
    void createEdges();

    // Create the edges starting at v0:
    // If direction is 0, move forward in the AssemblyGraph.
    // If direction is 1, move backward in the AssemblyGraph.
    void createEdges(vertex_descriptor v0, uint64_t direction);

    void writeGraphviz(const string& fileName) const;
};


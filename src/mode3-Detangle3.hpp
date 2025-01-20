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

    // Average coverage of the internal anchors.
    uint64_t coverage;

    // The total offset between the first and last internal anchors.
    uint64_t offset;

    Detangle3GraphVertex(
        AssemblyGraph::edge_descriptor e,
        AnchorId firstInternalAncorId,
        AnchorId lastInternalAncorId,
        uint64_t coverage,
        uint64_t offset
        ) :
        e(e),
        firstInternalAncorId(firstInternalAncorId),
        lastInternalAncorId(lastInternalAncorId),
        coverage(coverage),
        offset(offset)
        {}
};



class shasta::mode3::Detangle3GraphEdge {
public:

    // AnchorPairInfo between the last internal anchor of the source
    // vertex and the first internal anchor of the target vertex.
    AnchorPairInfo info;

    Detangle3GraphEdge(const AnchorPairInfo& info) : info(info) {}

    // hasMinimumOffset[0] gets set if this the edge with minimum
    // offset among all edges with the same source vertex.
    // hasMinimumOffset[1] gets set if this the edge with minimum
    // offset among all edges with the same target vertex.
    array<bool, 2> hasMinimumOffset = {false, false};

    // hasMaximumCommon[0] gets set if this the edge with maximum
    // number of common oriented reads among all edges with the same source vertex.
    // hasMaximumCommon[1] gets set if this the edge with maximum
    // number of common oriented reads among all edges with the same target vertex.
    array<bool, 2> hasMaximumCommon = {false, false};

    bool isStrong() const
    {
        return hasMaximumCommon[0] and hasMaximumCommon[1];
    }
};



class shasta::mode3::Detangle3Graph : public Detangle3GraphBaseClass {
public:

    Detangle3Graph(AssemblyGraph&);

private:
    AssemblyGraph& assemblyGraph;

    void createVertices(uint64_t minCoverage, uint64_t maxCoverage);
    void createVertex(AssemblyGraph::edge_descriptor e, uint64_t minCoverage, uint64_t maxCoverage);
    std::map<AssemblyGraph::edge_descriptor, vertex_descriptor> vertexMap;

    // Create all edges.
    void createEdges();

    // Graph cleanup functions.
    void removeStrongComponents();
    void removeWeakEdges();
    void removeIsolatedVertices();
    void transitiveReduction();

    // Create the edges starting at v0:
    // If direction is 0, move forward in the AssemblyGraph.
    // If direction is 1, move backward in the AssemblyGraph.
    void createEdges(vertex_descriptor v0, uint64_t direction);

    void writeGraphviz(const string& fileName) const;
};


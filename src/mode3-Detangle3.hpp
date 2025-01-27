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

    // The strong chain that this vertex is on, if any, and its position on the strong chain.
    uint64_t strongChainId = invalid<uint64_t>;
    uint64_t positionInStrongChain = invalid<uint64_t>;
    void resetStrongChain()
    {
        strongChainId = invalid<uint64_t>;
        positionInStrongChain = invalid<uint64_t>;
    }

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

    Detangle3GraphVertex() = default;
    Detangle3GraphVertex(const Detangle3GraphVertex&) = default;
};



class shasta::mode3::Detangle3GraphEdge {
public:

    // AnchorPairInfo between the last internal anchor of the source
    // vertex and the first internal anchor of the target vertex.
    AnchorPairInfo info;

    Detangle3GraphEdge(const AnchorPairInfo& info) : info(info) {}
    Detangle3GraphEdge() = default;
    Detangle3GraphEdge(const Detangle3GraphEdge&) = default;

    // hasMaximumCommon[0] gets set if this the edge with maximum
    // number of common oriented reads among all edges with the same source vertex.
    // hasMaximumCommon[1] gets set if this the edge with maximum
    // number of common oriented reads among all edges with the same target vertex.
    array<bool, 2> hasMaximumCommon = {false, false};

    bool isStrong() const
    {
        return hasMaximumCommon[0] and hasMaximumCommon[1];
    }

    bool isVeryWeak() const
    {
        return (not hasMaximumCommon[0]) and (not hasMaximumCommon[1]);
    }

    void resetHasMinimumOffsetFlags()
    {
        hasMaximumCommon = {false, false};
    }
};



class shasta::mode3::Detangle3Graph : public Detangle3GraphBaseClass {
public:

    Detangle3Graph(AssemblyGraph&);
    Detangle3Graph(const Detangle3Graph&);

private:
    AssemblyGraph& assemblyGraph;

    void createVertices(uint64_t minCoverage, uint64_t maxCoverage);
    void createVertex(AssemblyGraph::edge_descriptor e, uint64_t minCoverage, uint64_t maxCoverage);
    std::map<AssemblyGraph::edge_descriptor, vertex_descriptor> vertexMap;

    // Create all edges.
    void createEdges();

    // Set the hasMaximumCommon on all edges.
    void setHasMaximumCommonFlags();

    // Create the edges starting at v0:
    // If direction is 0, move forward in the AssemblyGraph.
    // If direction is 1, move backward in the AssemblyGraph.
    void createEdges(vertex_descriptor v0, uint64_t direction);

    // Graph cleanup functions.
    void removeStrongComponents();
    void removeWeakEdges();
    void removeVeryWeakEdges();
    void removeIsolatedVertices();
    void transitiveReduction();
    void prune(uint64_t maxLength);
    bool pruneIteration(uint64_t maxLength);

    void findLinearChains(vector< vector<vertex_descriptor> >&) const;
    void findLinearChains(vector< vector<edge_descriptor> >&) const;
    void findLinearChains(vector< vector<vertex_descriptor> >&, vector< vector<edge_descriptor> >&) const;

    // Find strong chains that can be used to stitch AssemblyGraph Chains together.
    vector< vector<vertex_descriptor> > strongChains;
    void findStrongChains(uint64_t maxPruneLength);

    // Use the strong chains to update the AssemblyGraph.
    void updateAssemblyGraph();
    void updateAssemblyGraph(const vector<vertex_descriptor>& strongChain);

    // Remove edges between vertices of the same strong chain,
    // except for the edges which form the strong chain itself.
    void removeInternalStrongChainEdges();

    // Remove edges incident on vertices that are internal to a strong chain,
    // except for the edges of strong chains themselves.
    void removeEdgesIncidentInsideStrongChains();

    string vertexStringId(vertex_descriptor) const;
    void writeGraphviz(const string& fileName) const;
};


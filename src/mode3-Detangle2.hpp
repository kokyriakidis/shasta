#pragma once

// Class mode3::Detangle2 and related classes used by mode3::AssemblyGraph::detangle2.

// Shasta.
#include "mode3-AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"

namespace shasta {
    namespace mode3 {

        class Detangle2Graph;
        class Detangle2GraphVertex;
        class Detangle2GraphEdge;

        using Detangle2GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Detangle2GraphVertex,
            Detangle2GraphEdge>;

    }
}



// A Detangle2GraphVertex corresponds to an AssemblyGraph which,
// by construction, is a BubbleChain consisting of a single, long Chain.
class shasta::mode3::Detangle2GraphVertex {
public:
    AssemblyGraph::edge_descriptor e;
    uint64_t offset;

    Detangle2GraphVertex(
        AssemblyGraph::edge_descriptor e,
        uint64_t offset) :
        e(e), offset(offset) {}
};



// A Detangle2GraphEdge v0->v1  (Chain0->Chain1) is created if we found:
// - A forward path starting at Chain0 and ending at Chain1.
// and/or
// - A backward path starting at Chain1 and ending at Chain0.
class shasta::mode3::Detangle2GraphEdge {
public:

    // Flags to indicate if we found the forward path and/or the backward path.
    array<bool, 2> found = {false, false};

    // The forward and backward paths found for this edge.
    array< vector<AnchorId>, 2> paths;
    const vector<AnchorId>& shortestPath() const
    {
        SHASTA_ASSERT(found[0] and found[1]);
        if(paths[0].size() < paths[1].size()) {
            return paths[0];
        } else {
            return paths[1];
        }
    }
};



class shasta::mode3::Detangle2Graph : public Detangle2GraphBaseClass {
public:

    // Construct the vertices from the long Chains of the AssemblyGraph.
    Detangle2Graph(
        AssemblyGraph&,
        uint64_t chainLengthThreshold);
    void addEdges();

    // Remove edges found in one direction only.
    void removeWeakEdges();

    void findLinearChains();
    vector< vector<edge_descriptor> > linearChains;

    // Use each linear chain to create a new Chain in the assembly graph.
    void createChains();
    void createChain(const vector<edge_descriptor>&);

    void writeGraphviz(const string& fileName) const;

private:
    AssemblyGraph& assemblyGraph;
    std::map<AssemblyGraph::edge_descriptor, vertex_descriptor> vertexMap;

    // Find forward and backward paths starting at the Chain corresponding
    // to a given AssemblyGraph edge, and generate edges when successful.
    void findForwardPath(vertex_descriptor);
    void findBackwardPath(vertex_descriptor);
    void findPath(vertex_descriptor, uint64_t direction);


    // This is used by findPath.
    class AnchorInfo {
    public:
        AnchorId anchorId;

        uint64_t totalOffset;   // Relative to the start anchors.

        AnchorInfo(
            AnchorId anchorId,
            uint64_t totalOffset) :
            anchorId(anchorId),
            totalOffset(totalOffset)
            {}

        // This controls the queueing order in findPath.
        // Not clear if < or > is better.
        bool operator<(const AnchorInfo& that) const {
            return totalOffset > that.totalOffset;
        }
    };
};


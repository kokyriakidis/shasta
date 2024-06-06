#pragma once

// Shasta.
#include "mode3-AssemblyGraph.hpp"


// The TangleGraph is a bipartite graph that stores
// the relationship between OrientedReadIds and MarkerGraphEdgeIds
// in a mode3::AssemblyGraph.

namespace shasta {
    namespace mode3 {

        class TangleGraph;
        class TangleGraphVertex;
        class TangleGraphEdge;

        using TangleGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            TangleGraphVertex,
            TangleGraphEdge>;

        enum class TangleGraphVertexType;
    }
}



enum class shasta::mode3::TangleGraphVertexType {
    OrientedRead,
    MarkerGraphEdge,
    Chain
};



class shasta::mode3::TangleGraphVertex {
public:
    TangleGraphVertexType type;

    // If type is OrientedRead, we store the index of the oriented read into
    // assemblyGraph.orientedReadIds;
    uint64_t orientedReadIndex;

    // If type is MarkerGraphEdge, we store the index of the marker graph edge into
    // assemblyGraph.markerGraphEdgeIds;
    uint64_t markerGraphEdgeIndex;

    // If type is Chain, we store the edge descriptor for the corresponding
    // AssemblyGraph edge.
    AssemblyGraph::edge_descriptor e;

    // In the first two constructors, the second argument is for dispatching only
    // andis ignored.
    TangleGraphVertex(uint64_t orientedReadIndex, OrientedReadId) :
        type(TangleGraphVertexType::OrientedRead),
        orientedReadIndex(orientedReadIndex) {}

    TangleGraphVertex(uint64_t markerGraphEdgeIndex, MarkerGraphEdgeId) :
        type(TangleGraphVertexType::MarkerGraphEdge),
        markerGraphEdgeIndex(markerGraphEdgeIndex) {}

    TangleGraphVertex(AssemblyGraph::edge_descriptor e) :
        type(TangleGraphVertexType::Chain),
        e(e) {}

};



class shasta::mode3::TangleGraphEdge {
public:
    uint64_t coverage = 0;;
};



class shasta::mode3::TangleGraph : public TangleGraphBaseClass {
public:
    TangleGraph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

};


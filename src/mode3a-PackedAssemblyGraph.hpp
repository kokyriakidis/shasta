#ifndef SHASTA_MODE3A_PACKED_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_PACKED_ASSEMBLY_GRAPH_HPP

// In the PackedAssemblyGraph, each vertex corresponds
// to a linear sequence of vertices (path) in the AssemblyGraph.
// A directed edge v0->v1 is created if the Jaccard similarity
// between the last assembly graph vertex of v0 and
// the first assembly graph vertex of v1 is high.

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3a {
        class PackedAssemblyGraph;
        class PackedAssemblyGraphVertex;
        class PackedAssemblyGraphEdge;
        class PackedAssemblyGraphJourneyEntry;

        using PackedAssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            PackedAssemblyGraphVertex, PackedAssemblyGraphEdge>;

        // Forward definitions to avoid including AssemblyGraph.hpp here.
        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            AssemblyGraphVertex, AssemblyGraphEdge>;
        class JourneyEntry;
    }
}



class shasta::mode3a::PackedAssemblyGraphVertex {
public:
    vector<AssemblyGraphBaseClass::vertex_descriptor> assemblyGraphVertices;
    vector<JourneyEntry> journeyEntries;
};



class shasta::mode3a::PackedAssemblyGraphEdge {
public:
    uint64_t coverage;
};



class shasta::mode3a::PackedAssemblyGraph : public PackedAssemblyGraphBaseClass {
public:
    PackedAssemblyGraph(
        const AssemblyGraph&,
        uint64_t minLinkCoverage1,
        uint64_t minLinkCoverage2,
        uint64_t minMarkerCount);
private:
    const AssemblyGraph& assemblyGraph;

    void createVertices(
        uint64_t minLinkCoverage1,
        uint64_t minMarkerCount);

    string vertexStringId(vertex_descriptor) const;

    // Oriented read journeys on the PackedAssemblyGraph.
    // Indexed by OrientedReadId.getValue();
    vector< vector<vertex_descriptor> > journeys;
    void computeJourneys();

    void createEdges(uint64_t minLinkCoverage2);

    void removeRoundTripEdges();

    void writeGraphviz() const;
};


#endif

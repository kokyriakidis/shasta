#ifndef SHASTA_MODE3A_JACCARD_GRAPH_HPP
#define SHASTA_MODE3A_JACCARD_GRAPH_HPP

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include <map>
#include "string.hpp"

namespace shasta {
    namespace mode3a {
        class JaccardGraph;
        class JaccardGraphVertex;
        class JaccardGraphEdge;

        using JaccardGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            JaccardGraphVertex, JaccardGraphEdge>;

        // Forward definition to avoid including AssemblyGraph.hpp here.
        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            AssemblyGraphVertex, AssemblyGraphEdge>;
    }
}



class shasta::mode3a::JaccardGraphVertex {
public:
    AssemblyGraphBaseClass::vertex_descriptor av;
};



class shasta::mode3a::JaccardGraphEdge {
public:
    double jaccard;
    bool keep;
};



class shasta::mode3a::JaccardGraph : public JaccardGraphBaseClass {
public:

    JaccardGraph(const AssemblyGraph& assemblyGraph);
    const AssemblyGraph& assemblyGraph;

    vertex_descriptor addVertex(AssemblyGraphBaseClass::vertex_descriptor);

    // This gives the JaccardGraph vertex corresponding to each
    // AssemblyGraph vertex.
    std::map<AssemblyGraphBaseClass::vertex_descriptor, vertex_descriptor> vertexMap;

    edge_descriptor addEdge(
        AssemblyGraphBaseClass::vertex_descriptor,
        AssemblyGraphBaseClass::vertex_descriptor,
        double jaccard);

    string vertexStringId(vertex_descriptor) const;

    void makeKnn(uint64_t m);

    void writeGraphviz(const string& fileName, double minJaccard) const;
    void writeGraphviz(ostream&, double minJaccard) const;
};

#endif

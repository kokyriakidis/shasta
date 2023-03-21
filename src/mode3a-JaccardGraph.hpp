#ifndef SHASTA_MODE3A_JACCARD_GRAPH_HPP
#define SHASTA_MODE3A_JACCARD_GRAPH_HPP

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include <map>
#include "memory.hpp"
#include "string.hpp"
#include "vector.hpp"

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

    // The length of the longest path ending here.
    uint64_t longestPathLength;

    bool isLongPathVertex = false;

    // This gets set if the vertex belongs to a non-trivial
    // strong connected component.
    bool isCyclic = false;
};



class shasta::mode3a::JaccardGraphEdge {
public:
    double jaccard;
    bool isLongPathEdge = false;

    // This is only used for displaying the JaccardGraph.
    bool display;

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


    // Compute large connected components.
    void computeConnectedComponents(
        uint64_t minComponentSize,
        vector< shared_ptr<JaccardGraph> >&
    );

    // Compute strongly connected components.
    void computeStronglyConnectedComponents(
        vector< vector<vertex_descriptor> >&
    );

    // Remove all vertices that belong to strogly connected components.
    void removeStronglyConnectedComponents();

    // Find long paths and mark edges and vertices on the long paths.
    bool findLongPaths(uint64_t minPathLength);

    // Mark the edges to be displayed.
    void markDisplayEdges();

    void writeGraphviz(const string& fileName, double minJaccard);
    void writeGraphviz(ostream&, double minJaccard);
};

#endif

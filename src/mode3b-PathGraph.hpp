#ifndef SHASTA_MODE3B_PATH_GRAPH_HPP
#define SHASTA_MODE3B_PATH_GRAPH_HPP

/*******************************************************************************

In the mode3b::PathGraph, each vertex corresponds to a primary edge of
of the marker graph, which is believed to correspond to a single copy
of sequence. It is characterized as follows:
- minPrimaryCoverage <= coverage <= maxPrimaryCoverage
- No duplicate oriented reads on the marker graph edge or its vertices.

A directed edge v0->v1 is generated if a sufficient number of oriented reads
visit the marker graph edge corresponding to v1 after
the marker graph edge corresponding to v0, without visiting other
primary marker graph edges in between.

*******************************************************************************/

// Shasta.
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Assembler;
    namespace mode3b {
        class PathGraph;
    }
}



class shasta::mode3b::PathGraph {
public:
    PathGraph(const Assembler&);
private:
    const Assembler& assembler;

    // Find out if a marker graph edge is a primary edge.
    uint64_t minPrimaryCoverage;
    uint64_t maxPrimaryCoverage;
    bool isPrimary(MarkerGraphEdgeId) const;

    // A table of all the primary marker graph edges.
    // Each entry in this table is a vertex of the PathGraph.
    // The index in this table is the vertexId.
    // The table is sorted.
    vector<MarkerGraphEdgeId> vertices;

    // Map MarkerGraphEdgeId to vertexIds.
    // This is indexed by the MarkerGraphEdgeId and contains invalid<uint64_t>
    // if that marker graph edge is not a primary marker graph edge.
    vector<uint64_t> vertexTable;

    // This fills in primaryEdges and the primaryEdgesTable.
    void findVertices();

    // Follow the reads to find edges of the PathGraph.
    // Each edge is stored as a pair of primaryIds (indexes into primaryEdges vector).
    vector< pair<uint64_t, uint64_t> > edges;
    vector<uint64_t> edgeCoverage;
    uint64_t minCoverage;
    void findEdges();

    // Write the entire PathGraph in graphviz format.
    void writeGraphviz() const;


    // Represent each connected component as a Graph.
    class Vertex {
    public:
        MarkerGraphEdgeId edgeId;
    };
    class Edge {
    public:
    };
    class Graph : public boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        Vertex,
        Edge> {
    public:
        std::map<MarkerGraphEdgeId, vertex_descriptor> vertexMap;
        void addVertex(MarkerGraphEdgeId);
        void addEdge(MarkerGraphEdgeId, MarkerGraphEdgeId);
    };
    vector<Graph> components;
    void createComponents();

    // Index the large components.
    // Sorted by decreasing size.
    uint64_t minComponentSize;
    vector< pair<uint64_t, uint64_t> > componentIndex; // (componentId, size)

};

#endif


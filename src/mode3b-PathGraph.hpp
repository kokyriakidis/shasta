#ifndef SHASTA_MODE3B_PATH_GRAPH_HPP
#define SHASTA_MODE3B_PATH_GRAPH_HPP

/*******************************************************************************

In the mode3b::PathGraph, each vertex corresponds to a primary edge of
of the marker graph, which is believed to correspond to a single copy
of sequence. It is characterized as follows:
- minPrimaryCoverage <= coverage <= maxPrimaryCoverage
- No duplicate oriented reads on the marker graph edge or its vertices.
- Is a branch edge (see isBranchEdge for details).

A directed edge v0->v1 is generated if a sufficient number of oriented reads
visit the marker graph edge corresponding to v1 after
the marker graph edge corresponding to v0, without visiting other
primary marker graph edges in between.

*******************************************************************************/

// Shasta.
#include "MarkerGraphEdgePairInfo.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
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

    uint64_t minCoverageA;  // Used by createEdges
    uint64_t minCoverageB;  // Used by recreate

    // A table of all the primary marker graph edges.
    // Each entry in this table is a vertex of the PathGraph.
    // The index in this table is the vertexId.
    // The table is sorted.
    // This cannot be names "vertices" due to naming conflicts caused by
    // boost graph macros.
    vector<MarkerGraphEdgeId> verticesVector;

    // This fills in the verticesVector.
    void createVertices();

    // The "journey" of each oriented read is the sequence of vertices it encounters.
    // It stores pairs (ordinal0, vertexId) for each oriented read, sorted by ordinal0.
    // The vertexId is the index in verticesVector.
    // Indexed by OrientedReadId::getValue.
    // Journeys are used to generate edges by "following the reads".
    vector < vector< pair<uint32_t, uint64_t> > > orientedReadJourneys;
    void computeOrientedReadJourneys();


    // Follow the reads to find edges of the PathGraph.
    // Each edge is stored as a pair of primaryIds (indexes into primaryEdges vector).
    // This cannot be names "edges" due to naming conflicts caused by
    // boost graph macros.
    vector< pair<uint64_t, uint64_t> > edgesVector;
    vector<uint64_t> edgeCoverage;
    void createEdges();

    // Write the entire PathGraph in graphviz format.
    void writeGraphviz() const;


    // Represent each connected component as a Graph.
    class Vertex {
    public:
        MarkerGraphEdgeId edgeId;
        bool wasProcessed = false;
    };
    class Edge {
    public:
        uint64_t coverage;
        MarkerGraphEdgePairInfo info;
        // Flat that is set if the marker graph edges corresponding to the two
        // vertice of this edge are adjacent in the mARKER GRAPH.
        bool adjacent;
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
        void addEdge(
            MarkerGraphEdgeId,
            MarkerGraphEdgeId,
            uint64_t coverage,
            const MarkerGraphEdgePairInfo&,
            bool adjacent);
        void writeGraphviz(uint64_t componentId, ostream&, uint64_t minCoverageA) const;
        void findLinearChains(
            uint64_t minChainLength,
            vector< vector<MarkerGraphEdgeId> >&);
    };
    vector<Graph> components;
    void createComponents();

    // Index the large components.
    // Sorted by decreasing size.
    uint64_t minComponentSize;
    vector< pair<uint64_t, uint64_t> > componentIndex; // (componentId, size)

    // The linear chains found in all connected components stored in the componentIndex.
    vector< vector<MarkerGraphEdgeId> > chains;
    void findChains(uint64_t minChainLength);

    // Find out if a marker graph edge is a branch edge.
    // A marker graph edge is a branch edge if:
    // - Its source vertex has more than one outgoing edge with coverage at least minPrimaryCoverage.
    // OR
    // - Its target vertex has more than one incoming edge with coverage at least minPrimaryCoverage.
    bool isBranchEdge(MarkerGraphEdgeId) const;

    // Create an updated version of the PathGraph from the chains, as follows:
    // - Only vertices that appear in chains are used.
    // - The oriented read journeys are recomputed.
    // - Edges are recreated as follows:
    //   1. Edges between successive vertices of each chain are created.
    //      These edges are guaranteed to have edge at least equal to minCoverageA.
    //   2. Edges bridging between two distinct chains are created if they have coverage
    //      at least equal to minCoverageB and join the end/begin of the two chains.
    void recreate();
};

#endif


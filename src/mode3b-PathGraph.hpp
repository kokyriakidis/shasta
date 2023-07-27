#ifndef SHASTA_MODE3B_PATH_GRAPH_HPP
#define SHASTA_MODE3B_PATH_GRAPH_HPP

/*******************************************************************************

In the mode3b::GlobalPathGraph, each vertex corresponds to a primary edge of
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

        // The global path graph.
        // Each vertex corresponds to a primary marker graph edge.
        class GlobalPathGraph;

        // A single connected component of the GlobalPathGraph.
        class PathGraphVertex;
        class PathGraphEdge;
        class PathGraph;
        using PathGraphGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            PathGraphVertex,
            PathGraphEdge>;

        // In the ChainGraph, each vertex represents a chain in a PathGraph.
        class ChainGraphVertex;
        class ChainGraphEdge;
        class ChainGraph;
        using ChainGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            ChainGraphVertex,
            ChainGraphEdge>;
    }
}



class shasta::mode3b::PathGraphVertex {
public:
    // The marker graph edge corresponding to this PathGraph vertex.
    MarkerGraphEdgeId edgeId;
};



class shasta::mode3b::PathGraphEdge {
public:
    uint64_t coverage;
    MarkerGraphEdgePairInfo info;

    // Flag that is set if the two marker graph edges corresponding to
    // the vertices of this PathGraph edge are adjacent.
    // This is only used to report this information in graphviz output.
    bool adjacent;

    // If this edge belongs to a chain, store the chainId and position.
    uint64_t chainId = invalid<uint64_t>;
    uint64_t positionInChain = invalid<uint64_t>;
};



class shasta::mode3b::PathGraph : public PathGraphGraphBaseClass {
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

    // Linear chains of edges.
    vector< vector<edge_descriptor> > chains;
    void findChains(
        double minCorrectedJaccard,
        uint64_t minTotalBaseOffset);
};



class shasta::mode3b::GlobalPathGraph {
public:
    GlobalPathGraph(const Assembler&);
private:
    const Assembler& assembler;

    // Find out if a marker graph edge is a primary edge.
    uint64_t minPrimaryCoverage;
    uint64_t maxPrimaryCoverage;
    bool isPrimary(MarkerGraphEdgeId) const;

    // The minimum coverage required to generate an edge when "following the reads".
    uint64_t minCoverage;

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

    // The connected components of the GlobalPathGraph.
    vector<PathGraph> components;
    void createComponents();

    // Index the large connected components.
    // Sorted by decreasing size.
    uint64_t minComponentSize;
    vector< pair<uint64_t, uint64_t> > componentIndex; // (componentId, size)

    // Find out if a marker graph edge is a branch edge.
    // A marker graph edge is a branch edge if:
    // - Its source vertex has more than one outgoing edge with coverage at least minPrimaryCoverage.
    // OR
    // - Its target vertex has more than one incoming edge with coverage at least minPrimaryCoverage.
    bool isBranchEdge(MarkerGraphEdgeId) const;
};



class shasta::mode3b::ChainGraphVertex {
public:
};



class shasta::mode3b::ChainGraphEdge {
public:
    MarkerGraphEdgePairInfo info;
};



class shasta::mode3b::ChainGraph : public ChainGraphBaseClass {
public:
    ChainGraph(
        const PathGraph&,
        const Assembler&,
        double minCorrectedJaccardForChain);
    void writeGraphviz(ostream&) const;
};



#endif


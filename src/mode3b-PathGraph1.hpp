#ifndef SHASTA_MODE3B_PATH_GRAPH1_HPP
#define SHASTA_MODE3B_PATH_GRAPH1_HPP

/*******************************************************************************

In the mode3b::GlobalPathGraph1, each vertex corresponds to a primary edge of
of the marker graph, which is believed to correspond to a single copy
of sequence. It is characterized as follows:
- minPrimaryCoverage <= coverage <= maxPrimaryCoverage
- No duplicate oriented reads on the marker graph edge or its vertices.

A directed edge v0->v1 is generated if:
- A sufficient number of oriented reads visit the marker graph edge
  corresponding to v1 shortly after the marker graph edge corresponding to v0.
- The corrected Jaccard similarity between v0 and v1 is sufficiently high.

*******************************************************************************/

// Shasta.
#include "MarkerGraphEdgePairInfo.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Assembler;
    namespace mode3b {

        // The global path graph.
        // Each vertex corresponds to a primary marker graph edge.
        class GlobalPathGraph1;

        // A single connected component of the GlobalPathGraph.
        class PathGraph1Vertex;
        class PathGraph1Edge;
        class PathGraph1;
        using PathGraph1GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            PathGraph1Vertex,
            PathGraph1Edge>;
    }
}



class shasta::mode3b::PathGraph1Vertex {
public:
    // The marker graph edge corresponding to this PathGraph1 vertex.
    MarkerGraphEdgeId edgeId;
};



class shasta::mode3b::PathGraph1Edge {
public:
    MarkerGraphEdgePairInfo info;
    bool keep;  // Used by PathGraph1::knn.
};



class shasta::mode3b::PathGraph1 : public PathGraph1GraphBaseClass {
public:

    std::map<MarkerGraphEdgeId, vertex_descriptor> vertexMap;
    void addVertex(MarkerGraphEdgeId);

    void addEdge(
        MarkerGraphEdgeId,
        MarkerGraphEdgeId,
        const MarkerGraphEdgePairInfo&);
    void writeGraphviz(uint64_t componentId, ostream&) const;

    // For each vertex, only keep the best k outgoing and k incoming edges.
    // "Best" as defined by correctedJaccard of the edges.
    void knn(uint64_t k);
};



class shasta::mode3b::GlobalPathGraph1 {
public:
    GlobalPathGraph1(const Assembler&);
private:
    const Assembler& assembler;

    // Find out if a marker graph edge is a primary edge.
    bool isPrimary(
        MarkerGraphEdgeId,
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage) const;

    // Find out if a marker graph edge is a branch edge.
    // A marker graph edge is a branch edge if:
    // - Its source vertex has more than one outgoing edge with coverage at least minEdgeCoverage.
    // OR
    // - Its target vertex has more than one incoming edge with coverage at least minEdgeCoverage.
    bool isBranchEdge(
        MarkerGraphEdgeId,
        uint64_t minEdgeCoverage) const;

    // Each vertex corresponds to a primary marker graph edge.
    // Store them here.
    // The index in this table is the vertexId.
    // The table is sorted.
    vector<MarkerGraphEdgeId> vertices;
    void createVertices(
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage);

    // The "journey" of each oriented read is the sequence of vertices it encounters.
    // It stores pairs (ordinal0, vertexId) for each oriented read, sorted by ordinal0.
    // The vertexId is the index in verticesVector.
    // Indexed by OrientedReadId::getValue.
    // Journeys are used to generate edges by "following the reads".
    vector < vector< pair<uint32_t, uint64_t> > > orientedReadJourneys;
    void computeOrientedReadJourneys();

    // Follow the reads to find edges of the PathGraph.
    // Each edge is stored as a pair of vertexIds (indexes into vertices vector).
    class Edge {
    public:
        uint64_t vertexId0;
        uint64_t vertexId1;
        MarkerGraphEdgePairInfo info;
    };
    vector<Edge> edges;
    void createEdges(
        uint64_t maxDistanceInJourney,
        uint64_t minEdgeCoverage,
        double minCorrectedJaccard);

    // Write the entire PathGraph in graphviz format.
    void writeGraphviz() const;

    // The connected components of the GlobalPathGraph1.
    // Stored sorted by decreasing size, as measured by number of vertices.
    vector< shared_ptr<PathGraph1> > components;
    void createComponents(
        double minCorrectedJaccard,
        uint64_t minComponentSize);

    // K-nn of each connected component:
    // for each vertex, keep only the k best outgoing and incoming edges,
    // as measured by correctedJaccard of each edge.
    // This can break contiguity of the connected component.
    void knn(uint64_t k);

    // Transitive reduction of each connected component.
    void transitiveReduction();

    // To create initial chains:
    // - Use only very strong edges (minCorrectedJaccard >= minCorrectedJaccard1).
    // - Compute connected components.
    // - Keep k best outgoing/incoming edges for each vertex.
    // - Transitive reduction.
    // - Longest path in each connected component.
    // This can cause contiguity breaks, which will be recovered later using
    // a more complete version of the GlobalPathGraph1.
    void createInitialChains(uint64_t minEstimatedLength);

    // Write each connected component in graphviz format.
    void writeGraphviz(const string& baseName) const;
};



#endif


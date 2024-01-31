#ifndef SHASTA_MODE3B_PATH_GRAPH_HPP
#define SHASTA_MODE3B_PATH_GRAPH_HPP

/*******************************************************************************

In the mode3b::GlobalPathGraph, each vertex corresponds to a primary edge of
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
#include "Base.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
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
        class GlobalPathGraph;
        class GlobalPathGraphVertex;
        class GlobalPathGraphEdge;

        // A subset of the GlobalPathGraph, for example
        // a single connected component, represented as a Boost graph.
        class PathGraphVertex;
        class PathGraphEdge;
        class PathGraph;
        using PathGraphGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            PathGraphVertex,
            PathGraphEdge>;

        class GlobalPathGraphDisplayOptions;

    }
}



// Class to control Graphviz output of GlobalPathGraph and PathGraph.
class shasta::mode3b::GlobalPathGraphDisplayOptions {
public:
    bool labels = true;
    bool tooltips = true;
    bool colorVertices = true;
    bool colorEdges = true;
    bool showNonTransitiveReductionEdges = true;

    // Thresholds for coloring by corrected Jaccard similarity J'.
    // If J' <= redJ, the edge is drawn red.
    // If J' >= greenJ, the edge is drawn green.
    // For values in between, the color is interpolated.
    double redJ;
    double greenJ;

    GlobalPathGraphDisplayOptions(double redJ = 0., double greenJ = 1.) :
        redJ(redJ), greenJ(greenJ) {}

    void makeCompact()
    {
        labels = false;
        tooltips = false;
        colorVertices = false;
        colorEdges = false;
    }
};



class shasta::mode3b::PathGraphVertex {
public:

    // The corresponding GlobalPathGraph vertexId.
    uint64_t vertexId;

    // The corresponding marker graph edgeId.
    MarkerGraphEdgeId edgeId;
};



class shasta::mode3b::PathGraphEdge {
public:
    MarkerGraphEdgePairInfo info;
    uint64_t coverage;
    bool isNonTransitiveReductionEdge = false;
};



// A subset of the GlobalPathGraph, for example
// a single connected component, represented as a Boost graph.
class shasta::mode3b::PathGraph : public PathGraphGraphBaseClass {
public:

    std::map<MarkerGraphEdgeId, vertex_descriptor> vertexMap;
    void addVertex(
        uint64_t vertexId,
        MarkerGraphEdgeId);

    void addEdge(
        MarkerGraphEdgeId,
        MarkerGraphEdgeId,
        const MarkerGraphEdgePairInfo&,
        uint64_t coverage);

    void writeGraphviz(
        const string& name,
        const GlobalPathGraphDisplayOptions&) const;

    // Create the connected components of this PathGraph,
    // without changing the PathGraph itself.
    vector< shared_ptr<PathGraph> > createConnectedComponents(uint64_t minComponentSize) const;

    void localTransitiveReduction(
        uint64_t distance,
        uint64_t maxCoverage);

    // Remove cross-edges.
    // This removes an edge v0->v1 if the following are all true:
    // - Its coverage is at most lowCoverageThreshold.
    // - Its estimated offset is at least minOffset.
    // - v0 has at least one out-edge with coverage at least highCoverageThreshold.
    // - v1 has at least one in-edge with coverage at least highCoverageThreshold.
    void removeCrossEdges(
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold,
        uint64_t minOffset);

};



class shasta::mode3b::GlobalPathGraphVertex {
public:
    MarkerGraphEdgeId edgeId;

    GlobalPathGraphVertex(MarkerGraphEdgeId edgeId) : edgeId(edgeId) {}

    // Compare by edgeId only.
    bool operator<(const GlobalPathGraphVertex& that) const
    {
        return edgeId < that.edgeId;
    }
};




class shasta::mode3b::GlobalPathGraphEdge {
public:
    uint64_t vertexId0;
    uint64_t vertexId1;
    MarkerGraphEdgePairInfo info;

    // The number of oriented reads for which the journey contains
    // a transition from vertexId0 to vertexId1.
    // This is less than or equal to info.common.
    uint64_t coverage = invalid<uint64_t>;
};



class shasta::mode3b::GlobalPathGraph {
public:
    static void assemble(
        const Assembler&,
        uint64_t threadCount0,
        uint64_t threadCount1);
    static void loadAndAssemble(
        const Assembler&,
        const string& fileName,
        uint64_t threadCount0,
        uint64_t threadCount1);
private:
    GlobalPathGraph(const Assembler&);
    const Assembler& assembler;

    // Each vertex corresponds to a primary marker graph edge.
    // Store them here.
    // The index in this table is the vertexId.
    // The table is sorted by MarkerGraphEdgeId.
    // It cannot be named "vertices" because of name conflicts with the Boost graph
    // iteration macros.
    vector<GlobalPathGraphVertex> verticesVector;
    void createVertices();

    // Return the vertexId corresponding to a given MarkerGraphEdgeId, or
    // invalid<MarkerGraphEdgeId> if no such a vertex exists.
    // This uses a binary search in the verticesVector (via std::lower_bound).
    uint64_t getVertexId(MarkerGraphEdgeId) const;

    vector<GlobalPathGraphEdge> edges;
    void createEdges();

    // The connected components of the GlobalPathGraph.
    // Stored sorted by decreasing size, as measured by number of vertices.
    vector< shared_ptr<PathGraph> > components;

    // Create connected components.
    // This only considers edges with corrected Jaccard at least equal to
    // minCorrectedJaccard, and only stores connected components with at
    // least minComponentSize vertices.
    void createComponents(
        double minCorrectedJaccard,
        uint64_t minComponentSize);

    // Local transitive reduction of each connected component.
    // Only edges with coverage up to maxCoverage are subject to flagging
    // as removed during transitive reduction.
    void localTransitiveReduction(
        uint64_t distance,
        uint64_t maxCoverage);

    // Write the entire PathGraph in graphviz format.
    void writeGraphviz() const;

    // Write each connected component in graphviz format.
    void writeComponentsGraphviz(
        const string& baseName,
        const GlobalPathGraphDisplayOptions&) const;

    void assembleComponent(
        uint64_t componentId,
        uint64_t transitiveReductionDistance,
        uint64_t transitiveReductionMaxCoverage,
        uint64_t crossEdgesLowCoverageThreshold,
        uint64_t crossEdgesHighCoverageThreshold,
        uint64_t crossEdgesMinOffset,
        uint64_t threadCount0,
        uint64_t threadCount1);
};



#endif


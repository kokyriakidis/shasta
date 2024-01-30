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
        class GlobalPathGraph1;
        class GlobalPathGraph1Vertex;
        class GlobalPathGraph1Edge;

        // A subset of the GlobalPathGraph1, for example
        // a single connected component, represented as a Boost graph.
        class PathGraph1Vertex;
        class PathGraph1Edge;
        class PathGraph1;
        using PathGraph1GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            PathGraph1Vertex,
            PathGraph1Edge>;

        class GlobalPathGraph1DisplayOptions;

        class CompressedPathGraph1A;
        class CompressedPathGraph1AVertex;
        class CompressedPathGraph1AEdge;
        using CompressedPathGraph1ABaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            CompressedPathGraph1AVertex,
            CompressedPathGraph1AEdge>;

        class AssemblyPath;
    }
}



// Class to control Graphviz output of GlobalPathGraph1 and PathGraph1.
class shasta::mode3b::GlobalPathGraph1DisplayOptions {
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

    GlobalPathGraph1DisplayOptions(double redJ = 0., double greenJ = 1.) :
        redJ(redJ), greenJ(greenJ) {}

    void makeCompact()
    {
        labels = false;
        tooltips = false;
        colorVertices = false;
        colorEdges = false;
    }
};



class shasta::mode3b::PathGraph1Vertex {
public:

    // The corresponding GlobalPathGraph1 vertexId.
    uint64_t vertexId;

    // The corresponding marker graph edgeId.
    MarkerGraphEdgeId edgeId;
};



class shasta::mode3b::PathGraph1Edge {
public:
    MarkerGraphEdgePairInfo info;
    uint64_t coverage;
    bool isNonTransitiveReductionEdge = false;
};



// A subset of the GlobalPathGraph1, for example
// a single connected component, represented as a Boost graph.
class shasta::mode3b::PathGraph1 : public PathGraph1GraphBaseClass {
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
        const vector<GlobalPathGraph1Vertex>& globalVertices,
        const string& name,
        const GlobalPathGraph1DisplayOptions&) const;

    // For each vertex, only keep the best k outgoing and k incoming edges.
    // "Best" as defined by correctedJaccard of the edges.
    void knn(uint64_t k);

    // Create the connected components of this PathGraph1,
    // without changing the PathGraph1 itself.
    vector< shared_ptr<PathGraph1> > createConnectedComponents(uint64_t minComponentSize) const;

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



class shasta::mode3b::GlobalPathGraph1Vertex {
public:
    MarkerGraphEdgeId edgeId;

    GlobalPathGraph1Vertex(MarkerGraphEdgeId edgeId) : edgeId(edgeId) {}

    // Information on the oriented reads that visit this vertex.
    class JourneyInfoItem {
    public:
        OrientedReadId orientedReadId;
        uint64_t positionInJourney;
    };
    vector<JourneyInfoItem> journeyInfoItems;

    // Compare by edgeId only.
    bool operator<(const GlobalPathGraph1Vertex& that) const
    {
        return edgeId < that.edgeId;
    }
};




class shasta::mode3b::GlobalPathGraph1Edge {
public:
    uint64_t vertexId0;
    uint64_t vertexId1;
    MarkerGraphEdgePairInfo info;

    // The number of oriented reads for which the journey contains
    // a transition from vertexId0 to vertexId1.
    // This is less than or equal to info.common.
    uint64_t coverage = invalid<uint64_t>;
};



class shasta::mode3b::GlobalPathGraph1 {
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
    GlobalPathGraph1(const Assembler&);
    const Assembler& assembler;

    // Find out if a marker graph edge is a primary edge.
    bool isPrimary(
        MarkerGraphEdgeId,
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage) const;

    // Each vertex corresponds to a primary marker graph edge.
    // Store them here.
    // The index in this table is the vertexId.
    // The table is sorted by MarkerGraphEdgeId.
    // It cannot be named "vertices" because of name conflicts with the Boost graph
    // iteration macros.
    vector<GlobalPathGraph1Vertex> verticesVector;
    void createVertices(
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage);

    // Return the vertexId corresponding to a given MarkerGraphEdgeId, or
    // invalid<MarkerGraphEdgeId> if no such a vertex exists.
    // This uses a binary search in the verticesVector (via std::lower_bound).
    uint64_t getVertexId(MarkerGraphEdgeId) const;

    // The "journey" of each oriented read is the sequence of vertices it encounters.
    // It stores pairs (ordinal0, vertexId) for each oriented read, sorted by ordinal0.
    // The vertexId is the index in verticesVector.
    // Indexed by OrientedReadId::getValue.
    // Journeys are used to generate edges by "following the reads".
    vector < vector< pair<uint32_t, uint64_t> > > orientedReadJourneys;
    void computeOrientedReadJourneys();

    vector<GlobalPathGraph1Edge> edges;
    void createEdges0(
        uint64_t maxDistanceInJourney,
        uint64_t minEdgeCoverage,
        double minCorrectedJaccard);

    // Find children edges of vertexId0.
    // The first element of each pair of the children vector
    // is the vertexId of the child vertex.
    void findChildren(
        uint64_t vertexId0,
        uint64_t minEdgeCoverage,
        double minCorrectedJaccard,
        vector< pair<uint64_t, MarkerGraphEdgePairInfo> >& children);
    void findParents(
        uint64_t vertexId0,
        uint64_t minEdgeCoverage,
        double minCorrectedJaccard,
        vector< pair<uint64_t, MarkerGraphEdgePairInfo> >& parents);

    // The connected components of the GlobalPathGraph1.
    // Stored sorted by decreasing size, as measured by number of vertices.
    vector< shared_ptr<PathGraph1> > components;

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
        const GlobalPathGraph1DisplayOptions&) const;

    static void assemble2(
        const Assembler&,
        uint64_t threadCount0,
        uint64_t threadCount1);
    void assemble2(
        uint64_t componentId,
        uint64_t transitiveReductionDistance,
        uint64_t transitiveReductionMaxCoverage,
        uint64_t crossEdgesLowCoverageThreshold,
        uint64_t crossEdgesHighCoverageThreshold,
        uint64_t crossEdgesMinOffset,
        uint64_t threadCount0,
        uint64_t threadCount1);
};



// A better compressed representation of a PathGraph1.
// Each linear sequence of edges, without branches, of the PathGraph1
// after transitive reduction
// is compressed to a single edge in the CompressedPathGraph1A.
// Therefore, initially:
// - Each edge of the CompressedPathGraph1A corresponds to a linear sequence of edges,
//   without branches, of the PathGraph1 after transitive reduction.
// - Each vertex of the CompressedPathGraph1A corresponds to a vertex of the PathGraph1.
// - But not all vertices of the PathGraph1 have a corresponding vertex in the CompressedPathGraph1A.
// Later, during detangling, edges can be merged, so an edge in the CompressedPathGraph1A
// no longer corresponds to a linear sequence of edges in the PathGraph1.
// For this reason a CompressedPathGraph1AEdge stores a sequence of MarkerGraphedgeIds.


class shasta::mode3b::CompressedPathGraph1AVertex {
public:
    PathGraph1::vertex_descriptor v;

    // The id of the ChokePointChain that this vertex belongs to, if any.
    uint64_t chokePointChainId = invalid<uint64_t>;
};



class shasta::mode3b::CompressedPathGraph1AEdge {
public:
    uint64_t id = invalid<uint64_t>;
    vector<MarkerGraphEdgeId> chain;

    // The id of the ChokePointChain that this edge belongs to, if any.
    uint64_t chokePointChainId = invalid<uint64_t>;
};



class shasta::mode3b::CompressedPathGraph1A :
    public CompressedPathGraph1ABaseClass,
    public MultithreadedObject<CompressedPathGraph1A> {
public:
    CompressedPathGraph1A(
        const PathGraph1&,
        uint64_t componentId,
        const Assembler&);

private:
    // Information stored by the constructor.
    const PathGraph1& graph;
    uint64_t componentId;
    const Assembler& assembler;

    // The CompressedPathGraph1A corresponding to a given PathGraph1 vertex.
    // Not all PathGraph1 vertices have a corresponding CompressedPathGraph1A vertex.
    std::map<PathGraph1::vertex_descriptor, vertex_descriptor> vertexMap;

    uint64_t nextEdgeId = 0;
    string edgeStringId(edge_descriptor) const;

    // Initial creation from the PathGraph1.
    void create();

    // Transitive reduction.
    // This removes an edge cv0->cv1 if:
    // - It corresponds to a single PathGraph1 edge with coverage at most maxCoverage.
    // - cv1 is reachable from cv0 with a path that:
    //   * Has length at most maxDistance, as measured by number of edges in the CompressedPathGraph1A.
    //   * Does not use the cv0->cv1.
    // Edges are considered in order of increasing coverage.
    void transitiveRedution(
        uint64_t maxCoverage,
        uint64_t maxDistance
    );

    // Detangling.
    void detangle(
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh,
        uint64_t pathLengthForChokePoints,
        uint64_t maxBubbleIndexDelta,
        uint64_t transitiveRedutionMaxCoverage,
        uint64_t transitiveRedutionMaxDistance);
    uint64_t detangleEdges(
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh
        );
    bool detangleEdge(
        edge_descriptor,
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh,
        vector<edge_descriptor>& removedEdges);
    /*
    uint64_t detangleBubbleChains(
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh
        );
    */
    void detangleUsingChokePoints(
        uint64_t pathLengthForChokePoints,
        uint64_t maxBubbleIndexDelta,
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh,
        bool debug
        );
    class ChokePointChain;
    uint64_t detangleChokePointChain(
        ChokePointChain&,
        uint64_t chokePointChainId,
        uint64_t maxBubbleIndexDelta,
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh,
        bool debug);



    // The tangle matrix defined by two vertices.
    class TangleMatrix {
    public:
        // The in-edges of cv0.
        vector<edge_descriptor> inEdges;

        // The out-edges of cv1.
        vector<edge_descriptor> outEdges;

        // m[i][j] contains the number of common oriented reads
        // between the i-th in-edge and the j-th out-edge,
        // computed using the second to last PathGraph1 vertex of each in-edge
        // and the second PathGraph1 vertex of each out-edge.
        vector< vector<uint64_t> > m;

        uint64_t inDegree() const;
        uint64_t outDegree() const;

        // Analyze the TangleMatrix.
        // A matrix element is considered negigible and treated as zero if it is <= lowThreshold.
        // A matrix element is considered significantly meaningful if it is >= highThreshold.
        // Otherwise a matrix element is considered ambiguous.
        void analyze(
            uint64_t lowThreshold,
            uint64_t highThreshold,
            int64_t& phase,
            uint64_t& minConcordant,
            uint64_t& maxDiscordant) const;
    };
    void computeTangleMatrix(
        vertex_descriptor,
        vertex_descriptor,
        TangleMatrix&,
        bool debug
        ) const;
    void computeTangleMatrix(
        const vector<edge_descriptor>& inEdges,
        const vector<edge_descriptor>& outEdges,
        TangleMatrix&,
        bool debug
        ) const;
    // This version only fills in m. The inEdges and out_edgesmust have already been filled in.
    void computeTangleMatrix(TangleMatrix&, bool debug) const;
    void writeTangleMatrix(const TangleMatrix&) const;


#if 0
    // Classes used by detangleBubbleChains.
    // Two successive Bubbles in a bubble chains are usually separated by an edge,
    // but the edge can be missing.
    class Bubble {
    public:
        vertex_descriptor source;
        vertex_descriptor target;
        vector<edge_descriptor> edges;
        uint64_t degree() const {
            return edges.size();
        }
        Bubble* previousBubble = 0;
        Bubble* nextBubble = 0;
    };
    class InterBubble {
    public:
        edge_descriptor edge;
        bool edgeExists;    // If false, edge is not valid.
    };
    class BubbleChain {
    public:
        vector<const Bubble*> bubbles;
        vector<InterBubble> interBubbles;   // One less that the number of bubbles.
    };
    void findBubbles(vector<Bubble>&) const;
    void findBubbleChains(vector<BubbleChain>&) const;
#endif



    // Classes and functions used by detangleUsingChokePoints.

    // A ChokePointChain describes a linear chain of "choke points" in the CompressedPathGraph1A
    // and the intervening "superbubbles" in-between.
    // A superbubble can have arbitrary complexity but will often consist of a single edge
    // or a bubble, most commonly a diploid bubble. Diploid bubbles are used for detangling.
    class Superbubble {
    public:
        vector<vertex_descriptor> internalVertices;
        vector<edge_descriptor> internalEdges;
        bool isDiploidBubble = false;
        vector<edge_descriptor> diploidEdges; // Only if isDiploidBubble is true;
    };
    class ChokePointChain {
    public:
        vector<vertex_descriptor> chokePoints;
        vector<Superbubble> superbubbles;       // Size is chokePoints.size() - 1

        // The indexes of the diploid bubbles in the superbubbles vector.
        // Indexes into this vector can be used as ids for diploid bubbles during phasing.
        vector<uint64_t> diploidBubblesIndexes;

        void getAllVertices(vector<vertex_descriptor>&) const;
        void getAllEdges(vector<edge_descriptor>&) const;

        // Flag that indicates this ChokePointChain overlaps another, larger ChokePointChain,
        // and should be discarded.
        bool discard = false;
    };
    void findChokePointChains(
        uint64_t pathLengthForChokePoints,
        vector<ChokePointChain>&,
        bool debug) const;
    // void analyzeChokePoints() const;
    void findVerticesAndEdgesBetweenChokePoints(
        vertex_descriptor,
        vertex_descriptor,
        vector<vertex_descriptor>&,
        vector<edge_descriptor>&) const;
    void flagOverlappingChokePointChains(
        vector<ChokePointChain>&,
        bool debug);



    // The PhasingGraph is used by detangleUsingChokePoints to phase
    // the diploid bubbles in a ChokePointChain.
    // It is an undirected graph in which each vertex represents a diploid bubble.
    // It uses vecS for the vertex container, so the vertex_descriptor
    // is the same as the bubble id in the ChokePointChain
    // (index in bubbleChainIndexes).
    class PhasingGraphVertex {
    public:
        int64_t phase = 0;  // +1 or -1 for phase vertices, 0 otherwise
        bool isInLargestComponent = false;
    };
    class PhasingGraphEdge {
    public:
        int64_t phase;          // +1 (in phase) or -1 (out of phase)

        // Tangle matrix metrics.
        // If phase = +1, minConcordant = min(m00, m11), maxDiscordant = max(m01, m10).
        // If phase = -1, minConcordant = min(m01, m10), maxDiscordant = max(m00, m11).
        uint64_t minConcordant;
        uint64_t maxDiscordant;

        bool operator<(const PhasingGraphEdge& that) const
        {
            if(maxDiscordant < that.maxDiscordant) {
                return true;
            }
            if(maxDiscordant > that.maxDiscordant) {
                return false;
            }
            return minConcordant > that.minConcordant;
        }
        bool isSpanningTreeEdge = false;
    };
    using PhasingGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::undirectedS,
        PhasingGraphVertex,
        PhasingGraphEdge>;
    class PhasingGraph : public PhasingGraphBaseClass {
    public:
        PhasingGraph(uint64_t bubbleCount);
        void addEdge(uint64_t bubbleId0, uint64_t bubbleId1, const PhasingGraphEdge& edge);

        // Assign a phase to a subset of the vertices.
        // Return the number of inconsistent edges.
        uint64_t phase();

        bool isConsistent(edge_descriptor) const;

        void writeGraphviz(ostream&, const string& name) const;
    };



    // Accessors.

    // Get the vertex_descriptor corresponding to a PathGraph1::vertex_descriptor,
    // adding a vertex if necessary.
    vertex_descriptor getCompressedVertex(PathGraph1::vertex_descriptor);

    // Get the MarkerGraphEdgeId corresponding to a given vertex.
    MarkerGraphEdgeId markerGraphEdgeId(vertex_descriptor) const;

    // Get MarkerGraphEdgeIds at the beginning and end of each edge.
    MarkerGraphEdgeId firstMarkerGraphEdgeId(edge_descriptor) const;
    MarkerGraphEdgeId lastMarkerGraphEdgeId(edge_descriptor) const;
    MarkerGraphEdgeId secondMarkerGraphEdgeId(edge_descriptor) const;
    MarkerGraphEdgeId secondToLastMarkerGraphEdgeId(edge_descriptor) const;

    uint64_t totalBaseOffset(edge_descriptor) const;



    // Output.
    void writeGraphviz(const string& fileNamePrefix) const;
    void writeGfa(const string& fileNamePrefix) const;
    void writeGfaAndGraphviz(const string& fileNamePrefix) const;
    // void writeBubble(const Bubble&, ostream&) const;
    void writeChokePointChain(const ChokePointChain&) const;
};


#endif


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

        // A compressed representation of a PathGraph1.
        // A vertex of the CompressedPathGraph1 corresponds to one of the following:
        // - A linear chain of PathGraph1 vertices with in-degree and out-degree 0 or 1.
        // OR
        // - A single PathGraph1 vertex with in-degree>1 and/or out-degree>1.
        // An edge of the CompressedPathGraph1 always corresponds to an edge of the PathGraph1.
        class CompressedPathGraph1;
        class CompressedPathGraph1Vertex;
        class CompressedPathGraph1Edge;
        using CompressedPathGraph1BaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            CompressedPathGraph1Vertex,
            CompressedPathGraph1Edge>;

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

    void localTransitiveReduction(uint64_t distance);

};



class shasta::mode3b::GlobalPathGraph1Vertex {
public:
    MarkerGraphEdgeId edgeId;
    uint64_t chainId = invalid<uint64_t>;
    uint64_t positionInChain = invalid<uint64_t>;
    bool isFirstInChain = false;
    bool isLastInChain = false;

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
    // Only filled in by createEdges0.
    uint64_t coverage = invalid<uint64_t>;
};



class shasta::mode3b::GlobalPathGraph1 {
public:
    static void assemble(
        const Assembler&,
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
    // The table is sorted by MarkerGraphEdgeId.
    // It cannot be names "vertices" because of name conflicts with the Boost graph
    // iteration macros.
    vector<GlobalPathGraph1Vertex> verticesVector;
    void createVertices(
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage);

    // Return the vertexId corresponding to a given MarkerGraphEdgeId, or
    // invalid<MarkerGraphEdgeId> if no such a vertex exists.
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
    void createEdges1(
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

    // K-nn of each connected component:
    // for each vertex, keep only the k best outgoing and incoming edges,
    // as measured by correctedJaccard of each edge.
    // This can break contiguity of the connected component.
    void knn(uint64_t k);

    // Local transitive reduction of each connected component.
    void localTransitiveReduction(uint64_t distance);



    // A chain is a sequence of vertices (not necessarily a path)
    // that can be used to generate an AssemblyPath.
    class Chain {
    public:

        // The vertexIds (indices in the verticesVector) of
        // the vertices of this chain.
        vector<uint64_t> vertexIds;

        // The corresponding MarkerGraphEdgeIds.
        vector<MarkerGraphEdgeId> edgeIds;

        // The MarkerGraphEdgePairInfos in between the vertices.
        // infos.size() is always vertexIds.size() - 1;
        vector<MarkerGraphEdgePairInfo> infos;

        uint64_t totalOffset() const;
    };



    // For each connected component, use the longest path
    // to create a Chain. Only keep the ones that are sufficiently long.
    // This also stores chain information in the vertices.
    void createChainsFromComponents(
        uint64_t minEstimatedLength,
        vector<Chain>&);

    // Seed chains.
    vector<Chain> seedChains;
    void writeSeedChains() const;
    void writeSeedChainsDetails() const;
    void writeSeedChainsStatistics() const;

    // Store in this PathGraph1 the seed chains found in another PathGraph1
    void storeSeedChains(const vector<Chain>&);

    // Connect seed chains by following reads.
    void connectSeedChains0();

    // Connect seed chains by walking the graph.
    class ChainConnector {
    public:
        uint64_t chainId0;
        uint64_t chainId1;
        vector<uint64_t> vertexIds;
        vector<MarkerGraphEdgePairInfo> infos;
        void reverse();
        uint64_t totalOffset() const;
    };
    void connectSeedChains1(
        uint64_t minEdgeCoverage,
        double minCorrectedJaccard,
        vector<ChainConnector>&
        );
    void connectSeedChain1(
        uint64_t chainId,
        uint64_t direction, // 0 = forward, 1 = backward
        uint64_t minEdgeCoverage,
        double minCorrectedJaccard,
        vector<ChainConnector>&,
        ostream& out
        );
    void connectSeedChains2(
        uint64_t minCommonCount,
        double minCorrectedJaccard,
        vector<ChainConnector>&);
    void connectSeedChains2(
        uint64_t componentId,
        const PathGraph1& component,
        uint64_t minCommonCount,
        double minCorrectedJaccard,
        ostream&,
        vector<ChainConnector>&);
    void writeConnectors(const vector<ChainConnector>&) const;

    // Extend a chain forward until we bump into another chain.
    // This returns the chainId of the chain we found, or invalid<uint64_t>
    // if none found.
    uint64_t extendChainForward(
        uint64_t chainId,
        Chain&,
        const PathGraph1& component,
        uint64_t minCommonCount,
        double minCorrectedJaccard) const;

    // Use the ChainConnectors to stitch together the seed chains.
    // The replaces the connected components with the connected
    // components of the stitched graph.
    void stitchSeedChains(
        const vector<ChainConnector>&,
        uint64_t minComponentSize);

    // This generates an AssemblyPath for each of the Chains passed in,
    // then assembles the AssemblyPath and writes assembled sequence to fasta.
    void assembleChains(
        const vector<Chain>&,
        ostream& fasta,
        const string& csvPrefix) const;

    // Write the entire PathGraph in graphviz format.
    void writeGraphviz() const;

    // Write each connected component in graphviz format.
    void writeComponentsGraphviz(
        const string& baseName,
        const GlobalPathGraph1DisplayOptions&) const;

    static void assemble0(const Assembler&);

    static void assemble1(
        const Assembler&,
        uint64_t threadCount0,
        uint64_t threadCount1);
    static void assemble1(
        GlobalPathGraph1&,
        uint64_t threadCount0,
        uint64_t threadCount1,
        uint64_t componentId,
        uint64_t transitiveReductionDistance);

};



// A compressed representation of a PathGraph1.
// A vertex of the CompressedPathGraph1 corresponds to one of the following:
// - A linear chain of PathGraph1 vertices with in-degree and out-degree 0 or 1.
// OR
// - A single PathGraph1 vertex with in-degree>1 and/or out-degree>1.
class shasta::mode3b::CompressedPathGraph1Vertex {
public:
    uint64_t id;
    vector<PathGraph1::vertex_descriptor> v;
    shared_ptr<AssemblyPath> assemblyPath;
};
class shasta::mode3b::CompressedPathGraph1Edge {
public:
    // MarkerGraphEdgePairInfo computed between the last marker graph
    // edge of the source vertex and the first marker graph edge of the target vertex.
    MarkerGraphEdgePairInfo info;
};



class shasta::mode3b::CompressedPathGraph1 :
    public CompressedPathGraph1BaseClass,
    public MultithreadedObject<CompressedPathGraph1> {
public:
    CompressedPathGraph1(
        const PathGraph1&,
        uint64_t componentId,
        const Assembler&);

    // Information stored by the constructor.
    const PathGraph1& graph;
    uint64_t componentId;
    const Assembler& assembler;

    // Operations to simplify the CompressedPathGraph1.
    bool mergeLinearChains();
    bool localTransitiveReduction(uint64_t distance);
    bool detangleVertices(uint64_t detangleTolerance);
    bool detangleLinearChains(uint64_t detangleTolerance);
    bool detangleVertex(vertex_descriptor, uint64_t detangleTolerance);
    // bool detangleSuperbubbles(uint64_t minReliableLength);
    bool detangleKnots();
    void detangleIteration(
        const string& name,     // For graphviz output
        uint64_t compressedTransitiveReductionDistance,
        uint64_t detangleTolerance);
    bool removeCrossEdges(uint64_t threshold1, uint64_t threshold2);
    void detangle();



    // A Knot of degree m>1 is a special kind of tangle in the CompressedPathGraph1,
    // loosely defined as the portion of the graph "between" a pair of vertices v0 and v1.
    // More precisely it is defined as follows:
    // - v0 has in-degree m and v1 has out-degree m (usually m=2).
    // - The m parents of v0 are called the sources of the Knot.
    // - The m children  v1 are called the targets of the Knot.
    // - All forward paths starting at v0 reach v1 before reaching any vertices outside the Knot.
    // - All backward paths ending at v1 reach v1 before reaching any vertices outside the Knot.
    // - The Knot internal vertices are the vertices forward accessible from v0
    //   and backward accessible from v1.
    // - There is a 1-1 mapping between the source vertices and the target vertices
    //   such that most/all reads "move" along that mapping.
    // These conditions imply that there are no edges to/from Knot vertices
    // to/from vertices outside the Knot.
    class Knot {
    public:
        vertex_descriptor cv0;
        vertex_descriptor cv1;
        vector<vertex_descriptor> sourceVertices;   // The parents of cv0
        vector<vertex_descriptor> targetVertices;
        vector<vertex_descriptor> internalVertices;

        uint64_t degree() const
        {
            const uint64_t m = sourceVertices.size();
            SHASTA_ASSERT(targetVertices.size() == m);
            return m;
        }

        // The tangle matrix stores the number of common oriented reads
        // between each pair of source/target vertices.
        // Indexed by [i][j] where i is an index into soruceVertices
        // and j is an index into targetVertices.
        vector< vector<uint64_t> > tangleMatrix;
    };
    void findKnots(vector<Knot>&) const;



    // The id of the next vertex to be created.
    uint64_t nextVertexId = 0;

    // Return a string to uniquely identify this vertex as componentId-id.
    string vertexIdString(vertex_descriptor) const;

    uint64_t totalBaseOffset(vertex_descriptor) const;

    uint64_t totalVertexBaseOffset() const;
    uint64_t totalEdgeBaseOffset() const;

    // Sequence assembly.
    void assembleVertices(
        uint64_t threadCount0,
        uint64_t threadCount1);
    class AssembleVerticesData {
    public:
        uint64_t threadCount1;
        vector<vertex_descriptor> allVertices;
    };
    void assembleVerticesThreadFunction(uint64_t threadId);
    AssembleVerticesData assembleVerticesData;
    void assembleVertex(
        vertex_descriptor,
        uint64_t threadCount1);


    // Output.
    void writeGraphviz(
        bool labels,
        const string& fileNamePrefix) const;
    void writeGraphviz(
        const string& fileNamePrefix) const;

    void writeVerticesCsv() const;
    void writeGfa(const string& fileNamePrefix) const;
    void writeGfaAndGraphviz(const string& fileNamePrefix) const;

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
};



class shasta::mode3b::CompressedPathGraph1AEdge {
public:
    uint64_t id;
    vector<MarkerGraphEdgeId> chain;
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

    // Detangling.
    void detangle(
        uint64_t detangleThresholdLow,
        uint64_t detangleThresholdHigh,
        uint64_t pathLengthForChokePoints,
        uint64_t maxBubbleIndexDelta
        );
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
        uint64_t maxBubbleIndexDelta);
    class ChokePointChain;
    void detangleChokePointChain(
        ChokePointChain&,
        uint64_t maxBubbleIndexDelta);



    // The tangle matrix defined by two vertices.
    class TangleMatrix {
    public:
        vertex_descriptor cv0;
        vertex_descriptor cv1;

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
    };
    void computeTangleMatrix(
        vertex_descriptor,
        vertex_descriptor,
        TangleMatrix&
        ) const;
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
    };
    void findChokePointChains(uint64_t pathLengthForChokePoints, vector<ChokePointChain>&) const;
    // void analyzeChokePoints() const;
    void findVerticesBetweenChokePoints(
        vertex_descriptor,
        vertex_descriptor,
        vector<vertex_descriptor>&
    ) const;



    // Accessors.

    // Get the vertex_descriptor corresponding to a PathGraph1::vertex_descriptor,
    // adding a vertex if necessary.
    vertex_descriptor getCompressedVertex(PathGraph1::vertex_descriptor);

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


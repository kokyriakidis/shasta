#ifndef SHASTA_MODE3B_PATH_FILLER2_HPP
#define SHASTA_MODE3B_PATH_FILLER2_HPP

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3b {
        class PathFiller2Vertex;
        class PathFiller2Edge;
        class PathFiller2;
        using PathFiller2BaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathFiller2Vertex,
            PathFiller2Edge
            >;
        class PathFiller2DisplayOptions;
    }
    class Assembler;
};



class shasta::mode3b::PathFiller2DisplayOptions {
public:

    // If this is not open, no output takes place.
    ostream& html;

    bool showGraph = false;
    bool showVertices = false;
    bool showVertexLabels = false;
    bool showEdgeLabels = false;
    bool showDebugInformation = false;

    PathFiller2DisplayOptions(ostream& html) : html(html) {}
};



class shasta::mode3b::PathFiller2Vertex {
public:

    // The corresponding vertex in the global marker graph.
    // There can be more than one PathFiller2Vertex corresponding to
    // a given marker graph vertex. The replicaIndex is used to
    // identify them.
    MarkerGraphVertexId vertexId;
    uint64_t replicaIndex = 0;
    string stringId() const;

    // The ordinals in this vertex for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one ordinal in a vertex.
    // This has the same size and is indexed in the same way
    // as the orientedReadInfos vector.
    vector< vector<int64_t> > ordinals;
    uint64_t coverage() const;

    // Required by approximateTopologicalSort and only used for display.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    PathFiller2Vertex(MarkerGraphVertexId, uint64_t orientedReadCount);
    PathFiller2Vertex() {}
};



class shasta::mode3b::PathFiller2Edge {
public:

    // The MarkerIntervals in this edge for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one MarkerInterval in an edge.
    // This vector has the same size as the PathFiller::orientedReadInfos vector
    // and is indexed in the same way.
    // This only stores the ordinals for the MarkerInterval as
    // the OrientedReadId is stored in the orientedReadInfos vector.
    vector< vector< pair<int64_t, int64_t> > > markerIntervals;

    uint64_t coverage() const;

    // Set by approximateTopologicalSort and only used for display.
    bool isDagEdge = false;

    // Sequence assembled using MSA of the sequences contributes by each
    // MarkerInterval.
    vector<Base> sequence;
};



class shasta::mode3b::PathFiller2 : public PathFiller2BaseClass {
public:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

    PathFiller2(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        const PathFiller2DisplayOptions&);

    void getSequence(
        bool includeFirstAndLastEdge,
        vector<Base>&) const;

private:

    // Store constructor arguments.
    const Assembler& assembler;
    MarkerGraphEdgeId edgeIdA;
    MarkerGraphEdgeId edgeIdB;
    const PathFiller2DisplayOptions& options;
    ostream& html;

    void checkAssumptions() const;



    // A class used to store an ordinal and the corresponding position
    // of a marker in an oriented read.
    // Store signed to facilitate manipulations that involve subtractions.
    class OrdinalAndPosition {
    public:
        int64_t ordinal = invalid<int64_t>;
        int64_t position  = invalid<int64_t>;
        bool isValid() const
        {
            return ordinal != invalid<int64_t>;
        }
    };



    // Information about the portion of an oriented read used in this assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        OrientedReadInfo(OrientedReadId orientedReadId) :
            orientedReadId(orientedReadId)
            {}

        // The position of the source and target vertex of edgeIdA in this oriented read.
        // If this oriented read does not appear in edgeIdA, these are left uninitialized.
        OrdinalAndPosition ordinalAndPositionA0;
        OrdinalAndPosition ordinalAndPositionA1;
        bool isOnA() const
        {
            return ordinalAndPositionA0.isValid();
        }

        // The position of the source and target vertex of edgeIdB in this oriented read.
        // If this oriented read does not appear in edgeIdB, this is left uninitialized.
        OrdinalAndPosition ordinalAndPositionB0;
        OrdinalAndPosition ordinalAndPositionB1;
        bool isOnB() const
        {
            return ordinalAndPositionB0.isValid();
        }

        // Order by OrientedReadId.
        bool operator<(const OrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }

        // The ordinal offset between the source vertex of edgeIdA
        // and the target vertex of edgeIdB.
        int64_t ordinalOffset() const
        {
            SHASTA_ASSERT(isOnA() and isOnB());
            return ordinalAndPositionB1.ordinal - ordinalAndPositionA0.ordinal;
        }

        // The base offset between the source vertex of edgeIdA
        // and the target vertex of edgeIdB.
        int64_t positionOffset() const
        {
            SHASTA_ASSERT(isOnA() and isOnB());
            return ordinalAndPositionB1.position - ordinalAndPositionA0.position;
        }

        // The first and last ordinals of this oriented read used for this assembly.
        int64_t firstOrdinal;
        int64_t lastOrdinal;

        // The vertices corresponding to each of the ordinals of
        // this oriented read between firstOrdinal and lastOrdinal included.
        vector<vertex_descriptor> vertices;
        vertex_descriptor& vertexAtOrdinal(uint64_t ordinal)
        {
            return vertices[ordinal - firstOrdinal];
        }
        const vertex_descriptor& vertexAtOrdinal(uint64_t ordinal) const
        {
            return vertices[ordinal - firstOrdinal];
        }
    };



    // For assembly, we use the union of the oriented reads
    // that appear in edgeIdA and edgeIdB.
    // In contrast, class PathFiller1 used the intersection of those oriented reads,
    // which resulted in low coverage and low accuracy when edgeIdA and edgeIdB are
    // very distant from each other, and there are not many oriented reads that
    // cover both of them.
    // OrientedReadInfos are stored sorted by OrientedReadId.
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    void writeOrientedReads() const;

    // The index of an OrientedReadId is its index in the orientedReadInfos vector.
    uint64_t getOrientedReadIndex(OrientedReadId) const;

    // Estimated offset in bases between the source vertex of edgeIdA
    // and the target vertex of edgeIdB.
    // The estimate is done using the oriented reads that appear
    // both in edgeIdA and edgeIdB.
    int64_t estimatedA0B1Offset;
    void estimateOffset();

    // Vertex creation.
    void createVertices(double estimatedOffsetRatio);
    void createVerticesHelper(
        uint64_t i,
        int64_t ordinal,
        std::map<MarkerGraphVertexId, vertex_descriptor>& vertexMap);
    void removeLowCoverageVertices(uint64_t minVertexCoverage);
    void removeVertex(vertex_descriptor);
    void splitVertices(int64_t maxBaseSkip);

    void createEdges();
    void removeAllEdges();

    // Remove all vertices in non-trivial strongly connected components.
    void removeStrongComponents();

    // The assembly path, including edgeIdA at the beginning
    // and edgeIdB at the end.
    vector<edge_descriptor> assemblyPath;
    void findAssemblyPath();

    // Assemble edges using MSA.
    // Sequences stored in the marker graph are not used.
    void assembleAssemblyPathEdges();
    void assembleEdge(edge_descriptor);

    // Output.
    void approximateTopologicalSort();
    void writeGraph() const;
    void writeGraph(const string& title);
    void writeGraphviz(ostream&) const;
    void vertexCoverageHistogram() const;
    void edgeCoverageHistogram() const;

};

#endif

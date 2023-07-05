#ifndef SHASTA_MODE3B_PATH_FILLER1_HPP
#define SHASTA_MODE3B_PATH_FILLER1_HPP

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include <map>
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"


namespace shasta {
    namespace mode3b {
        class PathFiller1Vertex;
        class PathFiller1Edge;
        class PathFiller1;
        using PathFiller1BaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathFiller1Vertex,
            PathFiller1Edge
            >;
    }

    class Assembler;
    class Base;
    class MarkerInterval;
};




class shasta::mode3b::PathFiller1Vertex {
public:
    // The corresponding vertex in the global marker graph.
    // There can be more than one PathFillerVertex corresponding to
    // a given marker graph vertex. The replicaIndex is used to
    // identify them.
    MarkerGraphVertexId vertexId;
    uint64_t replicaIndex = 0;
    string stringId() const;

    // The ordinals in this vertex for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one ordinal in a vertex.
    // This has the same size as the orientedReadInfos vector above
    // and is indexed in the same way.
    vector< vector<uint32_t> > ordinals;
    uint64_t coverage() const;

    // Required by approximateTopologicalSort and only used for display.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    PathFiller1Vertex(MarkerGraphVertexId, uint64_t orientedReadCount);
    PathFiller1Vertex() {}
};



class shasta::mode3b::PathFiller1Edge {
public:

    // The MarkerIntervals in this edge for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one MarkerInterval in an edge.
    // This vector has the same size as the PathFiller::orientedReadInfos vector
    // and is indexed in the same way.
    // This only stores the ordinals for the MarkerInterval as
    // the OrientedReadId is stored in the orientedReadInfos vector.
    vector< vector< pair<uint32_t, uint32_t> > > markerIntervals;

    // Sequence assembled using MSA of the sequences contributes by each
    // MarkerInterval.
    vector<Base> sequence;

    uint64_t coverage() const
    {
        uint64_t c = 0;
        for(const auto& v: markerIntervals) {
            c += v.size();
        }
        return c;
    }

    // Set by approximateTopologicalSort and only used for display.
    bool isDagEdge = false;
};



class shasta::mode3b::PathFiller1: public PathFiller1BaseClass {
public:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

    PathFiller1(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        ostream& html,
        bool showGraph,
        bool showVertices,
        bool showVertexLabels,
        bool showEdgeLabels,
        bool showDebugInformation);

    // Get the assembled sequence.
    // The sequences of edgeIdA and edgeIdB are only included if
    // includePrimary is true.
    void getSequence(vector<Base>&, bool includePrimary) const;

    uint64_t coverage() const
    {
        return orientedReadInfos.size();
    }

private:

    // The path secondary edges. This excludes the primary edges edgeIdA and edgeIdB.
    vector<MarkerGraphEdgeId> secondaryEdges;

    // Store constructor arguments.
    const Assembler& assembler;
    MarkerGraphEdgeId edgeIdA;
    MarkerGraphEdgeId edgeIdB;

    void checkAssumptions() const;



    // The OrientedReadIds we will be using for this assembly.
    // These are the ones that are common between edgeIdA and edgeIdB,
    // and that visit edgeIdB after visiting edgeIdA.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The ordinals of this oriented read on the start/end edges
        // edgeIdA and edgeIdB
        uint32_t ordinalA0;
        uint32_t ordinalA1;
        uint32_t ordinalB0;
        uint32_t ordinalB1;

        // The corresponding positions in the oriented read.
        uint32_t positionA0;
        uint32_t positionA1;
        uint32_t positionB0;
        uint32_t positionB1;

        // Ordinal and base offset between ordinal A0 and B1.
        uint32_t ordinalOffset;
        uint32_t baseOffset;

        // The vertices corresponding to each of the ordinals of
        // this oriented read between ordinalA0 and ordinalB1 included.
        vector<vertex_descriptor> vertices;
        vertex_descriptor getVertex(uint32_t ordinal) const;
    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    void writeOrientedReads(ostream& html) const;

    // Average ordinal offsets and base offsets.
    uint32_t ordinalOffset;
    uint32_t baseOffset;

    // Get the base sequence corresponding to an edge.
    // The edge can be virtual.
    void getEdgeSequence(edge_descriptor, vector<Base>&) const;
    void getEdgeSequence(MarkerGraphEdgeId, vector<Base>&) const;

    void createGraph(
        uint64_t maxBaseSkip,
        uint64_t minVertexCoverage,
        ostream& html,
        bool showGraph,
        bool showVertices,
        bool showVertexLabels,
        bool showEdgeLabels,
        bool showDebugInformation);
    void createVertices();
    void removeVertex(vertex_descriptor);
    void removeLowCoverageVertices(uint64_t minVertexCoverage);
    void splitVertices(uint64_t maxBaseSkip);
    void createEdges();
    void removeAllEdges();

    vertex_descriptor getEntrance() const;
    vertex_descriptor getExit() const;

    // Approximate topological sort is only used for better and
    // faster display in Graphviz. It sets rank and color in vertices
    // and isDagEdge in edges.
    void approximateTopologicalSort();

    void removeStrongComponents();

    // Linearize the graph by keeping only vertices on the dominator tree path
    // from edgeIdA to edgeIdB.
    void linearize();

    // Assemble edges using MSA.
    // Sequences stored in the marker graph are not used.
    bool edgesWereAssembled = false;
    void assembleEdges();
    void assembleEdge(edge_descriptor);

    // Get the edge sequence from the marker graph, for a regular edge,
    // or from the edge itself, for a virtual edge.
    span<const Base> getEdgeSequence(edge_descriptor) const;

    // The assembly path, including edgeIdA at the beginning
    // and edgeIdB at the end.
    vector<edge_descriptor> assemblyPath;
    void findAssemblyPath();

    // Output.
    void writeGraph(
        ostream& html,
        bool showVertices,
        bool showVertexLabels,
        bool showEdgeLabels) const;
    void writeVerticesCsv() const;
    void writeGraphviz(
        ostream&,
        bool showVertices,
        bool showVertexLabels,
        bool showEdgeLabels) const;

    void writeSequence(ostream& html) const;
    void writeSequenceFasta(ostream& html) const;
    void writeAssemblyDetails(ostream& csv) const;
};

#endif

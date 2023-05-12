#ifndef SHASTA_MODE3B_PATH_FILLER_HPP
#define SHASTA_MODE3B_PATH_FILLER_HPP

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
        class PathFillerVertex;
        class PathFillerEdge;
        class PathFiller;
        using PathFillerBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathFillerVertex,
            PathFillerEdge
            >;
    }

    class Assembler;
    class Base;
    class MarkerInterval;
};




class shasta::mode3b::PathFillerVertex {
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

    // Set if this vertex is part of a non-trivial strongly connected component.
    uint64_t strongComponentId = invalid<uint64_t>;
    bool isStrongComponentVertex() const {
        return strongComponentId != invalid<uint64_t>;
    }

    // Required by approximateTopologicalSort and only used for display.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    PathFillerVertex(MarkerGraphVertexId, uint64_t orientedReadCount);
    PathFillerVertex() {}
};


// If an edge was created using multiple sequence alignment,
// it is called a virtual edge.
// A virtual edge only stores a sequence of MarkerGraphEdgeIds.
class shasta::mode3b::PathFillerEdge {
public:
    // The corresponding edge in the global marker graph.
    // There can be more than one PathFillerEdge corresponding to
    // a given marker graph edge.
    MarkerGraphEdgeId edgeId = invalid<MarkerGraphEdgeId>;

    // The MarkerIntervals in this edges for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one MarkerInterval in an edge.
    // This vector has the same size as the PathFiller::orientedReadInfos vector
    // and is indexed in the same way.
    // This only stores the ordinals for the MarkerInterval as
    // the OrientedReadId is stored in the orientedReadInfos vector.
    // For non-virtual edges, for each marker interval,
    // the second ordinal equals the first ordinal plus one,
    // so we could store just the first one.
    vector< vector< pair<uint32_t, uint32_t> > > markerIntervals;

    // Fields only stored for a virtual edge.
    // The sequence of MarkerGraphEdges is obtained via MSA.
    vector<MarkerGraphEdgeId> edgeIds;
    uint64_t virtualCoverage = invalid<uint64_t>;
    bool isVirtual = false;

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



class shasta::mode3b::PathFiller: public PathFillerBaseClass {
public:

    PathFiller(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        ostream& html);

    // The path secondary edges. This excludes the primary edges edgeIdA and edgeIdB.
    vector<MarkerGraphEdgeId> secondaryEdges;


private:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

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

    void createGraph(uint64_t maxBaseSkip);
    void createVertices();
    void splitVertices(uint64_t maxBaseSkip);
    void createEdges();

    // Approximate topological sort is only used for better and
    // faster display in Graphviz. It sets rank and color in vertices
    // and isDagEdge in edges.
    void approximateTopologicalSort();

#if 0
    // Find the edge v0->v1 that contains the specified MarkerInterval
    // for the i-th oriented read.
    // If no such edge, the second field of the return value is false.
    pair<edge_descriptor, bool> findEdge(
        vertex_descriptor v0,
        vertex_descriptor v1,
        uint64_t i,
        uint32_t ordinal0,
        uint32_t ordinal1) const;
#endif

    // Strongly connected component.
    // Only store the non-trivial ones.
    // A non-trivial strong component has at least one internal edge.
    // This means that it either has more than one vertex,
    // or it consists of a single vertex with a self-edge.
    class StrongComponent {
    public:
        vector<vertex_descriptor> vertices;
        StrongComponent(const vector<vertex_descriptor>& vertices);
    };
    void computeStrongComponents();
    vector<StrongComponent> strongComponents;

    // Remove all strong component vertices and the edges that involve them,
    // and create virtual edges to replace the removed edges.
    void createVirtualEdges();

    // Assemble a virtual edge.
    // This uses MSA so compute optimal MarkerGraphEdgeIds.
    void assembleVirtualEdge(edge_descriptor);

    bool assemble(ostream& html);

    // Get the assembled sequence.
    // The sequences of edgeIdA and edgeIdB are only included if
    // includePrimary is true.
    void getSequence(vector<Base>&, bool includePrimary) const;

    // Output.
    void writeGraph(ostream& html) const;
    void writeVerticesCsv() const;
    void writeGraphviz(ostream&) const;
    void writeSequence(ostream& html) const;
};

#endif

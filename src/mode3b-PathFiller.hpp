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
    // a given marker graph vertex.
    MarkerGraphVertexId vertexId;
    uint64_t replicaIndex = 0;
    string stringId() const;

    // The ordinals in this vertex for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one ordinal in a vertex.
    // This has the same size as the orientedReadInfos vector above
    // and is indexed in the same way.
    vector< vector<uint32_t> > ordinals;

    // Return true if any oriented reads have more than one ordinal.
    bool hasDuplicateOrientedReads() const;

    // Required by approximateTopologicalSort.
    uint64_t color;
    uint64_t rank;

    PathFillerVertex(MarkerGraphVertexId, uint64_t orientedReadCount);
    PathFillerVertex() {}
};



class shasta::mode3b::PathFillerEdge {
public:
    // The corresponding edge in the global marker graph.
    // There can be more than one PathFillerEdge corresponding to
    // a given marker graph edge.
    MarkerGraphEdgeId edgeId = invalid<MarkerGraphEdgeId>;

    // The MarkerIntervals in this edges for each of the oriented reads.
    // If cycles are present within this local marker graph,
    // an oriented read can have more than one MarkerInterval in an edge.
    // This vector has the same size as the orientedReadInfos vector above
    // and is indexed in the same way.
    // This only stores the ordinals for the MarkerInterval as
    // the OrientedReadId is stored in the orientedReadInfos vector.
    vector< vector< pair<uint32_t, uint32_t> > > markerIntervals;

    // Return the total number of marker intervals.
    uint64_t coverage() const;

    // Return true if any oriented reads have more than one ordinal.
    bool hasDuplicateOrientedReads() const;

    // Set by approximateTopologicalSort.
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


    span<const shasta::Base> edgeSequence(edge_descriptor) const;
    span<const shasta::Base> edgeSequence(MarkerGraphEdgeId) const;

    void createGraph(uint64_t maxBaseSkip);
    void createVertices();
    void splitVertices(uint64_t maxBaseSkip);
    void createEdges();
    void approximateTopologicalSort();

    void writeGraph(ostream& html) const;
    void writeVerticesCsv() const;
    void writeGraphviz(ostream&) const;

    bool fillPathGreedy(ostream& html);

    // Get the sequence.
    // The sequences of edgeIdA and edgeIdB are only included if
    // includePrimary is true.
    void getSequence(vector<Base>&, bool includePrimary) const;

    void writeSequence(ostream& html) const;
};

#endif

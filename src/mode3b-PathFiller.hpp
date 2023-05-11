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
    uint64_t coverage() const;

    // Return true if any oriented reads have more than one ordinal.
    bool hasDuplicateOrientedReads() const;

    // Required by approximateTopologicalSort.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    // Set if this vertex is part of a non-trivial strongly connected component.
    uint64_t strongComponentId = invalid<uint64_t>;
    bool isStrongComponentVertex() const {
        return strongComponentId != invalid<uint64_t>;
    }
    bool isStrongComponentEntrance = false;
    bool isStrongComponentExit = false;

    // The indexes of the VirtualEdges that have this vertex as their source.
    vector<uint64_t> virtualEdgeIndexes;

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

    // Find the edge v0->v1 that contains the specified MarkerInterval
    // for the i-th oriented read.
    edge_descriptor findEdge(
        vertex_descriptor v0,
        vertex_descriptor v1,
        uint64_t i,
        uint32_t ordinal0,
        uint32_t ordinal1) const;

    class StrongComponent {
    public:
        vector<vertex_descriptor> vertices;
        vector<vertex_descriptor> entrances;
        vector<vertex_descriptor> exits;
        StrongComponent(const vector<vertex_descriptor>& vertices);
    };
    void computeStrongComponents();
    vector<StrongComponent> strongComponents;
    void writeStrongComponents(ostream&) const;
    void writeStrongComponent(ostream&, uint64_t strongComponentId) const;
    void writeStrongComponentGraphviz(ostream&, uint64_t strongComponentId) const;

    // This returns true if the specified edge is internal to a
    // strongly connected component.
    bool isStrongComponentEdge(edge_descriptor) const;



    // Virtual edges are used as a replacement for edges
    // internal to strongly connected components,
    // and as a way to avoid cycles during assembly.
    // A virtual edge joins an entrance of a strongly connected
    // component with an exit of the same strongly connected component.
    // It stores a sequence of marker graph edges obtained by MSA
    // of thew oriented reads that enter/exit the strongly connected component
    // at the entrance/exit joined by the virtual edge.
    // Virtual edges are not stored in the boost graph. They are stored here instead.
    class VirtualEdge {
    public:
        // The strong component, entrance, and exit that generated this virtual edge.
        // The entrance is the source vertex of the virtual edge.
        // The exit is the target vertex of the virtual edge.
        uint64_t strongComponentId;
        vertex_descriptor entrance;
        vertex_descriptor exit;

        uint64_t coverage;

        // The sequence of marker graph edges computed by the MSA
        // for this virtual edge.
        vector<MarkerGraphEdgeId> markerGraphEdges;
    };
    vector<VirtualEdge> virtualEdges;
    void createVirtualEdges();
    void createVirtualEdges(uint64_t strongComponentId);



    void writeGraph(ostream& html) const;
    void writeVerticesCsv() const;
    void writeGraphviz(ostream&, bool showVirtualEdges) const;

    bool fillPathGreedy(ostream& html);

    // Get the sequence.
    // The sequences of edgeIdA and edgeIdB are only included if
    // includePrimary is true.
    void getSequence(vector<Base>&, bool includePrimary) const;

    void writeSequence(ostream& html) const;
};

#endif

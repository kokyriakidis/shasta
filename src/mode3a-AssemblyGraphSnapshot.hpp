#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP

// Class mode3a::AssemblyGraphSnapshot creates a persistent
// snapshot of the mode3a::Assembly graph.
// It is stored using MemoryMapped classes.

// Shasta.
#include "mode3a-AssemblyGraph.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include "string.hpp"

namespace shasta {
    namespace mode3a {
        class AssemblyGraphSnapshot;

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;

        class PackedMarkerGraph;
    }

    class Base;
}



class shasta::mode3a::AssemblyGraphSnapshot :
    public MappedMemoryOwner {
public:

    // This creates a snapshot of the AssemblyGraph in the current state.
    AssemblyGraphSnapshot(
        const AssemblyGraph&,
        const string& name,
        const MappedMemoryOwner& mappedMemoryOwner);

    // This accesses an existing snapshot.
    AssemblyGraphSnapshot(
        const string& name,
        const MappedMemoryOwner& mappedMemoryOwner,
        const PackedMarkerGraph& packedMarkerGraph);

    const string name;

    // The MarkerGraph and PackedMarkerGraph that this AssemblyGraphSnapshot refers to.
    const PackedMarkerGraph& packedMarkerGraph;

    // Class Vertex stores the segmentId and segmentReplicaIndex from AssemblyGraphVertex.
    // For the purpose of the snapshot, the index into this vector
    // is used as a vertex id.
    class Vertex {
    public:
        uint64_t segmentId;
        uint64_t segmentReplicaIndex;
        Vertex(const AssemblyGraphVertex&);
        Vertex() {}
        string stringId() const
        {
            string s = to_string(segmentId);
            if(segmentReplicaIndex != 0) {
                s += "." + to_string(segmentReplicaIndex);
            }
            return s;
        }

        // A vertex can appear at most once as a primary vertex
        // in a TangledAssemblyPath.
        // Store that information here.
        uint64_t tangledPathId = invalid<uint64_t>;
        uint64_t positionInTangledPath = invalid<uint64_t>;

        // If this vertex is used somewhere in the PackedAssemblyGraph,
        // store that information here.
        uint64_t packedAssemblyGraphVertexId = invalid<uint64_t>;
        uint64_t positionInPackedAssemblyGraph = invalid<uint64_t>;
    };
    MemoryMapped::Vector<Vertex> vertexVector;  // Can't call it vertices due to boost graph macros.

    // Get the stringId for a given vertexId, or "None" if vertexId is invalid<uint64_t>.
    string vertexStringId(uint64_t vertexId) const;

    // Class Edge stores vertex ids (indexes into the vertices vector).
    class Edge {
    public:
        uint64_t vertexId0;
        uint64_t vertexId1;
    };
    MemoryMapped::Vector<Edge> edgeVector;  // Can't call it edges due to boost graph macros.

    // The journey of each oriented read, stored as a sequence of vertex ids
    // (that is, indexes into the segments vector above).
    // Indexed by OrientedRead::getValue().
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> journeys;

    // The journey entries for each vertex.
    MemoryMapped::VectorOfVectors<JourneyEntry, uint64_t> vertexJourneyEntries;

    // The transitions of each edge (link) don't need to be stored.
    // They can be computed from the vertexJourneyhEntries.
    class Transition {
    public:
        OrientedReadId orientedReadId;
        uint64_t position;  // The transition if between position and position+1
    };
    void getEdgeTransitions(uint64_t edgeId, vector<Transition>&) const;
    uint64_t getEdgeCoverage(uint64_t edgeId) const;

    // Find out if the segments of an edge are adjacent in the marker graph.
    bool segmentsAreAdjacent(uint64_t edgeId) const;

    // Connectivity.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> edgesBySource;
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> edgesByTarget;

    // Return true if the specified edge edgeId = v0->v1
    // satisfies both of the following two conditions:
    // - It has maximum coverage among the out-edges of v0.
    // - It has maximum coverage among the in-edges of v1.
    bool isMaximalCoverageEdge(uint64_t edgeId) const;

    // A data structure that allows to get a vertexId (index in vertexVector)
    // given a segmentId and segmentReplicaIndex.
    // Indexed by segmentId.
    // vertexTable[segmentId][segmentReplicaIndex] contains the vertexId
    // with the given segmentId and segmentReplicaIndex, or invalid<uint64_t>
    // if no such vertex.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> vertexTable;
    void createVertexTable(const PackedMarkerGraph&);

    // Use the vertex table to get a vertex id given segment id
    // and replica index. If not found, returns invalid<uint64_t>
    // and also fills in a message string.
    uint64_t getVertexId(
        uint64_t segmentId,
        uint64_t segmentReplicaIndex,
        string& message
    ) const;
    uint64_t getVertexId(
        uint64_t segmentId,
        uint64_t segmentReplicaIndex
    ) const;



    // Analyze the simple tangle at a given vertex vertexId1 by
    // "following the reads" one step forward or backward.
    // On return:
    // - The previousVertices and nextVertices vectors
    //   have size equal to the number of journey entries for vertexId1,
    //   and contain the vertices vertexId0 and vertexId2 each of
    //   the oriented reads visits immediately before/after visiting vertexId1.
    //   This can be invalid<uint64_t> if vertexId1 is the
    //   first or last vertex of a journey.
    // - inCoverage[vertexId0] gives the number of oriented reads that
    //   visit vertexId0 immediately before visiting vertexId1.
    // - outCoverage[vertexId2] gives the number of oriented reads that
    //   visit vertexId0 immediately after visiting vertexId1.
    // - tangle_matrix[make_pair(vertexId0, vertexId2)] is
    //   the number of oriented reads that visit vertexId0 immediately before vertexId1, and
    //   visit vertexId2 immediately after  vertexId1.
    // The three maps can include entries for which vertexId0 and/or vertexId2
    // are invalid<uint64_t>, if vertexId1 is the first or last
    // vertex of one or more oriented read journeys.
    void analyzeSimpleTangleAtVertex(
        uint64_t vertexId1,
        vector<uint64_t>& inVertices,
        vector<uint64_t>& outVertices,
        std::map<uint64_t, uint64_t>& inCoverage,
        std::map<uint64_t, uint64_t>& outCoverage,
        std::map< pair<uint64_t, uint64_t>, uint64_t>& tangleMatrix) const;

    // Compute the Jaccard similarity of the oriented read composition of two vertices.
    // Duplicate oriented reads in the path entries for the vertices are ignored.
    // This also computes the number of oriented reads in the union and intersection
    // of the two read compositions, as well as vectors containing the deduplicated
    // oriented reads for each of the two vertices.
    double jaccard(
        uint64_t vertexId0,
        uint64_t vertexId1,
        vector<OrientedReadId>& orientedReadIds0,   // Deduplicated
        vector<OrientedReadId>& orientedReadIds1,   // Deduplicated
        vector<OrientedReadId>& unionOrientedReads,
        vector<OrientedReadId>& intersectionOrientedReads) const;

    // Get deduplicated oriented read ids for a vertex.
    void getDeduplicatedOrientedReads(
        uint64_t vertexId,
        vector<OrientedReadId>&
        ) const;

    void write() const;
    void writeGfa(uint64_t minLinkCoverage) const;
    void writeJourneys() const;
    void writeJourneyEntries() const;
    void writeTransitions() const;
    void writeLinkTransitionsHtml(uint64_t linkId, ostream& html) const;

    // Compute the consensus sequence for a link.
    // Also compute the number of bases at the end
    // of the left segment and at the beginning of the
    // right segment that are overridden by the link
    // and should be ignored during assembly.
    void assembleLink(
        uint64_t linkId,
        ostream& html,
        vector<Base>& consensusSequence,
        uint64_t& leftOverride,
        uint64_t& rightOverride
        ) const;

    // Compute the MSA for a link using spoa.
    // Takes as input a vector of (sequence, frequency) containing
    // the MSA sequences for the oriented reads of the link
    // and the number of times each was found.
    static void linkMsaUsingSpoa(
        const vector< pair<vector<Base>, uint64_t> >& msaSequences,
        ostream& html,
        vector<Base>& consensusSequence
        );



    // The TangledAssemblyPaths from the AssemblyGraph.
    class TangledAssemblyPathEntry {
    public:
        uint64_t vertexId;
        bool isPrimary;
    };
    MemoryMapped::VectorOfVectors<TangledAssemblyPathEntry, uint64_t>
        tangledAssemblyPaths;



    // Create an assembly path and assemble its sequence.
    class AssemblyPath {
    public:
        void clear();

        class Vertex {
        public:
            uint64_t id;

            // The length of the complete assembled sequence for this segment.
            uint64_t sequenceLength;

            // The portion of the segment sequence that is used to assemble this path.
            uint64_t sequenceBegin;
            uint64_t sequenceEnd;

            // The position of this portion of sequence in the assembled path.
            uint64_t pathPosition;
        };
        vector<Vertex> vertices;

        class Edge {
        public:
            uint64_t edgeId;

            // Information obtained from assembleLink.
            uint64_t leftOverride;
            uint64_t rightOverride;
            vector<Base> sequence;

            // The portion of the link sequence that is used to assemble this path.
            uint64_t sequenceBegin;
            uint64_t sequenceEnd;

            // The position of this portion of sequence in the assembled path.
            uint64_t pathPosition;
        };
        vector<Edge> edges;
        void getAssembledSequence(
            const PackedMarkerGraph&,
            const AssemblyGraphSnapshot&,
            vector<Base>&) const;
    };
    void createAssemblyPath(
        const vector<uint64_t>& vertexIds,
        AssemblyPath&) const;
    void writeAssemblyPath(
        const AssemblyPath&,
        ostream& html
    ) const;
};

#endif

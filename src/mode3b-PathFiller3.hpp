#ifndef SHASTA_MODE3B_PATH_FILLER3_HPP
#define SHASTA_MODE3B_PATH_FILLER3_HPP

// PathFiller3 assembles the sequence between two primary marker graph edges.
// It uses a local marker graph.

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3b {
        class PathFiller3Vertex;
        class PathFiller3Edge;
        class PathFiller3;
        using PathFiller3BaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathFiller3Vertex,
            PathFiller3Edge
            >;
        class PathFiller3DisplayOptions;
        class PathFiller3MarkerIndexes;
    }
    class Assembler;
};



class shasta::mode3b::PathFiller3DisplayOptions {
public:

    // If this is not open, no output takes place.
    ostream& html;

    bool showGraph = false;
    bool showOrientedReads = false;
    bool showMarkers = false;
    bool showVertices = false;
    bool showVertexLabels = false;
    bool showEdgeLabels = false;
    bool showAssemblyDetails = false;
    bool showDebugInformation = false;

    PathFiller3DisplayOptions(ostream& html) : html(html) {}
};



// A way to identify a marker in PathFiller3, besides its id.
class shasta::mode3b::PathFiller3MarkerIndexes {
public:
    uint64_t i; // Index in orientedReadInfos
    uint64_t j; // Index in OrientedReadInfo::markerInfos;
};



class shasta::mode3b::PathFiller3Vertex {
public:
    uint64_t disjointSetId;
    bool isAccessibleA = false;
    bool isAccessibleB = false;
};



class shasta::mode3b::PathFiller3Edge {
public:

    // Each marker interval is identified by the two markers.
    vector< pair<PathFiller3MarkerIndexes, PathFiller3MarkerIndexes> > markerIntervals;

    uint64_t coverage() const
    {
        return markerIntervals.size();
    }

    // Consensus of the sequences contributes by each marker interval.
    vector<Base> consensusSequence;
    vector<uint64_t> consensusCoverage;
};



class shasta::mode3b::PathFiller3 : public PathFiller3BaseClass {
public:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

    PathFiller3(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        uint64_t minVertexCoverage, // 0 = automatic
        const PathFiller3DisplayOptions&);

    // Get the sequence between edgeIdA and edgeIdB.
    // This does not include the sequences of edgeIdA and edgeIdB themselves.
    void getSecondarySequence(
        vector<Base>&) const;

    // Get the complete sequence, including the sequences of edgeIdA and edgeIdB.
    void getCompleteSequence(
        vector<Base>&) const;

private:

    // Store constructor arguments.
    const Assembler& assembler;
    MarkerGraphEdgeId edgeIdA;
    MarkerGraphEdgeId edgeIdB;
    const PathFiller3DisplayOptions& options;
    ostream& html;

    MarkerGraphVertexId vertexIdA;  // The target vertex of marker graph edge edgeIdA.
    MarkerGraphVertexId vertexIdB;  // The target vertex of marker graph edge edgeIdA.

    void checkAssumptions() const;



    // A class used to store information about a marker of
    // an oriented read used in this assembly.
    // The ordinal and position are stored signed to facilitate manipulations
    // that involve subtractions.
    class MarkerInfo {
    public:
        int64_t ordinal;
        int64_t position;
        KmerId kmerId;

        // An id for this marker, global to the PathFiller3.
        // This is the index of this marker in the disjoint sets data structure.
        uint64_t id;

        // The id of the disjoint set this MarkerInfo belongs to.
        uint64_t disjointSetId;

    };



    // Information about the portion of an oriented read used in this assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        OrientedReadInfo(OrientedReadId orientedReadId) :
            orientedReadId(orientedReadId)
            {}

        // The ordinal of vertexIdA in this oriented read.
        // Only initialized for oriented reads that appear in edgeIdA.
        int64_t ordinalA = invalid<int64_t>;
        bool isOnA() const
        {
            return ordinalA != invalid<int64_t>;
        }

        // The ordinal of vertexIdB in this oriented read.
        // Only initialized for oriented reads that appear in edgeIdB.
        int64_t ordinalB = invalid<int64_t>;
        bool isOnB() const
        {
            return ordinalB != invalid<int64_t>;
        }

        // Note we are assuming that each oriented read appears once on edgeIdA, edgeIdB,
        // and their source and target vertices.

        // Order OrientedReadInfos by OrientedReadId.
        bool operator<(const OrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }


        // The ordinal offset between vertexIdA and vertexIdB.
        int64_t ordinalOffset() const
        {
            SHASTA_ASSERT(isOnA() and isOnB());
            return ordinalB - ordinalA;
        }

        // Information about the markers of this read we will use in this assembly.
        // The first one is at ordinal firstOrdinal.
        // The last one is a ordinal lastOrdinal.
        vector<MarkerInfo> markerInfos;

        // The first and last ordinals of this oriented read used for this assembly.
        // For reads on edgeIdA, firstOrdinal equals ordinalA.
        // For reads on edgeIdB, lastOrdinal  equals ordinalB.
        int64_t firstOrdinal()
        {
            SHASTA_ASSERT(not markerInfos.empty());
            return markerInfos.front().ordinal;
        }
        int64_t lastOrdinal()
        {
            SHASTA_ASSERT(not markerInfos.empty());
            return markerInfos.back().ordinal;
        }

    };

    // Get the base position of a marker in an oriented read
    // given the ordinal.
    int64_t basePosition(OrientedReadId, int64_t ordinal) const;

    // For assembly, we use the union of the oriented reads
    // that appear in edgeIdA and edgeIdB, and that have positive ordinal offset.
    // OrientedReadInfos are stored sorted by OrientedReadId.
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    void writeOrientedReads() const;

    // Estimated offset in bases between vertexIdA and vertexIdB.
    // The estimate is done using the oriented reads that appear
    // both in edgeIdA and edgeIdB.
    int64_t estimatedABOffset;
    void estimateOffset();

    // Fill in the markerInfos vector of each read.
    void gatherMarkers(double estimatedOffsetRatio);
    void writeMarkers();

    // Add the marker at given ordinal to the i-th oriented read.
    void addMarkerInfo(uint64_t i, int64_t ordinal);

    // Compute alignments and use them to create the disjoint set data structure,
    // from which the marker graph will be created.
    void alignAndDisjointSets(
        uint64_t matchScore,
        uint64_t mismatchScore,
        uint64_t gapScore,
        uint64_t maxSkipBases
        );

    // This stores the markers in each disjoint set.
    // Each marker is stored as pair(i, j)
    // where i is the index of the OrientedReadInfo in orientedReadInfos
    // and j is the index of the MarkerInfo in orientedReadInfo.markerInfos.
    // Keyed by the disjoint set id (the same also stored in each marker).
    std::map<uint64_t, vector<PathFiller3MarkerIndexes> > disjointSetsMap;

    vector<uint64_t> disjointSetsSizeHistogram;

    // Create vertices. Each disjoint set with at least minVertexCoverage markers
    // generates a vertex.
    void createVertices(
        uint64_t minVertexCoverage,
        double vertexSamplingRate);  // Only used if minVertexCoverage is 0;
    void removeVertex(vertex_descriptor);

    // The disjoint sets corresponding to vertexIdA and vertexIdB.
    // Those will always generate a vertex regardless of coverage.
    uint64_t disjointSetIdA = invalid<uint64_t>;
    uint64_t disjointSetIdB = invalid<uint64_t>;

    // Map that gives the vertex descriptor corresponding to a disjoint set id, if any.
    std::map<uint64_t, vertex_descriptor> vertexMap;

    // Create edges by following the reads.
    void createEdges();
    void removeAllEdges();

    void removeStrongComponents();

    // Remove vertices that are not accessible from vertexIdA
    // or from which vertexIdB is not accessible.
    // Returns the number of vertices that were removed.
    void removeInaccessibleVertices();

    // The assembly path, beginning at vertexIdA and ending at vertexIdB.
    // This means that the sequences of edgeIdA and edgeIdB are not included.
    vector<edge_descriptor> assemblyPath;
    void findAssemblyPath();
    bool assembleAssemblyPathEdges(uint64_t maxMsaLength);
    bool assembleEdge(uint64_t maxMsaLength, edge_descriptor);

    // Graphviz output.
    void writeGraph() const;
    void writeGraph(const string& title);
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    void writeCoverageCharacterToHtml(uint64_t coverage) const;
};

#endif

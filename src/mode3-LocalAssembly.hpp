#pragma once

// LocalAssembly assembles the sequence between two primary marker graph edges.
// It uses a local marker graph.

// Shasta.
#include "AssemblerOptions.hpp"
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
    namespace mode3 {
        class LocalAssemblyVertex;
        class LocalAssemblyEdge;
        class LocalAssembly;
        using LocalAssemblyBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            LocalAssemblyVertex,
            LocalAssemblyEdge
            >;
        class LocalAssemblyDisplayOptions;
        class LocalAssemblyMarkerIndexes;

        class Anchors;
        using AnchorId = uint64_t;
    }

    class CompressedMarker;
    class Reads;

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
};



class shasta::mode3::LocalAssemblyDisplayOptions {
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

    LocalAssemblyDisplayOptions(ostream& html) : html(html) {}
};



// A way to identify a marker in LocalAssembly, besides its id.
class shasta::mode3::LocalAssemblyMarkerIndexes {
public:
    uint64_t i; // Index in orientedReadInfos
    uint64_t j; // Index in OrientedReadInfo::markerInfos;
};



class shasta::mode3::LocalAssemblyVertex {
public:
    uint64_t disjointSetId;
    bool isAccessibleA = false;
    bool isAccessibleB = false;
};



class shasta::mode3::LocalAssemblyEdge {
public:

    // Each marker interval is identified by the two markers.
    vector< pair<LocalAssemblyMarkerIndexes, LocalAssemblyMarkerIndexes> > markerIntervals;

    uint64_t coverage() const
    {
        return markerIntervals.size();
    }

    // Consensus of the sequences contributes by each marker interval.
    vector<Base> consensusSequence;
    vector<uint64_t> consensusCoverage;
};



class shasta::mode3::LocalAssembly : public LocalAssemblyBaseClass {
public:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

    // The oriented reads common between anchorIdA and anchorIdB are always
    // used for assembly. The oriented reads that appear only
    // on anchorIdA or anchorIdB are used for assembly under control
    // of useA and useB.
    // So, if useA and useB are both true (the default), the assembly uses the
    // union of the oriented reads on anchorIdA and anchorIdB.
    // If they are both false, the assembly uses the
    // intersection of the oriented reads on anchorIdA and anchorIdB.
    // If useA is true and useB is false, the assembly uses the
    // oriented reads on anchorIdA, regardless of whether they appear on anchorIdB.
    // If useA is false and useB is true, the assembly uses the
    // oriented reads on anchorIdB, regardless of whether they appear on anchorIdA.
    LocalAssembly(
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const Anchors& anchors,
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        uint64_t minVertexCoverage, // 0 = automatic
        const LocalAssemblyDisplayOptions&,
        const Mode3AssemblyOptions::LocalAssemblyOptions&,
        bool useA = true,
        bool useB = true);

    // Get the sequence between edgeIdA and edgeIdB.
    // This does not include the sequences of edgeIdA and edgeIdB themselves.
    void getSecondarySequence(
        vector<Base>&) const;

    // Get the complete sequence, including the sequences of edgeIdA and edgeIdB.
    void getCompleteSequence(
        vector<Base>&) const;

private:

    // Store constructor arguments.
    uint64_t k;
    uint64_t kHalf;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const Anchors& anchors;
    AnchorId anchorIdA;
    AnchorId anchorIdB;
    const LocalAssemblyDisplayOptions& options;
    ostream& html;

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

        // An id for this marker, global to the LocalAssembly.
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
    void gatherOrientedReads(bool useA, bool useB);
    void writeOrientedReads() const;
    void writeOrientedReadsSequences() const;

    // Estimated offset in bases between vertexIdA and vertexIdB.
    // The estimate is done using the oriented reads that appear
    // both in edgeIdA and edgeIdB.
    // If the offset cannot be estimated because there are no
    // common oriented reads between egeIdA and edgeIdB,
    // it is set to invalid<int64_t>.
    // In that case, or if the offset is negative,
    // the assembly fails, which results in empty secondary sequence.
    int64_t estimatedABOffset;
    void estimateOffset();

    // Fill in the markerInfos vector of each read.
    void gatherMarkers(double estimatedOffsetRatio);
    void writeMarkers();

    // Add the marker at given ordinal to the i-th oriented read.
    void addMarkerInfo(uint64_t i, int64_t ordinal);

    // Compute alignments and use them to create the disjoint set data structure,
    // from which the marker graph will be created.
    // maxDrift is the maximum tolerated length drift of each read.
    // Used to compute the band for banded alignments.
    void alignAndDisjointSets(
        uint64_t matchScore,
        uint64_t mismatchScore,
        uint64_t gapScore,
        uint64_t maxSkipBases,
        double maxDrift,
        uint64_t minHalfBand,
        double minScoreRatio
        );

    // This stores the markers in each disjoint set.
    // Each marker is stored as pair(i, j)
    // where i is the index of the OrientedReadInfo in orientedReadInfos
    // and j is the index of the MarkerInfo in orientedReadInfo.markerInfos.
    // Keyed by the disjoint set id (the same also stored in each marker).
    std::map<uint64_t, vector<LocalAssemblyMarkerIndexes> > disjointSetsMap;

    vector<uint64_t> disjointSetsSizeHistogram;

    // Create vertices. Each disjoint set with at least minVertexCoverage markers
    // generates a vertex.
    // If minVertexCoverage is 0, a suitable value is computed.
    // This returns the value of minVertexCoverage actually used.
    uint64_t createVertices(
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

    // Remove strongly connected components.
    // Returns the number of vertices removed.
    uint64_t removeStrongComponents();

    // Remove vertices that are not accessible from vertexIdA
    // or from which vertexIdB is not accessible.
    // Returns the number of vertices that were removed.
    uint64_t removeInaccessibleVertices();

    // Possible courses of action when a long MSA is encountered.
    enum class LongMsaPolicy {
        throwException,
        assembleAtLowCoverage
    };

    // The assembly path, beginning at vertexIdA and ending at vertexIdB.
    // This means that the sequences of edgeIdA and edgeIdB are not included.
    vector<edge_descriptor> assemblyPath;
    void findAssemblyPath();
    void assembleAssemblyPathEdges(uint64_t maxMsaLength, LongMsaPolicy);
    void assembleEdge(
        uint64_t maxMsaLength,
        LongMsaPolicy,
        edge_descriptor);

    // Graphviz output.
    void writeGraph() const;
    void writeGraph(const string& title);
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    void writeCoverageCharacterToHtml(uint64_t coverage) const;

    // Remove all vertices and edges and clear the vertexMap and assemblyPath.
    // All other data are left alone.
    void clear();
};


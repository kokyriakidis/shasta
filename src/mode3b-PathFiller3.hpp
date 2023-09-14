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
    }
    class Assembler;
};



class shasta::mode3b::PathFiller3DisplayOptions {
public:

    // If this is not open, no output takes place.
    ostream& html;

    bool showGraph = false;
    bool showVertices = false;
    bool showVertexLabels = false;
    bool showEdgeLabels = false;
    bool showAssemblyDetails = false;
    bool showDebugInformation = false;

    PathFiller3DisplayOptions(ostream& html) : html(html) {}
};



class shasta::mode3b::PathFiller3Vertex {
public:

    // The ordinals in this vertex for each of the oriented reads used in this assembly.
    // This has one entry corresponding to each entry of the
    // PathFiller3::orientedInfos vector.
    vector< vector<int64_t> > ordinals;

};



class shasta::mode3b::PathFiller3Edge {
public:

    // The MarkerIntervals in this edge for each of the oriented reads used in this assembly.
    // This has one entry corresponding to each entry of the
    // PathFiller3::orientedInfos vector.
    // Each pair contains ordinals for the source and target vertex.
    vector< vector< pair<int64_t, int64_t> > > markerIntervals;

};



class shasta::mode3b::PathFiller3 : public PathFiller3BaseClass {
public:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

    PathFiller3(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        const PathFiller3DisplayOptions&);

    // Get the sequence between edgeIdA and edgeIdB.
    // This does not include the sequences od edgeIdA and edgeIdB themselves.
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
};

#endif

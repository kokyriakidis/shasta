#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "ReadId.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <list>
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3a {

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            AssemblyGraphVertex, AssemblyGraphEdge>;

        class JourneyEntry;
        class Transition;

        class PackedMarkerGraph;
    }
}



// An entry in the journey of an oriented read.
class shasta::mode3a::JourneyEntry {
public:
    OrientedReadId orientedReadId;

    // The position in the journey for this oriented read.
    uint64_t position;
};



// Each vertex represents a segment replica.
class shasta::mode3a::AssemblyGraphVertex {
public:

    // The PackedMarkerGraph segment that this vertex corresponds to.
    uint64_t segmentId;

    // Replica index among all vertices with the same segmentId.
    uint64_t segmentReplicaIndex;

    AssemblyGraphVertex(uint64_t segmentId, uint64_t segmentReplicaIndex=0) :
        segmentId(segmentId), segmentReplicaIndex(segmentReplicaIndex) {}

    string stringId() const
    {
        string s = to_string(segmentId);
        if(segmentReplicaIndex != 0) {
            s += "." + to_string(segmentReplicaIndex);
        }
        return s;
    }

    // The journey entries that go through this vertex.
    vector<JourneyEntry> journeyEntries;

    // Partial paths for a vertex are computed
    // by following the path entries in this vertex, forward or backward,
    // until we encounter a divergence or low coverage.
    // Both partial paths are stored with the closes vertex to
    // the start vertex first.
    vector<AssemblyGraphBaseClass::vertex_descriptor> forwardPartialPath;
    vector<AssemblyGraphBaseClass::vertex_descriptor> backwardPartialPath;
};



// A transition of an oriented read from a segment to the next.
class shasta::mode3a::Transition {
public:
    uint64_t position0;
    uint64_t position1;
    OrientedReadId orientedReadId;

    bool operator<(const Transition& that) const
    {
        return
            tie(orientedReadId, position0, position1) <
            tie(that.orientedReadId, that.position0, that.position1);
    }
    bool operator==(const Transition& that) const
    {
        return
            tie(orientedReadId, position0, position1) ==
            tie(that.orientedReadId, that.position0, that.position1);
    }
};



// Each edge represents a link.
class shasta::mode3a::AssemblyGraphEdge {
public:
};



class shasta::mode3a::AssemblyGraph :
    public AssemblyGraphBaseClass,
    public MultithreadedObject<shasta::mode3a::AssemblyGraph> {
public:
    AssemblyGraph(const PackedMarkerGraph&);
    void write(const string& name) const;

    // The sequence of vertices visited by each oriented read.
    // It has one entry for each entry of the corresponding journey in the PackedMarkerGraph.
    // An entry can be null_vertex(), if the corresponding vertex was removed
    // during detangling.
    // Indexed by OrientedReadId::getValue().
    vector< vector<vertex_descriptor> > journeys;

    // Each segment in the AssemblyGraph corresponds to a segment in
    // this PackedMarkerGraph.
    const PackedMarkerGraph& packedMarkerGraph;
private:

    // The vertices created so far for each segmentId.
    // Some of these may have been deleted and set to null_vertex().
    vector< vector<vertex_descriptor> > verticesBySegment;

    void createSegmentsAndJourneys();
    void createLinks();

    // Find out if two segments are adjacent in the marker graph.
    bool segmentsAreAdjacent(edge_descriptor e) const;
    bool segmentsAreAdjacent(vertex_descriptor, vertex_descriptor) const;

    // Get the stringId for a given vertex_descriptor, or "None" if v is null_vertex().
    string vertexStringId(vertex_descriptor v) const;

    // Get the transitions for an edge.
    void getEdgeTransitions(edge_descriptor, vector<Transition>&) const;
    uint64_t edgeCoverage(edge_descriptor) const;

    void writeGfa(const string& name, uint64_t minLinkCoverage) const;
    void writeLinkCoverageHistogram(const string& name) const;
    void writeJourneys(const string& name) const;

    // Simple detangling, one vertex at a time, looking only
    // at immediate parent and children.
public:
    void simpleDetangle(
        uint64_t minLinkCoverage,
        uint64_t minTangleCoverage);
private:
    void simpleDetangle(
        vertex_descriptor,
        uint64_t minLinkCoverage,
        uint64_t minTangleCoverage);

    // Find the previous and next vertex for each JourneyEntry in a given vertex.
    // On return, adjacentVertices contains a pair of vertex descriptors for
    // each JourneyEntry in vertex v, in the same order.
    // Those vertex descriptors are the previous and next vertex visited
    // by the oriented read for that JourneyEntry, and can be null_vertex()
    // if v is at the beginning or end of the journey of an oriented read.
    void findAdjacentVertices(
        vertex_descriptor v,
        vector< pair<vertex_descriptor, vertex_descriptor> >& adjacentVertices
    ) const;



    // Compute partial paths by following oriented reads in a starting vertex.
    // Each partial path is stored in its starting vertex.
    // See AssemblyGraphVertex for more information.
public:
    void computePartialPaths(
        uint64_t threadCount,
        uint64_t segmentCoverageThreshold1,
        uint64_t segmentCoverageThreshold2,
        uint64_t minLinkCoverage
        );
    void writePartialPaths() const;
private:
    void computePartialPathsThreadFunction(uint64_t threadId);

    // Simple computation of partial paths.
    void computePartialPath1(
        vertex_descriptor,
        uint64_t minLinkCoverage,
        ostream& debugOut
        );

    // More robust version that uses dominator trees.
    void computePartialPath2(
        vertex_descriptor,
        uint64_t segmentCoverageThreshold1,
        uint64_t segmentCoverageThreshold2,
        uint64_t minLinkCoverage,
        ostream& debugOut
        );

    // Data used by computePartialPaths multithreaded code.
    class ComputePartialPathsData {
    public:
        uint64_t segmentCoverageThreshold1;
        uint64_t segmentCoverageThreshold2;
        uint64_t minLinkCoverage;
    };
    ComputePartialPathsData computePartialPathsData;



    // Analyze partial paths to create assembly paths.
public:
    void analyzePartialPaths(uint64_t threadCount);
private:
    void analyzePartialPathsThreadFunction(uint64_t threadId);
    void analyzePartialPathsComponent(
        const vector<uint64_t>& component,
        const vector< pair<vertex_descriptor, vertex_descriptor> >& componentPairs,
        vector<vertex_descriptor>& longestPath,
        ostream& graphOut);
    class AnalyzePartialPathsData {
    public:

        // Contiguous numbering of the vertices of the AssemblyGraph
        // and the partial paths graph.
        vector<vertex_descriptor> vertexTable;
        std::map<vertex_descriptor, uint64_t> vertexMap;

        // This stores indexes in the vertexTable.
        vector< vector<uint64_t> > components;

        // The bidirectional pairs in each connected component.
        vector< vector< pair<vertex_descriptor, vertex_descriptor> > > componentPairs;

        // The computed longest paths for each connected component.
        // These paths are incomplete as they only include primary vertices.
        // Secondary vertices will be filled later by computeTangledAssemblyPaths.
        vector< vector<vertex_descriptor> > longestPaths;
    };
    AnalyzePartialPathsData analyzePartialPathsData;



    // Class TangledAssemblyPath describes a tangled
    // assembly path ready to be "ripped" from the AssemblyGraph.
    // Sequence assembly only happens after the "ripping" operation takes place.
    // At that time the assembly path becomes a linear sequence of vertices in
    // the AssemblyGraph.
    // Each TangledAssemblyPath corresponds to one of the longest paths
    // computed by analyzePartialPaths and stored in
    // analyzePartialPathsData.longestPaths.
    class TangledAssemblyPath {
    public:

        // The primary vertices are the ones that are unique to this path
        // (that is, they don't belong to any other assembly path).
        // They are computed by analyzePartialPaths.
        // All of the oriented reads in the primary vertices participate in the ripping.
        vector<vertex_descriptor> primaryVertices;

        // For each pair of consecutive primary vertices,
        // we store a vector of secondary vertices, in the order in which
        // they appear in the path. For each of the secondary vertices
        // we also store indices of the vertex journey entries that will
        // participate in the ripping.
        // The secondaryVertices vector has size one less than the
        // primaryVertices vector. Each entry corresponds to the interval
        // between two consecutive primaryVertices.
        class SecondaryVertexInfo {
        public:
            vertex_descriptor v;
            vector<uint64_t> journeyEntryIndexes;
        };
        vector< vector<SecondaryVertexInfo> > secondaryVertices;
    };
    vector<TangledAssemblyPath> tangledAssemblyPaths;
public:
    void computeTangledAssemblyPaths(uint64_t threadCount);
private:
    void computeTangledAssemblyPathsThreadFunction(uint64_t threadId);
    void computeTangledAssemblyPath(
        const vector<vertex_descriptor>& primaryVertices,
        TangledAssemblyPath&,
        ostream& debugOut
        );
    bool computeSecondaryVertices(
        vertex_descriptor v0,
        vertex_descriptor v1,
        vector<TangledAssemblyPath::SecondaryVertexInfo>&,
        ostream& debugOut);
};

#endif


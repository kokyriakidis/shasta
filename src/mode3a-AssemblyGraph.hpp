#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <list>
#include "memory.hpp"
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {

    class Base;
    namespace mode3a {

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            AssemblyGraphVertex, AssemblyGraphEdge>;
        class AssemblyGraphVertexPredicate;
        class AssemblyGraphEdgePredicate;

        class JourneyEntry;
        class Transition;
        class PathInformation;

        class PackedMarkerGraph;

        // Compute the MSA for a link using spoa.
        // Takes as input a vector of (sequence, frequency) containing
        // the MSA sequences for the oriented reads of the link
        // and the number of times each was found.
        void linkMsaUsingSpoa(
            const vector< pair<vector<Base>, uint64_t> >& msaSequences,
            ostream& html,
            vector<Base>& consensusSequence
            );
    }

}



// An entry in the journey of an oriented read.
class shasta::mode3a::JourneyEntry {
public:
    OrientedReadId orientedReadId;

    // The position in the journey for this oriented read.
    uint64_t position;
};



// Class used to store information about the paths that a vertex
// appears in.
class shasta::mode3a::PathInformation {
public:

    // Each vertex can apear at most once as a primary vertex in tangled paths.
    class PrimaryInfo {
    public:
        uint64_t pathId = invalid<uint64_t>;
        uint64_t positionInPath = invalid<uint64_t>;
    };
    PrimaryInfo primaryInfo;

    // Return the number of times this vertex appears as primary.
    uint64_t primaryCount() const
    {
        if(primaryInfo.pathId == invalid<uint64_t>) {
            return 0;
        } else {
            return 1;
        }
    }

    // Information about appearances of this vertex as a secondary vertex
    // in paths.
    class SecondaryInfo {
    public:
        uint64_t pathId;
        uint64_t positionInPath;
        uint64_t positionInLeg;
    };
    vector<SecondaryInfo> secondaryInfos;

    void clear()
    {
        primaryInfo.pathId = invalid<uint64_t>;
        primaryInfo.positionInPath = invalid<uint64_t>;
        secondaryInfos.clear();
    }
};



// Each vertex represents a segment replica.
class shasta::mode3a::AssemblyGraphVertex {
public:

    // The PackedMarkerGraph segment that this vertex corresponds to.
    uint64_t segmentId;

    // Replica index among all vertices with the same segmentId.
    uint64_t segmentReplicaIndex;

    // If this vertex is used somewhere in the PackedAssemblyGraph,
    // store that information here.
    uint64_t packedAssemblyGraphVertexId = invalid<uint64_t>;
    uint64_t positionInPackedAssemblyGraph = invalid<uint64_t>;

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
    // Because of the way they are constructed, the journey entries of a vertex
    // are ordered by OrientedReadId and then by position.
    vector<JourneyEntry> journeyEntries;

    // The deduplicated oriented reads in this vertex.
    vector<OrientedReadId> orientedReadIds;
    void computeOrientedReadIds();

    // This is used when we want to work with a filtered_graph.
    bool isActive = true;

    // Partial paths for a vertex are computed
    // by following the path entries in this vertex, forward or backward,
    // until we encounter a divergence or low coverage.
    // Both partial paths are stored with the closes vertex to
    // the start vertex first.
    vector<AssemblyGraphBaseClass::vertex_descriptor> forwardPartialPath;
    vector<AssemblyGraphBaseClass::vertex_descriptor> backwardPartialPath;

    // Information about the paths that this vertex appears in.
    PathInformation pathInformation;
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
    // This is used when we want to work with a filtered_graph.
    bool isActive = true;
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

public:
    // Get the stringId for a given vertex_descriptor, or "None" if v is null_vertex().
    string vertexStringId(vertex_descriptor v) const;

    // Get the transitions for an edge.
    void getEdgeTransitions(edge_descriptor, vector<Transition>&) const;
    uint64_t edgeCoverage(edge_descriptor) const;
private:

    void writeGfa(const string& name, uint64_t minLinkCoverage) const;
    void writeLinkCoverageHistogram(const string& name) const;
    void writeJourneys(const string& name) const;

    // Find descendants/ancestors of a given vertex up to a specified distance.
    // This is done by following the journeys.
    void findDescendants(
        vertex_descriptor, uint64_t distance,
        vector<vertex_descriptor>&) const;
    void findAncestors(
        vertex_descriptor, uint64_t distance,
        vector<vertex_descriptor>&) const;


    // Filename prefix for debug output.
    string debugOutputPrefix = "Mode3a-";
public:
    void setDebugOutputPrefix(const string& prefix) {
        debugOutputPrefix = prefix;
    }
private:

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


public:

    // Class AssemblyPath describes an assembly path in the AssemblyGraph.
    // An AssemblyPath consists of a sequence of primary vertices
    // interspeded with secondary vertices.
    // Successive primary vertices have similar read compositions
    // and represent the "native" sequence of the assembly paths.
    // The secondary vertices have extraneous reads
    // originating from unrelated sequence copies, and which
    // should not be used for assembly.
    // Each AssemblyPath corresponds to one of the longest paths
    // computed by analyzePartialPaths and stored in
    // analyzePartialPathsData.longestPaths.
    class AssemblyPath {
    public:

        vector<vertex_descriptor> primaryVertices;

        // For each pair of consecutive primary vertices we store a
        // SecondaryVertexInfo describing the secondary vertices in between.
        // The secondaryVerticesInfos vector has size one less than the
        // primaryVertices vector. Each entry corresponds to the interval
        // between two consecutive primaryVertices.
        class SecondaryVertexInfo {
        public:

            // The secondary vertices between these two primary vertices,
            // stored in the order in which they appear in the path.
            vector<vertex_descriptor> secondaryVertices;

            // The journey intervals for these secondary vertices.
            // This are the ones corresponding to the entries that will be
            // used to assemble this portion of the TangleAssemblyPath.
            class JourneyInterval {
            public:
                OrientedReadId orientedReadId;
                uint64_t begin;
                uint64_t end;
            };
            vector<JourneyInterval> journeyIntervals;

            // Given a vertex, find which journey entries in the vertex
            // are in one of the above journey intervals for this SecondaryVertexInfo.
            void getVertexJourneyEntries(
                const AssemblyGraphVertex&,
                vector<bool>&   // True of false for each of the journey entries in the vertex.
            ) const;

            // The fraction of graph vertices that become secondary vertices.
            double efficiency = 0.;

            bool failed = false;
        };
        vector<SecondaryVertexInfo> secondaryVerticesInfos;

        double efficiency;
        void computeEfficiency()
        {
            double sum = 0.;
            for(const SecondaryVertexInfo& secondaryVertexInfo : secondaryVerticesInfos) {
                sum += secondaryVertexInfo.efficiency;
            }
            efficiency = sum /double(secondaryVerticesInfos.size());
        }

        // Assembled sequence.
        vector<shasta::Base> assembledSequence;
    };
    vector<shared_ptr<AssemblyPath> > assemblyPaths;
    void computeAssemblyPaths(uint64_t threadCount);
private:

    // Class used to sort AssemblyPaths by decreasing efficiency.
    class OrderAssemblyPath {
    public:
        bool operator()(
            const shared_ptr<AssemblyPath>& x,
            const shared_ptr<AssemblyPath>& y
            ) const
        {
            return x->efficiency > y->efficiency;
        }

    };

    void computeAssemblyPathsThreadFunction(uint64_t threadId);

    // Compute an AssemblyPath given the primary vertices.
    void computeAssemblyPath(
        const vector<vertex_descriptor>& primaryVertices,
        AssemblyPath&,
        ostream& debugOut
        );

    // Find secondary vertices given a pair of consecutive primary vertices.
    void computeSecondaryVertices(
        vertex_descriptor v0,
        vertex_descriptor v1,
        AssemblyPath::SecondaryVertexInfo&,
        ostream& debugOut);

    // AssemblyPath debug output.
    void writeAssemblyPaths1() const;
    void writeAssemblyPaths2() const;
    void writeAssemblyPathsVertexSummary() const;
    void writeAssemblyPathsVertexInfo() const;
    void writeAssemblyPathsVertexHistogram() const;
    void writeAssemblyPathsJourneyInfo() const;
    void writeAssemblyPathsJourneyIntervals() const;
public:



    // A "flattened" version of the AssemblyPath,
    // in which the distinction between primary and secondary vertices
    // is lost, and for each vertex we store the indexes of the
    // JourneyEntries of the vertex that should be used for assembly.

    class FlattenedAssemblyPathSegment {
    public:
        vertex_descriptor v;
        vector<bool> useForAssembly;  // True of false for each of the journey entries in the vertex.

        // The length of the complete assembled sequence for this segment.
        uint64_t sequenceLength;

        // The portion of the segment sequence that is used to assemble this path.
        uint64_t sequenceBegin;
        uint64_t sequenceEnd;

        // The position of this portion of sequence in assembled sequence.
        uint64_t positionInAssembledSequence;

        // Constructor for a primary vertex.
        FlattenedAssemblyPathSegment(vertex_descriptor v, uint64_t n) :
            v(v), useForAssembly(vector<bool>(n, true)) {}

        // Constructor for a secondary vertex.
        FlattenedAssemblyPathSegment(vertex_descriptor v, const vector<bool>& useForAssembly) :
            v(v), useForAssembly(useForAssembly) {}
    };

    class FlattenedAssemblyPathLink {
    public:

        edge_descriptor e;

        // Information obtained from assembleLink.
        uint64_t leftOverride;
        uint64_t rightOverride;
        vector<shasta::Base> sequence;

        // The portion of the link sequence that is used to assemble this path.
        uint64_t sequenceBegin;
        uint64_t sequenceEnd;

        // The position of this portion of sequence in assembled sequence.
        uint64_t positionInAssembledSequence;

        FlattenedAssemblyPathLink(edge_descriptor e);
    };

    public:
    class FlattenedAssemblyPath {
    public:
        vector<FlattenedAssemblyPathSegment> segments;
        vector<FlattenedAssemblyPathLink> links;    // Size is segments.size() - 1
    };




    // Use the AssemblyPaths to assemble sequence.
    void assemble();
    void assemble(uint64_t assemblyPathId);
    void assemble(
        FlattenedAssemblyPath&,
        vector<shasta::Base>& sequence);
    void assembleLink(
        FlattenedAssemblyPathLink&,
        const FlattenedAssemblyPathSegment& previousSegment,
        const FlattenedAssemblyPathSegment& nextSegment
    ) const;

    // Find the Transitions to be used to assemble a link.
    vector<Transition> transitions;
    void getTransitionsForAssembly(
        const FlattenedAssemblyPathSegment& previousSegment,
        const FlattenedAssemblyPathSegment& nextSegment,
        vector<Transition>&) const;



    // Classes used to define detangling constructors.
    class DetangleUsingTangledAssemblyPaths {};
    class DetangleUsingTangleMatrices {};
    class DetangleUsingLocalClustering {};



    // Create a detangled AssemblyGraph using the AssemblyPaths
    // of another AssemblyGraph.
    AssemblyGraph(
        DetangleUsingTangledAssemblyPaths,
        const AssemblyGraph& oldAssemblyGraph);
private:
    void createFromAssemblyPaths(const AssemblyGraph& oldAssemblyGraph);



    // Create a detangled AssemblyGraph using tangle matrices to split vertices
    // of another AssemblyGraph.
    // Only tangle matrix entries that are at least equal to minCoverage as used.
    // As a result, the new AssemblyGraph can have missing journey entries.
    // That is, some journey entries will remain set to null_vertex().
public:
    AssemblyGraph(
        DetangleUsingTangleMatrices,
        const AssemblyGraph& oldAssemblyGraph,
        uint64_t minCoverage);
private:
    void createFromTangledMatrices(
        const AssemblyGraph& oldAssemblyGraph,
        uint64_t minCoverage);



    // Detangle by local clustering.
    // Create a detangled AssemblyGraph using clustering on another AssemblyGraph.
public:
    AssemblyGraph(
        DetangleUsingLocalClustering,
        const AssemblyGraph& oldAssemblyGraph);
private:
    void createByLocalClustering(
        const AssemblyGraph& oldAssemblyGraph);



    // Graph class used by computeSecondaryVertices.
    // The vertex stores vertexFrequency.
    // The edge stores transitionFrequency.
    using SecondaryVerticesGraphBaseClass =
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, uint64_t, uint64_t>;
    class SecondaryVerticesGraph : public SecondaryVerticesGraphBaseClass {
    public:
        SecondaryVerticesGraph(
            const AssemblyGraph& assemblyGraph,
            const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered,
            const vector<uint64_t>& vertexFrequency,
            const vector< pair<AssemblyGraph::vertex_descriptor, AssemblyGraph::vertex_descriptor> >&
                transitionsEncountered,
            const vector<uint64_t>& transitionFrequency,
            vertex_descriptor iv0,
            vertex_descriptor iv1);
        SecondaryVerticesGraph(
            const AssemblyGraph& assemblyGraph,
            const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered,
            const vector<uint64_t>& vertexFrequency,
            uint64_t n);
        bool isLinear(
            vertex_descriptor v0,
            vertex_descriptor v1
        ) const;
        string vertexStringId(vertex_descriptor v) const
        {
            return assemblyGraph.vertexStringId(verticesEncountered[v]);
        }
        void write(
            ostream&,
            const string& graphName
            ) const;

        const AssemblyGraph& assemblyGraph;
        const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered;
        const vector<uint64_t>& vertexFrequency;
        vertex_descriptor iv0;
        vertex_descriptor iv1;

        bool segmentsAreAdjacent(edge_descriptor e) const
        {
            const SecondaryVerticesGraph& graph = *this;
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const AssemblyGraph::vertex_descriptor u0 = verticesEncountered[v0];
            const AssemblyGraph::vertex_descriptor u1 = verticesEncountered[v1];
            return assemblyGraph.segmentsAreAdjacent(u0, u1);
        }

        // The path on the dominator tree from iv0 to iv1.
        vector<vertex_descriptor> dominatorTreePath;
        bool computeDominatorTreePath();

        // Compute the "best" path between iv0 and iv1.
        vector<edge_descriptor> bestPath;
        void computeBestPath(ostream& debugOut);

        // An edge predicate the retursn true for edges that are not dotted.
        class EdgePredicate {
        public:
            EdgePredicate(const SecondaryVerticesGraph* secondaryVerticesGraph = 0) :
                secondaryVerticesGraph(secondaryVerticesGraph) {}
            const SecondaryVerticesGraph* secondaryVerticesGraph;
            bool operator()(const edge_descriptor e) const
            {
                const vertex_descriptor u0 = source(e, *secondaryVerticesGraph);
                const vertex_descriptor u1 = target(e, *secondaryVerticesGraph);
                const AssemblyGraph::vertex_descriptor v0 = secondaryVerticesGraph->verticesEncountered[u0];
                const AssemblyGraph::vertex_descriptor v1 = secondaryVerticesGraph->verticesEncountered[u1];
                return secondaryVerticesGraph->assemblyGraph.segmentsAreAdjacent(v0, v1);
            }

        };

        // Handle dotted edges that "skip" a vertex.
        void handleDottedEdges1(ostream& debugOut);
        class HandleDottedEdges1Callback {
        public:
            const AssemblyGraph& assemblyGraph;
            const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered;
            const vector<uint64_t>& vertexFrequency;
            const SecondaryVerticesGraph& smallGraph;
            const SecondaryVerticesGraph& secondaryVerticesGraph;
            ostream& debugOut;

            HandleDottedEdges1Callback(
                const AssemblyGraph& assemblyGraph,
                const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered,
                const vector<uint64_t>& vertexFrequency,
                const SecondaryVerticesGraph& smallGraph,
                const SecondaryVerticesGraph& secondaryVerticesGraph,
                vector<edge_descriptor>& edgesToBeRemoved,
                ostream& debugOut) :
                assemblyGraph(assemblyGraph),
                verticesEncountered(verticesEncountered),
                vertexFrequency(vertexFrequency),
                smallGraph(smallGraph),
                secondaryVerticesGraph(secondaryVerticesGraph),
                debugOut(debugOut),
                edgesToBeRemoved(edgesToBeRemoved)
                {}

            template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
            bool operator()(
                const CorrespondenceMap1To2&,
                const CorrespondenceMap2To1&);

            // The edges that wil be removed.
            // These are edge descriptors in the secondaryVerticesGraph.
            vector<edge_descriptor>& edgesToBeRemoved;
        };
    };



    // Jaccard graph.
public:
    void computeJaccardGraph(
        uint64_t threadCount,
        double minJaccard,
        uint64_t minComponentSize,
        uint64_t minPathLength);
private:
    void findJaccardGraphCandidatePairs(uint64_t threadCount);
    void findJaccardGraphCandidatePairsThreadFunction(uint64_t threadId);
    void computeJaccardPairs(uint64_t threadCount, double minJaccard);
    void computeJaccardPairsThreadFunction(uint64_t threadId);
    class ComputeJaccardGraphData {
    public:

        using VertexPair = pair<vertex_descriptor, vertex_descriptor>;

        // The candidate pairs found by each thread.
        vector< vector<VertexPair> > threadCandidatePairs;

        // The consolidated and deduplicated candidate pairs.
        vector<VertexPair> candidatePairs;

        // The good pairs found by each thread.
        vector< vector<pair<VertexPair, double> > > threadGoodPairs;

        // The consolidated and deduplicated good pairs.
        vector< pair<VertexPair, double> > goodPairs;

        double minJaccard;
    };
    ComputeJaccardGraphData computeJaccardGraphData;

public:
    void computeVertexOrientedReadIds(uint64_t threadCount);
    void clearVertexOrientedReadIds();
private:
    void computeVertexOrientedReadIdsThreadFunction(uint64_t threadId);
    class ComputeVertexOrientedReadIdsData {
    public:
        vector<vertex_descriptor> allVertices;
    };
    ComputeVertexOrientedReadIdsData computeVertexOrientedReadIdsData;

    // Compute Jaccard similarity between two vertices.
    // This requires vertex oriented read ids to be available.
public:
    double computeJaccard(
        vertex_descriptor,
        vertex_descriptor,
        vector<OrientedReadId>& commonOrientedReadIds) const;
};



class shasta::mode3a::AssemblyGraphVertexPredicate {
public:
    AssemblyGraphVertexPredicate() : assemblyGraph(0) {}
    AssemblyGraphVertexPredicate(const AssemblyGraph& assemblyGraph) :
        assemblyGraph(&assemblyGraph) {}
    const AssemblyGraph* assemblyGraph;
    bool operator()(const AssemblyGraph::vertex_descriptor v) const
    {
        return (*assemblyGraph)[v].isActive;
    }
};



class shasta::mode3a::AssemblyGraphEdgePredicate {
public:
    AssemblyGraphEdgePredicate() : assemblyGraph(0) {}
    AssemblyGraphEdgePredicate(const AssemblyGraph& assemblyGraph) :
        assemblyGraph(&assemblyGraph) {}
    const AssemblyGraph* assemblyGraph;
    bool operator()(const AssemblyGraph::edge_descriptor e) const
    {
        return (*assemblyGraph)[e].isActive;
    }
};

#endif


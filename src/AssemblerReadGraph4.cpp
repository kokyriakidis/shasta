#include "Assembler.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
#include "orderPairs.hpp"
#include "Mode3Assembler.hpp"

using namespace shasta;

// Standard library.
#include "fstream.hpp"
#include "chrono.hpp"
#include "iterator.hpp"
#include <numeric>
#include <queue>
#include <random>
#include <stack>
#include <queue>
#include <set>

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graph_utility.hpp>


namespace shasta {
    class ReadGraph4;
    class ReadGraph4Vertex;
    class ReadGraph4Edge;

    using ReadGraph4BaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ReadGraph4Vertex,
        ReadGraph4Edge>;

    class ReadGraph4AllAlignments;
    class ReadGraph4AllAlignmentsVertex;
    class ReadGraph4AllAlignmentsEdge;

    using ReadGraph4AllAlignmentsBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ReadGraph4AllAlignmentsVertex,
        ReadGraph4AllAlignmentsEdge>;
}



class shasta::ReadGraph4Vertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

};

class shasta::ReadGraph4AllAlignmentsVertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

};


class shasta::ReadGraph4Edge {
public:
    uint64_t alignmentId;
    ReadGraph4Edge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};

class shasta::ReadGraph4AllAlignmentsEdge {
public:
    uint64_t alignmentId;
    ReadGraph4AllAlignmentsEdge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};

// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph4: public ReadGraph4BaseClass {
public:


    ReadGraph4(uint64_t n) : ReadGraph4BaseClass(n) {}
    
    uint64_t findPathWithPositiveOffset(
        OrientedReadId start,
        vector<vector<OrientedReadId>>& paths,
        vector<vector<double>>& pathsOffsets,
        vector<OrientedReadId>& currentPath,
        vector<double>& currentPathOffset,
        std::set<vertex_descriptor>& visited,
        uint64_t maxDistance,
        uint64_t currentDistance,
        MemoryMapped::Vector<AlignmentData>& alignmentData,
        ReadGraph4& readGraph);
    void findNeighborsUndirectedGraph(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideRight(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideLeft(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphBothSides(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsSkipSameComponentNodes(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachSameComponentNode(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};

// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph4AllAlignments: public ReadGraph4AllAlignmentsBaseClass {
public:

    ReadGraph4AllAlignments(uint64_t n) : ReadGraph4AllAlignmentsBaseClass(n) {}
    void computeShortPath(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        size_t maxDistance,
        vector<uint32_t>& path,
        vector<uint32_t>& distance,
        vector<OrientedReadId>& reachedVertices,
        vector<uint32_t>& parentEdges,
        MemoryMapped::Vector<AlignmentData>& alignmentData,
        ReadGraph4AllAlignments& readGraph);
    void findNeighborsUndirectedGraph(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideRight(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideLeft(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphBothSides(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsSkipSameComponentNodes(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachSameComponentNode(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachEndNode(OrientedReadId orientedReadId, vector<bool>& finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, uint64_t maxDistance, vector<OrientedReadId>& neighbors, vector<uint64_t>& neighborsAlignmentIds, ReadGraph4AllAlignments& readGraphAllAlignments);
    void findNeighborsEarlyStopWhenReachEndNode(OrientedReadId orientedReadId, vector<bool>& finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    uint64_t findAllPathsToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<vector<OrientedReadId>>& paths,
        vector<OrientedReadId>& currentPath,
        std::set<vertex_descriptor>& visited,
        uint64_t maxDistance,
        uint64_t currentDistance);
    uint64_t findSinglePathToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<vector<OrientedReadId>>& paths,
        vector<OrientedReadId>& currentPath,
        std::set<vertex_descriptor>& visited,
        uint64_t maxDistance,
        uint64_t currentDistance);
    uint64_t findShortestPathToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<OrientedReadId>& shortestPath);
    uint64_t findShortestPathToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<OrientedReadId>& shortestPath,
        uint64_t maxDistance);

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};








void ReadGraph4AllAlignments::computeShortPath(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxDistance,
    vector<uint32_t>& path,
    vector<uint32_t>& distance,
    vector<OrientedReadId>& reachedVertices,
    vector<uint32_t>& parentEdges,
    MemoryMapped::Vector<AlignmentData>& alignmentData,
    ReadGraph4AllAlignments& readGraph)
{
    const bool debug = false;
    path.clear();
    reachedVertices.clear();

    // Initialize all distances to infinite
    fill(distance.begin(), distance.end(), ReadGraph::infiniteDistance);

    // Initialize BFS
    std::queue<OrientedReadId> queuedVertices;
    queuedVertices.push(orientedReadId0);
    distance[orientedReadId0.getValue()] = 0;
    reachedVertices.push_back(orientedReadId0);

    // Do the BFS
    while(!queuedVertices.empty()) {
        const OrientedReadId vertex0 = queuedVertices.front();
        queuedVertices.pop();
        const uint32_t distance0 = distance[vertex0.getValue()];
        const uint32_t distance1 = distance0 + 1;

        if(distance1 > maxDistance) {
            continue;
        }


        BGL_FORALL_OUTEDGES(vertex0.getValue(), edge, *this, ReadGraph4AllAlignments) {
        
            vertex_descriptor targetVertex = target(edge, *this);
            const OrientedReadId vertex1 = OrientedReadId::fromValue(targetVertex);

            uint32_t alignmentId = readGraph[edge].alignmentId;
            AlignmentData alignment = alignmentData[alignmentId];

            //uint32_t alignmentId = ReadGraph4AllAlignments[edge].alignmentId;

            // Process new vertices
            if(distance[vertex1.getValue()] == ReadGraph::infiniteDistance) {
                distance[vertex1.getValue()] = distance1;
                reachedVertices.push_back(vertex1);
                parentEdges[vertex1.getValue()] = alignmentId;
                queuedVertices.push(vertex1);

                // Check if we reached the target
                if(vertex1 == orientedReadId1) {
                    // Reconstruct path
                    OrientedReadId vertex = vertex1;
                    while(vertex != orientedReadId0) {
                        const uint32_t alignmentId = parentEdges[vertex.getValue()];
                        path.push_back(alignmentId);
                        AlignmentData alignment = alignmentData[alignmentId];
                        const ReadId readId0 = alignment.readIds[0];
                        const ReadId readId1 = alignment.readIds[1];
                        const bool isSameStrand = alignment.isSameStrand;
                        SHASTA_ASSERT(readId0 < readId1);
                        const OrientedReadId A0 = OrientedReadId(readId0, 0);
                        const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);

                        if(vertex.getValue() == A0.getValue()) {
                            vertex = B0;
                        } else {
                            vertex = A0;
                        }
                    }
                    std::reverse(path.begin(), path.end());
                    return;
                }
            }
                

        }

        // for(const uint32_t edgeId: connectivity[vertex0.getValue()]) {
        //     const ReadGraphEdge& edge = edges[edgeId];
        //     if(edge.crossesStrands) {
        //         continue;
        //     }
        //     const OrientedReadId vertex1 = edge.getOther(vertex0);

        //     // Process new vertices
        //     if(distance[vertex1.getValue()] == ReadGraph::infiniteDistance) {
        //         distance[vertex1.getValue()] = distance1;
        //         reachedVertices.push_back(vertex1);
        //         parentEdges[vertex1.getValue()] = edgeId;
        //         queuedVertices.push(vertex1);

        //         // Check if we reached the target
        //         if(vertex1 == orientedReadId1) {
        //             // Reconstruct path
        //             OrientedReadId vertex = vertex1;
        //             while(vertex != orientedReadId0) {
        //                 const uint32_t edgeId = parentEdges[vertex.getValue()];
        //                 path.push_back(edgeId);
        //                 vertex = edges[edgeId].getOther(vertex);
        //             }
        //             std::reverse(path.begin(), path.end());
        //             return;
        //         }
        //     }
        // }
    }
}







void ReadGraph4::findNeighborsUndirectedGraph(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if(!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}

void ReadGraph4::findNeighborsDirectedGraphOneSideRight(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Process outgoing edges
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4::findNeighborsDirectedGraphOneSideLeft(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Process incoming edges to the current vertex
            BGL_FORALL_INEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor sourceVertex = source(edge, *this);
                if (!visitedVertices.contains(sourceVertex)) {
                    visitedVertices.insert(sourceVertex);
                    q.push(make_pair(sourceVertex, distance + 1));
                }
            }
        }
    }
}




void ReadGraph4::findNeighborsDirectedGraphBothSides(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Process outgoing edges
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }

            // Process incoming edges to the current vertex
            BGL_FORALL_INEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor sourceVertex = source(edge, *this);
                if (!visitedVertices.contains(sourceVertex)) {
                    visitedVertices.insert(sourceVertex);
                    q.push(make_pair(sourceVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighborsUndirectedGraph(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if(!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighborsDirectedGraphOneSideRight(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Process outgoing edges
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighborsDirectedGraphOneSideLeft(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Process incoming edges to the current vertex
            BGL_FORALL_INEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor sourceVertex = source(edge, *this);
                if (!visitedVertices.contains(sourceVertex)) {
                    visitedVertices.insert(sourceVertex);
                    q.push(make_pair(sourceVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighborsDirectedGraphBothSides(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Process outgoing edges
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }

            // Process incoming edges
            BGL_FORALL_INEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor sourceVertex = source(edge, *this);
                if (!visitedVertices.contains(sourceVertex)) {
                    visitedVertices.insert(sourceVertex);
                    q.push(make_pair(sourceVertex, distance + 1));
                }
            }
        }
    }
}












void ReadGraph4AllAlignments::findNeighborsEarlyStopWhenReachEndNode(
    OrientedReadId orientedReadId,
    vector<bool>& finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[OrientedReadId::fromValue(currentVertex).getValue()]) {
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[OrientedReadId::fromValue(targetVertex).getValue()]) {
                        neighbors.push_back(OrientedReadId::fromValue(targetVertex));
                        return;
                    }
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}




void ReadGraph4AllAlignments::findNeighborsEarlyStopWhenReachEndNode(
    OrientedReadId orientedReadId,
    vector<bool>& finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors,
    vector<uint64_t>& neighborsAlignmentIds,
    ReadGraph4AllAlignments& readGraphAllAlignments) 
{
    neighbors.clear();
    neighborsAlignmentIds.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            // if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[OrientedReadId::fromValue(currentVertex).getValue()]) {
            //     neighbors.push_back(OrientedReadId::fromValue(currentVertex));
            //     return;
            // }
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                uint32_t alignmentId = readGraphAllAlignments[edge].alignmentId;
                if (!visitedVertices.contains(targetVertex)) {
                    neighborsAlignmentIds.push_back(alignmentId);
                    if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[OrientedReadId::fromValue(targetVertex).getValue()]) {
                        neighbors.push_back(OrientedReadId::fromValue(targetVertex));
                        return;
                    }
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}













void ReadGraph4::findNeighborsSkipSameComponentNodes(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId != sourceComponentId) {
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
            }
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    // Only explore vertices in different components
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId != sourceComponentId) {
                        visitedVertices.insert(targetVertex);
                        q.push(make_pair(targetVertex, distance + 1));
                    }
                }
            }
        }
    }
}



void ReadGraph4AllAlignments::findNeighborsSkipSameComponentNodes(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId != sourceComponentId) {
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
            }
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    // Only explore vertices in different components
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId != sourceComponentId) {
                        visitedVertices.insert(targetVertex);
                        q.push(make_pair(targetVertex, distance + 1));
                    }
                }
            }
        }
    }
}









void ReadGraph4::findNeighborsEarlyStopWhenReachSameComponentNode(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId == sourceComponentId) {
                // If we find a vertex in the same component, add it to the neighbors and stop
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId == sourceComponentId) {
                        // If we find a vertex in the same component, add it to the neighbors and stop
                        neighbors.push_back(OrientedReadId::fromValue(targetVertex));
                        return;
                    }
                    // Only explore vertices in different components
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighborsEarlyStopWhenReachSameComponentNode(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId == sourceComponentId) {
                // If we find a vertex in the same component, add it to the neighbors and stop
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId == sourceComponentId) {
                        // If we find a vertex in the same component, add it to the neighbors and stop
                        neighbors.push_back(OrientedReadId::fromValue(targetVertex));
                        return;
                    }
                    // Only explore vertices in different components
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}








// Function to perform DFS to find a path with a positive offset
uint64_t ReadGraph4::findPathWithPositiveOffset(
    OrientedReadId start,
    vector<vector<OrientedReadId>>& paths,
    vector<vector<double>>& pathsOffsets,
    vector<OrientedReadId>& currentPath,
    vector<double>& currentPathOffset,
    std::set<vertex_descriptor>& visited,
    uint64_t maxDistance,
    uint64_t currentDistance,
    MemoryMapped::Vector<AlignmentData>& alignmentData,
    ReadGraph4& readGraph)
{

    if (currentPathOffset.empty()) {
        currentPathOffset.push_back(0.);
    }

    // Base case 1: We found a read path that gave us a positive offset
    // If the current path offset is positive, we have found a read path that gave us a positive offset
    // We add the current read to the path and mark it as visited.
    // We also add the path and the path offset to the paths and pathsOffsets vectors.
    // We return 1 to indicate that we found a read path that gave us a positive offset.
    if(currentPathOffset.back() > 0.) {
        if (currentPath.size() > 0) {
            currentPath.push_back(start);
            visited.insert(start.getValue());
            paths.push_back(currentPath);
            pathsOffsets.push_back(currentPathOffset); 
        }
        return 1;
    }
    
    // Base case 2: if current distance is greater than or equal to maxDistance, return
    // If we reach the max distance, we add the current path and the path offset to the paths and pathsOffsets vectors.
    // We return 0 to indicate that we reached the max distance without founding a read path that gave us a positive offset.
    if (currentDistance >= maxDistance) {
        if (currentPath.size() > 0) {
            paths.push_back(currentPath);
            pathsOffsets.push_back(currentPathOffset);
        }
        return 0;
    }


    // Add current vertex to path and mark as visited
    vertex_descriptor currentVertex = start.getValue();
    currentPath.push_back(start);
    visited.insert(start.getValue());

    // Process outgoing edges from the current vertex
    BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
        vertex_descriptor targetVertex = target(edge, *this);
        if (!visited.contains(targetVertex)) {

            uint32_t alignmentId = readGraph[edge].alignmentId;

            // Get the alignment from this alignmentID
            AlignmentData alignment = alignmentData[alignmentId];

            // double offset = alignment.info.offsetAtCenter();

            double offset = 0.;

            if (alignment.info.offsetAtCenter() > 0) {
                offset = alignment.info.offsetAtCenter();
            } else if (alignment.info.offsetAtCenter() < 0) {
                offset = - alignment.info.offsetAtCenter();

            }

            // Get the last value of currentOffset
            double lastOffset = currentPathOffset.back();

            double newPathOffset = lastOffset + offset;

            currentPathOffset.push_back(newPathOffset);
            
            uint64_t result = findPathWithPositiveOffset(OrientedReadId::fromValue(targetVertex), paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

            if(result == 1) {
                // cout << "Finished. Target Oriented ID: " << OrientedReadId::fromValue(targetVertex) << " Target vertex ID: " << targetVertex << " Source Oriented ID: " << OrientedReadId::fromValue(currentVertex) << " Source vertex ID: " << currentVertex << " New Offset: " << offset << " Old Offset: " << lastOffset << " Final Offset: " << newPathOffset << endl;
                // We found a read path that gave us a positive offset.
                // Propagate success up the call stack
                return 1; 
            }

            if(result == 0) {
                // cout << "maxDistance Exceeded. Target Oriented ID: " << OrientedReadId::fromValue(targetVertex) << " Target vertex ID: " << targetVertex << " Source Oriented ID: " << OrientedReadId::fromValue(currentVertex) << " Source vertex ID: " << currentVertex << " New Offset: " << offset << " Old Offset: " << lastOffset << " Final Offset: " << newPathOffset << endl;
                // If we reach the max distance, we remove the last offset from the current path offset.
                currentPathOffset.pop_back();
            }
        }
    }

    // Process incoming edges
    BGL_FORALL_INEDGES(currentVertex, edge, *this, ReadGraph4) {
        vertex_descriptor sourceVertex = source(edge, *this);
        if (!visited.contains(sourceVertex)) {
         
            uint32_t alignmentId = readGraph[edge].alignmentId;

            // get the alignment form the alignmentID
            AlignmentData alignment = alignmentData[alignmentId];

            // double offset = alignment.info.offsetAtCenter();

            double offset = 0.;

            if (alignment.info.offsetAtCenter() > 0) {
                offset = - alignment.info.offsetAtCenter();
            } else if (alignment.info.offsetAtCenter() < 0) {
                offset = alignment.info.offsetAtCenter();

            }

            // Get the last value of currentOffset
            double lastOffset = currentPathOffset.back();

            double newPathOffset = lastOffset + offset;

            currentPathOffset.push_back(newPathOffset);
            
            uint64_t result = findPathWithPositiveOffset(OrientedReadId::fromValue(sourceVertex), paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

            if(result == 1) {
                // cout << "Finished! Source Oriented ID: " << OrientedReadId::fromValue(sourceVertex) << " Source vertex ID: " << sourceVertex << " Target Oriented ID: " << OrientedReadId::fromValue(currentVertex) << " Target vertex ID: " << currentVertex << " New Offset: " << offset << " Old Offset: " << lastOffset << " Final Offset: " << newPathOffset << endl;
                // We found a read path that gave us a positive offset.
                // Propagate success up the call stack
                return 1; 
            }

            if(result == 0) {
                // cout << "maxDistance Exceeded. Source Oriented ID: " << OrientedReadId::fromValue(sourceVertex) << " Source vertex ID: " << sourceVertex << " Target Oriented ID: " << OrientedReadId::fromValue(currentVertex) << " Target vertex ID: " << currentVertex << " New Offset: " << offset << " Old Offset: " << lastOffset << " Final Offset: " << newPathOffset << endl;
                // If we reach the max distance, we remove the last offset from the current path offset.
                currentPathOffset.pop_back();
            }
        }
    }

    // Propagate failure up the call stack
    return 0;
}















uint64_t ReadGraph4AllAlignments::findShortestPathToNode(
    OrientedReadId start,
    OrientedReadId endNode,
    vector<OrientedReadId>& shortestPath,
    uint64_t maxDistance) 
{
    shortestPath.clear();
    
    // Regular queue for BFS (no need for priority queue)
    std::queue<pair<OrientedReadId, uint64_t>> q;  // Pair of (node, distance)
    
    // Keep track of visited nodes and parents for path reconstruction
    std::unordered_map<vertex_descriptor, vertex_descriptor> parent;
    std::set<vertex_descriptor> visited;
    
    // Start BFS with distance 0
    q.push({start, 0});
    visited.insert(start.getValue());
    
    //const uint64_t maxDistance = 10;  // Add maximum distance limit
    bool foundPath = false;
    
    while (!q.empty() && !foundPath) {
        OrientedReadId current = q.front().first;
        uint64_t currentDistance = q.front().second;
        q.pop();
        
        // If we reached the end node, we're done
        if (current == endNode) {
            foundPath = true;
            break;
        }
        
        // Stop exploring this path if we've exceeded maxDistance
        if (currentDistance >= maxDistance) {
            continue;
        }
        
        // Explore neighbors
        BGL_FORALL_OUTEDGES(current.getValue(), edge, *this, ReadGraph4AllAlignments) {
            vertex_descriptor nextVertex = target(edge, *this);
            if (visited.find(nextVertex) == visited.end()) {
                OrientedReadId next = OrientedReadId::fromValue(nextVertex);
                visited.insert(nextVertex);
                parent[nextVertex] = current.getValue();
                q.push({next, currentDistance + 1});
            }
        }
    }
    
    // If we found a path, reconstruct it
    if (foundPath) {
        // Reconstruct path from end to start
        OrientedReadId current = endNode;
        while (current != start) {
            shortestPath.push_back(current);
            current = OrientedReadId::fromValue(parent[current.getValue()]);
        }
        shortestPath.push_back(start);
        
        // Reverse to get path from start to end
        reverse(shortestPath.begin(), shortestPath.end());
        return 1;
    }
    
    // No path found
    return 0;
}



uint64_t ReadGraph4AllAlignments::findShortestPathToNode(
    OrientedReadId start,
    OrientedReadId endNode,
    vector<OrientedReadId>& shortestPath) 
{
    shortestPath.clear();
    
    // Regular queue for BFS (no need for priority queue)
    std::queue<OrientedReadId> q;
    
    // Keep track of visited nodes and parents for path reconstruction
    std::unordered_map<vertex_descriptor, vertex_descriptor> parent;
    std::set<vertex_descriptor> visited;
    
    // Start BFS
    q.push(start);
    visited.insert(start.getValue());
    
    bool foundPath = false;
    while (!q.empty() && !foundPath) {
        OrientedReadId current = q.front();
        q.pop();
        
        // If we reached the end node, we're done
        if (current == endNode) {
            foundPath = true;
            break;
        }
        
        // Explore neighbors
        BGL_FORALL_OUTEDGES(current.getValue(), edge, *this, ReadGraph4AllAlignments) {
            vertex_descriptor nextVertex = target(edge, *this);
            if (visited.find(nextVertex) == visited.end()) {
                OrientedReadId next = OrientedReadId::fromValue(nextVertex);
                visited.insert(nextVertex);
                parent[nextVertex] = current.getValue();
                q.push(next);
            }
        }
    }
    
    // If we found a path, reconstruct it
    if (foundPath) {
        // Reconstruct path from end to start
        OrientedReadId current = endNode;
        while (current != start) {
            shortestPath.push_back(current);
            current = OrientedReadId::fromValue(parent[current.getValue()]);
        }
        shortestPath.push_back(start);
        
        // Reverse to get path from start to end
        reverse(shortestPath.begin(), shortestPath.end());
        return 1;
    }
    
    // No path found
    return 0;
}



// Function to perform DFS from a start node to an endNode or up to maxDistance
uint64_t ReadGraph4AllAlignments::findSinglePathToNode(
    OrientedReadId start,
    OrientedReadId endNode,
    vector<vector<OrientedReadId>>& paths,
    vector<OrientedReadId>& currentPath,
    std::set<vertex_descriptor>& visited,
    uint64_t maxDistance,
    uint64_t currentDistance)
{
    // Base case: if current distance is greater than or equal to maxDistance, return
    if (currentDistance >= maxDistance) {
        return 0;
    }

    // Add current vertex to path and mark as visited
    vertex_descriptor currentVertex = start.getValue();
    currentPath.push_back(start);
    visited.insert(start.getValue());

    // If we reached the end node, add the path to the paths vector
    if (start == endNode) {
        paths.push_back(currentPath);
        return 1;
    }

    // Process outgoing edges from the current vertex
    BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
        vertex_descriptor targetVertex = target(edge, *this);
        if (!visited.contains(targetVertex)) {
            uint64_t result = findAllPathsToNode(OrientedReadId::fromValue(targetVertex), endNode, paths, currentPath, visited, maxDistance, currentDistance + 1);
            if (result == 1) {
                return 1; // Propagate success up the call stack
            }
        }
    }

    // Remove current vertex from path and visited set to backtrack
    visited.erase(start.getValue());
    currentPath.pop_back();

    return 0;
}







// Function to perform DFS to find all paths from start to end node
uint64_t ReadGraph4AllAlignments::findAllPathsToNode(
    OrientedReadId start,
    OrientedReadId endNode,
    vector<vector<OrientedReadId>>& paths,
    vector<OrientedReadId>& currentPath,
    std::set<vertex_descriptor>& visited,
    uint64_t maxDistance,
    uint64_t currentDistance)
{
    // Base case: if current distance is greater than or equal to maxDistance, return
    if (currentDistance >= maxDistance) {
        return 0;
    }

    // Add current vertex to path and mark as visited
    vertex_descriptor currentVertex = start.getValue();
    currentPath.push_back(start);
    visited.insert(start.getValue());

    // If we reached the end node, add the path to the paths vector
    // but continue searching for other paths
    uint64_t pathsFound = 0;
    if (start == endNode) {
        paths.push_back(currentPath);
        pathsFound = 1;
    }

    // Process outgoing edges from the current vertex
    BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
        vertex_descriptor targetVertex = target(edge, *this);
        if (!visited.contains(targetVertex)) {
            pathsFound += findAllPathsToNode(
                OrientedReadId::fromValue(targetVertex), 
                endNode, 
                paths, 
                currentPath, 
                visited, 
                maxDistance, 
                currentDistance + 1
            );
        }
    }

    // Remove current vertex from path and visited set to backtrack
    visited.erase(start.getValue());
    currentPath.pop_back();

    return pathsFound;
}























//class AlignmentStats{public: double errorRateRle; uint32_t alignedRange; uint32_t rightUnaligned; uint32_t leftUnaligned; uint32_t alignmentId;};


void Assembler::createReadGraph4withStrandSeparation(
    uint32_t maxAlignmentCount)
{
    cout << timestamp << "createReadGraph4 with strand separation begins" << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);


    const double minQRle60 = 60.;
    const double maxErrorRateRle60 = std::pow(10.0, - minQRle60 / 10.0);

    const double minQRle50 = 50.;
    const double maxErrorRateRle50 = std::pow(10.0, - minQRle50 / 10.0);

    const double minQRle40 = 40.;
    const double maxErrorRateRle40 = std::pow(10.0, - minQRle40 / 10.0);

    const double minQRle30 = 30.;
    const double maxErrorRateRle30 = std::pow(10.0, - minQRle30 / 10.0);



    const double WThreshold16 = 1e-16;
    const double logWThreshold16 = log(WThreshold16);

    const double WThreshold14 = 1e-14;
    const double logWThreshold14 = log(WThreshold14);

    const double WThreshold12 = 1e-12;
    const double logWThreshold12 = log(WThreshold12);

    const double WThreshold11 = 1e-11;
    const double logWThreshold11 = log(WThreshold11);

    const double WThreshold10 = 1e-10;
    const double logWThreshold10 = log(WThreshold10);

    const double WThreshold9 = 1e-9;
    const double logWThreshold9 = log(WThreshold9);

    const double WThreshold8 = 1e-8;
    const double logWThreshold8 = log(WThreshold8);

    const double WThreshold7 = 1e-7;
    const double logWThreshold7 = log(WThreshold7);

    const double WThreshold65 = 5e-6;
    const double logWThreshold65 = log(WThreshold65);

    const double WThreshold6 = 1e-6;
    const double logWThreshold6 = log(WThreshold6);

    const double WThreshold55 = 5e-5;
    const double logWThreshold55 = log(WThreshold55);

    const double WThreshold5 = 1e-5;
    const double logWThreshold5 = log(WThreshold5);

    const double WThreshold4 = 1e-4;
    const double logWThreshold4 = log(WThreshold4);

    const double WThreshold3 = 1e-3;
    const double logWThreshold3 = log(WThreshold3);

    const double WThreshold2 = 1e-2;
    const double logWThreshold2 = log(WThreshold2);


    const double WThresholdForBreaks = 1e+5;
    const double logWThresholdForBreaks = log(WThresholdForBreaks);


    const double WThresholdForBreaks2 = 1e+100000002;
    const double logWThresholdForBreaks2 = log(WThresholdForBreaks2);



    //
    // 1. Order alignments in order of increasing Q. 
    //


    // Gather in alignmentTable[alignmentID, Q]
    // alignments in order of increasing Q.
    // Q(n) = (1 + δ/2ε)^n * e-δL
    // ε = 1e-4, δ = 5e-4
    // logQ(n) = αn - δL


    vector< pair<uint64_t, double> > alignmentTable;
    vector< pair<uint64_t, double> > alignmentTableNotPassFilter;
    vector< pair<uint64_t, double> > QAlignments;
    const double epsilon = 1e-4;
    const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // Get stats about the reads
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    
    // Keep track of which readIds were used in alignments
    vector<bool> readUsed(readCount, false);

    // Flag alignments to be kept for break detection.
    vector<bool> keepAlignmentsForBreaks(alignmentCount, false);


    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if((alignmentId % 10000) == 0) {
            cout << timestamp << alignmentId << "/" << alignmentCount << endl;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // The alignment is stored as an alignment between readId0 on strand 0
        // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
        // The reverse complement alignment also exists, but is not stored explicitly.
        const ReadId readId0 = thisAlignmentData.readIds[0];
        const ReadId readId1 = thisAlignmentData.readIds[1];
        const bool isSameStrand = thisAlignmentData.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
        const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

        // Store this pair of edges in our edgeTable.
        const uint32_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
        const uint32_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = thisAlignmentData.info.mismatchCountRle;
        const double errorRateRle = thisAlignmentData.info.errorRateRle;
        const uint64_t nRLE = errorRateRle * 2 * L;
        const double markerCount = thisAlignmentData.info.markerCount;

        // logQ(n) = αn - δL
        const double logQ = alpha * double(nRLE) - delta * L;

        // if(errorRateRle<= maxErrorRateRle) {
        //     QAlignments.push_back(make_pair(alignmentId, errorRateRle));
        //     alignmentTable.push_back(make_pair(alignmentId, logQ));
        // }

        // if (logQ <= logQThreshold and errorRateRle <= maxErrorRateRle) {

        // if (logQ <= logQThreshold and errorRateRle <= maxErrorRateRle) {
        //     alignmentTable.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThresholdForBreaks){
        //     alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
        //     keepAlignmentsForBreaks[alignmentId] = true;
        // }


        // cout << "logQ: " << logQ << " L: " << L << " range0: " << range0 << " range1: " << range1 << " n: " << n << " nRLE: " << nRLE << " errorRateRle: " << errorRateRle<< endl;

        //if (logQ <= logWThreshold8) {
        if (logQ <= logWThreshold8) {
            alignmentTable.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logWThresholdForBreaks){
            alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
            keepAlignmentsForBreaks[alignmentId] = true;
        }

    }

    

    sort(alignmentTable.begin(), alignmentTable.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTableNotPassFilter.begin(), alignmentTableNotPassFilter.end(), OrderPairsBySecondOnly<uint64_t, double>());
    cout << "The alignmentTable has " << alignmentTable.size() << " entries." << endl; // 123863
    cout << "The alignmentTableNotPassFilter has " << alignmentTableNotPassFilter.size() << " entries." << endl;

    





    ///
    // 2. Process alignments in order of increasing Q. 
    //
    // i.   Start with no edges in the read graph. 
    // ii.  Process alignments in order of increasing Q. 
    // iii. If the alignment breaks strand separation, it is skipped. 
    // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped.  
    // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.
    
    // Maintain a vector containing the degree of each vertex
    // verticesDegree[vertexID] -> degree
    vector<uint64_t> verticesDegree(orientedReadCount, 0);
    vector<uint64_t> verticesDegreeGoodSupport(orientedReadCount, 0);
    vector<uint64_t> verticesDegree30to40(orientedReadCount, 0);


    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }
    
    cout << "Number of reads: " << readCount << endl;
    cout << "Number of oriented reads: " << orientedReadCount << endl;


    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);



    // Process alignments in order of increasing Q
    vector alignmentTablesToProcess({alignmentTable});
    
    uint64_t crossStrandEdgeCount = 0;

    for (auto alignmentTableToProcess : alignmentTablesToProcess) {
        for(auto it=alignmentTableToProcess.begin(); it!=alignmentTableToProcess.end(); ++it) {
            const pair<uint64_t, double>& p = *it;
            const uint64_t alignmentId = p.first;
            const double logQ = p.second;

            // Get the alignment data
            AlignmentData& alignment = alignmentData[alignmentId];
            const ReadId readId0 = alignment.readIds[0];
            const ReadId readId1 = alignment.readIds[1];
            const bool isSameStrand = alignment.isSameStrand;
            SHASTA_ASSERT(readId0 < readId1);
            const OrientedReadId A0 = OrientedReadId(readId0, 0);
            const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
            const OrientedReadId A1 = OrientedReadId(readId0, 1);
            const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

            SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
            SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
            SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
            SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

            // Get the connected components that these oriented reads are in.
            const uint32_t a0 = disjointSets.find_set(A0.getValue());
            const uint32_t b0 = disjointSets.find_set(B0.getValue());
            const uint32_t a1 = disjointSets.find_set(A1.getValue());
            const uint32_t b1 = disjointSets.find_set(B1.getValue());


            // If the alignment breaks strand separation, it is skipped.
            // If A0 and B1 are in the same connected component,
            // A1 and B0 also must be in the same connected component.
            // Adding this pair of edges would create a self-complementary
            // connected component containing A0, B0, A1, and B1,
            // and to ensure strand separation we don't want to do that.
            // So we mark these edges as cross-strand edges
            // and don't use them to update the disjoint set data structure.
            if(a0 == b1) {
                SHASTA_ASSERT(a1 == b0);
                crossStrandEdgeCount += 2;
                continue;
            }

            

            // If both vertices of the potential edge have at least the required minimum number 
            // of neighbors, the alignment is also skipped. 
            const uint64_t degreeA0 = verticesDegree[A0.getValue()];
            const uint64_t degreeB0 = verticesDegree[B0.getValue()];
            const uint64_t degreeA0GoodSupport = verticesDegreeGoodSupport[A0.getValue()];
            const uint64_t degreeB0GoodSupport = verticesDegreeGoodSupport[B0.getValue()];
            const uint64_t degreeA030to40 = verticesDegree30to40[A0.getValue()];
            const uint64_t degreeB030to40 = verticesDegree30to40[B0.getValue()];


            // if(degreeA030to40 >= 1 and degreeB030to40 >= 1) {
            //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
            //     continue;
            // }

            // if(degreeA0GoodSupport >= 10 and degreeB0GoodSupport >= 10) {
            //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
            //     continue;
            // }


            if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
                // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
                continue;
            }

            // Add the alignment to the read graph.
            keepAlignment[alignmentId] = true;
            alignment.info.isInReadGraph = 1;

            // Update vertex degrees
            verticesDegree[A0.getValue()]++;
            verticesDegree[B0.getValue()]++;
            verticesDegree[A1.getValue()]++;
            verticesDegree[B1.getValue()]++;
            
            // if(alignment.info.markerCount >= 1500) {
            //     verticesDegreeGoodSupport[A0.getValue()]++;
            //     verticesDegreeGoodSupport[B0.getValue()]++;
            //     verticesDegreeGoodSupport[A1.getValue()]++;
            //     verticesDegreeGoodSupport[B1.getValue()]++;
            // }

            // if(alignment.info.errorRateRle <= maxErrorRateRle30 and alignment.info.errorRateRle >= maxErrorRateRle40) {
            //     verticesDegree30to40[A0.getValue()]++;
            //     verticesDegree30to40[B0.getValue()]++;
            //     verticesDegree30to40[A1.getValue()]++;
            //     verticesDegree30to40[B1.getValue()]++;
            // }
            

            // Update disjoint sets
            disjointSets.union_set(a0, b0);
            disjointSets.union_set(a1, b1);

        }

        // Verify that for any read the two oriented reads are in distinct
        // connected components.
        for(ReadId readId=0; readId<readCount; readId++) {
            const OrientedReadId orientedReadId0(readId, 0);
            const OrientedReadId orientedReadId1(readId, 1);
            SHASTA_ASSERT(
                disjointSets.find_set(orientedReadId0.getValue()) !=
                disjointSets.find_set(orientedReadId1.getValue())
            );
        }
    }

    // Track the size of each set in the disjoint sets
    vector<uint64_t> setSizes(orientedReadCount, 0);
    for (int i = 0; i < orientedReadCount; ++i) {
        setSizes[disjointSets.find_set(i)]++;
    }


    // Print how many alignments were kept in this step
    const size_t keepCountR1 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Finding strict disjointSets step: Keeping " << keepCountR1 << " alignments of " << keepAlignment.size() << endl;















    //*
    //
    // Create the dynamically adjustable boost readGraph using the alignments we selected.
    //
    //*
    using boost::add_vertex;
    using boost::add_edge;

    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph4 readGraph(orientedReadCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // Record whether this alignment is used in the read graph.
        const bool keepThisAlignment = keepAlignment[alignmentId];
        const AlignmentData& alignment = alignmentData[alignmentId];

        // If this alignment is not used in the read graph, we are done.
        if(not keepThisAlignment) {
            continue;
        }

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(alignment.readIds[0], 0);
        OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

        // Swap them if necessary, depending on the average alignment offset at center.
        if(alignment.info.offsetAtCenter() < 0.) {
            swap(orientedReadId0, orientedReadId1);
        }

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    }


    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph4AllAlignments readGraphAllAlignments(orientedReadCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // Record whether this alignment is used in the read graph.
        const bool keepThisAlignment = keepAlignmentsForBreaks[alignmentId];
        const AlignmentData& alignment = alignmentData[alignmentId];

        // If this alignment is not used in the read graph, we are done.
        if(not keepThisAlignment) {
            continue;
        }

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(alignment.readIds[0], 0);
        OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

        // Swap them if necessary, depending on the average alignment offset at center.
        if(alignment.info.offsetAtCenter() < 0.) {
            swap(orientedReadId0, orientedReadId1);
        }

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
    }

    
    
    
    cout << "The read graph has " << num_vertices(readGraph) << " vertices and " << num_edges(readGraph) << " edges." << endl;





    //*
    //
    // Find possible dead end nodes on the directed graph.
    // We only keep the reads that have no outgoing neighbors.
    //
    //*
    vector<bool> potentialDeadEndReads(orientedReadCount, false);
    vector<bool> isolatedReads(orientedReadCount, false);

    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {
            OrientedReadId orientedReadId(readId, strand);
            
            // Find neighbors in the forward direction
            vector<OrientedReadId> forwardNeighbors;
            readGraph.findNeighborsDirectedGraphOneSideRight(orientedReadId, 1, forwardNeighbors);
            
            // Find neighbors in the backward direction 
            vector<OrientedReadId> leftNeighbors;
            readGraph.findNeighborsDirectedGraphOneSideLeft(orientedReadId, 1, leftNeighbors);

            if (forwardNeighbors.empty() && leftNeighbors.empty() ) {
                isolatedReads[orientedReadId.getValue()] = true;
                OrientedReadId reverseOrientedReadId = orientedReadId;
                reverseOrientedReadId.flipStrand();
                isolatedReads[reverseOrientedReadId.getValue()] = true;
                continue;
            }

            // If a read has neighbors only in the backward direction, it's a potential dead end
            if (forwardNeighbors.empty() && !leftNeighbors.empty()) {
                potentialDeadEndReads[orientedReadId.getValue()] = true;
            }
        }
    }

    // count the number of potential dead end reads
    uint64_t potentialDeadEndReadCount = count(potentialDeadEndReads.begin(), potentialDeadEndReads.end(), true);
    cout << "Found " << potentialDeadEndReadCount << " potential dead end reads." << endl;

    // // Print dead end Oriented reads
    // // iterate over all oriented reads
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     for (Strand strand = 0; strand < 2; strand++) {
    //         OrientedReadId orientedReadId(readId, strand);
    //         if (potentialDeadEndReads[orientedReadId.getValue()]) {
    //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a potential dead end read." << endl;
    //         }
    //     }
    // }






    //*
    //
    // Filter out the potential dead end nodes on the directed graph.
    // Look for a read path that lead to a positive offset.
    // If no such path is found, the orientedReadId is kept as a potential dead end node.
    //
    //*
    vector<bool> finalDeadEndReadsWithNoOutgoingNodes(orientedReadCount, false);
    vector<bool> finalDeadEndReadsWithNoIncomingNodes(orientedReadCount, false);
    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {

            OrientedReadId orientedReadId(readId, strand);

            OrientedReadId reverseOrientedReadId = orientedReadId;
            reverseOrientedReadId.flipStrand();

            // check if the orientedReadId is in potentialDeadEndReads
            if(!potentialDeadEndReads[orientedReadId.getValue()]) {
                continue;
            }

            // create the necessary variables for findAllPaths
            vector<vector<OrientedReadId>> paths;
            vector<vector<double>> pathsOffsets;
            vector<OrientedReadId> currentPath;
            vector<double> currentPathOffset;
            std::set<ReadGraph4BaseClass::vertex_descriptor> visited;
            uint64_t maxDistance = 5;
            uint64_t currentDistance = 0;

            uint64_t result = readGraph.findPathWithPositiveOffset(orientedReadId, paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

            // Check if we found a read path with positive offset.
            // If yes, the function findPathWithPositiveOffset will return 1, if not, it will return 0.
            if(result == 1) {
                // cout << "Found a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset" << endl;
                // //print the paths and then the pathsOffsets
                // for (size_t i = 0; i < paths.size(); i++) {
                //     cout << "Path " << i << ": ";
                //     for (size_t j = 0; j < paths[i].size(); j++) {
                //         cout << paths[i][j].getReadId() << " ";
                //     }
                //     cout << endl;
                //     cout << "PathOffsets " << i << ": ";
                //     for (size_t j = 0; j < pathsOffsets[i].size(); j++) {
                //         cout << pathsOffsets[i][j] << " ";
                //     }
                //     cout << endl;
                // }
            } else if (result == 0) {
                cout << "Did not find a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset. Keeping it as a potential dead end read." << endl;
                finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()] = true;
                finalDeadEndReadsWithNoIncomingNodes[reverseOrientedReadId.getValue()] = true;
            }
        }
    }


    // // print dead end Oriented reads
    // // iterate over all oriented reads
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     for (Strand strand = 0; strand < 2; strand++) {
    //         OrientedReadId orientedReadId(readId, strand);
    //         if (finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]) {
    //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
    //         } else if (finalDeadEndReadsWithNoIncomingNodes[orientedReadId.getValue()]) {
    //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
    //         }
    //     }
    // }











    //*
    //
    // Extend the potential dead end nodes list.
    // Add neighboring nodes of potential dead end nodes to the dead end node list.
    //
    //*
    vector<bool> finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors(orientedReadCount, false);
    vector<bool> finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors(orientedReadCount, false);
    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {
            OrientedReadId orientedReadId(readId, strand);

            if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){

                finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()] = true;
                
                OrientedReadId reverseOrientedReadId = orientedReadId;
                reverseOrientedReadId.flipStrand();

                finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;

                // Find distance 1 neighbors
                vector<OrientedReadId> distance1Neighbors;
                readGraph.findNeighborsDirectedGraphBothSides(orientedReadId, 1, distance1Neighbors);

                for(auto neighbor : distance1Neighbors) {
                    finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[neighbor.getValue()] = true;

                    OrientedReadId reverseOrientedReadId = neighbor;
                    reverseOrientedReadId.flipStrand();

                    finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;
                }
            }
        }
    }

    // // print dead end Oriented reads
    // // iterate over all oriented reads
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     for (Strand strand = 0; strand < 2; strand++) {
    //         OrientedReadId orientedReadId(readId, strand);
    //         if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
    //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
    //         } else if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
    //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
    //         }
    //     }
    // }










    



    //
    //
    // Find forbidden non used reads
    //
    //
    vector<bool> forbiddenReads(orientedReadCount, false);

    // // iterate over finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     for (Strand strand = 0; strand < 2; strand++) {
    //         // continue;
    //         OrientedReadId orientedReadId(readId, strand);
    //         // Check if the oriented read is a final dead end read with no outgoing nodes (nothing in its right side)
    //         if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {

    //             // Find neighbors in the forward direction
    //             vector<OrientedReadId> forwardNeighbors;
    //             readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, 5, forwardNeighbors);

    //             // cout << "Check 1 " << endl;

    //             // create a std::set of the forwardNeighbors for easy contain check
    //             std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
                

    //             if(forwardNeighbors.empty()) {
    //                 // cout << "Check 2 " << endl;
    //                 continue;
    //             }

    //             // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
    //             // or a non speficic node if we exceeded maxDistance
    //             OrientedReadId lastNode = forwardNeighbors.back();
                

    //             //
    //             // Ensure the last node does not belong to the other disjoint set
    //             //

    //             // First check if the last node is a potential dead end node with no INCOMING nodes
    //             if(finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[lastNode.getValue()]) {
                    

    //                 const OrientedReadId A0 = orientedReadId;
    //                 const OrientedReadId B0 = lastNode;
    //                 const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
    //                 const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

    //                 SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
    //                 SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
    //                 SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
    //                 SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

    //                 // Get the connected components that these oriented reads are in.
    //                 const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //                 const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //                 const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //                 const uint32_t b1 = disjointSets.find_set(B1.getValue());


    //                 // If the last node belongs to the disjoint set of the other strand
    //                 if(a0 != b0 and a0 == b1) {

    //                     cout << "Found bad actors " << endl;

    //                     // vector<vector<OrientedReadId>> paths;
    //                     // vector<OrientedReadId> currentPath;
    //                     // std::set<ReadGraph4AllAlignmentsBaseClass::vertex_descriptor> visited;

    //                     // Find the shortest path from A0 to B1
    //                     vector<OrientedReadId> shortestPath;
    //                     uint64_t maxDistance = 10;
    //                     readGraphAllAlignments.findShortestPathToNode(A0, B1, shortestPath, maxDistance);

    //                     // Print the shortest path
    //                     cout << "Shortest path from A0 to B1: ";
    //                     for (const OrientedReadId& orientedReadId : shortestPath) {
    //                         cout << "ReadID " << orientedReadId.getReadId() << " Strand " << orientedReadId.getStrand() << " ";
    //                     }
    //                     cout << endl;

    //                     // Iterate over the alignment ids and if the alignment has reads that are in the shortest path, set the keepAlignment to true
    //                     // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //                     for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {

    //                         const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
    //                         const uint64_t alignmentId = p.first;
    //                         const double logQ = p.second;

    //                         const bool keepThisAlignment = keepAlignment[alignmentId];

    //                         const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

    //                         if(keepThisAlignment) {
    //                             continue;
    //                         }

    //                         if(not keepThisBreaksAlignment) {
    //                             continue;
    //                         }

    //                         AlignmentData& alignment = alignmentData[alignmentId];
                    

    //                         // Get the alignment data
    //                         ReadId readId0v2 = alignment.readIds[0];
    //                         ReadId readId1v2 = alignment.readIds[1];
    //                         const bool isSameStrandv2 = alignment.isSameStrand;
    //                         SHASTA_ASSERT(readId0v2 < readId1v2);
    //                         OrientedReadId A0v2 = OrientedReadId(readId0v2, 0);
    //                         OrientedReadId B0v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 0 : 1);
    //                         OrientedReadId A1v2 = OrientedReadId(readId0v2, 1);
    //                         OrientedReadId B1v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 1 : 0);

    //                         SHASTA_ASSERT(A0v2.getReadId() == A1v2.getReadId());
    //                         SHASTA_ASSERT(B0v2.getReadId() == B1v2.getReadId());
    //                         SHASTA_ASSERT(A0v2.getStrand() == 1 - A1v2.getStrand());
    //                         SHASTA_ASSERT(B0v2.getStrand() == 1 - B1v2.getStrand());

    //                         // Get the connected components that these oriented reads are in.
    //                         const uint32_t a0v2 = disjointSets.find_set(A0v2.getValue());
    //                         const uint32_t b0v2 = disjointSets.find_set(B0v2.getValue());
    //                         const uint32_t a1v2 = disjointSets.find_set(A1v2.getValue());
    //                         const uint32_t b1v2 = disjointSets.find_set(B1v2.getValue());

    //                         // If the alignment breaks strand separation, it is skipped.
    //                         // If A0 and B1 are in the same connected component,
    //                         // A1 and B0 also must be in the same connected component.
    //                         // Adding this pair of edges would create a self-complementary
    //                         // connected component containing A0, B0, A1, and B1,
    //                         // and to ensure strand separation we don't want to do that.
    //                         // So we mark these edges as cross-strand edges
    //                         // and don't use them to update the disjoint set data structure.
    //                         if(a0v2 == b1v2) {
    //                             SHASTA_ASSERT(a1 == b0);
    //                             crossStrandEdgeCount += 2;
    //                             continue;
    //                         }

    //                         // Swap them if necessary, depending on the average alignment offset at center.
    //                         if(alignment.info.offsetAtCenter() < 0.) {
    //                             swap(A0v2, B0v2);
    //                         }

    //                         // Check if A0v2 is in the shortest path and check if B0v2 is exactly after A0v2 in the path
    //                         bool a0v2InPath = false;
    //                         bool b0v2InPath = false;
    //                         for(uint64_t i=0; i<shortestPath.size()-1; i++) {
    //                             if(shortestPath[i] == A0v2) {
    //                                 a0v2InPath = true;
    //                             }
    //                             if(shortestPath[i+1] == B0v2) {
    //                                 b0v2InPath = true;
    //                             }
    //                         }

    //                         if(a0v2InPath && b0v2InPath) {
    //                             // Add the alignment to the read graph.
    //                             keepAlignment[alignmentId] = true;
    //                             alignment.info.isInReadGraph = 1;

    //                             // Update vertex degrees
    //                             verticesDegree[A0v2.getValue()]++;
    //                             verticesDegree[B0v2.getValue()]++;
    //                             verticesDegree[A1v2.getValue()]++;
    //                             verticesDegree[B1v2.getValue()]++;

    //                             disjointSets.union_set(a0v2, b0v2);
    //                             disjointSets.union_set(a1v2, b1v2);

    //                             finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[A0.getValue()] = false;
    //                             finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[A1.getValue()] = false;
    //                             finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[B0.getValue()] = false;
    //                             finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[B1.getValue()] = false;
                                                               

    //                             // Update disjoint sets
    //                             //disjointSets.union_set(a0, a0v2);
    //                             // disjointSets.union_set(a0, b0v2);

    //                             //disjointSets.union_set(a1, a1v2);
    //                             // disjointSets.union_set(a1, b1v2);

    //                             cout << "Adding alignment " << alignmentId << " between " << A0v2.getReadId() << " and " << B0v2.getReadId() << endl;

                                
    //                         }
    //                     }

                        

    //                 }
    //             }

    //         }
    //     }
    // }



















    // iterate over finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors
    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {
            // continue;
            OrientedReadId orientedReadId(readId, strand);
            // Check if the oriented read is a final dead end read with no outgoing nodes (nothing in its right side)
            if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {

                // Find neighbors in the forward direction
                vector<OrientedReadId> forwardNeighbors;
                readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, 5, forwardNeighbors);

                // cout << "Check 1 " << endl;

                // create a std::set of the forwardNeighbors for easy contain check
                std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
                

                if(forwardNeighbors.empty()) {
                    // cout << "Check 2 " << endl;
                    continue;
                }

                // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
                // or a non speficic node if we exceeded maxDistance
                OrientedReadId lastNode = forwardNeighbors.back();
                

                //
                // Ensure the last node does not belong to the other disjoint set
                //

                // First check if the last node is a potential dead end node with no INCOMING nodes
                if(finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[lastNode.getValue()]) {
                    
                    const OrientedReadId A0 = orientedReadId;
                    const OrientedReadId B0 = lastNode;
                    const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
                    const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

                    SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
                    SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
                    SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
                    SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

                    // Get the connected components that these oriented reads are in.
                    const uint32_t a0 = disjointSets.find_set(A0.getValue());
                    const uint32_t b0 = disjointSets.find_set(B0.getValue());
                    const uint32_t a1 = disjointSets.find_set(A1.getValue());
                    const uint32_t b1 = disjointSets.find_set(B1.getValue());






                        // // Find distance 1 neighbors
                        // vector<OrientedReadId> startNodeDistance1Neighbors;
                        // readGraphAllAlignments.findNeighborsDirectedGraphBothSides(orientedReadId, 5, startNodeDistance1Neighbors);
                        
                        // // Find distance 1 neighbors
                        // vector<OrientedReadId> lastNodeDistance1Neighbors;
                        // readGraphAllAlignments.findNeighborsDirectedGraphBothSides(lastNode, 5, lastNodeDistance1Neighbors);



                        // // Find distinct neighbors of the start node
                        // // Vector to store neighbors unique to start node
                        // vector<OrientedReadId> distinctNeighborStartNode;

                        // // Iterate through all neighbors of start node
                        // for(auto neighbor : startNodeDistance1Neighbors) {
                        //     // Check if this neighbor is NOT in lastNodeDistance1Neighbors
                        //     if(std::find(lastNodeDistance1Neighbors.begin(), 
                        //                 lastNodeDistance1Neighbors.end(), 
                        //                 neighbor) == lastNodeDistance1Neighbors.end()) {
                        //         // If not found, add it to our distinct neighbors list
                        //         distinctNeighborStartNode.push_back(neighbor);
                        //     }
                        // }

                        // std::set<OrientedReadId> distinctNeighborStartNodeSet(distinctNeighborStartNode.begin(), distinctNeighborStartNode.end());

                        // distinctNeighborStartNodeSet.insert(orientedReadId);





                        // // Find distinct neighbors of the last node
                        // // Vector to store neighbors unique to last node
                        // vector<OrientedReadId> distinctNeighborLastNode;

                        // // Iterate through all neighbors of start node
                        // for(auto neighbor : startNodeDistance1Neighbors) {
                        //     // Check if this neighbor is NOT in lastNodeDistance1Neighbors
                        //     if(std::find(startNodeDistance1Neighbors.begin(), 
                        //                 startNodeDistance1Neighbors.end(), 
                        //                 neighbor) == startNodeDistance1Neighbors.end()) {
                        //         // If not found, add it to our distinct neighbors list
                        //         distinctNeighborLastNode.push_back(neighbor);
                        //     }
                        // }

                        // std::set<OrientedReadId> distinctNeighborLastNodeSet(distinctNeighborLastNode.begin(), distinctNeighborLastNode.end());

                        // distinctNeighborLastNodeSet.insert(lastNode);



                        // for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
                        // // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

                        //     const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
                        //     const uint64_t alignmentId = p.first;
                        //     const double logQ = p.second;

                        //     const bool keepThisAlignment = keepAlignment[alignmentId];

                        //     const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

                        //     if(keepThisAlignment) {
                        //         continue;
                        //     }

                        //     if(not keepThisBreaksAlignment) {
                        //         continue;
                        //     }

                        //     AlignmentData& alignment = alignmentData[alignmentId];
                        
                        //     // Get the OrientedReadIds.
                        //     OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
                        //     OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
                        //     SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

                        //     // Swap them if necessary, depending on the average alignment offset at center.
                        //     if(alignment.info.offsetAtCenter() < 0.) {
                        //         swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
                        //     }

                        //     if((alignmentOrientedReadId0.getValue() == orientedReadId.getValue() and distinctNeighborStartNodeSet.contains(alignmentOrientedReadId1)) 
                        //     ||  (distinctNeighborLastNodeSet.contains(alignmentOrientedReadId0) and alignmentOrientedReadId1.getValue() == lastNode.getValue()) 
                        //     || distinctNeighborLastNodeSet.contains(alignmentOrientedReadId0) and distinctNeighborLastNodeSet.contains(alignmentOrientedReadId1) 
                        //     || distinctNeighborStartNodeSet.contains(alignmentOrientedReadId0) and distinctNeighborStartNodeSet.contains(alignmentOrientedReadId1)) {
                                
                        //         // Get the alignment data
                        //         const ReadId readId0 = alignment.readIds[0];
                        //         const ReadId readId1 = alignment.readIds[1];
                        //         const bool isSameStrand = alignment.isSameStrand;
                        //         SHASTA_ASSERT(readId0 < readId1);
                        //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
                        //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
                        //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
                        //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

                        //         SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
                        //         SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
                        //         SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
                        //         SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

                        //         // Get the connected components that these oriented reads are in.
                        //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
                        //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
                        //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
                        //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

                        //         // If the alignment breaks strand separation, it is skipped.
                        //         // If A0 and B1 are in the same connected component,
                        //         // A1 and B0 also must be in the same connected component.
                        //         // Adding this pair of edges would create a self-complementary
                        //         // connected component containing A0, B0, A1, and B1,
                        //         // and to ensure strand separation we don't want to do that.
                        //         // So we mark these edges as cross-strand edges
                        //         // and don't use them to update the disjoint set data structure.
                        //         if(a0 == b1) {
                        //             SHASTA_ASSERT(a1 == b0);
                        //             crossStrandEdgeCount += 2;
                        //             continue;
                        //         }

                        //         // If both vertices of the potential edge have at least the required minimum number 
                        //         // of neighbors, the alignment is also skipped. 
                        //         const uint64_t degreeA0 = verticesDegree[A0.getValue()];
                        //         const uint64_t degreeB0 = verticesDegree[B0.getValue()];

                    
                        //         // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
                        //         //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
                        //         //     continue;
                        //         // }

                        //         // Add the alignment to the read graph.
                        //         keepAlignment[alignmentId] = true;
                        //         alignment.info.isInReadGraph = 1;

                        //         // Update vertex degrees
                        //         verticesDegree[A0.getValue()]++;
                        //         verticesDegree[B0.getValue()]++;
                        //         verticesDegree[A1.getValue()]++;
                        //         verticesDegree[B1.getValue()]++;
                                

                        //         // Update disjoint sets
                        //         disjointSets.union_set(a0, b0);
                        //         disjointSets.union_set(a1, b1);

                        //         // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

                        //         // Create the edge.
                        //         add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                        //         // Also create the reverse complemented edge.
                        //         alignmentOrientedReadId0.flipStrand();
                        //         alignmentOrientedReadId1.flipStrand();
                        //         add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                                
                        //     }
                        // }


                        // // The last node belongs to the disjoint set of the other strand
                        // finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[lastNode.getValue()] = false;

                        // // Remove them from the consideration
                        // for(auto neighbor : distinctNeighborLastNode) {
                        //     finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[neighbor.getValue()] = false;
                        // }

                        // // Find distance 1 neighbors
                        // vector<OrientedReadId> lastNodeReversedDistance1Neighbors;
                        // readGraph.findNeighborsDirectedGraphBothSides(B1, 1, lastNodeReversedDistance1Neighbors);

                        
                        // // Add them to the consideration
                        // for(auto neighbor : lastNodeReversedDistance1Neighbors) {
                        //     finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[neighbor.getValue()] = true;
                        // }
                        

                        // // Update the forwardNeighborsSet
                        // forwardNeighborsSet.clear();

                        // // Find neighbors in the forward direction
                        // readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, 5, forwardNeighbors);
                        
                        // forwardNeighborsSet.insert(forwardNeighbors.begin(), forwardNeighbors.end());


                        // if(forwardNeighbors.empty()) {
                        //     // cout << "Check 2 " << endl;
                        //     continue;
                        // }

                        // for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
                        // // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

                        //     const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
                        //     const uint64_t alignmentId = p.first;
                        //     const double logQ = p.second;

                        //     const bool keepThisAlignment = keepAlignment[alignmentId];

                        //     const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

                        //     if(keepThisAlignment) {
                        //         continue;
                        //     }

                        //     if(not keepThisBreaksAlignment) {
                        //         continue;
                        //     }

                        //     AlignmentData& alignment = alignmentData[alignmentId];
                        
                        //     // Get the OrientedReadIds.
                        //     OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
                        //     OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
                        //     SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

                        //     // Swap them if necessary, depending on the average alignment offset at center.
                        //     if(alignment.info.offsetAtCenter() < 0.) {
                        //         swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
                        //     }

                        //     // // convert the shortestPath vector to a set and check if it contains alignmentId
                        //     // std::set<uint32_t> shortestPathSet(shortestPath.begin(), shortestPath.end());
                        //     // if(shortestPathSet.contains(alignmentId)) {
                        //     //     keepAlignment[alignmentId] = true;
                        //     //     alignment.info.isInReadGraph = 1;
                        //     //     // Create the edge.
                        //     //     add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                        //     //     // Also create the reverse complemented edge.
                        //     //     alignmentOrientedReadId0.flipStrand();
                        //     //     alignmentOrientedReadId1.flipStrand();
                        //     //     add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                        //     //     continue;
                        //     // }

                        //     if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                                
                        //         // Get the alignment data
                        //         const ReadId readId0 = alignment.readIds[0];
                        //         const ReadId readId1 = alignment.readIds[1];
                        //         const bool isSameStrand = alignment.isSameStrand;
                        //         SHASTA_ASSERT(readId0 < readId1);
                        //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
                        //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
                        //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
                        //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

                        //         SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
                        //         SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
                        //         SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
                        //         SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

                        //         // Get the connected components that these oriented reads are in.
                        //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
                        //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
                        //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
                        //         const uint32_t b1 = disjointSets.find_set(B1.getValue());


                        //         // If the alignment breaks strand separation, it is skipped.
                        //         // If A0 and B1 are in the same connected component,
                        //         // A1 and B0 also must be in the same connected component.
                        //         // Adding this pair of edges would create a self-complementary
                        //         // connected component containing A0, B0, A1, and B1,
                        //         // and to ensure strand separation we don't want to do that.
                        //         // So we mark these edges as cross-strand edges
                        //         // and don't use them to update the disjoint set data structure.
                        //         if(a0 == b1) {
                        //             SHASTA_ASSERT(a1 == b0);
                        //             crossStrandEdgeCount += 2;
                        //             continue;
                        //         }

                                

                        //         // If both vertices of the potential edge have at least the required minimum number 
                        //         // of neighbors, the alignment is also skipped. 
                        //         const uint64_t degreeA0 = verticesDegree[A0.getValue()];
                        //         const uint64_t degreeB0 = verticesDegree[B0.getValue()];

                    
                        //         // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
                        //         //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
                        //         //     continue;
                        //         // }

                        //         // Add the alignment to the read graph.
                        //         keepAlignment[alignmentId] = true;
                        //         alignment.info.isInReadGraph = 1;

                        //         // Update vertex degrees
                        //         verticesDegree[A0.getValue()]++;
                        //         verticesDegree[B0.getValue()]++;
                        //         verticesDegree[A1.getValue()]++;
                        //         verticesDegree[B1.getValue()]++;
                                

                        //         // Update disjoint sets
                        //         disjointSets.union_set(a0, b0);
                        //         disjointSets.union_set(a1, b1);

                        //         // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

                        //         // Create the edge.
                        //         add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                        //         // Also create the reverse complemented edge.
                        //         alignmentOrientedReadId0.flipStrand();
                        //         alignmentOrientedReadId1.flipStrand();
                        //         add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


                        //     }

                        // }
                        
                        // continue;

                    

             



                    // for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
                    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

                        // const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
                        // const uint64_t alignmentId = p.first;
                        // const double logQ = p.second;

                        const bool keepThisAlignment = keepAlignment[alignmentId];

                        const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

                        if(keepThisAlignment) {
                            continue;
                        }

                        if(not keepThisBreaksAlignment) {
                            continue;
                        }

                        AlignmentData& alignment = alignmentData[alignmentId];
                    
                        // Get the OrientedReadIds.
                        OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
                        OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
                        SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

                        // Swap them if necessary, depending on the average alignment offset at center.
                        if(alignment.info.offsetAtCenter() < 0.) {
                            swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
                        }

                        // // convert the shortestPath vector to a set and check if it contains alignmentId
                        // std::set<uint32_t> shortestPathSet(shortestPath.begin(), shortestPath.end());
                        // if(shortestPathSet.contains(alignmentId)) {
                        //     keepAlignment[alignmentId] = true;
                        //     alignment.info.isInReadGraph = 1;
                        //     // Create the edge.
                        //     add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                        //     // Also create the reverse complemented edge.
                        //     alignmentOrientedReadId0.flipStrand();
                        //     alignmentOrientedReadId1.flipStrand();
                        //     add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                        //     continue;
                        // }

                        if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                            
                            // keepAlignment[alignmentId] = true;
                            // alignment.info.isInReadGraph = 1;

                            // continue;

                            // Get the alignment data
                            ReadId readId0v2 = alignment.readIds[0];
                            ReadId readId1v2 = alignment.readIds[1];
                            const bool isSameStrandv2 = alignment.isSameStrand;
                            SHASTA_ASSERT(readId0v2 < readId1v2);
                            OrientedReadId A0v2 = OrientedReadId(readId0v2, 0);
                            OrientedReadId B0v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 0 : 1);
                            OrientedReadId A1v2 = OrientedReadId(readId0v2, 1);
                            OrientedReadId B1v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 1 : 0);

                            SHASTA_ASSERT(A0v2.getReadId() == A1v2.getReadId());
                            SHASTA_ASSERT(B0v2.getReadId() == B1v2.getReadId());
                            SHASTA_ASSERT(A0v2.getStrand() == 1 - A1v2.getStrand());
                            SHASTA_ASSERT(B0v2.getStrand() == 1 - B1v2.getStrand());

                            // Get the connected components that these oriented reads are in.
                            const uint32_t a0v2 = disjointSets.find_set(A0v2.getValue());
                            const uint32_t b0v2 = disjointSets.find_set(B0v2.getValue());
                            const uint32_t a1v2 = disjointSets.find_set(A1v2.getValue());
                            const uint32_t b1v2 = disjointSets.find_set(B1v2.getValue());


                            // If the alignment breaks strand separation, it is skipped.
                            // If A0 and B1 are in the same connected component,
                            // A1 and B0 also must be in the same connected component.
                            // Adding this pair of edges would create a self-complementary
                            // connected component containing A0, B0, A1, and B1,
                            // and to ensure strand separation we don't want to do that.
                            // So we mark these edges as cross-strand edges
                            // and don't use them to update the disjoint set data structure.
                            if(a0v2 == b1v2) {
                                SHASTA_ASSERT(a1v2 == b0v2);
                                crossStrandEdgeCount += 2;
                                continue;
                            }

                            

                            
                                                            

                            // Update disjoint sets
                            //disjointSets.union_set(a0, a0v2);
                            // disjointSets.union_set(a0, b0v2);

                            //disjointSets.union_set(a1, a1v2);
                            // disjointSets.union_set(a1, b1v2);

                            // cout << "Adding alignment " << alignmentId << " between " << A0v2.getReadId() << " and " << B0v2.getReadId() << endl;

                            

                            // If both vertices of the potential edge have at least the required minimum number 
                            // of neighbors, the alignment is also skipped. 
                            // const uint64_t degreeA0 = verticesDegree[A0.getValue()];
                            // const uint64_t degreeB0 = verticesDegree[B0.getValue()];

                
                            // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
                            //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
                            //     continue;
                            // }

                            // Add the alignment to the read graph.
                            keepAlignment[alignmentId] = true;
                            alignment.info.isInReadGraph = 1;

                            // Update vertex degrees
                            verticesDegree[A0v2.getValue()]++;
                            verticesDegree[B0v2.getValue()]++;
                            verticesDegree[A1v2.getValue()]++;
                            verticesDegree[B1v2.getValue()]++;
                            

                            // Update disjoint sets
                            disjointSets.union_set(a0v2, b0v2);
                            disjointSets.union_set(a1v2, b1v2);

                            finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[A0.getValue()] = false;
                            finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[A1.getValue()] = false;
                            finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[B0.getValue()] = false;
                            finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[B1.getValue()] = false;


                            // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

                            // Create the edge.
                            add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                            // Also create the reverse complemented edge.
                            alignmentOrientedReadId0.flipStrand();
                            alignmentOrientedReadId1.flipStrand();
                            add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


                        }

                    }
                    
                }

            }
        }
    }
        
    



    






    // // iterate over finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     for (Strand strand = 0; strand < 2; strand++) {
    //         OrientedReadId orientedReadId(readId, strand);
    //         // Check if the oriented read is a final dead end read with no outgoing nodes (nothing in its right side)
    //         if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {

    //             // Find neighbors in the forward direction
    //             vector<OrientedReadId> forwardNeighbors;
    //             readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, 5, forwardNeighbors);

    //             // cout << "Check 1 " << endl;

    //             // create a std::set of the forwardNeighbors for easy contain check
    //             std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
                

    //             if(forwardNeighbors.empty()) {
    //                 // cout << "Check 2 " << endl;
    //                 continue;
    //             }

    //             // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
    //             // or a non speficic node if we exceeded maxDistance
    //             OrientedReadId lastNode = forwardNeighbors.back();

    //             // Check if this node is in the finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors
    //             if(finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[lastNode.getValue()]) {
    //                 // We found a path from a dead end node to another dead end node
    //                 // cout << "Check 3 " << endl;
    //                 // We can remove the dead end node
    //                 // finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()] = false;
    //                 // finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[lastNode.getValue()] = false;

    //                 // Find all alignments involving nodes in forwardNeighbors that have not been considered
                    
    //                 //for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
    //                 for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

    //                     // const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
    //                     // const uint64_t alignmentId = p.first;
    //                     // const double logQ = p.second;

    //                     const bool keepThisAlignment = keepAlignment[alignmentId];

    //                     const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

    //                     if(keepThisAlignment) {
    //                         continue;
    //                     }

    //                     if(not keepThisBreaksAlignment) {
    //                         continue;
    //                     }

    //                     AlignmentData& alignment = alignmentData[alignmentId];
                    
    //                     // Get the OrientedReadIds.
    //                     OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
    //                     OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //                     SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

    //                     // Swap them if necessary, depending on the average alignment offset at center.
    //                     if(alignment.info.offsetAtCenter() < 0.) {
    //                         swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
    //                     }

    //                     if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                            
    //                         // Get the alignment data
    //                         const ReadId readId0 = alignment.readIds[0];
    //                         const ReadId readId1 = alignment.readIds[1];
    //                         const bool isSameStrand = alignment.isSameStrand;
    //                         SHASTA_ASSERT(readId0 < readId1);
    //                         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //                         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //                         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //                         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //                         SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
    //                         SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
    //                         SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
    //                         SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

    //                         // Get the connected components that these oriented reads are in.
    //                         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //                         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //                         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //                         const uint32_t b1 = disjointSets.find_set(B1.getValue());


    //                         // If the alignment breaks strand separation, it is skipped.
    //                         // If A0 and B1 are in the same connected component,
    //                         // A1 and B0 also must be in the same connected component.
    //                         // Adding this pair of edges would create a self-complementary
    //                         // connected component containing A0, B0, A1, and B1,
    //                         // and to ensure strand separation we don't want to do that.
    //                         // So we mark these edges as cross-strand edges
    //                         // and don't use them to update the disjoint set data structure.
    //                         if(a0 == b1) {
    //                             SHASTA_ASSERT(a1 == b0);
    //                             crossStrandEdgeCount += 2;
    //                             continue;
    //                         }

                            

    //                         // If both vertices of the potential edge have at least the required minimum number 
    //                         // of neighbors, the alignment is also skipped. 
    //                         const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                         const uint64_t degreeB0 = verticesDegree[B0.getValue()];

                
    //                         // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                         //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
    //                         //     continue;
    //                         // }

    //                         // Add the alignment to the read graph.
    //                         keepAlignment[alignmentId] = true;
    //                         alignment.info.isInReadGraph = 1;

    //                         // Update vertex degrees
    //                         verticesDegree[A0.getValue()]++;
    //                         verticesDegree[B0.getValue()]++;
    //                         verticesDegree[A1.getValue()]++;
    //                         verticesDegree[B1.getValue()]++;
                            

    //                         // Update disjoint sets
    //                         disjointSets.union_set(a0, b0);
    //                         disjointSets.union_set(a1, b1);

    //                         cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

    //                         // Create the edge.
    //                         add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

    //                         // Also create the reverse complemented edge.
    //                         alignmentOrientedReadId0.flipStrand();
    //                         alignmentOrientedReadId1.flipStrand();
    //                         add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


    //                     }

    //                 }
                    

    //             }
                
    //         }
    //     }
    // }


    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }


    // Print how many alignments were kept
    const size_t keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Adding alignments for break bridging: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;
    

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);


    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;

    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    
    cout << "The read graph has " << componentMap.size() << " connected components." << endl;




















































































































    // //*
    // //
    // // Now it is time to find the breaks in the readGraph and add the
    // // alignments between two added reads that passed the initial strict 
    // // logQ filter and these reads are IN A DIFFERENT COMPONENT.
    // //
    // //*
    // bool weFoundAlignmentsInBreaksToAdd = true;
    // uint64_t allignmentsAdded = 0;
    // while (weFoundAlignmentsInBreaksToAdd) {
    //     weFoundAlignmentsInBreaksToAdd = false;
    //     cout << "We added " << allignmentsAdded << " Alignments" << endl;
    //     for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //         if((i % 100000) == 0) {
    //             cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //         }

    //         const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //         const double logQ = alignmentTableNotPassFilter[i].second;
    //         const AlignmentData& alignment = alignmentData[alignmentId];

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // Case 1a: There is an alignment between two added reads
    //         // that passed the initial strict logQ filter and these reads are 
    //         // IN THE SAME CONNECTED COMPONENT.
    //         // This case checks for direct connections between to nodes in a break
    //         // using alignments that did not pass the initial strict logQ filter.
    //         // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
    //         if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && a0 == b0 && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
    //             // Find neighbors of the orientedReadId in the filtered readGraph.
    //             // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
    //             // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
                
    //             uint64_t maxDistance = 10;
    //             vector<OrientedReadId> readGraphNeighbors;
    //             readGraphAllAlignments.findNeighborsEarlyStopWhenReachSameComponentNode(A0, disjointSets, maxDistance, readGraphNeighbors);

    //             // Get the last node in the neighbors, which contains the first node in the same component as A0
    //             OrientedReadId lastNode = readGraphNeighbors.back();

    //             // Verify it is in the same component as A0
    //             if(disjointSets.find_set(lastNode.getValue()) == a0) {
    //                 SHASTA_ASSERT(a1 == b0);
    //             }

    //             // Check if it is in the same component as A0
    //             if(disjointSets.find_set(lastNode.getValue()) != a0) {
    //                 // We did not found a node in the same component as A0
    //                 // Tha means that this is probably not a break situation
    //                 // Or that we did not search further enough
    //                 continue;

    //             }

    //             // Check if it is in the same component as A0
                
    //             endNodes[A0.getValue()] = true;
    //             endNodes[B0.getValue()] = true;
    //             add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //             keepAlignment[alignmentId] = true;
    //             weFoundAlignmentsInBreaksToAdd = true;
    //             allignmentsAdded++;
    //             alignmentTableNotPassFilterAlreadyConsidered[i] = true;
                

    //             // 


    //             // Check if B0 is NOT in the neighbors
    //             if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //                 endNodes[A0.getValue()] = true;
    //                 endNodes[B0.getValue()] = true;
    //                 add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //                 keepAlignment[alignmentId] = true;
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //                 allignmentsAdded++;
    //                 alignmentTableNotPassFilterAlreadyConsidered[i] = true;
    //             }
    //             continue;
    //         }
    //     }
    // }





































































        


    

    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

    //     // Record whether this alignment is used in the read graph.
    //     const bool keepThisAlignment = keepAlignmentsForBreaks[alignmentId];
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // If this alignment is not used in the read graph, we are done.
    //     if(not keepThisAlignment) {
    //         continue;
    //     }

    //     // Get the OrientedReadIds.
    //     OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //     OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //     SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //     // Swap them if necessary, depending on the average alignment offset at center.
    //     if(alignment.info.offsetAtCenter() < 0.) {
    //         swap(orientedReadId0, orientedReadId1);
    //     }

    //     // Create the edge.
    //     add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

    //     // Also create the reverse complemented edge.
    //     orientedReadId0.flipStrand();
    //     orientedReadId1.flipStrand();
    //     add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
    // }








    
    
    // // Find all reads that are dead ends (have neighbors only on one side)
    // vector<bool> deadEndNodes(readCount, false);

    // // Check each vertex in the read graph
    // vector<bool> readIdConsideredForBreakAndFailed(readCount, false);

    // for (ReadId readId = 0; readId < readCount; readId++) {

    //     OrientedReadId orientedReadId(readId, 0);

    //     if(readIdConsideredForBreakAndFailed[orientedReadId.getReadId()]) {
    //         continue;
    //     }

    //     vector<OrientedReadId> orientedReadsToCheckNeighbors;
    //     orientedReadsToCheckNeighbors.push_back(orientedReadId);

    //     vector<OrientedReadId> orientedReadsToCheckNeighborsTemp;

    //     std::map<OrientedReadId, double> offsetsNeighbors;
    //     offsetsNeighbors[orientedReadId] = 0.;

    //     vector<bool> readIdAlreadyConsidered(readCount, false);
    //     readIdAlreadyConsidered[orientedReadId.getReadId()] = true;

    //     bool isEndNode = true;
    //     bool isIsolatedNode = false;

    //     uint64_t iterations = 0;

    //     cout << "Check1 " << endl;

    //     while((isEndNode and iterations <= 4) and !isIsolatedNode) {

    //         iterations++;

    //         for (auto orientedReadIdToCheckNeighbors : orientedReadsToCheckNeighbors) {

    //             readIdAlreadyConsidered[orientedReadIdToCheckNeighbors.getReadId()] = true;

    //             // Find neighbors on distance 1
    //             vector<OrientedReadId> neighbors;
    //             readGraph.findNeighborsDirectedGraphBothSides(orientedReadIdToCheckNeighbors, 1, neighbors);

    //             if(neighbors.empty()) {
    //                 isIsolatedNode = true;
    //                 break;
    //             }

    //             cout << "Check2 " << endl;

    //             for(auto neighbor : neighbors) {

    //                 if(readIdAlreadyConsidered[neighbor.getReadId()]) {
    //                     continue;
    //                 }

    //                 // Get the alignmentId between two orientedReadIDs
    //                 std::pair<ReadGraph4::edge_descriptor, bool> edge_pair = boost::edge(
    //                     orientedReadIdToCheckNeighbors.getValue(), 
    //                     neighbor.getValue(), 
    //                     readGraph);

    //                 if (!edge_pair.second) {
    //                     continue;
    //                 }

    //                 cout << "Check3 " << endl;

    //                 uint32_t alignmentId = readGraph[edge_pair.first].alignmentId;

    //                 cout << "Check3.5 " << endl;

    //                 const AlignmentData& alignment = alignmentData[alignmentId];

    //                 double offset = alignment.info.offsetAtCenter();

    //                 offsetsNeighbors[neighbor] = offsetsNeighbors[orientedReadIdToCheckNeighbors] + offset;

    //                 cout << "Check4 " << endl;

    //                 if(offsetsNeighbors[neighbor] < 0.) {
    //                     readIdConsideredForBreakAndFailed[neighbor.getReadId()] = true;
    //                     orientedReadsToCheckNeighborsTemp.push_back(neighbor);
    //                 } else {
    //                     readIdConsideredForBreakAndFailed[orientedReadId.getReadId()] = true;
    //                     isEndNode = false;
    //                 }

    //                 cout << "Check5 " << endl;

                    
    //             }

    //             if(!isEndNode) {
    //                 break;
    //             }
    //         }

    //         if(isEndNode and !isIsolatedNode) {
    //             orientedReadsToCheckNeighbors.clear();
    //             orientedReadsToCheckNeighbors = orientedReadsToCheckNeighborsTemp;
    //             orientedReadsToCheckNeighborsTemp.clear();
    //             cout << "Check6 " << endl;
    //         }
            
    //     }

    //     if(isEndNode || isIsolatedNode) {
    //         deadEndNodes[orientedReadId.getReadId()] = true;
    //         cout << "Check7 " << endl;
    //     }

    // }


    // cout << "Found " << deadEndNodes.size() << " dead end nodes in the read graph." << endl;

    // // print dead end nodes
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     if(deadEndNodes[readId]) {
    //         cout << "Dead end node: " << readId << endl;
    //     }
    // }





    // vector<pair<uint64_t, size_t>> edgesToRemove;

    // // Keep only the top 10 neighbors with the highest markerCounts
    // BGL_FORALL_VERTICES(v, readGraph, ReadGraph4) {
    //     // Create a vector to store for each vertex the all of it's neighrbors markerCounts in distance 1
    //     vector<pair<uint64_t, double>> neighborsMarkerCount;
    //     vector<pair<uint64_t, double>> neighborsMarkerCount2;

    //     const ReadGraph4Vertex& vertex = readGraph[v];
    //     const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
        
    //     //cout << "Vertex " << orientedReadId << " has degree " << out_degree(v, readGraph) << endl;


    //     if (verticesDegree[orientedReadId.getValue()] <= 10) {
    //         continue;
    //     }

    //     uint64_t numberOfAlignmentsThatCanBeRemoved = verticesDegree[orientedReadId.getValue()] - 10;

    //     cout << "Vertex " << orientedReadId << " has degree " << out_degree(v, readGraph) << endl;
        

    //     BGL_FORALL_OUTEDGES(v, e, readGraph, ReadGraph4) {
    //         const uint32_t alignmentId = readGraph[e].alignmentId;
    //         const AlignmentData& alignment = alignmentData[alignmentId];
    //         const double errorRateRle = alignment.info.errorRateRle;
    //         const double markerCount = alignment.info.markerCount;

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         SHASTA_ASSERT(readId0 < readId1);
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         if (verticesDegree[B0.getValue()] <= 10) {
    //             continue;
    //         }

    //         neighborsMarkerCount.push_back(make_pair(target(e, readGraph), markerCount));
    //         neighborsMarkerCount2.push_back(make_pair(B0.getValue(), markerCount));
    //     }

    //     //cout << "Vertex " << orientedReadId << " has out degree " << out_degree(v, readGraph) << endl;
    //     //cout << "Vertex " << orientedReadId << " has in degree " << in_degree(v, readGraph) << endl;

    //     // Sort the neighbors by markerCount
    //     sort(neighborsMarkerCount.begin(), neighborsMarkerCount.end(), OrderPairsBySecondOnly<uint64_t, double>());
    //     sort(neighborsMarkerCount2.begin(), neighborsMarkerCount2.end(), OrderPairsBySecondOnly<uint64_t, double>());

    //     //cout << "Vertex " << orientedReadId << " has out degree " << out_degree(v, readGraph) << endl;
    //     //cout << "Vertex " << orientedReadId << " has in degree " << in_degree(v, readGraph) << endl;


    //     // Remove edges in the edgesToRemove vector
    //     for (uint32_t i = 0; i < numberOfAlignmentsThatCanBeRemoved; ++i) {
    //         edgesToRemove.push_back(make_pair(orientedReadId.getValue(), neighborsMarkerCount[i].first));
    //         //remove_edge(orientedReadId.getValue(), neighborsMarkerCount[i].first, readGraph);
    //         verticesDegree[orientedReadId.getValue()]--;
    //         verticesDegree[neighborsMarkerCount2[i].first]--;
    //         // Also decrease degree for reverse complement vertices
    //         OrientedReadId reverseOrientedReadId = orientedReadId;
    //         reverseOrientedReadId.flipStrand();
    //         OrientedReadId reverseNeighborId(OrientedReadId::fromValue(neighborsMarkerCount2[i].first));
    //         reverseNeighborId.flipStrand();
    //         verticesDegree[reverseOrientedReadId.getValue()]--;
    //         verticesDegree[reverseNeighborId.getValue()]--;
    //     }


    //     // cout << "Vertex " << orientedReadId << " has out degree " << out_degree(v, readGraph) << endl;
    //     //cout << "Vertex " << orientedReadId << " has in degree " << in_degree(v, readGraph) << endl;


    // }

    // cout << "Number of edges to remove: " << edgesToRemove.size() << endl;
    // // Remove them
    // for (uint32_t i = 0; i < edgesToRemove.size(); ++i) {
    //     //remove_edge(edgesToRemove[i], readGraph);
    //     remove_edge(edgesToRemove[i].first, edgesToRemove[i].second, readGraph);
    // }


    // // Initially, each alignment generates two edges.
    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

    //     // Record whether this alignment is used in the read graph.
    //     const bool keepThisAlignment = keepAlignment[alignmentId];
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // If this alignment is not used in the read graph, we are done.
    //     if(not keepThisAlignment) {
    //         continue;
    //     }

    //     // Get the OrientedReadIds.
    //     OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //     OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //     SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //     // Get the alignment data
    //     const ReadId readId0 = alignment.readIds[0];
    //     const ReadId readId1 = alignment.readIds[1];
    //     const bool isSameStrand = alignment.isSameStrand;
    //     SHASTA_ASSERT(readId0 < readId1);
    //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //     if (verticesDegree[A0.getValue()] < 10 or verticesDegree[B0.getValue()] < 10) {
    //         continue;
    //     }

        
        
    
    // }

    
    
    // // iterate over all alignments in readGraph and check if the alignment are in QAlignments
    
    // // Print the number of alignments in QAlignments
    // cout << "There are " << QAlignments.size() << " alignments in QAlignments." << endl;

    // // Now iterate over alignments in readGraph and check if they exist in QAlignments
    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        
    //     // Look for this alignmentId in QAlignments
    //     auto it = std::find_if(QAlignments.begin(), QAlignments.end(),
    //         [alignmentId](const pair<uint64_t, double>& p) { 
    //             return p.first == alignmentId; 
    //         });

    //     // If this alignment is not used in readGraph, and the alignment is in QAlignments, add the edge
    //     if(!alignmentData[alignmentId].info.isInReadGraph) {
    //         if(it != QAlignments.end()) {
    //             AlignmentData& alignment = alignmentData[alignmentId];
    //             // Get the OrientedReadIds.
    //             OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //             OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //             SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //             // Remove the edge.
    //             add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

    //             // Also remove the reverse complemented edge.
    //             orientedReadId0.flipStrand();
    //             orientedReadId1.flipStrand();
    //             add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

    //             alignment.info.isInReadGraph = 1;
    //             keepAlignment[alignmentId] = true;

    //             continue;
    //         }
    //     }
        
        
    //     // // If this alignment is used in readGraph, and the alignment is NOT in QAlignments, remove the edge
    //     // if(it == QAlignments.end()) {
    //     //     AlignmentData& alignment = alignmentData[alignmentId];
    //     //     // Get the OrientedReadIds.
    //     //     OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //     //     OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //     //     SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //     //     // Remove the edge.
    //     //     remove_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), readGraph);

    //     //     // Also remove the reverse complemented edge.
    //     //     orientedReadId0.flipStrand();
    //     //     orientedReadId1.flipStrand();
    //     //     remove_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), readGraph); 

    //     //     alignment.info.isInReadGraph = 0;
    //     //     keepAlignment[alignmentId] = false;
    //     // }
    // }


    
    







    // //*
    // //
    // // Now it is time to find the breaks in the readGraph
    // // In this case we add the alignments that did not pass the initial strict logQ filter
    // // involving reads that HAVE ALREADY BEEN USED (they are not isolated reads) AND are IN A DIFFERENT COMPONENT.
    // //
    // //*
    // vector<bool> endNodes(orientedReadCount, false);
    
    
    // bool weFoundAlignmentsInBreaksToAdd = true;

    // vector<bool> alignmentTableNotPassFilterAlreadyConsidered(alignmentTableNotPassFilter.size(), false);
    // const double QThresholdForAlignmentsInBreaks = 1e-10;
    // const double logQThresholdForAlignmentsInBreaks = log(QThresholdForAlignmentsInBreaks);

    // uint64_t allignmentsAdded = 0;
    // vector<vector<uint32_t>> betweenComponentsAlignments;
    // vector<pair<uint32_t, uint32_t>> componentsToUnite;


    // for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //     if((i % 100000) == 0) {
    //         cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //     }

    //     const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //     const double logQ = alignmentTableNotPassFilter[i].second;
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     const ReadId readId0 = alignment.readIds[0];
    //     const ReadId readId1 = alignment.readIds[1];
    //     const bool isSameStrand = alignment.isSameStrand;
    //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //     // Get the connected components that these oriented reads are in.
    //     const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //     const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //     const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //     const uint32_t b1 = disjointSets.find_set(B1.getValue());


    //     // Case 1: There is an alignment between two added reads
    //     // that passed the initial strict logQ filter and these reads are 
    //     // IN A DIFFERENT COMPONENT.
    //     // This case checks for direct connections between to nodes in a break
    //     // using alignments that did not pass the initial strict logQ filter.
    //     // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
    //     bool isInDifferentComponent = (a0 != b0) && (a0 != b1);
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && isInDifferentComponent && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
    //         // Find neighbors of the orientedReadId in the filtered readGraph.
    //         // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
    //         // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
    //         uint64_t maxDistance = 3;
    //         vector<OrientedReadId> readGraphNeighbors;
    //         readGraph.findNeighbors(A0, maxDistance, readGraphNeighbors);

    //         // Check if B0 is NOT in the neighbors
    //         if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //             cout << "Connecting components: b0 " << b0 << " and a0 " << a0 << endl;
    //             cout << "Components sizes are: b0 " << setSizes[b0] << " a0 " << setSizes[a0] << " b1 " << setSizes[b1] << " a1 " << setSizes[a1] << endl;
    //             cout << "Adding alignment of Read1: " << readId0 << " with the Read2: " << readId1 << endl;
    //             endNodes[A0.getValue()] = true;
    //             endNodes[B0.getValue()] = true;
    //             add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //             keepAlignment[alignmentId] = true;
    //             allignmentsAdded++;
    //             alignmentTableNotPassFilterAlreadyConsidered[i] = true;
    //             setSizes[a0] += setSizes[b0];
    //             setSizes[b0] = setSizes[a0];
    //             setSizes[a1] += setSizes[b1];
    //             setSizes[b1] = setSizes[a1];
    //             cout << "Components sizes are: b0 " << setSizes[b0] << " a0 " << setSizes[a0] << " b1 " << setSizes[b1] << " a1 " << setSizes[a1] << endl;
    //             componentsToUnite.push_back(make_pair(a0, b0));
    //             componentsToUnite.push_back(make_pair(a1, b1));
                
    //         }
    //         continue;
    //     }
    // }
    // for (auto& p : componentsToUnite) {
    //     disjointSets.union_set(p.first, p.second);
    // }

    // cout << "We added " << allignmentsAdded << " Alignments" << endl;







    // //*
    // //
    // // Now it is time to find the breaks in the readGraph and add the
    // // alignments between two added reads that passed the initial strict 
    // // logQ filter and these reads are IN A DIFFERENT COMPONENT.
    // //
    // //*
    // weFoundAlignmentsInBreaksToAdd = true;
    // while (weFoundAlignmentsInBreaksToAdd) {
    //     weFoundAlignmentsInBreaksToAdd = false;
    //     cout << "We added " << allignmentsAdded << " Alignments" << endl;
    //     for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //         if((i % 100000) == 0) {
    //             cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //         }

    //         const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //         const double logQ = alignmentTableNotPassFilter[i].second;
    //         const AlignmentData& alignment = alignmentData[alignmentId];

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // Case 1a: There is an alignment between two added reads
    //         // that passed the initial strict logQ filter and these reads are 
    //         // IN THE SAME CONNECTED COMPONENT.
    //         // This case checks for direct connections between to nodes in a break
    //         // using alignments that did not pass the initial strict logQ filter.
    //         // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
    //         if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && a0 == b0 && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
    //             // Find neighbors of the orientedReadId in the filtered readGraph.
    //             // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
    //             // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
                
    //             uint64_t maxDistance = 10;
    //             vector<OrientedReadId> readGraphNeighbors;
    //             readGraphAllAlignments.findNeighborsEarlyStopWhenReachSameComponentNode(A0, disjointSets, maxDistance, readGraphNeighbors);

    //             // Get the last node in the neighbors, which contains the first node in the same component as A0
    //             OrientedReadId lastNode = readGraphNeighbors.back();

    //             // Verify it is in the same component as A0
    //             if(disjointSets.find_set(lastNode.getValue()) == a0) {
    //                 SHASTA_ASSERT(a1 == b0);
    //             }

    //             // Check if it is in the same component as A0
    //             if(disjointSets.find_set(lastNode.getValue()) != a0) {
    //                 // We did not found a node in the same component as A0
    //                 // Tha means that this is probably not a break situation
    //                 // Or that we did not search further enough
    //                 continue;

    //             }

    //             // Check if it is in the same component as A0
                
    //             endNodes[A0.getValue()] = true;
    //             endNodes[B0.getValue()] = true;
    //             add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //             keepAlignment[alignmentId] = true;
    //             weFoundAlignmentsInBreaksToAdd = true;
    //             allignmentsAdded++;
    //             alignmentTableNotPassFilterAlreadyConsidered[i] = true;
                

    //             // 


    //             // Check if B0 is NOT in the neighbors
    //             if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //                 endNodes[A0.getValue()] = true;
    //                 endNodes[B0.getValue()] = true;
    //                 add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //                 keepAlignment[alignmentId] = true;
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //                 allignmentsAdded++;
    //                 alignmentTableNotPassFilterAlreadyConsidered[i] = true;
    //             }
    //             continue;
    //         }
    //     }
    // }






    //*
    //
    // Now it is time to find the breaks in the readGraph and add the
    // alignments between two added reads that passed the initial strict 
    // logQ filter and these reads are IN A DIFFERENT COMPONENT.
    //
    //*
    //vector<bool> endNodes(orientedReadCount, false);
    //vector<bool> alignmentTableNotPassFilterAlreadyConsidered(alignmentTableNotPassFilter.size(), false);
    //bool weFoundAlignmentsInBreaksToAdd = true;

    //const double QThresholdForAlignmentsInBreaks = 1e-10;
    //const double logQThresholdForAlignmentsInBreaks = log(QThresholdForAlignmentsInBreaks);

    //uint64_t allignmentsAdded = 0;

    // while (weFoundAlignmentsInBreaksToAdd) {
    //     weFoundAlignmentsInBreaksToAdd = false;
    //     cout << "We added " << allignmentsAdded << " Alignments" << endl;
    //     for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //         if((i % 100000) == 0) {
    //             cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //         }

    //         const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //         const double logQ = alignmentTableNotPassFilter[i].second;
    //         const AlignmentData& alignment = alignmentData[alignmentId];

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

            // // Case 1a: There is an alignment between two added reads
            // // that passed the initial strict logQ filter and these reads are 
            // // IN THE SAME CONNECTED COMPONENT.
            // // This case checks for direct connections between to nodes in a break
            // // using alignments that did not pass the initial strict logQ filter.
            // // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
            // if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && a0 == b0 && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
            //     // Find neighbors of the orientedReadId in the filtered readGraph.
            //     // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
            //     // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
            //     uint64_t maxDistance = 5;
            //     vector<OrientedReadId> readGraphNeighbors;
            //     readGraph.findNeighbors(A0, maxDistance, readGraphNeighbors);

            //     // Check if B0 is NOT in the neighbors
            //     if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
            //         endNodes[A0.getValue()] = true;
            //         endNodes[B0.getValue()] = true;
            //         add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
            //         keepAlignment[alignmentId] = true;
            //         weFoundAlignmentsInBreaksToAdd = true;
            //         allignmentsAdded++;
            //         alignmentTableNotPassFilterAlreadyConsidered[i] = true;
            //     }
            //     continue;
            // }

            

            // // Case 2: There is an alignment between an originaly used read that was in an alignment
            // // that passed the initial strict logQ filter, and an isolated read that was not used.
            // // if(readUsed[alignment.readIds[0]] && !readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks) {
            // if(readUsed[alignment.readIds[0]] && !readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
            //     // Find neighbors of the orientedReadId B0 (the isolated one) in the allAlignments readGraph.
            //     // Next, check which of these neighbors are in the same connected component as A0.
            //     // Next, check if these neighbors of B0 in the same connected component of A0 are reachable 
            //     // in a set distance from A0 in the filtered readGraph.
            //     // If some of them are not, we can add the alignment because we are in a break.
            //     // TODO: The maxDistances need to be adjusted probably.
            //     uint64_t maxDistance = 8;
            //     vector<OrientedReadId> readGraphAllAlignmentsNeighbors;
            //     readGraphAllAlignments.findNeighbors(B0, maxDistance, readGraphAllAlignmentsNeighbors);
                
            //     maxDistance = 12;
            //     vector<OrientedReadId> readGraphNeighbors;
            //     readGraph.findNeighbors(A0, maxDistance, readGraphNeighbors);

            //     // Check which of these neighbors of B0 are in the same connected component as A0.
            //     vector<OrientedReadId> neighborsInSameComponent;
            //     for(const OrientedReadId& neighbor: readGraphAllAlignmentsNeighbors) {
            //         if(disjointSets.find_set(neighbor.getValue()) == a0) {
            //             neighborsInSameComponent.push_back(neighbor);
            //         }
            //     }

            //     // Check if these neighbors of B0 in the same component of A0 are reachable from A0 in the filtered readGraph.
            //     // If any of them are not, we can add the alignment.
            //     bool weCanAddAlignment = false;
            //     for(const OrientedReadId& neighbor: neighborsInSameComponent) {
            //         if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), neighbor) == readGraphNeighbors.end()) {
            //             weCanAddAlignment = true;
            //             break;
            //         }
            //     }

            //     if(weCanAddAlignment) {
            //         endNodes[A0.getValue()] = true;
            //         endNodes[B0.getValue()] = true;
            //         add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
            //         keepAlignment[alignmentId] = true;
            //         readUsed[alignment.readIds[1]] = true;
            //         // Update disjoint sets
            //         disjointSets.union_set(a0, b0);
            //         disjointSets.union_set(a1, b1);
            //         weFoundAlignmentsInBreaksToAdd = true;
            //         allignmentsAdded++;
            //         alignmentTableNotPassFilterAlreadyConsidered[i] = true;
            //     }

            //     continue;
            // }
            
            // // Case 3: There is an alignment between an isolated read that was not used and
            // // an originaly used read that was in an alignment that passed the initial strict logQ filter.
            // // if(!readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks) {
            // if(!readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
            //     // Find neighbors of the orientedReadId A0 (the isolated one) in the allAlignments readGraph.
            //     // Next, check which of these neighbors are in the same connected component as B0.
            //     // Next, check if these neighbors of A0 in the same connected component of B0 are reachable 
            //     // in a set distance from B0 in the filtered readGraph.
            //     // If some of them are not, we can add the alignment because we are in a break.

            //     uint64_t maxDistance = 8;
            //     vector<OrientedReadId> readGraphAllAlignmentsNeighbors;
            //     readGraphAllAlignments.findNeighbors(A0, maxDistance, readGraphAllAlignmentsNeighbors);

            //     maxDistance = 12;
            //     vector<OrientedReadId> readGraphNeighbors;
            //     readGraph.findNeighbors(B0, maxDistance, readGraphNeighbors);

            //     // Check which of these neighbors of A0 are in the same connected component as B0.
            //     vector<OrientedReadId> neighborsInSameComponent;
            //     for(const OrientedReadId& neighbor: readGraphAllAlignmentsNeighbors) {
            //         if(disjointSets.find_set(neighbor.getValue()) == b0) {
            //             neighborsInSameComponent.push_back(neighbor);
            //         }
            //     }

            //     // Check if these neighbors of A0 in the same component of B0 are reachable from B0 in the filtered readGraph.
            //     // If any of them are not, we can add the alignment.
            //     bool weCanAddAlignment = false;
            //     for(const OrientedReadId& neighbor: neighborsInSameComponent) {
            //         if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), neighbor) == readGraphNeighbors.end()) {
            //             weCanAddAlignment = true;
            //             break;
            //         }
            //     }

            //     if(weCanAddAlignment) {
            //         endNodes[A0.getValue()] = true;
            //         endNodes[B0.getValue()] = true;
            //         add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
            //         keepAlignment[alignmentId] = true;
            //         readUsed[alignment.readIds[0]] = true;
            //         // Update disjoint sets
            //         disjointSets.union_set(a0, b0);
            //         disjointSets.union_set(a1, b1);
            //         weFoundAlignmentsInBreaksToAdd = true;
            //         allignmentsAdded++;
            //         alignmentTableNotPassFilterAlreadyConsidered[i] = true;
            //     }

            //     continue;
            // }
    //     }
    // }


    // Remove previously created graphs
    //readGraphAllAlignments.remove();












    // // In this step now we will add additional alignments that were filterd by
    // // the Q filter and that involve isolated oriented reads that are added (readAdded) 
    // // and they belong to the same connected component.
    // const double QThresholdForFinalRoundAddedAlignments = 1e-4;
    // const double logQThresholdForFinalRoundAddedAlignments = log(QThresholdForFinalRoundAddedAlignments);
    // bool weFoundAlignmentsInBreaksToAdd = true;
    // // while (weFoundAlignmentsInBreaksToAdd) {
    // //     weFoundAlignmentsInBreaksToAdd = false;
    // for(size_t i=0; i<alignmentTableNotPassFilterNotIncludedAlignments.size(); i++) {
    //     if((i % 10000) == 0) {
    //         cout << timestamp << i << "/" << alignmentTableNotPassFilterNotIncludedAlignments.size() << endl;
    //     }

    //     const uint64_t alignmentId = alignmentTableNotPassFilterNotIncludedAlignments[i].first;
    //     const double logQ = alignmentTableNotPassFilterNotIncludedAlignments[i].second;
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // Case 1: There is an alignment between two added reads
    //     if(readAdded[alignment.readIds[0]] && readAdded[alignment.readIds[1]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 //cout << "Case 1 happened" << endl;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //             }
    //         }
    //         continue;
    //     }

    //     // Case 2a: There is an alignment between an added read and a originaly used read
    //     if(readAdded[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && !readAdded[alignment.readIds[1]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 //cout << "Case 2a happened" << endl;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //             }
    //         }
    //         continue;
    //     }

    //     // Case 2b: There is an alignment between an added read and a originaly used read
    //     if(readAdded[alignment.readIds[1]] && readUsed[alignment.readIds[0]] && !readAdded[alignment.readIds[0]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 //cout << "Case 2b happened" << endl;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //             }
    //         }
    //         continue;
    //     }



    //     // Case 3: There is an alignment between two originaly used read
    //     // We want to add only alignments that connect breaks in the read graph.
    //     // To do that we will check if the two oriented reads are in the same connected component
    //     // and if their distance with BFS is greater than 1 given that we found an alignment that connects them
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 cout << "Case 3a happened" << endl;

    //                 // Find neighbors of the orientedReadId
    //                 uint64_t maxDistance = 3;

    //                 vector<OrientedReadId> readGraphNeighbors;
    //                 readGraphAllAlignments.findNeighbors(A0, maxDistance, readGraphNeighbors);

    //                 // Check if B0 is NOT in the neighbors
    //                 if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //                     cout << "Case 3b happened" << endl;
    //                     keepAlignment[alignmentId] = true;
    //                     // Update disjoint sets
    //                     disjointSets.union_set(a0, b0);
    //                     disjointSets.union_set(a1, b1);
    //                     weFoundAlignmentsInBreaksToAdd = true;
    //                 }
                    
    //             }
                
    //         }
            
    //     }
    // }



    
    


    // // In this step now we will add additional alignments that were filterd by
    // // the Q filter and that involve isolated oriented reads that are added (readAdded) 
    // // and they belong to the same connected component.
    // const double QThresholdForFinalRoundAddedAlignments = 1e-5;
    // const double logQThresholdForFinalRoundAddedAlignments = log(QThresholdForFinalRoundAddedAlignments);
    // for(size_t i=0; i<alignmentTableNotPassFilterNotIncludedAlignments.size(); i++) {
    //     const uint64_t alignmentId = alignmentTableNotPassFilterNotIncludedAlignments[i].first;
    //     const double logQ = alignmentTableNotPassFilterNotIncludedAlignments[i].second;
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // check if the reads were used in alignments
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
    //         // if(readAdded[alignment.readIds[0]] && readAdded[alignment.readIds[1]]) {
    //         //     const ReadId readId0 = alignment.readIds[0];
    //         //     const ReadId readId1 = alignment.readIds[1];
    //         //     const bool isSameStrand = alignment.isSameStrand;
    //         //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         //     // Get the connected components that these oriented reads are in.
    //         //     const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         //     const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         //     const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         //     const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         //     // check if the reads are in the same connected component
    //         //     if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //         //         keepAlignment[alignmentId] = true;
    //         //         // Update disjoint sets
    //         //         disjointSets.union_set(a0, b0);
    //         //         disjointSets.union_set(a1, b1);
    //         //     }
    //         // }
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //             }
                
    //         }
            
    //     }
    // }

    


    // // We iterate over each conected component we have constructed
    // // using strict disjointSets and we try to detect breaks.
    // // If there is a break in the readGraph and there is not a break in the AllAlignments graph
    // // we need to add the alignments that connect the two oriented reads to the readGraph.
    // for(auto it=componentMap.begin(); it!=componentMap.end(); ++it) {
    //     const vector<OrientedReadId>& orientedReadIds = it->second;
    //     for(const OrientedReadId& orientedReadId: orientedReadIds) {
    //         // Find neighbors of the orientedReadId
    //         uint64_t maxDistance = 5;

    //         vector<OrientedReadId> readGraphNeighbors;
    //         readGraph.findNeighbors(orientedReadId, maxDistance, readGraphNeighbors);

    //         vector<OrientedReadId> readGraphAllAlignmentsNeighbors;
    //         readGraphAllAlignments.findNeighbors(orientedReadId, maxDistance, readGraphAllAlignmentsNeighbors);

    //         // Find the neighbors that exist in the AllAlignments graph but not in the readGraph and also
    //         // they are in the disjointSets of the orientedReadId
    //         for(const OrientedReadId& orientedReadId1: readGraphAllAlignmentsNeighbors) {
    //             if(disjointSets.find_set(orientedReadId1.getValue()) == disjointSets.find_set(orientedReadId.getValue())) {

    //                 // Check if the distance between the OrientedReadIs are significant different in the readGraph and AllAlignments graph
    //                 // Find the shortest path between the two oriented reads
    //                 vector<uint32_t> distanceReadGraph(2*readCount, ReadGraph::infiniteDistance);
    //                 vector<OrientedReadId> reachedVerticesReadGraph;
    //                 vector<uint32_t> parentEdgesReadGraph(2*readCount);
    //                 vector<uint32_t> shortestPathReadGraph;
    //                 readGraph.computeShortPath(orientedReadId, orientedReadId1,
    //                     maxDistance, shortestPathReadGraph,
    //                     distanceReadGraph, reachedVerticesReadGraph, parentEdgesReadGraph);

    //                 vector<uint32_t> distanceReadGraphAllAlignments(2*readCount, ReadGraph::infiniteDistance);
    //                 vector<OrientedReadId> reachedVerticesReadGraphAllAlignments;
    //                 vector<uint32_t> parentEdgesReadGraphAllAlignments(2*readCount);
    //                 vector<uint32_t> shortestPathReadGraphAllAlignments;
    //                 readGraphAllAlignments.computeShortPath(orientedReadId, orientedReadId1,
    //                     maxDistance, shortestPathReadGraphAllAlignments,
    //                     distanceReadGraphAllAlignments, reachedVerticesReadGraphAllAlignments, parentEdgesReadGraphAllAlignments);

    //                 // If the distance between the OrientedReadIs are about the same then we can skip
    //                 // It means that there is not really a break in the readGraph. If there was a break
    //                 // we would have a significant difference in the distance between the two OrientedReadIs
    //                 // because a way bigger path would be needed to connect them.
    //                 uint64_t maxDistance = 5;
    //                 if (distanceReadGraph[orientedReadId1.getValue()] != ReadGraph::infiniteDistance && 
    //                     distanceReadGraphAllAlignments[orientedReadId1.getValue()] != ReadGraph::infiniteDistance) {
    //                     if(std::abs(int(distanceReadGraph[orientedReadId1.getValue()] - distanceReadGraphAllAlignments[orientedReadId1.getValue()])) < maxDistance) {
    //                         continue;
    //                     }

    //                 }

    //                 // Else, we need to add the alignments that connect the two oriented reads
    //                 // and add them to the readGraph if they do not break the disjoint sets introducing cross-strand edges
    //                 // Find all alignments that connect orientedReadId and orientedReadId1 even if they do not connect them directly
    //                 vector<uint32_t> alignmentsBetweenReads;

                    


    //                 if(std::find(readGraphNeighbors.begin(), readGraphNeighbors.end(), orientedReadId1) == readGraphNeighbors.end()) {
    //                     // Find all alignments that connect orientedReadId and orientedReadId1
    //                     vector<uint32_t> alignmentsBetweenReads;
    //                     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //                         const AlignmentData& alignment = alignmentData[alignmentId];
    //                         if(alignment.info.isInReadGraph) {
    //                             continue;  // Skip alignments already in read graph
    //                         }
                            
    //                         const ReadId readId0 = alignment.readIds[0];
    //                         const ReadId readId1 = alignment.readIds[1];
    //                         const bool isSameStrand = alignment.isSameStrand;
                            
    //                         OrientedReadId readId0Oriented(readId0, 0);
    //                         OrientedReadId readId1Oriented(readId1, isSameStrand ? 0 : 1);

    //                         // Check if this alignment connects our two reads
    //                         if((readId0Oriented == orientedReadId && readId1Oriented == orientedReadId1) ||
    //                            (readId0Oriented == orientedReadId1 && readId1Oriented == orientedReadId)) {
    //                             alignmentsBetweenReads.push_back(alignmentId);
    //                         }
    //                     }

    //                     // Add these alignments to the read graph
    //                     for(uint32_t alignmentId : alignmentsBetweenReads) {
    //                         AlignmentData& alignment = alignmentData[alignmentId];
    //                         alignment.info.isInReadGraph = true;
    //                         readGraph.addEdge(alignment.readIds[0], alignment.readIds[1], 
    //                                          alignment.isSameStrand, alignmentId);
    //                     }
    //                 }
    //             }
    //         }
            


    //         // iterate over the neighbors
    //         for(const OrientedReadId& orientedReadId1: readGraphNeighbors) {
    //             // Find the shortest path between the two oriented reads
    //             vector<uint32_t> distanceReadGraph(2*readCount, ReadGraph::infiniteDistance);
    //             vector<OrientedReadId> reachedVerticesReadGraph;
    //             vector<uint32_t> parentEdgesReadGraph(2*readCount);
    //             vector<uint32_t> shortestPathReadGraph;
    //             readGraph.computeShortPath(orientedReadId, orientedReadId1,
    //                 maxDistance, shortestPathReadGraph,
    //                 distanceReadGraph, reachedVerticesReadGraph, parentEdgesReadGraph);

    //             vector<uint32_t> distanceReadGraphAllAlignments(2*readCount, ReadGraph::infiniteDistance);
    //             vector<OrientedReadId> reachedVerticesReadGraphAllAlignments;
    //             vector<uint32_t> parentEdgesReadGraphAllAlignments(2*readCount);
    //             vector<uint32_t> shortestPathReadGraphAllAlignments;
    //             readGraph.computeShortPath(orientedReadId, orientedReadId1,
    //                 maxDistance, shortestPathReadGraphAllAlignments,
    //                 distanceReadGraphAllAlignments, reachedVerticesReadGraphAllAlignments, parentEdgesReadGraphAllAlignments);

    //             // If the sistance between the OrientedReadIs are about the same then we can skip
    //             if (distanceReadGraph[orientedReadId1.getValue()] != ReadGraph::infiniteDistance && 
    //                 distanceReadGraphAllAlignments[orientedReadId1.getValue()] != ReadGraph::infiniteDistance) {
    //                 if(std::abs(int(distanceReadGraph[orientedReadId1.getValue()] - distanceReadGraphAllAlignments[orientedReadId1.getValue()])) < 5) {
    //                     continue;
    //                 }

    //             }

    //             // Iterate over all reach vertices in the AllAlignments graph and check if they
    //             // are in the same set of the disjointSets
    //             for(const OrientedReadId& reachedVertex: reachedVerticesReadGraphAllAlignments) {
    //                 if(disjointSets.find_set(reachedVertex.getValue()) == disjointSets.find_set(orientedReadId.getValue())) {
    //                     // now check if it is reachable from the readGraph
    //                     if(distanceReadGraph[reachedVertex.getValue()] == ReadGraph::infiniteDistance) {
    //                         // It is not reachable from the readGraph
    //                         // This means we need to find all the alignments that connect the two oriented reads
    //                         // and add them to the readGraph.
    //                         // Find all alignments that connect orientedReadId and reachedVertex
    //                         vector<uint32_t> alignmentsBetweenReads;
    //                         for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //                             const AlignmentData& alignment = alignmentData[alignmentId];
    //                             if(alignment.info.isInReadGraph) {
    //                                 continue;  // Skip alignments already in read graph
    //                             }
                                
    //                             const ReadId readId0 = alignment.readIds[0];
    //                             const ReadId readId1 = alignment.readIds[1];
    //                             const bool isSameStrand = alignment.isSameStrand;
                                
    //                             OrientedReadId readId0Oriented(readId0, 0);
    //                             OrientedReadId readId1Oriented(readId1, isSameStrand ? 0 : 1);

    //                             // Check if this alignment connects our two reads
    //                             if((readId0Oriented == orientedReadId && readId1Oriented == reachedVertex) ||
    //                                (readId0Oriented == reachedVertex && readId1Oriented == orientedReadId)) {
    //                                 alignmentsBetweenReads.push_back(alignmentId);
    //                             }
    //                         }

    //                         // Add these alignments to the read graph
    //                         for(uint32_t alignmentId : alignmentsBetweenReads) {
    //                             AlignmentData& alignment = alignmentData[alignmentId];
    //                             alignment.info.isInReadGraph = true;
    //                             readGraph.addEdge(alignment.readIds[0], alignment.readIds[1], 
    //                                              alignment.isSameStrand, alignmentId);
    //                         }


    //                     }
                        
    //                 }
    //             }

    //             // Check if the distance is the same
    //             if(distanceReadGraph[orientedReadId1.getValue()] == ReadGraph::infiniteDistance && distanceReadGraphAllAlignments[orientedReadId1.getValue()]) {
    //             }

                
    //         }
            


    //     }

    // }













    // // Process notIncludedReadsToAlignments in order of increasing Q
    // // but now use the already constructed disjointSets as a "good starting point"
    // uint64_t crossStrandEdgeCountR2 = 0;
    // for(auto it=notIncludedReadsToAlignments.begin(); it!=notIncludedReadsToAlignments.end(); ++it) {
    //     const pair<uint64_t, double>& p = *it;
    //     const uint64_t alignmentId = p.first;

    //     // Get the alignment data
    //     AlignmentData& alignment = alignmentData[alignmentId];
    //     const ReadId readId0 = alignment.readIds[0];
    //     const ReadId readId1 = alignment.readIds[1];
    //     const bool isSameStrand = alignment.isSameStrand;
    //     SHASTA_ASSERT(readId0 < readId1);
    //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //     SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
    //     SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
    //     SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
    //     SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

    //     // Get the connected components that these oriented reads are in.
    //     const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //     const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //     const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //     const uint32_t b1 = disjointSets.find_set(B1.getValue());


    //     // If the alignment breaks strand separation, it is skipped.
    //     // If A0 and B1 are in the same connected component,
    //     // A1 and B0 also must be in the same connected component.
    //     // Adding this pair of edges would create a self-complementary
    //     // connected component containing A0, B0, A1, and B1,
    //     // and to ensure strand separation we don't want to do that.
    //     // So we mark these edges as cross-strand edges
    //     // and don't use them to update the disjoint set data structure.
    //     if(a0 == b1) {
    //         SHASTA_ASSERT(a1 == b0);
    //         crossStrandEdgeCountR2 += 2;
    //         continue;
    //     }

    //     // If both vertices of the potential edge have at least the required minimum number 
    //     // of neighbors, the alignment is also skipped. 
    //     const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //     const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //     if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //         // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
    //         continue;
    //     }

    //     // Add the alignment to the read graph.
    //     keepAlignment[alignmentId] = true;
    //     alignment.info.isInReadGraph = 1;

    //     // Update vertex degrees
    //     verticesDegree[A0.getValue()]++;
    //     verticesDegree[B0.getValue()]++;
    //     verticesDegree[A1.getValue()]++;
    //     verticesDegree[B1.getValue()]++;
        

    //     // Update disjoint sets
    //     disjointSets.union_set(a0, b0);
    //     disjointSets.union_set(a1, b1);
    // }

    // // Verify that for any read the two oriented reads are in distinct
    // // connected components.
    // for(ReadId readId=0; readId<readCount; readId++) {
    //     const OrientedReadId orientedReadId0(readId, 0);
    //     const OrientedReadId orientedReadId1(readId, 1);
    //     SHASTA_ASSERT(
    //         disjointSets.find_set(orientedReadId0.getValue()) !=
    //         disjointSets.find_set(orientedReadId1.getValue())
    //     );
    // }


    // // Print how many alignments were kept
    // const size_t keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
    // cout << "Adding isolated reads step: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;

    

    cout << timestamp << "Done processing alignments." << endl;


    

    cout << timestamp << "createReadGraph4 with strand separation ends." << endl;

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << num_edges(readGraph) << " total in round 1." << endl;
    
    // cout << "Strand separation flagged " << crossStrandEdgeCountR2 <<
    //     " read graph edges out of " << readGraph.edges.size() << " total in round 2." << endl;


    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // Stop writing when we reach connected components
        // consisting of a single isolated read.
        if(component.size() == 1) {
            break;
        }

        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        // The above process of strand separation should have removed
        // all self-complementary components.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        SHASTA_ASSERT(not isSelfComplementary);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // For Mode 2 and Mode 3 assembly, we will only assemble one connected component
    // of each pair. In each pair, we choose the component in the pair
    // that has the lowest numbered read on strand 0.
    // Then, for each read we store in its ReadFlags the strand
    // that the read appears in in this component.
    // That flag will be used in Mode 2 assembly to
    // select portions of the marker graph that should be assembled.
    uint64_t n = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If the lowest numbered read is on strand 1, this is not one of
        // the connected components we want to use.
        if(component.front().getStrand() == 1) {
            continue;
        }

        // Store the strand for each read in this component.
        for(const OrientedReadId orientedReadId: component) {
            reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
        }
        n += component.size();
    }
    SHASTA_ASSERT(n == readCount);


}
























void Assembler::createReadGraph4(
    uint32_t maxAlignmentCount)
{
    const bool debug = false;

    bool assemblyIsAvailable = false;
    try {
        accessMode3Assembler();
        SHASTA_ASSERT(mode3Assembler);
        const mode3::Anchors& anchors = mode3Assembler->anchors();
        SHASTA_ASSERT(anchors.journeys.isOpen());
        SHASTA_ASSERT(anchors.journeys.size() == markers.size());
        assemblyIsAvailable= true;
    } catch (...) {
    }


    if(assemblyIsAvailable) { 
        // Skip if assembly is available
    } else {
        // QRle threshold to use an alignment in the read graph.
        const double minQRle = 26;

        const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

        cout << timestamp << "createReadGraph4 begins, minQRle " << minQRle << endl;

        // Get the total number of stored alignments.
        const uint64_t alignmentCount = alignmentData.size();
        SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

        // Flag all alignments as not to be kept.
        vector<bool> keepAlignment(alignmentCount, false);

        // Loop over all alignments.
        for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
            if((alignmentId % 1000) == 0) {
                cout << timestamp << alignmentId << "/" << alignmentCount << endl;
            }

            // Get information for this alignment.
            AlignmentData& thisAlignmentData = alignmentData[alignmentId];
            // keepAlignment[alignmentId] = true;
            // thisAlignmentData.info.isInReadGraph = 1;

            // The alignment is stored as an alignment between readId0 on strand 0
            // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
            // The reverse complement alignment also exists, but is not stored explicitly.
            const ReadId readId0 = thisAlignmentData.readIds[0];
            const ReadId readId1 = thisAlignmentData.readIds[1];
            const bool isSameStrand = thisAlignmentData.isSameStrand;
            SHASTA_ASSERT(readId0 < readId1);
            const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
            const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

            const uint32_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
            const uint32_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
            const double L = (range0 + range1)/2;
            const uint64_t n = thisAlignmentData.info.mismatchCountRle;
            const double errorRateRle = double(n)/(2.0*L);;

            // If the RLE Q is large enough, flag thus alignment as to be kept.
            if((errorRateRle <= maxErrorRateRle)) {
                keepAlignment[alignmentId] = true;
                thisAlignmentData.info.isInReadGraph = 1;
            } else {
                keepAlignment[alignmentId] = false;
                thisAlignmentData.info.isInReadGraph = 0;
            }

        }

        cout << timestamp << "Done processing alignments." << endl;

        const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
        cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

        // Create the read graph using the alignments we selected.
        createReadGraphUsingSelectedAlignments(keepAlignment);

        cout << timestamp << "Initial createReadGraph4 ends." << endl;

        // Remove bridges from the read graph.
        removeReadGraphBridges(5);

        cout << timestamp << "createReadGraph4 ends." << endl;
    }

    

}




// Strict strand separation in the read graph based on using Poisson CDF.
// This guarantees that the read graph contains no self-complementary
// connected components.
// In other words, for any ReadId x, the two oriented reads
// x-0 and x-1 are guaranteed to be in distinct components of the
// read graph.
void Assembler::flagCrossStrandReadGraphEdges4()
{
    // Each alignment used in the read graph generates a pair of
    // consecutively numbered edges in the read graph
    // which are the reverse complement of each other.
    SHASTA_ASSERT((readGraph.edges.size() % 2) == 0);

    // Below, each pair is identified by the (even) index of
    // the first edge in the pair.

    // For each number of aligned markers alignedMarkerCount,
    // Gather in edgeTable[alignedMarkerCount] pairs
    // with that number of aligned markers.
    vector< pair<uint64_t, double> > edgeTable;

    // To loop over pairs of edges, increment by 2.
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId+=2) {
        const ReadGraphEdge& edge = readGraph.edges[edgeId];
        SHASTA_ASSERT(not edge.crossesStrands);

        // Skip edges flagged as inconsistent.
        if(edge.hasInconsistentAlignment) {
            continue;
        }

        const uint64_t alignmentId = edge.alignmentId;
        const AlignmentData& alignment = alignmentData[alignmentId];

        // Skip edges involving reads classified as chimeric.
        if(getReads().getFlags(alignment.readIds[0]).isChimeric) {
            continue;
        }
        if(getReads().getFlags(alignment.readIds[1]).isChimeric) {
            continue;
        }

        // Sanity check.
        SHASTA_ASSERT(alignmentData[alignmentId].info.isInReadGraph);

        // Check that the next edge is the reverse complement of
        // this edge.
        {
            const ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
            SHASTA_ASSERT(not nextEdge.crossesStrands);
            array<OrientedReadId, 2> nextEdgeOrientedReadIds = nextEdge.orientedReadIds;
            nextEdgeOrientedReadIds[0].flipStrand();
            nextEdgeOrientedReadIds[1].flipStrand();
            SHASTA_ASSERT(nextEdgeOrientedReadIds == edge.orientedReadIds);
            SHASTA_ASSERT(nextEdge.alignmentId == alignmentId);
        }


        // // Store this pair of edges in our edgeTable.
        // const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[0], 0, markers);
        // const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[1], 1, markers);
        // const double L = (range0 + range1) / 2;
        // const uint64_t n = alignment.info.mismatchCountRle;
        // const double m = L * 0.0001;

        // // Calculate Q using Poisson CDF
        // // Helper function to calculate gamma_q(k+1, λ)
        // auto calculate_gamma_q = [](uint64_t k, double lambda) -> double
        // {
        //     // gamma_q(a,z) is the normalized upper incomplete gamma function
        //     // For Poisson CDF, we need gamma_q(k+1, λ)
        //     double a = static_cast<double>(k + 1);
        //     return boost::math::gamma_q(a, lambda);
        // };

        // const double poissonCDF = calculate_gamma_q(n, m);


        // edgeTable.push_back(make_pair(edgeId, poissonCDF));

        // Store this pair of edges in our edgeTable.
        const double errorRateRle = alignment.info.errorRateRle;
        const double L = (alignment.info.range(0) + alignment.info.range(1)) / 0.04;
        const double m = L * 0.0001;
        const uint64_t n = static_cast<uint64_t>(std::round(errorRateRle * L));
        // Calculate Q using Poisson CDF
        // Helper function to calculate gamma_q(k+1, λ)
        auto calculate_gamma_q = [](uint64_t k, double lambda) -> double
        {
            // gamma_q(a,z) is the normalized upper incomplete gamma function
            // For Poisson CDF, we need gamma_q(k+1, λ)
            double a = static_cast<double>(k + 1);
            return boost::math::gamma_q(a, lambda);
        };
        const double poissonCDF = calculate_gamma_q(n, m);

        edgeTable.push_back(make_pair(edgeId, poissonCDF));
        
    }

    // Sort by increasing Q
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<uint64_t, double>());

    // Print out top 100 alignments by Q value
    cout << "Top 100 alignments by Q value:" << endl;
    cout << "EdgeId\tQ\tMarkerCount\tErrorRateRleNew\tErrorRateRleOld\tL\tn" << endl;
    for(size_t i = 0; i < min(size_t(100), edgeTable.size()); i++) {
        const auto& p = edgeTable[i];
        const uint64_t edgeId = p.first;
        const double q = p.second;
        const AlignmentData& alignment = alignmentData[readGraph.edges[edgeId].alignmentId];
        const uint64_t markerCount = alignment.info.markerCount;
        const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[0], 0, markers);
        const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[1], 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = alignment.info.mismatchCountRle;
        const double errorRateRleNew = double(n)/(2.0*L);;
        const double errorRateRleOld = alignment.info.errorRateRle;
        const double m = L * 0.0001;
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        cout << edgeId << "(" << readId0 << "," << readId1 << ")" << "\t" << q << "\t" << markerCount << "\t" << errorRateRleNew << "\t" << errorRateRleOld << "\t" << L << "\t" << n << endl;
    }


    // Create and initialize the disjoint sets data structure needed below.
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    SHASTA_ASSERT(readGraph.connectivity.size() == orientedReadCount);
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }


    // Now process our edge table in order of decreasing alignedMarkerCount.
    // Each pair consists of two edges:
    // A0--B0
    // A1--B1
    // where A1 is the reverse complement of A0 and B1 is the reverse complement of B0.
    // Maintain a disjoint set data structure for read graph vertices.
    // For each pair, compute the current connected components
    // for A0 B0 A1 B1 from the disjoint set data structure, call them a0 b0 a1 b1.
    // If a1==b0 (in which case also b1==a0), adding these edges
    // would create a self-complementary connected component,
    // so we mark them as cross-strand edges and don't add them to the
    // disjoint set data structure.
    // Otherwise, the two edges are added to the disjoint set data structure.
    uint64_t crossStrandEdgeCount = 0;
    for(auto it=edgeTable.begin(); it!=edgeTable.end(); ++it) {
        const pair<uint64_t, double>& p = *it;
        const uint64_t edgeId = p.first;
        ReadGraphEdge& edge = readGraph.edges[edgeId];
        ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
        const OrientedReadId A0 = edge.orientedReadIds[0];
        const OrientedReadId B0 = edge.orientedReadIds[1];
        const OrientedReadId A1 = nextEdge.orientedReadIds[0];
        const OrientedReadId B1 = nextEdge.orientedReadIds[1];

        SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
        SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
        SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
        SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

        // Get the connected components that these oriented reads are in.
        const uint32_t a0 = disjointSets.find_set(A0.getValue());
        const uint32_t b0 = disjointSets.find_set(B0.getValue());
        const uint32_t a1 = disjointSets.find_set(A1.getValue());
        const uint32_t b1 = disjointSets.find_set(B1.getValue());

        // If A0 and B0 are in the same connected component,
        // A1 and B1 also must be in the same connected component.
        // There is nothing to do as this pair of edges
        // does not affect the connected components.
        if (a0 == b0) {
            SHASTA_ASSERT(a1 == b1);
            continue;
        }

        // If A0 and B1 are in the same connected component,
        // A1 and B0 also must be in the same connected component.
        // Adding this pair of edges would create a self-complementary
        // connected component containing A0, B0, A1, and B1,
        // and to ensure strand separation we don't want to do that.
        // So we mark these edges as cross-strand edges
        // and don't use them to update the disjoint set data structure.
        if(a0 == b1) {
            SHASTA_ASSERT(a1 == b0);
            edge.crossesStrands = 1;
            nextEdge.crossesStrands = 1;
            alignmentData[edge.alignmentId].info.isInReadGraph = false;
            crossStrandEdgeCount += 2;
            continue;
        }

        // Otherwise, just update the disjoint sets data structure
        // with these two edges.
        disjointSets.union_set(a0, b0);
        disjointSets.union_set(a1, b1);
    }

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << readGraph.edges.size() << " total." << endl;



    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }


    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    // cout << "The read graph has " << componentMap.size() << " connected components." << endl;



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // Stop writing when we reach connected components
        // consisting of a single isolated read.
        if(component.size() == 1) {
            break;
        }

        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        // The above process of strand separation should have removed
        // all self-complementary components.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        SHASTA_ASSERT(not isSelfComplementary);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // For Mode 2 assembly, we will only assemble one connected component
    // of each pair. In each pair, we choose the component in the pair
    // that has the lowest numbered read on strand 0.
    // Then, for each read we store in its ReadFlags the strand
    // that the read appears in in this component.
    // That flag will be used in Mode 2 assembly to
    // select portions of the marker graph that should be assembled.
    uint64_t n = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If the lowest numbered read is on strand 1, this is not one of
        // the connected components we want to use.
        if(component.front().getStrand() == 1) {
            continue;
        }

        // Store the strand for each read in this component.
        for(const OrientedReadId orientedReadId: component) {
            reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
        }
        n += component.size();
    }
    SHASTA_ASSERT(n == readCount);

}









// Strict strand separation in the read graph based on using Poisson CDF.
// This guarantees that the read graph contains no self-complementary
// connected components.
// In other words, for any ReadId x, the two oriented reads
// x-0 and x-1 are guaranteed to be in distinct components of the
// read graph.
void Assembler::flagCrossStrandReadGraphEdges5()
{
    // Each alignment used in the read graph generates a pair of
    // consecutively numbered edges in the read graph
    // which are the reverse complement of each other.
    SHASTA_ASSERT((readGraph.edges.size() % 2) == 0);

    // Below, each pair is identified by the (even) index of
    // the first edge in the pair.

    // For each number of aligned markers alignedMarkerCount,
    // Gather in edgeTable[alignedMarkerCount] pairs
    // with that number of aligned markers.
    vector< pair<uint64_t, double> > edgeTable;
    const double epsilon = 1e-4;
    const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // To loop over pairs of edges, increment by 2.
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId+=2) {
        const ReadGraphEdge& edge = readGraph.edges[edgeId];
        SHASTA_ASSERT(not edge.crossesStrands);

        // Skip edges flagged as inconsistent.
        if(edge.hasInconsistentAlignment) {
            continue;
        }

        const uint64_t alignmentId = edge.alignmentId;
        const AlignmentData& alignment = alignmentData[alignmentId];

        // Skip edges involving reads classified as chimeric.
        if(getReads().getFlags(alignment.readIds[0]).isChimeric) {
            continue;
        }
        if(getReads().getFlags(alignment.readIds[1]).isChimeric) {
            continue;
        }

        // Sanity check.
        SHASTA_ASSERT(alignmentData[alignmentId].info.isInReadGraph);

        // Check that the next edge is the reverse complement of
        // this edge.
        {
            const ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
            SHASTA_ASSERT(not nextEdge.crossesStrands);
            array<OrientedReadId, 2> nextEdgeOrientedReadIds = nextEdge.orientedReadIds;
            nextEdgeOrientedReadIds[0].flipStrand();
            nextEdgeOrientedReadIds[1].flipStrand();
            SHASTA_ASSERT(nextEdgeOrientedReadIds == edge.orientedReadIds);
            SHASTA_ASSERT(nextEdge.alignmentId == alignmentId);
        }
        

        // Q(n) = (1 + δ/2ε)^n * e-δL
        // ε = 1e-4, δ = 5e-4
        // Store this pair of edges in our edgeTable.
        const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[0], 0, markers);
        const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[1], 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = alignment.info.mismatchCountRle;     

        // logQ(n) = αn - δL
        const double logQ_n = alpha * double(n) - delta * L;

        edgeTable.push_back(make_pair(edgeId, logQ_n));
        
    }

    // Sort by increasing Q
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<uint64_t, double>());

    // Print out top 100 alignments by Q value
    cout << "Top 100 alignments by logQ value:" << endl;
    cout << "EdgeId\tlogQ\tMarkerCount\tErrorRateRle\tL\tn" << endl;
    for(size_t i = 0; i < min(size_t(100), edgeTable.size()); i++) {
        const auto& p = edgeTable[i];
        const uint64_t edgeId = p.first;
        const double logQ = p.second;
        const AlignmentData& alignment = alignmentData[readGraph.edges[edgeId].alignmentId];
        const uint64_t markerCount = alignment.info.markerCount;
        const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[0], 0, markers);
        const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[1], 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = alignment.info.mismatchCountRle;
        const double errorRateRle = double(n)/(2.0*L);;
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        cout << edgeId << " ( " << readId0 << "," << readId1 << " )" << "\t" << logQ << "\t" << markerCount << "\t" << errorRateRle << "\t" << L << "\t" << n << endl;
    }


    // Create and initialize the disjoint sets data structure needed below.
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    SHASTA_ASSERT(readGraph.connectivity.size() == orientedReadCount);
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }


    // Now process our edge table in order of decreasing alignedMarkerCount.
    // Each pair consists of two edges:
    // A0--B0
    // A1--B1
    // where A1 is the reverse complement of A0 and B1 is the reverse complement of B0.
    // Maintain a disjoint set data structure for read graph vertices.
    // For each pair, compute the current connected components
    // for A0 B0 A1 B1 from the disjoint set data structure, call them a0 b0 a1 b1.
    // If a1==b0 (in which case also b1==a0), adding these edges
    // would create a self-complementary connected component,
    // so we mark them as cross-strand edges and don't add them to the
    // disjoint set data structure.
    // Otherwise, the two edges are added to the disjoint set data structure.
    uint64_t crossStrandEdgeCount = 0;
    for(auto it=edgeTable.begin(); it!=edgeTable.end(); ++it) {
        const pair<uint64_t, double>& p = *it;
        const uint64_t edgeId = p.first;
        ReadGraphEdge& edge = readGraph.edges[edgeId];
        ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
        const OrientedReadId A0 = edge.orientedReadIds[0];
        const OrientedReadId B0 = edge.orientedReadIds[1];
        const OrientedReadId A1 = nextEdge.orientedReadIds[0];
        const OrientedReadId B1 = nextEdge.orientedReadIds[1];

        SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
        SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
        SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
        SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

        // Get the connected components that these oriented reads are in.
        const uint32_t a0 = disjointSets.find_set(A0.getValue());
        const uint32_t b0 = disjointSets.find_set(B0.getValue());
        const uint32_t a1 = disjointSets.find_set(A1.getValue());
        const uint32_t b1 = disjointSets.find_set(B1.getValue());

        // If A0 and B0 are in the same connected component,
        // A1 and B1 also must be in the same connected component.
        // There is nothing to do as this pair of edges
        // does not affect the connected components.
        if (a0 == b0) {
            SHASTA_ASSERT(a1 == b1);
            continue;
        }

        // If A0 and B1 are in the same connected component,
        // A1 and B0 also must be in the same connected component.
        // Adding this pair of edges would create a self-complementary
        // connected component containing A0, B0, A1, and B1,
        // and to ensure strand separation we don't want to do that.
        // So we mark these edges as cross-strand edges
        // and don't use them to update the disjoint set data structure.
        if(a0 == b1) {
            SHASTA_ASSERT(a1 == b0);
            edge.crossesStrands = 1;
            nextEdge.crossesStrands = 1;
            alignmentData[edge.alignmentId].info.isInReadGraph = false;
            crossStrandEdgeCount += 2;
            continue;
        }

        // Otherwise, just update the disjoint sets data structure
        // with these two edges.
        disjointSets.union_set(a0, b0);
        disjointSets.union_set(a1, b1);
    }

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << readGraph.edges.size() << " total." << endl;



    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }


    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    // cout << "The read graph has " << componentMap.size() << " connected components." << endl;



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // Stop writing when we reach connected components
        // consisting of a single isolated read.
        if(component.size() == 1) {
            break;
        }

        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        // The above process of strand separation should have removed
        // all self-complementary components.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        SHASTA_ASSERT(not isSelfComplementary);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // For Mode 2 assembly, we will only assemble one connected component
    // of each pair. In each pair, we choose the component in the pair
    // that has the lowest numbered read on strand 0.
    // Then, for each read we store in its ReadFlags the strand
    // that the read appears in in this component.
    // That flag will be used in Mode 2 assembly to
    // select portions of the marker graph that should be assembled.
    uint64_t n = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If the lowest numbered read is on strand 1, this is not one of
        // the connected components we want to use.
        if(component.front().getStrand() == 1) {
            continue;
        }

        // Store the strand for each read in this component.
        for(const OrientedReadId orientedReadId: component) {
            reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
        }
        n += component.size();
    }
    SHASTA_ASSERT(n == readCount);

}











// Add this helper method to properly remove the read graph
void Assembler::removeReadGraph()
{
    // Close existing read graph if it's open
    if (readGraph.edges.isOpen) {
        readGraph.edges.close();
    }
    if (readGraph.connectivity.isOpen()) {
        readGraph.connectivity.close();
    }

    cout << timestamp << "Removed previous read graph files" << endl;
}






















































































// void Assembler::createReadGraph4(
//     uint32_t maxAlignmentCount)
// {
//     const bool debug = false;

//     // QRle threshold to use an alignment in the read graph.
//     const double minQRle = 100000.;

//     const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

//     cout << timestamp << "createReadGraph4 begins, maxAlignmentCount skata " << maxAlignmentCount << endl;

//     // Get the total number of stored alignments.
//     const uint64_t alignmentCount = alignmentData.size();
//     SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

//     // Flag all alignments as not to be kept.
//     vector<bool> keepAlignment(alignmentCount, false);

//     // The computation of projected alignments is expensive, and so it should be done once for each alignment 
//     // in an initial step and store the Error Rates of the projected alignment
//     vector<double> alignmentErrorRateRle(alignmentCount);

//     // This will hold the decomepressed Alignment.
//     // Defined here to reduce memory allocation activity.
//     Alignment alignment;

//     // Loop over all alignments.
//     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         if((alignmentId % 1000) == 0) {
//             cout << timestamp << alignmentId << "/" << alignmentCount << endl;
//         }

//         // Get information for this alignment.
//         AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         thisAlignmentData.info.isInReadGraph = 0;

//         // The alignment is stored as an alignment between readId0 on strand 0
//         // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         // The reverse complement alignment also exists, but is not stored explicitly.
//         const ReadId readId0 = thisAlignmentData.readIds[0];
//         const ReadId readId1 = thisAlignmentData.readIds[1];
//         const bool isSameStrand = thisAlignmentData.isSameStrand;
//         SHASTA_ASSERT(readId0 < readId1);
//         const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
//         const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

//         // The alignment is stored in compressed form as a string,
//         // so we have to decompress it.
//         span<const char> compressedAlignment = compressedAlignments[alignmentId];
//         shasta::decompress(compressedAlignment, alignment);

//         // Project this alignment to base space.
//         const ProjectedAlignment projectedAlignment(
//             *this,
//             {orientedReadId0, orientedReadId1},
//             alignment,
//             true);

//         const double errorRateRle = projectedAlignment.errorRateRle();

//         alignmentErrorRateRle[alignmentId] = errorRateRle;

//     }


//     // Vector to keep the alignments for each read,
//     // with their number of markers.
//     // Contains pairs(errorRateRle, alignment id).
//     vector<AlignmentStats> readAlignments;

    
//     // Find the number of reads and oriented reads.
//     const ReadId orientedReadCount = uint32_t(markers.size());
//     SHASTA_ASSERT((orientedReadCount % 2) == 0);
//     const ReadId readCount = orientedReadCount / 2;

//     // Loop over reads.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         if(debug) {
//             cout << "Working on read " << readId << endl;
//         }

//         OrientedReadId alignmentOrientedReadId0(readId, 0);

//         // Loop over the alignments that this oriented read is involved in, with the proper orientation.
//         const vector< pair<OrientedReadId, AlignmentInfo> > correctOrientedAlignments =
//             findOrientedAlignments(alignmentOrientedReadId0, false);



//         // Gather the alignments for this read, considering alignment range and right unaligned portion.
//         readAlignments.clear();

//         for(uint32_t i=0; i<correctOrientedAlignments.size(); i++) {

//             const uint32_t alignmentId = alignmentTable[alignmentOrientedReadId0.getValue()][i];

//             const double errorRateRle = alignmentErrorRateRle[alignmentId];

//             // Calculate alignment range
//             const uint32_t alignmentRange = correctOrientedAlignments[i].second.markerCount;

//             // Calculate right unaligned portion
//             const uint32_t rightUnaligned0 = correctOrientedAlignments[i].second.rightTrim(0);

//             // Calculate left unaligned portion
//             const uint32_t leftUnaligned0 = correctOrientedAlignments[i].second.leftTrim(0);

//             // Calculate right unaligned portion
//             const uint32_t rightUnaligned1 = correctOrientedAlignments[i].second.rightTrim(1);

//             // Calculate left unaligned portion
//             const uint32_t leftUnaligned1 = correctOrientedAlignments[i].second.leftTrim(1);

//             // if(rightUnaligned1 != 0 and leftUnaligned1 != 0 and errorRateRle<= maxErrorRateRle) {
//             //     readAlignments.push_back(AlignmentStats{errorRateRle, alignmentRange, rightUnaligned1, leftUnaligned1, alignmentId});
//             // }

//             if(errorRateRle<= maxErrorRateRle) {
//                 readAlignments.push_back(AlignmentStats{errorRateRle, alignmentRange, rightUnaligned1, leftUnaligned1, alignmentId});
//             }

//         }
        
//         cout << "Working on read " << readId << endl;


//         // // Find the alignment with the highest alignedRange
//         // auto bestAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return a.alignedRange < b.alignedRange;
//         //     });

//         // if (bestAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestAlignment = *bestAlignmentIt;

//         //     // Clear existing alignments and keep only the best one
//         //     readAlignments.clear();
//         //     readAlignments.push_back(bestAlignment);

//         //     if(debug) {
//         //         cout << "Kept alignment with the highest alignedRange: " << bestAlignment.alignedRange << endl;
//         //     }
//         // } else {
//         //     if(debug) {
//         //         cout << "No alignments found for this read." << endl;
//         //     }
//         // }


//         // // Find the 10 alignments with the lowest errorRateRle
//         // std::partial_sort(readAlignments.begin(), 
//         //                   readAlignments.begin() + std::min(10UL, readAlignments.size()), 
//         //                   readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return a.errorRateRle < b.errorRateRle;
//         //     });

//         // // Find the 10 alignments with the highest sum of rightUnaligned, leftUnaligned, and alignedRange
//         // std::partial_sort(readAlignments.begin(), 
//         //                   readAlignments.begin() + std::min(10UL, readAlignments.size()), 
//         //                   readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return (a.rightUnaligned + a.leftUnaligned + a.alignedRange) >
//         //                (b.rightUnaligned + b.leftUnaligned + b.alignedRange);
//         //     });

//         cout << "read " << readId << " has that many alignments skata:" << readAlignments.size() << endl;

//         // // Keep only the top 10 alignments (or fewer if there are less than 10)
//         // readAlignments.resize(std::min(10UL, readAlignments.size()));

//         // cout << "read " << readId << " has that many alignments:" << readAlignments.size() << endl;

//         // if(debug) {
//         //     cout << "Top 5 alignments (or fewer) based on sum of unaligned portions and aligned range:" << endl;
//         //     for(const auto& alignment : readAlignments) {
//         //         cout << "AlignmentId: " << alignment.alignmentId 
//         //              << ", Sum: " << (alignment.rightUnaligned + alignment.leftUnaligned + alignment.alignedRange)
//         //              << " (Right: " << alignment.rightUnaligned 
//         //              << ", Left: " << alignment.leftUnaligned 
//         //              << ", Aligned: " << alignment.alignedRange << ")" << endl;
//         //     }
//         // }


//         // // Find the alignment with the highest alignedRange which also has the highest rightUnaligned
//         // auto bestRightAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         // if (a.alignedRange == b.alignedRange) {
//         //         //     return a.rightUnaligned > b.rightUnaligned;
//         //         // }
//         //         // return a.alignedRange < b.alignedRange;
//         //         return a.rightUnaligned+a.leftUnaligned+a.alignedRange < b.rightUnaligned+b.leftUnaligned+b.alignedRange;
//         //     });

//         // if (bestRightAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestRightAlignment = *bestRightAlignmentIt;

//         //     readAlignmentsRightTrimmed.clear();
//         //     readAlignmentsRightTrimmed.push_back(bestRightAlignment);
//         // }


//         // // Find the alignment with the highest alignedRange which also has the highest leftUnaligned
//         // auto bestLeftAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         // if (a.alignedRange == b.alignedRange) {
//         //         //     return a.leftUnaligned > b.leftUnaligned;
//         //         // }
//         //         // return a.alignedRange < b.alignedRange;
//         //         return a.leftUnaligned < b.leftUnaligned;
//         //     });

//         // if (bestLeftAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestLeftAlignment = *bestLeftAlignmentIt;

//         //     readAlignmentsLeftTrimmed.clear();
//         //     readAlignmentsLeftTrimmed.push_back(bestLeftAlignment);
//         // }

        

//         // if (!readAlignmentsRightTrimmed.empty()){
//         // cout << "The best rightTrimm alignment is " << readAlignmentsRightTrimmed[0].alignmentId << " alignmentId." << endl;
//         // cout << "The best rightTrimm has " << readAlignmentsRightTrimmed[0].rightUnaligned << " rightUnaligned." << endl;
//         // }

//         // if (!readAlignmentsLeftTrimmed.empty()){
//         //     cout << "The best leftTrimm alignment is " << readAlignmentsLeftTrimmed[0].alignmentId << " alignmentId." << endl;
//         //     cout << "The best leftTrimm has " << readAlignmentsRightTrimmed[0].leftUnaligned << " leftUnaligned." << endl;
//         // }
        

//         // if(debug) {
//         //     cout << "Found " << readAlignments.size() << " alignments." << endl;
//         // }

//         // // Keep the best maxAlignmentCount.
//         // if(readAlignments.size() > maxAlignmentCount) {
//         //     std::nth_element(
//         //         readAlignments.begin(),
//         //         readAlignments.begin() + maxAlignmentCount,
//         //         readAlignments.end(),
//         //         std::less< pair<double, uint32_t> >());
//         //     readAlignments.resize(maxAlignmentCount);
//         // }
//         // if(debug) {
//         //     cout << "Kept " << readAlignments.size() << " alignments." << endl;
//         // }

//         // Mark the surviving alignments as to be kept.
//         for(const auto& p: readAlignments) {
//             // Get information for this alignment.
//             const uint32_t alignmentId = p.alignmentId;
//             AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//             keepAlignment[alignmentId] = true;
//             thisAlignmentData.info.isInReadGraph = 1;
//             if(debug) {
//                 const AlignmentData& alignment = alignmentData[alignmentId];
//                 cout << "Marked alignment " << alignment.readIds[0] << " " <<
//                     alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
//             }
//         }

//         // // Mark the surviving alignments as to be kept.
//         // for(const auto& p: readAlignmentsLeftTrimmed) {
//         //     // Get information for this alignment.
//         //     const uint32_t alignmentId = p.alignmentId;
//         //     AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         //     keepAlignment[alignmentId] = true;
//         //     thisAlignmentData.info.isInReadGraph = 1;
//         //     if(debug) {
//         //         const AlignmentData& alignment = alignmentData[alignmentId];
//         //         cout << "Marked alignment " << alignment.readIds[0] << " " <<
//         //             alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
//         //     }
//         // }
//     }


//     cout << timestamp << "Done processing alignments." << endl;

//     const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

//     // Create the read graph using the alignments we selected.
//     createReadGraphUsingSelectedAlignments(keepAlignment);

//     cout << timestamp << "createReadGraph4 ends." << endl;
// }

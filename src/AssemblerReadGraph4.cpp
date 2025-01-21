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
#include <unordered_map>

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








// void ReadGraph4AllAlignments::computeShortPath(
//     OrientedReadId orientedReadId0,
//     OrientedReadId orientedReadId1,
//     size_t maxDistance,
//     vector<uint32_t>& path,
//     vector<uint32_t>& distance,
//     vector<OrientedReadId>& reachedVertices,
//     vector<uint32_t>& parentEdges,
//     MemoryMapped::Vector<AlignmentData>& alignmentData,
//     ReadGraph4AllAlignments& readGraph)
// {
//     const bool debug = false;
//     path.clear();
//     reachedVertices.clear();

//     // Initialize all distances to infinite
//     fill(distance.begin(), distance.end(), ReadGraph::infiniteDistance);

//     // Initialize BFS
//     std::queue<OrientedReadId> queuedVertices;
//     queuedVertices.push(orientedReadId0);
//     distance[orientedReadId0.getValue()] = 0;
//     reachedVertices.push_back(orientedReadId0);

//     // Do the BFS
//     while(!queuedVertices.empty()) {
//         const OrientedReadId vertex0 = queuedVertices.front();
//         queuedVertices.pop();
//         const uint32_t distance0 = distance[vertex0.getValue()];
//         const uint32_t distance1 = distance0 + 1;

//         if(distance1 > maxDistance) {
//             continue;
//         }


//         BGL_FORALL_OUTEDGES(vertex0.getValue(), edge, *this, ReadGraph4AllAlignments) {
        
//             vertex_descriptor targetVertex = target(edge, *this);
//             const OrientedReadId vertex1 = OrientedReadId::fromValue(targetVertex);

//             uint32_t alignmentId = readGraph[edge].alignmentId;
//             AlignmentData alignment = alignmentData[alignmentId];

//             //uint32_t alignmentId = ReadGraph4AllAlignments[edge].alignmentId;

//             // Process new vertices
//             if(distance[vertex1.getValue()] == ReadGraph::infiniteDistance) {
//                 distance[vertex1.getValue()] = distance1;
//                 reachedVertices.push_back(vertex1);
//                 parentEdges[vertex1.getValue()] = alignmentId;
//                 queuedVertices.push(vertex1);

//                 // Check if we reached the target
//                 if(vertex1 == orientedReadId1) {
//                     // Reconstruct path
//                     OrientedReadId vertex = vertex1;
//                     while(vertex != orientedReadId0) {
//                         const uint32_t alignmentId = parentEdges[vertex.getValue()];
//                         path.push_back(alignmentId);
//                         AlignmentData alignment = alignmentData[alignmentId];
//                         const ReadId readId0 = alignment.readIds[0];
//                         const ReadId readId1 = alignment.readIds[1];
//                         const bool isSameStrand = alignment.isSameStrand;
//                         SHASTA_ASSERT(readId0 < readId1);
//                         const OrientedReadId A0 = OrientedReadId(readId0, 0);
//                         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);

//                         if(vertex.getValue() == A0.getValue()) {
//                             vertex = B0;
//                         } else {
//                             vertex = A0;
//                         }
//                     }
//                     std::reverse(path.begin(), path.end());
//                     return;
//                 }
//             }
                

//         }

//         // for(const uint32_t edgeId: connectivity[vertex0.getValue()]) {
//         //     const ReadGraphEdge& edge = edges[edgeId];
//         //     if(edge.crossesStrands) {
//         //         continue;
//         //     }
//         //     const OrientedReadId vertex1 = edge.getOther(vertex0);

//         //     // Process new vertices
//         //     if(distance[vertex1.getValue()] == ReadGraph::infiniteDistance) {
//         //         distance[vertex1.getValue()] = distance1;
//         //         reachedVertices.push_back(vertex1);
//         //         parentEdges[vertex1.getValue()] = edgeId;
//         //         queuedVertices.push(vertex1);

//         //         // Check if we reached the target
//         //         if(vertex1 == orientedReadId1) {
//         //             // Reconstruct path
//         //             OrientedReadId vertex = vertex1;
//         //             while(vertex != orientedReadId0) {
//         //                 const uint32_t edgeId = parentEdges[vertex.getValue()];
//         //                 path.push_back(edgeId);
//         //                 vertex = edges[edgeId].getOther(vertex);
//         //             }
//         //             std::reverse(path.begin(), path.end());
//         //             return;
//         //         }
//         //     }
//         // }
//     }
// }







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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
            // OrientedReadId::fromValue(currentVertex).getValue()
            if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[currentVertex]) {
                neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    // OrientedReadId::fromValue(targetVertex).getValue()
                    if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[targetVertex]) {
                        neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)));
                        return;
                    }
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }

            // Process incoming edges
            BGL_FORALL_INEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor sourceVertex = source(edge, *this);
                if (!visitedVertices.contains(sourceVertex)) {
                    if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[sourceVertex]) {
                        neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(sourceVertex)));
                        return;
                    }
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
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                uint64_t alignmentId = readGraphAllAlignments[edge].alignmentId;
                if (!visitedVertices.contains(targetVertex)) {
                    neighborsAlignmentIds.push_back(alignmentId);
                    // OrientedReadId::fromValue(targetVertex).getValue()
                    if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[targetVertex]) {
                        neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
                        neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(currentVertex)));
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
                        neighbors.push_back(OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)));
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

            uint64_t alignmentId = readGraph[edge].alignmentId;

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
            
            uint64_t result = findPathWithPositiveOffset(OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)), paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

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
         
            uint64_t alignmentId = readGraph[edge].alignmentId;

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
            
            uint64_t result = findPathWithPositiveOffset(OrientedReadId::fromValue(static_cast<ReadId>(sourceVertex)), paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

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
                OrientedReadId next = OrientedReadId::fromValue(static_cast<ReadId>(nextVertex));
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
            current = OrientedReadId::fromValue(static_cast<ReadId>(parent[current.getValue()]));
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
                OrientedReadId next = OrientedReadId::fromValue(static_cast<ReadId>(nextVertex));
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
            current = OrientedReadId::fromValue(static_cast<ReadId>(parent[current.getValue()]));
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
            uint64_t result = findAllPathsToNode(OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)), endNode, paths, currentPath, visited, maxDistance, currentDistance + 1);
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
                OrientedReadId::fromValue(static_cast<ReadId>(targetVertex)), 
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
    uint32_t maxAlignmentCount,
    double epsilon,
    double delta,
    double WThreshold,
    double WThresholdForBreaks
    )
{
    cout << timestamp << "createReadGraph4 with strand separation begins" << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);


    //*
    //
    // Order alignments in order of increasing Q. 
    //
    // Gather in alignmentTable[alignmentID, Q]
    // alignments in order of increasing Q.
    // Q(n) = (1 + δ/2ε)^n * e-δL
    // ε = 1e-4, δ = 5e-4
    // logQ(n) = αn - δL
    //
    //*

    // const double epsilon = 1e-4;
    // const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // const double WThreshold = 1e-7;
    const double logWThreshold = log(WThreshold);

    // const double WThresholdForBreaks = 1e+15;
    const double logWThresholdForBreaks = log(WThresholdForBreaks);



    vector< pair<uint64_t, double> > alignmentTable;
    vector< pair<uint64_t, double> > alignmentTableNotPassFilter;
    vector< pair<uint64_t, double> > QAlignments;
    

    // Get stats about the reads
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    
    // Keep track of which readIds were used in alignments
    vector<bool> readUsed(readCount, false);

    // Flag alignments to be kept for break detection.
    vector<bool> keepAlignmentsForBreaks(alignmentCount, false);


    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if((alignmentId % 100000) == 0) {
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
        // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
        const double errorRateRle = thisAlignmentData.info.errorRateRle;
        const double nRLE = errorRateRle * 2 * L;
        // const double markerCount = thisAlignmentData.info.markerCount;

        // logQ(n) = αn - δL
        const double logQ = alpha * double(nRLE) - delta * L;

        if (logQ <= logWThreshold) {
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
    cout << "The alignmentTable has " << alignmentTable.size() << " entries." << endl;
    cout << "The alignmentTableNotPassFilter has " << alignmentTableNotPassFilter.size() << " entries." << endl;
    





    ///
    // Process alignments in order of increasing Q. 
    //
    // i.   Start with no edges in the read graph. 
    // ii.  Process alignments in order of increasing Q. 
    // iii. If the alignment breaks strand separation, it is skipped. 
    // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped.  
    // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.
    
    // Maintain a vector containing the degree of each vertex
    // verticesDegree[vertexID] -> degree
    vector<uint64_t> verticesDegree(orientedReadCount, 0);

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
            // const double logQ = p.second;

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
    for (std::size_t i = 0; i < orientedReadCount; ++i) {
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

    cout << "The read graph has " << num_vertices(readGraph) << " vertices and " << num_edges(readGraph) << " edges." << endl;



    //*
    //
    // Create the dynamically adjustable boost readGraph using the alignments that did not pass the strict filter.
    // These alignments are used to create the read graph that will aid in the break detection.
    //
    //*
    
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
    
    cout << "The read graph for break detection has " << num_vertices(readGraphAllAlignments) << " vertices and " << num_edges(readGraphAllAlignments) << " edges." << endl;





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
    uint64_t finalNumberOfPotentialDeadEndNodes = 0;
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
            uint64_t maxDistance = 4;
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
                finalNumberOfPotentialDeadEndNodes++;
                // cout << "Did not find a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset. Keeping it as a potential dead end read." << endl;
                finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()] = true;
                finalDeadEndReadsWithNoIncomingNodes[reverseOrientedReadId.getValue()] = true;
            }
        }
    }

    cout << "After filtering we are left with " << finalNumberOfPotentialDeadEndNodes << " potential dead end reads." << endl;


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











    // //*
    // //
    // // Extend the potential dead end nodes list.
    // // Add neighboring nodes of potential dead end nodes to the dead end node list.
    // //
    // //*
    // vector<bool> finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors(orientedReadCount, false);
    // vector<bool> finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors(orientedReadCount, false);
    // for (ReadId readId = 0; readId < readCount; readId++) {
    //     for (Strand strand = 0; strand < 2; strand++) {
    //         OrientedReadId orientedReadId(readId, strand);

    //         if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){

    //             finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()] = true;
                
    //             OrientedReadId reverseOrientedReadId = orientedReadId;
    //             reverseOrientedReadId.flipStrand();

    //             finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;

    //             // Find distance 1 neighbors
    //             vector<OrientedReadId> distance1Neighbors;
    //             readGraph.findNeighborsDirectedGraphBothSides(orientedReadId, 1, distance1Neighbors);

    //             for(auto neighbor : distance1Neighbors) {
    //                 finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[neighbor.getValue()] = true;

    //                 OrientedReadId reverseOrientedReadId = neighbor;
    //                 reverseOrientedReadId.flipStrand();

    //                 finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;
    //             }
    //         }
    //     }
    // }

    // // // print dead end Oriented reads
    // // // iterate over all oriented reads
    // // for (ReadId readId = 0; readId < readCount; readId++) {
    // //     for (Strand strand = 0; strand < 2; strand++) {
    // //         OrientedReadId orientedReadId(readId, strand);
    // //         if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
    // //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
    // //         } else if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
    // //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
    // //         }
    // //     }
    // // }





    //*
    //
    // Create a map of endNodesWithNoOutgoingNodes to other endNodesWithNoIncomingNodes that are possible to connect to.
    //
    //*

    // Case 1: the NoOut deadEnd node will be mapped to at least 3 NoIn deadEnd nodes 
    // IN THE SAME disjointSet.
    //
    //  NoIn - - - - | - - - - -  - - - - - - | - - - - -
    //  NoIn - - - - | - - NoOut     NoIn - - | - - - - -

    std::unordered_map<uint32_t, vector<uint32_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet;
    vector<bool> endNodesWithNoOutgoingNodesConsidered(orientedReadCount, false);

    // First, we try to find EndNodesWithNoIncomingNodes that are in the same connected component
    // as EndNodesWithNoOutgoingNodes. These have priority over other EndNodesWithNoIncomingNodes in other connected components.
    // These will also include telomeric nodes because they do not have incoming nodes!
    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {
            OrientedReadId orientedReadId(readId, strand);

            // check if the orientedReadId is an endNode with no outgoing nodes
            if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
                // Get it's disjointSet
                uint32_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

                // Find all EndNodesWithNoIncomingNodes that are in the same disjointSet
                for(uint32_t id=0; id<orientedReadCount; id++) {
                    // check if the id is an endNode with no incoming nodes
                    if(finalDeadEndReadsWithNoIncomingNodes[id]){
                        // get the disjointSet of the id
                        uint32_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
                        // check if the disjointSet of the id is the same as the disjointSet of the endNode with no outgoing nodes
                        if(deadEndReadWithNoOutgoingNodesDisjointSetId == endNodeWithNoIncomingNodesDisjointSetId){
                            endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet[orientedReadId.getValue()].push_back(id);
                            endNodesWithNoOutgoingNodesConsidered[orientedReadId.getValue()] = true;
                        }
                    }
                }

            }
        }
    }


    for(auto& p : endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet) {
        uint32_t value = p.first;
        OrientedReadId orientedReadId = OrientedReadId::fromValue(value);
        SHASTA_ASSERT(orientedReadId.getValue() == value);
        // const ReadId readId = orientedReadId.getReadId();
        // const Strand strand = orientedReadId.getStrand();

        // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in the same disjointSet:" << endl;
        vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
        for(auto& node : p.second) {
            OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(node);
            // cout << "ReadID " << nodeOrientedReadId.getReadId() << " strand " << nodeOrientedReadId.getStrand() << endl;
            SHASTA_ASSERT(nodeOrientedReadId.getValue() == node);
            endNodesWithNoIncomingNodes[node] = true;
        }

        
        // Find neighbors in the forward direction of the ALL ALIGNMENTS read graph 
        // starting from orientedReadId which is an endNode with no outgoing nodes in the filtered read graph. 
        // Early stop when we reach an endNode with no incoming nodes (a node with endNodesWithNoIncomingNodes set to true).
        vector<OrientedReadId> forwardNeighbors;
        readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, endNodesWithNoIncomingNodes, 5, forwardNeighbors);

        // // print the forward neighbors
        // cout << "Forward neighbors of the endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " are:" << endl;
        // for(auto& neighbor : forwardNeighbors) {
        //     cout << "ReadID " << neighbor.getReadId() << " strand " << neighbor.getStrand() << endl;
        // }


        // create a std::set of the forwardNeighbors for easy contain check
        std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
        

        if(forwardNeighbors.empty()) {
            // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
            // cout << "No forward neighbors found" << endl;
            continue;
        }


        // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
        // OR a non speficic node if we exceeded maxDistance
        OrientedReadId lastNode = forwardNeighbors.back();


        bool success = false;


        // First check if the last node is a potential dead end node with no INCOMING nodes
        if(endNodesWithNoIncomingNodes[lastNode.getValue()]) {

            // cout << "The last node is: " << lastNode.getReadId() << " strand " << lastNode.getStrand() << endl;
            
            const OrientedReadId A0 = orientedReadId;
            const OrientedReadId B0 = lastNode;
            const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
            const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

            SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
            SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
            SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
            SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

            // Get the connected components that these oriented reads are in.
            // const uint32_t a0 = disjointSets.find_set(A0.getValue());
            // const uint32_t b0 = disjointSets.find_set(B0.getValue());
            // const uint32_t a1 = disjointSets.find_set(A1.getValue());
            // const uint32_t b1 = disjointSets.find_set(B1.getValue());


            for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
            // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

                const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
                const uint64_t alignmentId = p.first;
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


                if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                    
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

                    // Make sure all alignments added are not considered as dead ends anymore
                    finalDeadEndReadsWithNoOutgoingNodes[A0v2.getValue()] = false;
                    finalDeadEndReadsWithNoOutgoingNodes[A1v2.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B0v2.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B1v2.getValue()] = false;

                    // Make sure the start and last nodes are not considered as dead ends anymore
                    finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] = false;
                    finalDeadEndReadsWithNoOutgoingNodes[A1.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B1.getValue()] = false;

                    success = true;

                    // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

                    // Create the edge.
                    add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                    // Also create the reverse complemented edge.
                    alignmentOrientedReadId0.flipStrand();
                    alignmentOrientedReadId1.flipStrand();
                    add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


                }

            }

            // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
            // it means that we have connected them with alignments.
            if(success) {
                // cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
            }

            // // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
            // // it means that we have connected them with alignments.
            // if(finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] == false and finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] == false) {
            //     cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
            // } else {
            //     cout << "Did not connect endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with any other endNode with no incoming nodes" << endl;
            // }

        }

        if(!success) {
           // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
           // cout << "No alignments added" << endl;
        }

            
    }


    // Print how many alignments were kept
    const size_t keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Adding alignments for break bridging to connect endNodes in the same disjointSet: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;




    
    // Case 2: the NoOut deadEnd node will be mapped to at least 2 NoIn deadEnd nodes in the same disjointSet.
    // Other NoIn deadEnd nodes will be in a different disjointSet.
    //
    //  NoIn - - - - | - - - - -  - - - - - - - - - - -
    //  NoIn - - - - | - - NoOut     
    //                                     NoIn - - - - - - - (Other disjointSet)

    // Case 3: Breaks in a haploid chromosome (chrX chrY).
    //
    //  NoIn - - - - - - NoOut
    //                             NoIn - - - - - - - (Other disjointSet) 


    std::unordered_map<uint32_t, vector<uint32_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet;

    // Now, we try to find EndNodesWithNoIncomingNodes that are in different disjointSets than EndNodesWithNoOutgoingNodes.
    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {
            OrientedReadId orientedReadId(readId, strand);

            // check if the orientedReadId is an endNode with no outgoing nodes
            if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
                // Get it's disjointSet
                uint32_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

                // Find all EndNodesWithNoIncomingNodes that are in different disjointSet
                for(uint32_t id=0; id<orientedReadCount; id++) {
                    // check if the id is an endNode with no incoming nodes
                    if(finalDeadEndReadsWithNoIncomingNodes[id]){
                        // get the disjointSet of the id
                        uint32_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
                        // check if the disjointSet of the id is different from the disjointSet of the endNode with no outgoing nodes
                        if(deadEndReadWithNoOutgoingNodesDisjointSetId != endNodeWithNoIncomingNodesDisjointSetId){
                            endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet[orientedReadId.getValue()].push_back(id);
                            endNodesWithNoOutgoingNodesConsidered[orientedReadId.getValue()] = true;
                        }
                    }
                }

            }
        }
    }


    for(auto& p : endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet) {
        uint32_t value = p.first;
        OrientedReadId orientedReadId = OrientedReadId::fromValue(value);
        SHASTA_ASSERT(orientedReadId.getValue() == value);
        // const ReadId readId = orientedReadId.getReadId();
        // const Strand strand = orientedReadId.getStrand();

        // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in a different disjointSet: " << endl;
        vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
        for(auto& node : p.second) {
            OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(node);
            // cout << "ReadID " << nodeOrientedReadId.getReadId() << " strand " << nodeOrientedReadId.getStrand() << endl;
            SHASTA_ASSERT(nodeOrientedReadId.getValue() == node);
            endNodesWithNoIncomingNodes[node] = true;
        }
        
        // Find neighbors in the forward direction of the ALL ALIGNMENTS read graph 
        // starting from orientedReadId which is an endNode with no outgoing nodes in the filtered read graph. 
        // Early stop when we reach an endNode with no incoming nodes (a node with endNodesWithNoIncomingNodes set to true).
        vector<OrientedReadId> forwardNeighbors;
        readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, endNodesWithNoIncomingNodes, 5, forwardNeighbors);

        // create a std::set of the forwardNeighbors for easy contain check
        std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
        

        if(forwardNeighbors.empty()) {
            // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
            // cout << "No forward neighbors found" << endl;
            continue;
        }

        // // print the forward neighbors
        // cout << "Forward neighbors of the endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " are:" << endl;
        // for(auto& neighbor : forwardNeighbors) {
        //     cout << "ReadID " << neighbor.getReadId() << " strand " << neighbor.getStrand() << endl;
        // }

        // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
        // OR a non speficic node if we exceeded maxDistance
        OrientedReadId lastNode = forwardNeighbors.back();

        bool success = false;

        // First check if the last node is a potential dead end node with no INCOMING nodes
        if(endNodesWithNoIncomingNodes[lastNode.getValue()]) {
            
            const OrientedReadId A0 = orientedReadId;
            const OrientedReadId B0 = lastNode;
            const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
            const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

            SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
            SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
            SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
            SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

            // Get the connected components that these oriented reads are in.
            // const uint32_t a0 = disjointSets.find_set(A0.getValue());
            // const uint32_t b0 = disjointSets.find_set(B0.getValue());
            // const uint32_t a1 = disjointSets.find_set(A1.getValue());
            // const uint32_t b1 = disjointSets.find_set(B1.getValue());


            for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
            // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

                const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
                const uint64_t alignmentId = p.first;
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


                if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                    
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

                    // Make sure all alignments added are not considered as dead ends anymore
                    finalDeadEndReadsWithNoOutgoingNodes[A0v2.getValue()] = false;
                    finalDeadEndReadsWithNoOutgoingNodes[A1v2.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B0v2.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B1v2.getValue()] = false;

                    // Make sure the start and last nodes are not considered as dead ends anymore
                    finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] = false;
                    finalDeadEndReadsWithNoOutgoingNodes[A1.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] = false;
                    finalDeadEndReadsWithNoIncomingNodes[B1.getValue()] = false;

                    success = true;


                    // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

                    // Create the edge.
                    add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

                    // Also create the reverse complemented edge.
                    alignmentOrientedReadId0.flipStrand();
                    alignmentOrientedReadId1.flipStrand();
                    add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


                }

            }

            // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
            // it means that we have connected them with alignments.
            if(success) {
                // cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
            }


        }

        if(!success) {
            // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
            // cout << "No alignments added" << endl;
        }

            
    }


    // Print how many alignments were kept
    const size_t keepCountR3 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Adding alignments for break bridging to connect endNodes in different disjointSets: Keeping " << keepCountR3 << " alignments of " << keepAlignment.size() << endl;



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















    










    //
    //
    // Find forbidden non used reads
    //
    //
    // vector<bool> forbiddenReads(orientedReadCount, false);

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
    uint32_t /* maxAlignmentCount */)
{
    // const bool debug = false;

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
        // const double m = L * 0.0001;
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






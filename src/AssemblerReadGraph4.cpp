// Ignore some warnings in this file for now.
// But they should eventually be fixed.
// #pragma GCC diagnostic ignored "-Wunused-variable"
// #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
// #pragma GCC diagnostic ignored "-Wunused-parameter"
// #pragma GCC diagnostic ignored "-Wconversion"
// #pragma GCC diagnostic ignored "-Wfloat-conversion"
// #pragma GCC diagnostic ignored "-Wsign-compare"

#include "Assembler.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "diploidBayesianPhase.hpp"
#include "extractKmer.hpp"
#include "performanceLog.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "shastaTypes.hpp"
#include "timestamp.hpp"
#include "orderPairs.hpp"
#include "Mode3Assembler.hpp"
#include <iostream>
#include <vector>
#include "deduplicate.hpp"
#include <unordered_set>
#include <boost/icl/interval_map.hpp> // Include Boost Interval Container Library header
#include <set> // Include set for the interval map value

using namespace shasta;

// Standard library.
#include "fstream.hpp"
#include <queue>
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
    
    bool findPathWithPositiveOffset(
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
        uint64_t maxDistance,
        vector<uint64_t>& path,
        vector<uint64_t>& distance,
        vector<OrientedReadId>& reachedVertices,
        vector<uint64_t>& parentEdges,
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







/*
void ReadGraph4AllAlignments::computeShortPath(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    uint64_t maxDistance,
    vector<uint64_t>& path,
    vector<uint64_t>& distance,
    vector<OrientedReadId>& reachedVertices,
    vector<uint64_t>& parentEdges,
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
        const uint64_t distance0 = distance[vertex0.getValue()];
        const uint64_t distance1 = distance0 + 1;

        if(distance1 > maxDistance) {
            continue;
        }


        BGL_FORALL_OUTEDGES(vertex0.getValue(), edge, *this, ReadGraph4AllAlignments) {
        
            vertex_descriptor targetVertex = target(edge, *this);
            const OrientedReadId vertex1 = OrientedReadId::fromValue(targetVertex);

            uint64_t alignmentId = readGraph[edge].alignmentId;
            AlignmentData alignment = alignmentData[alignmentId];

            //uint64_t alignmentId = ReadGraph4AllAlignments[edge].alignmentId;

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
                        const uint64_t alignmentId = parentEdges[vertex.getValue()];
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

        // for(const uint64_t edgeId: connectivity[vertex0.getValue()]) {
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
        //                 const uint64_t edgeId = parentEdges[vertex.getValue()];
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
*/







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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    // OrientedReadId::fromValue(targetVertex).getValue()
                    if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[targetVertex]) {
                        neighbors.push_back(OrientedReadId::fromValue(ReadId(targetVertex)));
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
                        neighbors.push_back(OrientedReadId::fromValue(ReadId(sourceVertex)));
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
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
                        neighbors.push_back(OrientedReadId::fromValue(ReadId(targetVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
                        neighbors.push_back(OrientedReadId::fromValue(ReadId(targetVertex)));
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
                neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(ReadId(currentVertex)));
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
                        neighbors.push_back(OrientedReadId::fromValue(ReadId(targetVertex)));
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
bool ReadGraph4::findPathWithPositiveOffset(
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
            
            bool result = findPathWithPositiveOffset(OrientedReadId::fromValue(ReadId(targetVertex)), paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

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
            
            bool result = findPathWithPositiveOffset(OrientedReadId::fromValue(ReadId(sourceVertex)), paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

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
                OrientedReadId next = OrientedReadId::fromValue(ReadId(nextVertex));
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
            current = OrientedReadId::fromValue(ReadId(parent[current.getValue()]));
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
                OrientedReadId next = OrientedReadId::fromValue(ReadId(nextVertex));
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
            current = OrientedReadId::fromValue(ReadId(parent[current.getValue()]));
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
            uint64_t result = findAllPathsToNode(OrientedReadId::fromValue(ReadId(targetVertex)), endNode, paths, currentPath, visited, maxDistance, currentDistance + 1);
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
                OrientedReadId::fromValue(ReadId(targetVertex)), 
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





































// Find in the alignment table the alignments involving
// a given oriented read, and return them with the correct
// orientation (this may involve a swap and/or reverse complement
// of the AlignmentInfo stored in the alignmentTable).
vector< pair<OrientedReadId, uint64_t> >
    Assembler::findOrientedAlignmentsPlusIds(
        OrientedReadId orientedReadId0Argument) const
{
    const ReadId readId0 = orientedReadId0Argument.getReadId();
    const ReadId strand0 = orientedReadId0Argument.getStrand();

    vector< pair<OrientedReadId, uint64_t> > result;

    // Loop over alignment involving this read, as stored in the
    // alignment table.
    const auto alignmentTable0 = alignmentTable[orientedReadId0Argument.getValue()];
    for(const auto alignmentId: alignmentTable0) {
        const AlignmentData& ad = alignmentData[alignmentId];

        // Get the oriented read ids that the AlignmentData refers to.
        OrientedReadId orientedReadId0(ad.readIds[0], 0);
        OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        AlignmentInfo alignmentInfo = ad.info;

        // Swap oriented reads, if necessary.
        if(orientedReadId0.getReadId() != readId0) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }
        SHASTA_ASSERT(orientedReadId0.getReadId() == readId0);

        // Reverse complement, if necessary.
        if(orientedReadId0.getStrand() != strand0) {
            orientedReadId0.flipStrand();
            orientedReadId1.flipStrand();
            alignmentInfo.reverseComplement();
        }
        SHASTA_ASSERT(orientedReadId0.getStrand() == strand0);
        SHASTA_ASSERT(orientedReadId0 == orientedReadId0Argument);

        result.push_back(make_pair(orientedReadId1, alignmentId));
    }
    return result;
}





















bool isSiteInHomopolymerRegion(
    uint64_t sitePositionInRead,
    const shasta::LongBaseSequenceView& readSequence,
    uint64_t homopolymerThreshold = 3)
{
    const uint64_t readLength = readSequence.baseCount;
    SHASTA_ASSERT(sitePositionInRead < readLength);

    const shasta::Base siteBase = readSequence[sitePositionInRead];
    const uint8_t siteBaseValue = siteBase.value; // Assuming Base::value exists

    // Determine scan boundaries for characters to the left/right of the site.
    // These are the furthest characters (inclusive) that the scan loops will consider.
    uint64_t leftScanLoopEndPos = 0;
    if (sitePositionInRead - homopolymerThreshold < 0) {
        leftScanLoopEndPos = 0;
    } else {
        leftScanLoopEndPos = sitePositionInRead - homopolymerThreshold;
    }

    uint64_t rightScanLoopEndPos = 0;
    if (sitePositionInRead + homopolymerThreshold >= readLength) {
        rightScanLoopEndPos = readLength - 1;
    } else {
        rightScanLoopEndPos = sitePositionInRead + homopolymerThreshold;
    }


    // --- Forward scan (strictly to the right of sitePositionInRead) ---
    uint8_t forwardRunBaseVal = 0; // Will hold the base value of the forward run
    uint64_t forwardRunLength = 0;
    bool forwardBaseFound = false;

    // Start scanning from the base immediately to the right of the site
    if (sitePositionInRead + 1 < readLength) { // Check if there is any base to the right
        for (uint64_t i = sitePositionInRead + 1; i <= rightScanLoopEndPos; ++i) {    
            const shasta::Base currentForwardBase = readSequence[i];
            if (!forwardBaseFound) {
                forwardRunBaseVal = currentForwardBase.value;
                forwardRunLength = 1;
                forwardBaseFound = true;
            } else {
                if (currentForwardBase.value == forwardRunBaseVal) {
                    forwardRunLength++;
                } else {
                    break; // Homopolymer run ended
                }
            }
        }
    }

    // --- Backward scan (strictly to the left of sitePositionInRead) ---
    uint8_t backwardRunBaseVal = 0; // Will hold the base value of the backward run
    uint64_t backwardRunLength = 0;
    bool backwardBaseFound = false;

    // Start scanning from the base immediately to the left of the site
    if (sitePositionInRead > 0) { // Check if there is any base to the left
        for (int64_t i = static_cast<int64_t>(sitePositionInRead) - 1; i >= leftScanLoopEndPos; --i) {
            const shasta::Base currentBackwardBase = readSequence[static_cast<uint64_t>(i)];
            if (!backwardBaseFound) {
                backwardRunBaseVal = currentBackwardBase.value;
                backwardRunLength = 1;
                backwardBaseFound = true;
            } else {
                if (currentBackwardBase.value == backwardRunBaseVal) {
                    backwardRunLength++;
                } else {
                    break; // Homopolymer run ended
                }
            }
        }
    }

    // --- Combine results and check conditions ---
    // 'effectiveForwardLength' and 'effectiveBackwardLength' will store the length of runs
    // including the site base, if the site base extends them.
    uint64_t effectiveForwardLength = forwardRunLength;
    uint64_t effectiveBackwardLength = backwardRunLength;

    if (forwardBaseFound && forwardRunBaseVal == siteBaseValue) {
        effectiveForwardLength++; // Site extends the forward run
    } else if (backwardBaseFound && backwardRunBaseVal == siteBaseValue) {
        // Site extends the backward run (and not the forward one, due to 'else if')
        effectiveBackwardLength++;
    }

    if (effectiveForwardLength >= homopolymerThreshold || effectiveBackwardLength >= homopolymerThreshold) {
        return true;
    }

    // If the site base, the forward run base, and the backward run base are all identical,
    // and the total length of (run_strictly_left + site + run_strictly_right) meets threshold.
    if (forwardBaseFound && backwardBaseFound &&
        siteBaseValue == forwardRunBaseVal && 
        backwardRunBaseVal == forwardRunBaseVal) {
        // All three (site, char of forward run, char of backward run) are identical.
        // The length to check is the sum of the original strict run lengths plus 1 for the site itself.
        if ((forwardRunLength + backwardRunLength + 1) >= homopolymerThreshold) {
            return true;
        }
    }

    return false;
}





// Function to check for a specific strand bias pattern
// Returns true if the strand bias is detected, and populates outDominantStrandHet1 and outDominantStrandHet2.
// dominantStrand values: 0 (all on strand 0), 1 (all on strand 1), 2 (no dominance or not enough reads).
bool hasSiteStrandBias(
    const std::set<OrientedReadId>& hetBase1OrientedReads,
    const std::set<OrientedReadId>& hetBase2OrientedReads,
    uint64_t& outDominantStrandHet1,
    uint64_t& outDominantStrandHet2)
{
    // For this specific bias pattern, we need at least one read supporting each allele's strand pattern.
    const uint64_t minReadsPerAlleleForThisBias = 1;

    outDominantStrandHet1 = 2; // Initialize to no dominance
    if (hetBase1OrientedReads.size() >= minReadsPerAlleleForThisBias) {
        uint64_t hetBase1CountStrand0 = 0;
        uint64_t hetBase1CountStrand1 = 0;
        for (const auto& hetBase1OrientedRead : hetBase1OrientedReads) {
            if (hetBase1OrientedRead.getStrand() == 0) {
                hetBase1CountStrand0++;
            } else {
                hetBase1CountStrand1++;
            }
        }
        if (hetBase1CountStrand0 > 0 && hetBase1CountStrand1 == 0) {
            outDominantStrandHet1 = 0; // All reads for hetBase1 are on strand 0
        } else if (hetBase1CountStrand1 > 0 && hetBase1CountStrand0 == 0) {
            outDominantStrandHet1 = 1; // All reads for hetBase1 are on strand 1
        }
    }

    outDominantStrandHet2 = 2; // Initialize to no dominance
    if (hetBase2OrientedReads.size() >= minReadsPerAlleleForThisBias) {
        uint64_t hetBase2CountStrand0 = 0;
        uint64_t hetBase2CountStrand1 = 0;
        for (const auto& hetBase2OrientedRead : hetBase2OrientedReads) {
            if (hetBase2OrientedRead.getStrand() == 0) {
                hetBase2CountStrand0++;
            } else {
                hetBase2CountStrand1++;
            }
        }
        if (hetBase2CountStrand0 > 0 && hetBase2CountStrand1 == 0) {
            outDominantStrandHet2 = 0; // All reads for hetBase2 are on strand 0
        } else if (hetBase2CountStrand1 > 0 && hetBase2CountStrand0 == 0) {
            outDominantStrandHet2 = 1; // All reads for hetBase2 are on strand 1
        }
    }

    // Check for the specific bias pattern:
    // Both alleles must show complete strand dominance, and on opposite strands.
    if (outDominantStrandHet1 != 2 && outDominantStrandHet2 != 2 && outDominantStrandHet1 != outDominantStrandHet2) {
        return true; // Bias detected
    }

    return false; // No bias detected
}




class AlignmentPositionBaseStats{
    public:
        //
        // Fields for potential heterozygous site analysis
        //
        uint64_t positionInRead0 = 0;
        uint64_t baseOfReadId0 = 0; // Base value (0=A, 1=C, 2=G, 3=T, 4=Gap)

        uint64_t totalNumberOfA = 0;
        uint64_t totalNumberOfC = 0;
        uint64_t totalNumberOfG = 0;
        uint64_t totalNumberOfT = 0;
        uint64_t totalNumberOfGap = 0;
        uint64_t totalNumberOfAlignments = 0;

        std::set<OrientedReadId> orientedReadIdsWithA;
        std::set<OrientedReadId> orientedReadIdsWithC;
        std::set<OrientedReadId> orientedReadIdsWithG;
        std::set<OrientedReadId> orientedReadIdsWithT;
        std::set<OrientedReadId> orientedReadIdsWithGap;

        std::set<ReadId> readIdsWithA;
        std::set<ReadId> readIdsWithC;
        std::set<ReadId> readIdsWithG;
        std::set<ReadId> readIdsWithT;
        std::set<ReadId> readIdsWithGap;

        std::set<uint64_t> alignmentIdsWithA;
        std::set<uint64_t> alignmentIdsWithC;
        std::set<uint64_t> alignmentIdsWithG;
        std::set<uint64_t> alignmentIdsWithT;
        std::set<uint64_t> alignmentIdsWithGap;

        //
        // Fields for heterozygous site analysis
        //
        uint64_t hetBase1 = 100; // Invalid base value initially
        std::set<OrientedReadId> hetBase1OrientedReadIds;
        std::set<ReadId> hetBase1ReadIds;
        std::set<uint64_t> hetBase1AlignmentIds;
        uint64_t totalNumberOfHetBase1 = 0;
        double percentageOfHetBase1 = 0;

        uint64_t hetBase2 = 100; // Invalid base value initially
        std::set<OrientedReadId> hetBase2OrientedReadIds;
        std::set<ReadId> hetBase2ReadIds;
        std::set<uint64_t> hetBase2AlignmentIds;
        uint64_t totalNumberOfHetBase2 = 0;
        double percentageOfHetBase2 = 0;
};


class Site {
    public:
        std::set<OrientedReadId> orientedReads;
        OrientedReadId targetOrientedReadId;
        std::set<OrientedReadId> excludedOrientedReads;
    };

// Structure to hold data for the phasing threads
struct PhasingThreadData {
    // Pointers to assembler data (const access needed in thread)
    const Assembler *assembler; // Use const Assembler*
    uint64_t alignmentCount;
    uint64_t readCount; // Add readCount for checks
    uint64_t orientedReadCount; // Add orientedReadCount for checks

    // Parameters needed by the thread function
    // (Add any other parameters from createReadGraph4withStrandSeparation if needed inside the loop)

    // Thread-local storage for results
    // Each outer vector is indexed by threadId
    // Each inner vector<bool> is indexed by alignmentId
    vector<vector<bool> > threadForbiddenAlignments;
    vector<vector<bool> > threadFirstPassHetAlignments;
    vector<vector<bool> > threadAlignmentsAlreadyConsidered;
    vector<vector<Site> > threadSites;
    vector<bool> isReadIdContained;

    // Mutex for thread-safe cout if debugging is needed inside threads
    // std::mutex coutMutex; // Uncomment if needed

    PhasingThreadData(const Assembler *asmPtr, size_t threadCount) : assembler(asmPtr) {
        // --- Check asmPtr validity early ---
        SHASTA_ASSERT(assembler != nullptr);

        alignmentCount = assembler->alignmentData.size();
        readCount = assembler->getReads().readCount();
        orientedReadCount = readCount * 2;

        // --- Check threadCount validity ---
        SHASTA_ASSERT(threadCount > 0);

        threadForbiddenAlignments.resize(threadCount);
        threadFirstPassHetAlignments.resize(threadCount);
        threadAlignmentsAlreadyConsidered.resize(threadCount);
        threadSites.resize(threadCount);

        for (size_t i = 0; i < threadCount; ++i) {
            // Resize inner vectors
            threadForbiddenAlignments[i].resize(alignmentCount, false);
            threadFirstPassHetAlignments[i].resize(alignmentCount, false);
            threadAlignmentsAlreadyConsidered[i].resize(alignmentCount, false);
        }
    }
};

// Global pointer for thread access (similar pattern to computeAlignmentsData)
PhasingThreadData *phasingThreadData = nullptr;


// Thread function for the phasing analysis part
void Assembler::createReadGraph4PhasingThreadFunction(size_t threadId) {

    ofstream debugOut("Debug-" + to_string(threadId) + ".txt");
    debugOut << "Thread ID: " << threadId << " is starting read phasing" << endl;

    // --- Ensure phasingThreadData is valid ---
    SHASTA_ASSERT(phasingThreadData != nullptr);

    // Get reference to thread-local storage
    SHASTA_ASSERT(threadId < phasingThreadData->threadForbiddenAlignments.size()); // Check threadId validity
    vector<bool>& forbiddenAlignments = phasingThreadData->threadForbiddenAlignments[threadId];
    vector<bool>& firstPassHetAlignments = phasingThreadData->threadFirstPassHetAlignments[threadId];
    vector<bool>& alignmentsAlreadyConsidered = phasingThreadData->threadAlignmentsAlreadyConsidered[threadId];
    vector<Site>& sites = phasingThreadData->threadSites[threadId];
    const uint64_t alignmentCount = phasingThreadData->alignmentCount;
    const uint64_t orientedReadCount = phasingThreadData->orientedReadCount; // Get from shared data

    // Define types for the interval map
    using Interval = boost::icl::discrete_interval<uint32_t>;
    // Store pairs of (alignmentId, the other OrientedReadId relative to orientedReadId0)
    using AlignmentInfoPair = std::pair<uint64_t, OrientedReadId>;
    // Define a custom comparison for the set to handle pairs
    struct CompareAlignmentInfoPair {
        bool operator()(const AlignmentInfoPair& a, const AlignmentInfoPair& b) const {
            if (a.first != b.first) {
                return a.first < b.first; // Compare alignmentId first
            }
            return a.second < b.second; // Then compare OrientedReadId
        }
    };
    using AlignmentInfoSet = std::set<AlignmentInfoPair, CompareAlignmentInfoPair>;
    using AlignmentIntervalMap = boost::icl::interval_map<uint32_t, AlignmentInfoSet>;


    // Access assembler data via the pointer
    const Assembler& assembler = *(phasingThreadData->assembler);

    // Get batches of read IDs to process
    uint64_t readIdBegin;
    uint64_t readIdEnd; // variable to receive the end of the batch
    while (getNextBatch(readIdBegin, readIdEnd)) {

        debugOut << "Thread ID: " << threadId << " is processing read IDs from " << readIdBegin << " to " << readIdEnd << endl;

        // --- Check batch validity ---
        SHASTA_ASSERT(readIdBegin < phasingThreadData->readCount);

        // Ensure readIdEnd doesn't wrap around or go below readIdBegin after clamping
        SHASTA_ASSERT(readIdEnd >= readIdBegin);

        // Loop over read IDs in this batch
        for (ReadId readId = readIdBegin; readId < readIdEnd; readId++) {

            // --- Start of per-read analysis code ---
            // if ((readId != 2228) && (readId != 2229) && (readId != 2230) && (readId != 2231) && (readId != 2232) && (readId != 2233) && (readId != 2234) && (readId != 2235) && (readId != 2236) && (readId != 2237)) { // Keep this for debugging specific reads if needed
            //     continue;
            // }

            // if (((readId < 19300) || (readId > 19420))) { // Keep this for debugging specific reads if needed
            //     if ((readId < 1980) || (readId > 2050)) {
            //         continue;
            //     }
            // }

            // // --- Start of per-read analysis code ---
            // if (readId != 0) {
            //     continue;
            // }

            debugOut << "Thread ID: " << threadId << " is processing read ID: " << readId << endl;
            debugOut << timestamp << endl;

            const ReadId readId0 = readId;
            const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
            const OrientedReadId orientedReadId0(readId0, strand0);


            // // We do not perform the analysis on contained reads because the longer reads that contain them will
            // // make the decisions for them (where they best belong) given that they will have more informative
            // // sites to consider.
            // if (phasingThreadData->isReadIdContained[readId0]) {
            //     // Skip contained reads
            //     continue;
            // }



            // XXX
            // --- Build Interval Tree for alignments relative to orientedReadId0 ---
            //
            AlignmentIntervalMap alignmentIntervals; // Create the interval map for this readId0
            
            debugOut << timestamp << "Thread ID: " << threadId << " Building interval tree for read ID: " << readId0 << endl;

            const auto alignmentTableForIntervalTreeConstruction = alignmentTable[orientedReadId0.getValue()];
            for (const auto alignmentId : alignmentTableForIntervalTreeConstruction) {

                const AlignmentData& ad = alignmentData[alignmentId];

                OrientedReadId currentOrientedReadId0(ad.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = ad.info;

                Alignment alignment;
                span<const char> compressedAlignment = assembler.compressedAlignments[alignmentId];
                shasta::decompress(compressedAlignment, alignment);
            
                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != readId0) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                    alignment.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

                uint32_t markerCount0 = uint32_t(markers[currentOrientedReadId0.getValue()].size());
                uint32_t markerCount1 = uint32_t(markers[currentOrientedReadId1.getValue()].size());

                // Reverse complement if necessary
                if (currentOrientedReadId0.getStrand() != strand0) {
                    // Ensure marker counts are valid before calling reverseComplement
                    alignment.reverseComplement(markerCount0, markerCount1);
                    // alignment.checkStrictlyIncreasing(); // Can be expensive, maybe remove in production
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);
                SHASTA_ASSERT(currentOrientedReadId0 == orientedReadId0); // Check consistency

                // if ((alignmentInfo.data[0].leftTrim() <= 50) || (alignmentInfo.data[1].rightTrim() <= 50)) {
                //     continue;
                // }

                const auto orientedReadMarkers = markers[currentOrientedReadId0.getValue()];
                const uint32_t position0 = orientedReadMarkers[alignmentInfo.data[0].firstOrdinal].position;
                const uint32_t position1 = orientedReadMarkers[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assembler.assemblerInfo->k);

                // --- Add interval to the map ---
                if (position1 > position0) { // Ensure a valid interval
                    Interval interval = Interval::closed(position0, position1); // Create a [position0, position1) interval
                    // Create a pair with alignmentId and the correctly oriented other read
                    AlignmentInfoPair infoPair = {alignmentId, currentOrientedReadId1};
                    AlignmentInfoSet currentAlignmentInfoSet = {infoPair}; // Create a set with the pair
                    alignmentIntervals.add({interval, currentAlignmentInfoSet}); // Add to the map (handles overlaps)
                    debugOut << "Thread ID: " << threadId << " Added interval [" << position0 << ", " << position1 << ") for alignment " << alignmentId << " (Other Read: " << currentOrientedReadId1 << ")" << endl;
                } else {
                    // debugOut << "Thread ID: " << threadId << " Skipping invalid interval [" << position0 << ", " << position1 << "] for alignment " << alignmentId << endl;
                }
            }

            debugOut << timestamp << "Thread ID: " << threadId << " Finished building interval tree for read ID: " << readId0 << endl;

            // XXX
            // --- END OF: Build Interval Tree for alignments relative to orientedReadId0 ---
            //
            



            // // --- DEBUG: Query interval tree for position 35102 ---
            // uint32_t debugPosition = 35102;
            // auto it = alignmentIntervals.find(debugPosition);

            // if (it != alignmentIntervals.end()) {
            //     debugOut << "--- DEBUG INFO: Alignments covering position " << debugPosition << " in read " << orientedReadId0 << " ---" << endl;
            //     // The value type is now AlignmentInfoSet
            //     const AlignmentInfoSet& coveringAlignmentInfo = it->second;

            //     // Get the count directly from the size of the set
            //     uint64_t alignmentCountAtPosition = coveringAlignmentInfo.size();
            //     debugOut << "  Number of alignments covering position: " << alignmentCountAtPosition << endl;

            //     debugOut << "  Alignment Info (ID, OtherRead): " << endl;
            //     // Iterate through the pairs in the set
            //     for (const AlignmentInfoPair& infoPair : coveringAlignmentInfo) {
            //         uint64_t alignmentId = infoPair.first;
            //         OrientedReadId storedOrientedReadId1 = infoPair.second; // Get the stored OrientedReadId

            //         // Print the alignment ID and the stored OrientedReadId1
            //         debugOut << "(" << alignmentId << ", " << storedOrientedReadId1 << ") " << endl;
            //     }
            //     debugOut << endl;
            //     debugOut << "--- END DEBUG INFO for position " << debugPosition << " ---" << endl;
            // } else {
            //     debugOut << "--- DEBUG INFO: No alignments found covering position " << debugPosition << " in read " << orientedReadId0 << " ---" << endl;
            // }
            // // --- END DEBUG ---










            debugOut << timestamp << "Thread ID: " << threadId << " Starting looping over alignments for read ID: " << readId0 << endl;

            std::map<uint64_t, AlignmentPositionBaseStats> positionStatsOnOrientedReadId0;
            std::map<uint64_t, AlignmentPositionBaseStats> potentialHetSitesOnOrientedReadId0;




            // Loop over alignments involving the target read (readId0)
            const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];

            
            // XXX
            // --- Find and forbid alignments between the same read on different strands.
            //
            vector <bool> readIdsUsed(phasingThreadData->readCount, false);
            vector <bool> forbiddenReadIds(phasingThreadData->readCount, false);

            for (const auto alignmentId : alignmentTable0) {
                
                const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
                OrientedReadId orientedReadId0(thisAlignmentData.readIds[0], 0);
                OrientedReadId orientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thisAlignmentData.info;

            
                const ReadId currentReadId0 = thisAlignmentData.readIds[0];
                const ReadId currentReadId1 = thisAlignmentData.readIds[1];
                if (currentReadId0 == readId0) {
                    if(readIdsUsed[currentReadId1]) {
                        forbiddenAlignments[alignmentId] = true;
                        forbiddenReadIds[currentReadId1] = true;
                        continue;
                    }
                    readIdsUsed[currentReadId1] = true;
                } else if (currentReadId1 == readId0) {
                    if(readIdsUsed[currentReadId0]) {
                        forbiddenAlignments[alignmentId] = true;
                        forbiddenReadIds[currentReadId0] = true;
                        continue;
                    }
                    readIdsUsed[currentReadId0] = true;                
                }
                
            }

            // We need a second loop over the alignments to forbid the second alignment instance in some cases
            for (const auto alignmentId : alignmentTable0) {
                
                const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
                const ReadId currentReadId0 = thisAlignmentData.readIds[0];
                const ReadId currentReadId1 = thisAlignmentData.readIds[1];
                if (currentReadId0 == readId0) {
                    if(forbiddenReadIds[currentReadId1]) {
                        forbiddenAlignments[alignmentId] = true;
                    }
                } else if (currentReadId1 == readId0) {
                    if(forbiddenReadIds[currentReadId0]) {
                        forbiddenAlignments[alignmentId] = true;
                    }               
                }

            }

            // XXX
            // --- END OF: Find and forbid alignments between the same read on different strands.
            //

            

            // XXX
            // --- We found and forbid alignments between the same read on different strands.
            //     Now we can proceed with the analysis of the rest of alignments.

            for (const auto alignmentId : alignmentTable0) {

                if (forbiddenAlignments[alignmentId]) {
                    // Skip alignments that are forbidden due to strand issues
                    continue;
                }

                const AlignmentData& ad = alignmentData[alignmentId];
                const size_t alignedMarkerCount = ad.info.markerCount;

                OrientedReadId currentOrientedReadId0(ad.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = ad.info;

                Alignment alignment;
                span<const char> compressedAlignment = assembler.compressedAlignments[alignmentId];
                shasta::decompress(compressedAlignment, alignment);

                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != readId0) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                    alignment.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

                uint32_t markerCount0 = uint32_t(markers[currentOrientedReadId0.getValue()].size());
                uint32_t markerCount1 = uint32_t(markers[currentOrientedReadId1.getValue()].size());

                // Reverse complement if necessary
                if (currentOrientedReadId0.getStrand() != strand0) {
                    // Ensure marker counts are valid before calling reverseComplement
                    alignment.reverseComplement(markerCount0, markerCount1);
                    // alignment.checkStrictlyIncreasing(); // Can be expensive, maybe remove in production
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);





                // if ((alignmentInfo.data[0].leftTrim() <= 50) || (alignmentInfo.data[1].rightTrim() <= 50)) {
                //     continue;
                // }





                // // Explore only alignments that are not hanging on the left.
                // // The target read is always in the most "left" position
                // if ((ad.info.data[0].leftTrim() <= 200)) {
                //     continue;
                // }


                // if (phasingThreadData->isReadIdContained[currentOrientedReadId1.getReadId()]) {
                //     // Skip contained reads
                //     continue;
                // }


                // if (ad.info.data[1].rightTrim() <= 50) {
                //     // Skip alignments that are not right trimmed
                //     continue;
                // }

                // if ((ad.info.data[1].rightTrim() <= 50) || (ad.info.offsetAtCenter() <= 0.)) {
                //      // Skip alignments that are not right trimmed
                //      continue;
                // }




                const auto orientedReadMarkers = markers[currentOrientedReadId0.getValue()];
                const uint32_t position0 = orientedReadMarkers[alignmentInfo.data[0].firstOrdinal].position;
                const uint32_t position1 = orientedReadMarkers[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assembler.assemblerInfo->k);

                
                // XXX
                // --- Project this alignment to base space.
                //

                ProjectedAlignment projectedAlignment(
                    assembler,
                    {currentOrientedReadId0, currentOrientedReadId1},
                    alignment,
                    ProjectedAlignment::Method::QuickRaw
                );

                // XXX
                // --- END OF: Project this alignment to base space.
                //


                // ofstream html("a.html");
                // projectedAlignment.writeHtml(html, false);
                // SHASTA_ASSERT(0);

                // debugOut << "Projected alignment between ReadId: " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and ReadId: " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " completed." << endl;
                
                // Loop over the RAW segments of the projectedAlignment.
                for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
                    
                    // Get the RAW sequences
                    const vector<Base>& sequence0 = segment.sequences[0];
                    const vector<Base>& sequence1 = segment.sequences[1];

                    // Align them base by base to get the right sequence alignment representation.
                    uint64_t position0 = 0;
                    uint64_t position1 = 0;
                    vector<AlignedBase> rawAlignmentSequence0;
                    vector<AlignedBase> rawAlignmentSequence1;

                    // Retreive the RAW alignment sequences.
                    for (const pair<bool, bool>& p: segment.alignment) {
                        const bool hasBase0 = p.first;
                        const bool hasBase1 = p.second;
                        
                        AlignedBase baseAtPosition0;
                        AlignedBase baseAtPosition1;
                        uint64_t positionWithDetectedChange0 = 0;
                        uint64_t positionWithDetectedChange1 = 0;
                        if(hasBase0) {
                            SHASTA_ASSERT(position0 < sequence0.size());
                            baseAtPosition0 = AlignedBase(sequence0[position0]);
                            rawAlignmentSequence0.push_back(AlignedBase(sequence0[position0]));
                            positionWithDetectedChange0 = position0;
                            position0++;
                        } else {
                            rawAlignmentSequence0.push_back(AlignedBase::gap());
                            baseAtPosition0 = AlignedBase::gap();
                        }

                        if(hasBase1) {
                            SHASTA_ASSERT(position1 < sequence1.size());
                            baseAtPosition1 = AlignedBase(sequence1[position1]);
                            rawAlignmentSequence1.push_back(AlignedBase(sequence1[position1]));
                            positionWithDetectedChange1 = position1;
                            position1++;
                        } else {
                            rawAlignmentSequence1.push_back(AlignedBase::gap());
                            baseAtPosition1 = AlignedBase::gap();
                        }

                        // We need to keep track of the gaps even though we do not rely on them for phasing.
                        // The number of alignments with gaps in this position is important to infer statistics
                        // about this position later (like the total alignments in that position).
                        if (hasBase0 && (baseAtPosition0 != baseAtPosition1)) {
                            


                            // //
                            // // --- DEBUG: Print the RAW alignment sequences that have differences ---
                            // //            AND the differences are not gaps
                            // //
                            // if(hasBase1) {
                            //     debugOut << endl;
                            //     debugOut << "rawAlignmentSequence0 marker positionsA start: " << segment.positionsA[0] << endl;
                            //     debugOut << "rawAlignmentSequence1 marker positionsA start: " << segment.positionsA[1] << endl;
                            //     debugOut << "sequence0 base change positionsA start: " << segment.positionsA[0] + positionWithDetectedChange0 << endl;
                            //     debugOut << "sequence1 base change positionsA start: " << segment.positionsA[1] + positionWithDetectedChange1 << endl;
                            //     debugOut << "sequence0 base at base change positionsA: " << baseAtPosition0 << endl;
                            //     debugOut << "sequence1 base at base change positionsA: " << baseAtPosition1 << endl;
                            //     debugOut << "We found " << segment.editDistance << " edit distance differences between the RAW sequences" << endl;
                            //     debugOut << "rawAlignmentSequence0: ";
                            //     for(uint64_t i=0; i<rawAlignmentSequence0.size(); i++) {
                            //         const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
                            //         if(isDifferent) {
                            //             debugOut << "[";
                            //         }
                            //         debugOut << rawAlignmentSequence0[i].character();
                            //         if(isDifferent) {
                            //             debugOut << "]";
                            //         }
                            //     }
                            //     debugOut << endl;

                            //     debugOut << "rawAlignmentSequence1: ";
                            //     for(uint64_t i=0; i<rawAlignmentSequence1.size(); i++) {
                            //         const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
                            //         if(isDifferent) {
                            //             debugOut << "[";
                            //         }
                            //         debugOut << rawAlignmentSequence1[i].character();
                            //         if(isDifferent) {
                            //             debugOut << "]";
                            //         }
                            //     }
                            //     debugOut << endl;
                            // }
                            // //
                            // // --- END DEBUG ---
                            // //


                            // XXX
                            // --- Update statistics for positionInRead0 ---
                            //

                            uint64_t positionInRead0 = segment.positionsA[0] + positionWithDetectedChange0;
                            uint64_t positionInRead1 = segment.positionsA[1] + positionWithDetectedChange1;
                            // Initialize stats for this position if not already present
                            if(not positionStatsOnOrientedReadId0.contains(positionInRead0)) {
                                AlignmentPositionBaseStats thisPositionStats;
                                thisPositionStats.positionInRead0 = positionInRead0;
                                thisPositionStats.baseOfReadId0 = baseAtPosition0.value; 
                                positionStatsOnOrientedReadId0.emplace(positionInRead0, thisPositionStats);
                            }

                            // Update counts based on rawAlignmentSequence1[i]
                            auto& positionStatsInRead0 = positionStatsOnOrientedReadId0[positionInRead0];

                            // // Increment alignment count for this position
                            // // if (hasBase0 && hasBase1 && (baseAtPosition0 != baseAtPosition1)) {
                            // if (hasBase0 && hasBase1 && (baseAtPosition0 != baseAtPosition1)) {
                            //     positionStatsInRead0.totalNumberOfAlignments++;
                            // }
                                           
                            if(baseAtPosition1.value == 0) { // A in sequence1
                                positionStatsInRead0.totalNumberOfA++;
                                positionStatsInRead0.orientedReadIdsWithA.insert(currentOrientedReadId1);
                                positionStatsInRead0.readIdsWithA.insert(currentOrientedReadId1.getReadId());
                                positionStatsInRead0.alignmentIdsWithA.insert(alignmentId);
                            } else if(baseAtPosition1.value == 1) { // C in sequence1
                                positionStatsInRead0.totalNumberOfC++;
                                positionStatsInRead0.orientedReadIdsWithC.insert(currentOrientedReadId1);
                                positionStatsInRead0.readIdsWithC.insert(currentOrientedReadId1.getReadId());
                                positionStatsInRead0.alignmentIdsWithC.insert(alignmentId);
                            } else if(baseAtPosition1.value == 2) { // G in sequence1
                                positionStatsInRead0.totalNumberOfG++;
                                positionStatsInRead0.orientedReadIdsWithG.insert(currentOrientedReadId1);
                                positionStatsInRead0.readIdsWithG.insert(currentOrientedReadId1.getReadId());
                                positionStatsInRead0.alignmentIdsWithG.insert(alignmentId);
                            } else if(baseAtPosition1.value == 3) { // T in sequence1
                                positionStatsInRead0.totalNumberOfT++;
                                positionStatsInRead0.orientedReadIdsWithT.insert(currentOrientedReadId1);
                                positionStatsInRead0.readIdsWithT.insert(currentOrientedReadId1.getReadId());
                                positionStatsInRead0.alignmentIdsWithT.insert(alignmentId);
                            } else if(baseAtPosition1.value == 4) { // Gap in sequence1
                                positionStatsInRead0.totalNumberOfGap++;
                                positionStatsInRead0.orientedReadIdsWithGap.insert(currentOrientedReadId1);
                                positionStatsInRead0.readIdsWithGap.insert(currentOrientedReadId1.getReadId());
                                positionStatsInRead0.alignmentIdsWithGap.insert(alignmentId);
                            }

                        } // --- End of printing the sequences that have differences and updating the stats of this position ---

                    } // --- End of loop over the segment ---
                    SHASTA_ASSERT(rawAlignmentSequence0.size() == rawAlignmentSequence1.size());

                } // --- End loop over segments ---

                // debugOut << "Loop over the RAW segments of the Projected alignment between ReadId: " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and ReadId: " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " completed." << endl;  
            
            } // --- End loop over alignments for readId0 ---

            debugOut << timestamp << "Thread ID: " << threadId << " Finished looping over alignments for read ID: " << readId0 << endl;

            
            
            // XXX
            // --- Analyze potential heterozygous sites ---
            //

            // We finished analyzing all alignments for the target read (readId0).
            // Now we need to check each potential site in readId0 (in positionStatsOnOrientedReadId0) 
            // to see if it involves a heterozygous site.
            uint64_t sitesSkippedDueToInsufficientCoverage = 0;
            for (auto const& [positionInRead0, positionStatsInRead0] : positionStatsOnOrientedReadId0) {

                // if (positionInRead0 == 35102) {
                //     cout << "positionInRead0: " << positionInRead0 << endl;
                //     cout << "positionStatsInRead0.totalNumberOfAlignments: " << positionStatsInRead0.totalNumberOfAlignments << endl;
                //     cout << "positionStatsInRead0.totalNumberOfA: " << positionStatsInRead0.totalNumberOfA << endl;
                //     cout << "positionStatsInRead0.totalNumberOfC: " << positionStatsInRead0.totalNumberOfC << endl;
                //     cout << "positionStatsInRead0.totalNumberOfG: " << positionStatsInRead0.totalNumberOfG << endl;
                //     cout << "positionStatsInRead0.totalNumberOfT: " << positionStatsInRead0.totalNumberOfT << endl;
                //     cout << "positionStatsInRead0.totalNumberOfGap: " << positionStatsInRead0.totalNumberOfGap << endl;
                //     cout << "positionStatsInRead0.baseOfReadId0: " << positionStatsInRead0.baseOfReadId0 << endl;
                //     cout << "positionStatsInRead0.orientedReadIdsWithA.size(): " << positionStatsInRead0.orientedReadIdsWithA.size() << endl;
                //     cout << "positionStatsInRead0.orientedReadIdsWithC.size(): " << positionStatsInRead0.orientedReadIdsWithC.size() << endl;
                //     cout << "positionStatsInRead0.orientedReadIdsWithG.size(): " << positionStatsInRead0.orientedReadIdsWithG.size() << endl;
                //     cout << "positionStatsInRead0.orientedReadIdsWithT.size(): " << positionStatsInRead0.orientedReadIdsWithT.size() << endl;
                //     cout << "positionStatsInRead0.orientedReadIdsWithGap.size(): " << positionStatsInRead0.orientedReadIdsWithGap.size() << endl;
                //     cout << "positionStatsInRead0.readIdsWithA.size(): " << positionStatsInRead0.readIdsWithA.size() << endl;
                //     cout << "positionStatsInRead0.readIdsWithC.size(): " << positionStatsInRead0.readIdsWithC.size() << endl;
                //     cout << "positionStatsInRead0.readIdsWithG.size(): " << positionStatsInRead0.readIdsWithG.size() << endl;
                //     cout << "positionStatsInRead0.readIdsWithT.size(): " << positionStatsInRead0.readIdsWithT.size() << endl;
                //     cout << "positionStatsInRead0.readIdsWithGap.size(): " << positionStatsInRead0.readIdsWithGap.size() << endl;
                //     cout << "positionStatsInRead0.alignmentIdsWithA.size(): " << positionStatsInRead0.alignmentIdsWithA.size() << endl;
                //     cout << "positionStatsInRead0.alignmentIdsWithC.size(): " << positionStatsInRead0.alignmentIdsWithC.size() << endl;
                //     cout << "positionStatsInRead0.alignmentIdsWithG.size(): " << positionStatsInRead0.alignmentIdsWithG.size() << endl;
                //     cout << "positionStatsInRead0.alignmentIdsWithT.size(): " << positionStatsInRead0.alignmentIdsWithT.size() << endl;
                //     cout << "positionStatsInRead0.alignmentIdsWithGap.size(): " << positionStatsInRead0.alignmentIdsWithGap.size() << endl;
                // }






                // XXX
                // --- Check if this potential het site in this position ---
                //     is in a homopolymer region.
                //

                ShortBaseSequence8 kmer;

                const LongBaseSequenceView readId0Sequence = reads->getRead(readId0);

                // Extract the kmer from the readId0 sequence that includes
                // the 2 preceding and the 2 following bases on that position.
                extractKmer(readId0Sequence, positionInRead0-3, 7, kmer);
                
                uint64_t homopolymerThreshold = 3;
                bool siteIsInHomopolymerRegion = isSiteInHomopolymerRegion(positionInRead0, readId0Sequence, homopolymerThreshold);
                if (siteIsInHomopolymerRegion) {
                    // Skip this site because it is in a homopolymer region
                    sitesSkippedDueToInsufficientCoverage++;
                    // debugOut << "Skipping position " << positionInRead0 << " due to homopolymer region." << std::endl;
                    // kmer.write(debugOut, 7) << std::endl;
                    continue;
                }

                // XXX
                // --- END OF: Check if this potential het site in this position ---
                //             is in a homopolymer region.
                //


                


                // Get the number of alignment supporting each base.
                const uint64_t totalNumberOfA = positionStatsInRead0.totalNumberOfA;
                const uint64_t totalNumberOfC = positionStatsInRead0.totalNumberOfC;
                const uint64_t totalNumberOfG = positionStatsInRead0.totalNumberOfG;
                const uint64_t totalNumberOfT = positionStatsInRead0.totalNumberOfT;
                const uint64_t totalNumberOfGap = positionStatsInRead0.totalNumberOfGap;

                uint64_t totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead =
                    totalNumberOfA +
                    totalNumberOfC +
                    totalNumberOfG +
                    totalNumberOfT +
                    totalNumberOfGap;

                uint64_t totalNumberOfAlignmentsWithoutGapsSupportingAChangeFromTheTargetRead = totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead - totalNumberOfGap;
                
                // This is the minimum needed number of alignments that
                // support a change from the base of the target read in that position
                uint64_t minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead = 3;

                // Check if there is a supporting base change with sufficient coverage
                if ((totalNumberOfA < minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead) &&
                    (totalNumberOfC < minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead) &&
                    (totalNumberOfG < minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead) &&
                    (totalNumberOfT < minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead)) {
                    sitesSkippedDueToInsufficientCoverage++;
                    // if (positionInRead0 == 35102) {
                    //     cout << "THIS SHOULD NOT HAPPEN. Skipping position " << positionInRead0 << " due to insufficient coverage: " << totalNumberOfAlignments << endl;
                    // }
                    // debugOut << "Skipping position " << positionInRead0 << " due to insufficient coverage supporting a change from the target read: " << totalNumberOfAlignmentsWithoutGapsSupportingAChangeFromTheTargetRead << endl;
                    // debugOut << "Position: " << positionInRead0 << " A: " << totalNumberOfA << " C: " << totalNumberOfC << " G: " << totalNumberOfG << " T: " << totalNumberOfT << " Gap: " << totalNumberOfGap << endl;
                    continue;
                }
                
                // Sort the base counts in descending order
                vector<pair<uint64_t, uint64_t>> baseCounts = {
                    {totalNumberOfA, 0},
                    {totalNumberOfC, 1},
                    {totalNumberOfG, 2},
                    {totalNumberOfT, 3},
                    {totalNumberOfGap, 4}
                };
                std::sort(baseCounts.begin(), baseCounts.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

                
                // Check if this is a potential heterozygous site
                // Criteria: The 2 alleles with the highest coverage are above 4 coverage each 
                // and the target read supports the base of one of them
                const uint64_t highestCoverageBaseCounts = baseCounts[0].first;
                const uint64_t secondHighestCoverageBaseCounts = baseCounts[1].first;
                
                // If the highest coverage base has coverage 3 or more (minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead) 
                // AND the base is not a gap, then proceede.
                if ((highestCoverageBaseCounts >= minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead) and (baseCounts[0].second != 4)) {

                    // If the second highest coverage base has coverage more than 0 then skip this site for now.
                    // This site might be a complex one with a second allele base difference
                    // or it might include alignments with gaps (which make the site suspicious)
                    // TODO: If there is a second allele base difference (complex site)
                    // then skip this site for now
                    if ((secondHighestCoverageBaseCounts > 0) and (baseCounts[1].second != 4)) {
                        // debugOut << "Skipping position " << positionInRead0 << " due to having a second highest coverage base that is not a gap: " << secondHighestCoverageBaseCounts << " (complex site)." << endl;
                        // debugOut << "Position: " << positionInRead0 << " A: " << totalNumberOfA << " C: " << totalNumberOfC << " G: " << totalNumberOfG << " T: " << totalNumberOfT << " Gap: " << totalNumberOfGap << endl;
                        continue;
                    }

                    //
                    // This is a potential heterozygous site
                    //

                    // We don't have the number of alignments that agree with the base of the target read
                    // in that specific position. We only have the number of alignments that support a change.
                    // So, we need to find the total coverage of that position from the intervalTree.
                    uint64_t totalCoverageInThisPosition = 0;
                    std::set<OrientedReadId> coveringOrientedReadIds;
                    std::set<uint64_t> coveringAlignmentIds;
                    
                    auto it = alignmentIntervals.find(positionInRead0);
                    if (it != alignmentIntervals.end()) {
                        
                        // Retrieve the covering alignment info pairs of this position,
                        // which include the pairs (alignmentId, orientedReadId of the other read).
                        const AlignmentInfoSet& coveringAlignmentInfo = it->second;
                        
                        // Calculate the total coverage in this position. 
                        totalCoverageInThisPosition = coveringAlignmentInfo.size();

                        // Iterate through the pairs in coveringAlignmentInfo to
                        // extract the alignmentId and the orientedReadId of the other read.
                        for (const AlignmentInfoPair& infoPair : coveringAlignmentInfo) {
                            uint64_t alignmentId = infoPair.first;
                            OrientedReadId storedOrientedReadId1 = infoPair.second;
                            if (positionInRead0 == 513) {
                                // debugOut << "Found covering alignment ID: " << alignmentId << " for position " << positionInRead0 << endl;
                                // debugOut << "OrientedReadId1: " << storedOrientedReadId1 << endl;
                            }
                            // if (phasingThreadData->isReadIdContained[storedOrientedReadId1.getReadId()]) {
                            //     // Skip contained reads
                            //     continue;
                            // }
                            coveringAlignmentIds.insert(alignmentId);
                            coveringOrientedReadIds.insert(storedOrientedReadId1);
                        }

                    }

                    // Calculate the total coverage in this position without including alignments with gaps
                    uint64_t totalCoverageInThisPositionWithoutGaps = totalCoverageInThisPosition - positionStatsInRead0.totalNumberOfGap;

                    // This is the minimum needed number of alignments that
                    // support the target read in that position. 
                    // This is 2 + 1(the target read) = 3
                    uint64_t minimumNumberOfAlignmentsSupportingTheTargetRead = 2;

                    // Skip this position if the coverage supporting the target read is less than 2
                    if (totalCoverageInThisPosition - totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead < minimumNumberOfAlignmentsSupportingTheTargetRead) {
                        sitesSkippedDueToInsufficientCoverage++;
                        // debugOut << "Skipping position " << positionInRead0 << " due to insufficient coverage supporting the target read: " << totalCoverageInThisPosition - totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead << endl;
                        // debugOut << "Position: " << positionInRead0 << " A: " << totalNumberOfA << " C: " << totalNumberOfC << " G: " << totalNumberOfG << " T: " << totalNumberOfT << " Gap: " << totalNumberOfGap << " Total Alignments that support a change: " << totalNumberOfAlignmentsWithoutGapsSupportingAChangeFromTheTargetRead << ". Total alignments: " << totalCoverageInThisPosition << ". Total alignments without gaps: " << totalCoverageInThisPositionWithoutGaps << "." << endl;
                        continue;
                    }

                    // Print position and base counts
                    // debugOut << "Position: " << positionInRead0 << " A: " << totalNumberOfA << " C: " << totalNumberOfC << " G: " << totalNumberOfG << " T: " << totalNumberOfT << " Gap: " << totalNumberOfGap << " Total Alignments that support a change: " << totalNumberOfAlignmentsWithoutGapsSupportingAChangeFromTheTargetRead << ". Total alignments: " << totalCoverageInThisPosition << ". Total alignments without gaps: " << totalCoverageInThisPositionWithoutGaps << "." << endl;
                    

                    // if (positionInRead0 == 35102) {
                    //     cout << "Position " << positionInRead0 << " is a potential heterozygous site with " << totalCoverageInThisPosition << " total alignments and " << positionStatsInRead0.totalNumberOfAlignments << " alignments supporting the other allele." << endl;
                    // }

                    
                    const double supportingTheTargetReadBaseCountsPercentage = double(totalCoverageInThisPositionWithoutGaps - highestCoverageBaseCounts) / double(totalCoverageInThisPositionWithoutGaps);
                    const double supportingTheHetReadBaseCountsPercentage = double(highestCoverageBaseCounts) / double(totalCoverageInThisPositionWithoutGaps);


                    potentialHetSitesOnOrientedReadId0[positionInRead0] = positionStatsInRead0;
                    // hetBase1 will always be the base of the target read in that position
                    potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1 = positionStatsInRead0.baseOfReadId0;
                    potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2 = baseCounts[0].second;

                    // if (positionInRead0 == 35102) {
                    //     cout << "In position " << positionInRead0 << ", the base of the readId0 is one of the top two het bases." << endl;
                    //     cout << "hetBase1 (The base of the target read): " << potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1 << endl;
                    //     cout << "hetBase2: " << potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2 << endl;
                    // }

                    auto& orientedReadIdsSupportingTheTargetReadBase = potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds;
                    auto& orientedReadIdsSupportingTheOtherBase = potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds;
                    if (baseCounts[0].second == 0) { // A is the 'other' base
                        orientedReadIdsSupportingTheOtherBase = positionStatsInRead0.orientedReadIdsWithA;

                        // Calculate the difference: coveringOrientedReadIds - orientedReadIdsSupportingTheOtherBase
                        std::set_difference(
                            coveringOrientedReadIds.begin(), coveringOrientedReadIds.end(),
                            orientedReadIdsSupportingTheOtherBase.begin(), orientedReadIdsSupportingTheOtherBase.end(),
                            std::inserter(orientedReadIdsSupportingTheTargetReadBase, orientedReadIdsSupportingTheTargetReadBase.begin())
                        );

                    } else if (baseCounts[0].second == 1) { // C is the 'other' base
                        orientedReadIdsSupportingTheOtherBase = positionStatsInRead0.orientedReadIdsWithC;
                    
                        std::set_difference(
                            coveringOrientedReadIds.begin(), coveringOrientedReadIds.end(),
                            orientedReadIdsSupportingTheOtherBase.begin(), orientedReadIdsSupportingTheOtherBase.end(),
                            std::inserter(orientedReadIdsSupportingTheTargetReadBase, orientedReadIdsSupportingTheTargetReadBase.begin())
                        );
                    
                    } else if (baseCounts[0].second == 2) { // G is the 'other' base
                        orientedReadIdsSupportingTheOtherBase = positionStatsInRead0.orientedReadIdsWithG;
                    
                        std::set_difference(
                            coveringOrientedReadIds.begin(), coveringOrientedReadIds.end(),
                            orientedReadIdsSupportingTheOtherBase.begin(), orientedReadIdsSupportingTheOtherBase.end(),
                            std::inserter(orientedReadIdsSupportingTheTargetReadBase, orientedReadIdsSupportingTheTargetReadBase.begin())
                        );
                    
                    } else if (baseCounts[0].second == 3) { // T is the 'other' base
                        orientedReadIdsSupportingTheOtherBase = positionStatsInRead0.orientedReadIdsWithT;
                    
                        std::set_difference(
                            coveringOrientedReadIds.begin(), coveringOrientedReadIds.end(),
                            orientedReadIdsSupportingTheOtherBase.begin(), orientedReadIdsSupportingTheOtherBase.end(),
                            std::inserter(orientedReadIdsSupportingTheTargetReadBase, orientedReadIdsSupportingTheTargetReadBase.begin())
                        );
                    
                    } else if (baseCounts[0].second == 4) { // Gap is the 'other' base
                        orientedReadIdsSupportingTheOtherBase = positionStatsInRead0.orientedReadIdsWithGap;
                    
                        std::set_difference(
                            coveringOrientedReadIds.begin(), coveringOrientedReadIds.end(),
                            orientedReadIdsSupportingTheOtherBase.begin(), orientedReadIdsSupportingTheOtherBase.end(),
                            std::inserter(orientedReadIdsSupportingTheTargetReadBase, orientedReadIdsSupportingTheTargetReadBase.begin())
                        );
                    }

                    // There are cases where the coveringOrientedReadIds also include orientedReadIds that support a gap.
                    // We should remove those alignments to get the accurate set of orientedReadIds that support the target read base.
                    std::set<OrientedReadId> tempResult;
                    std::set_difference(
                        orientedReadIdsSupportingTheTargetReadBase.begin(), orientedReadIdsSupportingTheTargetReadBase.end(),
                        positionStatsInRead0.orientedReadIdsWithGap.begin(), positionStatsInRead0.orientedReadIdsWithGap.end(),
                        std::inserter(tempResult, tempResult.begin())
                    );
                    orientedReadIdsSupportingTheTargetReadBase.swap(tempResult);

                    // if (positionInRead0 == 35102) {
                    //     // print which orientedReadIds are supporting the target read base
                    //     cout << "OrientedReadIds supporting the target read base: " << endl;
                    //     for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheTargetReadBase) {
                    //         cout << orientedReadId << endl;
                    //     }
                    //     // print which orientedReadIds are supporting the other base
                    //     cout << "OrientedReadIds supporting the other base: " << endl;
                    //     for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheOtherBase) {
                    //         cout << orientedReadId << endl;
                    //     }
                    // }
                    
                    potentialHetSitesOnOrientedReadId0[positionInRead0].totalNumberOfHetBase1 = (totalCoverageInThisPositionWithoutGaps - highestCoverageBaseCounts);
                    potentialHetSitesOnOrientedReadId0[positionInRead0].totalNumberOfHetBase2 = highestCoverageBaseCounts;
                    potentialHetSitesOnOrientedReadId0[positionInRead0].percentageOfHetBase1 = supportingTheTargetReadBaseCountsPercentage;
                    potentialHetSitesOnOrientedReadId0[positionInRead0].percentageOfHetBase2 = supportingTheHetReadBaseCountsPercentage;
                    
                    // // // Optional: Log information about this site
                    // char bases[] = {'A', 'C', 'G', 'T', '-'};
                    // cout << "Potential heterozygous site at position " << positionInRead0 
                    //     << " in readId " << readId0 << ":" << endl
                    //     << "  Base of readId0: " << bases[positionStatsInRead0.baseOfReadId0] << endl
                    //     << "  First variant: " << bases[baseCounts[0].second] << " (" << firstBasePercentage * 100 << "%)" << endl
                    //     << "  Second variant: " << bases[baseCounts[1].second] << " (" << secondBasePercentage * 100 << "%)" << endl;
                    




                    // XXX
                    // --- Check for strand bias ---
                    // This site is strand biased if:
                    // 1. All reads supporting the target read's base (hetBase1) are on one specific strand.
                    // 2. AND all reads supporting the other heterozygous base (hetBase2) are on the *opposite* strand.
                    //

                    uint64_t dominantStrandHet1Value;
                    uint64_t dominantStrandHet2Value;

                    bool strandBiasDetected = hasSiteStrandBias(
                        potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds,
                        potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds,
                        dominantStrandHet1Value,
                        dominantStrandHet2Value
                    );

                    if (strandBiasDetected) {
                        // debugOut << "Strand bias detected at position " << positionInRead0 
                        //            << " for readId " << readId0 << ". Reads for target base on strand " << dominantStrandHet1Value
                        //            << ", reads for other base on strand " << dominantStrandHet2Value
                        //            << ". Skipping this site." << endl;
                        sitesSkippedDueToInsufficientCoverage++; // Using existing counter, consider a new one for clarity if needed
                        potentialHetSitesOnOrientedReadId0.erase(positionInRead0);
                        continue; // Skip to the next positionStatsInRead0
                    }

                    // XXX
                    // --- END OF: Check if we have strand bias issues in this position ---
                    //

                } // --- End of loop over a specific potential heterozygous site ---

            } // --- End of loop over all potential heterozygous sites ---
            



            // XXX
            // --- Now we need to analyze the potential heterozygous sites and try to find sets of reads that belong to the same haplotype ---
            //

            if (potentialHetSitesOnOrientedReadId0.size() > 0) {
                debugOut << endl;
                debugOut << "Found " << potentialHetSitesOnOrientedReadId0.size() << " potential heterozygous sites in readId " << readId0 << endl;
                debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
            } else {
                debugOut << endl;
                debugOut << "Found no potential heterozygous sites in readId " << readId0 << endl;
                debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
                continue;
            }

            // Loop over the potential heterozygous sites and print them
            char bases[] = {'A', 'C', 'G', 'T', '-'};
            for (const auto& [positionInRead0, positionStatsInRead0] : potentialHetSitesOnOrientedReadId0) {
                if (positionStatsInRead0.hetBase1 == 4 || positionStatsInRead0.hetBase2 == 4) {
                    // Skip gaps
                    continue;
                }
                debugOut << "Potential heterozygous site at position " << positionInRead0
                    << " in readId " << readId0 << " in strand " << strand0 << ":" << endl
                    << "  Base of readId0: " << bases[positionStatsInRead0.baseOfReadId0] << endl
                    << "  First variant: " << bases[positionStatsInRead0.hetBase1] << " (" << positionStatsInRead0.totalNumberOfHetBase1 << " "<< positionStatsInRead0.percentageOfHetBase1 * 100 << "%)" << endl
                    << "  Second variant: " << bases[positionStatsInRead0.hetBase2] << " (" << positionStatsInRead0.totalNumberOfHetBase2 << " " << positionStatsInRead0.percentageOfHetBase2 * 100 << "%)" << endl;
                
                debugOut << "ReadIds that support the first variant: " << endl;
                if (positionStatsInRead0.hetBase1 == 0) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase1OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase1 == 1) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase1OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase1 == 2) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase1OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase1 == 3) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase1OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase1 == 4) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase1OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                }
                

                debugOut << "ReadIds that support the second variant: " << endl;
                if (positionStatsInRead0.hetBase2 == 0) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase2OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase2 == 1) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase2OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase2 == 2) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase2OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase2 == 3) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase2OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                } else if (positionStatsInRead0.hetBase2 == 4) {
                    for (const auto& orientedReadId : positionStatsInRead0.hetBase2OrientedReadIds) {
                        debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                    }
                }
            }


            

            // XXX
            // --- Dynamic Programming for Grouping Compatible Sites and Phasing ---
            //
            // debugOut << timestamp << "Starting DP for grouping compatible sites for read " << orientedReadId0 << endl;

            // 0. Prepare sorted list of sites
            std::vector<std::pair<uint64_t, AlignmentPositionBaseStats>> sortedSites;
            for (const auto& pair : potentialHetSitesOnOrientedReadId0) {
                sortedSites.push_back(pair);
            }
            // Ensure sites are sorted by position (map iteration order is usually sorted, but explicit sort is safer)
            std::sort(sortedSites.begin(), sortedSites.end(),
                    [](const auto& a, const auto& b) { return a.first < b.first; });

            uint64_t N = sortedSites.size();
            if (N == 0) {
                // debugOut << "No potential het sites found after filtering for read " << orientedReadId0 << ", skipping DP and subsequent phasing." << endl;
                // debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
                continue; // Skip to the next readId if no sites
            } else {
                // debugOut << "Found " << N << " potential heterozygous sites for DP in readId " << readId0 << endl;
                // debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
            }

            // Define compatibility check function locally (or move to class scope)
            auto areCompatibleSites =
                [&](const std::pair<uint64_t, AlignmentPositionBaseStats>& site_pair_i,
                    const std::pair<uint64_t, AlignmentPositionBaseStats>& site_pair_j,
                    const OrientedReadId& targetReadId) -> bool
            {
                const AlignmentPositionBaseStats& stats_i = site_pair_i.second;
                const AlignmentPositionBaseStats& stats_j = site_pair_j.second;

                // --- Check if target read covers both sites informatively ---
                // Target read must have one of the two identified heterozygous bases at each site
                // for a meaningful phase comparison relative to the target.
                bool target_covers_i = (stats_i.baseOfReadId0 == stats_i.hetBase1 || stats_i.baseOfReadId0 == stats_i.hetBase2);
                bool target_covers_j = (stats_j.baseOfReadId0 == stats_j.hetBase1 || stats_j.baseOfReadId0 == stats_j.hetBase2);

                if (!target_covers_i || !target_covers_j) {
                    // Target read doesn't provide a consistent reference phase across both sites.
                    // cout << "Target read " << targetReadId << " does not cover sites " << site_pair_i.first << " and " << site_pair_j.first << " informatively for relative phasing." << endl;
                    return false;
                }
                // --- Target read covers both sites informatively ---

                // 1a. Get all reads overlapping site i
                std::set<OrientedReadId> reads_at_i = stats_i.hetBase1OrientedReadIds;
                reads_at_i.insert(stats_i.hetBase2OrientedReadIds.begin(), stats_i.hetBase2OrientedReadIds.end());
                // Add target read since we know it covers informatively
                // reads_at_i.insert(targetReadId);

                // 1b. Get all reads overlapping site j
                std::set<OrientedReadId> reads_at_j = stats_j.hetBase1OrientedReadIds;
                reads_at_j.insert(stats_j.hetBase2OrientedReadIds.begin(), stats_j.hetBase2OrientedReadIds.end());
                // Add target read since we know it covers informatively
                // reads_at_j.insert(targetReadId);

                // 2. Find common overlapping reads
                std::vector<OrientedReadId> commonOverlappingReads;
                std::set_intersection(reads_at_i.begin(), reads_at_i.end(),
                                    reads_at_j.begin(), reads_at_j.end(),
                                    std::back_inserter(commonOverlappingReads));

                bool found_phase_00 = false; // Flag for reads matching target at both sites (Phase 0-0)
                bool found_phase_11 = false; // Flag for reads differing from target at both sites (Phase 1-1)
                uint64_t found_phase_00_count = 0; // Count of reads matching target at both sites
                uint64_t found_phase_11_count = 0; // Count of reads differing from target at both sites
                // 3. Check phase consistency for each common read relative to target 
                //    and also evidence for both phase patterns
                for (const auto& read : commonOverlappingReads) {
                    // Determine phase at site i (0: matches target, 1: differs from target, -1: unknown/not informative)
                    int phase_i = -1;
                    uint64_t allele_at_i = 100; // Invalid allele
                    if (stats_i.hetBase1OrientedReadIds.count(read)) {
                        allele_at_i = stats_i.hetBase1;
                    } else if (stats_i.hetBase2OrientedReadIds.count(read)) {
                        allele_at_i = stats_i.hetBase2;
                    }

                    if (allele_at_i <= 4) { // Is the allele valid (A,C,G,T,-)?
                        if (allele_at_i == stats_i.baseOfReadId0) {
                            phase_i = 0; // Matches the base of the target read
                        } else {
                            phase_i = 1; // Differs from the base of the target read
                        }
                    } // else phase_i remains -1 (read doesn't support either het allele at this site)

                    // debugOut << "Read " << read << " at site " << site_pair_i.first << " has phase " << phase_i << endl;

                    // Determine phase at site j (0: matches target, 1: differs from target, -1: unknown/not informative)
                    int phase_j = -1;
                    uint64_t allele_at_j = 100;
                    if (stats_j.hetBase1OrientedReadIds.count(read)) {
                        allele_at_j = stats_j.hetBase1;
                    } else if (stats_j.hetBase2OrientedReadIds.count(read)) {
                        allele_at_j = stats_j.hetBase2;
                    }

                    if (allele_at_j <= 4) {
                        if (allele_at_j == stats_j.baseOfReadId0) {
                            phase_j = 0; // Matches target
                        } else {
                            phase_j = 1; // Differs from target
                        } 
                    } // else phase_j remains -1

                    // debugOut << "Read " << read << " at site " << site_pair_j.first << " has phase " << phase_j << endl;

                    // Check for inconsistency or update flags if consistent and informative
                    if (phase_i != -1 && phase_j != -1) { // Both phases defined relative to target for this read
                        if (phase_i != phase_j) {
                            // Found a read with inconsistent phasing relative to the target across the two sites.
                            // debugOut << "Incompatibility: Read " << read << " phase " << phase_i << " at site " << site_pair_i.first << ", phase " << phase_j << " at site " << site_pair_j.first << " relative to target." << endl;
                            return false; // Inconsistent phasing detected
                        } else if (phase_i == 0) { // phase_i == 0 && phase_j == 0
                            found_phase_00 = true; // Found evidence for 0-0 linkage relative to target
                            found_phase_00_count++;
                        } else { // phase_i == 1 && phase_j == 1
                            found_phase_11 = true; // Found evidence for 1-1 linkage relative to target
                            found_phase_11_count++;
                        }
                    } else if (phase_i == -1 || phase_j == -1) {
                        // debugOut << "Read " << read << " is not informative for phasing relative to target at site " << site_pair_i.first << " or site " << site_pair_j.first << "." << endl;
                        return false; // Inconsistent phasing detected
                    }

                }

                // Sites are compatible *relative to the target* only if:
                // 1. No common reads showed inconsistent phasing between positions i and j (phase_i != phase_j).
                // 2. There was at least one read supporting the 0-0 linkage relative to target.
                // 3. There was at least one read supporting the 1-1 linkage relative to target.
                // cout << "Compatibility check for sites " << site_pair_i.first << " and " << site_pair_j.first << " relative to target " << targetReadId << ": found_phase_00=" << found_phase_00 << ", found_phase_11=" << found_phase_11 << endl;
                // return found_phase_00 && found_phase_11;
                return found_phase_00 && found_phase_11 && (found_phase_00_count >= 3) && (found_phase_11_count >= 3);

            };

            // 1. Initialize DP table and parent pointers
            std::vector<uint64_t> LCG(N, 1); // LCG[i] = size of largest compatible group ending at i
            std::vector<int64_t> parent(N, -1); // parent[i] = index j that gave the max LCG[i]

            // 2. Fill DP table using recurrence relation: LCG(i) = max_{j<i, S[j]<->S[i]} {LCG(j)} + 1
            for (uint64_t i = 1; i < N; ++i) {
                for (uint64_t j = 0; j < i; ++j) {
                    // Check compatibility S[j] <-> S[i]
                    if (areCompatibleSites(sortedSites[i], sortedSites[j], orientedReadId0)) {
                        if (LCG[j] + 1 > LCG[i]) {
                            LCG[i] = LCG[j] + 1;
                            parent[i] = j;
                        }
                    }
                }
            }

            // --- Traceback and Phasing ---
            std::vector<bool> isAssigned(N, false);

            // Create pairs of (LCG value, index) to sort by LCG value descending
            std::vector<std::pair<uint64_t, size_t>> sortedLCGIndices;
            std::vector<std::pair<uint64_t, size_t>> isolatedSitesIndices;
            for(size_t i = 0; i < N; ++i) {
                if (LCG[i] > 1) { // Groupped sites (they have LCG > 1)
                    sortedLCGIndices.push_back({LCG[i], i});
                } else if (LCG[i] == 1) { // Isolated sites
                    isolatedSitesIndices.push_back({LCG[i], i});
                }
            }
            // Sort descending by LCG value
            std::sort(sortedLCGIndices.rbegin(), sortedLCGIndices.rend()); 

            // debugOut << "DP results for read " << orientedReadId0 << ": LCG values > 1:" << endl;
            for(const auto& p : sortedLCGIndices) {
                // debugOut << "  Site Index: " << p.second << " (Pos: " << sortedSites[p.second].first << "), LCG: " << p.first << endl;
            }

            std::set<OrientedReadId> tempInPhaseOrientedReads;
            std::set<OrientedReadId> excludedOutOfPhaseOrientedReads;
            std::set<OrientedReadId> involvedOrientedReadsInInformativeSites;
            std::set<OrientedReadId> finalInPhaseOrientedReads;

            for (const auto& lcgPair : sortedLCGIndices) {
                size_t current_i = lcgPair.second;
                uint64_t current_LCG = lcgPair.first;

                if (!isAssigned[current_i]) {
                    std::vector<uint64_t> newClusterPositions; // Keep this for informativeSites logic later if needed
                    std::vector<size_t> newClusterIndices;    // Store indices for phasing logic
                    int traceIndex = current_i;
                    // debugOut << "Starting traceback for cluster ending at index " << current_i << " (Pos: " << sortedSites[current_i].first << ") with LCG " << LCG[current_i] << endl;
    
                    while (traceIndex != -1 && !isAssigned[traceIndex]) {
    
                        const AlignmentPositionBaseStats& currentIndexStats = sortedSites[traceIndex].second;
                        uint64_t targetAllele = currentIndexStats.baseOfReadId0;
                        uint64_t otherAllele = (targetAllele == currentIndexStats.hetBase1) ? currentIndexStats.hetBase2 : currentIndexStats.hetBase1;
    
                        // 1. Populate involved reads for this site
                        // Add all reads supporting either heterozygous allele at this site
                        involvedOrientedReadsInInformativeSites.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
                        involvedOrientedReadsInInformativeSites.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
    
                        // 2. Get in-phase and out-of-phase reads *at this specific site index* relative to the target read
                        // Populate in-phase set (reads with the same allele as the target read at this site)
                        if (targetAllele == currentIndexStats.hetBase1) {
                            tempInPhaseOrientedReads.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
                        } else { // targetAllele must be stats.hetBase2
                            tempInPhaseOrientedReads.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
                        }
                        tempInPhaseOrientedReads.insert(orientedReadId0); // The target read is always in phase with itself
    
                        // Populate out-of-phase set (reads with the other heterozygous allele)
                        if (otherAllele == currentIndexStats.hetBase1) {
                            excludedOutOfPhaseOrientedReads.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
                        } else { // otherAllele must be stats.hetBase2
                            excludedOutOfPhaseOrientedReads.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
                        }
    
                        newClusterPositions.push_back(sortedSites[traceIndex].first); // Store site position
                        newClusterIndices.push_back(traceIndex); // Store site index
                        isAssigned[traceIndex] = true;
                        // debugOut << "  Adding site index " << traceIndex << " (Pos: " << sortedSites[traceIndex].first << ") to compatible group." << endl;
                        traceIndex = parent[traceIndex];
    
                    }
    
                    std::reverse(newClusterPositions.begin(), newClusterPositions.end()); // Store cluster in positional order
                    std::reverse(newClusterIndices.begin(), newClusterIndices.end());   // Store indices in positional order
                    //siteClusters.push_back(newClusterPositions); // Keep the original structure if needed elsewhere
                    //siteIndexClusters.push_back(newClusterIndices); // Use this for phasing
                    // debugOut << "  Finished constructing of a compatible group with " << newClusterPositions.size() << " sites." << endl;
                }
            }


            for (const auto& lcgPair : isolatedSitesIndices) {
                size_t current_i = lcgPair.second;
                uint64_t current_LCG = lcgPair.first;

                // An isolated site is considered informative only if it is supported by a
                // sufficient high number of reads.
                if ((current_LCG == 1) && (!isAssigned[current_i])) {
                    const AlignmentPositionBaseStats& currentIndexStats = sortedSites[current_i].second;
                    uint64_t targetAllele = currentIndexStats.baseOfReadId0;
                    uint64_t otherAllele = (targetAllele == currentIndexStats.hetBase1) ? currentIndexStats.hetBase2 : currentIndexStats.hetBase1;

                    uint64_t hetBase1Coverage = currentIndexStats.hetBase1OrientedReadIds.size();
                    uint64_t hetBase2Coverage = currentIndexStats.hetBase2OrientedReadIds.size();

                    double threshold1 = 0.35;
                    double threshold2 = 0.24;

                    uint64_t minCoverageAmongHetBases = std::min(hetBase1Coverage, hetBase2Coverage);
                    double available = minCoverageAmongHetBases;
                    uint64_t totalCoverageAmongHetBases = hetBase1Coverage + hetBase2Coverage;
                    available = available/((double)(totalCoverageAmongHetBases));

                    if(minCoverageAmongHetBases >= 5) {
                        if(available < threshold2 || totalCoverageAmongHetBases < 10) {
                            continue;
                        }
                    } else if (available < threshold1 || hetBase1Coverage < 4 || totalCoverageAmongHetBases < 10) {
                        continue;
                    }
                    
                    // 1. Populate involved reads for this site
                    // Add all reads supporting either heterozygous allele at this site
                    involvedOrientedReadsInInformativeSites.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
                    involvedOrientedReadsInInformativeSites.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());

                    // 2. Get in-phase and out-of-phase reads *at this specific site index* relative to the target read
                    // Populate in-phase set (reads with the same allele as the target read at this site)
                    if (targetAllele == currentIndexStats.hetBase1) {
                        tempInPhaseOrientedReads.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
                    } else { // targetAllele must be stats.hetBase2
                        tempInPhaseOrientedReads.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
                    }
                    tempInPhaseOrientedReads.insert(orientedReadId0); // The target read is always in phase with itself

                    // Populate out-of-phase set (reads with the other heterozygous allele)
                    if (otherAllele == currentIndexStats.hetBase1) {
                        excludedOutOfPhaseOrientedReads.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
                    } else { // otherAllele must be stats.hetBase2
                        excludedOutOfPhaseOrientedReads.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
                    }

                    isAssigned[current_i] = true;
                    // debugOut << "  Adding site index " << current_i << " (Pos: " << sortedSites[current_i].first << ") with LCG " << current_LCG << "." << endl;
                }
            }



            
            // Find which reads exist only in the tempInPhaseOrientedReads set and not in the excludedOutOfPhaseOrientedReads set and add them to finalInPhaseOrientedReads
            for (const auto& read : tempInPhaseOrientedReads) {
                if (excludedOutOfPhaseOrientedReads.find(read) == excludedOutOfPhaseOrientedReads.end()) {
                    finalInPhaseOrientedReads.insert(read);
                }
            }
            debugOut << "Final in-phase reads for read " << orientedReadId0 << ": " << endl;
            for (const auto& read : finalInPhaseOrientedReads) {
                debugOut << "  ReadId: " << read.getReadId() << " Strand: " << read.getStrand() << endl;
            }

            // XXX
            // --- END OF: Dynamic Programming for Grouping Compatible Sites and Phasing ---
            //





            // XXX
            // --- Keep the correct alignments and forbid the bad alignments based on the compatible site phasing ---
            //

            uint64_t numberOfFirstPassHetAlignments = 0;
            // --- Update thread-local boolean vectors ---
            // Loop over all alignments and mark those involving only reads within finalInPhaseOrientedReads
            // Also forbid alignments involving reads in excludedOutOfPhaseOrientedReads

            // Loop over alignments involving the target read (readId0)
            for (const auto alignmentId : alignmentTable0) {

                const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
                ReadId currentReadId0 = thiAlignmentData.readIds[0];
                ReadId currentReadId1 = thiAlignmentData.readIds[1];
                OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thiAlignmentData.info;

                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != readId0) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

                // Flip strands, if necessary
                if (currentOrientedReadId0.getStrand() != strand0) {
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

                // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
                bool involvesfinalInPhaseOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && finalInPhaseOrientedReads.count(currentOrientedReadId1));
                bool involvesExcludedOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && excludedOutOfPhaseOrientedReads.count(currentOrientedReadId1));
                
                if (involvesfinalInPhaseOrientedReads) {
                    if (phasingThreadData->isReadIdContained[currentOrientedReadId0.getReadId()]) {
                        continue;
                    }
                    // Mark firstPassHet in thread-local vector
                    // debugOut << "Marking alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " as first pass (intra-phase)" << endl;
                    firstPassHetAlignments[alignmentId] = true;
                    numberOfFirstPassHetAlignments++;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

                if (involvesExcludedOrientedReads) {
                    if (phasingThreadData->isReadIdContained[currentOrientedReadId0.getReadId()]) {
                        continue;
                    }
                    // Mark forbidden in thread-local vector
                    // debugOut << "Forbidding alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
                    forbiddenAlignments[alignmentId] = true;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

            }

            // // Loop over all oriented reads in finalInPhaseOrientedReads 
            // for (const auto& orientedRead : finalInPhaseOrientedReads) {
            //     const auto alignmentTableForOrientedRead = alignmentTable[orientedRead.getValue()];
            //     for (const auto alignmentId : alignmentTableForOrientedRead) {

            //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
            //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
            //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            //         AlignmentInfo alignmentInfo = thiAlignmentData.info;

            //         // Swap oriented reads, if necessary.
            //         if (currentOrientedReadId0.getReadId() != orientedRead.getReadId()) {
            //             swap(currentOrientedReadId0, currentOrientedReadId1);
            //             alignmentInfo.swap();
            //         }
            //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedRead.getReadId());

            //         // Flip strands, if necessary
            //         if (currentOrientedReadId0.getStrand() != orientedRead.getStrand()) {
            //             currentOrientedReadId0.flipStrand();
            //             currentOrientedReadId1.flipStrand();
            //             alignmentInfo.reverseComplement();
            //         }
            //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedRead.getStrand());

            //         // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
            //         bool involvesfinalInPhaseOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && finalInPhaseOrientedReads.count(currentOrientedReadId1));
            //         bool involvesExcludedOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && excludedOutOfPhaseOrientedReads.count(currentOrientedReadId1));
                    
            //         if (involvesfinalInPhaseOrientedReads) {
            //             // Mark firstPassHet in thread-local vector
            //             // debugOut << "Marking alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " as first pass (intra-phase)" << endl;
            //             firstPassHetAlignments[alignmentId] = true;
            //             numberOfFirstPassHetAlignments++;
            //             alignmentsAlreadyConsidered[alignmentId] = true;
            //         }

            //         if (involvesExcludedOrientedReads) {
            //             // Mark forbidden in thread-local vector
            //             // debugOut << "Forbidding alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
            //             forbiddenAlignments[alignmentId] = true;
            //             alignmentsAlreadyConsidered[alignmentId] = true;
            //         }
                    
            //     }
            // }
            






            // // We do not perform the analysis on contained reads because the longer reads will
            // // make the decision for them (where they best belong) given that they will have more informative
            // // sites to consider. However, these contained reads should be assigned to one cluster of reads only.
            // // Otherwise, we observe between-haplotypes alignments that connect to
            // // these contained reads from both sides of a het region.
            // // So, we also need to forbid alignments between contained reads and excludedOutOfPhaseOrientedReads.
            // for (const auto& orientedRead : finalInPhaseOrientedReads) {
            //     if (phasingThreadData->isReadIdContained[orientedRead.getReadId()]) {
            //         const auto alignmentTableForContainedRead = alignmentTable[orientedRead.getValue()];
            //         // Loop over alignments involving contained read.
            //         for (const auto alignmentId : alignmentTableForContainedRead) {
            //
            //             const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            //             ReadId currentReadId0 = thiAlignmentData.readIds[0];
            //             ReadId currentReadId1 = thiAlignmentData.readIds[1];
            //             OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            //             OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            //
            //             // Swap oriented reads, if necessary.
            //             if (currentOrientedReadId0.getReadId() != orientedRead.getReadId()) {
            //                 swap(currentOrientedReadId0, currentOrientedReadId1);
            //             }
            //             SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedRead.getReadId());
            //
            //             // Flip strands, if necessary
            //             if (currentOrientedReadId0.getStrand() != orientedRead.getStrand()) {
            //                 currentOrientedReadId0.flipStrand();
            //                 currentOrientedReadId1.flipStrand();
            //             }
            //             SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedRead.getStrand());
            //
            //             // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
            //             bool involvesExcludedOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && excludedOutOfPhaseOrientedReads.count(currentOrientedReadId1));
            //
            //             if (involvesExcludedOrientedReads) {
            //                 debugOut << "Forbidding alignment " << alignmentId << " involving contained read: " << currentReadId0 << " and read: " << currentReadId1 << " because the other read is in excludedOutOfPhaseOrientedReads." << endl;
            //                 forbiddenAlignments[alignmentId] = true;
            //                 alignmentsAlreadyConsidered[alignmentId] = true;
            //             }
            //
            //         }
            //
            //     }
            // }

            // cout << "Number of first pass het alignments for readId " << readId0 << ": " << numberOfFirstPassHetAlignments << endl;

            // debugOut << timestamp << "Finished DP for grouping compatible sites for read " << orientedReadId0 << endl;

            // XXX
            // --- END OF: Keep the correct alignments and forbid the bad alignments based on the compatible site phasing ---
            //


            
            // XXX
            // --- Populate the sites vector for this read's haplotype block ---
            //

            // Create a new Site object for the identified orientedReadId cluster of the target read
            Site newSite;
            newSite.orientedReads.insert(finalInPhaseOrientedReads.begin(), finalInPhaseOrientedReads.end());
            newSite.targetOrientedReadId = orientedReadId0;
            newSite.excludedOrientedReads.insert(excludedOutOfPhaseOrientedReads.begin(), excludedOutOfPhaseOrientedReads.end());

            if (!newSite.orientedReads.empty()) {
                sites.push_back(newSite);
            }

            // Create a site of the reverse strand
            Site newSiteReverseStrand;
            for (const auto& orientedRead : finalInPhaseOrientedReads) {
                OrientedReadId reverseOrientedRead(orientedRead.getReadId(), orientedRead.getStrand() == 0 ? 1 : 0);
                newSiteReverseStrand.orientedReads.insert(reverseOrientedRead);
            }
            for (const auto& orientedRead : excludedOutOfPhaseOrientedReads) {
                OrientedReadId reverseOrientedRead(orientedRead.getReadId(), orientedRead.getStrand() == 0 ? 1 : 0);
                newSiteReverseStrand.excludedOrientedReads.insert(reverseOrientedRead);
            }
            newSiteReverseStrand.targetOrientedReadId = OrientedReadId(orientedReadId0.getReadId(), orientedReadId0.getStrand() == 0 ? 1 : 0);

            if (!newSiteReverseStrand.orientedReads.empty()) {
                sites.push_back(newSiteReverseStrand);
            }


            // XXX
            // --- END OF: Populate the sites vector for this read's haplotype block ---
            //


        } // End loop over read IDs in batch

    } // End while getNextBatch

}



void Assembler::createReadGraph4withStrandSeparation(
    uint64_t maxAlignmentCount,
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

    // Get stats about the reads
    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2*readCount;

    // Initialize result vectors
    vector<bool> forbiddenAlignments(alignmentCount, false);
    vector<bool> firstPassHetAlignments(alignmentCount, false);
    vector<bool> alignmentsAlreadyConsidered(alignmentCount, false);
    vector<Site> sites;

    // --- Parallel Phasing Section ---
    cout << timestamp << "Starting parallel phasing analysis." << endl;
    uint64_t threadCount = assemblerInfo->threadCount;
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Initialize shared data structure for threads
    phasingThreadData = new PhasingThreadData(this, threadCount);

    phasingThreadData->isReadIdContained.resize(readCount, false);

    // Loop over read IDs to find contained reads
    for (ReadId readId = 0; readId < readCount; readId++) {

        const ReadId readId0 = readId;
        const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
        const OrientedReadId orientedReadId0(readId0, strand0);

        // Loop over alignments involving this read
        const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
        for (const auto alignmentId : alignmentTable0) {

            const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            ReadId currentReadId0 = thiAlignmentData.readIds[0];
            ReadId currentReadId1 = thiAlignmentData.readIds[1];
            OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = thiAlignmentData.info;

            // Swap oriented reads, if necessary.
            if (currentOrientedReadId0.getReadId() != readId0) {
                swap(currentOrientedReadId0, currentOrientedReadId1);
                alignmentInfo.swap();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

            // Flip strands, if necessary
            if (currentOrientedReadId0.getStrand() != strand0) {
                currentOrientedReadId0.flipStrand();
                currentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

            const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();

            const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            // TODO: this is wrong
            // const uint64_t maxTrim = assemblerInfo->actualMaxTrim;
            // assemblerOptions.alignOptions.maxTrim
            const uint64_t maxTrim = 50;

            // Check for containment.
            const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
            const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);

            if (isContained0 && !isContained1) {
                // 0 is unambiguously contained in 1.
                phasingThreadData->isReadIdContained[currentReadId0] = true;
            }
            if(isContained1 && !isContained0) {
                // 1 is unambiguously contained in 0.
                phasingThreadData->isReadIdContained[currentReadId1] = true;
            }
            if(isContained0 && isContained1) {
                // Near complete overlap.
                continue;
            }

        }
    }

    // Print the contained reads
    uint64_t containedReads = 0;
    cout << "Contained reads: " << endl;
    for (ReadId readId = 0; readId < readCount; readId++) {
        if (phasingThreadData->isReadIdContained[readId]) {
            cout << "  ReadId: " << readId << endl;
            containedReads++;
        }
    }
    cout << "Finished printing contained reads." << endl;
    cout << "Total number of contained reads is: " << containedReads << endl;

    // Copy the isReadIdContained from threads
    vector<bool> isReadIdContained;
    isReadIdContained.resize(readCount);
    for (ReadId readId = 0; readId < readCount; readId++) {
        isReadIdContained[readId] = phasingThreadData->isReadIdContained[readId];
    }


    // Initialize batch processing for reads
    const uint64_t requestedBatchSize = 1; // Or some other suitable value
    setupLoadBalancing(readCount, requestedBatchSize);


    runThreads(&Assembler::createReadGraph4PhasingThreadFunction, threadCount);


    cout << timestamp << "Parallel phasing analysis completed." << endl;

    // Aggregate results from threads
    cout << timestamp << "Aggregating results from threads." << endl;
    for(uint64_t threadId=0; threadId<threadCount; ++threadId) {
        for(uint64_t alignmentId=0; alignmentId<alignmentCount; ++alignmentId) {
            if(phasingThreadData->threadForbiddenAlignments[threadId][alignmentId]) {
                forbiddenAlignments[alignmentId] = true;
            }
            if(phasingThreadData->threadAlignmentsAlreadyConsidered[threadId][alignmentId]) {
                alignmentsAlreadyConsidered[alignmentId] = true;
            }
            if(phasingThreadData->threadFirstPassHetAlignments[threadId][alignmentId]) {
                firstPassHetAlignments[alignmentId] = true;
            }
        }
        // Aggregate the sites from each thread
        for(const Site& site : phasingThreadData->threadSites[threadId]) {
            sites.push_back(site);
        }
    }

    // Clean up shared data
    delete phasingThreadData;
    phasingThreadData = nullptr;
    cout << timestamp << "Finished aggregating results." << endl;
    // --- End Parallel Phasing Section ---

    // Count results
    const long forbiddenCount = count(forbiddenAlignments.begin(), forbiddenAlignments.end(), true);
    const long firstPassCount = count(firstPassHetAlignments.begin(), firstPassHetAlignments.end(), true);
    const uint64_t sitesCount = sites.size();
    cout << timestamp << "Total forbidden alignments after phasing: " << forbiddenCount << endl;
    cout << timestamp << "Total first pass het alignments after phasing: " << firstPassCount << endl;
    cout << timestamp << "Total sites after phasing: " << sitesCount << endl;




    vector <bool> isReadIdContained2;
    isReadIdContained2.resize(readCount);
    // Loop over read IDs to find contained reads
    for (ReadId readId = 0; readId < readCount; readId++) {

        const ReadId readId0 = readId;
        const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
        const OrientedReadId orientedReadId0(readId0, strand0);

        // Loop over alignments involving this read
        const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
        for (const auto alignmentId : alignmentTable0) {

            if (!firstPassHetAlignments[alignmentId]) {
                continue;
            }
            
            if (forbiddenAlignments[alignmentId]) {
                continue;
            }

            const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            ReadId currentReadId0 = thiAlignmentData.readIds[0];
            ReadId currentReadId1 = thiAlignmentData.readIds[1];
            OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = thiAlignmentData.info;

            // Swap oriented reads, if necessary.
            if (currentOrientedReadId0.getReadId() != readId0) {
                swap(currentOrientedReadId0, currentOrientedReadId1);
                alignmentInfo.swap();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

            // Flip strands, if necessary
            if (currentOrientedReadId0.getStrand() != strand0) {
                currentOrientedReadId0.flipStrand();
                currentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

            const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();

            const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            // TODO: this is wrong
            // const uint64_t maxTrim = assemblerInfo->actualMaxTrim;
            // assemblerOptions.alignOptions.maxTrim
            const uint64_t maxTrim = 50;

            // Check for containment.
            const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
            const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);

            if (isContained0 && !isContained1) {
                // 0 is unambiguously contained in 1.
                isReadIdContained2[currentReadId0] = true;
            }
            if(isContained1 && !isContained0) {
                // 1 is unambiguously contained in 0.
                isReadIdContained2[currentReadId1] = true;
            }
            if(isContained0 && isContained1) {
                // Near complete overlap.
                continue;
            }

        }
    }

    // Print the contained reads
    uint64_t newContainedReads = 0;
    cout << "Contained reads: " << endl;
    for (ReadId readId = 0; readId < readCount; readId++) {
        if (isReadIdContained2[readId]) {
            cout << "  ReadId: " << readId << endl;
            newContainedReads++;
        }
    }
    cout << "Finished printing contained reads." << endl;
    cout << "New total number of contained reads is: " << newContainedReads << endl;









    // // XXX
    // // --- Keep the correct alignments and forbid the bad alignments based on the compatible site phasing ---
    // //

    // uint64_t numberOfFirstPassHetAlignments = 0;
    // // --- Update thread-local boolean vectors ---
    // // Loop over all alignments and mark those involving only reads within finalInPhaseOrientedReads
    // // Also forbid alignments involving reads in excludedOutOfPhaseOrientedReads

    // // Loop over alignments involving the target read (readId0)
    // for (const auto alignmentId : alignmentTable0) {

    //     const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //     ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //     ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //     OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //     OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //     AlignmentInfo alignmentInfo = thiAlignmentData.info;

    //     // Swap oriented reads, if necessary.
    //     if (currentOrientedReadId0.getReadId() != readId0) {
    //         swap(currentOrientedReadId0, currentOrientedReadId1);
    //         alignmentInfo.swap();
    //     }
    //     SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

    //     // Flip strands, if necessary
    //     if (currentOrientedReadId0.getStrand() != strand0) {
    //         currentOrientedReadId0.flipStrand();
    //         currentOrientedReadId1.flipStrand();
    //         alignmentInfo.reverseComplement();
    //     }
    //     SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

    //     // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
    //     bool involvesfinalInPhaseOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && finalInPhaseOrientedReads.count(currentOrientedReadId1));
    //     bool involvesExcludedOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && excludedOutOfPhaseOrientedReads.count(currentOrientedReadId1));
        
    //     if (involvesfinalInPhaseOrientedReads) {
    //         // if (phasingThreadData->isReadIdContained[currentOrientedReadId0.getReadId()]) {
    //         //     continue;
    //         // }
    //         // Mark firstPassHet in thread-local vector
    //         // debugOut << "Marking alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " as first pass (intra-phase)" << endl;
    //         firstPassHetAlignments[alignmentId] = true;
    //         numberOfFirstPassHetAlignments++;
    //         alignmentsAlreadyConsidered[alignmentId] = true;
    //     }

    //     if (involvesExcludedOrientedReads) {
    //         if (phasingThreadData->isReadIdContained[currentOrientedReadId0.getReadId()]) {
    //             continue;
    //         }
    //         // Mark forbidden in thread-local vector
    //         // debugOut << "Forbidding alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
    //         forbiddenAlignments[alignmentId] = true;
    //         alignmentsAlreadyConsidered[alignmentId] = true;
    //     }

    // }











    // for (ReadId readId = 0; readId < readCount; readId++) {
    //
    //     if (readId != 8156 && readId != 8157 && readId != 16706) {
    //         continue;
    //     }
    //
    //     const ReadId readId0 = readId;
    //     const Strand strand0 = 0;
    //     const OrientedReadId orientedReadId0(readId0, strand0);
    //     const auto alignmentTableForOrientedReadId0 = alignmentTable[orientedReadId0.getValue()];
    //     std::set <OrientedReadId> rightHangingOrientedReadIds;
    //     std::set <OrientedReadId> leftHangingOrientedReadIds;
    //     uint64_t numberOfRightHangingReads = 0;
    //     uint64_t numberOfLeftHangingReads = 0;
    //     for (const auto alignmentId : alignmentTableForOrientedReadId0) {
    //
    //         if (firstPassHetAlignments[alignmentId] != true) {
    //             continue;
    //         }
    //
    //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //
    //         // Swap oriented reads, if necessary.
    //         if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
    //             swap(currentOrientedReadId0, currentOrientedReadId1);
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
    //
    //         // Flip strands, if necessary
    //         if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
    //             currentOrientedReadId0.flipStrand();
    //             currentOrientedReadId1.flipStrand();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());
    //
    //         if (thiAlignmentData.info.rightTrim(1) >= 50) {
    //             // This currentOrientedReadId1 is hanging on the right of currentOrientedReadId0
    //             rightHangingOrientedReadIds.insert(currentOrientedReadId1);
    //             continue;
    //         }
    //
    //         if (thiAlignmentData.info.leftTrim(1) >= 50) {
    //             // This currentOrientedReadId1 is hanging on the right of currentOrientedReadId0
    //             leftHangingOrientedReadIds.insert(currentOrientedReadId1);
    //             continue;
    //         }
    //
    //
    //     }
    //     numberOfRightHangingReads = rightHangingOrientedReadIds.size();
    //     numberOfLeftHangingReads = leftHangingOrientedReadIds.size();
    //
    //
    //     bool foundRightHangingOrientedReadsThatDoNotOverlap = false;
    //
    //     // Loop over right-Hanging-OrientedReadIds
    //     for (const auto rightHangingOrientedRead: rightHangingOrientedReadIds) {
    //         uint64_t numberOfRightHangingReadsOnNeighboringOrientedRead = 0;
    //         const auto alignmentTableForRightHangingOrientedRead = alignmentTable[rightHangingOrientedRead.getValue()];
    //         for (const auto alignmentId : alignmentTableForRightHangingOrientedRead) {
    //
    //             const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //             ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //             ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //             OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //             OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //
    //             // Swap oriented reads, if necessary.
    //             if (currentOrientedReadId0.getReadId() != rightHangingOrientedRead.getReadId()) {
    //                 swap(currentOrientedReadId0, currentOrientedReadId1);
    //             }
    //             SHASTA_ASSERT(currentOrientedReadId0.getReadId() == rightHangingOrientedRead.getReadId());
    //
    //             // Flip strands, if necessary
    //             if (currentOrientedReadId0.getStrand() != rightHangingOrientedRead.getStrand()) {
    //                 currentOrientedReadId0.flipStrand();
    //                 currentOrientedReadId1.flipStrand();
    //             }
    //             SHASTA_ASSERT(currentOrientedReadId0.getStrand() == rightHangingOrientedRead.getStrand());
    //
    //             // check if currentOrientedReadId1 is in rightHangingOrientedReadIds
    //             if (rightHangingOrientedReadIds.count(currentOrientedReadId1)) {
    //                 numberOfRightHangingReadsOnNeighboringOrientedRead++;
    //                 continue;
    //             }
    //         }
    //
    //         cout << "The number of Right Hanging Oriented Reads of " << readId0 << " is:" << rightHangingOrientedReadIds.size() << endl;
    //         cout << "The number of Left Hanging Oriented Reads of " << readId0 << " is:" << leftHangingOrientedReadIds.size() << endl;
    //         cout << "The number Of Right Hanging Oriented Reads of " << rightHangingOrientedRead << " is:" << numberOfRightHangingReadsOnNeighboringOrientedRead << endl;
    //         cout << endl;
    //
    //         if (numberOfRightHangingReadsOnNeighboringOrientedRead < rightHangingOrientedReadIds.size() * 0.50) {
    //             foundRightHangingOrientedReadsThatDoNotOverlap = true;
    //             break;
    //         }
    //
    //         if (foundRightHangingOrientedReadsThatDoNotOverlap == true) {
    //             break;
    //         }
    //
    //     }
    //
    //     if (foundRightHangingOrientedReadsThatDoNotOverlap == true) {
    //         // find all alignmentId involving readId0 and forbid them
    //
    //     }
    //
    //
    //
    // }
    //
    //
    //
    //
    //
    //
    
    
    
    
    // Loop over all sites and check if a read is included in both strands in the site
    vector<bool> sitesThatHaveStrandIssuesForbidden(sites.size(), false);
    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        const Site& site = sites[siteId];
        std::set<ReadId> readIdsInSite;
        for(const OrientedReadId orientedReadId: site.orientedReads) {
            readIdsInSite.insert(orientedReadId.getReadId());
        }
        // Check if both strands of the same read are present in the site
        for(const ReadId readId: readIdsInSite) {
            OrientedReadId orientedReadId0(readId, 0);
            OrientedReadId orientedReadId1(readId, 1);
            if (site.orientedReads.count(orientedReadId0) && site.orientedReads.count(orientedReadId1)) {
                cout << timestamp << "Found both strands of read " << readId << " in site " << siteId << "." << endl;
                // print all reads in the site
                cout << "Reads in site " << siteId << ": ";
                for(const OrientedReadId orientedReadId: site.orientedReads) {
                    cout << orientedReadId.getReadId() << " ";
                }
                cout << endl;
                // Mark the site as having strand issues
                sitesThatHaveStrandIssuesForbidden[siteId] = true;
            }
        }
    }
    
    
    // We need to find, for each orientedReadId, the sites that are associated with it
    // and fill in the orientedReadSites
    vector< vector<uint64_t> > orientedReadSites(readCount * 2); // Indexed via OrientedReadId::getValue()
    
    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        // if (sitesThatHaveStrandIssuesForbidden[siteId]) {
        //     continue; // Skip sites that have strand issues
        // }
        const Site& site = sites[siteId];
        OrientedReadId targetOrientedReadId = site.targetOrientedReadId;
        // for(const OrientedReadId orientedReadId: site.orientedReads) {
        //     orientedReadSites[orientedReadId.getValue()].push_back(siteId);
        // }
        orientedReadSites[targetOrientedReadId.getValue()].push_back(siteId);
    }
    
    // for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
    //     // if (sitesThatHaveStrandIssuesForbidden[siteId]) {
    //     //     continue; // Skip sites that have strand issues
    //     // }
    //     const Site& site = sites[siteId];
    //     for(const OrientedReadId orientedReadId: site.orientedReads) {
    //         orientedReadSites[orientedReadId.getValue()].push_back(siteId);
    //     }
    // }
    
    cout << timestamp << "Finished finding for each orientedReadId, the sites that are associated with it." << endl;



    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        const Site& site = sites[siteId];
        OrientedReadId targetOrientedReadId = site.targetOrientedReadId;
        if (isReadIdContained[targetOrientedReadId.getReadId()] && !isReadIdContained2[targetOrientedReadId.getReadId()]) {
            cout << "Trying to analyze readId: " << targetOrientedReadId.getReadId() << endl;
            const auto alignmentTable0 = alignmentTable[targetOrientedReadId.getValue()];
            for (const auto alignmentId : alignmentTable0) {

                const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
                ReadId currentReadId0 = thiAlignmentData.readIds[0];
                ReadId currentReadId1 = thiAlignmentData.readIds[1];
                OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thiAlignmentData.info;

                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != targetOrientedReadId.getReadId()) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetOrientedReadId.getReadId());

                // Flip strands, if necessary
                if (currentOrientedReadId0.getStrand() != targetOrientedReadId.getStrand()) {
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == targetOrientedReadId.getStrand());

                // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
                bool involvesfinalInPhaseOrientedReads = (site.orientedReads.count(currentOrientedReadId0) && site.orientedReads.count(currentOrientedReadId1));
                bool involvesExcludedOrientedReads = (site.excludedOrientedReads.count(currentOrientedReadId0) && site.excludedOrientedReads.count(currentOrientedReadId1));
                
                if (involvesfinalInPhaseOrientedReads) {
                    cout << "Added firstPassHetAlignments alignmentId: " << alignmentId << endl;
                    // if (phasingThreadData->isReadIdContained[currentOrientedReadId0.getReadId()]) {
                    //     continue;
                    // }
                    // Mark firstPassHet in thread-local vector
                    // debugOut << "Marking alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " as first pass (intra-phase)" << endl;
                    firstPassHetAlignments[alignmentId] = true;
                    // numberOfFirstPassHetAlignments++;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

                if (involvesExcludedOrientedReads) {
                    cout << "Added forbiddenAlignments alignmentId: " << alignmentId << endl;
                    // Mark forbidden in thread-local vector
                    // debugOut << "Forbidding alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
                    forbiddenAlignments[alignmentId] = true;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

            }
        }
        
    }

    




    // // Loop over read IDs to find contained reads
    // for (ReadId readId = 0; readId < readCount; readId++) {

    //     const ReadId readId0 = readId;
    //     const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
    //     const OrientedReadId orientedReadId0(readId0, strand0);

    //     // Loop over alignments involving this read
    //     const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
    //     for (const auto alignmentId : alignmentTable0) {

    //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //         AlignmentInfo alignmentInfo = thiAlignmentData.info;

    //         // Swap oriented reads, if necessary.
    //         if (currentOrientedReadId0.getReadId() != readId0) {
    //             swap(currentOrientedReadId0, currentOrientedReadId1);
    //             alignmentInfo.swap();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

    //         // Flip strands, if necessary
    //         if (currentOrientedReadId0.getStrand() != strand0) {
    //             currentOrientedReadId0.flipStrand();
    //             currentOrientedReadId1.flipStrand();
    //             alignmentInfo.reverseComplement();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);
    //     }
    // }


    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //     const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //     ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //     ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //     OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //     OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //     AlignmentInfo alignmentInfo = thiAlignmentData.info;

    //     // Procede only if both reads are contained.
    //     if (!isReadIdContained[currentReadId0] || !isReadIdContained[currentReadId1]) {
    //         continue;
    //     }

    //     bool foundNonContainedReadThatIsPartOfTheSiteReadSet = false;
    //     for (const uint64_t siteId: orientedReadSites[currentOrientedReadId0.getValue()]) {
    //         const Site &site = sites[siteId];
    //         for (const OrientedReadId &orientedReadId: site.orientedReads) {
    //             if (!isReadIdContained[orientedReadId.getReadId()]) {
    //                 foundNonContainedReadThatIsPartOfTheSiteReadSet = true;
    //                 break;
    //                 // forbiddenAlignments[alignmentId] = true;
    //                 // break;
    //             }
    //         }
    //     }

    //     if (foundNonContainedReadThatIsPartOfTheSiteReadSet) {
    //         for (const uint64_t siteId: orientedReadSites[currentOrientedReadId0.getValue()]) {
    //             const Site &site = sites[siteId];
    //             for (const OrientedReadId &orientedReadId: site.orientedReads) {
    //                 if (isReadIdContained[orientedReadId.getReadId()] && (orientedReadId.getReadId() == currentReadId1)) {
    //                     forbiddenAlignments[alignmentId] = true;
    //                     break;
    //                 }
    //             }
    //         }
    //     }

    // }





























    
    // // print all the reads in each site of the orientedReadId 19364-1
    // OrientedReadId testOrientedReadId(3094, 0);
    // cout << "OrientedReadId: " << testOrientedReadId << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId2(3117, 0);
    // cout << "OrientedReadId: " << testOrientedReadId2 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId2 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId2.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId3(3119, 0);
    // cout << "OrientedReadId: " << testOrientedReadId3 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId3 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId3.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId4(2949, 1);
    // cout << "OrientedReadId: " << testOrientedReadId4 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId4 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId4.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId5(2957, 0);
    // cout << "OrientedReadId: " << testOrientedReadId5 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId5 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId5.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId6(2935, 0);
    // cout << "OrientedReadId: " << testOrientedReadId6 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId6 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId6.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId7(1180, 0);
    // cout << "OrientedReadId: " << testOrientedReadId7 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId7 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId7.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId8(1185, 0);
    // cout << "OrientedReadId: " << testOrientedReadId8 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId8 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId8.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId9(1191, 1);
    // cout << "OrientedReadId: " << testOrientedReadId9 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId9 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId9.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId10(1194, 1);
    // cout << "OrientedReadId: " << testOrientedReadId10 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId10 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId10.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId11(1200, 0);
    // cout << "OrientedReadId: " << testOrientedReadId11 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId11 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId11.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId12(3118, 0);
    // cout << "OrientedReadId: " << testOrientedReadId12 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId12 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId12.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;



    // // print all the reads in each site of the orientedReadId 19364-1
    // OrientedReadId testOrientedReadId(3094, 0);
    // cout << "OrientedReadId: " << testOrientedReadId << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId2(3117, 0);
    // cout << "OrientedReadId: " << testOrientedReadId2 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId2 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId2.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId3(3119, 0);
    // cout << "OrientedReadId: " << testOrientedReadId3 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId3 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId3.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId4(2949, 1);
    // cout << "OrientedReadId: " << testOrientedReadId4 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId4 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId4.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId5(1055, 1);
    // cout << "OrientedReadId: " << testOrientedReadId5 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId5 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId5.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId6(1045, 0);
    // cout << "OrientedReadId: " << testOrientedReadId6 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId6 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId6.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId7(1057, 0);
    // cout << "OrientedReadId: " << testOrientedReadId7 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId7 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId7.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId8(1043, 1);
    // cout << "OrientedReadId: " << testOrientedReadId8 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId8 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId8.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId9(1498, 1);
    // cout << "OrientedReadId: " << testOrientedReadId9 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId9 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId9.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId10(1494, 1);
    // cout << "OrientedReadId: " << testOrientedReadId10 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId10 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId10.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId11(1489, 1);
    // cout << "OrientedReadId: " << testOrientedReadId11 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId11 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId11.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId12(1488, 0);
    // cout << "OrientedReadId: " << testOrientedReadId12 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId12 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId12.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;



    //*
    //
    // Create the dynamically adjustable boost readGraph using 
    // all the alignments that are not forbidden.
    //
    //*
    
    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph4AllAlignments readGraphAllAlignments(orientedReadCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // if (!firstPassHetAlignments[alignmentId]) {
        //     continue;
        // }

        if (forbiddenAlignments[alignmentId]) {
            continue;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(thisAlignmentData.readIds[0], 0);
        OrientedReadId orientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);
        AlignmentInfo alignmentInfo = thisAlignmentData.info;

        // Swap them if necessary, depending on the average alignment offset at center.
        if(alignmentInfo.offsetAtCenter() < 0.) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
    }
    
    cout << "The read graph with all alignments that are not forbidden has " << num_vertices(readGraphAllAlignments) << " vertices and " << num_edges(readGraphAllAlignments) << " edges." << endl;




    

    // for (ReadId readId = 0; readId < readCount; readId++) {

    //     if (readId != 3094) {
    //         continue;
    //     }

    //     std::set <OrientedReadId> orientedReadsThatShouldBeForbidden;

    //     const ReadId readId0 = readId;
    //     const Strand strand0 = 0;
    //     const OrientedReadId orientedReadId0(readId0, strand0);
    //     const auto alignmentTableForOrientedReadId0 = alignmentTable[orientedReadId0.getValue()];
        
    //     std::set <OrientedReadId> rightHangingOrientedReadIds;
    //     std::set <OrientedReadId> leftHangingOrientedReadIds;
    //     uint64_t numberOfRightHangingReads = 0;
    //     uint64_t numberOfLeftHangingReads = 0;


    //     for (const auto alignmentId : alignmentTableForOrientedReadId0) {

    //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //         AlignmentInfo alignmentInfo = thiAlignmentData.info;
    
    //         // Swap oriented reads, if necessary.
    //         if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
    //             swap(currentOrientedReadId0, currentOrientedReadId1);
    //             alignmentInfo.swap();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
    
    //         // Flip strands, if necessary
    //         if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
    //             currentOrientedReadId0.flipStrand();
    //             currentOrientedReadId1.flipStrand();
    //             alignmentInfo.reverseComplement();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());

    //         if (firstPassHetAlignments[alignmentId] && !forbiddenAlignments[alignmentId]) {
    //             // cout << "First pass het alignment between " << orientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " is not forbidden." << endl;
    //             continue;
    //         }
    
    //         if (alignmentInfo.rightTrim(1) >= 50 && alignmentInfo.offsetAtCenter() > 0) {
    //             orientedReadsThatShouldBeForbidden.insert(currentOrientedReadId1);
    //         }

    //     }


    //     for (const auto alignmentId : alignmentTableForOrientedReadId0) {

    //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //         AlignmentInfo alignmentInfo = thiAlignmentData.info;
    
    //         // Swap oriented reads, if necessary.
    //         if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
    //             swap(currentOrientedReadId0, currentOrientedReadId1);
    //             alignmentInfo.swap();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
    
    //         // Flip strands, if necessary
    //         if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
    //             currentOrientedReadId0.flipStrand();
    //             currentOrientedReadId1.flipStrand();
    //             alignmentInfo.reverseComplement();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());

    //         if (!firstPassHetAlignments[alignmentId] || forbiddenAlignments[alignmentId]) {
    //             continue;
    //         }
    
    //         if (alignmentInfo.rightTrim(1) >= 50 && alignmentInfo.offsetAtCenter() > 0) {
    //             for (const uint64_t siteId: orientedReadSites[currentOrientedReadId1.getValue()]) {
    //                 const Site &site = sites[siteId];
    //                 for (const OrientedReadId &orientedReadId: site.orientedReads) {
    //                     // check if the orientedReadId is in the orientedReadsThatShouldBeForbidden
    //                     if (orientedReadsThatShouldBeForbidden.find(orientedReadId) != orientedReadsThatShouldBeForbidden.end()) {
    //                         forbiddenAlignments[alignmentId] = true;
    //                         break;
    //                     }
    //                 }
    //             }
    //         }

    //     }






        // for (const auto alignmentId : alignmentTableForOrientedReadId0) {

        //     const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
        //     ReadId currentReadId0 = thiAlignmentData.readIds[0];
        //     ReadId currentReadId1 = thiAlignmentData.readIds[1];
        //     OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
        //     OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
        //     AlignmentInfo alignmentInfo = thiAlignmentData.info;
    
        //     // Swap oriented reads, if necessary.
        //     if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
        //         swap(currentOrientedReadId0, currentOrientedReadId1);
        //         alignmentInfo.swap();
        //     }
        //     SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
    
        //     // Flip strands, if necessary
        //     if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
        //         currentOrientedReadId0.flipStrand();
        //         currentOrientedReadId1.flipStrand();
        //         alignmentInfo.reverseComplement();
        //     }
        //     SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());

        //     if (firstPassHetAlignments[alignmentId] && !forbiddenAlignments[alignmentId]) {
        //         // cout << "First pass het alignment between " << orientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " is not forbidden." << endl;
        //         continue;
        //     }
    
        //     if (alignmentInfo.rightTrim(1) >= 50 && alignmentInfo.offsetAtCenter() > 0) {
        //         for (const uint64_t siteId: orientedReadSites[currentOrientedReadId1.getValue()]) {
        //             const Site &site = sites[siteId];
        //             for (const OrientedReadId &readId: site.orientedReads) {
        //                 orientedReadsThatShouldBeForbidden.insert(readId);
        //             }
        //         }
        //     }
        // }

        // // // print the orientedReadsThatShouldBeForbidden
        // // cout << "Oriented reads that should be forbidden:" << endl;
        // // for (const OrientedReadId &readId: orientedReadsThatShouldBeForbidden) {
        // //     cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
        // // }
        // // cout << endl;

        // if (orientedReadsThatShouldBeForbidden.size() > 0) {
            
        //     for (const auto alignmentId : alignmentTableForOrientedReadId0) {

        //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
        //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
        //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
        //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
        //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
        //         AlignmentInfo alignmentInfo = thiAlignmentData.info;
        
        //         // Swap oriented reads, if necessary.
        //         if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
        //             swap(currentOrientedReadId0, currentOrientedReadId1);
        //             alignmentInfo.swap();
        //         }
        //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
        
        //         // Flip strands, if necessary
        //         if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
        //             currentOrientedReadId0.flipStrand();
        //             currentOrientedReadId1.flipStrand();
        //             alignmentInfo.reverseComplement();
        //         }
        //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());
        
        //         if (alignmentInfo.rightTrim(1) >= 50 && alignmentInfo.offsetAtCenter() > 0) {
        //             cout << "Forbidden alignment between " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " is not forbidden." << endl;
        //             if (orientedReadsThatShouldBeForbidden.find(currentOrientedReadId1) != orientedReadsThatShouldBeForbidden.end()) {
        //                 cout << "Forbidden alignment between " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " is forbidden." << endl;
        //                 forbiddenAlignments[alignmentId] = true;
        //             }
        //         }
        //     }
            
        // }



        // for (const auto alignmentId : alignmentTableForOrientedReadId0) {

        //     const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
        //     ReadId currentReadId0 = thiAlignmentData.readIds[0];
        //     ReadId currentReadId1 = thiAlignmentData.readIds[1];
        //     OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
        //     OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
        //     AlignmentInfo alignmentInfo = thiAlignmentData.info;
    
        //     // Swap oriented reads, if necessary.
        //     if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
        //         swap(currentOrientedReadId0, currentOrientedReadId1);
        //         alignmentInfo.swap();
        //     }
        //     SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
    
        //     // Flip strands, if necessary
        //     if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
        //         currentOrientedReadId0.flipStrand();
        //         currentOrientedReadId1.flipStrand();
        //         alignmentInfo.reverseComplement();
        //     }
        //     SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());

        //     if (!firstPassHetAlignments[alignmentId] || forbiddenAlignments[alignmentId]) {
        //         continue;
        //     }
    
        //     if (alignmentInfo.rightTrim(1) >= 50 && alignmentInfo.offsetAtCenter() > 0) {
        //         for (const uint64_t siteId: orientedReadSites[currentOrientedReadId1.getValue()]) {
        //             const Site &site = sites[siteId];
        //             for (const OrientedReadId &orientedReadId: site.orientedReads) {
        //                 // check if the orientedReadId is in the orientedReadsThatShouldBeForbidden
        //                 if (orientedReadsThatShouldBeForbidden.find(orientedReadId) != orientedReadsThatShouldBeForbidden.end()) {
        //                     forbiddenAlignments[alignmentId] = true;
        //                     break;
        //                 }
        //             }
        //         }
        //     }
        // }



    
    // }




    // for (ReadId readId = 0; readId < readCount; readId++) {

    //     if (isReadIdContained[readId]) {
    //         continue;
    //     }

    //     const ReadId readId0 = readId;
    //     const Strand strand0 = 0;
    //     const OrientedReadId orientedReadId0(readId0, strand0);
    //     const auto alignmentTableForOrientedReadId0 = alignmentTable[orientedReadId0.getValue()];
    //     std::set <OrientedReadId> rightHangingOrientedReadIds;
    //     std::set <OrientedReadId> leftHangingOrientedReadIds;
    //     uint64_t numberOfRightHangingReads = 0;
    //     uint64_t numberOfLeftHangingReads = 0;

    //     bool foundRightHangingOrientedReadsThatAreNotContained = false;
    //     for (const auto alignmentId : alignmentTableForOrientedReadId0) {

    //         if (!firstPassHetAlignments[alignmentId]) {
    //             continue;
    //         }

    //         if (forbiddenAlignments[alignmentId]) {
    //             continue;
    //         }

    //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //         ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //         ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //         OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //         OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //         AlignmentInfo alignmentInfo = thiAlignmentData.info;
    
    //         // Swap oriented reads, if necessary.
    //         if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
    //             swap(currentOrientedReadId0, currentOrientedReadId1);
    //             alignmentInfo.swap();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
    
    //         // Flip strands, if necessary
    //         if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
    //             currentOrientedReadId0.flipStrand();
    //             currentOrientedReadId1.flipStrand();
    //             alignmentInfo.reverseComplement();
    //         }
    //         SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());
    
    //         if (alignmentInfo.rightTrim(1) >= 50) {
    //             if (!isReadIdContained[currentOrientedReadId1.getReadId()]) {
    //                 foundRightHangingOrientedReadsThatAreNotContained = true;
    //                 break;
    //             }
    //         }
    //     }

    //     if (foundRightHangingOrientedReadsThatAreNotContained) {
            
    //         for (const auto alignmentId : alignmentTableForOrientedReadId0) {

    //             if (!firstPassHetAlignments[alignmentId]) {
    //                 continue;
    //             }

    //             if (forbiddenAlignments[alignmentId]) {
    //                 continue;
    //             }

    //             const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //             ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //             ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //             OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //             OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //             AlignmentInfo alignmentInfo = thiAlignmentData.info;
        
    //             // Swap oriented reads, if necessary.
    //             if (currentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
    //                 swap(currentOrientedReadId0, currentOrientedReadId1);
    //                 alignmentInfo.swap();
    //             }
    //             SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedReadId0.getReadId());
        
    //             // Flip strands, if necessary
    //             if (currentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
    //                 currentOrientedReadId0.flipStrand();
    //                 currentOrientedReadId1.flipStrand();
    //                 alignmentInfo.reverseComplement();
    //             }
    //             SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedReadId0.getStrand());
        
    //             if (alignmentInfo.rightTrim(1) >= 50) {
    //                 if (isReadIdContained[currentOrientedReadId1.getReadId()]) {
    //                     forbiddenAlignments[alignmentId] = true;
    //                 }
    //             }
    //         }
            
    //     }
    
    // }
    
    
    
    
    
    















    
    
    



    // // loop over the indices of the orientedReadSites and check if the orientedReadId0 is present in the sites
    // vector< pair<uint64_t, uint64_t> > pairsOfSites;
    // for (uint64_t readId = 0; readId < readCount; readId++) {
    //
    //     if (readId != 1112) {
    //         continue; // Skip if not the specific readId we are interested in
    //     }
    //
    //     OrientedReadId orientedReadId0(readId, 0);
    //
    //     vector<uint64_t> v = orientedReadSites[orientedReadId0.getValue()];
    //
    //     // Add a check to prevent accessing v.size()-1 when v is empty or has only one element.
    //     if (v.size() == 0) {
    //         cout << timestamp << "Skipping site pair generation for orientedReadId with less than 2 associated sites." << endl;
    //         continue; // Skip if there are not enough elements to form a pair
    //     }
    //
    //     if (v.size() == 1) {
    //         const uint64_t siteId0 = v[0];
    //         pairsOfSites.push_back(make_pair(siteId0, siteId0));
    //         // cout << timestamp << "Skipping site pair generation for orientedReadId with less than 2 associated sites." << endl;
    //         continue; // Skip if there are not enough elements to form a pair
    //     }
    //
    //     for(uint64_t i0=0; i0<v.size()-1; i0++) {
    //         const uint64_t siteId0 = v[i0];
    //         if (sitesThatHaveStrandIssuesForbidden[siteId0]) {
    //             continue; // Skip sites that have strand issues
    //         }
    //         for(uint64_t i1=i0+1; i1<v.size(); i1++) {
    //             const uint64_t siteId1 = v[i1];
    //             if (sitesThatHaveStrandIssuesForbidden[siteId1]) {
    //                 continue; // Skip sites that have strand issues
    //             }
    //             pairsOfSites.push_back(make_pair(siteId0, siteId1));
    //         }
    //     }
    //
    // }
    //
    //
    // // We need to generate all possible pairs of sites from the orientedReadSites.
    // // The siteId of the first site in each pair will always be less than the siteId of the second site in each pair
    // vector< pair<uint64_t, uint64_t> > pairsOfSites;
    // for(const vector<uint64_t>& v: orientedReadSites) {
    //
    //     // Add a check to prevent accessing v.size()-1 when v is empty or has only one element.
    //     if (v.size() == 0) {
    //         cout << timestamp << "Skipping site pair generation for orientedReadId with less than 2 associated sites." << endl;
    //         continue; // Skip if there are not enough elements to form a pair
    //     }
    //
    //     if (v.size() == 1) {
    //         const uint64_t siteId0 = v[0];
    //         pairsOfSites.push_back(make_pair(siteId0, siteId0));
    //         // cout << timestamp << "Skipping site pair generation for orientedReadId with less than 2 associated sites." << endl;
    //         continue; // Skip if there are not enough elements to form a pair
    //     }
    //
    //     for(uint64_t i0=0; i0<v.size()-1; i0++) {
    //         const uint64_t siteId0 = v[i0];
    //         if (sitesThatHaveStrandIssuesForbidden[siteId0]) {
    //             continue; // Skip sites that have strand issues
    //         }
    //         for(uint64_t i1=i0+1; i1<v.size(); i1++) {
    //             const uint64_t siteId1 = v[i1];
    //             if (sitesThatHaveStrandIssuesForbidden[siteId1]) {
    //                 continue; // Skip sites that have strand issues
    //             }
    //             pairsOfSites.push_back(make_pair(siteId0, siteId1));
    //         }
    //     }
    // }
    //
    // cout << timestamp << "Finished generating all possible pairs of sites from the orientedReadSites." << endl;
    //
    //
    // // Because the siteId of the first site in each pair will always be less than the siteId of
    // // the second site in each pair, we will have to deduplicate the pairs.
    // // --- Deduplicate pairs and count common reads ---
    // vector<uint64_t> commonReadsCount; // Stores the count of common reads for each unique pair
    // deduplicateAndCount(pairsOfSites, commonReadsCount);
    // SHASTA_ASSERT(commonReadsCount.size() == pairsOfSites.size());
    //
    // cout << timestamp << "Finished deduplicating pairs." << endl;


    // // Set up the disjoint sets data structure. ---
    // // --- Each vertex is a heterozygous phased site ---
    // const uint64_t vertexCount = sites.size();
    // vector<uint64_t> rankHetSites(vertexCount);
    // vector<uint64_t> parentHetSites(vertexCount);
    // boost::disjoint_sets<uint64_t*, uint64_t*> disjointSetsHetSites(&rankHetSites[0], &parentHetSites[0]);
    // for(uint64_t i=0; i<vertexCount; i++) {
    //     disjointSetsHetSites.make_set(i);
    // }

    // cout << timestamp << "Finished fifth pass of phasing." << endl;

    // // --- Loop over all pairs of sites that share at least m OrientedReadIds. ---
    // vector<bool> sitesThatHaveStrandIssues(sites.size(), false);
    // const uint64_t minCommonReadsForMerging = 2;
    // uint64_t mergeCount = 0;
    // uint64_t pairsConsidered = 0;
    // for(uint64_t i=0; i<pairsOfSites.size(); i++) {
    //     pairsConsidered++;
    //     const pair<uint64_t, uint64_t>& p = pairsOfSites[i];
    //     const uint64_t commonOrientedReadCount = commonReadsCount[i];
    //     const uint64_t siteId0 = p.first;
    //     const uint64_t siteId1 = p.second;
    //     const Site& site0 = sites[siteId0];
    //     const Site& site1 = sites[siteId1];

    //     if (sitesThatHaveStrandIssuesForbidden[siteId1] || sitesThatHaveStrandIssuesForbidden[siteId0]) {
    //         continue; // Skip sites that have strand issues
    //     }
        
    //     // // Decide if these two should be merged based on the number of common reads
    //     // bool merge = (commonOrientedReadCount >= minCommonReadsForMerging);


    //     // find the number of common reads between the two sites without using commonOrientedReadCount
    //     std::set<OrientedReadId> commonReadsSet;
    //     std::set_intersection(
    //         site0.orientedReads.begin(), site0.orientedReads.end(),
    //         site1.orientedReads.begin(), site1.orientedReads.end(),
    //         std::inserter(commonReadsSet, commonReadsSet.begin())
    //     );
    //     const uint64_t calculatedCommonReadCount = commonReadsSet.size();

    //     // Decide if these two should be merged based on the number of common reads
    //     bool merge = (calculatedCommonReadCount >= minCommonReadsForMerging);

    //     // Check if the sites contain the same read in different strands
    //     if (merge) {
    //         for (const OrientedReadId& readId0 : site0.orientedReads) {
    //             for (const OrientedReadId& readId1 : site1.orientedReads) {
    //                 if (readId0.getReadId() == readId1.getReadId() && 
    //                     readId0.getStrand() != readId1.getStrand()) {
    //                     // If the same read is present in both sites but on different strands,
    //                     // we should not merge them.
    //                     cout << timestamp << "Not merging sites " << siteId0 << " and " << siteId1 
    //                          << " because they contain the same read " 
    //                          << readId0.getReadId() << " on different strands." << endl;
    //                     sitesThatHaveStrandIssues[siteId0] = true;
    //                     sitesThatHaveStrandIssues[siteId1] = true;
    //                 }
    //             }
    //         }
    //     }
        
    //     if(merge && !sitesThatHaveStrandIssues[siteId0] && !sitesThatHaveStrandIssues[siteId1]) {
    //         // Merge the sets (connectedComponents) containing siteId0 and siteId1
    //         disjointSetsHetSites.union_set(siteId0, siteId1);
    //         mergeCount++;
    //     }
    // }
    // cout << timestamp << "Finished site merging. Considered " << pairsConsidered << " pairs, merged " << mergeCount << " pairs based on common read count >= " << minCommonReadsForMerging << "." << endl;


    // // --- At this point, disjointSetsHetSites.find_set(siteId) gives the id of the merged set
    // // --- that siteId is part of. This id is in [0, vertexCount).
    // vector< vector<uint64_t> > connectedComponents;
    // connectedComponents.resize(vertexCount);
    // for(uint64_t siteId=0; siteId<vertexCount; siteId++) {
    //     if (sitesThatHaveStrandIssuesForbidden[siteId]) {
    //         continue; // Skip sites that have strand issues
    //     }
    //     const uint64_t componentId = disjointSetsHetSites.find_set(siteId);
    //     connectedComponents[componentId].push_back(siteId);
    // }


    // // loop over the connected components and print the reads participating in this component
    // cout << timestamp << "Connected components after merging:" << endl;
    // for (uint64_t componentId = 0; componentId < connectedComponents.size(); ++componentId) {
    //     const vector<uint64_t>& siteIdsInComponent = connectedComponents[componentId];

    //     // Skip components that became empty after merging (their sites were moved)
    //     if (siteIdsInComponent.empty()) {
    //         continue;
    //     }

    //     cout << "Component " << componentId << " (Size: " << siteIdsInComponent.size() << " sites):" << endl;

    //     // Collect all unique OrientedReadIds from all sites in this component
    //     std::set<OrientedReadId> componentOrientedReads;
    //     for (const uint64_t siteId : siteIdsInComponent) {
            
    //         const Site& site = sites[siteId];
    //         componentOrientedReads.insert(site.orientedReads.begin(), site.orientedReads.end());
    //         cout << "  Site " << siteId << " has " << site.orientedReads.size() << " oriented reads." << endl;
    //         // print the reads in this site
    //         for (const OrientedReadId& orientedReadId : site.orientedReads) {
    //             cout << "    Oriented Read: " << orientedReadId.getReadId() << "-" << orientedReadId.getStrand() << endl;
    //         }
            
    //     }

    //     // Print the unique OrientedReadIds for this component
    //     cout << "  Oriented Reads (" << componentOrientedReads.size() << "): " << endl;
    //     for (const OrientedReadId& orientedReadId : componentOrientedReads) {
    //         cout << orientedReadId << endl;
    //     }
    //     cout << endl;
    // }
    



    //
    //
    // We finished looping over the reads and found the haplotype sets for each one
    // Now we need to create the first pass read graph that involves only those haplotype specific 
    // alignments we found in the het sites
    //
    //


    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);
    // vector<bool> keepAlignment(alignmentCount, true);
    // createReadGraphUsingSelectedAlignments(keepAlignment);
    // return;

    // for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {

    //     AlignmentData& alignment = alignmentData[alignmentId];

    //     if (forbiddenAlignments[alignmentId]) {
    //         continue;
    //     }
        
    //     if (firstPassHetAlignments[alignmentId]) {
    //         // Add the alignment to the read graph.
    //         keepAlignment[alignmentId] = true;
    //         alignment.info.isInReadGraph = 1;
    //     }
    // }

    // createReadGraphUsingSelectedAlignments(keepAlignment);
    // return;

    
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

    // const double WThreshold = 1e-8;
    const double logWThreshold = log(WThreshold);

    // const double WThresholdForBreaks = 1e+15;
    const double logWThresholdForBreaks = log(WThresholdForBreaks);

    vector< pair<uint64_t, double> > alignmentTableHetSites;
    vector< pair<uint64_t, double> > alignmentTableHetSitesPlusNotForbidden;
    
    // Keep track of which readIds were used in alignments
    vector<bool> readUsed(readCount, false);



    //
    //
    //
    //
    // Do a first pass in which we allow only het loci in order of increasing Q, 
    // then do another pass to fill the breaks where you also allow all the other alignments
    //
    //
    //
    //



    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {

        cout << "Processing alignment " << alignmentId << endl;
        
        if (!firstPassHetAlignments[alignmentId]) {
            continue;
        }

        if (forbiddenAlignments[alignmentId]) {
            continue;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // The alignment is stored as an alignment between readId0 on strand 0
        // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
        // The reverse complement alignment also exists, but is not stored explicitly.
        
        ReadId readId0 = thisAlignmentData.readIds[0];
        ReadId readId1 = thisAlignmentData.readIds[1];
        bool isSameStrand = thisAlignmentData.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
        OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.
        AlignmentInfo alignmentInfo = thisAlignmentData.info;

        // // Swap them if necessary, depending on the average alignment offset at center.
        // if(alignmentInfo.offsetAtCenter() < 0.) {
        //     swap(orientedReadId0, orientedReadId1);
        //     alignmentInfo.swap();
        // }

        // if (readId0 == 3118 || readId1 == 3118 || readId0 == 3251 || readId1 == 3251) {
        //     continue;
        // }

        // if (bridgingReads.count(readId0) || bridgingReads.count(readId1)) {
        //     continue;
        // }

        // // Swap them if necessary, depending on the average alignment offset at center.
        // if(thisAlignmentData.info.offsetAtCenter() < 0.) {
        //     swap(orientedReadId0, orientedReadId1);
        // }

        // if (isReadIdContained2[readId0] || isReadIdContained2[readId1]) {
        //     continue;
        // }

        // if (isReadIdContained2[readId1]) {
        //     continue;
        // }

        // Store this pair of edges in our edgeTable.
        const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
        const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
        const double L = double(range0 + range1)/2.;
        // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
        const double errorRateRle = thisAlignmentData.info.errorRateRle;
        const double nRLE = errorRateRle * 2 * L;
        // const double markerCount = thisAlignmentData.info.markerCount;

        // logQ(n) = αn - δL
        const double logQ = alpha * double(nRLE) - delta * L;

        // Add the alignment to the table along with its logQ.
        alignmentTableHetSites.push_back(make_pair(alignmentId, logQ));
        readUsed[readId0] = true;
        readUsed[readId1] = true;

        keepAlignment[alignmentId] = true;
        thisAlignmentData.info.isInReadGraph = 1;
        cout << "Added alignment " << alignmentId << " to the read graph." << endl;

        // // This time use the regular Bayesian filtering
        // if (logQ <= logWThreshold) {
        //     alignmentTableHetSites.push_back(make_pair(alignmentId, logQ));
        //     alignmentsAlreadyConsidered[alignmentId] = true;
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        //     keepAlignment[alignmentId] = true;
        //     thisAlignmentData.info.isInReadGraph = 1;
        // }
        // } else if(logQ <= logWThresholdForBreaks){
        //     alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
        //     alignmentsAlreadyConsidered[alignmentId] = true;
        //     keepAlignmentsForBreaks[alignmentId] = true;
        // }

    }

    sort(alignmentTableHetSites.begin(), alignmentTableHetSites.end(), OrderPairsBySecondOnly<uint64_t, double>());
    cout << "The alignmentTableHetSites has " << alignmentTableHetSites.size() << " entries." << endl;



    createReadGraphUsingSelectedAlignments(keepAlignment);
    return;























    


    // Maintain a vector containing the degree of each vertex
    // verticesDegree[vertexID] -> degree
    vector<uint64_t> verticesDegree(orientedReadCount, 0);

    cout << "Number of reads: " << readCount << endl;
    cout << "Number of oriented reads: " << orientedReadCount << endl;


    ///
    // Process the HET SITES ALIGNMENTS in order of increasing Q. 
    //
    // i.   Start with no edges in the read graph. 
    // ii.  Process alignments spanning the het sites in order of increasing Q. 
    // iii. If the alignment breaks strand separation, it is skipped. 
    // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped. (this step is skipped)  
    // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.

    // Initiallize disjoint sets for HET SITES ALIGNMENTS
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }

    // // Flag all alignments as not to be kept.
    // vector<bool> keepAlignment(alignmentCount, false);

    // Process alignments in order of increasing Q
    vector alignmentTablesToProcess({alignmentTableHetSites});
    
    uint64_t crossStrandEdgeCountHet = 0;

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


            if (isReadIdContained[readId0] || isReadIdContained[readId1]) {
                continue;
            }

            // if (readId0 == 8156 || readId0 == 8156 || readId1 == 8157 || readId1 == 8157) {
            //     continue;
            // }


            // keepAlignment[alignmentId] = true;
            // alignment.info.isInReadGraph = 1;
            // continue;

            // Get the connected components that these oriented reads are in.
            const uint64_t a0 = disjointSets.find_set(A0.getValue());
            const uint64_t b0 = disjointSets.find_set(B0.getValue());
            const uint64_t a1 = disjointSets.find_set(A1.getValue());
            const uint64_t b1 = disjointSets.find_set(B1.getValue());


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
                crossStrandEdgeCountHet += 2;
                continue;
            }

            

            // // If both vertices of the potential edge have at least the required minimum number 
            // // of neighbors, the alignment is also skipped. 
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

    // Print how many alignments were kept in this step
    const long keepCountR1 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Finding strict disjointSets step involving alignments in HET SITES: Keeping " << keepCountR1 << " alignments in HET sites out of " << keepAlignment.size() << " alignments in general."<< endl;
    cout << "Found " << crossStrandEdgeCountHet << " cross strand edges in HET SITES during disjointSets construction." << endl;










    // Process alignments in order of increasing Q
    // vector alignmentTablesToProcess({alignmentTableHetSites});
    
    // uint64_t crossStrandEdgeCountHet = 0;

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


            // if (isReadIdContained[readId0] || isReadIdContained[readId1]) {
            //     continue;
            // }

            // if (readId0 == 8156 || readId0 == 8156 || readId1 == 8157 || readId1 == 8157) {
            //     continue;
            // }


            // keepAlignment[alignmentId] = true;
            // alignment.info.isInReadGraph = 1;
            // continue;

            // Get the connected components that these oriented reads are in.
            const uint64_t a0 = disjointSets.find_set(A0.getValue());
            const uint64_t b0 = disjointSets.find_set(B0.getValue());
            const uint64_t a1 = disjointSets.find_set(A1.getValue());
            const uint64_t b1 = disjointSets.find_set(B1.getValue());


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
                crossStrandEdgeCountHet += 2;
                continue;
            }

            

            // // If both vertices of the potential edge have at least the required minimum number 
            // // of neighbors, the alignment is also skipped. 
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




    





































    // // const double epsilon = 1e-4;
    // // const double delta = 5e-4;
    // const double alpha = log(1 + delta/(2*epsilon));

    // // const double WThreshold = 1e-8;
    // const double logWThreshold = log(WThreshold);

    // // const double WThresholdForBreaks = 1e+15;
    // const double logWThresholdForBreaks = log(WThresholdForBreaks);

    // vector< pair<uint64_t, double> > alignmentTableHetSites;
    
    // // Keep track of which readIds were used in alignments
    // vector<bool> readUsed(readCount, false);

    // // Maintain a vector containing the degree of each vertex
    // // verticesDegree[vertexID] -> degree
    // vector<uint64_t> verticesDegree(orientedReadCount, 0);

    // // Initiallize disjoint sets for HET SITES ALIGNMENTS
    // vector<ReadId> rank(orientedReadCount);
    // vector<ReadId> parent(orientedReadCount);
    // boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    // for(ReadId readId=0; readId<readCount; readId++) {
    //     for(Strand strand=0; strand<2; strand++) {
    //         disjointSets.make_set(OrientedReadId(readId, strand).getValue());
    //     }
    // }

    // // Flag all alignments as not to be kept.
    // vector<bool> keepAlignment(alignmentCount, false);

    // // Process alignments in order of increasing Q
    // vector alignmentTablesToProcess({alignmentTableHetSites});
    
    // uint64_t crossStrandEdgeCountHet = 0;

    // // loop over orientedReadSites
    // for (OrientedReadId::Int orientedReadIdValue = 0; orientedReadIdValue < orientedReadSites.size(); ++orientedReadIdValue) {
    //     const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
    //     const ReadId readId = orientedReadId.getReadId();

    //     if (isReadIdContained[readId]) {
    //         continue; // Skip contained reads
    //     }

    //     const vector<uint64_t>& sitesForThisRead = orientedReadSites[orientedReadIdValue];

    //     // loop over the sites for this orientedReadId
    //     std::set<ReadId> siteReadIds;
    //     for (const uint64_t siteId : sitesForThisRead) {
    //         const Site& site = sites[siteId];
            
    //         // loop over the oriented reads in this site           
    //         for (const OrientedReadId& orientedReadIdInSite : site.orientedReads) {
    //             siteReadIds.insert(orientedReadIdInSite.getReadId());
    //         }
            
    //     }

    //     // Loop over all alignments.
    //     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {

    //         if (forbiddenAlignments[alignmentId]) {
    //             continue;
    //         }

    //         // Get the alignment data
    //         AlignmentData& alignment = alignmentData[alignmentId];
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

    //         if (siteReadIds.count(readId0) == 0 || siteReadIds.count(readId1) == 0) {
    //             // Neither read of the alignment is in the current site.
    //             continue;
    //         }


    //         // Get the connected components that these oriented reads are in.
    //         const uint64_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint64_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint64_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint64_t b1 = disjointSets.find_set(B1.getValue());


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
    //             crossStrandEdgeCountHet += 2;
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

    //     }

    //     // Verify that for any read the two oriented reads are in distinct
    //     // connected components.
    //     for(ReadId readId=0; readId<readCount; readId++) {
    //         const OrientedReadId orientedReadId0(readId, 0);
    //         const OrientedReadId orientedReadId1(readId, 1);
    //         SHASTA_ASSERT(
    //             disjointSets.find_set(orientedReadId0.getValue()) !=
    //             disjointSets.find_set(orientedReadId1.getValue())
    //         );
    //     }


    // }







































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
        if(!keepThisAlignment) {
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

    cout << "The read graph that used the alignments in HET sites has " << num_vertices(readGraph) << " vertices and " << num_edges(readGraph) << " edges." << endl;



    // UNCOMMENT
    vector< pair<uint64_t, double> > alignmentTable;
    vector< pair<uint64_t, double> > alignmentTableNotPassFilter;
    // Flag alignments to be kept for break detection.
    vector<bool> keepAlignmentsForBreaks(alignmentCount, false);

    //
    //
    //
    //
    // Do a second pass in which we allow all the other alignments in order of increasing Q
    // after we filter them out using the Bayesian filtering
    //
    //
    //
    //

    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if((alignmentId % 100000) == 0) {
            cout << timestamp << alignmentId << "/" << alignmentCount << endl;
        }

        if (alignmentsAlreadyConsidered[alignmentId]) {
            continue;
        }

        if (forbiddenAlignments[alignmentId]) {
            continue;
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
        const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
        const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
        const double L = double(range0 + range1)/2.;
        // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
        const double errorRateRle = thisAlignmentData.info.errorRateRle;
        const double nRLE = errorRateRle * 2 * L;
        // const double markerCount = thisAlignmentData.info.markerCount;

        // logQ(n) = αn - δL
        const double logQ = alpha * double(nRLE) - delta * L;

        // This time use the regular Bayesian filtering
        if (logQ <= logWThreshold) {
            alignmentTable.push_back(make_pair(alignmentId, logQ));
            alignmentsAlreadyConsidered[alignmentId] = true;
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logWThresholdForBreaks){
            alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
            alignmentsAlreadyConsidered[alignmentId] = true;
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
    

    // boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    // for(ReadId readId=0; readId<readCount; readId++) {
    //     for(Strand strand=0; strand<2; strand++) {
    //         disjointSets.make_set(OrientedReadId(readId, strand).getValue());
    //     }
    // }

    // // Flag all alignments as not to be kept.
    // vector<bool> keepAlignment(alignmentCount, false);

    // Process alignments in order of increasing Q
    alignmentTablesToProcess = {alignmentTable};

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
            const uint64_t a0 = disjointSets.find_set(A0.getValue());
            const uint64_t b0 = disjointSets.find_set(B0.getValue());
            const uint64_t a1 = disjointSets.find_set(A1.getValue());
            const uint64_t b1 = disjointSets.find_set(B1.getValue());


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



            // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
            //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
            //     continue;
            // }

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

    // Print how many alignments were kept in this step
    const long keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Finding strict disjointSets step: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;






    //*
    //
    // Update the dynamically adjustable boost readGraph using the alignments we selected not involving HET sites.
    //
    //*


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
    // ReadGraph4AllAlignments readGraphAllAlignments(orientedReadCount);

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
    long potentialDeadEndReadCount = count(potentialDeadEndReads.begin(), potentialDeadEndReads.end(), true);
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

            bool result = readGraph.findPathWithPositiveOffset(orientedReadId, paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

            // Check if we found a read path with positive offset.
            // If yes, the function findPathWithPositiveOffset will return 1, if not, it will return 0.
            if(result == 1) {
                // cout << "Found a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset" << endl;
                // //print the paths and then the pathsOffsets
                // for (uint64_t i = 0; i < paths.size(); i++) {
                //     cout << "Path " << i << ": ";
                //     for (uint64_t j = 0; j < paths[i].size(); j++) {
                //         cout << paths[i][j].getReadId() << " ";
                //     }
                //     cout << endl;
                //     cout << "PathOffsets " << i << ": ";
                //     for (uint64_t j = 0; j < pathsOffsets[i].size(); j++) {
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

    std::unordered_map<uint64_t, vector<uint64_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet;
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
                uint64_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

                // Find all EndNodesWithNoIncomingNodes that are in the same disjointSet
                for(uint64_t id=0; id<orientedReadCount; id++) {
                    // check if the id is an endNode with no incoming nodes
                    if(finalDeadEndReadsWithNoIncomingNodes[id]){
                        // get the disjointSet of the id
                        uint64_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
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
        uint64_t value = p.first;
        OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(value));
        SHASTA_ASSERT(orientedReadId.getValue() == value);
        // const ReadId readId = orientedReadId.getReadId();
        // const Strand strand = orientedReadId.getStrand();

        // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in the same disjointSet:" << endl;
        vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
        for(auto& node : p.second) {
            OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(ReadId(node));
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
            // const uint64_t a0 = disjointSets.find_set(A0.getValue());
            // const uint64_t b0 = disjointSets.find_set(B0.getValue());
            // const uint64_t a1 = disjointSets.find_set(A1.getValue());
            // const uint64_t b1 = disjointSets.find_set(B1.getValue());


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
                    const uint64_t a0v2 = disjointSets.find_set(A0v2.getValue());
                    const uint64_t b0v2 = disjointSets.find_set(B0v2.getValue());
                    const uint64_t a1v2 = disjointSets.find_set(A1v2.getValue());
                    const uint64_t b1v2 = disjointSets.find_set(B1v2.getValue());


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
    const long keepCountR3 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Adding alignments for break bridging to connect endNodes in the same disjointSet: Keeping " << keepCountR3 << " alignments of " << keepAlignment.size() << endl;




    
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


    std::unordered_map<uint64_t, vector<uint64_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet;

    // Now, we try to find EndNodesWithNoIncomingNodes that are in different disjointSets than EndNodesWithNoOutgoingNodes.
    for (ReadId readId = 0; readId < readCount; readId++) {
        for (Strand strand = 0; strand < 2; strand++) {
            OrientedReadId orientedReadId(readId, strand);

            // check if the orientedReadId is an endNode with no outgoing nodes
            if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
                // Get it's disjointSet
                uint64_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

                // Find all EndNodesWithNoIncomingNodes that are in different disjointSet
                for(uint64_t id=0; id<orientedReadCount; id++) {
                    // check if the id is an endNode with no incoming nodes
                    if(finalDeadEndReadsWithNoIncomingNodes[id]){
                        // get the disjointSet of the id
                        uint64_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
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
        uint64_t value = p.first;
        OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(value));
        SHASTA_ASSERT(orientedReadId.getValue() == value);
        // const ReadId readId = orientedReadId.getReadId();
        // const Strand strand = orientedReadId.getStrand();

        // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in a different disjointSet: " << endl;
        vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
        for(auto& node : p.second) {
            OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(ReadId(node));
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
            // const uint64_t a0 = disjointSets.find_set(A0.getValue());
            // const uint64_t b0 = disjointSets.find_set(B0.getValue());
            // const uint64_t a1 = disjointSets.find_set(A1.getValue());
            // const uint64_t b1 = disjointSets.find_set(B1.getValue());


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
                    const uint64_t a0v2 = disjointSets.find_set(A0v2.getValue());
                    const uint64_t b0v2 = disjointSets.find_set(B0v2.getValue());
                    const uint64_t a1v2 = disjointSets.find_set(A1v2.getValue());
                    const uint64_t b1v2 = disjointSets.find_set(B1v2.getValue());


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
    const long keepCountR4 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Adding alignments for break bridging to connect endNodes in different disjointSets: Keeping " << keepCountR4 << " alignments of " << keepAlignment.size() << endl;



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




    //*
    //
    // Create the read graph using the alignments we selected.
    //
    //*
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
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<uint64_t, uint64_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[ReadId(p.second)]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    uint64_t accumulatedOrientedReadCount = 0;
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










// void Assembler::createReadGraph4withStrandSeparation(
//     uint64_t maxAlignmentCount,
//     double epsilon,
//     double delta,
//     double WThreshold,
//     double WThresholdForBreaks
//     )
// {
//     cout << timestamp << "createReadGraph4 with strand separation begins" << endl;

//     // Get the total number of stored alignments.
//     const uint64_t alignmentCount = alignmentData.size();
//     SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

//     // Get stats about the reads
//     const uint64_t readCount = reads->readCount();
//     const uint64_t orientedReadCount = 2*readCount;

    


//     vector<bool> forbiddenAlignments(alignmentCount, false);
//     vector<bool> forbiddenReads(readCount, false);
//     vector<uint64_t> firstPassHetAlignments(alignmentCount, false);
//     std::set<ReadId> bridgingReads;
//     vector<uint64_t> alignmentsAlreadyConsidered(alignmentCount, false);
//     std::set<ReadId> readIdsPseudoHap1;
//     std::set<ReadId> readIdsPseudoHap2;
//     vector<Site> sites;



//     // Loop over reads.
//     for(ReadId readId=0; readId<readCount; readId++) {

//         if (readId != 0) {
//             continue;
//         }

//         cout << "Working on read " << readId << endl;

//         const ReadId readId0 = readId;
//         const Strand strand0 = 0;
//         const OrientedReadId orientedReadId0(readId0, strand0);

//         std::map<uint64_t, AlignmentPositionBaseStats> positionStatsOnOrientedReadId0;
//         std::map<uint64_t, AlignmentPositionBaseStats> potentialHetSitesOnOrientedReadId0;

//         // // Find the alignments that this oriented read is involved in, with the proper orientation.
//         // const vector< pair<OrientedReadId, uint64_t> > alignments =
//         //     findOrientedAlignmentsPlusIds(orientedReadId0);

//         const auto startingTime = std::chrono::steady_clock::now();
//         std::chrono::steady_clock::duration compressedAlignmentsTime{};
//         std::chrono::steady_clock::duration projectedAlignmentsTime{};
//         std::chrono::steady_clock::duration rleLoopsTime{};

//         std::chrono::steady_clock::time_point t0;
//         std::chrono::steady_clock::time_point t1;
//         std::chrono::steady_clock::time_point t2;
//         std::chrono::steady_clock::time_point t3;

//         // 1.e-9* double((std:chrono:duration_cast<std::chrono::nanoseconds>(t1-t0).count())) << "s.";

//         // Loop over alignment involving this read, as stored in the
//         // alignment table.
//         const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
//         for(const auto alignmentId: alignmentTable0) {
//             const AlignmentData& ad = alignmentData[alignmentId];

//             // Get the oriented read ids that the AlignmentData refers to.
//             OrientedReadId orientedReadId0(ad.readIds[0], 0);
//             OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
//             AlignmentInfo alignmentInfo = ad.info;

//             cout << "Working on alignment between the reads: " << endl;
//             cout << "ReadId0: " << orientedReadId0.getReadId() << " Strand: " << orientedReadId0.getStrand() << endl;
//             cout << "ReadId1: " << orientedReadId1.getReadId() << " Strand: " << orientedReadId1.getStrand() << endl;

//             // Now we need to check if the compressed alignment needs modifications to
//             // represent the alignments between the correct oriented reads.
//             Alignment alignment;

//             t0 = std::chrono::steady_clock::now();

//             // The alignment is stored in compressed form as a string,
//             // so we have to decompress it.
//             span<const char> compressedAlignment = compressedAlignments[alignmentId];
//             shasta::decompress(compressedAlignment, alignment);

//             t1 = std::chrono::steady_clock::now();

//             compressedAlignmentsTime += t1 - t0;

//             // Swap oriented reads, if necessary.
//             if(orientedReadId0.getReadId() != readId0) {
//                 swap(orientedReadId0, orientedReadId1);
//                 alignmentInfo.swap();
//                 alignment.swap();
//             }
//             SHASTA_ASSERT(orientedReadId0.getReadId() == readId0);

//             // Get the markerCount info from the alignment
//             uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());
//             uint32_t markerCount1 = uint32_t(markers[orientedReadId1.getValue()].size());

//             // We need to reverseComplement the ordinals of the bases in the alignment
//             if (orientedReadId0.getStrand() != strand0) {
//                 alignment.reverseComplement(markerCount0, markerCount1);
//                 alignment.checkStrictlyIncreasing();
//                 orientedReadId0.flipStrand();
//                 orientedReadId1.flipStrand();
//                 alignmentInfo.reverseComplement();
//             }
//             SHASTA_ASSERT(orientedReadId0.getStrand() == strand0);


//             // Project this alignment to base space.
//             const ProjectedAlignment projectedAlignment(
//                 *this,
//                 {orientedReadId0, orientedReadId1},
//                 alignment,
//                 false);

//             cout << "Projected alignment between ReadId: " << orientedReadId0.getReadId() << "-" << orientedReadId0.getStrand() << " and ReadId: " << orientedReadId1.getReadId() << "-" << orientedReadId1.getStrand() << " completed." << endl;
            
//             t2 = std::chrono::steady_clock::now();

//             projectedAlignmentsTime += t2 - t1;

//             // Loop over the RAW and RLE segments of the projected alignment.
//             for(const ProjectedAlignmentSegment& segment: projectedAlignment.segments) {
                
//                 // Get the RLE sequences
//                 const vector<Base>& rleSequence0 = segment.rleSequences[0];
//                 const vector<Base>& rleSequence1 = segment.rleSequences[1];

//                 // Get the RAW sequences
//                 const vector<Base>& sequence0 = segment.sequences[0];
//                 const vector<Base>& sequence1 = segment.sequences[1];
                
//                 // Align them base by base to get the right sequence alignment representation
//                 uint64_t position0 = 0;
//                 uint64_t position1 = 0;
//                 vector<AlignedBase> rawAlignmentSequence0;
//                 vector<AlignedBase> rawAlignmentSequence1;
                
//                 // Retreive the raw alignment sequences
//                 for(const pair<bool, bool>& p: segment.alignment) {
//                     const bool hasBase0 = p.first;
//                     const bool hasBase1 = p.second;
    
//                     if(hasBase0) {
//                         rawAlignmentSequence0.push_back(AlignedBase(sequence0[position0++]));
//                     } else {
//                         rawAlignmentSequence0.push_back(AlignedBase::gap());
//                     }
    
//                     if(hasBase1) {
//                         rawAlignmentSequence1.push_back(AlignedBase(sequence1[position1++]));
//                     } else {
//                         rawAlignmentSequence1.push_back(AlignedBase::gap());
//                     }
    
//                 }
//                 SHASTA_ASSERT(rawAlignmentSequence0.size() == rawAlignmentSequence1.size());



//                 position0 = 0;
//                 position1 = 0;
//                 vector<AlignedBase> rleAlignmentSequence0;
//                 vector<AlignedBase> rleAlignmentSequence1;

//                 // Retreive the RLE alignment sequences
//                 for(const pair<bool, bool>& p: segment.rleAlignment) {
//                     const bool hasBase0 = p.first;
//                     const bool hasBase1 = p.second;
    
//                     if(hasBase0) {
//                         rleAlignmentSequence0.push_back(AlignedBase(rleSequence0[position0++]));
//                     } else {
//                         rleAlignmentSequence0.push_back(AlignedBase::gap());
//                     }
    
//                     if(hasBase1) {
//                         rleAlignmentSequence1.push_back(AlignedBase(rleSequence1[position1++]));
//                     } else {
//                         rleAlignmentSequence1.push_back(AlignedBase::gap());
//                     }
    
//                 }
//                 SHASTA_ASSERT(rleAlignmentSequence0.size() == rleAlignmentSequence1.size());


                
//                 // Print the RAW and RLE alignment sequences for debugging
//                 // if the edit distance in RLE is not 0
//                 if(segment.rleEditDistance != 0) {
//                     cout << "We found " << segment.editDistance << " edit distance differences between the RAW sequences" << endl;
//                     cout << "We found " << segment.rleEditDistance << " edit distance differences between the RLE sequences" << endl;
//                     cout << "rawAlignmentSequence0: ";
//                     for(uint64_t i=0; i<rawAlignmentSequence0.size(); i++) {
//                         const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
//                         if(isDifferent) {
//                             cout << "[";
//                         }
//                         cout << rawAlignmentSequence0[i].character();
//                         if(isDifferent) {
//                             cout << "]";
//                         }
//                     }
//                     cout << endl;
                
//                     cout << "rawAlignmentSequence1: ";
//                     for(uint64_t i=0; i<rawAlignmentSequence1.size(); i++) {
//                         const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
//                         if(isDifferent) {
//                             cout << "[";
//                         }
//                         cout << rawAlignmentSequence1[i].character();
//                         if(isDifferent) {
//                             cout << "]";
//                         }
//                     }
//                     cout << endl;

//                     // cout << "rleAlignmentSequence0: ";
//                     // for(uint64_t i=0; i<rleAlignmentSequence0.size(); i++) {
//                     //     const bool isDifferent = (rleAlignmentSequence0[i] != rleAlignmentSequence1[i]);
//                     //     if(isDifferent) {
//                     //         cout << "[";
//                     //     }
//                     //     cout << rleAlignmentSequence0[i].character();
//                     //     if(isDifferent) {
//                     //         cout << "]";
//                     //     }
//                     // }
//                     // cout << endl;
                
//                     // cout << "rleAlignmentSequence1: ";
//                     // for(uint64_t i=0; i<rleAlignmentSequence1.size(); i++) {
//                     //     const bool isDifferent = (rleAlignmentSequence0[i] != rleAlignmentSequence1[i]);
//                     //     if(isDifferent) {
//                     //         cout << "[";
//                     //     }
//                     //     cout << rleAlignmentSequence1[i].character();
//                     //     if(isDifferent) {
//                     //         cout << "]";
//                     //     }
//                     // }
//                     cout << endl;
//                     cout << "rawAlignmentSequence0 positionsA start: " << segment.positionsA[0] << endl;
//                     cout << "rawAlignmentSequence1 positionsA start: " << segment.positionsA[1] << endl;
//                 }




//                 //
//                 //
//                 //
//                 // UNCOMMENT IF WE WANT TO WORK IN RLE
//                 //
//                 //
//                 //

//                 // // VERY IMPORTANT. We need to map the positions of the bases in the RLE alignment 
//                 // // to the positions of the bases in the RAW alignment.
//                 // vector<uint64_t> positionsOfRleSegmentInRawCoordinates;
//                 // uint64_t iRLE = 0;
//                 // for(uint64_t iRAW=0; iRAW<rawAlignmentSequence0.size(); iRAW++) {
//                 //     if(rawAlignmentSequence0[iRAW] == rleAlignmentSequence0[iRLE]) {
//                 //         uint64_t adjustedPositionCoordinate = segment.positionsA[0] + iRAW;
//                 //         positionsOfRleSegmentInRawCoordinates.push_back(adjustedPositionCoordinate);
//                 //         iRLE++;
//                 //     } else {
//                 //         continue;
//                 //     }
//                 // }


//                 // for(uint64_t i=0; i<rleAlignmentSequence0.size(); i++) {

//                 //     // const bool isDifferent = (rleAlignmentSequence0[i] != rleAlignmentSequence1[i]);
//                 //     // if(not isDifferent) {
//                 //     //     continue;
//                 //     // }
                
//                 //     uint64_t positionInRead0 = positionsOfRleSegmentInRawCoordinates[i];

//                 //     if(not positionStatsOnOrientedReadId0.contains(positionInRead0)) {
//                 //         AlignmentPositionBaseStats thisPositionStats;
//                 //         thisPositionStats.positionInRead0 = positionInRead0;
//                 //         thisPositionStats.baseOfReadId0 = static_cast<uint64_t>(rleAlignmentSequence0[i].value);
//                 //         positionStatsOnOrientedReadId0.emplace(positionInRead0, thisPositionStats);

//                 //         // cout << "Position in read0: " << positionInRead0 << " Base of readId0: " << thisPositionStats.baseOfReadId0 << endl;
//                 //     }

//                 //     // Now update the AlignmentPositionBaseStats for this position.
//                 //     // The byte value is always one of 0, 1, 2, 3, 4.
//                 //     // The byte value is always one of A, C, G, T, -.
//                 //     if(static_cast<uint64_t>(rleAlignmentSequence1[i].value) == 0) {
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfA++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithA.push_back(orientedReadId1);
//                 //         // positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithA.push_back(alignmentId);
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].percentageOfA = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfA) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                 //     } else if(static_cast<uint64_t>(rleAlignmentSequence1[i].value) == 1) {
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfC++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithC.push_back(orientedReadId1);
//                 //         // positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithC.push_back(alignmentId);
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].percentageOfC = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfC) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                 //     } else if(static_cast<uint64_t>(rleAlignmentSequence1[i].value) == 2) {
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfG++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithG.push_back(orientedReadId1);
//                 //         // positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithG.push_back(alignmentId);
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].percentageOfG = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfG) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                 //     } else if(static_cast<uint64_t>(rleAlignmentSequence1[i].value) == 3) {
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfT++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithT.push_back(orientedReadId1);
//                 //         // positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithT.push_back(alignmentId);
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].percentageOfT = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfT) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                 //     } else if(static_cast<uint64_t>(rleAlignmentSequence1[i].value) == 4) {
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfGap++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithGap.push_back(orientedReadId1);
//                 //         // positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithGap.push_back(alignmentId);
//                 //         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                 //         positionStatsOnOrientedReadId0[positionInRead0].percentageOfGap = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfGap) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                 //     }

//                 //     // cout << "Position in read0: " << positionInRead0 << " Base of readId0: " << static_cast<uint64_t>(rleAlignmentSequence0[i].value) << " Base of readId1: " << static_cast<uint64_t>(rleAlignmentSequence1[i].value) << endl;
//                 //     // cout << "We found a editDistance variance in Position in read0: " << positionInRead0 << endl;
//                 // }

//                 //
//                 //
//                 //
//                 // UNCOMMENT IF WE WANT TO WORK IN RLE
//                 //
//                 //
//                 //



//                 // When there is an insertion on the orientedReadId1, it will appear as a deletion on the
//                 // orientedReadId0. This will mess with the statisticts of the true position on orientedReadId0.
//                 // So we need to check if the base in the rawAlignmentSequence0 is a gap.
//                 uint64_t actualPositionIndex = 0;
//                 uint64_t positionsConsideredInRead0 = 0;
//                 for(uint64_t i=0; i<rawAlignmentSequence0.size(); i++) {

//                     // const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
//                     // if(not isDifferent) {
//                     //     continue;
//                     // }

//                     // Check if the base in the rawAlignmentSequence0 is a gap
//                     if(rawAlignmentSequence0[i].isGap()) {
//                         continue;
//                     }
                    
//                     uint64_t positionInRead0;

//                     // Edge case where the deletion starts at the beginning position of the segment
//                     if (positionsConsideredInRead0 == 0) {
//                         positionInRead0 = segment.positionsA[0] + actualPositionIndex;
//                         positionsConsideredInRead0++;
//                     } else {
//                         actualPositionIndex++;
//                         positionInRead0 = segment.positionsA[0] + actualPositionIndex;
//                         positionsConsideredInRead0++;
//                     }
//                     SHASTA_ASSERT(positionsConsideredInRead0 == actualPositionIndex + 1);
                
                    

//                     if(not positionStatsOnOrientedReadId0.contains(positionInRead0)) {
//                         AlignmentPositionBaseStats thisPositionStats;
//                         thisPositionStats.positionInRead0 = positionInRead0;
//                         thisPositionStats.baseOfReadId0 = static_cast<uint64_t>(rawAlignmentSequence0[i].value);
//                         positionStatsOnOrientedReadId0.emplace(positionInRead0, thisPositionStats);

//                         // cout << "Position in read0: " << positionInRead0 << " Base of readId0: " << thisPositionStats.baseOfReadId0 << endl;
//                     }

//                     // Now update the AlignmentPositionBaseStats for this position.
//                     // The byte value is always one of 0, 1, 2, 3, 4.
//                     // The byte value is always one of A, C, G, T, -.
//                     if(rawAlignmentSequence1[i].value == 0) {
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfA++;
//                         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithA.insert(orientedReadId1);
//                         positionStatsOnOrientedReadId0[positionInRead0].readIdsWithA.insert(orientedReadId1.getReadId());
//                         positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithA.insert(alignmentId);
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                         positionStatsOnOrientedReadId0[positionInRead0].percentageOfA = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfA) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                     } else if(rawAlignmentSequence1[i].value == 1) {
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfC++;
//                         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithC.insert(orientedReadId1);
//                         positionStatsOnOrientedReadId0[positionInRead0].readIdsWithC.insert(orientedReadId1.getReadId());
//                         positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithC.insert(alignmentId);
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                         positionStatsOnOrientedReadId0[positionInRead0].percentageOfC = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfC) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                     } else if(rawAlignmentSequence1[i].value == 2) {
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfG++;
//                         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithG.insert(orientedReadId1);
//                         positionStatsOnOrientedReadId0[positionInRead0].readIdsWithG.insert(orientedReadId1.getReadId());
//                         positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithG.insert(alignmentId);
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                         positionStatsOnOrientedReadId0[positionInRead0].percentageOfG = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfG) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                     } else if(rawAlignmentSequence1[i].value == 3) {
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfT++;
//                         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithT.insert(orientedReadId1);
//                         positionStatsOnOrientedReadId0[positionInRead0].readIdsWithT.insert(orientedReadId1.getReadId());
//                         positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithT.insert(alignmentId);
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                         positionStatsOnOrientedReadId0[positionInRead0].percentageOfT = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfT) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                     } else if(rawAlignmentSequence1[i].value == 4) {
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfGap++;
//                         positionStatsOnOrientedReadId0[positionInRead0].orientedReadIdsWithGap.insert(orientedReadId1);
//                         positionStatsOnOrientedReadId0[positionInRead0].readIdsWithGap.insert(orientedReadId1.getReadId());
//                         positionStatsOnOrientedReadId0[positionInRead0].alignmentIdsWithGap.insert(alignmentId);
//                         positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments++;
//                         positionStatsOnOrientedReadId0[positionInRead0].percentageOfGap = double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfGap) / double(positionStatsOnOrientedReadId0[positionInRead0].totalNumberOfAlignments);
//                     }

//                     // cout << "Position in read0: " << positionInRead0 << " Base of readId0: " << static_cast<uint64_t>(rawAlignmentSequence0[i].value) << " Base of readId1: " << static_cast<uint64_t>(rawAlignmentSequence1[i].value) << endl;
//                     // cout << "We found a editDistance variance in Position in read0: " << positionInRead0 << endl;
//                 }


    
//             }

//             t3 = std::chrono::steady_clock::now();

//             rleLoopsTime += t3 - t2;

//             // cout << "Loop over the RLE segments of the Projected alignment between ReadId: " << orientedReadId0.getReadId() << "-" << orientedReadId0.getStrand() << " and ReadId: " << orientedReadId1.getReadId() << "-" << orientedReadId1.getStrand() << " completed." << endl;
//             cout << "Loop over the RAW segments of the Projected alignment between ReadId: " << orientedReadId0.getReadId() << "-" << orientedReadId0.getStrand() << " and ReadId: " << orientedReadId1.getReadId() << "-" << orientedReadId1.getStrand() << " completed." << endl;           

//         }

//         cout << "Time taken to decompress alignment: " << std::chrono::duration_cast<std::chrono::nanoseconds>(compressedAlignmentsTime).count() << " ns" << endl;
//         cout << "Time taken to loop do the projectedAlignments: " << std::chrono::duration_cast<std::chrono::nanoseconds>(projectedAlignmentsTime).count() << " ns" << endl;
//         cout << "Time taken to loop over segments: " << std::chrono::duration_cast<std::chrono::nanoseconds>(rleLoopsTime).count() << " ns" << endl;
//         //return;













//         // We finized analyzing all alignments for this read. Now we can check for potential het sites.
//         // Now we need to check each potential site in readId0 (in positionStatsOnOrientedReadId0) and check if it contains a heterozygous site
//         uint64_t sitesSkippedDueToInsufficientCoverage = 0;
//         for (const auto& [positionInRead0, positionStats] : positionStatsOnOrientedReadId0) {
//             // Skip positions with insufficient coverage
//             // At least ~5x per haplotype
//             if (positionStats.totalNumberOfAlignments < 8) {
//                 sitesSkippedDueToInsufficientCoverage++;
//                 continue;
//             }

//             // Get the number of each base
//             const uint64_t numberOfA = positionStats.totalNumberOfA;
//             const uint64_t numberOfC = positionStats.totalNumberOfC;
//             const uint64_t numberOfG = positionStats.totalNumberOfG;
//             const uint64_t numberOfT = positionStats.totalNumberOfT;
//             const uint64_t numberOfGap = positionStats.totalNumberOfGap;

//             // Sort the base counts in descending order
//             vector<pair<uint64_t, uint64_t>> baseCounts = {
//                 {numberOfA, 0},
//                 {numberOfC, 1},
//                 {numberOfG, 2},
//                 {numberOfT, 3},
//                 {numberOfGap, 4}
//             };
//             std::sort(baseCounts.begin(), baseCounts.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
            
//             // Check if this is a potential heterozygous site
//             // Criteria: top two bases (with highest counts) each represent at least 20% of reads
//             // and together they represent at least 70% of reads
//             if (baseCounts[0].first > 0 && baseCounts[1].first > 0) {
//                 const double firstBasePercentage = double(baseCounts[0].first) / double(positionStats.totalNumberOfAlignments);
//                 const uint64_t firstBaseCounts = baseCounts[0].first;
//                 const double secondBasePercentage = double(baseCounts[1].first) / double(positionStats.totalNumberOfAlignments);
//                 const uint64_t secondBaseCounts = baseCounts[1].first;
//                 const double combinedPercentage = firstBasePercentage + secondBasePercentage;
                
//                 //if (firstBasePercentage >= 0.2 && firstBaseCounts >=3 && secondBasePercentage >= 0.2 && secondBaseCounts >=3 && combinedPercentage >= 0.7) {
//                 if (firstBaseCounts >=4 && secondBaseCounts >=4) {
//                     // This is a potential heterozygous site

//                     // If the base of the readId0 is not one of the top two het bases then we don't want to include it as potential het site
//                     if (positionStats.baseOfReadId0 != baseCounts[0].second && positionStats.baseOfReadId0 != baseCounts[1].second) {
//                         continue;
//                     }

//                     potentialHetSitesOnOrientedReadId0[positionInRead0] = positionStats;
//                     potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1 = baseCounts[0].second;
//                     potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2 = baseCounts[1].second;
//                     if (baseCounts[0].second == 0) {
//                         // Base A
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds = positionStats.orientedReadIdsWithA;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1ReadIds = positionStats.readIdsWithA;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1Alignments = positionStats.alignmentIdsWithA;
//                     } else if (baseCounts[0].second == 1) {
//                         // Base C
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds = positionStats.orientedReadIdsWithC;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1ReadIds = positionStats.readIdsWithC;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1Alignments = positionStats.alignmentIdsWithC;
//                     } else if (baseCounts[0].second == 2) {
//                         // Base G
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds = positionStats.orientedReadIdsWithG;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1ReadIds = positionStats.readIdsWithG;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1Alignments = positionStats.alignmentIdsWithG;
//                     } else if (baseCounts[0].second == 3) {
//                         // Base T
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds = positionStats.orientedReadIdsWithT;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1ReadIds = positionStats.readIdsWithT;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1Alignments = positionStats.alignmentIdsWithT;
//                     } else if (baseCounts[0].second == 4) {
//                         // Gap
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1OrientedReadIds = positionStats.orientedReadIdsWithGap;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1ReadIds = positionStats.readIdsWithGap;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase1Alignments = positionStats.alignmentIdsWithGap;
//                     }

//                     if (baseCounts[1].second == 0) {
//                         // Base A
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds = positionStats.orientedReadIdsWithA;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2ReadIds = positionStats.readIdsWithA;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2Alignments = positionStats.alignmentIdsWithA;
//                     } else if (baseCounts[1].second == 1) {
//                         // Base C
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds = positionStats.orientedReadIdsWithC;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2ReadIds = positionStats.readIdsWithC;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2Alignments = positionStats.alignmentIdsWithC;
//                     } else if (baseCounts[1].second == 2) {
//                         // Base G
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds = positionStats.orientedReadIdsWithG;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2ReadIds = positionStats.readIdsWithG;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2Alignments = positionStats.alignmentIdsWithG;
//                     } else if (baseCounts[1].second == 3) {
//                         // Base T
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds = positionStats.orientedReadIdsWithT;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2ReadIds = positionStats.readIdsWithT;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2Alignments = positionStats.alignmentIdsWithT;
//                     } else if (baseCounts[1].second == 4) {
//                         // Gap
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2OrientedReadIds = positionStats.orientedReadIdsWithGap;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2ReadIds = positionStats.readIdsWithGap;
//                         potentialHetSitesOnOrientedReadId0[positionInRead0].hetBase2Alignments = positionStats.alignmentIdsWithGap;
//                     }
//                     potentialHetSitesOnOrientedReadId0[positionInRead0].totalNumberOfHetBase1 = baseCounts[0].first;
//                     potentialHetSitesOnOrientedReadId0[positionInRead0].totalNumberOfHetBase2 = baseCounts[1].first;
//                     potentialHetSitesOnOrientedReadId0[positionInRead0].percentageOfHetBase1 = firstBasePercentage;
//                     potentialHetSitesOnOrientedReadId0[positionInRead0].percentageOfHetBase2 = secondBasePercentage;
                    
//                     // // // Optional: Log information about this site
//                     // char bases[] = {'A', 'C', 'G', 'T', '-'};
//                     // cout << "Potential heterozygous site at position " << positionInRead0 
//                     //     << " in readId " << readId0 << ":" << endl
//                     //     << "  Base of readId0: " << bases[positionStats.baseOfReadId0] << endl
//                     //     << "  First variant: " << bases[baseCounts[0].second] << " (" << firstBasePercentage * 100 << "%)" << endl
//                     //     << "  Second variant: " << bases[baseCounts[1].second] << " (" << secondBasePercentage * 100 << "%)" << endl;
//                 }
//             }
//         }




//         // Now we need to analyze the potential heterozygous sites and try to find sets of reads that belong to the same haplotype

//         if (potentialHetSitesOnOrientedReadId0.size() > 0) {
//             cout << "Found " << potentialHetSitesOnOrientedReadId0.size() << " potential heterozygous sites in readId " << readId0 << endl;
//             cout << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
//         } else {
//             cout << "Found no potential heterozygous sites in readId " << readId0 << endl;
//             cout << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
//             continue;
//         }

//         // Loop over the potential heterozygous sites and print them
//         char bases[] = {'A', 'C', 'G', 'T', '-'};
//         for (const auto& [positionInRead0, positionStats] : potentialHetSitesOnOrientedReadId0) {
//             if (positionStats.hetBase1 == 4 || positionStats.hetBase2 == 4) {
//                 // Skip gaps
//                 continue;
//             }
//             cout << "Potential heterozygous site at position " << positionInRead0 
//                 << " in readId " << readId0 << " in strand " << strand0 << ":" << endl
//                 << "  Base of readId0: " << bases[positionStats.baseOfReadId0] << endl
//                 << "  First variant: " << bases[positionStats.hetBase1] << " (" << positionStats.totalNumberOfHetBase1 << " "<< positionStats.percentageOfHetBase1 * 100 << "%)" << endl
//                 << "  Second variant: " << bases[positionStats.hetBase2] << " (" << positionStats.totalNumberOfHetBase2 << " " << positionStats.percentageOfHetBase2 * 100 << "%)" << endl;
            
//             cout << "ReadIds that support the first variant: " << endl;
//             std::set<OrientedReadId> orientedReadIdsHetBase1;
//             if (positionStats.hetBase1 == 0) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithA) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase1.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase1 == 1) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithC) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase1.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase1 == 2) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithG) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase1.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase1 == 3) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithT) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase1.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase1 == 4) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithGap) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase1.insert(orientedReadId);
//                 }
//             }
            

//             cout << "ReadIds that support the second variant: " << endl;
//             std::set<OrientedReadId> orientedReadIdsHetBase2;
//             if (positionStats.hetBase2 == 0) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithA) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase2.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase2 == 1) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithC) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase2.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase2 == 2) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithG) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase2.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase2 == 3) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithT) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase2.insert(orientedReadId);
//                 }
//             } else if (positionStats.hetBase2 == 4) {
//                 for (const auto& orientedReadId : positionStats.orientedReadIdsWithGap) {
//                     cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//                     orientedReadIdsHetBase2.insert(orientedReadId);
//                 }
//             }


//             // Check if orientedReadIdsHetBase1 and orientedReadIdsHetBase2 have any common reads
//             // There should NOT be any common reads between the two sets
//             // Helper function to check for common elements between sets
//             auto hasCommonElements = [](const std::set<OrientedReadId>& set1, 
//                 const std::set<OrientedReadId>& set2) -> bool
//             {
//                 // Use the smaller set for iteration
//                 const auto& smaller = (set1.size() < set2.size()) ? set1 : set2;
//                 const auto& larger = (set1.size() < set2.size()) ? set2 : set1;

//                 for (const auto& elem : smaller) {
//                     if (larger.find(elem) != larger.end()) {
//                         return true;  // Found common element
//                     }
//                 }
//                 return false;  // No common elements
//             };

//             SHASTA_ASSERT(hasCommonElements(orientedReadIdsHetBase1, orientedReadIdsHetBase2) == false);


//         }

        

//         // return;


//         // Remove potential het sites that have a base deletion in either allele 
//         // if there is a deletion in the previous or next base.
//         // Most of them correspond to homopolymer run errors
//         vector<uint64_t> positionsToRemove;
//         for (const auto& [positionInRead0, positionStats] : potentialHetSitesOnOrientedReadId0) {
//             if ((positionStats.hetBase1 == 4) or (positionStats.hetBase2 == 4)) {
//                 positionsToRemove.push_back(positionInRead0);
//             }
//         }

//         // for (const auto& [positionInRead0, positionStats] : potentialHetSitesOnOrientedReadId0) {
//         //     if ((positionStats.hetBase1 == 4) || (positionStats.hetBase2 == 4)) {
//         //         // Check if the next position also has a gap
//         //         auto nextPositionIter = potentialHetSitesOnOrientedReadId0.find(positionInRead0 + 2);
//         //         auto previousPositionIter = potentialHetSitesOnOrientedReadId0.find(positionInRead0 - 2);
//         //         if (nextPositionIter != potentialHetSitesOnOrientedReadId0.end() and 
//         //             (nextPositionIter->second.hetBase1 == 4 or nextPositionIter->second.hetBase2 == 4)) {
//         //             // Next position has a gap, so we remove this position
//         //             positionsToRemove.push_back(positionInRead0);
//         //         }
//         //         if (previousPositionIter != potentialHetSitesOnOrientedReadId0.end() and 
//         //             (previousPositionIter->second.hetBase1 == 4 or previousPositionIter->second.hetBase2 == 4)) {
//         //             // Previous position has a gap, so we remove this position
//         //             positionsToRemove.push_back(positionInRead0);
//         //         }
//         //     }
//         // }

//         for (const auto& positionInRead0 : positionsToRemove) {
//             potentialHetSitesOnOrientedReadId0.erase(positionInRead0);
//         }

//         // Print to verify that we removed the right positions
//         if (potentialHetSitesOnOrientedReadId0.size() > 0) {
//             cout << "Found " << potentialHetSitesOnOrientedReadId0.size() << " potential heterozygous sites in readId " << readId0 << ", after removing sites with deletions" << endl;
//             cout << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
//         } else {
//             cout << "Found no potential heterozygous sites in readId " << readId0 << ", after removing sites with deletions" << endl;
//             cout << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
//             continue;
//         }





//         // --- Dynamic Programming for Grouping Compatible Sites ---
//         cout << timestamp << "Starting DP for grouping compatible sites for read " << orientedReadId0 << endl;

//         // 0. Prepare sorted list of sites
//         std::vector<std::pair<uint64_t, AlignmentPositionBaseStats>> sortedSites;
//         for (const auto& pair : potentialHetSitesOnOrientedReadId0) {
//             sortedSites.push_back(pair);
//         }
//         // Ensure sites are sorted by position (map iteration order is usually sorted, but explicit sort is safer)
//         std::sort(sortedSites.begin(), sortedSites.end(),
//                     [](const auto& a, const auto& b) { return a.first < b.first; });

//         uint64_t N = sortedSites.size();
//         if (N == 0) {
//             cout << "No potential het sites found after filtering for read " << orientedReadId0 << ", skipping DP and subsequent phasing." << endl;
//             cout << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
//             continue; // Skip to the next readId if no sites
//         } else {
//             cout << "Found " << N << " potential heterozygous sites for DP in readId " << readId0 << endl;
//             cout << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
//         }


//         // Define compatibility check function locally (or move to class scope)
//         auto areCompatibleSites =
//             [&](const std::pair<uint64_t, AlignmentPositionBaseStats>& site_pair_i,
//                 const std::pair<uint64_t, AlignmentPositionBaseStats>& site_pair_j,
//                 const OrientedReadId& targetReadId) -> bool
//         {
//             const AlignmentPositionBaseStats& stats_i = site_pair_i.second;
//             const AlignmentPositionBaseStats& stats_j = site_pair_j.second;

//             // --- Check if target read covers both sites informatively ---
//             // Target read must have one of the two identified heterozygous bases at each site
//             // for a meaningful phase comparison relative to the target.
//             bool target_covers_i = (stats_i.baseOfReadId0 == stats_i.hetBase1 || stats_i.baseOfReadId0 == stats_i.hetBase2);
//             bool target_covers_j = (stats_j.baseOfReadId0 == stats_j.hetBase1 || stats_j.baseOfReadId0 == stats_j.hetBase2);

//             if (!target_covers_i || !target_covers_j) {
//                 // Target read doesn't provide a consistent reference phase across both sites.
//                 cout << "Target read " << targetReadId << " does not cover sites " << site_pair_i.first << " and " << site_pair_j.first << " informatively for relative phasing." << endl;
//                 return false;
//             }
//             // --- Target read covers both sites informatively ---

//             // 1a. Get all reads overlapping site i
//             std::set<OrientedReadId> reads_at_i = stats_i.hetBase1OrientedReadIds;
//             reads_at_i.insert(stats_i.hetBase2OrientedReadIds.begin(), stats_i.hetBase2OrientedReadIds.end());
//             // Add target read since we know it covers informatively
//             //reads_at_i.insert(targetReadId);

//             // 1b. Get all reads overlapping site j
//             std::set<OrientedReadId> reads_at_j = stats_j.hetBase1OrientedReadIds;
//             reads_at_j.insert(stats_j.hetBase2OrientedReadIds.begin(), stats_j.hetBase2OrientedReadIds.end());
//             // Add target read since we know it covers informatively
//             //reads_at_j.insert(targetReadId);

//             // 2. Find common overlapping reads
//             std::vector<OrientedReadId> commonOverlappingReads;
//             std::set_intersection(reads_at_i.begin(), reads_at_i.end(),
//                                     reads_at_j.begin(), reads_at_j.end(),
//                                     std::back_inserter(commonOverlappingReads));

//             bool found_phase_00 = false; // Flag for reads matching target at both sites (Phase 0-0)
//             bool found_phase_11 = false; // Flag for reads differing from target at both sites (Phase 1-1)

//             // 3. Check phase consistency for each common read relative to target 
//             //    and also evidence for both phase patterns
//             for (const auto& read : commonOverlappingReads) {
//                 // Determine phase at site i (0: matches target, 1: differs from target, -1: unknown/not informative)
//                 int phase_i = -1;
//                 uint64_t allele_at_i = 100; // Invalid allele
//                 if (stats_i.hetBase1OrientedReadIds.count(read)) {
//                     allele_at_i = stats_i.hetBase1;
//                 } else if (stats_i.hetBase2OrientedReadIds.count(read)) {
//                     allele_at_i = stats_i.hetBase2;
//                 }

//                 if (allele_at_i <= 4) { // Is the allele valid (A,C,G,T,-)?
//                     if (allele_at_i == stats_i.baseOfReadId0) {
//                         phase_i = 0; // Matches the base of the target read
//                     } else {
//                         phase_i = 1; // Differs from the base of the target read
//                     }
//                 } // else phase_i remains -1 (read doesn't support either het allele at this site)

//                 cout << "Read " << read << " at site " << site_pair_i.first << " has phase " << phase_i << endl;

//                 // Determine phase at site j (0: matches target, 1: differs from target, -1: unknown/not informative)
//                 int phase_j = -1;
//                 uint64_t allele_at_j = 100;
//                 if (stats_j.hetBase1OrientedReadIds.count(read)) {
//                     allele_at_j = stats_j.hetBase1;
//                 } else if (stats_j.hetBase2OrientedReadIds.count(read)) {
//                     allele_at_j = stats_j.hetBase2;
//                 }

//                 if (allele_at_j <= 4) {
//                     if (allele_at_j == stats_j.baseOfReadId0) {
//                         phase_j = 0; // Matches target
//                     } else {
//                         phase_j = 1; // Differs from target
//                     } 
//                 } // else phase_j remains -1

//                 cout << "Read " << read << " at site " << site_pair_j.first << " has phase " << phase_j << endl;

//                 // Check for inconsistency or update flags if consistent and informative
//                 if (phase_i != -1 && phase_j != -1) { // Both phases defined relative to target for this read
//                     if (phase_i != phase_j) {
//                         // Found a read with inconsistent phasing relative to the target across the two sites.
//                         cout << "Incompatibility: Read " << read << " phase " << phase_i << " at site " << site_pair_i.first << ", phase " << phase_j << " at site " << site_pair_j.first << " relative to target." << endl;
//                         return false; // Inconsistent phasing detected
//                     } else if (phase_i == 0) { // phase_i == 0 && phase_j == 0
//                         found_phase_00 = true; // Found evidence for 0-0 linkage relative to target
//                     } else { // phase_i == 1 && phase_j == 1
//                         found_phase_11 = true; // Found evidence for 1-1 linkage relative to target
//                     }
//                 } else if (phase_i == -1 || phase_j == -1) {
//                     cout << "Read " << read << " is not informative for phasing relative to target at site " << site_pair_i.first << " or site " << site_pair_j.first << "." << endl;
//                     return false; // Inconsistent phasing detected
//                 }

//             }

//             // Sites are compatible *relative to the target* only if:
//             // 1. No common reads showed inconsistent phasing between positions i and j (phase_i != phase_j).
//             // 2. There was at least one read supporting the 0-0 linkage relative to target.
//             // 3. There was at least one read supporting the 1-1 linkage relative to target.
//             // cout << "Compatibility check for sites " << site_pair_i.first << " and " << site_pair_j.first << " relative to target " << targetReadId << ": found_phase_00=" << found_phase_00 << ", found_phase_11=" << found_phase_11 << endl;
//             return found_phase_00 && found_phase_11;

//         };


//         // 1. Initialize DP table and parent pointers
//         std::vector<uint64_t> LCG(N, 1); // LCG[i] = size of largest compatible group ending at i
//         std::vector<int64_t> parent(N, -1); // parent[i] = index j that gave the max LCG[i]

//         // 2. Fill DP table using recurrence relation: LCG(i) = max_{j<i, S[j]<->S[i]} {LCG(j)} + 1
//         for (uint64_t i = 1; i < N; ++i) {
//             for (uint64_t j = 0; j < i; ++j) {
//                 // Check compatibility S[j] <-> S[i]
//                 if (areCompatibleSites(sortedSites[i], sortedSites[j], orientedReadId0)) {
//                     if (LCG[j] + 1 > LCG[i]) {
//                         LCG[i] = LCG[j] + 1;
//                         parent[i] = j;
//                     }
//                 }
//             }
//         }

//         // 3. Traceback for grouping sites
//         std::vector<bool> isAssigned(N, false);

//         // Create pairs of (LCG value, index) to sort by LCG value descending
//         std::vector<std::pair<int, size_t>> sortedLCGIndices;
//         for(size_t i = 0; i < N; ++i) {
//             if (LCG[i] > 1) { // Only consider sites part of a group > 1
//                 sortedLCGIndices.push_back({LCG[i], i});
//             }
//         }
//         // Sort descending by LCG value
//         std::sort(sortedLCGIndices.rbegin(), sortedLCGIndices.rend()); 

//         cout << "DP results for read " << orientedReadId0 << ": LCG values > 1:" << endl;
//         for(const auto& p : sortedLCGIndices) {
//             cout << "  Site Index: " << p.second << " (Pos: " << sortedSites[p.second].first << "), LCG: " << p.first << endl;
//         }

//         std::set<OrientedReadId> tempInPhaseOrientedReads;
//         std::set<OrientedReadId> excludedOutOfPhaseOrientedReads;
//         std::set<OrientedReadId> involvedOrientedReadsInInformativeSites;
//         std::set<OrientedReadId> finalInPhaseOrientedReads;

//         for (const auto& lcgPair : sortedLCGIndices) {
//             size_t current_i = lcgPair.second;

//             if (!isAssigned[current_i]) {
//                 std::vector<uint64_t> newClusterPositions; // Keep this for informativeSites logic later if needed
//                 std::vector<size_t> newClusterIndices;    // Store indices for phasing logic
//                 int traceIndex = current_i;
//                 cout << "Starting traceback for cluster ending at index " << current_i << " (Pos: " << sortedSites[current_i].first << ") with LCG " << LCG[current_i] << endl;

//                 while (traceIndex != -1 && !isAssigned[traceIndex]) {

//                     const AlignmentPositionBaseStats& currentIndexStats = sortedSites[traceIndex].second;
//                     uint64_t targetAllele = currentIndexStats.baseOfReadId0;
//                     uint64_t otherAllele = (targetAllele == currentIndexStats.hetBase1) ? currentIndexStats.hetBase2 : currentIndexStats.hetBase1;

//                     // 1. Populate involved reads for this site
//                     // Add all reads supporting either heterozygous allele at this site
//                     involvedOrientedReadsInInformativeSites.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
//                     involvedOrientedReadsInInformativeSites.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());

//                     // 2. Get in-phase and out-of-phase reads *at this specific site index* relative to the target read
//                     // Populate in-phase set (reads with the same allele as the target read at this site)
//                     if (targetAllele == currentIndexStats.hetBase1) {
//                         tempInPhaseOrientedReads.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
//                     } else { // targetAllele must be stats.hetBase2
//                         tempInPhaseOrientedReads.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
//                     }
//                     tempInPhaseOrientedReads.insert(orientedReadId0); // The target read is always in phase with itself

//                     // Populate out-of-phase set (reads with the other heterozygous allele)
//                     if (otherAllele == currentIndexStats.hetBase1) {
//                         excludedOutOfPhaseOrientedReads.insert(currentIndexStats.hetBase1OrientedReadIds.begin(), currentIndexStats.hetBase1OrientedReadIds.end());
//                     } else { // otherAllele must be stats.hetBase2
//                         excludedOutOfPhaseOrientedReads.insert(currentIndexStats.hetBase2OrientedReadIds.begin(), currentIndexStats.hetBase2OrientedReadIds.end());
//                     }

//                     newClusterPositions.push_back(sortedSites[traceIndex].first); // Store site position
//                     newClusterIndices.push_back(traceIndex); // Store site index
//                     isAssigned[traceIndex] = true;
//                     cout << "  Adding site index " << traceIndex << " (Pos: " << sortedSites[traceIndex].first << ") to compatible group." << endl;
//                     traceIndex = parent[traceIndex];

//                 }

//                 std::reverse(newClusterPositions.begin(), newClusterPositions.end()); // Store cluster in positional order
//                 std::reverse(newClusterIndices.begin(), newClusterIndices.end());   // Store indices in positional order
//                 //siteClusters.push_back(newClusterPositions); // Keep the original structure if needed elsewhere
//                 //siteIndexClusters.push_back(newClusterIndices); // Use this for phasing
//                 cout << "  Finished constructing of a compatible group with " << newClusterPositions.size() << " sites." << endl;
//             }
//         }


//         // Find which reads exist only in the tempInPhaseOrientedReads set and not in the excludedOutOfPhaseOrientedReads set and add them to finalInPhaseOrientedReads
//         for (const auto& read : tempInPhaseOrientedReads) {
//             if (excludedOutOfPhaseOrientedReads.find(read) == excludedOutOfPhaseOrientedReads.end()) {
//                 finalInPhaseOrientedReads.insert(read);
//             }
//         }
//         cout << "Final in-phase reads for read " << orientedReadId0 << ": " << endl;
//         for (const auto& read : finalInPhaseOrientedReads) {
//             cout << "  ReadId: " << read.getReadId() << " Strand: " << read.getStrand() << endl;
//         }






//         // Loop over all alignments and mark those involving only reads within finalInPhaseOrientedReads
//         // Also forbid alignments involving reads in excludedOutOfPhaseOrientedReads
//         for(uint64_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
//             // Skip if already considered to avoid redundant checks or overwriting previous decisions
//             if (alignmentsAlreadyConsidered[alignmentId]) {
//                 continue;
//             }

//             const AlignmentData& ad = alignmentData[alignmentId];
//             ReadId alnReadId0 = ad.readIds[0];
//             ReadId alnReadId1 = ad.readIds[1];
//             bool isSameStrand = ad.isSameStrand;

//             // Create the primary oriented read pair for this alignment
//             OrientedReadId primaryOrientedRead0(alnReadId0, 0);
//             OrientedReadId primaryOrientedRead1(alnReadId1, isSameStrand ? 0 : 1);

//             // Create the reverse complement oriented read pair
//             OrientedReadId rcOrientedRead0(alnReadId0, 1);
//             OrientedReadId rcOrientedRead1(alnReadId1, isSameStrand ? 1 : 0);

//             // --- Check if any read involved in the alignment is in excludedOutOfPhaseOrientedReads ---
//             bool involvesExcludedRead = (finalInPhaseOrientedReads.count(primaryOrientedRead0) && excludedOutOfPhaseOrientedReads.count(primaryOrientedRead1)) ||
//                                         (finalInPhaseOrientedReads.count(primaryOrientedRead1) && excludedOutOfPhaseOrientedReads.count(primaryOrientedRead0)) ||
//                                         (finalInPhaseOrientedReads.count(rcOrientedRead0) && excludedOutOfPhaseOrientedReads.count(rcOrientedRead1)) ||
//                                         (finalInPhaseOrientedReads.count(rcOrientedRead1) && excludedOutOfPhaseOrientedReads.count(rcOrientedRead0));

//             if (involvesExcludedRead) {
//                 cout << "Forbidding alignment " << alignmentId << " involving reads: " << alnReadId0 << " and " << alnReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
//                 forbiddenAlignments[alignmentId] = true;
//                 alignmentsAlreadyConsidered[alignmentId] = true; // Mark as considered
//                 continue; // Skip further checks for this alignment
//             }
//             // --- End check for excluded reads ---


//             // --- Check if both reads are in finalInPhaseOrientedReads (existing logic) ---
//             // Check if both reads in the primary orientation are in finalInPhaseOrientedReads
//             bool bothInPhasePrimary = finalInPhaseOrientedReads.count(primaryOrientedRead0) &&
//                                         finalInPhaseOrientedReads.count(primaryOrientedRead1);

//             // Check if both reads in the reverse complement orientation are in finalInPhaseOrientedReads
//             bool bothInPhaseRC = finalInPhaseOrientedReads.count(rcOrientedRead0) &&
//                                     finalInPhaseOrientedReads.count(rcOrientedRead1);

//             // If both reads of the alignment (in either orientation pair consistent with the alignment)
//             // are present in the finalInPhaseOrientedReads set, mark the alignment.
//             if (bothInPhasePrimary || bothInPhaseRC) {
//                 cout << "Marking alignment " << alignmentId << " involving reads: " << alnReadId0 << " and " << alnReadId1 << " as first pass (intra-phase)" << endl;
//                 firstPassHetAlignments[alignmentId] = true;
//                 alignmentsAlreadyConsidered[alignmentId] = true; // Mark as considered
//             }
//             // --- End check for in-phase reads ---
//         }



//         const long firstPassHetAlignmentsCount = count(firstPassHetAlignments.begin(), firstPassHetAlignments.end(), true);
//         cout << "Found " << firstPassHetAlignmentsCount << " first pass alignments in readId " << readId0 << endl;


        

//         // return;

        


//         // // Handle isolated sites (LCG[i] == 1) - criteria (i) from text
//         // std::vector<uint64_t> informativeIsolatedSites;
//         // // Define the threshold for sufficiently high number of reads for an isolated site
//         // // This value might need tuning based on data characteristics.
//         // const uint64_t minReadsPerHetBaseForIsolatedSite = 4; // Example threshold
//         // for(size_t i = 0; i < N; ++i) {
//         //     if (!isAssigned[i]) { // If not part of any cluster found above (i.e., LCG[i] == 1 or part of a smaller cluster already processed)
//         //         // Check criteria (i): "An isolated site is considered informative only if it is supported by a sufficiently high number of reads"
//         //         const AlignmentPositionBaseStats& stats = sortedSites[i].second;
//         //         // Total reads supporting the two major alleles at this site
//         //         uint64_t supportingReadsHet1 = stats.totalNumberOfHetBase1;
//         //         uint64_t supportingReadsHet2 = stats.totalNumberOfHetBase2;
//         //         if (supportingReadsHet1 >= minReadsPerHetBaseForIsolatedSite && supportingReadsHet2 >= minReadsPerHetBaseForIsolatedSite) {
//         //             informativeIsolatedSites.push_back(sortedSites[i].first);
//         //             // cout << "Site index " << i << " (Pos: " << sortedSites[i].first << ") is isolated but informative (support=" << supportingReads << " >= " << minReadsForIsolatedSite << ")" << endl;
//         //         } else {
//         //             // cout << "Site index " << i << " (Pos: " << sortedSites[i].first << ") is isolated and non-informative (support=" << supportingReads << " < " << minReadsForIsolatedSite << ")" << endl;
//         //         }
//         //     }
//         // }


        


        


























//         // // --- Global Phasing Implementation Starts ---

//         // // 1. Identify Involved Reads
//         // std::set<OrientedReadId> involvedReads;
//         // bool read0Involved = false;
//         // for (const auto& [pos, stats] : potentialHetSitesOnOrientedReadId0) {
//         //     involvedReads.insert(stats.hetBase1OrientedReadIds.begin(), stats.hetBase1OrientedReadIds.end());
//         //     involvedReads.insert(stats.hetBase2OrientedReadIds.begin(), stats.hetBase2OrientedReadIds.end());
//         //     if (stats.baseOfReadId0 == stats.hetBase1 || stats.baseOfReadId0 == stats.hetBase2) {
//         //         read0Involved = true;
//         //     }
//         // }
//         // if (read0Involved) {
//         //     // Check if orientedReadId0 actually exists in the alignment data (should always be true here)
//         //     // and add it if it contributed to any het site definition.
//         //     involvedReads.insert(orientedReadId0);
//         // }

//         // std::set<OrientedReadId> finalHaplotype1Reads;
//         // std::set<OrientedReadId> finalHaplotype2Reads;

//         // if (involvedReads.empty() || potentialHetSitesOnOrientedReadId0.empty()) {
//         //     cout << "No involved reads or sites for global phasing for read " << orientedReadId0 << endl;
//         //     // If no reads/sites, haplotypes remain empty. Skip graph building/partitioning.
//         // } else {
//         //     cout << "Starting global phasing for read " << orientedReadId0 << " with " << involvedReads.size() << " involved reads and " << potentialHetSitesOnOrientedReadId0.size() << " sites." << endl;

//         //     // Create mapping for graph vertices
//         //     std::map<OrientedReadId, size_t> readToIndex;
//         //     std::vector<OrientedReadId> indexToRead;
//         //     size_t currentIndex = 0;
//         //     for(const auto& read : involvedReads) {
//         //         readToIndex[read] = currentIndex;
//         //         indexToRead.push_back(read);
//         //         currentIndex++;
//         //     }
//         //     size_t numInvolvedReads = involvedReads.size();

//         //     // 2. Build Phasing Graph
//         //     // Use boost::adjacency_list. Edge property stores the weight.
//         //     using PhasingGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_weight_t, double>>;
//         //     PhasingGraph phasingGraph(numInvolvedReads);
//         //     // Get the property map for weights
//         //     auto weightMap = boost::get(boost::edge_weight, phasingGraph);

//         //     // Use a map to accumulate weights between pairs of nodes before adding edges
//         //     // Key: pair of indices (u, v) where u < v. Value: accumulated weight.
//         //     std::map<std::pair<size_t, size_t>, double> edgeWeightSums;

//         //     for (const auto& [positionInRead0, stats] : potentialHetSitesOnOrientedReadId0) {
//         //         const auto& reads1 = stats.hetBase1OrientedReadIds;
//         //         const auto& reads2 = stats.hetBase2OrientedReadIds;
//         //         uint64_t allele1 = stats.hetBase1;
//         //         uint64_t allele2 = stats.hetBase2;

//         //         // Combine reads covering this site (including read0 if applicable and involved)
//         //         std::set<OrientedReadId> readsAtSite;
//         //         // Add reads supporting allele 1 if they are in the overall involved set
//         //         for(const auto& r : reads1) {
//         //             if (involvedReads.count(r)) readsAtSite.insert(r);
//         //         }
//         //         // Add reads supporting allele 2 if they are in the overall involved set
//         //         for(const auto& r : reads2) {
//         //             if (involvedReads.count(r)) readsAtSite.insert(r);
//         //         }
//         //         // Add read0 if it supports an allele and is in the overall involved set
//         //         if ((stats.baseOfReadId0 == allele1 || stats.baseOfReadId0 == allele2) && involvedReads.count(orientedReadId0)) {
//         //             readsAtSite.insert(orientedReadId0);
//         //         }


//         //         // Iterate through pairs of reads covering this site
//         //         for (auto itA = readsAtSite.begin(); itA != readsAtSite.end(); ++itA) {
//         //             auto itB = itA;
//         //             std::advance(itB, 1); // Start B after A
//         //             for (; itB != readsAtSite.end(); ++itB) {
//         //                 const OrientedReadId& readA = *itA;
//         //                 const OrientedReadId& readB = *itB;

//         //                 // Determine allele for readA at this site
//         //                 uint64_t alleleA = 100; // Use an invalid value initially
//         //                 if (readA == orientedReadId0) alleleA = stats.baseOfReadId0;
//         //                 else if (reads1.count(readA)) alleleA = allele1;
//         //                 else if (reads2.count(readA)) alleleA = allele2;

//         //                 // Determine allele for readB at this site
//         //                 uint64_t alleleB = 100;
//         //                 if (readB == orientedReadId0) alleleB = stats.baseOfReadId0;
//         //                 else if (reads1.count(readB)) alleleB = allele1;
//         //                 else if (reads2.count(readB)) alleleB = allele2;

//         //                 // Ensure both reads have a valid allele assigned for this site
//         //                 if (alleleA > 4 || alleleB > 4) {
//         //                     // This might happen if a read is in involvedReads but doesn't support either allele at this specific site
//         //                     continue;
//         //                 }

//         //                 double weightIncrement = 0;
//         //                 if (alleleA == alleleB) {
//         //                     weightIncrement = +1.0; // Agreement weight
//         //                 } else {
//         //                     weightIncrement = -1.0; // Conflict weight
//         //                 }

//         //                 // Get indices and update the weight sum map
//         //                 size_t u = readToIndex[readA];
//         //                 size_t v = readToIndex[readB];
//         //                 if (u > v) std::swap(u, v); // Ensure consistent key order (u < v)
//         //                 edgeWeightSums[{u, v}] += weightIncrement;
//         //             }
//         //         }
//         //     }

//         //     // Add edges to the actual BGL graph from the summed weights
//         //     for(const auto& [edgePair, totalWeight] : edgeWeightSums) {
//         //         if (std::abs(totalWeight) > 0) { // Only add edges with non-zero net weight
//         //             auto edge_desc = boost::add_edge(edgePair.first, edgePair.second, phasingGraph);
//         //             // Check if edge was added successfully before setting weight
//         //             if(edge_desc.second) {
//         //                 weightMap[edge_desc.first] = totalWeight;
//         //             } else {
//         //                  // This case should ideally not happen with vecS/vecS if indices are valid
//         //                  cout << "Warning: Failed to add edge between " << edgePair.first << " and " << edgePair.second << endl;
//         //             }
//         //         }
//         //     }

//         //     // 3. Partition the Graph (Greedy BFS/Queue-based Heuristic)
//         //     std::set<OrientedReadId> visitedPhasing; // Keep track of assigned reads across components
//         //     std::queue<OrientedReadId> q;

//         //     for (const auto& startNodeRead : involvedReads) {
//         //         if (visitedPhasing.count(startNodeRead)) {
//         //             continue; // Already processed in a previous component
//         //         }

//         //         // Start phasing a new connected component
//         //         std::set<OrientedReadId> currentHaplotype1;
//         //         std::set<OrientedReadId> currentHaplotype2;

//         //         // Start the component with the current unvisited read in Haplotype 1
//         //         currentHaplotype1.insert(startNodeRead);
//         //         visitedPhasing.insert(startNodeRead);
//         //         q.push(startNodeRead);

//         //         while(!q.empty()) {
//         //             OrientedReadId currentRead = q.front();
//         //             q.pop();
//         //             size_t u = readToIndex[currentRead];
//         //             bool isInHap1 = currentHaplotype1.count(currentRead); // Check which haplotype currentRead belongs to within this component

//         //             // Iterate over neighbors in the phasing graph
//         //             PhasingGraph::adjacency_iterator neighborIt, neighborEnd;
//         //             for (boost::tie(neighborIt, neighborEnd) = boost::adjacent_vertices(u, phasingGraph); neighborIt != neighborEnd; ++neighborIt) {
//         //                 size_t v_idx = *neighborIt;
//         //                 OrientedReadId neighborRead = indexToRead[v_idx];

//         //                 if (visitedPhasing.count(neighborRead)) {
//         //                     // Optional: Check for consistency if already assigned in this component run?
//         //                     // bool neighborInHap1 = currentHaplotype1.count(neighborRead);
//         //                     // bool neighborInHap2 = currentHaplotype2.count(neighborRead);
//         //                     // if (neighborInHap1 || neighborInHap2) { ... check consistency ... }
//         //                     continue; // Already assigned globally
//         //                 }

//         //                 // Get edge weight
//         //                 auto edge_desc = boost::edge(u, v_idx, phasingGraph);
//         //                 if (!edge_desc.second) continue; // Should exist if adjacent_vertices found it
//         //                 double weight = weightMap[edge_desc.first];

//         //                 // Assign neighbor based on connection to currentRead
//         //                 if (weight > 0) { // Agreement: Assign to the same haplotype
//         //                     if (isInHap1) currentHaplotype1.insert(neighborRead);
//         //                     else currentHaplotype2.insert(neighborRead);
//         //                 } else { // Conflict (weight < 0): Assign to the opposite haplotype
//         //                     if (isInHap1) currentHaplotype2.insert(neighborRead);
//         //                     else currentHaplotype1.insert(neighborRead);
//         //                 }
//         //                 visitedPhasing.insert(neighborRead);
//         //                 q.push(neighborRead);
//         //             }
//         //         }

//         //         // Merge current component's haplotypes into final sets
//         //         // Simple merge: Add all reads. Assumes first component sets the phase.
//         //         // A more robust merge could check overlaps if final sets are non-empty.
//         //         finalHaplotype1Reads.insert(currentHaplotype1.begin(), currentHaplotype1.end());
//         //         finalHaplotype2Reads.insert(currentHaplotype2.begin(), currentHaplotype2.end());
//         //     }

//         //     // Handle reads not connected in the phasing graph (e.g., singletons or pairs with zero net weight)
//         //     // These reads could not be confidently assigned to either haplotype based on the graph structure.
//         //     for(const auto& read : involvedReads) {
//         //         if (!visitedPhasing.count(read)) {
//         //             // Add disconnected/unphased reads to a separate set instead of arbitrarily assigning
//         //             forbiddenReads[read.getReadId()] = true; // Mark as forbidden
//         //             cout << "Warning: Adding disconnected/unphased read " << read.getReadId() << " to forbiddenReads set." << endl;
//         //         }
//         //     }
//         // }

//         // // --- Global Phasing Implementation Ends ---

//         // // print the finalHaplotype1Reads and finalHaplotype2Reads
//         // cout << "Working on ReadId " << readId0 << " in strand " << strand0 << endl;
//         // cout << "Final Haplotype 1 Reads: " << endl;
//         // for (const auto& read : finalHaplotype1Reads) {
//         //     cout << "ReadId: " << read.getReadId() << " Strand: " << read.getStrand() << endl;
//         // }
//         // cout << "Final Haplotype 2 Reads: " << endl;
//         // for (const auto& read : finalHaplotype2Reads) {
//         //     cout << "ReadId: " << read.getReadId() << " Strand: " << read.getStrand() << endl;
//         // }
//         // // Check if the finalHaplotype1Reads and finalHaplotype2Reads have any common reads
//         // // There should NOT be any common reads between the two sets
//         // // Helper function to check for common elements between sets
//         // auto hasCommonElements = [](const std::set<OrientedReadId>& set1, 
//         //     const std::set<OrientedReadId>& set2) -> bool
//         // {
//         //     // Use the smaller set for iteration
//         //     const auto& smaller = (set1.size() < set2.size()) ? set1 : set2;
//         //     const auto& larger = (set1.size() < set2.size()) ? set2 : set1;

//         //     for (const auto& elem : smaller) {
//         //         if (larger.find(elem) != larger.end()) {
//         //             return true;  // Found common element
//         //         }
//         //     }
//         //     return false;  // No common elements
//         // };
//         // SHASTA_ASSERT(hasCommonElements(finalHaplotype1Reads, finalHaplotype2Reads) == false);

//         // vector< std::pair< std::set<OrientedReadId>, std::set<OrientedReadId> > > readHaplotypeSetsToKeep;
//         // if (!finalHaplotype1Reads.empty() || !finalHaplotype2Reads.empty()) {
//         //     readHaplotypeSetsToKeep.push_back(std::make_pair(finalHaplotype1Reads, finalHaplotype2Reads));
//         // }

//         // // --- Populate the sites vector for this read's haplotype block ---
//         // for (const auto& haplotypeSetPair : readHaplotypeSetsToKeep) {
//         //     // Create a new Site object for the identified haplotype pair
//         //     Site newSite;
//         //     newSite.orientedReads[0] = haplotypeSetPair.first;  // Assign haplotype 1 reads
//         //     newSite.orientedReads[1] = haplotypeSetPair.second; // Assign haplotype 2 reads

//         //     // Add the newly created site to the global sites vector
//         //     // Note: This might add duplicate or overlapping site information if multiple reads
//         //     // belong to the same underlying haplotype block. A deduplication step might
//         //     // be needed later depending on how 'sites' is used.
//         //     sites.push_back(newSite);
//         // }
//         // // --- End of populating sites vector ---













//         //return;

//         // vector< std::pair< std::set<OrientedReadId>, std::set<OrientedReadId> > > readHaplotypeSetsToKeep;
//         // std::set<OrientedReadId> unassignedOrientedReads;
//         // double numberOfTimesReadId0DoesNotHaveOneOfTheBasesInHetSites = 0;
//         // // Loop over the potential heterozygous sites starting from left to right position in readId0
//         // // that's why the potentialHetSitesOnOrientedReadId0 is modelled as a std::map
//         // for (const auto& [positionInRead0, positionStats] : potentialHetSitesOnOrientedReadId0) {
//         //     // Get the number of reads that support this position
//         //     // const uint64_t numberOfReads = positionStats.totalNumberOfAlignments;

//         //     // Get the two bases that are present at this position
//         //     const uint64_t currectHetBase1 = positionStats.hetBase1;
//         //     const uint64_t currectHetBase2 = positionStats.hetBase2;

//         //     // Get the reads that support these bases
//         //     std::set<OrientedReadId> currectHetReads1 = positionStats.hetBase1OrientedReadIds;
//         //     std::set<OrientedReadId> currectHetReads2 = positionStats.hetBase2OrientedReadIds;
//         //     // std::set<ReadId> currectHetReads1 = positionStats.hetBase1ReadIds;
//         //     // std::set<ReadId> currectHetReads2 = positionStats.hetBase2ReadIds;

//         //     // Get the alignments that support these bases
//         //     std::set<uint64_t> currectHetReads1Alignments = positionStats.hetBase1Alignments;
//         //     std::set<uint64_t> currectHetReads2Alignments = positionStats.hetBase2Alignments;


//         //     // Check if the reference read belongs to the current read haplotype set pair
//         //     if (positionStats.baseOfReadId0 == currectHetBase1) {
//         //         // The reference read belongs to the current read haplotype set pair
//         //         // We need to add the reference read to the current read haplotype set pair
//         //         currectHetReads1.insert(orientedReadId0);
//         //     } else if (positionStats.baseOfReadId0 == currectHetBase2) {
//         //         // The reference read belongs to the current read haplotype set pair
//         //         // We need to add the reference read to the current read haplotype set pair
//         //         currectHetReads2.insert(orientedReadId0);
//         //     } else {
//         //         // The reference read does not belong to the current read haplotype set pair
//         //         numberOfTimesReadId0DoesNotHaveOneOfTheBasesInHetSites++;
//         //         // TODO: Check if there are other reads supporting that base
//         //         // TODO: Those reads might belong to another similar copy of the sequence
//         //     }


            
//         //     // readHaplotypeSetsToKeep.push_back(std::make_pair(currectHetReads1, currectHetReads2));
//         //     // continue;




//         //     // If this is the first het site we meet, we need to add the current read haplotype set pair
//         //     // to the readHaplotypeSetsToKeep and move on to the next het site
//         //     if (readHaplotypeSetsToKeep.empty()) {
//         //         readHaplotypeSetsToKeep.push_back(std::make_pair(currectHetReads1, currectHetReads2));
//         //         continue;
//         //     }

//         //     // Get the last read haplotype set pair from the readHaplotypeSetsToKeep
//         //     auto lastPair = readHaplotypeSetsToKeep.back();

//         //     const auto& lastHetReads1 = lastPair.first;
//         //     const auto& lastHetReads2 = lastPair.second;

//         //     // Remove the last element
//         //     // We removed it because if we manage to find common reads between the last and current
//         //     // read haplotype sets, we will merge the sets and we will add the new merged set
//         //     readHaplotypeSetsToKeep.pop_back();

//         //     //
//         //     // Check if there are common reads between the lastHetReads and currentHetReads
//         //     //

//         //     //
//         //     // TODO: Need to check what is happening with those not common reads
//         //     // TODO: Check if they are not common because they are not lengthy enough?
//         //     // TODO: or because they do not share the same base in the adjacent het site?
//         //     //

//         //     // Case1: Common reads between lastHetReads1 and currectHetReads1
//         //     vector<OrientedReadId> commonReadsBetweenLastHetReads1AndCurrentHetReads1;
//         //     vector<OrientedReadId> notCommonReadsBetweenLastHetReads1AndCurrentHetReads1;
//         //     for (const auto& lastHetRead1 : lastHetReads1) {
//         //         if (std::find(currectHetReads1.begin(), currectHetReads1.end(), lastHetRead1) != currectHetReads1.end()) {
//         //             // We have a common read
//         //             // Add the common orientedReadId to the common reads vectors
//         //             commonReadsBetweenLastHetReads1AndCurrentHetReads1.push_back(lastHetRead1);
//         //         } else {
//         //             notCommonReadsBetweenLastHetReads1AndCurrentHetReads1.push_back(lastHetRead1);
//         //         }
//         //     }

//         //     // Case2: Common reads between lastHetReads1 and currectHetReads2
//         //     vector<OrientedReadId> commonReadsBetweenLastHetReads1AndCurrentHetReads2;
//         //     vector<OrientedReadId> notCommonReadsBetweenLastHetReads1AndCurrentHetReads2;
//         //     for (const auto& lastHetRead1 : lastHetReads1) {
//         //         if (std::find(currectHetReads2.begin(), currectHetReads2.end(), lastHetRead1) != currectHetReads2.end()) {
//         //             // We have a common read
//         //             // Add the common orientedReadId to the common reads vectors
//         //             commonReadsBetweenLastHetReads1AndCurrentHetReads2.push_back(lastHetRead1);
//         //         } else {
//         //             notCommonReadsBetweenLastHetReads1AndCurrentHetReads2.push_back(lastHetRead1);
//         //         }
//         //     }

//         //     // Case3: Common reads between lastHetReads2 and nextHetBase1
//         //     vector<OrientedReadId> commonReadsBetweenLastHetReads2AndCurrentHetReads1;
//         //     vector<OrientedReadId> notCommonReadsBetweenLastHetReads2AndCurrentHetReads1;
//         //     for (const auto& lastHetRead2 : lastHetReads2) {
//         //         if (std::find(currectHetReads1.begin(), currectHetReads1.end(), lastHetRead2) != currectHetReads1.end()) {
//         //             // We have a common read
//         //             // Add the common orientedReadId to the common reads vectors
//         //             commonReadsBetweenLastHetReads2AndCurrentHetReads1.push_back(lastHetRead2);
//         //         } else {
//         //             notCommonReadsBetweenLastHetReads2AndCurrentHetReads1.push_back(lastHetRead2);
//         //         }
//         //     }

//         //     // Case4: Common reads between lastHetReads2 and nextHetBase2
//         //     vector<OrientedReadId> commonReadsBetweenLastHetReads2AndCurrentHetReads2;
//         //     vector<OrientedReadId> notCommonReadsBetweenLastHetReads2AndCurrentHetReads2;
//         //     for (const auto& lastHetRead2 : lastHetReads2) {
//         //         if (std::find(currectHetReads2.begin(), currectHetReads2.end(), lastHetRead2) != currectHetReads2.end()) {
//         //             // We have a common read
//         //             // Add the common orientedReadId to the common reads vectors
//         //             commonReadsBetweenLastHetReads2AndCurrentHetReads2.push_back(lastHetRead2);
//         //         } else {
//         //             notCommonReadsBetweenLastHetReads2AndCurrentHetReads2.push_back(lastHetRead2);
//         //         }
//         //     }

//         //     uint64_t numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads1 = commonReadsBetweenLastHetReads1AndCurrentHetReads1.size();
//         //     uint64_t numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads2 = commonReadsBetweenLastHetReads1AndCurrentHetReads2.size();
//         //     uint64_t numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads1 = commonReadsBetweenLastHetReads2AndCurrentHetReads1.size();
//         //     uint64_t numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads2 = commonReadsBetweenLastHetReads2AndCurrentHetReads2.size();
//         //     cout << "Stats for the common reads between the last and current het sites on current position " << positionInRead0 << " of ReadId " << readId << " :" << endl;
//         //     cout << "Number of common reads between lastHetReads1 and currectHetReads1: " << numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads1 << endl;
//         //     cout << "Number of common reads between lastHetReads1 and currectHetReads2: " << numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads2 << endl;
//         //     cout << "Number of common reads between lastHetReads2 and currectHetReads1: " << numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads1 << endl;
//         //     cout << "Number of common reads between lastHetReads2 and currectHetReads2: " << numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads2 << endl;



//         //     // If we have no common reads between the 2 adjacent het sites
//         //     // we should add lastHetRead and currentHetRead sets without modifying them
//         //     // and then move to the next het site
//         //     bool noCommonReads = (numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads1 == 0 && numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads2 == 0 && numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads1 == 0 && numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads2 == 0);
//         //     if (noCommonReads) {
//         //         readHaplotypeSetsToKeep.push_back(std::make_pair(lastHetReads1, lastHetReads2));
//         //         readHaplotypeSetsToKeep.push_back(std::make_pair(currectHetReads1, currectHetReads2));
//         //         continue;
//         //     }

//         //     //
//         //     // If we reached this part it means that 
//         //     // we have common reads between the 2 adjacent het sites
//         //     // and we need to merge the sets of reads
//         //     // 

//         //     // Create a new set for the merged reads
//         //     std::set<OrientedReadId> mergedReadsHap1;
//         //     std::set<OrientedReadId> mergedReadsHap2;

//         //     // Find the combination with the most common reads
//         //     bool managedToPhaseHetSites = false;

//         //     // diploidBayesianPhase uses a Bayesian model to evaluate the phasing of two bubbles relative to each other.
//         //     uint64_t m00 = numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads1;
//         //     uint64_t m01 = numberOfCommonReadsBetweenLastHetReads1AndCurrentHetReads2;
//         //     uint64_t m10 = numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads1;
//         //     uint64_t m11 = numberOfCommonReadsBetweenLastHetReads2AndCurrentHetReads2;
//         //     const array<array<uint64_t, 2>, 2> matrix = {m00, m01, m10, m11};
//         //     double logPin, logPout;
//         //     double epsilon = 0.1;
//         //     tie(logPin, logPout) = diploidBayesianPhase(matrix, epsilon);

//         //     cout << "logPin " << logPin << " logPout " << logPout << endl;
            
//         //     double minLogP = 20;
//         //     const bool isInPhase    = (logPin - logPout) >= minLogP;
//         //     const bool isOutOfPhase = (logPout - logPin) >= minLogP;

//         //     if (isInPhase) {
//         //         // loop over commonReadsBetweenLastHetReads1AndCurrentHetReads1
//         //         for (const auto& commonRead : commonReadsBetweenLastHetReads1AndCurrentHetReads1) {
//         //             mergedReadsHap1.insert(commonRead);
//         //         }
//         //         // loop over commonReadsBetweenLastHetReads2AndCurrentHetReads2
//         //         for (const auto& commonRead : commonReadsBetweenLastHetReads2AndCurrentHetReads2) {
//         //             mergedReadsHap2.insert(commonRead);
//         //         }

//         //         // mergedReadsHap1.insert(lastHetReads1.begin(), lastHetReads1.end());
//         //         // mergedReadsHap1.insert(currectHetReads1.begin(), currectHetReads1.end());

//         //         // mergedReadsHap2.insert(lastHetReads2.begin(), lastHetReads2.end());
//         //         // mergedReadsHap2.insert(currectHetReads2.begin(), currectHetReads2.end());

//         //         if (m00 >=3 and m11 >=3) {
//         //             managedToPhaseHetSites = true;
//         //         }
                
//         //     } else if (isOutOfPhase) {
//         //         // loop over commonReadsBetweenLastHetReads1AndCurrentHetReads1
//         //         for (const auto& commonRead : commonReadsBetweenLastHetReads1AndCurrentHetReads2) {
//         //             mergedReadsHap1.insert(commonRead);
//         //         }
//         //         // loop over commonReadsBetweenLastHetReads2AndCurrentHetReads2
//         //         for (const auto& commonRead : commonReadsBetweenLastHetReads2AndCurrentHetReads1) {
//         //             mergedReadsHap2.insert(commonRead);
//         //         }

//         //         // mergedReadsHap1.insert(lastHetReads1.begin(), lastHetReads1.end());
//         //         // mergedReadsHap1.insert(currectHetReads2.begin(), currectHetReads2.end());

//         //         // mergedReadsHap2.insert(lastHetReads2.begin(), lastHetReads2.end());
//         //         // mergedReadsHap2.insert(currectHetReads1.begin(), currectHetReads1.end());

//         //         if (m01 >=3 and m10 >=3) {
//         //             managedToPhaseHetSites = true;
//         //         }
//         //     }

//         //     if (managedToPhaseHetSites) {
//         //         // Check if orientedReadIdsHetBase1 and orientedReadIdsHetBase2 have any common reads
//         //         // Helper function to check for common elements between sets
//         //         auto hasCommonElements = [](const std::set<OrientedReadId>& set1, 
//         //             const std::set<OrientedReadId>& set2) -> bool
//         //         {
//         //             // Use the smaller set for iteration
//         //             const auto& smaller = (set1.size() < set2.size()) ? set1 : set2;
//         //             const auto& larger = (set1.size() < set2.size()) ? set2 : set1;

//         //             for (const auto& elem : smaller) {
//         //                 if (larger.find(elem) != larger.end()) {
//         //                     return true;  // Found common element
//         //                 }
//         //             }
//         //             return false;  // No common elements
//         //         };

//         //         SHASTA_ASSERT(hasCommonElements(mergedReadsHap1, mergedReadsHap2) == false);
//         //     }

            


//         //     // We managed to phase the het sites
//         //     // Modify the read sets of both sites so that they have the up to date merged read sets
//         //     if (managedToPhaseHetSites) {
//         //         readHaplotypeSetsToKeep.push_back(std::make_pair(mergedReadsHap1, mergedReadsHap2));
//         //         continue;
//         //     }

//         //     // We did not manage to phase the het sites
//         //     if (!managedToPhaseHetSites) {
//         //         readHaplotypeSetsToKeep.push_back(std::make_pair(lastHetReads1, lastHetReads2));
//         //         readHaplotypeSetsToKeep.push_back(std::make_pair(currectHetReads1, currectHetReads2));
//         //         continue;
//         //     }


//         // }

//         // TO DELETE
//         // // Loop over the readHaplotypeSetsToKeep and add the reads to readIdsPseudoHap1 and readIdsPseudoHap2
//         // // readIdsPseudoHap1 and readIdsPseudoHap2 are the reads that belong to the first and second haplotype respectively
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {

//         //     const auto& haplotype1Reads = readHaplotypeSet.first;
//         //     const auto& haplotype2Reads = readHaplotypeSet.second;

//         //     // check how many reads in haplotype1Reads are already in readIdsPseudoHap1
//         //     uint64_t numberOfHaplotype1ReadsInReadIdsPseudoHap1 = 0;
//         //     for (const auto& haplotype1Read : haplotype1Reads) {
//         //         if (readIdsPseudoHap1.find(haplotype1Read.getReadId()) == readIdsPseudoHap1.end()) {
//         //             numberOfHaplotype1ReadsInReadIdsPseudoHap1++;
//         //         }
//         //     }
//         //     // check how many reads in haplotype1Reads are already in readIdsPseudoHap2
//         //     uint64_t numberOfHaplotype1ReadsInReadIdsPseudoHap2 = 0;
//         //     for (const auto& haplotype1Read : haplotype1Reads) {
//         //         if (readIdsPseudoHap2.find(haplotype1Read.getReadId()) == readIdsPseudoHap2.end()) {
//         //             numberOfHaplotype1ReadsInReadIdsPseudoHap2++;
//         //         }
//         //     }
//         //     // check how many reads in haplotype2Reads are already in readIdsPseudoHap1
//         //     uint64_t numberOfHaplotype2ReadsInReadIdsPseudoHap1 = 0;
//         //     for (const auto& haplotype2Read : haplotype2Reads) {
//         //         if (readIdsPseudoHap1.find(haplotype2Read.getReadId()) == readIdsPseudoHap1.end()) {
//         //             numberOfHaplotype2ReadsInReadIdsPseudoHap1++;
//         //         }
//         //     }
//         //     // check how many reads in haplotype2Reads are already in readIdsPseudoHap2
//         //     uint64_t numberOfHaplotype2ReadsInReadIdsPseudoHap2 = 0;
//         //     for (const auto& haplotype2Read : haplotype2Reads) {
//         //         if (readIdsPseudoHap2.find(haplotype2Read.getReadId()) == readIdsPseudoHap2.end()) {
//         //             numberOfHaplotype2ReadsInReadIdsPseudoHap2++;
//         //         }
//         //     }

//         //     if (numberOfHaplotype1ReadsInReadIdsPseudoHap1 == 0 and numberOfHaplotype1ReadsInReadIdsPseudoHap2 == 0 and numberOfHaplotype2ReadsInReadIdsPseudoHap1 == 0 and numberOfHaplotype2ReadsInReadIdsPseudoHap2 == 0) {
//         //         // haplotype1Reads and haplotype2Reads are not in readIdsPseudoHap1 or readIdsPseudoHap2
//         //         // add them to readIdsPseudoHap1 and readIdsPseudoHap2 respectively
//         //         for (const auto& haplotype1Read : haplotype1Reads) {
//         //             readIdsPseudoHap1.insert(haplotype1Read.getReadId());
//         //         }
//         //         for (const auto& haplotype2Read : haplotype2Reads) {
//         //             readIdsPseudoHap2.insert(haplotype2Read.getReadId());
//         //         }
//         //         continue;
//         //     } 

//         //     if (numberOfHaplotype1ReadsInReadIdsPseudoHap1 > numberOfHaplotype1ReadsInReadIdsPseudoHap2) {
//         //         for (const auto& haplotype1Read : haplotype1Reads) {
//         //             readIdsPseudoHap1.insert(haplotype1Read.getReadId());
//         //         }
//         //         for (const auto& haplotype2Read : haplotype2Reads) {
//         //             readIdsPseudoHap2.insert(haplotype2Read.getReadId());
//         //         }
//         //         continue;
//         //     } else if (numberOfHaplotype1ReadsInReadIdsPseudoHap1 == 0 and numberOfHaplotype1ReadsInReadIdsPseudoHap2 != 0 and numberOfHaplotype1ReadsInReadIdsPseudoHap2 > numberOfHaplotype2ReadsInReadIdsPseudoHap1) {
//         //         for (const auto& haplotype1Read : haplotype1Reads) {
//         //             readIdsPseudoHap2.insert(haplotype1Read.getReadId());
//         //         }
//         //         for (const auto& haplotype2Read : haplotype2Reads) {
//         //             readIdsPseudoHap1.insert(haplotype2Read.getReadId());
//         //         }
//         //         continue;
//         //     }

//         // }




//         // XXX
//         //
//         // We finished looping over the potential heterozygous sites
//         // and finished trying to phase all the reads in the haplotype block of readId0
//         //

//         // Find alignments between reads in different haplotype sets
//         // and mark them as forbidden
//         // This is to ensure that we do not use alignments between reads in different haplotype sets
//         // as they are not reliable


        

//         // // Iterate over the determined haplotype sets for this block
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {
//         //     const std::set<OrientedReadId>& haplotype1OrientedReadIdSet = readHaplotypeSet.first;
//         //     const std::set<OrientedReadId>& haplotype2OrientedReadIdSet = readHaplotypeSet.second;

//         //     // Helper function to extract ReadIds from a set of OrientedReadIds
//         //     auto getReadIdsFromOrientedReadIds =
//         //         [](const std::set<OrientedReadId>& orientedReadIds) -> std::set<ReadId> {
//         //         std::set<ReadId> readIds;
//         //         for (const auto& orientedReadId : orientedReadIds) {
//         //             readIds.insert(orientedReadId.getReadId());
//         //         }
//         //         return readIds;
//         //     };

//         //     // Get the ReadIds for each haplotype set
//         //     std::set<ReadId> haplotype1ReadIdSet = getReadIdsFromOrientedReadIds(haplotype1OrientedReadIdSet);
//         //     std::set<ReadId> haplotype2ReadIdSet = getReadIdsFromOrientedReadIds(haplotype2OrientedReadIdSet);

//         //     // Iterate through all alignments to find those connecting reads across haplotypes
//         //     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         //         // Skip if already marked forbidden or considered
//         //         if (forbiddenAlignments[alignmentId] || alignmentsAlreadyConsidered[alignmentId]) {
//         //             continue;
//         //         }

//         //         const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         //         // Get the readIds that the AlignmentData refers to.
//         //         ReadId thisAlignmentReadId0 = thisAlignmentData.readIds[0];
//         //         ReadId thisAlignmentReadId1 = thisAlignmentData.readIds[1];

//         //         // Skip if either read involved in the alignment is forbidden
//         //         if (forbiddenReads[thisAlignmentReadId0] || forbiddenReads[thisAlignmentReadId1]) {
//         //             // Optionally mark this alignment as forbidden too, if desired
//         //             forbiddenAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true;
//         //             continue;
//         //         }

//         //         // Check if the alignment connects a read from haplotype 1 to a read from haplotype 2
//         //         bool crossesHaplotypes =
//         //             (haplotype1ReadIdSet.count(thisAlignmentReadId0) && haplotype2ReadIdSet.count(thisAlignmentReadId1)) ||
//         //             (haplotype1ReadIdSet.count(thisAlignmentReadId1) && haplotype2ReadIdSet.count(thisAlignmentReadId0));

//         //         // Check if the alignment connects two reads within the same haplotype
//         //         bool withinHaplotype1 = haplotype1ReadIdSet.count(thisAlignmentReadId0) && haplotype1ReadIdSet.count(thisAlignmentReadId1);
//         //         bool withinHaplotype2 = haplotype2ReadIdSet.count(thisAlignmentReadId0) && haplotype2ReadIdSet.count(thisAlignmentReadId1);
//         //         bool withinSameHaplotype = withinHaplotype1 || withinHaplotype2;


//         //         if (crossesHaplotypes) {
//         //             // Mark this alignment as forbidden
//         //             cout << "Forbidding alignment " << alignmentId << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << " (cross-haplotype)" << endl;
//         //             forbiddenAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true; // Mark as considered to avoid redundant checks
//         //         } else if (withinSameHaplotype) {
//         //             // Mark this alignment as a first pass het alignment (intra-haplotype)
//         //             // cout << "Marking alignment " << alignmentId << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << " as first pass (intra-haplotype)" << endl;
//         //             firstPassHetAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true; // Mark as considered
//         //         }

//         //         // if (
//         //         //     (haplotype1ReadIdSet.count(thisAlignmentReadId0) && !haplotype1ReadIdSet.count(thisAlignmentReadId1)) ||
//         //         //     (haplotype1ReadIdSet.count(thisAlignmentReadId1) && !haplotype1ReadIdSet.count(thisAlignmentReadId0))
//         //         // ) {
//         //         //     firstPassHetAlignments[alignmentId] = true;
//         //         //     alignmentsAlreadyConsidered[alignmentId] = true;
//         //         //     haplotype1ReadIdSet.insert(thisAlignmentReadId0);
//         //         //     haplotype1ReadIdSet.insert(thisAlignmentReadId1);
//         //         //     continue;
//         //         // }

//         //         // if (
//         //         //     (haplotype2ReadIdSet.count(thisAlignmentReadId0) && !haplotype2ReadIdSet.count(thisAlignmentReadId1)) ||
//         //         //     (haplotype2ReadIdSet.count(thisAlignmentReadId1) && !haplotype2ReadIdSet.count(thisAlignmentReadId0))
//         //         // ) {
//         //         //     firstPassHetAlignments[alignmentId] = true;
//         //         //     alignmentsAlreadyConsidered[alignmentId] = true;
//         //         //     haplotype2ReadIdSet.insert(thisAlignmentReadId0);
//         //         //     haplotype2ReadIdSet.insert(thisAlignmentReadId1);
//         //         //     continue;
//         //         // }



//         //     }
//         // }



        
//         // // Iterate over the determined haplotype sets for this block
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {
//         //     const std::set<OrientedReadId>& haplotype1OrientedReadIdSet = readHaplotypeSet.first;
//         //     const std::set<OrientedReadId>& haplotype2OrientedReadIdSet = readHaplotypeSet.second;

//         //     // Helper function to extract ReadIds from a set of OrientedReadIds
//         //     auto getReadIdsFromOrientedReadIds =
//         //         [](const std::set<OrientedReadId>& orientedReadIds) -> std::set<ReadId> {
//         //         std::set<ReadId> readIds;
//         //         for (const auto& orientedReadId : orientedReadIds) {
//         //             readIds.insert(orientedReadId.getReadId());
//         //         }
//         //         return readIds;
//         //     };

//         //     // Get the ReadIds for each haplotype set
//         //     std::set<ReadId> haplotype1ReadIdSet = getReadIdsFromOrientedReadIds(haplotype1OrientedReadIdSet);
//         //     std::set<ReadId> haplotype2ReadIdSet = getReadIdsFromOrientedReadIds(haplotype2OrientedReadIdSet);

//         //     std::set<ReadId> connectsToHap1;
//         //     std::set<ReadId> connectsToHap2;

//         //     std::set<ReadId> addedToHap1;
//         //     std::set<ReadId> addedToHap2;

//         //     for (uint64_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
//         //         // Optional: Skip forbidden alignments
//         //         if (!forbiddenAlignments.empty() && forbiddenAlignments[alignmentId]) {
//         //              continue;
//         //         }

//         //         if (alignmentsAlreadyConsidered[alignmentId]) {
//         //             continue;
//         //         }
        
//         //         const AlignmentData& aln = alignmentData[alignmentId];
//         //         ReadId readA = aln.readIds[0];
//         //         ReadId readB = aln.readIds[1];
        
//         //         // Optional: Skip if reads themselves are forbidden
//         //         if (!forbiddenReads.empty() && (forbiddenReads[readA] || forbiddenReads[readB])) {
//         //             continue;
//         //         }
        
//         //         bool aInHap1 = haplotype1ReadIdSet.count(readA);
//         //         bool bInHap1 = haplotype1ReadIdSet.count(readB);
//         //         bool aInHap2 = haplotype2ReadIdSet.count(readA);
//         //         bool bInHap2 = haplotype2ReadIdSet.count(readB);

//         //         if(aInHap1 && !bInHap1 && !bInHap2 && !addedToHap1.count(readB) && !addedToHap2.count(readB)) {
//         //             // firstPassHetAlignments[alignmentId] = true;
//         //             // alignmentsAlreadyConsidered[alignmentId] = true;
//         //             addedToHap1.insert(readB);
//         //         } 
//         //         else if(bInHap1 && !aInHap1 && !aInHap2 && !addedToHap1.count(readA) && !addedToHap2.count(readA)) {
//         //             // firstPassHetAlignments[alignmentId] = true;
//         //             // alignmentsAlreadyConsidered[alignmentId] = true;
//         //             addedToHap1.insert(readA);
//         //         }
//         //         else if(aInHap2 && !bInHap1 && !bInHap2 && !addedToHap1.count(readB) && !addedToHap2.count(readB)) {
//         //             // firstPassHetAlignments[alignmentId] = true;
//         //             // alignmentsAlreadyConsidered[alignmentId] = true;
//         //             addedToHap2.insert(readB);
//         //         } 
//         //         else if(bInHap2 && !aInHap1 && !aInHap2 && !addedToHap1.count(readA) && !addedToHap2.count(readA)) {
//         //             // firstPassHetAlignments[alignmentId] = true;
//         //             // alignmentsAlreadyConsidered[alignmentId] = true;
//         //             addedToHap2.insert(readA);
//         //         }
        
//         //         // Record which reads connect to which haplotype set via this alignment
//         //         if (aInHap1) connectsToHap1.insert(readB);
//         //         if (bInHap1) connectsToHap1.insert(readA);
//         //         if (aInHap2) connectsToHap2.insert(readB);
//         //         if (bInHap2) connectsToHap2.insert(readA);
//         //     }

//         //     // Find the intersection: reads connected to BOTH Hap1 and Hap2
//         //     std::set_intersection(
//         //         connectsToHap1.begin(), connectsToHap1.end(),
//         //         connectsToHap2.begin(), connectsToHap2.end(),
//         //         std::inserter(bridgingReads, bridgingReads.begin())
//         //     );


//         //     // --- New loop to forbid alignments based on addedToHap1/addedToHap2 ---
//         //     // Iterate through all alignments again to enforce new forbidden connections
//         //     for (uint64_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
//         //         // Skip if already forbidden
//         //         if (!forbiddenAlignments.empty() && forbiddenAlignments[alignmentId]) {
//         //             continue;
//         //         }

//         //         if (alignmentsAlreadyConsidered[alignmentId]) {
//         //             continue;
//         //         }

//         //         const AlignmentData& aln = alignmentData[alignmentId];
//         //         ReadId readA = aln.readIds[0];
//         //         ReadId readB = aln.readIds[1];

//         //         // Check if alignment connects a read tentatively added to Hap1 with a read in original Hap2
//         //         bool added1_vs_orig2 = (addedToHap1.count(readA) && haplotype2ReadIdSet.count(readB)) ||
//         //                                (addedToHap1.count(readB) && haplotype2ReadIdSet.count(readA));

//         //         // Check if alignment connects a read tentatively added to Hap2 with a read in original Hap1
//         //         bool added2_vs_orig1 = (addedToHap2.count(readA) && haplotype1ReadIdSet.count(readB)) ||
//         //                                (addedToHap2.count(readB) && haplotype1ReadIdSet.count(readA));

//         //         // Check if alignment connects a read tentatively added to Hap1 with a read in original Hap1
//         //         bool added1_vs_orig1 = (addedToHap1.count(readA) && haplotype1ReadIdSet.count(readB)) ||
//         //                                (addedToHap1.count(readB) && haplotype1ReadIdSet.count(readA));

//         //         // Check if alignment connects a read tentatively added to Hap2 with a read in original Hap2
//         //         bool added2_vs_orig2 = (addedToHap2.count(readA) && haplotype2ReadIdSet.count(readB)) ||
//         //                                (addedToHap2.count(readB) && haplotype2ReadIdSet.count(readA));

//         //         if (added1_vs_orig2) {
//         //             // This alignment connects a read tentatively assigned to Hap1 extension
//         //             // with a read from the original Hap2 set. Forbid it.
//         //             cout << "Forbidding alignment " << alignmentId << " involving reads: " << readA << " and " << readB << " (addedToHap1 vs haplotype2ReadIdSet)" << endl;
//         //             forbiddenAlignments[alignmentId] = true;
//         //             // If this alignment was previously marked as firstPassHet, unmark it.
//         //             if (firstPassHetAlignments[alignmentId]) {
//         //                 firstPassHetAlignments[alignmentId] = false;
//         //             }
//         //         }

//         //         if (added2_vs_orig1) {
//         //             // This alignment connects a read tentatively assigned to Hap2 extension
//         //             // with a read from the original Hap1 set. Forbid it.
//         //             cout << "Forbidding alignment " << alignmentId << " involving reads: " << readA << " and " << readB << " (addedToHap2 vs haplotype1ReadIdSet)" << endl;
//         //             forbiddenAlignments[alignmentId] = true;
//         //              // If this alignment was previously marked as firstPassHet, unmark it.
//         //             if (firstPassHetAlignments[alignmentId]) {
//         //                 firstPassHetAlignments[alignmentId] = false;
//         //             }
//         //         }

//         //         // if (added1_vs_orig1) {
//         //         //     // This alignment connects a read tentatively assigned to Hap1 extension
//         //         //     // with a read from the original Hap1 set. Allow it.
//         //         //     cout << "Adding alignment " << alignmentId << " involving reads: " << readA << " and " << readB << " (addedToHap1 vs haplotype1ReadIdSet)" << endl;
//         //         //     firstPassHetAlignments[alignmentId] = true;
//         //         //     alignmentsAlreadyConsidered[alignmentId] = true;
//         //         // }

//         //         // if (added2_vs_orig2) {
//         //         //     // This alignment connects a read tentatively assigned to Hap2 extension
//         //         //     // with a read from the original Hap2 set. Allow it.
//         //         //     cout << "Adding alignment " << alignmentId << " involving reads: " << readA << " and " << readB << " (addedToHap2 vs haplotype2ReadIdSet)" << endl;
//         //         //     firstPassHetAlignments[alignmentId] = true;
//         //         //     alignmentsAlreadyConsidered[alignmentId] = true;
//         //         // }
//         //     }
//         //     // --- End of new loop ---
//         // }
















//         // // Gather together all the reads that managed to get phased in the haplotype block
//         // std::set<OrientedReadId> orientedReadIdsThatManagedToGetPhasedInTheHaplotypeBlock;
//         // std::set<ReadId> readIdsThatManagedToGetPhasedInTheHaplotypeBlock;
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {
//         //         const std::set<OrientedReadId>& haplotype1OrientedReadIdSet = readHaplotypeSet.first;
//         //         const std::set<OrientedReadId>& haplotype2OrientedReadIdSet = readHaplotypeSet.second;
//         //         for (const auto& haplotype1OrientedReadId : haplotype1OrientedReadIdSet) {
//         //             readIdsThatManagedToGetPhasedInTheHaplotypeBlock.insert(haplotype1OrientedReadId.getReadId());
//         //         }
//         //         for (const auto& haplotype2OrientedReadId : haplotype2OrientedReadIdSet) {
//         //             readIdsThatManagedToGetPhasedInTheHaplotypeBlock.insert(haplotype2OrientedReadId.getReadId());
//         //         }
//         //         orientedReadIdsThatManagedToGetPhasedInTheHaplotypeBlock.insert(haplotype1OrientedReadIdSet.begin(), haplotype1OrientedReadIdSet.end());
//         //         orientedReadIdsThatManagedToGetPhasedInTheHaplotypeBlock.insert(haplotype2OrientedReadIdSet.begin(), haplotype2OrientedReadIdSet.end());
//         // }

//         // // find the alignments in readHaplotypeSetsToKeep and forbid them if they are between reads in opposite sets.
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {
//         //     const std::set<OrientedReadId>& haplotype1OrientedReadIdSet = readHaplotypeSet.first;
//         //     const std::set<OrientedReadId>& haplotype2OrientedReadIdSet = readHaplotypeSet.second;
//         //     // Helper function to extract ReadIds from a set of OrientedReadIds
//         //     auto getReadIdsFromOrientedReadIds =
//         //         [](const std::set<OrientedReadId>& orientedReadIds) -> std::set<ReadId> {
//         //         std::set<ReadId> readIds;
//         //         for (const auto& orientedReadId : orientedReadIds) {
//         //             readIds.insert(orientedReadId.getReadId());
//         //         }
//         //         return readIds;
//         //     };
//         //     const std::set<ReadId> haplotype1ReadIdSet = getReadIdsFromOrientedReadIds(haplotype1OrientedReadIdSet);
//         //     const std::set<ReadId> haplotype2ReadIdSet = getReadIdsFromOrientedReadIds(haplotype2OrientedReadIdSet);
//         //     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         //         const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         //         // Get the readIds that the AlignmentData refers to.
//         //         ReadId thisAlignmentReadId0 = thisAlignmentData.readIds[0];
//         //         ReadId thisAlignmentReadId1 = thisAlignmentData.readIds[1];

//         //         // Check if the readIds are in the readIdsThatManagedToGetPhasedInTheHaplotypeBlock
//         //         if (haplotype1ReadIdSet.contains(thisAlignmentReadId0) and haplotype2ReadIdSet.contains(thisAlignmentReadId1)) {
//         //             cout << "Removing alignment " << alignmentId << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << endl;
//         //             forbiddenAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true;
//         //             continue;
//         //         }

//         //         if (haplotype1ReadIdSet.contains(thisAlignmentReadId1) and haplotype2ReadIdSet.contains(thisAlignmentReadId0)) {
//         //             cout << "Removing alignment " << alignmentId << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << endl;
//         //             forbiddenAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true;
//         //             continue;
//         //         }

//         //     }
//         // }














        

//         // // Now we need to check if the reference read is involved in the het sites
//         // // If numberOfTimesReadId0DoesNotHaveOneOfTheBasesInHetSites is greater than 0
//         // // then the ref read should not align with the reads involved in the het heplotypes.
//         // // WE TREAT THE HET SITES AS A STRONG SIGNAL FOR WHETHER THE ALIGNMENT IS GOOD OR NOT.
//         // bool isRefReadInvolvedInHetSites = false;

//         // if (potentialHetSitesOnOrientedReadId0.size() > 3) {
//         //     if(numberOfTimesReadId0DoesNotHaveOneOfTheBasesInHetSites / double(potentialHetSitesOnOrientedReadId0.size()) <= 0.25) {
//         //         isRefReadInvolvedInHetSites = true;
//         //     } else {
//         //         isRefReadInvolvedInHetSites = false;
//         //     }
//         // } else if (potentialHetSitesOnOrientedReadId0.size() <= 3) {
//         //     if(numberOfTimesReadId0DoesNotHaveOneOfTheBasesInHetSites == 0) {
//         //         isRefReadInvolvedInHetSites = true;
//         //     } else {
//         //         isRefReadInvolvedInHetSites = false;
//         //     }
//         // }




//         // // If the ref read is NOT involved in the het sites
//         // // Forbid the alignments of ref read with the reads in the het sites
//         // // TODO: Also try to do the same for the other reads that do not belong with the het reads
//         // if (not isRefReadInvolvedInHetSites) {
//         //     cout << "The reference read does not have the bases of the het sites." << endl;

//         //     for (const auto& [positionInRead0, positionStats] : potentialHetSitesOnOrientedReadId0) {
    
//         //         // Get the alignments that support these bases
//         //         std::set<uint64_t> currectHetReads1Alignments = positionStats.hetBase1Alignments;
//         //         std::set<uint64_t> currectHetReads2Alignments = positionStats.hetBase2Alignments;

//         //         for (const auto& currectHetReads1Alignment : currectHetReads1Alignments) {

//         //             const AlignmentData& thisAlignmentData = alignmentData[currectHetReads1Alignment];
//         //             // Get the readIds that the AlignmentData refers to.
//         //             ReadId thisAlignmentReadId0 = thisAlignmentData.readIds[0];
//         //             ReadId thisAlignmentReadId1 = thisAlignmentData.readIds[1];

//         //             if (readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId0) and readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId1)) {
//         //                 if (alignmentsAlreadyConsidered[currectHetReads1Alignment]) {
//         //                     continue;
//         //                 }
//         //                 cout << "Removing alignment " << currectHetReads1Alignment << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << endl;
//         //                 forbiddenAlignments[currectHetReads1Alignment] = true;
//         //                 alignmentsAlreadyConsidered[currectHetReads1Alignment] = true;
//         //             }
                    
//         //         }
//         //         for (const auto& currectHetReads2Alignment : currectHetReads2Alignments) {
                    
//         //             const AlignmentData& thisAlignmentData = alignmentData[currectHetReads2Alignment];
//         //             // Get the readIds that the AlignmentData refers to.
//         //             ReadId thisAlignmentReadId0 = thisAlignmentData.readIds[0];
//         //             ReadId thisAlignmentReadId1 = thisAlignmentData.readIds[1];

//         //             if (readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId0) and readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId1)) {
//         //                 if (alignmentsAlreadyConsidered[currectHetReads2Alignment]) {
//         //                     continue;
//         //                 }
//         //                 cout << "Removing alignment " << currectHetReads2Alignment << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << endl;
//         //                 forbiddenAlignments[currectHetReads2Alignment] = true;
//         //                 alignmentsAlreadyConsidered[currectHetReads2Alignment] = true;
//         //             }
                    
//         //         }
//         //     }
//         // }

//         // // If the ref read IS involved in the het sites
//         // // Use the alignments of ref read with the reads in the het sites
//         // if (isRefReadInvolvedInHetSites) {
//         //     cout << "The reference read have the bases of the het sites." << endl;
//         //     for (const auto& [positionInRead0, positionStats] : potentialHetSitesOnOrientedReadId0) {
    
//         //         // Get the alignments that support these bases
//         //         std::set<uint64_t> currectHetReads1Alignments = positionStats.hetBase1Alignments;
//         //         std::set<uint64_t> currectHetReads2Alignments = positionStats.hetBase2Alignments;

//         //         for (const auto& currectHetReads1Alignment : currectHetReads1Alignments) {

//         //             const AlignmentData& thisAlignmentData = alignmentData[currectHetReads1Alignment];
//         //             // Get the readIds that the AlignmentData refers to.
//         //             ReadId thisAlignmentReadId0 = thisAlignmentData.readIds[0];
//         //             ReadId thisAlignmentReadId1 = thisAlignmentData.readIds[1];
//         //             if (readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId0) and readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId1)) {
//         //                 if (alignmentsAlreadyConsidered[currectHetReads1Alignment]) {
//         //                     continue;
//         //                 }
//         //                 // cout << "Using alignment " << currectHetReads1Alignment << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << endl;
//         //                 firstPassHetAlignments[currectHetReads1Alignment] = true;
//         //                 alignmentsAlreadyConsidered[currectHetReads1Alignment] = true;
//         //             }
                    
//         //         }
//         //         for (const auto& currectHetReads2Alignment : currectHetReads2Alignments) {

//         //             const AlignmentData& thisAlignmentData = alignmentData[currectHetReads2Alignment];
//         //             // Get the readIds that the AlignmentData refers to.
//         //             ReadId thisAlignmentReadId0 = thisAlignmentData.readIds[0];
//         //             ReadId thisAlignmentReadId1 = thisAlignmentData.readIds[1];
//         //             if (readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId0) and readIdsThatManagedToGetPhasedInTheHaplotypeBlock.contains(thisAlignmentReadId1)) {
//         //                 if (alignmentsAlreadyConsidered[currectHetReads2Alignment]) {
//         //                     continue;
//         //                 }
//         //                 // cout << "Using alignment " << currectHetReads2Alignment << " involving reads: " << thisAlignmentReadId0 << " and " << thisAlignmentReadId1 << endl;
//         //                 firstPassHetAlignments[currectHetReads2Alignment] = true;
//         //                 alignmentsAlreadyConsidered[currectHetReads2Alignment] = true;
//         //             }
                    
//         //         }
//         //     }
//         // }


//         // cout << "Found " << readHaplotypeSetsToKeep.size() << " haplotype sets in readId " << readId0 << endl;
//         // // Print the readHaplotypeSetsToKeep
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {
//         //     cout << "Haplotype set: " << endl;
//         //     cout << "Haplotype 1: " << endl;
//         //     for (const auto& orientedReadId : readHaplotypeSet.first) {
//         //         cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//         //     }
//         //     cout << "Haplotype 2: " << endl;
//         //     for (const auto& orientedReadId : readHaplotypeSet.second) {
//         //         cout << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
//         //     }
//         // }


//         // vector< std::pair< std::set<ReadId>, std::set<ReadId> > > readHaplotypeSetsToKeepReadIdsOnly;
//         // for (const auto& readHaplotypeSet : readHaplotypeSetsToKeep) {
//         //     const std::set<OrientedReadId>& haplotype1OrientedReadIdSet = readHaplotypeSet.first;
//         //     const std::set<OrientedReadId>& haplotype2OrientedReadIdSet = readHaplotypeSet.second;
//         //     std::set<ReadId> haplotype1ReadIdSet;
//         //     std::set<ReadId> haplotype2ReadIdSet;
//         //     for (const auto& orientedReadId : haplotype1OrientedReadIdSet) {
//         //         haplotype1ReadIdSet.insert(orientedReadId.getReadId());
//         //     }
//         //     for (const auto& orientedReadId : haplotype2OrientedReadIdSet) {
//         //         haplotype2ReadIdSet.insert(orientedReadId.getReadId());
//         //     }
//         //     readHaplotypeSetsToKeepReadIdsOnly.push_back(std::make_pair(haplotype1ReadIdSet, haplotype2ReadIdSet));
            
//         // }



//         // //
//         // //
//         // // Finally, we need to forbid the alignments between the reads in the haplotype sets
//         // // that belong to different haplotypes.
//         // // We have the alignmentIds between the reference read and the reads in the haplotype sets.
//         // // We don't have the alignmentIds between the reads in the haplotype sets.
//         // // So, we need to find them first.
//         // //
//         // //

//         // for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {

//         //     if (forbiddenAlignments[alignmentId]) {
//         //         continue;
//         //     }

//         //     // Get information for this alignment.
//         //     AlignmentData& thisAlignmentData = alignmentData[alignmentId];

//         //     // The alignment is stored as an alignment between readId0 on strand 0
//         //     // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         //     // The reverse complement alignment also exists, but is not stored explicitly.
            
//         //     const ReadId readId0 = thisAlignmentData.readIds[0];
//         //     const ReadId readId1 = thisAlignmentData.readIds[1];
//         //     const bool isSameStrand = thisAlignmentData.isSameStrand;
//         //     SHASTA_ASSERT(readId0 < readId1);
//         //     const OrientedReadId A0(readId0, 0);
//         //     const OrientedReadId B0(readId1, isSameStrand ? 0 : 1);
//         //     const OrientedReadId A1(readId0, 1);
//         //     const OrientedReadId B1(readId1, isSameStrand ? 1 : 0);

//         //     // Check if the alignment involves a read in the haplotype sets
//         //     for (const auto& readHaplotypeSet : readHaplotypeSetsToKeepReadIdsOnly) {
//         //         const std::set<ReadId>& haplotype1ReadIdSet = readHaplotypeSet.first;
//         //         const std::set<ReadId>& haplotype2ReadIdSet = readHaplotypeSet.second;
                
//         //         // We have an alignment that involves two reads that belong to different haplotypes
//         //         // We need to forbid the alignment
//         //         if ((haplotype1ReadIdSet.contains(readId0) && haplotype2ReadIdSet.contains(readId1)) 
//         //             || (haplotype1ReadIdSet.contains(readId1) && haplotype2ReadIdSet.contains(readId0))) {
//         //             cout << "Forbidding alignment " << alignmentId << " involving reads: " << readId0 << " and " << readId1 << " because they belong to different haplotypes." << endl;
//         //             forbiddenAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true;
//         //         }

//         //         // We have an alignment that involves two reads that belong to the same haplotype
//         //         // We need to use the alignment
//         //         if ((haplotype1ReadIdSet.contains(readId0) && haplotype1ReadIdSet.contains(readId1)) 
//         //             || (haplotype2ReadIdSet.contains(readId0) && haplotype2ReadIdSet.contains(readId1))) {
//         //             // cout << "Using alignment " << alignmentId << " involving reads: " << readId0 << " and " << readId1 << " because they belong to the same haplotype." << endl;
//         //             firstPassHetAlignments[alignmentId] = true;
//         //             alignmentsAlreadyConsidered[alignmentId] = true;
//         //         }

//         //     }

//         // }

//         // const long firstPassHetAlignmentsCount = count(firstPassHetAlignments.begin(), firstPassHetAlignments.end(), true);
//         // cout << "Found " << firstPassHetAlignmentsCount << " first pass alignments in readId " << readId0 << endl;


//     }

//     // print the current timestamp
//     cout << timestamp << "Finished first pass of phasing." << endl;

//     // return;

//     // cout << timestamp << "Finished first pass of phasing." << endl;

//     // // --- We need to find, for each orientedReadId, the sites that are associated with it ---
//     // // --- and fill in the orientedReadSites ---
//     // vector< vector<uint64_t> > orientedReadSites(readCount * 2);    // Indexed via OrientedReadId::getValue()

//     // for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
//     //     const Site& site = sites[siteId];
//     //     for(const std::set<OrientedReadId>& s: site.orientedReads) {
//     //         for(const OrientedReadId orientedReadId: s) {
//     //             orientedReadSites[orientedReadId.getValue()].push_back(siteId);
//     //         }
//     //     }
//     // }

//     // cout << timestamp << "Finished second pass of phasing." << endl;


//     // // --- We need to generate all possible pairs of sites from the orientedReadSites ---
//     // // --- The siteId of the first site in each pair will always be less than the siteId of the second site in each pair ---
//     // vector< pair<uint64_t, uint64_t> > pairsOfSites;
//     // for(const vector<uint64_t>& v: orientedReadSites) {
//     //     // Add a check to prevent accessing v.size()-1 when v is empty or has only one element.
//     //     if (v.size() < 2) {
//     //         continue; // Skip if there are not enough elements to form a pair
//     //     }
//     //     for(uint64_t i0=0; i0<v.size()-1; i0++) {
//     //         const uint64_t siteId0 = v[i0];
//     //         for(uint64_t i1=i0+1; i1<v.size(); i1++) {
//     //             const uint64_t siteId1 = v[i1];
//     //             pairsOfSites.push_back(make_pair(siteId0, siteId1));
//     //         }
//     //     }
//     // }

//     // cout << timestamp << "Finished third pass of phasing." << endl;

//     // // --- Because the siteId of the first site in each pair will always be less than the siteId of 
//     // // the second site in each pair, we will have to deduplicate the pairs ---
//     // // --- Deduplicate pairs and count common reads ---
//     // vector<uint64_t> commonReadsCount; // Stores the count of common reads for each unique pair
//     // deduplicateAndCount(pairsOfSites, commonReadsCount);
//     // SHASTA_ASSERT(commonReadsCount.size() == pairsOfSites.size());

//     // cout << timestamp << "Finished fourth pass of phasing." << endl;

//     // // --- Set up the disjoint sets data structure. ---
//     // // --- Each vertex is a heterozygous phased site ---
//     // const uint64_t vertexCount = sites.size();
//     // vector<uint64_t> rankHetSites(vertexCount);
//     // vector<uint64_t> parentHetSites(vertexCount);
//     // boost::disjoint_sets<uint64_t*, uint64_t*> disjointSetsHetSites(&rankHetSites[0], &parentHetSites[0]);
//     // for(uint64_t i=0; i<vertexCount; i++) {
//     //     disjointSetsHetSites.make_set(i);
//     // }

//     // cout << timestamp << "Finished fifth pass of phasing." << endl;

//     // // --- Loop over all pairs of sites that share at least m OrientedReadIds. ---
//     // const uint64_t minCommonReadsForMerging = 2;
//     // uint64_t mergeCount = 0;
//     // uint64_t pairsConsidered = 0;
//     // for(uint64_t i=0; i<pairsOfSites.size(); i++) {
//     //     pairsConsidered++;
//     //     const pair<uint64_t, uint64_t>& p = pairsOfSites[i];
//     //     const uint64_t commonOrientedReadCount = commonReadsCount[i];
//     //     const uint64_t siteId0 = p.first;
//     //     const uint64_t siteId1 = p.second;
//     //     const Site& site0 = sites[siteId0];
//     //     const Site& site1 = sites[siteId1];
        
//     //     // Decide if these two should be merged based on the number of common reads
//     //     const bool merge = (commonOrientedReadCount >= minCommonReadsForMerging);
        
//     //     if(merge) {
//     //         // Merge the sets (connectedComponents) containing siteId0 and siteId1
//     //         disjointSetsHetSites.union_set(siteId0, siteId1);
//     //         mergeCount++;
//     //     }
//     // }
//     // cout << timestamp << "Finished site merging. Considered " << pairsConsidered << " pairs, merged " << mergeCount << " pairs based on common read count >= " << minCommonReadsForMerging << "." << endl;


//     // // --- At this point, disjointSetsHetSites.find_set(siteId) gives the id of the merged set
//     // // --- that siteId is part of. This id is in [0, vertexCount).
//     // vector< vector<uint64_t> > connectedComponents;
//     // for(uint64_t siteId=0; siteId<vertexCount; siteId++) {
//     //     const uint64_t componentId = disjointSetsHetSites.find_set(siteId);
//     //     connectedComponents[componentId].push_back(siteId);
//     // }

//     // cout << timestamp << "Finished sixth pass of phasing." << endl;


//     // // --- Global Phasing Implementation Starts Here ---
//     // cout << timestamp << "Starting global phasing of connectedComponents." << endl;
//     // uint64_t totalPhasedReads = 0;
//     // uint64_t totalComponentsPhased = 0;

//     // // Loop over the connected connectedComponents.
//     // // Each component is a set of Sites that were recursively merged together.
//     // // We only process compoconnectedComponentsnents that actually contain sites.
//     // for(uint64_t componentId=0; componentId<connectedComponents.size(); componentId++) {
//     //     const vector<uint64_t>& componentSiteIds = connectedComponents[componentId];

//     //     cout << timestamp << "Finished seventh pass of phasing." << endl;

//     //     // Skip empty connectedComponents (these represent siteIds that were not part of any merged site pair)
//     //     if (componentSiteIds.empty()) {
//     //         continue;
//     //     }

//     //     // Skip connectedComponents with only one site, as phasing requires at least two linked sites.
//     //     if (componentSiteIds.size() < 2) {
//     //          cout << "Skipping component " << componentId << " with only " << componentSiteIds.size() << " site(s)." << endl;
//     //          continue;
//     //     }

//     //     cout << "Phasing component " << componentId << " with " << componentSiteIds.size() << " sites." << endl;
//     //     totalComponentsPhased++;

//     //     // 1. Identify Involved Reads within this component
//     //     std::set<OrientedReadId> involvedReads;
//     //     for(const uint64_t siteId : componentSiteIds) {
//     //         const Site& site = sites[siteId];
//     //         involvedReads.insert(site.orientedReads[0].begin(), site.orientedReads[0].end());
//     //         involvedReads.insert(site.orientedReads[1].begin(), site.orientedReads[1].end());
//     //     }

//     //     if (involvedReads.empty()) {
//     //         cout << "Component " << componentId << " has no involved reads. Skipping." << endl;
//     //         continue;
//     //     }
//     //     cout << "Component " << componentId << " involves " << involvedReads.size() << " reads." << endl;


//     //     // Create mapping for graph vertices
//     //     std::map<OrientedReadId, size_t> readToIndex;
//     //     std::vector<OrientedReadId> indexToRead;
//     //     size_t currentIndex = 0;
//     //     for(const auto& read : involvedReads) {
//     //         readToIndex[read] = currentIndex;
//     //         indexToRead.push_back(read);
//     //         currentIndex++;
//     //     }
//     //     size_t numInvolvedReads = involvedReads.size();

//     //     // 2. Build Phasing Graph for the component
//     //     using PhasingGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_weight_t, double>>;
//     //     PhasingGraph phasingGraph(numInvolvedReads);
//     //     auto weightMap = boost::get(boost::edge_weight, phasingGraph);
//     //     std::map<std::pair<size_t, size_t>, double> edgeWeightSums;

//     //     // Iterate over all pairs of involved reads
//     //     for (auto itA = involvedReads.begin(); itA != involvedReads.end(); ++itA) {
//     //         auto itB = itA;
//     //         std::advance(itB, 1);
//     //         for (; itB != involvedReads.end(); ++itB) {
//     //             const OrientedReadId& readA = *itA;
//     //             const OrientedReadId& readB = *itB;
//     //             double currentWeightSum = 0;

//     //             // Check all sites in the component if they contain both readA and readB
//     //             for (const uint64_t siteId : componentSiteIds) {
//     //                  const Site& site = sites[siteId];
//     //                  bool aInHap0 = site.orientedReads[0].count(readA);
//     //                  bool aInHap1 = site.orientedReads[1].count(readA);
//     //                  bool bInHap0 = site.orientedReads[0].count(readB);
//     //                  bool bInHap1 = site.orientedReads[1].count(readB);

//     //                  // Only contribute weight if both reads are present in this site
//     //                  if ((aInHap0 || aInHap1) && (bInHap0 || bInHap1)) {
//     //                      if ((aInHap0 && bInHap0) || (aInHap1 && bInHap1)) {
//     //                          currentWeightSum += 1.0; // Agreement
//     //                      } else {
//     //                          currentWeightSum -= 1.0; // Conflict
//     //                      }
//     //                  }
//     //             }

//     //             // If there's a non-zero net weight, record it
//     //             if (std::abs(currentWeightSum) > 0) {
//     //                 size_t u = readToIndex[readA];
//     //                 size_t v = readToIndex[readB];
//     //                 if (u > v) std::swap(u, v);
//     //                 edgeWeightSums[{u, v}] = currentWeightSum; // Use direct assignment, assuming weight comes only from sites within this component
//     //             }
//     //         }
//     //     }

//     //     // Add edges to the BGL graph
//     //     for(const auto& [edgePair, totalWeight] : edgeWeightSums) {
//     //         auto edge_desc = boost::add_edge(edgePair.first, edgePair.second, phasingGraph);
//     //         if(edge_desc.second) {
//     //             weightMap[edge_desc.first] = totalWeight;
//     //         } else {
//     //              cout << "Warning: Failed to add edge between " << edgePair.first << " and " << edgePair.second << " for component " << componentId << endl;
//     //         }
//     //     }

//     //     // 3. Partition the Graph (Greedy BFS/Queue-based Heuristic)
//     //     std::set<OrientedReadId> componentHaplotype1;
//     //     std::set<OrientedReadId> componentHaplotype2;
//     //     std::set<OrientedReadId> componentUnphased; // Reads that couldn't be confidently assigned
//     //     std::set<OrientedReadId> visitedInComponent;
//     //     std::queue<OrientedReadId> q;

//     //     for (const auto& startNodeRead : involvedReads) {
//     //         if (visitedInComponent.count(startNodeRead)) {
//     //             continue; // Already processed in this component
//     //         }

//     //         // Start a new sub-component phasing run
//     //         std::set<OrientedReadId> currentSubHaplotype1;
//     //         std::set<OrientedReadId> currentSubHaplotype2;

//     //         // Start with the current unvisited read in Haplotype 1 of the sub-component
//     //         currentSubHaplotype1.insert(startNodeRead);
//     //         visitedInComponent.insert(startNodeRead);
//     //         q.push(startNodeRead);

//     //         while(!q.empty()) {
//     //             OrientedReadId currentRead = q.front();
//     //             q.pop();
//     //             size_t u = readToIndex[currentRead];
//     //             bool isInSubHap1 = currentSubHaplotype1.count(currentRead);

//     //             PhasingGraph::adjacency_iterator neighborIt, neighborEnd;
//     //             for (boost::tie(neighborIt, neighborEnd) = boost::adjacent_vertices(u, phasingGraph); neighborIt != neighborEnd; ++neighborIt) {
//     //                 size_t v_idx = *neighborIt;
//     //                 OrientedReadId neighborRead = indexToRead[v_idx];

//     //                 if (visitedInComponent.count(neighborRead)) {
//     //                     // Consistency Check (Optional but recommended)
//     //                     bool neighborInSubHap1 = currentSubHaplotype1.count(neighborRead);
//     //                     bool neighborInSubHap2 = currentSubHaplotype2.count(neighborRead);
//     //                     if (neighborInSubHap1 || neighborInSubHap2) {
//     //                          auto edge_desc_check = boost::edge(u, v_idx, phasingGraph);
//     //                          if(edge_desc_check.second) {
//     //                              double weight_check = weightMap[edge_desc_check.first];
//     //                              bool consistent = (weight_check > 0 && ((isInSubHap1 && neighborInSubHap1) || (!isInSubHap1 && neighborInSubHap2))) ||
//     //                                                (weight_check < 0 && ((isInSubHap1 && neighborInSubHap2) || (!isInSubHap1 && neighborInSubHap1)));
//     //                              if (!consistent) {
//     //                                  cout << "Warning: Phasing inconsistency detected in component " << componentId << " between " << currentRead << " and " << neighborRead << ". Weight: " << weight_check << endl;
//     //                                  // Handle inconsistency: e.g., mark both reads as unphased
//     //                                  componentUnphased.insert(currentRead);
//     //                                  componentUnphased.insert(neighborRead);
//     //                                  // Remove from current sub-haplotypes if present
//     //                                  currentSubHaplotype1.erase(currentRead);
//     //                                  currentSubHaplotype2.erase(currentRead);
//     //                                  currentSubHaplotype1.erase(neighborRead);
//     //                                  currentSubHaplotype2.erase(neighborRead);
//     //                              }
//     //                          }
//     //                     }
//     //                     continue; // Already assigned or handled
//     //                 }

//     //                 auto edge_desc = boost::edge(u, v_idx, phasingGraph);
//     //                 if (!edge_desc.second) continue;
//     //                 double weight = weightMap[edge_desc.first];

//     //                 // Assign neighbor based on connection
//     //                 if (weight > 0) { // Agreement
//     //                     if (isInSubHap1) currentSubHaplotype1.insert(neighborRead);
//     //                     else currentSubHaplotype2.insert(neighborRead);
//     //                 } else { // Conflict
//     //                     if (isInSubHap1) currentSubHaplotype2.insert(neighborRead);
//     //                     else currentSubHaplotype1.insert(neighborRead);
//     //                 }
//     //                 visitedInComponent.insert(neighborRead);
//     //                 q.push(neighborRead);
//     //             }
//     //         }
//     //         // Merge sub-component results into component results
//     //         componentHaplotype1.insert(currentSubHaplotype1.begin(), currentSubHaplotype1.end());
//     //         componentHaplotype2.insert(currentSubHaplotype2.begin(), currentSubHaplotype2.end());
//     //     }

//     //     // Handle reads that were involved but never visited (singletons in phasing graph)
//     //     for(const auto& read : involvedReads) {
//     //         if (!visitedInComponent.count(read)) {
//     //             componentUnphased.insert(read);
//     //         }
//     //     }

//     //     // Remove any reads marked as unphased due to inconsistencies from final haplotypes
//     //     for(const auto& unphasedRead : componentUnphased) {
//     //         componentHaplotype1.erase(unphasedRead);
//     //         componentHaplotype2.erase(unphasedRead);
//     //     }


//     //     cout << "Component " << componentId << " phased into: Hap1=" << componentHaplotype1.size()
//     //          << ", Hap2=" << componentHaplotype2.size() << ", Unphased=" << componentUnphased.size() << endl;
//     //     totalPhasedReads += componentHaplotype1.size() + componentHaplotype2.size();


//     //     // 4. Apply Phasing Results (Update forbiddenAlignments)
//     //     // Get ReadIds for each haplotype set
//     //     auto getReadIds = [](const std::set<OrientedReadId>& orientedReads) {
//     //         std::set<ReadId> readIds;
//     //         for (const auto& orId : orientedReads) readIds.insert(orId.getReadId());
//     //         return readIds;
//     //     };
//     //     std::set<ReadId> hap1ReadIds = getReadIds(componentHaplotype1);
//     //     std::set<ReadId> hap2ReadIds = getReadIds(componentHaplotype2);
//     //     std::set<ReadId> unphasedReadIds = getReadIds(componentUnphased);


//     //     // Iterate through all alignments to forbid cross-haplotype ones within this component
//     //     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//     //         // Skip if already forbidden
//     //         if (forbiddenAlignments[alignmentId]) {
//     //             continue;
//     //         }

//     //         const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//     //         ReadId alnReadId0 = thisAlignmentData.readIds[0];
//     //         ReadId alnReadId1 = thisAlignmentData.readIds[1];

//     //         // Check if both reads belong to this component's phased sets
//     //         bool read0InComp = hap1ReadIds.count(alnReadId0) || hap2ReadIds.count(alnReadId0) || unphasedReadIds.count(alnReadId0);
//     //         bool read1InComp = hap1ReadIds.count(alnReadId1) || hap2ReadIds.count(alnReadId1) || unphasedReadIds.count(alnReadId1);

//     //         if (read0InComp && read1InComp) {
//     //             // Both reads are in this component's scope. Check phasing.
//     //             bool read0Hap1 = hap1ReadIds.count(alnReadId0);
//     //             bool read0Hap2 = hap2ReadIds.count(alnReadId0);
//     //             bool read1Hap1 = hap1ReadIds.count(alnReadId1);
//     //             bool read1Hap2 = hap2ReadIds.count(alnReadId1);
//     //             bool read0Unphased = unphasedReadIds.count(alnReadId0);
//     //             bool read1Unphased = unphasedReadIds.count(alnReadId1);

//     //             // Forbid if reads are in different haplotypes
//     //             if ((read0Hap1 && read1Hap2) || (read0Hap2 && read1Hap1)) {
//     //                 // cout << "Forbidding alignment " << alignmentId << " (reads " << alnReadId0 << ", " << alnReadId1 << ") due to cross-haplotype connection in component " << componentId << endl;
//     //                 forbiddenAlignments[alignmentId] = true;
//     //                 // Optionally unmark from firstPassHetAlignments if it was marked
//     //                 if (firstPassHetAlignments[alignmentId]) {
//     //                     firstPassHetAlignments[alignmentId] = false;
//     //                 }
//     //             }
//     //             // Forbid if either read is unphased (conservative approach)
//     //             else if (read0Unphased || read1Unphased) {
//     //                 // cout << "Forbidding alignment " << alignmentId << " (reads " << alnReadId0 << ", " << alnReadId1 << ") due to unphased read in component " << componentId << endl;
//     //                 forbiddenAlignments[alignmentId] = true;
//     //                 if (firstPassHetAlignments[alignmentId]) {
//     //                     firstPassHetAlignments[alignmentId] = false;
//     //                 }
//     //             }
                    
//     //             // If reads are in the same haplotype, confirm marking as first pass.
//     //             else if ((read0Hap1 && read1Hap1) || (read0Hap2 && read1Hap2)) {
//     //                 // This alignment is intra-haplotype within the component according to global phasing.
//     //                 // Mark or re-confirm it as a first pass alignment.
//     //                 // This might override previous decisions if the global phasing provides a clearer picture.
//     //                 if (!forbiddenAlignments[alignmentId]) { // Only mark if not already forbidden for other reasons
//     //                     firstPassHetAlignments[alignmentId] = true;
//     //                     // Optionally add a log message if needed:
//     //                     // cout << "Confirming alignment " << alignmentId << " (reads " << alnReadId0 << ", " << alnReadId1 << ") as intra-haplotype in component " << componentId << endl;
//     //                 }
//     //             }
//     //         }
//     //     }
//     //     // --- End Apply Phasing Results ---

//     // }

//     // cout << timestamp << "Finished global phasing of " << totalComponentsPhased << " connectedComponents. Total reads phased: " << totalPhasedReads << endl;
//     // // --- Global Phasing Implementation Ends ---






//     //
//     //
//     // We finished looping over the reads and found the haplotype sets for each one
//     // Now we need to create the first pass read graph that involves only those haplotype specific 
//     // alignments we found in the het sites
//     //
//     //

    
//     //*
//     //
//     // Order alignments in order of increasing Q. 
//     //
//     // Gather in alignmentTable[alignmentID, Q]
//     // alignments in order of increasing Q.
//     // Q(n) = (1 + δ/2ε)^n * e-δL
//     // ε = 1e-4, δ = 5e-4
//     // logQ(n) = αn - δL
//     //
//     //*

//     // const double epsilon = 1e-4;
//     // const double delta = 5e-4;
//     const double alpha = log(1 + delta/(2*epsilon));

//     // const double WThreshold = 1e-8;
//     const double logWThreshold = log(WThreshold);

//     // const double WThresholdForBreaks = 1e+15;
//     const double logWThresholdForBreaks = log(WThresholdForBreaks);

//     vector< pair<uint64_t, double> > alignmentTableHetSites;
    
//     // Keep track of which readIds were used in alignments
//     vector<bool> readUsed(readCount, false);



//     //
//     //
//     //
//     //
//     // Do a first pass in which we allow only het loci in order of increasing Q, 
//     // then do another pass to fill the breaks where you also allow all the other alignments
//     //
//     //
//     //
//     //



//     // Loop over all alignments.
//     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        
//         if (!firstPassHetAlignments[alignmentId]) {
//             continue;
//         }

//         if (forbiddenAlignments[alignmentId]) {
//             continue;
//         }

//         // Get information for this alignment.
//         AlignmentData& thisAlignmentData = alignmentData[alignmentId];

//         // The alignment is stored as an alignment between readId0 on strand 0
//         // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         // The reverse complement alignment also exists, but is not stored explicitly.
        
//         const ReadId readId0 = thisAlignmentData.readIds[0];
//         const ReadId readId1 = thisAlignmentData.readIds[1];
//         const bool isSameStrand = thisAlignmentData.isSameStrand;
//         SHASTA_ASSERT(readId0 < readId1);
//         const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
//         const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

//         // Skip if either read involved in the alignment is forbidden
//         if (forbiddenReads[readId0] || forbiddenReads[readId1]) {
//             continue;
//         }

//         // if (bridgingReads.count(readId0) || bridgingReads.count(readId1)) {
//         //     continue;
//         // }

//         // Store this pair of edges in our edgeTable.
//         const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
//         const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
//         const double L = double(range0 + range1)/2.;
//         // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
//         const double errorRateRle = thisAlignmentData.info.errorRateRle;
//         const double nRLE = errorRateRle * 2 * L;
//         // const double markerCount = thisAlignmentData.info.markerCount;

//         // logQ(n) = αn - δL
//         const double logQ = alpha * double(nRLE) - delta * L;

//         // Add the alignment to the table along with its logQ.
//         alignmentTableHetSites.push_back(make_pair(alignmentId, logQ));
//         readUsed[readId0] = true;
//         readUsed[readId1] = true;

//     }

//     sort(alignmentTableHetSites.begin(), alignmentTableHetSites.end(), OrderPairsBySecondOnly<uint64_t, double>());
//     cout << "The alignmentTableHetSites has " << alignmentTableHetSites.size() << " entries." << endl;




//     // Maintain a vector containing the degree of each vertex
//     // verticesDegree[vertexID] -> degree
//     vector<uint64_t> verticesDegree(orientedReadCount, 0);

//     cout << "Number of reads: " << readCount << endl;
//     cout << "Number of oriented reads: " << orientedReadCount << endl;


//     ///
//     // Process the HET SITES ALIGNMENTS in order of increasing Q. 
//     //
//     // i.   Start with no edges in the read graph. 
//     // ii.  Process alignments spanning the het sites in order of increasing Q. 
//     // iii. If the alignment breaks strand separation, it is skipped. 
//     // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped. (this step is skipped)  
//     // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.

//     // Initiallize disjoint sets for HET SITES ALIGNMENTS
//     vector<ReadId> rank(orientedReadCount);
//     vector<ReadId> parent(orientedReadCount);
//     boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
//     for(ReadId readId=0; readId<readCount; readId++) {
//         for(Strand strand=0; strand<2; strand++) {
//             disjointSets.make_set(OrientedReadId(readId, strand).getValue());
//         }
//     }

//     // Flag all alignments as not to be kept.
//     vector<bool> keepAlignment(alignmentCount, false);

//     // Process alignments in order of increasing Q
//     vector alignmentTablesToProcess({alignmentTableHetSites});
    
//     uint64_t crossStrandEdgeCountHet = 0;

//     for (auto alignmentTableToProcess : alignmentTablesToProcess) {
//         for(auto it=alignmentTableToProcess.begin(); it!=alignmentTableToProcess.end(); ++it) {
//             const pair<uint64_t, double>& p = *it;
//             const uint64_t alignmentId = p.first;
//             // const double logQ = p.second;

//             // Get the alignment data
//             AlignmentData& alignment = alignmentData[alignmentId];
//             const ReadId readId0 = alignment.readIds[0];
//             const ReadId readId1 = alignment.readIds[1];
//             const bool isSameStrand = alignment.isSameStrand;
//             SHASTA_ASSERT(readId0 < readId1);
//             const OrientedReadId A0 = OrientedReadId(readId0, 0);
//             const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
//             const OrientedReadId A1 = OrientedReadId(readId0, 1);
//             const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             // If the alignment breaks strand separation, it is skipped.
//             // If A0 and B1 are in the same connected component,
//             // A1 and B0 also must be in the same connected component.
//             // Adding this pair of edges would create a self-complementary
//             // connected component containing A0, B0, A1, and B1,
//             // and to ensure strand separation we don't want to do that.
//             // So we mark these edges as cross-strand edges
//             // and don't use them to update the disjoint set data structure.
//             if(a0 == b1) {
//                 SHASTA_ASSERT(a1 == b0);
//                 crossStrandEdgeCountHet += 2;
//                 continue;
//             }

            

//             // If both vertices of the potential edge have at least the required minimum number 
//             // of neighbors, the alignment is also skipped. 
//             const uint64_t degreeA0 = verticesDegree[A0.getValue()];
//             const uint64_t degreeB0 = verticesDegree[B0.getValue()];



//             if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
//                 // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
//                 continue;
//             }

//             // Add the alignment to the read graph.
//             keepAlignment[alignmentId] = true;
//             alignment.info.isInReadGraph = 1;
            

//             // Update vertex degrees
//             verticesDegree[A0.getValue()]++;
//             verticesDegree[B0.getValue()]++;
//             verticesDegree[A1.getValue()]++;
//             verticesDegree[B1.getValue()]++;
            

//             // Update disjoint sets
//             disjointSets.union_set(a0, b0);
//             disjointSets.union_set(a1, b1);

//         }

//         // Verify that for any read the two oriented reads are in distinct
//         // connected components.
//         for(ReadId readId=0; readId<readCount; readId++) {
//             const OrientedReadId orientedReadId0(readId, 0);
//             const OrientedReadId orientedReadId1(readId, 1);
//             SHASTA_ASSERT(
//                 disjointSets.find_set(orientedReadId0.getValue()) !=
//                 disjointSets.find_set(orientedReadId1.getValue())
//             );
//         }
//     }

//     // Print how many alignments were kept in this step
//     const long keepCountR1 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Finding strict disjointSets step involving alignments in HET SITES: Keeping " << keepCountR1 << " alignments in HET sites out of " << keepAlignment.size() << " alignments in general."<< endl;
//     cout << "Found " << crossStrandEdgeCountHet << " cross strand edges in HET SITES during disjointSets construction." << endl;


//     //*
//     //
//     // Create the dynamically adjustable boost readGraph using the alignments we selected.
//     //
//     //*
//     using boost::add_vertex;
//     using boost::add_edge;

//     // The vertex_descriptor is OrientedReadId::getValue().
//     ReadGraph4 readGraph(orientedReadCount);

//     // Initially, each alignment generates two edges.
//     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//         // Record whether this alignment is used in the read graph.
//         const bool keepThisAlignment = keepAlignment[alignmentId];
//         const AlignmentData& alignment = alignmentData[alignmentId];

//         // If this alignment is not used in the read graph, we are done.
//         if(!keepThisAlignment) {
//             continue;
//         }

//         // Get the OrientedReadIds.
//         OrientedReadId orientedReadId0(alignment.readIds[0], 0);
//         OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//         SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

//         // Swap them if necessary, depending on the average alignment offset at center.
//         if(alignment.info.offsetAtCenter() < 0.) {
//             swap(orientedReadId0, orientedReadId1);
//         }

//         // Create the edge.
//         add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//         // Also create the reverse complemented edge.
//         orientedReadId0.flipStrand();
//         orientedReadId1.flipStrand();
//         add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
//     }

//     cout << "The read graph that used the alignments in HET sites has " << num_vertices(readGraph) << " vertices and " << num_edges(readGraph) << " edges." << endl;



//     // UNCOMMENT
//     vector< pair<uint64_t, double> > alignmentTable;
//     vector< pair<uint64_t, double> > alignmentTableNotPassFilter;
//     // Flag alignments to be kept for break detection.
//     vector<bool> keepAlignmentsForBreaks(alignmentCount, false);

//     //
//     //
//     //
//     //
//     // Do a second pass in which we allow all the other alignments in order of increasing Q
//     // after we filter them out using the Bayesian filtering
//     //
//     //
//     //
//     //

//     // Loop over all alignments.
//     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         if((alignmentId % 100000) == 0) {
//             cout << timestamp << alignmentId << "/" << alignmentCount << endl;
//         }

//         if (alignmentsAlreadyConsidered[alignmentId]) {
//             continue;
//         }

//         if (forbiddenAlignments[alignmentId]) {
//             continue;
//         }

//         // Get information for this alignment.
//         AlignmentData& thisAlignmentData = alignmentData[alignmentId];

//         // The alignment is stored as an alignment between readId0 on strand 0
//         // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         // The reverse complement alignment also exists, but is not stored explicitly.
//         const ReadId readId0 = thisAlignmentData.readIds[0];
//         const ReadId readId1 = thisAlignmentData.readIds[1];
//         const bool isSameStrand = thisAlignmentData.isSameStrand;
//         SHASTA_ASSERT(readId0 < readId1);
//         const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
//         const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

//         // Skip if either read involved in the alignment is forbidden
//         if (forbiddenReads[readId0] || forbiddenReads[readId1]) {
//             continue;
//         }

//         // Store this pair of edges in our edgeTable.
//         const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
//         const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
//         const double L = double(range0 + range1)/2.;
//         // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
//         const double errorRateRle = thisAlignmentData.info.errorRateRle;
//         const double nRLE = errorRateRle * 2 * L;
//         // const double markerCount = thisAlignmentData.info.markerCount;

//         // logQ(n) = αn - δL
//         const double logQ = alpha * double(nRLE) - delta * L;

//         // This time use the regular Bayesian filtering
//         if (logQ <= logWThreshold) {
//             alignmentTable.push_back(make_pair(alignmentId, logQ));
//             alignmentsAlreadyConsidered[alignmentId] = true;
//             readUsed[readId0] = true;
//             readUsed[readId1] = true;
//         } else if(logQ <= logWThresholdForBreaks){
//             alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
//             alignmentsAlreadyConsidered[alignmentId] = true;
//             keepAlignmentsForBreaks[alignmentId] = true;
//         }

//     }

//     sort(alignmentTable.begin(), alignmentTable.end(), OrderPairsBySecondOnly<uint64_t, double>());
//     sort(alignmentTableNotPassFilter.begin(), alignmentTableNotPassFilter.end(), OrderPairsBySecondOnly<uint64_t, double>());
//     cout << "The alignmentTable has " << alignmentTable.size() << " entries." << endl;
//     cout << "The alignmentTableNotPassFilter has " << alignmentTableNotPassFilter.size() << " entries." << endl;





//     ///
//     // Process alignments in order of increasing Q. 
//     //
//     // i.   Start with no edges in the read graph. 
//     // ii.  Process alignments in order of increasing Q. 
//     // iii. If the alignment breaks strand separation, it is skipped. 
//     // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped.  
//     // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.
    

//     // boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
//     // for(ReadId readId=0; readId<readCount; readId++) {
//     //     for(Strand strand=0; strand<2; strand++) {
//     //         disjointSets.make_set(OrientedReadId(readId, strand).getValue());
//     //     }
//     // }

//     // // Flag all alignments as not to be kept.
//     // vector<bool> keepAlignment(alignmentCount, false);

//     // Process alignments in order of increasing Q
//     alignmentTablesToProcess = {alignmentTable};

//     uint64_t crossStrandEdgeCount = 0;

//     for (auto alignmentTableToProcess : alignmentTablesToProcess) {
//         for(auto it=alignmentTableToProcess.begin(); it!=alignmentTableToProcess.end(); ++it) {
//             const pair<uint64_t, double>& p = *it;
//             const uint64_t alignmentId = p.first;
//             // const double logQ = p.second;

//             // Get the alignment data
//             AlignmentData& alignment = alignmentData[alignmentId];
//             const ReadId readId0 = alignment.readIds[0];
//             const ReadId readId1 = alignment.readIds[1];
//             const bool isSameStrand = alignment.isSameStrand;
//             SHASTA_ASSERT(readId0 < readId1);
//             const OrientedReadId A0 = OrientedReadId(readId0, 0);
//             const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
//             const OrientedReadId A1 = OrientedReadId(readId0, 1);
//             const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             // If the alignment breaks strand separation, it is skipped.
//             // If A0 and B1 are in the same connected component,
//             // A1 and B0 also must be in the same connected component.
//             // Adding this pair of edges would create a self-complementary
//             // connected component containing A0, B0, A1, and B1,
//             // and to ensure strand separation we don't want to do that.
//             // So we mark these edges as cross-strand edges
//             // and don't use them to update the disjoint set data structure.
//             if(a0 == b1) {
//                 SHASTA_ASSERT(a1 == b0);
//                 crossStrandEdgeCount += 2;
//                 continue;
//             }

            

//             // If both vertices of the potential edge have at least the required minimum number 
//             // of neighbors, the alignment is also skipped. 
//             const uint64_t degreeA0 = verticesDegree[A0.getValue()];
//             const uint64_t degreeB0 = verticesDegree[B0.getValue()];



//             // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
//             //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
//             //     continue;
//             // }

//             // Add the alignment to the read graph.
//             keepAlignment[alignmentId] = true;
//             alignment.info.isInReadGraph = 1;

//             // Update vertex degrees
//             verticesDegree[A0.getValue()]++;
//             verticesDegree[B0.getValue()]++;
//             verticesDegree[A1.getValue()]++;
//             verticesDegree[B1.getValue()]++;
            

//             // Update disjoint sets
//             disjointSets.union_set(a0, b0);
//             disjointSets.union_set(a1, b1);

//         }

//         // Verify that for any read the two oriented reads are in distinct
//         // connected components.
//         for(ReadId readId=0; readId<readCount; readId++) {
//             const OrientedReadId orientedReadId0(readId, 0);
//             const OrientedReadId orientedReadId1(readId, 1);
//             SHASTA_ASSERT(
//                 disjointSets.find_set(orientedReadId0.getValue()) !=
//                 disjointSets.find_set(orientedReadId1.getValue())
//             );
//         }
//     }

//     // Print how many alignments were kept in this step
//     const long keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Finding strict disjointSets step: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;






//     //*
//     //
//     // Update the dynamically adjustable boost readGraph using the alignments we selected not involving HET sites.
//     //
//     //*


//     // Initially, each alignment generates two edges.
//     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//         // Record whether this alignment is used in the read graph.
//         const bool keepThisAlignment = keepAlignment[alignmentId];
//         const AlignmentData& alignment = alignmentData[alignmentId];

//         // If this alignment is not used in the read graph, we are done.
//         if(not keepThisAlignment) {
//             continue;
//         }

//         // Get the OrientedReadIds.
//         OrientedReadId orientedReadId0(alignment.readIds[0], 0);
//         OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//         SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

//         // Swap them if necessary, depending on the average alignment offset at center.
//         if(alignment.info.offsetAtCenter() < 0.) {
//             swap(orientedReadId0, orientedReadId1);
//         }

//         // Create the edge.
//         add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//         // Also create the reverse complemented edge.
//         orientedReadId0.flipStrand();
//         orientedReadId1.flipStrand();
//         add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
//     }

//     cout << "The read graph has " << num_vertices(readGraph) << " vertices and " << num_edges(readGraph) << " edges." << endl;



//     //*
//     //
//     // Create the dynamically adjustable boost readGraph using the alignments that did not pass the strict filter.
//     // These alignments are used to create the read graph that will aid in the break detection.
//     //
//     //*
    
//     // The vertex_descriptor is OrientedReadId::getValue().
//     ReadGraph4AllAlignments readGraphAllAlignments(orientedReadCount);

//     // Initially, each alignment generates two edges.
//     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//         // Record whether this alignment is used in the read graph.
//         const bool keepThisAlignment = keepAlignmentsForBreaks[alignmentId];
//         const AlignmentData& alignment = alignmentData[alignmentId];

//         // If this alignment is not used in the read graph, we are done.
//         if(not keepThisAlignment) {
//             continue;
//         }

//         // Get the OrientedReadIds.
//         OrientedReadId orientedReadId0(alignment.readIds[0], 0);
//         OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//         SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

//         // Swap them if necessary, depending on the average alignment offset at center.
//         if(alignment.info.offsetAtCenter() < 0.) {
//             swap(orientedReadId0, orientedReadId1);
//         }

//         // Create the edge.
//         add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

//         // Also create the reverse complemented edge.
//         orientedReadId0.flipStrand();
//         orientedReadId1.flipStrand();
//         add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
//     }
    
//     cout << "The read graph for break detection has " << num_vertices(readGraphAllAlignments) << " vertices and " << num_edges(readGraphAllAlignments) << " edges." << endl;


//     //*
//     //
//     // Find possible dead end nodes on the directed graph.
//     // We only keep the reads that have no outgoing neighbors.
//     //
//     //*
//     vector<bool> potentialDeadEndReads(orientedReadCount, false);
//     vector<bool> isolatedReads(orientedReadCount, false);

//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {
//             OrientedReadId orientedReadId(readId, strand);
            
//             // Find neighbors in the forward direction
//             vector<OrientedReadId> forwardNeighbors;
//             readGraph.findNeighborsDirectedGraphOneSideRight(orientedReadId, 1, forwardNeighbors);
            
//             // Find neighbors in the backward direction 
//             vector<OrientedReadId> leftNeighbors;
//             readGraph.findNeighborsDirectedGraphOneSideLeft(orientedReadId, 1, leftNeighbors);

//             if (forwardNeighbors.empty() && leftNeighbors.empty() ) {
//                 isolatedReads[orientedReadId.getValue()] = true;
//                 OrientedReadId reverseOrientedReadId = orientedReadId;
//                 reverseOrientedReadId.flipStrand();
//                 isolatedReads[reverseOrientedReadId.getValue()] = true;
//                 continue;
//             }

//             // If a read has neighbors only in the backward direction, it's a potential dead end
//             if (forwardNeighbors.empty() && !leftNeighbors.empty()) {
//                 potentialDeadEndReads[orientedReadId.getValue()] = true;
//             }
//         }
//     }

//     // count the number of potential dead end reads
//     long potentialDeadEndReadCount = count(potentialDeadEndReads.begin(), potentialDeadEndReads.end(), true);
//     cout << "Found " << potentialDeadEndReadCount << " potential dead end reads." << endl;

//     // // Print dead end Oriented reads
//     // // iterate over all oriented reads
//     // for (ReadId readId = 0; readId < readCount; readId++) {
//     //     for (Strand strand = 0; strand < 2; strand++) {
//     //         OrientedReadId orientedReadId(readId, strand);
//     //         if (potentialDeadEndReads[orientedReadId.getValue()]) {
//     //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a potential dead end read." << endl;
//     //         }
//     //     }
//     // }






//     //*
//     //
//     // Filter out the potential dead end nodes on the directed graph.
//     // Look for a read path that lead to a positive offset.
//     // If no such path is found, the orientedReadId is kept as a potential dead end node.
//     //
//     //*
//     uint64_t finalNumberOfPotentialDeadEndNodes = 0;
//     vector<bool> finalDeadEndReadsWithNoOutgoingNodes(orientedReadCount, false);
//     vector<bool> finalDeadEndReadsWithNoIncomingNodes(orientedReadCount, false);
//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {

//             OrientedReadId orientedReadId(readId, strand);

//             OrientedReadId reverseOrientedReadId = orientedReadId;
//             reverseOrientedReadId.flipStrand();

//             // check if the orientedReadId is in potentialDeadEndReads
//             if(!potentialDeadEndReads[orientedReadId.getValue()]) {
//                 continue;
//             }

//             // create the necessary variables for findAllPaths
//             vector<vector<OrientedReadId>> paths;
//             vector<vector<double>> pathsOffsets;
//             vector<OrientedReadId> currentPath;
//             vector<double> currentPathOffset;
//             std::set<ReadGraph4BaseClass::vertex_descriptor> visited;
//             uint64_t maxDistance = 4;
//             uint64_t currentDistance = 0;

//             bool result = readGraph.findPathWithPositiveOffset(orientedReadId, paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

//             // Check if we found a read path with positive offset.
//             // If yes, the function findPathWithPositiveOffset will return 1, if not, it will return 0.
//             if(result == 1) {
//                 // cout << "Found a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset" << endl;
//                 // //print the paths and then the pathsOffsets
//                 // for (uint64_t i = 0; i < paths.size(); i++) {
//                 //     cout << "Path " << i << ": ";
//                 //     for (uint64_t j = 0; j < paths[i].size(); j++) {
//                 //         cout << paths[i][j].getReadId() << " ";
//                 //     }
//                 //     cout << endl;
//                 //     cout << "PathOffsets " << i << ": ";
//                 //     for (uint64_t j = 0; j < pathsOffsets[i].size(); j++) {
//                 //         cout << pathsOffsets[i][j] << " ";
//                 //     }
//                 //     cout << endl;
//                 // }
//             } else if (result == 0) {
//                 finalNumberOfPotentialDeadEndNodes++;
//                 // cout << "Did not find a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset. Keeping it as a potential dead end read." << endl;
//                 finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()] = true;
//                 finalDeadEndReadsWithNoIncomingNodes[reverseOrientedReadId.getValue()] = true;
//             }
//         }
//     }

//     cout << "After filtering we are left with " << finalNumberOfPotentialDeadEndNodes << " potential dead end reads." << endl;


//     // // print dead end Oriented reads
//     // // iterate over all oriented reads
//     // for (ReadId readId = 0; readId < readCount; readId++) {
//     //     for (Strand strand = 0; strand < 2; strand++) {
//     //         OrientedReadId orientedReadId(readId, strand);
//     //         if (finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]) {
//     //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
//     //         } else if (finalDeadEndReadsWithNoIncomingNodes[orientedReadId.getValue()]) {
//     //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
//     //         }
//     //     }
//     // }











//     // //*
//     // //
//     // // Extend the potential dead end nodes list.
//     // // Add neighboring nodes of potential dead end nodes to the dead end node list.
//     // //
//     // //*
//     // vector<bool> finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors(orientedReadCount, false);
//     // vector<bool> finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors(orientedReadCount, false);
//     // for (ReadId readId = 0; readId < readCount; readId++) {
//     //     for (Strand strand = 0; strand < 2; strand++) {
//     //         OrientedReadId orientedReadId(readId, strand);

//     //         if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){

//     //             finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()] = true;
                
//     //             OrientedReadId reverseOrientedReadId = orientedReadId;
//     //             reverseOrientedReadId.flipStrand();

//     //             finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;

//     //             // Find distance 1 neighbors
//     //             vector<OrientedReadId> distance1Neighbors;
//     //             readGraph.findNeighborsDirectedGraphBothSides(orientedReadId, 1, distance1Neighbors);

//     //             for(auto neighbor : distance1Neighbors) {
//     //                 finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[neighbor.getValue()] = true;

//     //                 OrientedReadId reverseOrientedReadId = neighbor;
//     //                 reverseOrientedReadId.flipStrand();

//     //                 finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;
//     //             }
//     //         }
//     //     }
//     // }

//     // // // print dead end Oriented reads
//     // // // iterate over all oriented reads
//     // // for (ReadId readId = 0; readId < readCount; readId++) {
//     // //     for (Strand strand = 0; strand < 2; strand++) {
//     // //         OrientedReadId orientedReadId(readId, strand);
//     // //         if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
//     // //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
//     // //         } else if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
//     // //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
//     // //         }
//     // //     }
//     // // }





//     //*
//     //
//     // Create a map of endNodesWithNoOutgoingNodes to other endNodesWithNoIncomingNodes that are possible to connect to.
//     //
//     //*

//     // Case 1: the NoOut deadEnd node will be mapped to at least 3 NoIn deadEnd nodes 
//     // IN THE SAME disjointSet.
//     //
//     //  NoIn - - - - | - - - - -  - - - - - - | - - - - -
//     //  NoIn - - - - | - - NoOut     NoIn - - | - - - - -

//     std::unordered_map<uint64_t, vector<uint64_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet;
//     vector<bool> endNodesWithNoOutgoingNodesConsidered(orientedReadCount, false);

//     // First, we try to find EndNodesWithNoIncomingNodes that are in the same connected component
//     // as EndNodesWithNoOutgoingNodes. These have priority over other EndNodesWithNoIncomingNodes in other connected components.
//     // These will also include telomeric nodes because they do not have incoming nodes!
//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {
//             OrientedReadId orientedReadId(readId, strand);

//             // check if the orientedReadId is an endNode with no outgoing nodes
//             if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
//                 // Get it's disjointSet
//                 uint64_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

//                 // Find all EndNodesWithNoIncomingNodes that are in the same disjointSet
//                 for(uint64_t id=0; id<orientedReadCount; id++) {
//                     // check if the id is an endNode with no incoming nodes
//                     if(finalDeadEndReadsWithNoIncomingNodes[id]){
//                         // get the disjointSet of the id
//                         uint64_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
//                         // check if the disjointSet of the id is the same as the disjointSet of the endNode with no outgoing nodes
//                         if(deadEndReadWithNoOutgoingNodesDisjointSetId == endNodeWithNoIncomingNodesDisjointSetId){
//                             endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet[orientedReadId.getValue()].push_back(id);
//                             endNodesWithNoOutgoingNodesConsidered[orientedReadId.getValue()] = true;
//                         }
//                     }
//                 }

//             }
//         }
//     }


//     for(auto& p : endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet) {
//         uint64_t value = p.first;
//         OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(value));
//         SHASTA_ASSERT(orientedReadId.getValue() == value);
//         // const ReadId readId = orientedReadId.getReadId();
//         // const Strand strand = orientedReadId.getStrand();

//         // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in the same disjointSet:" << endl;
//         vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
//         for(auto& node : p.second) {
//             OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(ReadId(node));
//             // cout << "ReadID " << nodeOrientedReadId.getReadId() << " strand " << nodeOrientedReadId.getStrand() << endl;
//             SHASTA_ASSERT(nodeOrientedReadId.getValue() == node);
//             endNodesWithNoIncomingNodes[node] = true;
//         }

        
//         // Find neighbors in the forward direction of the ALL ALIGNMENTS read graph 
//         // starting from orientedReadId which is an endNode with no outgoing nodes in the filtered read graph. 
//         // Early stop when we reach an endNode with no incoming nodes (a node with endNodesWithNoIncomingNodes set to true).
//         vector<OrientedReadId> forwardNeighbors;
//         readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, endNodesWithNoIncomingNodes, 5, forwardNeighbors);

//         // // print the forward neighbors
//         // cout << "Forward neighbors of the endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " are:" << endl;
//         // for(auto& neighbor : forwardNeighbors) {
//         //     cout << "ReadID " << neighbor.getReadId() << " strand " << neighbor.getStrand() << endl;
//         // }


//         // create a std::set of the forwardNeighbors for easy contain check
//         std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
        

//         if(forwardNeighbors.empty()) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No forward neighbors found" << endl;
//             continue;
//         }


//         // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
//         // OR a non speficic node if we exceeded maxDistance
//         OrientedReadId lastNode = forwardNeighbors.back();


//         bool success = false;


//         // First check if the last node is a potential dead end node with no INCOMING nodes
//         if(endNodesWithNoIncomingNodes[lastNode.getValue()]) {

//             // cout << "The last node is: " << lastNode.getReadId() << " strand " << lastNode.getStrand() << endl;
            
//             const OrientedReadId A0 = orientedReadId;
//             const OrientedReadId B0 = lastNode;
//             const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
//             const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             // const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             // const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             // const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             // const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
//             // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//                 const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
//                 const uint64_t alignmentId = p.first;
//                 // const double logQ = p.second;

//                 const bool keepThisAlignment = keepAlignment[alignmentId];

//                 const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

//                 if(keepThisAlignment) {
//                     continue;
//                 }

//                 if(not keepThisBreaksAlignment) {
//                     continue;
//                 }

//                 AlignmentData& alignment = alignmentData[alignmentId];
            
//                 // Get the OrientedReadIds.
//                 OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
//                 OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//                 SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

//                 // Swap them if necessary, depending on the average alignment offset at center.
//                 if(alignment.info.offsetAtCenter() < 0.) {
//                     swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
//                 }


//                 if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                    
//                     // Get the alignment data
//                     ReadId readId0v2 = alignment.readIds[0];
//                     ReadId readId1v2 = alignment.readIds[1];
//                     const bool isSameStrandv2 = alignment.isSameStrand;
//                     SHASTA_ASSERT(readId0v2 < readId1v2);
//                     OrientedReadId A0v2 = OrientedReadId(readId0v2, 0);
//                     OrientedReadId B0v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 0 : 1);
//                     OrientedReadId A1v2 = OrientedReadId(readId0v2, 1);
//                     OrientedReadId B1v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 1 : 0);

//                     SHASTA_ASSERT(A0v2.getReadId() == A1v2.getReadId());
//                     SHASTA_ASSERT(B0v2.getReadId() == B1v2.getReadId());
//                     SHASTA_ASSERT(A0v2.getStrand() == 1 - A1v2.getStrand());
//                     SHASTA_ASSERT(B0v2.getStrand() == 1 - B1v2.getStrand());

//                     // Get the connected components that these oriented reads are in.
//                     const uint64_t a0v2 = disjointSets.find_set(A0v2.getValue());
//                     const uint64_t b0v2 = disjointSets.find_set(B0v2.getValue());
//                     const uint64_t a1v2 = disjointSets.find_set(A1v2.getValue());
//                     const uint64_t b1v2 = disjointSets.find_set(B1v2.getValue());


//                     // If the alignment breaks strand separation, it is skipped.
//                     // If A0 and B1 are in the same connected component,
//                     // A1 and B0 also must be in the same connected component.
//                     // Adding this pair of edges would create a self-complementary
//                     // connected component containing A0, B0, A1, and B1,
//                     // and to ensure strand separation we don't want to do that.
//                     // So we mark these edges as cross-strand edges
//                     // and don't use them to update the disjoint set data structure.
//                     if(a0v2 == b1v2) {
//                         SHASTA_ASSERT(a1v2 == b0v2);
//                         crossStrandEdgeCount += 2;
//                         continue;
//                     }

//                     // Add the alignment to the read graph.
//                     keepAlignment[alignmentId] = true;
//                     alignment.info.isInReadGraph = 1;

//                     // Update vertex degrees
//                     verticesDegree[A0v2.getValue()]++;
//                     verticesDegree[B0v2.getValue()]++;
//                     verticesDegree[A1v2.getValue()]++;
//                     verticesDegree[B1v2.getValue()]++;
                    

//                     // Update disjoint sets
//                     disjointSets.union_set(a0v2, b0v2);
//                     disjointSets.union_set(a1v2, b1v2);

//                     // Make sure all alignments added are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1v2.getValue()] = false;

//                     // Make sure the start and last nodes are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1.getValue()] = false;

//                     success = true;

//                     // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

//                     // Create the edge.
//                     add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//                     // Also create the reverse complemented edge.
//                     alignmentOrientedReadId0.flipStrand();
//                     alignmentOrientedReadId1.flipStrand();
//                     add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


//                 }

//             }

//             // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
//             // it means that we have connected them with alignments.
//             if(success) {
//                 // cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
//             }

//             // // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
//             // // it means that we have connected them with alignments.
//             // if(finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] == false and finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] == false) {
//             //     cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
//             // } else {
//             //     cout << "Did not connect endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // }

//         }

//         if(!success) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No alignments added" << endl;
//         }

            
//     }


//     // Print how many alignments were kept
//     const long keepCountR3 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Adding alignments for break bridging to connect endNodes in the same disjointSet: Keeping " << keepCountR3 << " alignments of " << keepAlignment.size() << endl;




    
//     // Case 2: the NoOut deadEnd node will be mapped to at least 2 NoIn deadEnd nodes in the same disjointSet.
//     // Other NoIn deadEnd nodes will be in a different disjointSet.
//     //
//     //  NoIn - - - - | - - - - -  - - - - - - - - - - -
//     //  NoIn - - - - | - - NoOut     
//     //                                     NoIn - - - - - - - (Other disjointSet)

//     // Case 3: Breaks in a haploid chromosome (chrX chrY).
//     //
//     //  NoIn - - - - - - NoOut
//     //                             NoIn - - - - - - - (Other disjointSet) 


//     std::unordered_map<uint64_t, vector<uint64_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet;

//     // Now, we try to find EndNodesWithNoIncomingNodes that are in different disjointSets than EndNodesWithNoOutgoingNodes.
//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {
//             OrientedReadId orientedReadId(readId, strand);

//             // check if the orientedReadId is an endNode with no outgoing nodes
//             if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
//                 // Get it's disjointSet
//                 uint64_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

//                 // Find all EndNodesWithNoIncomingNodes that are in different disjointSet
//                 for(uint64_t id=0; id<orientedReadCount; id++) {
//                     // check if the id is an endNode with no incoming nodes
//                     if(finalDeadEndReadsWithNoIncomingNodes[id]){
//                         // get the disjointSet of the id
//                         uint64_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
//                         // check if the disjointSet of the id is different from the disjointSet of the endNode with no outgoing nodes
//                         if(deadEndReadWithNoOutgoingNodesDisjointSetId != endNodeWithNoIncomingNodesDisjointSetId){
//                             endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet[orientedReadId.getValue()].push_back(id);
//                             endNodesWithNoOutgoingNodesConsidered[orientedReadId.getValue()] = true;
//                         }
//                     }
//                 }

//             }
//         }
//     }


//     for(auto& p : endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet) {
//         uint64_t value = p.first;
//         OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(value));
//         SHASTA_ASSERT(orientedReadId.getValue() == value);
//         // const ReadId readId = orientedReadId.getReadId();
//         // const Strand strand = orientedReadId.getStrand();

//         // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in a different disjointSet: " << endl;
//         vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
//         for(auto& node : p.second) {
//             OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(ReadId(node));
//             // cout << "ReadID " << nodeOrientedReadId.getReadId() << " strand " << nodeOrientedReadId.getStrand() << endl;
//             SHASTA_ASSERT(nodeOrientedReadId.getValue() == node);
//             endNodesWithNoIncomingNodes[node] = true;
//         }
        
//         // Find neighbors in the forward direction of the ALL ALIGNMENTS read graph 
//         // starting from orientedReadId which is an endNode with no outgoing nodes in the filtered read graph. 
//         // Early stop when we reach an endNode with no incoming nodes (a node with endNodesWithNoIncomingNodes set to true).
//         vector<OrientedReadId> forwardNeighbors;
//         readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, endNodesWithNoIncomingNodes, 5, forwardNeighbors);

//         // create a std::set of the forwardNeighbors for easy contain check
//         std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
        

//         if(forwardNeighbors.empty()) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No forward neighbors found" << endl;
//             continue;
//         }

//         // // print the forward neighbors
//         // cout << "Forward neighbors of the endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " are:" << endl;
//         // for(auto& neighbor : forwardNeighbors) {
//         //     cout << "ReadID " << neighbor.getReadId() << " strand " << neighbor.getStrand() << endl;
//         // }

//         // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
//         // OR a non speficic node if we exceeded maxDistance
//         OrientedReadId lastNode = forwardNeighbors.back();

//         bool success = false;

//         // First check if the last node is a potential dead end node with no INCOMING nodes
//         if(endNodesWithNoIncomingNodes[lastNode.getValue()]) {
            
//             const OrientedReadId A0 = orientedReadId;
//             const OrientedReadId B0 = lastNode;
//             const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
//             const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             // const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             // const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             // const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             // const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
//             // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//                 const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
//                 const uint64_t alignmentId = p.first;
//                 // const double logQ = p.second;

//                 const bool keepThisAlignment = keepAlignment[alignmentId];

//                 const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

//                 if(keepThisAlignment) {
//                     continue;
//                 }

//                 if(not keepThisBreaksAlignment) {
//                     continue;
//                 }

//                 AlignmentData& alignment = alignmentData[alignmentId];
            
//                 // Get the OrientedReadIds.
//                 OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
//                 OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//                 SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

//                 // Swap them if necessary, depending on the average alignment offset at center.
//                 if(alignment.info.offsetAtCenter() < 0.) {
//                     swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
//                 }


//                 if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                    
//                     // Get the alignment data
//                     ReadId readId0v2 = alignment.readIds[0];
//                     ReadId readId1v2 = alignment.readIds[1];
//                     const bool isSameStrandv2 = alignment.isSameStrand;
//                     SHASTA_ASSERT(readId0v2 < readId1v2);
//                     OrientedReadId A0v2 = OrientedReadId(readId0v2, 0);
//                     OrientedReadId B0v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 0 : 1);
//                     OrientedReadId A1v2 = OrientedReadId(readId0v2, 1);
//                     OrientedReadId B1v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 1 : 0);

//                     SHASTA_ASSERT(A0v2.getReadId() == A1v2.getReadId());
//                     SHASTA_ASSERT(B0v2.getReadId() == B1v2.getReadId());
//                     SHASTA_ASSERT(A0v2.getStrand() == 1 - A1v2.getStrand());
//                     SHASTA_ASSERT(B0v2.getStrand() == 1 - B1v2.getStrand());

//                     // Get the connected components that these oriented reads are in.
//                     const uint64_t a0v2 = disjointSets.find_set(A0v2.getValue());
//                     const uint64_t b0v2 = disjointSets.find_set(B0v2.getValue());
//                     const uint64_t a1v2 = disjointSets.find_set(A1v2.getValue());
//                     const uint64_t b1v2 = disjointSets.find_set(B1v2.getValue());


//                     // If the alignment breaks strand separation, it is skipped.
//                     // If A0 and B1 are in the same connected component,
//                     // A1 and B0 also must be in the same connected component.
//                     // Adding this pair of edges would create a self-complementary
//                     // connected component containing A0, B0, A1, and B1,
//                     // and to ensure strand separation we don't want to do that.
//                     // So we mark these edges as cross-strand edges
//                     // and don't use them to update the disjoint set data structure.
//                     if(a0v2 == b1v2) {
//                         SHASTA_ASSERT(a1v2 == b0v2);
//                         crossStrandEdgeCount += 2;
//                         continue;
//                     }

//                     // Add the alignment to the read graph.
//                     keepAlignment[alignmentId] = true;
//                     alignment.info.isInReadGraph = 1;

//                     // Update vertex degrees
//                     verticesDegree[A0v2.getValue()]++;
//                     verticesDegree[B0v2.getValue()]++;
//                     verticesDegree[A1v2.getValue()]++;
//                     verticesDegree[B1v2.getValue()]++;
                    

//                     // Update disjoint sets
//                     disjointSets.union_set(a0v2, b0v2);
//                     disjointSets.union_set(a1v2, b1v2);

//                     // Make sure all alignments added are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1v2.getValue()] = false;

//                     // Make sure the start and last nodes are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1.getValue()] = false;

//                     success = true;


//                     // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

//                     // Create the edge.
//                     add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//                     // Also create the reverse complemented edge.
//                     alignmentOrientedReadId0.flipStrand();
//                     alignmentOrientedReadId1.flipStrand();
//                     add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


//                 }

//             }

//             // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
//             // it means that we have connected them with alignments.
//             if(success) {
//                 // cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
//             }


//         }

//         if(!success) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No alignments added" << endl;
//         }

            
//     }


//     // Print how many alignments were kept
//     const long keepCountR4 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Adding alignments for break bridging to connect endNodes in different disjointSets: Keeping " << keepCountR4 << " alignments of " << keepAlignment.size() << endl;



//     // Verify that for any read the two oriented reads are in distinct
//     // connected components.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         const OrientedReadId orientedReadId0(readId, 0);
//         const OrientedReadId orientedReadId1(readId, 1);
//         SHASTA_ASSERT(
//             disjointSets.find_set(orientedReadId0.getValue()) !=
//             disjointSets.find_set(orientedReadId1.getValue())
//         );
//     }




//     //*
//     //
//     // Create the read graph using the alignments we selected.
//     //
//     //*
//     createReadGraphUsingSelectedAlignments(keepAlignment);


//     // Gather the vertices of each component.
//     std::map<ReadId, vector<OrientedReadId> > componentMap;

//     for(ReadId readId=0; readId<readCount; readId++) {
//         for(Strand strand=0; strand<2; strand++) {
//             const OrientedReadId orientedReadId(readId, strand);
//             const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
//             componentMap[componentId].push_back(orientedReadId);
//         }
//     }
    
//     cout << "The read graph has " << componentMap.size() << " connected components." << endl;

    

//     cout << timestamp << "Done processing alignments." << endl;


    

//     cout << timestamp << "createReadGraph4 with strand separation ends." << endl;

//     cout << "Strand separation flagged " << crossStrandEdgeCount <<
//         " read graph edges out of " << num_edges(readGraph) << " total in round 1." << endl;
    
//     // cout << "Strand separation flagged " << crossStrandEdgeCountR2 <<
//     //     " read graph edges out of " << readGraph.edges.size() << " total in round 2." << endl;


//     // Verify that for any read the two oriented reads are in distinct
//     // connected components.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         const OrientedReadId orientedReadId0(readId, 0);
//         const OrientedReadId orientedReadId1(readId, 1);
//         SHASTA_ASSERT(
//             disjointSets.find_set(orientedReadId0.getValue()) !=
//             disjointSets.find_set(orientedReadId1.getValue())
//         );
//     }



//     // Sort the components by decreasing size (number of reads).
//     // componentTable contains pairs(size, componentId as key in componentMap).
//     vector< pair<uint64_t, uint64_t> > componentTable;
//     for(const auto& p: componentMap) {
//         const vector<OrientedReadId>& component = p.second;
//         componentTable.push_back(make_pair(component.size(), p.first));
//     }
//     sort(componentTable.begin(), componentTable.end(), std::greater<pair<uint64_t, uint64_t>>());



//     // Store components in this order of decreasing size.
//     vector< vector<OrientedReadId> > components;
//     for(const auto& p: componentTable) {
//         components.push_back(componentMap[ReadId(p.second)]);
//     }
//     performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



//     // Write information for each component.
//     ofstream csv("ReadGraphComponents.csv");
//     csv << "Component,RepresentingRead,OrientedReadCount,"
//         "AccumulatedOrientedReadCount,"
//         "AccumulatedOrientedReadCountFraction\n";
//     uint64_t accumulatedOrientedReadCount = 0;
//     for(ReadId componentId=0; componentId<components.size(); componentId++) {
//         const vector<OrientedReadId>& component = components[componentId];

//         // Stop writing when we reach connected components
//         // consisting of a single isolated read.
//         if(component.size() == 1) {
//             break;
//         }

//         accumulatedOrientedReadCount += component.size();
//         const double accumulatedOrientedReadCountFraction =
//             double(accumulatedOrientedReadCount)/double(orientedReadCount);

//         // The above process of strand separation should have removed
//         // all self-complementary components.
//         const bool isSelfComplementary =
//             component.size() > 1 &&
//             (component[0].getReadId() == component[1].getReadId());
//         SHASTA_ASSERT(not isSelfComplementary);


//         // Write out.
//         csv << componentId << ",";
//         csv << component.front() << ",";
//         csv << component.size() << ",";
//         csv << accumulatedOrientedReadCount << ",";
//         csv << accumulatedOrientedReadCountFraction << "\n";
//     }



//     // For Mode 2 and Mode 3 assembly, we will only assemble one connected component
//     // of each pair. In each pair, we choose the component in the pair
//     // that has the lowest numbered read on strand 0.
//     // Then, for each read we store in its ReadFlags the strand
//     // that the read appears in in this component.
//     // That flag will be used in Mode 2 assembly to
//     // select portions of the marker graph that should be assembled.
//     uint64_t n = 0;
//     for(ReadId componentId=0; componentId<components.size(); componentId++) {
//         const vector<OrientedReadId>& component = components[componentId];

//         // If the lowest numbered read is on strand 1, this is not one of
//         // the connected components we want to use.
//         if(component.front().getStrand() == 1) {
//             continue;
//         }

//         // Store the strand for each read in this component.
//         for(const OrientedReadId orientedReadId: component) {
//             reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
//         }
//         n += component.size();
//     }
//     SHASTA_ASSERT(n == readCount);

// }




























































// //
// //
// //
// // UNCOMMENT THIS TO MAKE IT WORK
// //
// //
// //



// //class AlignmentStats{public: double errorRateRle; uint64_t alignedRange; uint64_t rightUnaligned; uint64_t leftUnaligned; uint64_t alignmentId;};


// void Assembler::createReadGraph4withStrandSeparation(
//     uint64_t maxAlignmentCount,
//     double epsilon,
//     double delta,
//     double WThreshold,
//     double WThresholdForBreaks
//     )
// {
//     cout << timestamp << "createReadGraph4 with strand separation begins" << endl;

//     // Get the total number of stored alignments.
//     const uint64_t alignmentCount = alignmentData.size();
//     SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);


//     //*
//     //
//     // Order alignments in order of increasing Q. 
//     //
//     // Gather in alignmentTablePassFilter[alignmentID, Q]
//     // alignments in order of increasing Q.
//     // Q(n) = (1 + δ/2ε)^n * e-δL
//     // ε = 1e-4, δ = 5e-4
//     // logQ(n) = αn - δL
//     //
//     //*

//     // const double epsilon = 1e-4;
//     // const double delta = 5e-4;
//     const double alpha = log(1 + delta/(2*epsilon));

//     // const double WThreshold = 1e-8;
//     const double logWThreshold = log(WThreshold);

//     // const double WThresholdForBreaks = 1e+15;
//     const double logWThresholdForBreaks = log(WThresholdForBreaks);



//     vector< pair<uint64_t, double> > alignmentTablePassFilter;
//     vector< pair<uint64_t, double> > alignmentTableNotPassFilter;
    

//     // Get stats about the reads
//     const uint64_t readCount = reads->readCount();
//     const uint64_t orientedReadCount = 2*readCount;
    
//     // Keep track of which readIds were used in alignments
//     vector<bool> readIDsPassFilter(readCount, false);

//     // Flag alignments to be kept for break detection.
//     vector<bool> keepAlignmentsForBreaks(alignmentCount, false);


//     // Loop over all alignments.
//     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         if((alignmentId % 100000) == 0) {
//             cout << timestamp << alignmentId << "/" << alignmentCount << endl;
//         }

//         // Get information for this alignment.
//         AlignmentData& thisAlignmentData = alignmentData[alignmentId];

//         // The alignment is stored as an alignment between readId0 on strand 0
//         // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         // The reverse complement alignment also exists, but is not stored explicitly.
//         const ReadId readId0 = thisAlignmentData.readIds[0];
//         const ReadId readId1 = thisAlignmentData.readIds[1];
//         const bool isSameStrand = thisAlignmentData.isSameStrand;
//         SHASTA_ASSERT(readId0 < readId1);
//         const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
//         const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

//         // Store this pair of edges in our edgeTable.
//         const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
//         const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
//         const double L = double(range0 + range1)/2.;
//         // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
//         const double errorRateRle = thisAlignmentData.info.errorRateRle;
//         const double nRLE = errorRateRle * 2 * L;
//         // const double markerCount = thisAlignmentData.info.markerCount;

//         // logQ(n) = αn - δL
//         const double logQ = alpha * double(nRLE) - delta * L;

//         // const uint64_t thisAlignmentMarkerCount = thisAlignmentData.info.markerCount;
        
//         //if (logQ <= logWThreshold)
//         if (logQ <= logWThreshold) {
//             alignmentTablePassFilter.push_back(make_pair(alignmentId, logQ));
//             // alignmentTablePassFilter.push_back(make_pair(alignmentId, thisAlignmentMarkerCount));
//             readIDsPassFilter[readId0] = true;
//             readIDsPassFilter[readId1] = true;
//         } else if(logQ <= logWThresholdForBreaks){
//             alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
//             // alignmentTableNotPassFilter.push_back(make_pair(alignmentId, thisAlignmentMarkerCount));
//             keepAlignmentsForBreaks[alignmentId] = true;
//         }

//     }

    

//     sort(alignmentTablePassFilter.begin(), alignmentTablePassFilter.end(), OrderPairsBySecondOnly<uint64_t, double>());
//     sort(alignmentTableNotPassFilter.begin(), alignmentTableNotPassFilter.end(), OrderPairsBySecondOnly<uint64_t, double>());
//     cout << "The alignmentTable has " << alignmentTablePassFilter.size() << " entries." << endl;
//     cout << "The alignmentTableNotPassFilter has " << alignmentTableNotPassFilter.size() << " entries." << endl;
    





//     ///
//     // Process alignments in order of increasing Q. 
//     //
//     // i.   Start with no edges in the read graph. 
//     // ii.  Process alignments in order of increasing Q. 
//     // iii. If the alignment breaks strand separation, it is skipped. 
//     // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped.  
//     // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.
    
//     // Maintain a vector containing the degree of each vertex
//     // verticesDegree[vertexID] -> degree
//     vector<uint64_t> verticesDegree(orientedReadCount, 0);

//     vector<ReadId> rank(orientedReadCount);
//     vector<ReadId> parent(orientedReadCount);
//     boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
//     for(ReadId readId=0; readId<readCount; readId++) {
//         for(Strand strand=0; strand<2; strand++) {
//             disjointSets.make_set(OrientedReadId(readId, strand).getValue());
//         }
//     }

//     cout << "Number of reads: " << readCount << endl;
//     cout << "Number of oriented reads: " << orientedReadCount << endl;

//     // Flag all alignments as not to be kept.
//     vector<bool> keepAlignment(alignmentCount, false);

//     // Process alignments in order of increasing Q
//     vector alignmentTablesToProcess({alignmentTablePassFilter});
    
//     uint64_t crossStrandEdgeCount = 0;

//     for (auto alignmentTableToProcess : alignmentTablesToProcess) {
//         for(auto it=alignmentTableToProcess.begin(); it!=alignmentTableToProcess.end(); ++it) {
//             const pair<uint64_t, double>& p = *it;
//             const uint64_t alignmentId = p.first;
//             // const double logQ = p.second;

//             // Get the alignment data
//             AlignmentData& alignment = alignmentData[alignmentId];
//             const ReadId readId0 = alignment.readIds[0];
//             const ReadId readId1 = alignment.readIds[1];
//             const bool isSameStrand = alignment.isSameStrand;
//             SHASTA_ASSERT(readId0 < readId1);
//             const OrientedReadId A0 = OrientedReadId(readId0, 0);
//             const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
//             const OrientedReadId A1 = OrientedReadId(readId0, 1);
//             const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             // If the alignment breaks strand separation, it is skipped.
//             // If A0 and B1 are in the same connected component,
//             // A1 and B0 also must be in the same connected component.
//             // Adding this pair of edges would create a self-complementary
//             // connected component containing A0, B0, A1, and B1,
//             // and to ensure strand separation we don't want to do that.
//             // So we mark these edges as cross-strand edges
//             // and don't use them to update the disjoint set data structure.
//             if(a0 == b1) {
//                 SHASTA_ASSERT(a1 == b0);
//                 crossStrandEdgeCount += 2;
//                 continue;
//             }

            

//             // // If both vertices of the potential edge have at least the required minimum number 
//             // // of neighbors, the alignment is also skipped. 
//             // const uint64_t degreeA0 = verticesDegree[A0.getValue()];
//             // const uint64_t degreeB0 = verticesDegree[B0.getValue()];



//             // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
//             //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
//             //     continue;
//             // }

//             // Add the alignment to the read graph.
//             keepAlignment[alignmentId] = true;
//             alignment.info.isInReadGraph = 1;

//             // Update vertex degrees
//             verticesDegree[A0.getValue()]++;
//             verticesDegree[B0.getValue()]++;
//             verticesDegree[A1.getValue()]++;
//             verticesDegree[B1.getValue()]++;
            

//             // Update disjoint sets
//             disjointSets.union_set(a0, b0);
//             disjointSets.union_set(a1, b1);

//         }

//         // Verify that for any read the two oriented reads are in distinct
//         // connected components.
//         for(ReadId readId=0; readId<readCount; readId++) {
//             const OrientedReadId orientedReadId0(readId, 0);
//             const OrientedReadId orientedReadId1(readId, 1);
//             SHASTA_ASSERT(
//                 disjointSets.find_set(orientedReadId0.getValue()) !=
//                 disjointSets.find_set(orientedReadId1.getValue())
//             );
//         }
//     }

//     // Track the size of each set in the disjoint sets
//     vector<uint64_t> setSizes(orientedReadCount, 0);
//     for (std::uint64_t i = 0; i < orientedReadCount; ++i) {
//         setSizes[disjointSets.find_set(i)]++;
//     }


//     // Print how many alignments were kept in this step
//     const long keepCountR1 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Finding strict disjointSets step: Keeping " << keepCountR1 << " alignments of " << keepAlignment.size() << endl;





//     //*
//     //
//     // Create the dynamically adjustable boost readGraph using the alignments we selected.
//     //
//     //*
//     using boost::add_vertex;
//     using boost::add_edge;

//     // The vertex_descriptor is OrientedReadId::getValue().
//     ReadGraph4 readGraph(orientedReadCount);

//     // Initially, each alignment generates two edges.
//     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//         // Record whether this alignment is used in the read graph.
//         const bool keepThisAlignment = keepAlignment[alignmentId];
//         const AlignmentData& alignment = alignmentData[alignmentId];

//         // If this alignment is not used in the read graph, we are done.
//         if(not keepThisAlignment) {
//             continue;
//         }

//         // Get the OrientedReadIds.
//         OrientedReadId orientedReadId0(alignment.readIds[0], 0);
//         OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//         SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

//         // Swap them if necessary, depending on the average alignment offset at center.
//         if(alignment.info.offsetAtCenter() < 0.) {
//             swap(orientedReadId0, orientedReadId1);
//         }

//         // Create the edge.
//         add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//         // Also create the reverse complemented edge.
//         orientedReadId0.flipStrand();
//         orientedReadId1.flipStrand();
//         add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
//     }

//     cout << "The read graph has " << num_vertices(readGraph) << " vertices and " << num_edges(readGraph) << " edges." << endl;



//     //*
//     //
//     // Create the dynamically adjustable boost readGraph using the alignments that did not pass the strict filter.
//     // These alignments are used to create the read graph that will aid in the break detection.
//     //
//     //*
    
//     // The vertex_descriptor is OrientedReadId::getValue().
//     ReadGraph4AllAlignments readGraphAllAlignments(orientedReadCount);

//     // Initially, each alignment generates two edges.
//     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//         // Record whether this alignment is used in the read graph.
//         const bool keepThisAlignment = keepAlignmentsForBreaks[alignmentId];
//         const AlignmentData& alignment = alignmentData[alignmentId];

//         // If this alignment is not used in the read graph, we are done.
//         if(not keepThisAlignment) {
//             continue;
//         }

//         // Get the OrientedReadIds.
//         OrientedReadId orientedReadId0(alignment.readIds[0], 0);
//         OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//         SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

//         // Swap them if necessary, depending on the average alignment offset at center.
//         if(alignment.info.offsetAtCenter() < 0.) {
//             swap(orientedReadId0, orientedReadId1);
//         }

//         // Create the edge.
//         add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

//         // Also create the reverse complemented edge.
//         orientedReadId0.flipStrand();
//         orientedReadId1.flipStrand();
//         add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
//     }
    
//     cout << "The read graph for break detection has " << num_vertices(readGraphAllAlignments) << " vertices and " << num_edges(readGraphAllAlignments) << " edges." << endl;





//     //*
//     //
//     // Find possible dead end nodes on the directed graph.
//     // We only keep the reads that have no outgoing neighbors.
//     //
//     //*
//     vector<bool> potentialDeadEndReads(orientedReadCount, false);
//     vector<bool> isolatedReads(orientedReadCount, false);

//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {
//             OrientedReadId orientedReadId(readId, strand);
            
//             // Find neighbors in the forward direction
//             vector<OrientedReadId> forwardNeighbors;
//             readGraph.findNeighborsDirectedGraphOneSideRight(orientedReadId, 1, forwardNeighbors);
            
//             // Find neighbors in the backward direction 
//             vector<OrientedReadId> leftNeighbors;
//             readGraph.findNeighborsDirectedGraphOneSideLeft(orientedReadId, 1, leftNeighbors);

//             if (forwardNeighbors.empty() && leftNeighbors.empty() ) {
//                 isolatedReads[orientedReadId.getValue()] = true;
//                 OrientedReadId reverseOrientedReadId = orientedReadId;
//                 reverseOrientedReadId.flipStrand();
//                 isolatedReads[reverseOrientedReadId.getValue()] = true;
//                 continue;
//             }

//             // If a read has neighbors only in the backward direction, it's a potential dead end
//             if (forwardNeighbors.empty() && !leftNeighbors.empty()) {
//                 potentialDeadEndReads[orientedReadId.getValue()] = true;
//             }
//         }
//     }

//     // count the number of potential dead end reads
//     long potentialDeadEndReadCount = count(potentialDeadEndReads.begin(), potentialDeadEndReads.end(), true);
//     cout << "Found " << potentialDeadEndReadCount << " potential dead end reads." << endl;

//     // // Print dead end Oriented reads
//     // // iterate over all oriented reads
//     // for (ReadId readId = 0; readId < readCount; readId++) {
//     //     for (Strand strand = 0; strand < 2; strand++) {
//     //         OrientedReadId orientedReadId(readId, strand);
//     //         if (potentialDeadEndReads[orientedReadId.getValue()]) {
//     //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a potential dead end read." << endl;
//     //         }
//     //     }
//     // }






//     //*
//     //
//     // Filter out the potential dead end nodes on the directed graph.
//     // Look for a read path that lead to a positive offset.
//     // If no such path is found, the orientedReadId is kept as a potential dead end node.
//     //
//     //*
//     uint64_t finalNumberOfPotentialDeadEndNodes = 0;
//     vector<bool> finalDeadEndReadsWithNoOutgoingNodes(orientedReadCount, false);
//     vector<bool> finalDeadEndReadsWithNoIncomingNodes(orientedReadCount, false);
//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {

//             OrientedReadId orientedReadId(readId, strand);

//             OrientedReadId reverseOrientedReadId = orientedReadId;
//             reverseOrientedReadId.flipStrand();

//             // check if the orientedReadId is in potentialDeadEndReads
//             if(!potentialDeadEndReads[orientedReadId.getValue()]) {
//                 continue;
//             }

//             // create the necessary variables for findAllPaths
//             vector<vector<OrientedReadId>> paths;
//             vector<vector<double>> pathsOffsets;
//             vector<OrientedReadId> currentPath;
//             vector<double> currentPathOffset;
//             std::set<ReadGraph4BaseClass::vertex_descriptor> visited;
//             uint64_t maxDistance = 4;
//             uint64_t currentDistance = 0;

//             bool result = readGraph.findPathWithPositiveOffset(orientedReadId, paths, pathsOffsets, currentPath, currentPathOffset, visited, maxDistance, currentDistance + 1, alignmentData, readGraph);

//             // Check if we found a read path with positive offset.
//             // If yes, the function findPathWithPositiveOffset will return 1, if not, it will return 0.
//             if(result == 1) {
//                 // cout << "Found a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset" << endl;
//                 // //print the paths and then the pathsOffsets
//                 // for (uint64_t i = 0; i < paths.size(); i++) {
//                 //     cout << "Path " << i << ": ";
//                 //     for (uint64_t j = 0; j < paths[i].size(); j++) {
//                 //         cout << paths[i][j].getReadId() << " ";
//                 //     }
//                 //     cout << endl;
//                 //     cout << "PathOffsets " << i << ": ";
//                 //     for (uint64_t j = 0; j < pathsOffsets[i].size(); j++) {
//                 //         cout << pathsOffsets[i][j] << " ";
//                 //     }
//                 //     cout << endl;
//                 // }
//             } else if (result == 0) {
//                 finalNumberOfPotentialDeadEndNodes++;
//                 // cout << "Did not find a path for the orientedRead with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " with positive offset. Keeping it as a potential dead end read." << endl;
//                 finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()] = true;
//                 finalDeadEndReadsWithNoIncomingNodes[reverseOrientedReadId.getValue()] = true;
//             }
//         }
//     }

//     cout << "After filtering we are left with " << finalNumberOfPotentialDeadEndNodes << " potential dead end reads." << endl;


//     // // print dead end Oriented reads
//     // // iterate over all oriented reads
//     // for (ReadId readId = 0; readId < readCount; readId++) {
//     //     for (Strand strand = 0; strand < 2; strand++) {
//     //         OrientedReadId orientedReadId(readId, strand);
//     //         if (finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]) {
//     //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
//     //         } else if (finalDeadEndReadsWithNoIncomingNodes[orientedReadId.getValue()]) {
//     //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
//     //         }
//     //     }
//     // }











//     // //*
//     // //
//     // // Extend the potential dead end nodes list.
//     // // Add neighboring nodes of potential dead end nodes to the dead end node list.
//     // //
//     // //*
//     // vector<bool> finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors(orientedReadCount, false);
//     // vector<bool> finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors(orientedReadCount, false);
//     // for (ReadId readId = 0; readId < readCount; readId++) {
//     //     for (Strand strand = 0; strand < 2; strand++) {
//     //         OrientedReadId orientedReadId(readId, strand);

//     //         if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){

//     //             finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()] = true;
                
//     //             OrientedReadId reverseOrientedReadId = orientedReadId;
//     //             reverseOrientedReadId.flipStrand();

//     //             finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;

//     //             // Find distance 1 neighbors
//     //             vector<OrientedReadId> distance1Neighbors;
//     //             readGraph.findNeighborsDirectedGraphBothSides(orientedReadId, 1, distance1Neighbors);

//     //             for(auto neighbor : distance1Neighbors) {
//     //                 finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[neighbor.getValue()] = true;

//     //                 OrientedReadId reverseOrientedReadId = neighbor;
//     //                 reverseOrientedReadId.flipStrand();

//     //                 finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[reverseOrientedReadId.getValue()] = true;
//     //             }
//     //         }
//     //     }
//     // }

//     // // // print dead end Oriented reads
//     // // // iterate over all oriented reads
//     // // for (ReadId readId = 0; readId < readCount; readId++) {
//     // //     for (Strand strand = 0; strand < 2; strand++) {
//     // //         OrientedReadId orientedReadId(readId, strand);
//     // //         if (finalDeadEndReadsWithNoOutgoingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
//     // //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no outgoing nodes." << endl;
//     // //         } else if (finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors[orientedReadId.getValue()]) {
//     // //             cout << "OrientedReadId with ReadID " << orientedReadId.getReadId() << " and strand " << orientedReadId.getStrand() << " is a final dead end read with no incoming nodes." << endl;
//     // //         }
//     // //     }
//     // // }





//     //*
//     //
//     // Create a map of endNodesWithNoOutgoingNodes to other endNodesWithNoIncomingNodes that are possible to connect to.
//     //
//     //*

//     // Case 1: the NoOut deadEnd node will be mapped to at least 3 NoIn deadEnd nodes 
//     // IN THE SAME disjointSet.
//     //
//     //  NoIn - - - - | - - - - -  - - - - - - | - - - - -
//     //  NoIn - - - - | - - NoOut     NoIn - - | - - - - -

//     std::unordered_map<uint64_t, vector<uint64_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet;
//     vector<bool> endNodesWithNoOutgoingNodesConsidered(orientedReadCount, false);

//     // First, we try to find EndNodesWithNoIncomingNodes that are in the same connected component
//     // as EndNodesWithNoOutgoingNodes. These have priority over other EndNodesWithNoIncomingNodes in other connected components.
//     // These will also include telomeric nodes because they do not have incoming nodes!
//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {
//             OrientedReadId orientedReadId(readId, strand);

//             // check if the orientedReadId is an endNode with no outgoing nodes
//             if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
//                 // Get it's disjointSet
//                 uint64_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

//                 // Find all EndNodesWithNoIncomingNodes that are in the same disjointSet
//                 for(uint64_t id=0; id<orientedReadCount; id++) {
//                     // check if the id is an endNode with no incoming nodes
//                     if(finalDeadEndReadsWithNoIncomingNodes[id]){
//                         // get the disjointSet of the id
//                         uint64_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
//                         // check if the disjointSet of the id is the same as the disjointSet of the endNode with no outgoing nodes
//                         if(deadEndReadWithNoOutgoingNodesDisjointSetId == endNodeWithNoIncomingNodesDisjointSetId){
//                             endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet[orientedReadId.getValue()].push_back(id);
//                             endNodesWithNoOutgoingNodesConsidered[orientedReadId.getValue()] = true;
//                         }
//                     }
//                 }

//             }
//         }
//     }


//     for(auto& p : endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesSameDisjointSet) {
//         uint64_t value = p.first;
//         OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(value));
//         SHASTA_ASSERT(orientedReadId.getValue() == value);
//         // const ReadId readId = orientedReadId.getReadId();
//         // const Strand strand = orientedReadId.getStrand();

//         // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in the same disjointSet:" << endl;
//         vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
//         for(auto& node : p.second) {
//             OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(ReadId(node));
//             // cout << "ReadID " << nodeOrientedReadId.getReadId() << " strand " << nodeOrientedReadId.getStrand() << endl;
//             SHASTA_ASSERT(nodeOrientedReadId.getValue() == node);
//             endNodesWithNoIncomingNodes[node] = true;
//         }

        
//         // Find neighbors in the forward direction of the ALL ALIGNMENTS read graph 
//         // starting from orientedReadId which is an endNode with no outgoing nodes in the filtered read graph. 
//         // Early stop when we reach an endNode with no incoming nodes (a node with endNodesWithNoIncomingNodes set to true).
//         vector<OrientedReadId> forwardNeighbors;
//         readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, endNodesWithNoIncomingNodes, 5, forwardNeighbors);

//         // // print the forward neighbors
//         // cout << "Forward neighbors of the endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " are:" << endl;
//         // for(auto& neighbor : forwardNeighbors) {
//         //     cout << "ReadID " << neighbor.getReadId() << " strand " << neighbor.getStrand() << endl;
//         // }


//         // create a std::set of the forwardNeighbors for easy contain check
//         std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
        

//         if(forwardNeighbors.empty()) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No forward neighbors found" << endl;
//             continue;
//         }


//         // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
//         // OR a non speficic node if we exceeded maxDistance
//         OrientedReadId lastNode = forwardNeighbors.back();


//         bool success = false;


//         // First check if the last node is a potential dead end node with no INCOMING nodes
//         if(endNodesWithNoIncomingNodes[lastNode.getValue()]) {

//             // cout << "The last node is: " << lastNode.getReadId() << " strand " << lastNode.getStrand() << endl;
            
//             const OrientedReadId A0 = orientedReadId;
//             const OrientedReadId B0 = lastNode;
//             const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
//             const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             // const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             // const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             // const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             // const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
//             // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//                 const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
//                 const uint64_t alignmentId = p.first;
//                 // const double logQ = p.second;

//                 const bool keepThisAlignment = keepAlignment[alignmentId];

//                 const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

//                 if(keepThisAlignment) {
//                     continue;
//                 }

//                 if(not keepThisBreaksAlignment) {
//                     continue;
//                 }

//                 AlignmentData& alignment = alignmentData[alignmentId];
            
//                 // Get the OrientedReadIds.
//                 OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
//                 OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//                 SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

//                 // Swap them if necessary, depending on the average alignment offset at center.
//                 if(alignment.info.offsetAtCenter() < 0.) {
//                     swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
//                 }


//                 if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                    
//                     // Get the alignment data
//                     ReadId readId0v2 = alignment.readIds[0];
//                     ReadId readId1v2 = alignment.readIds[1];
//                     const bool isSameStrandv2 = alignment.isSameStrand;
//                     SHASTA_ASSERT(readId0v2 < readId1v2);
//                     OrientedReadId A0v2 = OrientedReadId(readId0v2, 0);
//                     OrientedReadId B0v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 0 : 1);
//                     OrientedReadId A1v2 = OrientedReadId(readId0v2, 1);
//                     OrientedReadId B1v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 1 : 0);

//                     SHASTA_ASSERT(A0v2.getReadId() == A1v2.getReadId());
//                     SHASTA_ASSERT(B0v2.getReadId() == B1v2.getReadId());
//                     SHASTA_ASSERT(A0v2.getStrand() == 1 - A1v2.getStrand());
//                     SHASTA_ASSERT(B0v2.getStrand() == 1 - B1v2.getStrand());

//                     // Get the connected components that these oriented reads are in.
//                     const uint64_t a0v2 = disjointSets.find_set(A0v2.getValue());
//                     const uint64_t b0v2 = disjointSets.find_set(B0v2.getValue());
//                     const uint64_t a1v2 = disjointSets.find_set(A1v2.getValue());
//                     const uint64_t b1v2 = disjointSets.find_set(B1v2.getValue());


//                     // If the alignment breaks strand separation, it is skipped.
//                     // If A0 and B1 are in the same connected component,
//                     // A1 and B0 also must be in the same connected component.
//                     // Adding this pair of edges would create a self-complementary
//                     // connected component containing A0, B0, A1, and B1,
//                     // and to ensure strand separation we don't want to do that.
//                     // So we mark these edges as cross-strand edges
//                     // and don't use them to update the disjoint set data structure.
//                     if(a0v2 == b1v2) {
//                         SHASTA_ASSERT(a1v2 == b0v2);
//                         crossStrandEdgeCount += 2;
//                         continue;
//                     }

//                     // Add the alignment to the read graph.
//                     keepAlignment[alignmentId] = true;
//                     alignment.info.isInReadGraph = 1;

//                     // Update vertex degrees
//                     verticesDegree[A0v2.getValue()]++;
//                     verticesDegree[B0v2.getValue()]++;
//                     verticesDegree[A1v2.getValue()]++;
//                     verticesDegree[B1v2.getValue()]++;
                    

//                     // Update disjoint sets
//                     disjointSets.union_set(a0v2, b0v2);
//                     disjointSets.union_set(a1v2, b1v2);

//                     // Make sure all alignments added are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1v2.getValue()] = false;

//                     // Make sure the start and last nodes are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1.getValue()] = false;

//                     success = true;

//                     // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

//                     // Create the edge.
//                     add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//                     // Also create the reverse complemented edge.
//                     alignmentOrientedReadId0.flipStrand();
//                     alignmentOrientedReadId1.flipStrand();
//                     add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


//                 }

//             }

//             // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
//             // it means that we have connected them with alignments.
//             if(success) {
//                 // cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
//             }

//             // // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
//             // // it means that we have connected them with alignments.
//             // if(finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] == false and finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] == false) {
//             //     cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
//             // } else {
//             //     cout << "Did not connect endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // }

//         }

//         if(!success) {
//            // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//            // cout << "No alignments added" << endl;
//         }

            
//     }


//     // Print how many alignments were kept
//     const long keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Adding alignments for break bridging to connect endNodes in the same disjointSet: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;




    
//     // Case 2: the NoOut deadEnd node will be mapped to at least 2 NoIn deadEnd nodes in the same disjointSet.
//     // Other NoIn deadEnd nodes will be in a different disjointSet.
//     //
//     //  NoIn - - - - | - - - - -  - - - - - - - - - - -
//     //  NoIn - - - - | - - NoOut     
//     //                                     NoIn - - - - - - - (Other disjointSet)

//     // Case 3: Breaks in a haploid chromosome (chrX chrY).
//     //
//     //  NoIn - - - - - - NoOut
//     //                             NoIn - - - - - - - (Other disjointSet) 


//     std::unordered_map<uint64_t, vector<uint64_t>> endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet;

//     // Now, we try to find EndNodesWithNoIncomingNodes that are in different disjointSets than EndNodesWithNoOutgoingNodes.
//     for (ReadId readId = 0; readId < readCount; readId++) {
//         for (Strand strand = 0; strand < 2; strand++) {
//             OrientedReadId orientedReadId(readId, strand);

//             // check if the orientedReadId is an endNode with no outgoing nodes
//             if(finalDeadEndReadsWithNoOutgoingNodes[orientedReadId.getValue()]){
//                 // Get it's disjointSet
//                 uint64_t deadEndReadWithNoOutgoingNodesDisjointSetId = disjointSets.find_set(orientedReadId.getValue());

//                 // Find all EndNodesWithNoIncomingNodes that are in different disjointSet
//                 for(uint64_t id=0; id<orientedReadCount; id++) {
//                     // check if the id is an endNode with no incoming nodes
//                     if(finalDeadEndReadsWithNoIncomingNodes[id]){
//                         // get the disjointSet of the id
//                         uint64_t endNodeWithNoIncomingNodesDisjointSetId = disjointSets.find_set(id);
//                         // check if the disjointSet of the id is different from the disjointSet of the endNode with no outgoing nodes
//                         if(deadEndReadWithNoOutgoingNodesDisjointSetId != endNodeWithNoIncomingNodesDisjointSetId){
//                             endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet[orientedReadId.getValue()].push_back(id);
//                             endNodesWithNoOutgoingNodesConsidered[orientedReadId.getValue()] = true;
//                         }
//                     }
//                 }

//             }
//         }
//     }


//     for(auto& p : endNodesWithNoOutgoingNodesToEndNodesWithNoIncomingNodesDifferentDisjointSet) {
//         uint64_t value = p.first;
//         OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(value));
//         SHASTA_ASSERT(orientedReadId.getValue() == value);
//         // const ReadId readId = orientedReadId.getReadId();
//         // const Strand strand = orientedReadId.getStrand();

//         // cout << "EndNodeWithNoOutgoingNodes ReadID " << readId << " and strand " << strand << " is mapped to these NoIn deadEnd nodes in a different disjointSet: " << endl;
//         vector<bool> endNodesWithNoIncomingNodes(orientedReadCount, false);
//         for(auto& node : p.second) {
//             OrientedReadId nodeOrientedReadId = OrientedReadId::fromValue(ReadId(node));
//             // cout << "ReadID " << nodeOrientedReadId.getReadId() << " strand " << nodeOrientedReadId.getStrand() << endl;
//             SHASTA_ASSERT(nodeOrientedReadId.getValue() == node);
//             endNodesWithNoIncomingNodes[node] = true;
//         }
        
//         // Find neighbors in the forward direction of the ALL ALIGNMENTS read graph 
//         // starting from orientedReadId which is an endNode with no outgoing nodes in the filtered read graph. 
//         // Early stop when we reach an endNode with no incoming nodes (a node with endNodesWithNoIncomingNodes set to true).
//         vector<OrientedReadId> forwardNeighbors;
//         readGraphAllAlignments.findNeighborsEarlyStopWhenReachEndNode(orientedReadId, endNodesWithNoIncomingNodes, 5, forwardNeighbors);

//         // create a std::set of the forwardNeighbors for easy contain check
//         std::set<OrientedReadId> forwardNeighborsSet(forwardNeighbors.begin(), forwardNeighbors.end());
        

//         if(forwardNeighbors.empty()) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No forward neighbors found" << endl;
//             continue;
//         }

//         // // print the forward neighbors
//         // cout << "Forward neighbors of the endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " are:" << endl;
//         // for(auto& neighbor : forwardNeighbors) {
//         //     cout << "ReadID " << neighbor.getReadId() << " strand " << neighbor.getStrand() << endl;
//         // }

//         // Get the last item from forwardNeighbors. It contains the first encountered dead end node with no INCOMING nodes
//         // OR a non speficic node if we exceeded maxDistance
//         OrientedReadId lastNode = forwardNeighbors.back();

//         bool success = false;

//         // First check if the last node is a potential dead end node with no INCOMING nodes
//         if(endNodesWithNoIncomingNodes[lastNode.getValue()]) {
            
//             const OrientedReadId A0 = orientedReadId;
//             const OrientedReadId B0 = lastNode;
//             const OrientedReadId A1 = OrientedReadId(A0.getReadId(), A0.getStrand() == 0 ? 1 : 0);
//             const OrientedReadId B1 = OrientedReadId(B0.getReadId(), B0.getStrand() == 0 ? 1 : 0);

//             SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
//             SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
//             SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
//             SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

//             // Get the connected components that these oriented reads are in.
//             // const uint64_t a0 = disjointSets.find_set(A0.getValue());
//             // const uint64_t b0 = disjointSets.find_set(B0.getValue());
//             // const uint64_t a1 = disjointSets.find_set(A1.getValue());
//             // const uint64_t b1 = disjointSets.find_set(B1.getValue());


//             for(uint64_t index=0; index<alignmentTableNotPassFilter.size(); index++) {
//             // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

//                 const pair<uint64_t, double>& p = alignmentTableNotPassFilter[index];
//                 const uint64_t alignmentId = p.first;
//                 // const double logQ = p.second;

//                 const bool keepThisAlignment = keepAlignment[alignmentId];

//                 const bool keepThisBreaksAlignment = keepAlignmentsForBreaks[alignmentId];

//                 if(keepThisAlignment) {
//                     continue;
//                 }

//                 if(not keepThisBreaksAlignment) {
//                     continue;
//                 }

//                 AlignmentData& alignment = alignmentData[alignmentId];
            
//                 // Get the OrientedReadIds.
//                 OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
//                 OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
//                 SHASTA_ASSERT(alignmentOrientedReadId0 < alignmentOrientedReadId1);

//                 // Swap them if necessary, depending on the average alignment offset at center.
//                 if(alignment.info.offsetAtCenter() < 0.) {
//                     swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
//                 }


//                 if(alignmentOrientedReadId0.getValue() == orientedReadId.getValue() || forwardNeighborsSet.contains(alignmentOrientedReadId0) || forwardNeighborsSet.contains(alignmentOrientedReadId1) ){
                    
//                     // Get the alignment data
//                     ReadId readId0v2 = alignment.readIds[0];
//                     ReadId readId1v2 = alignment.readIds[1];
//                     const bool isSameStrandv2 = alignment.isSameStrand;
//                     SHASTA_ASSERT(readId0v2 < readId1v2);
//                     OrientedReadId A0v2 = OrientedReadId(readId0v2, 0);
//                     OrientedReadId B0v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 0 : 1);
//                     OrientedReadId A1v2 = OrientedReadId(readId0v2, 1);
//                     OrientedReadId B1v2 = OrientedReadId(readId1v2, isSameStrandv2 ? 1 : 0);

//                     SHASTA_ASSERT(A0v2.getReadId() == A1v2.getReadId());
//                     SHASTA_ASSERT(B0v2.getReadId() == B1v2.getReadId());
//                     SHASTA_ASSERT(A0v2.getStrand() == 1 - A1v2.getStrand());
//                     SHASTA_ASSERT(B0v2.getStrand() == 1 - B1v2.getStrand());

//                     // Get the connected components that these oriented reads are in.
//                     const uint64_t a0v2 = disjointSets.find_set(A0v2.getValue());
//                     const uint64_t b0v2 = disjointSets.find_set(B0v2.getValue());
//                     const uint64_t a1v2 = disjointSets.find_set(A1v2.getValue());
//                     const uint64_t b1v2 = disjointSets.find_set(B1v2.getValue());


//                     // If the alignment breaks strand separation, it is skipped.
//                     // If A0 and B1 are in the same connected component,
//                     // A1 and B0 also must be in the same connected component.
//                     // Adding this pair of edges would create a self-complementary
//                     // connected component containing A0, B0, A1, and B1,
//                     // and to ensure strand separation we don't want to do that.
//                     // So we mark these edges as cross-strand edges
//                     // and don't use them to update the disjoint set data structure.
//                     if(a0v2 == b1v2) {
//                         SHASTA_ASSERT(a1v2 == b0v2);
//                         crossStrandEdgeCount += 2;
//                         continue;
//                     }

//                     // Add the alignment to the read graph.
//                     keepAlignment[alignmentId] = true;
//                     alignment.info.isInReadGraph = 1;

//                     // Update vertex degrees
//                     verticesDegree[A0v2.getValue()]++;
//                     verticesDegree[B0v2.getValue()]++;
//                     verticesDegree[A1v2.getValue()]++;
//                     verticesDegree[B1v2.getValue()]++;
                    

//                     // Update disjoint sets
//                     disjointSets.union_set(a0v2, b0v2);
//                     disjointSets.union_set(a1v2, b1v2);

//                     // Make sure all alignments added are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0v2.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1v2.getValue()] = false;

//                     // Make sure the start and last nodes are not considered as dead ends anymore
//                     finalDeadEndReadsWithNoOutgoingNodes[A0.getValue()] = false;
//                     finalDeadEndReadsWithNoOutgoingNodes[A1.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B0.getValue()] = false;
//                     finalDeadEndReadsWithNoIncomingNodes[B1.getValue()] = false;

//                     success = true;


//                     // cout << "Adding alignment " << alignmentId << " between " << alignmentOrientedReadId0.getReadId() << " and " << alignmentOrientedReadId1.getReadId() << endl;

//                     // Create the edge.
//                     add_edge(alignmentOrientedReadId0.getValue(), alignmentOrientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

//                     // Also create the reverse complemented edge.
//                     alignmentOrientedReadId0.flipStrand();
//                     alignmentOrientedReadId1.flipStrand();
//                     add_edge(alignmentOrientedReadId1.getValue(), alignmentOrientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);


//                 }

//             }

//             // If the endNode with no outgoing nodes and the endNode with no incoming nodes are not set to false,
//             // it means that we have connected them with alignments.
//             if(success) {
//                 // cout << "Connected endNode with no outgoing nodes ReadID " << A0.getReadId() << " strand " << A0.getStrand() << " with endNode with no incoming nodes ReadID " << B0.getReadId() << " strand " << B0.getStrand() << endl;
//             }


//         }

//         if(!success) {
//             // cout << "Did not connect endNode with no outgoing nodes ReadID " << orientedReadId.getReadId() << " strand " << orientedReadId.getStrand() << " with any other endNode with no incoming nodes" << endl;
//             // cout << "No alignments added" << endl;
//         }

            
//     }


//     // Print how many alignments were kept
//     const long keepCountR3 = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Adding alignments for break bridging to connect endNodes in different disjointSets: Keeping " << keepCountR3 << " alignments of " << keepAlignment.size() << endl;



//     // Verify that for any read the two oriented reads are in distinct
//     // connected components.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         const OrientedReadId orientedReadId0(readId, 0);
//         const OrientedReadId orientedReadId1(readId, 1);
//         SHASTA_ASSERT(
//             disjointSets.find_set(orientedReadId0.getValue()) !=
//             disjointSets.find_set(orientedReadId1.getValue())
//         );
//     }




//     //*
//     //
//     // Create the read graph using the alignments we selected.
//     //
//     //*
//     createReadGraphUsingSelectedAlignments(keepAlignment);


//     // Gather the vertices of each component.
//     std::map<ReadId, vector<OrientedReadId> > componentMap;

//     for(ReadId readId=0; readId<readCount; readId++) {
//         for(Strand strand=0; strand<2; strand++) {
//             const OrientedReadId orientedReadId(readId, strand);
//             const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
//             componentMap[componentId].push_back(orientedReadId);
//         }
//     }
    
//     cout << "The read graph has " << componentMap.size() << " connected components." << endl;

    

//     cout << timestamp << "Done processing alignments." << endl;


    

//     cout << timestamp << "createReadGraph4 with strand separation ends." << endl;

//     cout << "Strand separation flagged " << crossStrandEdgeCount <<
//         " read graph edges out of " << num_edges(readGraph) << " total in round 1." << endl;
    
//     // cout << "Strand separation flagged " << crossStrandEdgeCountR2 <<
//     //     " read graph edges out of " << readGraph.edges.size() << " total in round 2." << endl;


//     // Verify that for any read the two oriented reads are in distinct
//     // connected components.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         const OrientedReadId orientedReadId0(readId, 0);
//         const OrientedReadId orientedReadId1(readId, 1);
//         SHASTA_ASSERT(
//             disjointSets.find_set(orientedReadId0.getValue()) !=
//             disjointSets.find_set(orientedReadId1.getValue())
//         );
//     }



//     // Sort the components by decreasing size (number of reads).
//     // componentTable contains pairs(size, componentId as key in componentMap).
//     vector< pair<uint64_t, uint64_t> > componentTable;
//     for(const auto& p: componentMap) {
//         const vector<OrientedReadId>& component = p.second;
//         componentTable.push_back(make_pair(component.size(), p.first));
//     }
//     sort(componentTable.begin(), componentTable.end(), std::greater<pair<uint64_t, uint64_t>>());



//     // Store components in this order of decreasing size.
//     vector< vector<OrientedReadId> > components;
//     for(const auto& p: componentTable) {
//         components.push_back(componentMap[ReadId(p.second)]);
//     }
//     performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



//     // Write information for each component.
//     ofstream csv("ReadGraphComponents.csv");
//     csv << "Component,RepresentingRead,OrientedReadCount,"
//         "AccumulatedOrientedReadCount,"
//         "AccumulatedOrientedReadCountFraction\n";
//     uint64_t accumulatedOrientedReadCount = 0;
//     for(ReadId componentId=0; componentId<components.size(); componentId++) {
//         const vector<OrientedReadId>& component = components[componentId];

//         // Stop writing when we reach connected components
//         // consisting of a single isolated read.
//         if(component.size() == 1) {
//             break;
//         }

//         accumulatedOrientedReadCount += component.size();
//         const double accumulatedOrientedReadCountFraction =
//             double(accumulatedOrientedReadCount)/double(orientedReadCount);

//         // The above process of strand separation should have removed
//         // all self-complementary components.
//         const bool isSelfComplementary =
//             component.size() > 1 &&
//             (component[0].getReadId() == component[1].getReadId());
//         SHASTA_ASSERT(not isSelfComplementary);


//         // Write out.
//         csv << componentId << ",";
//         csv << component.front() << ",";
//         csv << component.size() << ",";
//         csv << accumulatedOrientedReadCount << ",";
//         csv << accumulatedOrientedReadCountFraction << "\n";
//     }



//     // For Mode 2 and Mode 3 assembly, we will only assemble one connected component
//     // of each pair. In each pair, we choose the component in the pair
//     // that has the lowest numbered read on strand 0.
//     // Then, for each read we store in its ReadFlags the strand
//     // that the read appears in in this component.
//     // That flag will be used in Mode 2 assembly to
//     // select portions of the marker graph that should be assembled.
//     uint64_t n = 0;
//     for(ReadId componentId=0; componentId<components.size(); componentId++) {
//         const vector<OrientedReadId>& component = components[componentId];

//         // If the lowest numbered read is on strand 1, this is not one of
//         // the connected components we want to use.
//         if(component.front().getStrand() == 1) {
//             continue;
//         }

//         // Store the strand for each read in this component.
//         for(const OrientedReadId orientedReadId: component) {
//             reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
//         }
//         n += component.size();
//     }
//     SHASTA_ASSERT(n == readCount);


// }






// //
// //
// //
// // UNCOMMENT THIS TO MAKE IT WORK
// //
// //
// //













































































































































































































































































void Assembler::createReadGraph4(
    uint64_t /* maxAlignmentCount */)
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

            const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
            const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
            const double L = double(range0 + range1)/2.;
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

        const long keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
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
            std::array<OrientedReadId, 2> nextEdgeOrientedReadIds = nextEdge.orientedReadIds;
            nextEdgeOrientedReadIds[0].flipStrand();
            nextEdgeOrientedReadIds[1].flipStrand();
            SHASTA_ASSERT(nextEdgeOrientedReadIds == edge.orientedReadIds);
            SHASTA_ASSERT(nextEdge.alignmentId == alignmentId);
        }


        // // Store this pair of edges in our edgeTable.
        // const uint64_t range0 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[0], 0, markers);
        // const uint64_t range1 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[1], 1, markers);
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
    for(uint64_t i = 0; i < min(uint64_t(100), edgeTable.size()); i++) {
        const auto& p = edgeTable[i];
        const uint64_t edgeId = p.first;
        const double q = p.second;
        const AlignmentData& alignment = alignmentData[readGraph.edges[edgeId].alignmentId];
        const uint64_t markerCount = alignment.info.markerCount;
        const uint64_t range0 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[0], 0, markers);
        const uint64_t range1 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[1], 1, markers);
        const double L = double(range0 + range1)/2.;
        const uint64_t n = alignment.info.mismatchCountRle;
        const double errorRateRleNew = double(n)/(2.0*L);;
        const double errorRateRleOld = alignment.info.errorRateRle;
        // const double m = L * 0.0001;
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        cout << edgeId << "(" << readId0 << "," << readId1 << ")" << "\t" << q << "\t" << markerCount << "\t" << errorRateRleNew << "\t" << errorRateRleOld << "\t" << L << "\t" << n << endl;
    }


    // Create and initialize the disjoint sets data structure needed below.
    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2*readCount;
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
        const uint64_t a0 = disjointSets.find_set(A0.getValue());
        const uint64_t b0 = disjointSets.find_set(B0.getValue());
        const uint64_t a1 = disjointSets.find_set(A1.getValue());
        const uint64_t b1 = disjointSets.find_set(B1.getValue());

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
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<uint64_t, uint64_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[ReadId(p.second)]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    uint64_t accumulatedOrientedReadCount = 0;
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
            std::array<OrientedReadId, 2> nextEdgeOrientedReadIds = nextEdge.orientedReadIds;
            nextEdgeOrientedReadIds[0].flipStrand();
            nextEdgeOrientedReadIds[1].flipStrand();
            SHASTA_ASSERT(nextEdgeOrientedReadIds == edge.orientedReadIds);
            SHASTA_ASSERT(nextEdge.alignmentId == alignmentId);
        }
        

        // Q(n) = (1 + δ/2ε)^n * e-δL
        // ε = 1e-4, δ = 5e-4
        // Store this pair of edges in our edgeTable.
        const uint64_t range0 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[0], 0, markers);
        const uint64_t range1 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[1], 1, markers);
        const double L = double(range0 + range1)/2.;
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
    for(uint64_t i = 0; i < min(uint64_t(100), edgeTable.size()); i++) {
        const auto& p = edgeTable[i];
        const uint64_t edgeId = p.first;
        const double logQ = p.second;
        const AlignmentData& alignment = alignmentData[readGraph.edges[edgeId].alignmentId];
        const uint64_t markerCount = alignment.info.markerCount;
        const uint64_t range0 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[0], 0, markers);
        const uint64_t range1 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[1], 1, markers);
        const double L = double(range0 + range1)/2.;
        const uint64_t n = alignment.info.mismatchCountRle;
        const double errorRateRle = double(n)/(2.0*L);;
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        cout << edgeId << " ( " << readId0 << "," << readId1 << " )" << "\t" << logQ << "\t" << markerCount << "\t" << errorRateRle << "\t" << L << "\t" << n << endl;
    }


    // Create and initialize the disjoint sets data structure needed below.
    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2*readCount;
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
        const uint64_t a0 = disjointSets.find_set(A0.getValue());
        const uint64_t b0 = disjointSets.find_set(B0.getValue());
        const uint64_t a1 = disjointSets.find_set(A1.getValue());
        const uint64_t b1 = disjointSets.find_set(B1.getValue());

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
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<uint64_t, uint64_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[ReadId(p.second)]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    uint64_t accumulatedOrientedReadCount = 0;
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






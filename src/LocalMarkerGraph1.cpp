// Shasta.
#include "LocalMarkerGraph1.hpp"
#include "MarkerGraph.hpp"
using namespace shasta;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"

// Standard library.
#include "fstream.hpp"
#include <queue>



LocalMarkerGraph1::LocalMarkerGraph1(
    const MarkerGraph& markerGraph,
    MarkerGraphVertexId startVertexId,
    uint64_t maxDistance,
    uint64_t minVertexCoverage,
    uint64_t minEdgeCoverage) :
    markerGraph(markerGraph)
{
    LocalMarkerGraph1& graph = *this;

    // Do a BFS to generate the vertices.
    // Edges will be created later.
    const vertex_descriptor vStart = addVertex(startVertexId, 0);
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    while(!q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalMarkerGraph1Vertex& vertex0 = graph[v0];
        const MarkerGraphVertexId vertexId0 = vertex0.vertexId;
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over outgoing edges.
        for(uint64_t edgeId: markerGraph.edgesBySource[vertexId0]) {
            const auto& edge = markerGraph.edges[edgeId];

            // If coverage is too low, skip it.
            if(markerGraph.edgeCoverage(edgeId) < minEdgeCoverage) {
                continue;
            }

            // Get the target vertex.
            const MarkerGraph::VertexId vertexId1 = edge.target;
            SHASTA_ASSERT(edge.source == vertexId0);
            SHASTA_ASSERT(vertexId1 < markerGraph.vertexCount());

            // If vertex coverage is too low, skip it.
            if(markerGraph.vertexCoverage(vertexId1) < minVertexCoverage) {
                continue;
            }

            // Add this vertex, if we don't already have it.
            if(not vertexMap.contains(vertexId1)) {
                const vertex_descriptor v1 = graph.addVertex(vertexId1, distance1);

                // Also enqueue it, unless it is at maximum distance.
                if(distance1 < maxDistance) {
                    q.push(v1);
                }
            }
        }

        // Loop over incoming edges.
        for(uint64_t edgeId: markerGraph.edgesByTarget[vertexId0]) {
            const auto& edge = markerGraph.edges[edgeId];

            // If coverage is too low, skip it.
            if(markerGraph.edgeCoverage(edgeId) < minEdgeCoverage) {
                continue;
            }

            // Get the source vertex.
            const MarkerGraph::VertexId vertexId1 = edge.source;
            SHASTA_ASSERT(edge.target == vertexId0);
            SHASTA_ASSERT(vertexId1 < markerGraph.vertexCount());

            // If vertex coverage is too low, skip it.
            if(markerGraph.vertexCoverage(vertexId1) < minVertexCoverage) {
                continue;
            }

            // Add this vertex, if we don't already have it.
            if(not vertexMap.contains(vertexId1)) {
                const vertex_descriptor v1 = graph.addVertex(vertexId1, distance1);

                // Also enqueue it, unless it is at maximum distance.
                if(distance1 < maxDistance) {
                    q.push(v1);
                }
            }
        }
    }



    // Create edges.
    BGL_FORALL_VERTICES(v0, graph, LocalMarkerGraph1) {
        const LocalMarkerGraph1Vertex& vertex0 = graph[v0];
        const MarkerGraphVertexId vertexId0 = vertex0.vertexId;

        for(uint64_t edgeId: markerGraph.edgesBySource[vertexId0]) {

            // If coverage is too low, skip it.
            if(markerGraph.edgeCoverage(edgeId) < minEdgeCoverage) {
                continue;
            }
            const auto& edge = markerGraph.edges[edgeId];

            const MarkerGraph::VertexId vertexId1 = edge.target;
            SHASTA_ASSERT(edge.source == vertexId0);
            SHASTA_ASSERT(vertexId1 < markerGraph.vertexCount());

            // If vertexId1 is in the local marker graph, add this edge.
            auto it = vertexMap.find(vertexId1);
            if(it != vertexMap.end()) {
                const vertex_descriptor v1 = it->second;
                add_edge(v0, v1, LocalMarkerGraph1Edge(edgeId), graph);
            } else {
            }
        }
    }
}



LocalMarkerGraph1::vertex_descriptor LocalMarkerGraph1::addVertex(
    MarkerGraphVertexId vertexId,
    uint64_t distance)
{
    LocalMarkerGraph1& graph = *this;

    SHASTA_ASSERT(not vertexMap.contains(vertexId));
    const vertex_descriptor v = add_vertex(LocalMarkerGraph1Vertex(vertexId, distance), graph);
    vertexMap.insert(make_pair(vertexId, v));

    return v;
}



void LocalMarkerGraph1::writeGfa(const string& fileName) const
{
    const LocalMarkerGraph1& graph = *this;
    ofstream gfa(fileName);

    // Write the header.
    gfa << "H\tVN:Z:1.0\n";

    // Write one segment for each edge.
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph1) {
        const MarkerGraphEdgeId edgeId = graph[e].edgeId;
        gfa <<
            "S\t" << edgeId << "\t" << "*"
            // "\tLN:i:" << path.size() <<
            "\n";
    }



    // Write the links.
    // For each vertex, we write links between all pairs of incomint/outgoing edges.
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph1) {
        BGL_FORALL_INEDGES(v, e0, graph, LocalMarkerGraph1) {
            const MarkerGraphEdgeId edgeId0 = graph[e0].edgeId;
            BGL_FORALL_OUTEDGES(v, e1, graph, LocalMarkerGraph1) {
                const MarkerGraphEdgeId edgeId1 = graph[e1].edgeId;
                gfa << "L\t" <<
                    edgeId0 << "\t+\t" <<
                    edgeId1 << "\t+\t0M\n";
            }
        }
    }
}

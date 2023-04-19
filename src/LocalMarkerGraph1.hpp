#ifndef SHASTA_LOCAL_MARKER_GRAPH1_HPP
#define SHASTA_LOCAL_MARKER_GRAPH1_HPP

// Shasta.
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include "string.hpp"

namespace shasta {

    class LocalMarkerGraph1Vertex;
    class LocalMarkerGraph1Edge;
    class LocalMarkerGraph1;
    using LocalMarkerGraph1BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        LocalMarkerGraph1Vertex,
        LocalMarkerGraph1Edge
        >;

    class MarkerGraph;
}


class shasta::LocalMarkerGraph1Vertex {
public:

    // The id of the corresponding marker graph vertex.
    MarkerGraphVertexId vertexId;

    // The distance from the start vertex.
    uint64_t distance;

    LocalMarkerGraph1Vertex(
        MarkerGraphVertexId vertexId,
        uint64_t distance) :
        vertexId(vertexId),
        distance(distance)
    {
    }

};



class shasta::LocalMarkerGraph1Edge {
public:

    // The id of the corresponding marker graph edge.
    MarkerGraphEdgeId edgeId;

    LocalMarkerGraph1Edge(MarkerGraphEdgeId edgeId) :
        edgeId(edgeId)
    {
    }

};



class shasta::LocalMarkerGraph1 :
    public LocalMarkerGraph1BaseClass {
public:

    LocalMarkerGraph1(
        const MarkerGraph&,
        MarkerGraphVertexId,
        uint64_t maxDistance,
        uint64_t minVertexCoverage,
        uint64_t minEdgeCoverage
    );

    const MarkerGraph& markerGraph;

    std::map<MarkerGraphVertexId, vertex_descriptor> vertexMap;
    vertex_descriptor addVertex(MarkerGraphVertexId, uint64_t distance);

    void writeGfa(const string& fileName) const;
};

#endif

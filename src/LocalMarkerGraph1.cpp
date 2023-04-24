// Shasta.
#include "LocalMarkerGraph1.hpp"
#include "Base.hpp"
#include "computeLayout.hpp"
#include "MarkerGraph.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include "fstream.hpp"
#include <queue>

#include <chrono>
#include <thread>



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
            "S\t" << edgeId << "\t";

        auto sequence = markerGraph.edgeSequence[edgeId];
        copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(gfa));

        // RC is multiplied by sequence length so reports the number of reads
        // (edge coverage) as depth.
        gfa <<
            "\tLN:i:" << sequence.size() <<
            "\tRC:i:" << sequence.size() * markerGraph.edgeCoverage(edgeId) <<
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



void LocalMarkerGraph1::writeHtml0(
    ostream& html,
    uint64_t sizePixels,
    uint64_t quality,
    double timeout,
    bool useSvg) const
{
    const LocalMarkerGraph1& graph = *this;

    // Compute the layout.
    std::map<edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph1) {
        edgeLengthMap.insert(make_pair(e, 1.));
    }
    std::map<vertex_descriptor, array<double, 2> > positionMap;
    const auto t0 = steady_clock::now();
    const ComputeLayoutReturnCode returnCode = computeLayoutCustom(
        graph, edgeLengthMap, positionMap, quality, timeout);
    const auto t1 = steady_clock::now();
    html << "<br>Graph layout computation took " << seconds(t1 - t0) << "s.";
    if(returnCode == ComputeLayoutReturnCode::Timeout) {
        throw runtime_error("Graph layout took too long. "
            "Increase the timeout or decrease the maximum distance.");
    }
    if(returnCode != ComputeLayoutReturnCode::Success) {
        throw runtime_error("Graph layout failed.");
    }

    // Find minimum and maximum of x and y.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = xMin;
    double yMax = xMax;
    for(const auto& p: positionMap) {
        const auto& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        xMin = min(xMin, x);
        xMax = max(xMax, x);
        yMin = min(yMin, y);
        yMax = max(yMax, y);
    }
    const double range = max(xMax - xMin, yMax - yMin);
    const double factor = double(sizePixels) / range;



    // Gather positions, discretized to integers.
    // Each of these will generate a pixel.
    class PixelInfo {
    public:
        uint64_t maxCoverage;
        MarkerGraphVertexId vertexId;
    };
    std::map< pair<int64_t, int64_t>, PixelInfo> pixels;
    for(const auto& p: positionMap) {
        const vertex_descriptor v = p.first;
        const auto& xy = p.second;
        const MarkerGraphVertexId vertexId = graph[v].vertexId;
        const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
        const double x = xy[0];
        const double y = xy[1];
        const uint64_t ix = int64_t(x * factor);
        const uint64_t iy = int64_t(y * factor);
        auto it = pixels.find({ix, iy});
        if(it == pixels.end()) {
            pixels.insert(make_pair(make_pair(ix, iy), PixelInfo({coverage, vertexId})));
        } else {
            if(coverage > it->second.maxCoverage) {
                it->second.maxCoverage = coverage;
                it->second.vertexId = vertexId;
            }
        }
    }



    // Find minimum and maximum ix, iy.
    int64_t ixMin = std::numeric_limits<int64_t>::max();
    int64_t ixMax = std::numeric_limits<int64_t>::min();
    int64_t iyMin = ixMin;
    int64_t iyMax = ixMax;
    for(const auto& pixel :pixels) {
        const auto& ixy = pixel.first;
        ixMin = min(ixMin, ixy.first);
        ixMax = max(ixMax, ixy.first);
        iyMin = min(iyMin, ixy.second);
        iyMax = max(iyMax, ixy.second);
    }

    const int64_t width = ixMax - ixMin + 1;
    const int64_t height = iyMax - iyMin + 1;



    if(useSvg) {

        // Display using svg.
        html << "\n<br><svg width=" << width << " height=" << height << ">";
        const string coverage1Color = "red";
        const string coverage2Color = "yellow";
        const string highCoverageColor = "black";

        for(const auto& pixel :pixels) {
            const auto& ixy = pixel.first;
            const uint64_t coverage = pixel.second.maxCoverage;
            const MarkerGraphVertexId vertexId = pixel.second.vertexId;
            const int64_t ix = ixy.first - ixMin;
            SHASTA_ASSERT(ix >= 0);
            SHASTA_ASSERT(ix < width);
            const int64_t iy = ixy.second - iyMin;
            SHASTA_ASSERT(iy >= 0);
            SHASTA_ASSERT(iy < height);

            string color;
            if(coverage == 1) {
                color = coverage1Color;
            } else if(coverage == 2) {
                color = coverage2Color;
            } else {
                color = highCoverageColor;
            }

            html <<
                "\n<a href='"
                "exploreMarkerGraph1?vertexId=" << vertexId << "&outputType=createAndOpenGfa"
                "'>"
                "<line x1=" << ix << " y1=" << iy << " x2=" << ix << " y2=" << iy <<
                " stroke=" << color << " stroke-width=1px stroke-linecap=square />"
                "</a>";

        }


        html << "</svg>";



    } else {

        // Display using canvas
        const array<uint8_t, 3> coverage1Color = {255, 0, 0};
        const array<uint8_t, 3> coverage2Color = {255, 255, 0};
        const array<uint8_t, 3> highCoverageColor = {0, 0, 0};
        html <<
            "\n<br><canvas id=canvas width=" << width << " height=" << height <<
            ">"
            "\n <script>"
            "\n var canvas = document.getElementById('canvas');"
            "\n var ctx = canvas.getContext('2d');"
            "\n var i = ctx.createImageData(" << width << "," << height << ");\n";
        for(const auto& pixel :pixels) {
            const auto& ixy = pixel.first;
            const uint64_t coverage = pixel.second.maxCoverage;
            const int64_t ix = ixy.first - ixMin;
            SHASTA_ASSERT(ix >= 0);
            SHASTA_ASSERT(ix < width);
            const int64_t iy = ixy.second - iyMin;
            SHASTA_ASSERT(iy >= 0);
            SHASTA_ASSERT(iy < height);
            const uint64_t index = (4 * width) * iy + 4 * ix;
            if(coverage == 1) {
                for(uint64_t k=0; k<3; k++) {
                    html << "i.data[" << index+k << "]=" << int(coverage1Color[k]) << ";";
                }
            } else if(coverage == 2) {
                for(uint64_t k=0; k<3; k++) {
                    html << "i.data[" << index+k << "]=" << int(coverage2Color[k]) << ";";
                }
            } else {
                for(uint64_t k=0; k<3; k++) {
                    html << "i.data[" << index+k << "]=" << int(highCoverageColor[k]) << ";";
                }
            }
            html << "i.data[" << index+3 << "]=255;";
        }
        html <<
            "\n ctx.putImageData(i, 0, 0);"
            "\n </script>";
    }

}



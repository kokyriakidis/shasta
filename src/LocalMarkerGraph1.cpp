// Shasta.
#include "LocalMarkerGraph1.hpp"
#include "Base.hpp"
#include "computeLayout.hpp"
#include "findLinearChains.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MurmurHash2.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;

// Boost libraries.
#include "boost/graph/filtered_graph.hpp"
#include "boost/graph/iteration_macros.hpp"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include "fstream.hpp"
#include <queue>
#include <stack>



LocalMarkerGraph1::LocalMarkerGraph1(
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    MarkerGraphVertexId startVertexId,
    uint64_t maxDistance,
    uint64_t minVertexCoverage,
    uint64_t minEdgeCoverage) :
    markers(markers),
    markerGraph(markerGraph),
    maxDistance(maxDistance)
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
                edge_descriptor e;
                tie(e, ignore) = add_edge(v0, v1, LocalMarkerGraph1Edge(edgeId), graph);
                edgeMap.insert({edgeId, e});
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
    // const auto t0 = steady_clock::now();
    const ComputeLayoutReturnCode returnCode = computeLayoutCustom(
        graph, edgeLengthMap, positionMap, quality, timeout);
    // const auto t1 = steady_clock::now();
    // html << "<br>Graph layout computation took " << seconds(t1 - t0) << "s.";
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



void LocalMarkerGraph1::writeHtml1(
    ostream& html,
    uint64_t sizePixels,
    double thicknessScaling,
    uint64_t quality,
    double edgeResolution,
    const string& coloring,
    uint64_t redCoverage,
    uint64_t greenCoverage,
    MarkerGraphEdgeId readFollowingStartEdgeId,
    int64_t firstMarkerOffset,
    int64_t lastMarkerOffset,
    bool showLabels,
    double timeout) const
{
    const LocalMarkerGraph1& graph = *this;



    // To compute the layout, use an auxiliary graph with a vertex
    // for each vertex of the LocalMarkerGraph1 plus zero or more vertices
    // for each edge of the LocalMarkerGraph1.
    // In this initial implementation we divide each LocalMarkerGraph1 edge into a number
    // of AuxiliaryGraph edges equal to the number of bases in its sequence.
    using AuxiliaryGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
    AuxiliaryGraph auxiliaryGraph;

    // The auxiliary graph vertex corresponding to each vertex of the LocalMarkerGraph1.
    std::map<vertex_descriptor, AuxiliaryGraph::vertex_descriptor> auxiliaryVertexMap;

    // The auxiliary graph vertices corresponding to each edge of the LocalMarkerGraph1.
    std::map<edge_descriptor, vector<AuxiliaryGraph::vertex_descriptor> > auxiliaryEdgeMap;

    // The desired length of each edge of the auxiliary graph.
    std::map<AuxiliaryGraph::edge_descriptor, double> auxiliaryEdgeLengthMap;

    // Create vertices and edges of the AuxiliaryGraph.
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph1) {
        auxiliaryVertexMap.insert({v, add_vertex(auxiliaryGraph)});
    }
    vector<AuxiliaryGraph::vertex_descriptor> auxiliaryVertices;
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph1) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const MarkerGraphEdgeId edgeId = graph[e].edgeId;
        const uint64_t sequenceLength = markerGraph.edgeSequence[edgeId].size();
        const uint64_t auxiliaryVertexCount = max(1UL, uint64_t(edgeResolution * double(sequenceLength)));
        const double edgeLength = double(sequenceLength) / double(auxiliaryVertexCount + 1);
        auxiliaryVertices.clear();
        for(uint64_t i=0; i<auxiliaryVertexCount; i++) {
            auxiliaryVertices.push_back(add_vertex(auxiliaryGraph));
        }
        auxiliaryEdgeMap.insert({e, auxiliaryVertices});

        // Add the necessary auxiliary graph edges.
        AuxiliaryGraph::edge_descriptor ae;
        if(auxiliaryVertexCount == 0) {
            tie(ae, ignore) = add_edge(auxiliaryVertexMap[v0], auxiliaryVertexMap[v1], auxiliaryGraph);
            auxiliaryEdgeLengthMap.insert({ae, edgeLength});
        } else {
            tie(ae, ignore) = add_edge(auxiliaryVertexMap[v0], auxiliaryVertices.front(), auxiliaryGraph);
            auxiliaryEdgeLengthMap.insert({ae, edgeLength});
            for(uint64_t i=1; i<auxiliaryVertexCount; i++) {
                tie(ae, ignore) = add_edge(auxiliaryVertices[i-1], auxiliaryVertices[i], auxiliaryGraph);
                auxiliaryEdgeLengthMap.insert({ae, edgeLength});
            }
            tie(ae, ignore) = add_edge(auxiliaryVertices.back(), auxiliaryVertexMap[v1], auxiliaryGraph);
            auxiliaryEdgeLengthMap.insert({ae, edgeLength});
        }
    }

    // Compute the layout of the auxiliary graph.
    std::map<AuxiliaryGraph::vertex_descriptor, array<double, 2> > positionMap;
    computeLayoutCustom(auxiliaryGraph, auxiliaryEdgeLengthMap, positionMap, quality, timeout);



    // If we are doing read following, we need to compute
    // followed read coverage for each edge.
    std::map<edge_descriptor, uint64_t> readFollowingCoverageMap;
    uint64_t readFollowingStartEdgeCoverage = 0;
    if(coloring == "readFollowing") {
        readFollowingStartEdgeCoverage = markerGraph.edgeCoverage(readFollowingStartEdgeId);

        // Loop over the MarkerIntervals of the start edge for read following.
        for(const MarkerInterval& startMarkerInterval:
            markerGraph.edgeMarkerIntervals[readFollowingStartEdgeId]) {
            const OrientedReadId orientedReadId = startMarkerInterval.orientedReadId;
            const int64_t startOrdinal0 = int64_t(startMarkerInterval.ordinals[0]);

            // The number of markers in this oriented read.
            const int64_t orientedReadMarkerCount = int64_t(markers.size(orientedReadId.getValue()));

            // Get the MarkerId of the first marker of this oriented read.
            // We can use this later to easily get the MarkerId corresponding to any
            // marker in the same oriented read.
            const MarkerId firstOrientedReadMarkerId =
                markers.begin(orientedReadId.getValue()) - markers.begin();

            // Loop over the requested range of offsets.
            for(int64_t offset=firstMarkerOffset; offset<=lastMarkerOffset; offset++) {
                const int64_t ordinal0 = startOrdinal0 + offset;
                if(ordinal0 < 0) {
                    // This offset takes us before the beginning of this oriented read.
                    continue;
                }
                const int64_t ordinal1 = ordinal0 + 1;
                if(ordinal1 > orientedReadMarkerCount-1) {
                    // This offset takes us past the end of this oriented read.
                    continue;
                }

                // Find the MarkerIds corresponding to these two ordinals.
                const MarkerId markerId0 = firstOrientedReadMarkerId + ordinal0;
                const MarkerId markerId1 = firstOrientedReadMarkerId + ordinal1;

                // Find the corresponding marker graph vertices.
                // We are using the complete marker graph, so the vertices must exist.
                const MarkerGraph::CompressedVertexId compressedVertexId0 = markerGraph.vertexTable[markerId0];
                const MarkerGraph::CompressedVertexId compressedVertexId1 = markerGraph.vertexTable[markerId1];
                SHASTA_ASSERT(compressedVertexId0 != MarkerGraph::invalidCompressedVertexId);
                SHASTA_ASSERT(compressedVertexId1 != MarkerGraph::invalidCompressedVertexId);
                const MarkerGraphVertexId vertexId0 = compressedVertexId0;
                // const MarkerGraphVertexId vertexId1 = compressedVertexId1;

                // Find the edge vertexId0->vertexId1 that contains the MarkerInterval
                // with these oriented read and ordinals.
                MarkerInterval targetMarkerInterval(orientedReadId, uint32_t(ordinal0), uint32_t(ordinal1));
                MarkerGraphEdgeId edgeId = invalid<MarkerGraphEdgeId>;
                for(const MarkerGraphEdgeId candidateEdgeId: markerGraph.edgesBySource[vertexId0]) {
                    const auto edgeMarkerIntervals = markerGraph.edgeMarkerIntervals[candidateEdgeId];
                    if(find(edgeMarkerIntervals.begin(), edgeMarkerIntervals.end(), targetMarkerInterval)
                        != edgeMarkerIntervals.end()) {
                        edgeId = candidateEdgeId;
                        break;
                    }
                }
                SHASTA_ASSERT(edgeId != invalid<MarkerGraphEdgeId>);

                // cout << orientedReadId << " at offset " << offset << endl;

                // If this edge is in the LocalMarkerGraph1, increment its read following coverage.
                auto it = edgeMap.find(edgeId);
                if(it != edgeMap.end()) {
                    const edge_descriptor e = it->second;
                    auto jt = readFollowingCoverageMap.find(e);
                    if(jt == readFollowingCoverageMap.end()){
                        readFollowingCoverageMap.insert({e, 1});
                        // cout << "Added a new entry in the readFollowingCoverageMap." << endl;
                    } else {
                        ++jt->second;
                        // cout << "Incremented readFollowingCoverageMap to " << jt->second << endl;
                    }
                } else {
                    // cout << "Not found in the edge map." << endl;
                }
            }
        }

        /*
        for(const auto& p: readFollowingCoverageMap) {
            const edge_descriptor e = p.first;
            const uint64_t coverage = p.second;
            cout << graph[e].edgeId << " " << coverage << endl;
        }
        */
    }



    // Compute the view box.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = xMin;
    double yMax = xMax;
    for(const auto& p: positionMap) {
        const array<double, 2>& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        xMin = min(xMin, x);
        xMax = max(xMax, x);
        yMin = min(yMin, y);
        yMax = max(yMax, y);
    }
    const double extend = thicknessScaling;
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;
    const double fontSize = 16. * max(xMax-xMin, yMax-yMin) / double(sizePixels);

    // Make the "arrow" length equal to the desired length of 1 base.
    const double arrowLength = 1.;

    // Begin the svg.
    const string svgId = "LocalMarkerGraph1";
    html << "\n<div style='display: inline-block; vertical-align:top'>"
        "<br><svg id='" << svgId <<
        "' width='" <<  sizePixels <<
        "' height='" << sizePixels <<
        "' viewbox='" << xMin << " " << yMin << " " <<
        xMax - xMin << " " <<
        yMax - yMin << "'"
        " font-size='" << fontSize << "' style='border-style:solid;border-color:Black;stroke-linecap:round'"
        " font-family=monospace"
        ">\n";



    // Write the edges.
    html << "\n<g id=edges stroke-width='" << thicknessScaling << "'>";
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph1) {
        const MarkerGraphEdgeId edgeId = graph[e].edgeId;
        const uint64_t coverage = markerGraph.edgeCoverage(edgeId);
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const auto& p0 = positionMap[auxiliaryVertexMap[v0]];
        const auto& p1 = positionMap[auxiliaryVertexMap[v1]];
        const vector<AuxiliaryGraph::vertex_descriptor>& auxiliaryVertices =  auxiliaryEdgeMap[e];

        string color;
        uint64_t readFollowingCoverage = 0;
        if(coloring == "random") {
            const uint32_t hue = MurmurHash2(&edgeId, sizeof(edgeId), 231) % 360;
            color = "hsl(" + to_string(hue) + ",50%,50%)";
        } else if(coloring == "byCoverage") {
            if(coverage <= redCoverage) {
                color = "Red";
            } else if(coverage >= greenCoverage) {
                color = "Green";
            } else {
                const uint32_t hue = uint32_t(120. *
                    (double(coverage) - double(redCoverage)) / (double(greenCoverage) - double(redCoverage)));
                color = "hsl(" + to_string(hue) + ",50%,50%)";
            }
        } else if(coloring == "readFollowing") {
            auto it = readFollowingCoverageMap.find(e);
            if(it == readFollowingCoverageMap.end()) {
                color = "rgba(192,192,192,0.4)";    // Transparent grey
            } else {
                const uint64_t coverage = it->second;
                readFollowingCoverage = coverage;
                if(coverage <= redCoverage) {
                    color = "Red";
                } else if(coverage >= greenCoverage) {
                    color = "Green";
                } else {
                    const uint32_t hue = uint32_t(120. *
                        (double(coverage) - double(redCoverage)) / (double(greenCoverage) - double(redCoverage)));
                    color = "hsl(" + to_string(hue) + ",50%,50%)";
                }
            }
        } else {
            SHASTA_ASSERT(0);
        }
        const string properties = "stroke='" + color + "'";

        SHASTA_ASSERT(not auxiliaryVertices.empty());

        // Create a group for this edge.
        const auto sequence = markerGraph.edgeSequence[edgeId];
        html << "<g>";

        // Add a title.
        html <<
            "<title>Edge " << edgeId << ", coverage " << coverage <<
            ", " << sequence.size() << " bases: ";
        copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(html));
        if(coloring == "readFollowing") {
            html << ", read following coverage " << readFollowingCoverage << "/" <<
                readFollowingStartEdgeCoverage;
        }
        html << "</title>";

        // Add a hyperlink.
        html << "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId << "'>";

        // Line from p0 to the first auxiliary vertex.
        const auto& xyFirst = positionMap[auxiliaryVertices.front()];
        html << "\n<line x1=" << p0[0] << " y1=" << p0[1] <<
            " x2=" << xyFirst[0] << " y2=" << xyFirst[1] << " " << properties << " />";

        // Lines between auxiliary vertices.
        for(uint64_t i=1; i<auxiliaryVertices.size(); i++) {
            const auto& xyA = positionMap[auxiliaryVertices[i-1]];
            const auto& xyB = positionMap[auxiliaryVertices[i]];
            html << "\n<line x1=" << xyA[0] << " y1=" << xyA[1] <<
                " x2=" << xyB[0] << " y2=" << xyB[1] << " " << properties << " />";
        }

        // Line from the last auxiliary vertex to p1.
        const auto& xyLast = positionMap[auxiliaryVertices.back()];
        html << "\n<line x1=" << xyLast[0] << " y1=" << xyLast[1] <<
            " x2=" << p1[0] << " y2=" << p1[1] << " " << properties << " />";
        html << "</a>";

        // Label.
        if(showLabels) {
            double x, y;
            if((auxiliaryVertices.size() %2) == 0) {
                const auto positionA = positionMap[auxiliaryVertices[auxiliaryVertices.size()/2 -1]];
                const auto positionB = positionMap[auxiliaryVertices[auxiliaryVertices.size()/2]];
                x = (positionA[0] + positionB[0]) / 2;
                y = (positionA[1] + positionB[1]) / 2;
            } else {
                const auto position = positionMap[auxiliaryVertices[auxiliaryVertices.size()/2]];
                x = position[0];
                y = position[1];
            }
            html << "<text x='" << x << "' << y='" << y << "' dominant-baseline=middle text-anchor=middle>";
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(html));
            html << "</text>";
        }

        // End the group for this edge.
        html << "</g>";
    }
    html << "\n</g>";



    // Write the "arrows".
    html << "\n<g id=arrows stroke-width='" << thicknessScaling/3. << "'>";
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph1) {
        const vertex_descriptor v1 = target(e, graph);
        const auto& p1 = positionMap[auxiliaryVertexMap[v1]];
        const vector<AuxiliaryGraph::vertex_descriptor>& auxiliaryVertices =  auxiliaryEdgeMap[e];
        SHASTA_ASSERT(not auxiliaryVertices.empty());

        // Position of the last auxiliary vertex.
        const auto& xyLast = positionMap[auxiliaryVertices.back()];

        // Draw the "arrow".
        // We need to compute a unit vector in the direction (p1, xyLast).
        const double vx = xyLast[0] - p1[0];
        const double vy = xyLast[1] - p1[1];
        const double v = sqrt(vx*vx + vy * vy);
        if(v < 1.e-3) {
            // Trouble. This can happen if two vertices are very close. Skip the arrow.
            continue;
        }
        const double ux = vx / v;
        const double uy = vy / v;
        const double xArrow = p1[0] + ux * arrowLength;
        const double yArrow = p1[1] + uy * arrowLength;
        html << "\n<line x1=" << xArrow << " y1=" << yArrow <<
            " x2=" << p1[0] << " y2=" << p1[1] << " stroke=Black />";
    }
    html << "\n</g>";


    // Write the vertices.
    html << "\n<g id=vertices stroke-width='" << thicknessScaling << "'>";
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph1) {
        const auto& p = positionMap[auxiliaryVertexMap[v]];
        const double x = p[0];
        const double y = p[1];
        const string color = (graph[v].distance == maxDistance ? "Grey" : "Black");

        // Create a group for this edge.
        const MarkerGraphVertexId vertexId = graph[v].vertexId;
        const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
        html << "<g><title>Vertex " << vertexId << ", coverage " << coverage;
        html << "</title>";
        html << "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId << "'>";

        // Write the vertex.
        html << "\n<line x1=" << x << " y1=" << y <<
            " x2=" << x << " y2=" << y << " stroke=" << color << " />";

        // End the group.
        html << "</a></g>";
    }
    html << "\n</g>";

    // Finish the svg.
    html << "\n</svg></div>";

    // Add drag and zoom.
    addSvgDragAndZoom(html);

    // Side panel.
    html << "<div style='display: inline-block'>";

    // Change thickness
    html << R"stringDelimiter(
    <p><table>
    <tr><th class=left>Thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>
        <script>
        function changeThickness(factor)
        {
            edges = document.getElementById('edges');
            edges.setAttribute('stroke-width', factor * edges.getAttribute('stroke-width'));
            vertices = document.getElementById('vertices');
            vertices.setAttribute('stroke-width', factor * vertices.getAttribute('stroke-width'));
            arrows = document.getElementById('arrows');
            arrows.setAttribute('stroke-width', factor * arrows.getAttribute('stroke-width'));
        }
        </script>
        )stringDelimiter";



    // Zoom buttons.
    html << R"stringDelimiter(
        <tr title='Or use the mouse wheel.'><th class=left>Zoom<td>
        <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
        <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
        <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
        <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
        <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
        <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>
    )stringDelimiter";



    // End the side panel.
    html << "</table></div>";
}



void LocalMarkerGraph1::pruneLowCoverageLeaves(uint64_t maxPruneEdgeCoverage)
{
    if(maxPruneEdgeCoverage == 0) {
        return;
    }

    pruneLowCoverageForwardLeaves(maxPruneEdgeCoverage);
    pruneLowCoverageBackwardLeaves(maxPruneEdgeCoverage);

}



void LocalMarkerGraph1::pruneLowCoverageForwardLeaves(uint64_t maxPruneCoverage)
{
    LocalMarkerGraph1& graph = *this;

    // Start will all vertices with out-degree 0 and low coverage.
    std::stack<vertex_descriptor> leaves;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph1) {
        const MarkerGraphVertexId vertexId = graph[v].vertexId;
        const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
        if(coverage > maxPruneCoverage) {
            continue;
        }
        if(out_degree(v, graph) == 0) {
            leaves.push(v);
        }
    }

    // Main loop. At each iteration we remove a leaf, and add others as required.
    while(not leaves.empty()) {
        const vertex_descriptor leaf = leaves.top();
        leaves.pop();

        // If any parent has out-degree 1 and low coverage,
        // it becomes a leaf to be removed when we remove this one.
        BGL_FORALL_INEDGES(leaf, e, graph, LocalMarkerGraph1) {
            const vertex_descriptor parent = source(e, graph);
            if(parent == leaf) {
                continue;
            }
            const MarkerGraphVertexId vertexId = graph[parent].vertexId;
            const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
            if(coverage > maxPruneCoverage) {
                continue;
            }
            if(out_degree(parent, graph) == 1) {
                leaves.push(parent);
            }
        }

        clear_vertex(leaf, graph);
        remove_vertex(leaf, graph);
    }
}



void LocalMarkerGraph1::pruneLowCoverageBackwardLeaves(uint64_t maxPruneCoverage)
{
    LocalMarkerGraph1& graph = *this;

    // Start will all vertices with in-degree 0 and low coverage.
    std::stack<vertex_descriptor> leaves;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph1) {
        const MarkerGraphVertexId vertexId = graph[v].vertexId;
        const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
        if(coverage > maxPruneCoverage) {
            continue;
        }
        if(in_degree(v, graph)==0) {
            leaves.push(v);
        }
    }

    // Main loop. At each iteration we remove a leaf, and add others as required.
    while(not leaves.empty()) {
        const vertex_descriptor leaf = leaves.top();
        leaves.pop();

        // If any child has in-degree 1 and low coverage,
        // it becomes a leaf to be removed when we remove this one.
        BGL_FORALL_OUTEDGES(leaf, e, graph, LocalMarkerGraph1) {
            const vertex_descriptor child = target(e, graph);
            if(child == leaf) {
                continue;
            }
            const MarkerGraphVertexId vertexId = graph[child].vertexId;
            const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
            if(coverage > maxPruneCoverage) {
                continue;
            }
            if(in_degree(child, graph)==1) {
                leaves.push(child);
            }
        }

        clear_vertex(leaf, graph);
        remove_vertex(leaf, graph);
    }

}



void LocalMarkerGraph1::findLowCoverageChains(
    uint64_t maxChainCoverage,
    vector< vector<vertex_descriptor> >& chains
    ) const
{
    const LocalMarkerGraph1& graph = *this;

    // Create a filtered graph containing only the vertices
    // with coverage up to maxChainCoverage.
    class VertexPredicate {
    public:
        VertexPredicate() : graph(0), maxChainCoverage(invalid<uint64_t>) {}
        VertexPredicate(
            const LocalMarkerGraph1& graph,
            uint64_t maxChainCoverage) :
            graph(&graph),
            maxChainCoverage(maxChainCoverage)
            {}
        const LocalMarkerGraph1* graph;
        uint64_t maxChainCoverage;
        bool operator()(const vertex_descriptor v) const
        {
            const MarkerGraphVertexId vertexId = (*graph)[v].vertexId;
            const uint64_t coverage = graph->markerGraph.vertexCoverage(vertexId);
            return coverage <= maxChainCoverage;
        }
    };
    boost::filtered_graph<LocalMarkerGraph1, boost::keep_all, VertexPredicate>
        filteredGraph(graph,  boost::keep_all(), VertexPredicate(graph, maxChainCoverage));

    // Find linear chains in this filtered graph.
    findLinearVertexChains(filteredGraph, chains);
}



void LocalMarkerGraph1::removeLongLowCoverageChains(
    uint64_t maxChainCoverage,
    uint64_t minLength)
{
    LocalMarkerGraph1& graph = *this;

    // Find low coverage chains.
    vector< vector<LocalMarkerGraph1::vertex_descriptor> > lowCoverageChains;
    findLowCoverageChains(1, lowCoverageChains);

    // Remove the long ones.
    for(const auto& chain: lowCoverageChains) {
        if(chain.size() >= minLength) {
            for(const vertex_descriptor v: chain) {
                clear_vertex(v, graph);
                remove_vertex(v, graph);
            }
        }
    }

}


// Shasta.
#include "Assembler.hpp"
#include "computeLayout.hpp"
#include "html.hpp"
#include "mode3b-PathFinder.hpp"
#include "orderPairs.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>
#include <queue>



void Assembler::exploreMode3bPathGraph(const vector<string>& request, ostream& html)
{
#if 0
    // Get these from a form.
    const uint64_t maxMarkerOffset = 1000;
    const uint64_t minCoverage = 15;
    const uint64_t maxCoverage = 35;
    const uint64_t minCommonCount = 6;
    const double minCorrectedJaccard = 0.8;
    const uint64_t maxEdgeCount = 4;
    const uint64_t maxDistance = 6;
#endif


    // Get the request parameters.
    uint64_t edgeId = invalid<uint64_t>;
    getParameterValue(request, "edgeId", edgeId);

    uint64_t maxMarkerOffset = 10000;
    getParameterValue(request, "maxMarkerOffset", maxMarkerOffset);

    uint64_t minCoverage = 10;
    getParameterValue(request, "minCoverage", minCoverage);

    uint64_t maxCoverage = 30;
    getParameterValue(request, "maxCoverage", maxCoverage);

    uint64_t minCommonCount = 6;
    getParameterValue(request, "minCommonCount", minCommonCount);

    double minCorrectedJaccard = 0.8;
    getParameterValue(request, "minCorrectedJaccard", minCorrectedJaccard);

    uint64_t maxEdgeCount = 4;
    getParameterValue(request, "maxEdgeCount", maxEdgeCount);

    uint64_t maxDistance = 4;
    getParameterValue(request, "maxDistance", maxDistance);

    // The path direction can be forward, backward, or bidirectional.
    string graphDirection = "bidirectional";
    HttpServer::getParameterValue(request, "graphDirection", graphDirection);

    // Other quantities to add to the request.
    const double sizePixels = 900.;
    const double initialVertexRadiusPixels = 2.;
    const double initialEdgeThicknessPixels = 0.5;
    const double timeout = 600.;


    // Write the form.
    html <<
        "<h2>Path graph</h2>"
        "<form>"

        "<table>"

        "<tr>"
        "<td>Marker graph edge to start from"
        "<td class=centered><input type=text required name=edgeId size=8 style='text-align:center'"
        << ((edgeId == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeId) + "'")) <<
        ">"

        "<tr>"
        "<td>maxMarkerOffset"
        "<td class=centered>"
        "<input type=text required name=maxMarkerOffset size=8 style='text-align:center'" <<
        "value='" << maxMarkerOffset << "'>"

        "<tr>"
        "<td>minCoverage"
        "<td class=centered>"
        "<input type=text required name=minCoverage size=8 style='text-align:center'" <<
        "value='" << minCoverage << "'>"

        "<tr>"
        "<td>maxCoverage"
        "<td class=centered>"
        "<input type=text required name=maxCoverage size=8 style='text-align:center'" <<
        "value='" << maxCoverage << "'>"

        "<tr>"
        "<td>minCommonCount"
        "<td class=centered>"
        "<input type=text required name=minCommonCount size=8 style='text-align:center'" <<
        "value='" << minCommonCount << "'>"

        "<tr>"
        "<td>minCorrectedJaccard"
        "<td class=centered>"
        "<input type=text required name=minCorrectedJaccard size=8 style='text-align:center'" <<
        "value='" << minCorrectedJaccard << "'>"

        "<tr>"
        "<td>maxEdgeCount"
        "<td class=centered>"
        "<input type=text required name=maxEdgeCount size=8 style='text-align:center'" <<
        "value='" << maxEdgeCount << "'>"

        "<tr>"
        "<td>maxDistance"
        "<td class=centered>"
        "<input type=text required name=maxDistance size=8 style='text-align:center'" <<
        "value='" << maxDistance << "'>"

        "<tr>"
        "<td>Direction"
        "<td>"
        "<input type=radio name=graphDirection value=forward" <<
        (graphDirection=="forward" ? " checked=checked" : "") << "> Forward"
        "<br><input type=radio name=graphDirection value=backward" <<
        (graphDirection=="backward" ? " checked=checked" : "") << "> Backward"
        "<br><input type=radio name=graphDirection value=bidirectional" <<
        (graphDirection=="bidirectional" ? " checked=checked" : "") << "> Both" <<

        "</table>"
        "<br><input type=submit value='Do it'>"
        "</form>";



    // If the start edge is missing, don't do anything.
    if(edgeId == invalid<uint64_t>) {
        return;
    }



    // The local path graph.
    class PathGraphVertex {
    public:
        MarkerGraphEdgeId edgeId;
        uint64_t distance;
    };
    class PathGraphEdge {
    public:
        MarkerGraphEdgePairInfo info;
    };
    using PathGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        PathGraphVertex,
        PathGraphEdge
        >;
    class PathGraph : public PathGraphBaseClass {
    public:
        std::map<MarkerGraphEdgeId, vertex_descriptor> vertexMap;
        vertex_descriptor addVertex(MarkerGraphEdgeId edgeId, uint64_t distance)
        {
            SHASTA_ASSERT(not vertexMap.contains(edgeId));
            const vertex_descriptor v = add_vertex(PathGraphVertex({edgeId, distance}), *this);
            vertexMap.insert(make_pair(edgeId, v));
            return v;
        }
        vertex_descriptor getVertex(MarkerGraphEdgeId edgeId)
        {
            auto it = vertexMap.find(edgeId);
            if(it == vertexMap.end()) {
                return null_vertex();
            } else {
                return it->second;
            }
        }
    };
    using vertex_descriptor = PathGraph::vertex_descriptor;
    using edge_descriptor = PathGraph::edge_descriptor;
    PathGraph pathGraph;

    // Create the PathFinder, which knows how to find edges with similar read composition.
    PathFinder pathFinder(*this);
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > nextPrimaryEdges;


    // BFS-like loop to create the local path graph.
    std::queue<vertex_descriptor> q;
    const vertex_descriptor v = pathGraph.addVertex(edgeId, 0);
    q.push(v);
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const PathGraphVertex& vertex0 = pathGraph[v0];
        const MarkerGraphEdgeId edgeId0 = vertex0.edgeId;
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            if(direction == 0 and graphDirection == "backward") {
                continue;
            }
            if(direction == 1 and graphDirection == "forward") {
                continue;
            }

            pathFinder.findNextPrimaryEdges(
                edgeId0,
                direction,
                minCoverage,
                maxCoverage,
                maxEdgeCount,
                maxMarkerOffset,
                minCommonCount,
                minCorrectedJaccard,
                nextPrimaryEdges
                );

            // Loop over the next edges we found.
            for(auto& p: nextPrimaryEdges) {
                const MarkerGraphEdgeId edgeId1 = p.first;
                MarkerGraphEdgePairInfo& info = p.second;
                if(direction == 1) {
                    info.reverse();
                }

                vertex_descriptor v1 = pathGraph.getVertex(edgeId1);
                if(v1 == PathGraph::null_vertex()) {
                    v1 = pathGraph.addVertex(edgeId1, distance0 + 1);
                    if(distance1 < maxDistance) {
                        q.push(v1);
                    }
                }



                if(direction == 0) {
                    edge_descriptor e01;
                    bool edgeExists;
                    tie(e01, edgeExists) = edge(v0, v1, pathGraph);
                    if(not edgeExists) {
                        add_edge(v0, v1, PathGraphEdge({info}), pathGraph);
                    }
                } else {
                    edge_descriptor e10;
                    bool edgeExists;
                    tie(e10, edgeExists) = edge(v1, v0, pathGraph);
                    if(not edgeExists) {
                        add_edge(v1, v0, PathGraphEdge({info}), pathGraph);
                    }
                }

            }
        }
    }
    html << "<p>The local path graph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges.";



    // Compute the layout.
    std::map<edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        edgeLengthMap.insert(make_pair(e, pathGraph[e].info.offsetInBases));
    }
    std::map<vertex_descriptor, array<double, 2> > positionMap;
    const auto returnCode = computeLayoutCustom(pathGraph, edgeLengthMap, positionMap, 2, timeout);
    if(returnCode == ComputeLayoutReturnCode::Timeout) {
        throw runtime_error("Graph layout took too long. "
            "Increase the timeout or decrease the maximum distance.");
    }
    if(returnCode != ComputeLayoutReturnCode::Success) {
        throw runtime_error("Graph layout failed.");
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
    const double fontSize = 16. * max(xMax-xMin, yMax-yMin) / double(sizePixels);

    // Begin the svg.
    const string svgId = "PathGraph";
    html << "\n<div style='display: inline-block; vertical-align:top'>"
        "<br><svg id='" << svgId <<
        "' width='" <<  sizePixels <<
        "' height='" << sizePixels <<
        "' viewbox='" << xMin << " " << yMin << " " <<
        xMax - xMin << " " <<
        yMax - yMin << "'"
        " font-size='" << fontSize << "' style='border-style:solid;border-color:Black;stroke-linecap:round'"
        ">\n";



    // Write edges first, to avoid obscuring the vertices.
    const double edgeThickness = initialEdgeThicknessPixels * max(xMax-xMin, yMax-yMin) / sizePixels;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const PathGraphEdge& edge = pathGraph[e];

        const vertex_descriptor v1 = source(e, pathGraph);
        const vertex_descriptor v2 = target(e, pathGraph);

        auto it1 = positionMap.find(v1);
        SHASTA_ASSERT(it1 != positionMap.end());
        const auto& p1 = it1->second;
        const double x1 = p1[0];
        const double y1 = p1[1];

        auto it2 = positionMap.find(v2);
        SHASTA_ASSERT(it2 != positionMap.end());
        const auto& p2 = it2->second;
        const double x2 = p2[0];
        const double y2 = p2[1];

        // To show the direction, the two halves of the edge are
        // displayed in different colors.
        const double xm = 0.5 * (x1 + x2);
        const double ym = 0.5 * (y1 + y2);

        html <<
            "\n<g><title>" << pathGraph[v1].edgeId << "->" << pathGraph[v2].edgeId << ", offset " <<
            edge.info.offsetInBases << " bases</title>"
            "<line x1=" << x1 << " y1=" << y1 <<
            " x2=" << xm << " y2=" << ym <<
            " stroke=Black"
            " stroke-width=" << edgeThickness <<
            " />"
            "<line x1=" << xm << " y1=" << ym <<
            " x2=" << x2 << " y2=" << y2 <<
            " stroke=Green"
            " stroke-width=" << edgeThickness <<
            " /></g>";
    }


    // Write the vertices in order of decreasing distance
    // to avoid obscuring the ones at low distance.
    vector< pair<vertex_descriptor, uint64_t> > sortedVertices;
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        sortedVertices.push_back({v, pathGraph[v].distance});
    }
    sort(sortedVertices.begin(), sortedVertices.end(),
        OrderPairsBySecondOnlyGreater<vertex_descriptor, uint64_t>());


    // Write the vertices (each corresponding to a marker graph edge).
    const double vertexRadius = initialVertexRadiusPixels * max(xMax-xMin, yMax-yMin) / sizePixels;
    for(const auto& q: sortedVertices) {
        const vertex_descriptor v = q.first;
        const PathGraphVertex& vertex = pathGraph[v];
        auto it = positionMap.find(v);
        SHASTA_ASSERT(it != positionMap.end());
        const auto& p = it->second;
        const double x = p[0];
        const double y = p[1];

        string color = "Black";
        if(vertex.distance == 0) {
            color = "Purple";
        } else if(vertex.distance == maxDistance) {
            color = "Cyan";
        }

        html <<
            "<g>"
            "<title>" << vertex.edgeId << " distance " << vertex.distance << "</title>"
            "\n<circle cx=" << x << " cy=" << y << " r=" << vertexRadius <<
            " fill=" << color << " />"
            "</g>";

    }


    // Finish the svg.
    html << "\n</svg></div>";

    // Add drag and zoom.
    addSvgDragAndZoom(html);
}


// Shasta.
#include "mode3-LocalAnchorGraph.hpp"
#include "HttpServer.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>



LocalAnchorGraph::LocalAnchorGraph(
    const Anchors& anchors,
    const vector<AnchorId>& anchorIds,
    uint64_t maxDistance,
    bool filterEdgesByCoverageLoss,
    double maxCoverageLoss) :
    anchors(anchors),
    maxDistance(maxDistance)
{
    LocalAnchorGraph& graph = *this;

    // Initialize a BFS from these AnchorIds.
    std::queue<vertex_descriptor> q;
    for(const AnchorId anchorId: anchorIds) {
        SHASTA_ASSERT(not vertexMap.contains(anchorId));
        const vertex_descriptor v = boost::add_vertex(LocalAnchorGraphVertex(anchorId, 0), graph);
        vertexMap.insert({anchorId, v});
        q.push(v);
    }

    // BFS to find the vertices. We will add the edges later.
    vector<AnchorId> neighbors;
    vector<uint64_t> coverage;
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const AnchorId anchorId0 = vertex0.anchorId;
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        anchors.findChildren(anchorId0, neighbors, coverage);
        for(uint64_t i=0; i<neighbors.size(); i++) {
            const AnchorId anchorId1 = neighbors[i];
            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                continue;
            }

            // Filter by coverage loss, if requested.
            if(filterEdgesByCoverageLoss) {
                LocalAnchorGraphEdge edge;
                edge.coverage = coverage[i];
                anchors.analyzeAnchorPair(anchorId0, anchorId1, edge.info);
                if(edge.coverageLoss() > maxCoverageLoss) {
                    continue;
                }
            }

            const vertex_descriptor v1 = boost::add_vertex(LocalAnchorGraphVertex(anchorId1, distance1), graph);
            vertexMap.insert({anchorId1, v1});
            if(distance1 < maxDistance) {
                q.push(v1);
            }
        }

        anchors.findParents(anchorId0, neighbors, coverage);
        for(uint64_t i=0; i<neighbors.size(); i++) {
            const AnchorId anchorId1 = neighbors[i];
            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                continue;
            }

            // Filter by coverage loss, if requested.
            if(filterEdgesByCoverageLoss) {
                LocalAnchorGraphEdge edge;
                edge.coverage = coverage[i];
                anchors.analyzeAnchorPair(anchorId0, anchorId1, edge.info);
                if(edge.coverageLoss() > maxCoverageLoss) {
                    continue;
                }
            }

            const vertex_descriptor v1 = boost::add_vertex(LocalAnchorGraphVertex(anchorId1, distance1), graph);
            vertexMap.insert({anchorId1, v1});
            if(distance1 < maxDistance) {
                q.push(v1);
            }
        }
    }



    // Now add the edges.
    BGL_FORALL_VERTICES(v0, graph, LocalAnchorGraph) {
        const AnchorId anchorId0 = graph[v0].anchorId;
        anchors.findChildren(anchorId0, neighbors, coverage);
        for(uint64_t i=0; i<neighbors.size(); i++) {
            const AnchorId& anchorId1 = neighbors[i];
            auto it1 = vertexMap.find(anchorId1);
            if(it1 == vertexMap.end()) {
                continue;
            }
            const vertex_descriptor v1 = it1->second;

            // Create the tentative edge.
            LocalAnchorGraphEdge edge;
            edge.coverage = coverage[i];
            anchors.analyzeAnchorPair(anchorId0, anchorId1, edge.info);

            // Add it if requested.
            if((not filterEdgesByCoverageLoss) or
                (edge.coverageLoss() <= maxCoverageLoss)) {
                edge_descriptor e;
                tie(e, ignore) = add_edge(v0, v1, edge, graph);
            }
        }
    }
}



void LocalAnchorGraph::writeGraphviz(
    const string& fileName,
    double vertexSize,
    bool vertexSizeByCoverage,
    const string& edgeColoring,
    double edgeThickness,
    bool edgeThicknessByCoverage,
    double minimumEdgeLength,
    double additionalEdgeLengthPerKb,
    double arrowSize) const
{
    ofstream file(fileName);
    writeGraphviz(
        file,
        vertexSize, vertexSizeByCoverage,
        edgeColoring,
        edgeThickness, edgeThicknessByCoverage,
        minimumEdgeLength, additionalEdgeLengthPerKb,
        arrowSize);
}



void LocalAnchorGraph::writeGraphviz(
    ostream& s,
    double vertexSize,
    bool vertexSizeByCoverage,
    const string& edgeColoring,
    double edgeThickness,
    bool edgeThicknessByCoverage,
    double minimumEdgeLength,
    double additionalEdgeLengthPerKb,
    double arrowSize) const
{
    const LocalAnchorGraph& graph = *this;

    s << "digraph LocalAnchorGraph {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
        const LocalAnchorGraphVertex& vertex = graph[v];
        const AnchorId anchorId = vertex.anchorId;
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = anchors[anchorId].coverage();

        // Vertex name.
        s << "\"" << anchorIdString << "\"";

        // Begin vertex attributes.
        s << "[";

        // URL
        s << "URL=\"exploreAnchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "\"";

        // Tooltip.
        s << " tooltip=\"" << anchorIdString << " " << coverage << "\"";

        // Color.
        if(vertex.distance == 0) {
            s << " color=blue";
        } else if(vertex.distance == maxDistance) {
            s << " color=cyan";
        }

        // Size.
        const double displaySize =
            (vertexSizeByCoverage ?
            vertexSize * sqrt(0.1 * double(coverage)) :
            vertexSize
            ) / 72.;
        s << " width=" << displaySize ;
        s << " penwidth=" << 0.5 * displaySize;

        // End vertex attributes.
        s << "]";

        // End the line for this vertex.
        s << ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        const LocalAnchorGraphEdge& edge = graph[e];
        const double loss = edge.coverageLoss();

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const LocalAnchorGraphVertex& vertex1 = graph[v1];

        const AnchorId anchorId0 = vertex0.anchorId;
        const AnchorId anchorId1 = vertex1.anchorId;

        const string anchorId0String = anchorIdToString(anchorId0);
        const string anchorId1String = anchorIdToString(anchorId1);

        s << "\"" << anchorId0String << "\"->";
        s << "\"" << anchorId1String << "\"";

        // Begin edge attributes.
        s << " [";

        // URL
        s << "URL=\"exploreAnchorPair?"
            "anchorIdAString=" << HttpServer::urlEncode(anchorId0String) << "&"
            "anchorIdBString=" << HttpServer::urlEncode(anchorId1String) << "\"";

        // Tooltip.
        s << " tooltip="
            "\"" << anchorId0String << " to "
            << anchorId1String <<
            " " << edge.coverage << "/" << edge.info.common <<
            " loss " << std::fixed << std::setprecision(2) << loss <<
            " offset " << edge.info.offsetInBases << "\"";

        // Color.
        if(edgeColoring == "byCoverageLoss") {
            const double hue = (1. - loss) / 3.;
            s << " color=\"" << std::fixed << std::setprecision(2) << hue << " 1. 1.\"";
        }

        // Thickness.
        if(edgeThicknessByCoverage) {
            s << " penwidth=" << 0.1 * edgeThickness * double(edge.coverage);
        } else {
            s << " penwidth=" << edgeThickness;
        }

        // Arrow size.
        s << " arrowsize=" << arrowSize;

        // Length. Only use by fdp and neato layouts.
        const double displayLength =
            (minimumEdgeLength +
            additionalEdgeLengthPerKb * 0.001 * double(edge.info.offsetInBases)) / 72.;
        s << " len=" << displayLength;



        // End edge attributes.
        s << "]";

        // End the line for this edge.
        s << ";\n";
    }


    s << "}\n";
}

#include "mode3b-PathGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3b;

#include "fstream.hpp"


PathGraph::PathGraph(const Assembler& assembler) :
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    minPrimaryCoverage = 10;
    maxPrimaryCoverage = 25;
    minCoverage = 4;

    findVertices();
    cout << "The path graph has " << vertices.size() << " vertices. "
        "Each vertex corresponds to a primary edge of the marker graph. " <<
        "The marker graph has a total " <<
        assembler.markerGraph.edges.size() << " edges." << endl;

    findEdges();
    cout << "The path graph has " << edges.size() << " edges." << endl;

    writeGraphviz();

}



// Find out if a marker graph edge is a primary edge.
bool PathGraph::isPrimary(MarkerGraphEdgeId edgeId) const
{
    // Check coverage.
    const MarkerGraph& markerGraph = assembler.markerGraph;
    const uint64_t coverage = markerGraph.edgeCoverage(edgeId);
    if(coverage < minPrimaryCoverage) {
        return false;
    }
    if(coverage > maxPrimaryCoverage) {
        return false;
    }

    // Check for duplicate oriented reads on the edge.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeId)) {
        return false;
    }

    // Check for duplicate oriented reads on its vertices.
    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
    if(
        const auto& markers = assembler.markers;
        markerGraph.vertexHasDuplicateOrientedReadIds(edge.source, markers) or
        markerGraph.vertexHasDuplicateOrientedReadIds(edge.target, markers)) {
        return false;
    }

    // If all above checks passed, this is a primary edge.
    return true;
}



void PathGraph::findVertices()
{
    const MarkerGraph& markerGraph = assembler.markerGraph;

    vertices.clear();
    vertexTable.resize(markerGraph.edges.size());
    fill(vertexTable.begin(), vertexTable.end(), invalid<MarkerGraphEdgeId>);

    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        if(isPrimary(edgeId)) {
            vertexTable[edgeId] = vertices.size();
            vertices.push_back(edgeId);
        }
    }
}


void PathGraph::findEdges()
{
    // Store pairs (ordinal0, vertexId) for each oriented read.
    // This is indexed by OrientedReadId::getValue.
    vector < vector< pair<uint32_t, uint64_t> > > orientedReadPaths(assembler.markers.size());
    for(uint64_t vertexId=0; vertexId<vertices.size(); vertexId++) {
        const MarkerGraphEdgeId edgeId = vertices[vertexId];

        // Loop over MarkerIntervals of this primary marker graph edge.
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            orientedReadPaths[orientedReadId.getValue()].push_back(make_pair(ordinal0, vertexId));
        }
    }

    // Now follow the path of each oriented read.
    edges.clear();
    for(uint64_t i=0; i<orientedReadPaths.size(); i++) {
        auto& orientedReadPath = orientedReadPaths[i];
        sort(orientedReadPath.begin(), orientedReadPath.end(),
            OrderPairsByFirstOnly<uint32_t, uint64_t>());

        for(uint64_t j=1; j<orientedReadPath.size(); j++) {
            edges.push_back({orientedReadPath[j-1].second, orientedReadPath[j].second});
        }
    }

    // Find coverage for each edge and remove the ones with low coverage;
    edgeCoverage.clear();
    deduplicateAndCountWithThreshold(edges, edgeCoverage, minCoverage);
    edges.shrink_to_fit();

    // Write a coverage histogram.
    {
        vector<uint64_t> histogram;
        for(const uint64_t c: edgeCoverage) {
            if(c >= histogram.size()) {
                histogram.resize(c+1, 0);
            }
            ++histogram[c];
        }

        ofstream csv("PathGraphCoverageHistogram.csv");
        csv << "Coverage,Frequency\n";
        for(uint64_t c=0; c<histogram.size(); c++) {
            const uint64_t frequency = histogram[c];
            if(frequency) {
                csv << c << ",";
                csv << frequency << "\n";
            }
        }
    }

}



// Write the entire PathGraph in graphviz format.
void PathGraph::writeGraphviz() const
{
    ofstream out("PathGraph.dot");
    out << "digraph PathGraph {\n";

    for(uint64_t i=0; i<edges.size(); i++) {
        const pair<uint64_t, uint64_t>& edge = edges[i];
        const uint64_t vertexId0 = edge.first;
        const uint64_t vertexId1 = edge.second;
        const MarkerGraphEdgeId edgeId0 = vertices[vertexId0];
        const MarkerGraphEdgeId edgeId1 = vertices[vertexId1];
        out << edgeId0 << "->";
        out << edgeId1;
        out << "[tooltip=\"" << edgeId0 << "->" << edgeId1 << " " << edgeCoverage[i] << "\"];\n";
    }

    out << "}\n";
}


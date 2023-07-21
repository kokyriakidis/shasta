// Shasta.
#include "mode3b-PathGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"


PathGraph::PathGraph(const Assembler& assembler) :
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    minPrimaryCoverage = 10;
    maxPrimaryCoverage = 25;
    minCoverage = 4;
    minComponentSize = 100;

    findVertices();
    cout << "The path graph has " << verticesVector.size() << " vertices. "
        "Each vertex corresponds to a primary edge of the marker graph. " <<
        "The marker graph has a total " <<
        assembler.markerGraph.edges.size() << " edges." << endl;

    findEdges();
    cout << "The path graph has " << edgesVector.size() << " edges." << endl;

    createComponents();

    // Graphviz output.
    // writeGraphviz();    // Entire graph.
    for(uint64_t componentRank=0; componentRank<componentIndex.size(); componentRank++) {
        const uint64_t componentId = componentIndex[componentRank].first;
        const Graph& component = components[componentId];
        ofstream out("PathGraphComponent" + to_string(componentRank) + ".dot");
        component.writeGraphviz(componentRank, out);
    }
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

    verticesVector.clear();
    vertexTable.resize(markerGraph.edges.size());
    fill(vertexTable.begin(), vertexTable.end(), invalid<MarkerGraphEdgeId>);

    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        if(isPrimary(edgeId)) {
            vertexTable[edgeId] = verticesVector.size();
            verticesVector.push_back(edgeId);
        }
    }
}


void PathGraph::findEdges()
{
    // Store pairs (ordinal0, vertexId) for each oriented read.
    // This is indexed by OrientedReadId::getValue.
    vector < vector< pair<uint32_t, uint64_t> > > orientedReadPaths(assembler.markers.size());
    for(uint64_t vertexId=0; vertexId<verticesVector.size(); vertexId++) {
        const MarkerGraphEdgeId edgeId = verticesVector[vertexId];

        // Loop over MarkerIntervals of this primary marker graph edge.
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            orientedReadPaths[orientedReadId.getValue()].push_back(make_pair(ordinal0, vertexId));
        }
    }

    // Now follow the path of each oriented read.
    edgesVector.clear();
    for(uint64_t i=0; i<orientedReadPaths.size(); i++) {
        auto& orientedReadPath = orientedReadPaths[i];
        sort(orientedReadPath.begin(), orientedReadPath.end(),
            OrderPairsByFirstOnly<uint32_t, uint64_t>());

        for(uint64_t j=1; j<orientedReadPath.size(); j++) {
            edgesVector.push_back({orientedReadPath[j-1].second, orientedReadPath[j].second});
        }
    }

    // Find coverage for each edge and remove the ones with low coverage;
    edgeCoverage.clear();
    deduplicateAndCountWithThreshold(edgesVector, edgeCoverage, minCoverage);
    edgesVector.shrink_to_fit();

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

    for(uint64_t i=0; i<edgesVector.size(); i++) {
        const pair<uint64_t, uint64_t>& edge = edgesVector[i];
        const uint64_t vertexId0 = edge.first;
        const uint64_t vertexId1 = edge.second;
        const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0];
        const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1];
        out << edgeId0 << "->";
        out << edgeId1;
        out << "[tooltip=\"" << edgeId0 << "->" << edgeId1 << " " << edgeCoverage[i] << "\"];\n";
    }

    out << "}\n";
}



// Write a component of the PathGraph in graphviz format.
void PathGraph::Graph::writeGraphviz(uint64_t componentId, ostream& out) const
{
    const Graph& graph = *this;
    out << "digraph PathGraphComponent" << componentId << " {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        out << graph[v].edgeId << ";\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        out << graph[v0].edgeId << "->";
        out << graph[v1].edgeId << ";\n";
    }

    out << "}\n";
}



void PathGraph::createComponents()
{
    // Compute connected components.
    const uint64_t n = verticesVector.size();
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
    for(const auto& p: edgesVector) {
        disjointSets.union_set(p.first, p.second);
    }

    // Generate vertices of each connected component.
    components.clear();
    components.resize(n);
    for(uint64_t i=0; i<n; i++) {
        const uint64_t componentId = disjointSets.find_set(i);
        components[componentId].addVertex(verticesVector[i]);
    }

    // Create edges of each connected component.
    for(const auto& p: edgesVector) {
        const uint64_t vertexId0 = p.first;
        const uint64_t vertexId1 = p.second;
        const uint64_t componentId = disjointSets.find_set(vertexId0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(vertexId1));
        const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0];
        const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1];
        components[componentId].addEdge(edgeId0, edgeId1);
    }

    // Transitive reduction of each connected component.
    for(uint64_t componentId=0; componentId<n; componentId++) {
        Graph& component = components[componentId];
        if(num_vertices(component) > 2) {
            try {
                transitiveReduction(component);
            } catch(const exception& e) {
                cout << "Transitive reduction failed for component " <<
                    componentId << endl;
            }
        }
    }

    // Create the componentIndex.
    componentIndex.clear();
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const uint64_t componentSize = num_vertices(components[componentId]);
        if(componentSize >= minComponentSize) {
            componentIndex.push_back({componentId, componentSize});
        }
    }
    sort(componentIndex.begin(), componentIndex.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());
    cout << "Found " << componentIndex.size() << " connected components of the path graph "
        "with at least " << minComponentSize << " vertices." << endl;
}



void PathGraph::Graph::addVertex(MarkerGraphEdgeId edgeId)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    vertexMap.insert({edgeId, add_vertex(Vertex{edgeId}, *this)});
}



void PathGraph::Graph::addEdge(
    MarkerGraphEdgeId edgeId0,
    MarkerGraphEdgeId edgeId1)
{
    auto it0 = vertexMap.find(edgeId0);
    auto it1 = vertexMap.find(edgeId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    add_edge(v0, v1, *this);
}

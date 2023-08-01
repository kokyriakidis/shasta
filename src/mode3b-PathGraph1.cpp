// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "orderPairs.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>



GlobalPathGraph1::GlobalPathGraph1(const Assembler& assembler) :
    assembler(assembler)
{
    // PARAMETERS TO BE EXPOSED WHEN CODE STABILIZES.

    // The coverage range for primary vertices.
    // This controls the generation of the GlobalPathGraph1 vertices.
    const uint64_t minPrimaryCoverage = 8;
    const uint64_t maxPrimaryCoverage = 25;

    // Parameters that control generation of the GlobalPathGraph1 edges.
    const uint64_t maxDistanceInJourney = 20;
    const uint64_t minEdgeCoverage = 2;
    const double minCorrectedJaccard0 = 0.5;

    // Parameters that control the initial creation of chains.
    const double minCorrectedJaccard1 = 0.8;
    const uint64_t minComponentSize = 3;
    const uint64_t k = 1;   // For k-nn
    const uint64_t minAssembledLength = 10000;  // Minimum estimated length in bases

    // Create the GlobalPathGraph1.
    createVertices(minPrimaryCoverage, maxPrimaryCoverage);
    computeOrientedReadJourneys();
    createEdges(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard0);

    // To create initial chains, use only the best edges and
    // do K-nn and transitive reduction, then compute the longest path
    // in each connected component.
    createComponents(minCorrectedJaccard1, minComponentSize);
    knn(k);
    transitiveReduction();
    createInitialChains(minAssembledLength);
}



// Write each connected component in graphviz format.
void GlobalPathGraph1::writeGraphviz(const string& baseName) const
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        const PathGraph1& component = *components[componentRank];
        ofstream out(baseName + "-" + to_string(componentRank) + ".dot");
        component.writeGraphviz(componentRank, out);
    }
}



// Find out if a marker graph edge is a primary edge.
bool GlobalPathGraph1::isPrimary(
    MarkerGraphEdgeId edgeId,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage) const
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

#if 0
    // Check that is also is a branch edge.
    if(not isBranchEdge(edgeId)) {
        return false;
    }
#endif

    // If all above checks passed, this is a primary edge.
    return true;
}



// Find out if a marker graph edge is a branch edge.
// A marker graph edge is a branch edge if:
// - Its source vertex has more than one outgoing edge with coverage at least minPrimaryCoverage.
// OR
// - Its target vertex has more than one incoming edge with coverage at least minPrimaryCoverage.
bool GlobalPathGraph1::isBranchEdge(
    MarkerGraphEdgeId edgeId,
    uint64_t minEdgeCoverage) const
{
    // Access this marker graph edge and its vertices.
    const MarkerGraph::Edge& edge = assembler.markerGraph.edges[edgeId];
    const MarkerGraphVertexId vertexId0 = edge.source;
    const MarkerGraphVertexId vertexId1 = edge.target;

    // Check outgoing edges of vertexId0.
    const auto outgoingEdges0 = assembler.markerGraph.edgesBySource[vertexId0];
    uint64_t count0 = 0;
    for(const MarkerGraphEdgeId edgeId0: outgoingEdges0) {
        if(assembler.markerGraph.edgeCoverage(edgeId0) >= minEdgeCoverage) {
            ++count0;
        }
    }
    if(count0 > 1) {
        return true;
    }

    // Check incoming edges of vertexId1.
    const auto incomingEdges1 = assembler.markerGraph.edgesByTarget[vertexId1];
    uint64_t count1 = 0;
    for(const MarkerGraphEdgeId edgeId1: incomingEdges1) {
        if(assembler.markerGraph.edgeCoverage(edgeId1) >= minEdgeCoverage) {
            ++count1;
        }
    }
    if(count1 > 1) {
        return true;
    }

    return false;
}



void GlobalPathGraph1::createVertices(
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage)
{
    const MarkerGraph& markerGraph = assembler.markerGraph;

    vertices.clear();

    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        if(isPrimary(edgeId, minPrimaryCoverage, maxPrimaryCoverage)) {
            vertices.push_back(edgeId);
        }
    }
}



// The "journey" of each oriented read is the sequence of vertices it encounters.
// It stores pairs (ordinal0, vertexId) for each oriented read, sorted by ordinal0.
// The vertexId is the index in verticesVector.
// Indexed by OrientedReadId::getValue.
// Journeys are used to generate edges by "following the reads".
void GlobalPathGraph1::computeOrientedReadJourneys()
{
    orientedReadJourneys.clear();
    orientedReadJourneys.resize(assembler.markers.size());

    for(uint64_t vertexId=0; vertexId<vertices.size(); vertexId++) {
        const MarkerGraphEdgeId edgeId = vertices[vertexId];

        // Loop over MarkerIntervals of this primary marker graph edge.
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            orientedReadJourneys[orientedReadId.getValue()].push_back(make_pair(ordinal0, vertexId));
        }
    }

    // Now sort them, for each oriented read.
    for(auto& orientedReadJourney: orientedReadJourneys) {
        sort(orientedReadJourney.begin(), orientedReadJourney.end(),
            OrderPairsByFirstOnly<uint32_t, uint64_t>());
    }

    // Write the journeys to csv.
    ofstream csv("GlobalPathGraphJourneys.csv");
    for(ReadId readId=0; readId<assembler.markers.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = orientedReadJourneys[orientedReadId.getValue()];
            csv << orientedReadId;
            for(const auto& p: journey) {
                const uint64_t vertexId = p.second;
                const MarkerGraphEdgeId edgeId = vertices[vertexId];
                csv << "," << edgeId;
            }
            csv << "\n";
        }
    }
}



void GlobalPathGraph1::createEdges(
    uint64_t maxDistanceInJourney,
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard)
{
    // Candidate edges are pairs of vertices that appear near each other
    // in oriented read journeys.
    vector< pair<uint64_t, uint64_t> > candidateEdges;
    for(uint64_t i=0; i<orientedReadJourneys.size(); i++) {
        const auto& journey = orientedReadJourneys[i];

        for(uint64_t position0=0; position0 < journey.size(); position0++) {
            for(uint64_t distance = 1; distance <= maxDistanceInJourney; distance++) {
                const uint64_t position1 = position0 + distance;
                if(position1 >= journey.size()) {
                    break;
                }
                candidateEdges.push_back({journey[position0].second, journey[position1].second});
            }
        }
    }

    // Deduplicate the candidate edges and count the number of times
    // each of them was found. Keep only the ones that occurred at least
    // minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateEdges, coverage, minEdgeCoverage);

    // For each candidate edge, compute correctedJaccard, and if high enough
    // generate an edge.
    edges.clear();
    for(const auto& p: candidateEdges) {
        const uint64_t vertexId0 = p.first;
        const uint64_t vertexId1 = p.second;
        Edge edge;
        const MarkerGraphEdgeId edgeId0 = vertices[vertexId0];
        const MarkerGraphEdgeId edgeId1 = vertices[vertexId1];
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, edge.info));
        if(edge.info.correctedJaccard() >= minCorrectedJaccard) {
            edge.vertexId0 = vertexId0;
            edge.vertexId1 = vertexId1;
            edges.push_back(edge);
        }
    }

}



// Write the entire GlobalPathGraph in graphviz format.
void GlobalPathGraph1::writeGraphviz() const
{
    ofstream out("GlobalPathGraph.dot");
    out << "digraph GlobalPathGraph {\n";

    for(const Edge& edge: edges) {
        const MarkerGraphEdgeId edgeId0 = vertices[edge.vertexId0];
        const MarkerGraphEdgeId edgeId1 = vertices[edge.vertexId1];
        out << edgeId0 << "->";
        out << edgeId1 << ";\n";
    }

    out << "}\n";
}



void GlobalPathGraph1::createComponents(
    double minCorrectedJaccard,
    uint64_t minComponentSize)
{
    // Compute connected components.
    const uint64_t n = vertices.size();
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    for(const Edge& edge: edges) {
        if(edge.info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        disjointSets.union_set(edge.vertexId0, edge.vertexId1);
    }

    // Generate vertices of each connected component.
    vector< shared_ptr<PathGraph1> > allComponents;
    for(uint64_t componentId=0; componentId<n; componentId++) {
        allComponents.push_back(make_shared<PathGraph1>());
    }
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        const uint64_t componentId = disjointSets.find_set(vertexId);
        allComponents[componentId]->addVertex(vertexId, vertices[vertexId]);
    }

    // Create edges of each connected component.
    for(const Edge& edge: edges) {
        if(edge.info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        const uint64_t vertexId0 = edge.vertexId0;
        const uint64_t vertexId1 = edge.vertexId1;
        const uint64_t componentId = disjointSets.find_set(vertexId0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(vertexId1));
        const MarkerGraphEdgeId edgeId0 = vertices[vertexId0];
        const MarkerGraphEdgeId edgeId1 = vertices[vertexId1];
        allComponents[componentId]->addEdge(edgeId0, edgeId1, edge.info);
    }

    // Only keep the connected components that have at least
    // minComponentSize vertices.
    components.clear();
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const shared_ptr<PathGraph1> componentPointer = allComponents[componentId];
        const PathGraph1& component = *componentPointer;
        if(num_vertices(component) >= minComponentSize) {
            components.push_back(componentPointer);
        }
    }



    // Sort the pointers to connected components by decreasing size.
    class ComponentSorter {
    public:
        bool operator()(
            const shared_ptr<PathGraph1>& x,
            const shared_ptr<PathGraph1>& y)
        {
            return num_vertices(*x) > num_vertices(*y);
        }
    };
    sort(components.begin(), components.end(), ComponentSorter());
}



// K-nn of each connected component:
// for each vertex, keep only the k best outgoing and incoming edges,
// as measured by correctedJaccard of each edge.
// This can break contiguity of the connected component.
void GlobalPathGraph1::knn(uint64_t k)
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        component.knn(k);
    }
}



// Transitive reduction of each connected component.
// Ths can fail for connected components that contain cycles.
void GlobalPathGraph1::transitiveReduction()
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        shasta::transitiveReduction(component);
    }
}



void PathGraph1::addVertex(
    uint64_t vertexId,
    MarkerGraphEdgeId edgeId)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    vertexMap.insert({edgeId, add_vertex({vertexId, edgeId}, *this)});
}



void PathGraph1::addEdge(
    MarkerGraphEdgeId edgeId0,
    MarkerGraphEdgeId edgeId1,
    const MarkerGraphEdgePairInfo& info)
{
    auto it0 = vertexMap.find(edgeId0);
    auto it1 = vertexMap.find(edgeId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    add_edge(v0, v1, {info}, *this);
}



// Write a component of the GlobalPathGraph in graphviz format.
void PathGraph1::writeGraphviz(
    uint64_t componentId,
    ostream& out) const
{
    const PathGraph1& graph = *this;
    out << "digraph PathGraphComponent" << componentId << " {\n";

    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        out << graph[v].edgeId << ";\n";
    }

    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Color based on corrected Jaccard.
        const double correctedJaccard = edge.info.correctedJaccard();
        const double hue = correctedJaccard / 3.;   // 0=red, 0.5=yellow, 1=green

        out <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId <<
            " [tooltip=\"" <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId << " " <<
            std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << " " <<
            edge.info.offsetInBases;
        out << "\" color=\"" << hue << ",1,1\"];\n";
    }

    out << "}\n";
}



// For each vertex, only keep the best k outgoing and k incoming edges.
// "Best" as defined by correctedJaccard of the edges.
void PathGraph1::knn(uint64_t k)
{
    PathGraph1& graph = *this;

    // First mark all edges as not to be kept.
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        graph[e].keep = false;
    }



    // For each vertex, keep as to be kept the best k outgoing
    // and incoming edges.
    vector< pair<edge_descriptor, double> > adjacentEdges;  // With correctedJaccard.
    BGL_FORALL_VERTICES(v, graph, PathGraph1) {

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            adjacentEdges.clear();

            if(direction == 0) {
                BGL_FORALL_OUTEDGES(v, e, graph, PathGraph1) {
                    adjacentEdges.push_back({e, graph[e].info.correctedJaccard()});
                }
            } else {
                BGL_FORALL_INEDGES(v, e, graph, PathGraph1) {
                    adjacentEdges.push_back({e, graph[e].info.correctedJaccard()});
                }
            }

            // Only keep the k best.
            if(adjacentEdges.size() > k) {
                std::nth_element(
                    adjacentEdges.begin(),
                    adjacentEdges.begin() + k,
                    adjacentEdges.end(),
                    OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
                adjacentEdges.resize(k);
            }

            // Mark them as to be kept.
            for(const auto& p:adjacentEdges) {
                const edge_descriptor e = p.first;
                graph[e].keep = true;
            }
        }
    }

    // Remove edges not marked as to be kept.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        if(not graph[e].keep) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// To create initial chains:
// - Use only very strong edges (minCorrectedJaccard >= minCorrectedJaccard1).
// - Compute connected components.
// - Keep k best outgoing/incoming edges for each vertex.
// - Transitive reduction.
// This can cause contiguity breaks, which will be recovered later using
// a more complete version of the GlobalPathGraph1.
void GlobalPathGraph1::createInitialChains(uint64_t minAssembledLength)
{
    ofstream fasta("InitialChains.fasta");
    ofstream csv("InitialChains.csv");
    csv << "Rank,Vertices,Edges,Longest path length,Estimated length,Assembled length\n";

    // Compute the longest path in each component.
    vector<PathGraph1::vertex_descriptor> chain;
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        longestPath(component, chain);

        // Compute total base offset for this chain.
        uint64_t totalBaseOffset = 0;
        for(uint64_t i=1; i<chain.size(); i++) {
            const PathGraph1::vertex_descriptor v0 = chain[i-1];
            const PathGraph1::vertex_descriptor v1 = chain[i];
            PathGraph1::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, component);
            SHASTA_ASSERT(edgeExists);
            totalBaseOffset += component[e].info.offsetInBases;
        }

        if(totalBaseOffset < minAssembledLength) {
            continue;
        }



        // Generate an AssemblyPath using this chain.
        vector<MarkerGraphEdgeId> markerGraphEdgeIds;
        vector<MarkerGraphEdgePairInfo> infos;
        for(uint64_t i=0; i<chain.size(); i++) {
            markerGraphEdgeIds.push_back(component[chain[i]].edgeId);
        }
        for(uint64_t i=1; i<chain.size(); i++) {
            const PathGraph1::vertex_descriptor v0 = chain[i-1];
            const PathGraph1::vertex_descriptor v1 = chain[i];
            PathGraph1::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, component);
            infos.push_back(component[e].info);
        }
        AssemblyPath assemblyPath(assembler, markerGraphEdgeIds, infos);
        assemblyPath.writeFasta(fasta, to_string(componentRank));

        vector<Base> sequence;
        assemblyPath.getSequence(sequence);

        if(sequence.size() < minAssembledLength) {
            continue;
        }

        csv << componentRank << ",";
        csv << num_vertices(component) << ",";
        csv << num_edges(component) << ",";
        csv << chain.size() << ",";
        csv << totalBaseOffset << ",";
        csv << sequence.size() << "\n";
    }
}

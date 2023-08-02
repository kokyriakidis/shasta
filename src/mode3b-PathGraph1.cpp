// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "MurmurHash2.hpp"
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

    // Parameters that control the creation of seed chains.
    const double minCorrectedJaccard1 = 0.8;
    const uint64_t minComponentSize = 3;
    const uint64_t k = 1;   // For k-nn
    const uint64_t minEstimatedLength = 10000;  // Minimum estimated length in bases

    // Create the GlobalPathGraph1.
    createVertices(minPrimaryCoverage, maxPrimaryCoverage);
    computeOrientedReadJourneys();
    createEdges(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard0);

    // To create seed chains, use only the best edges and
    // do K-nn and transitive reduction, then compute the longest path
    // in each connected component.
    createComponents(minCorrectedJaccard1, minComponentSize);
    knn(k);
    transitiveReduction();
    createSeedChains(minEstimatedLength);
    // assembleSeedChains();
    connectSeedChains();

#if 0
    // Attempt to display the seed chains.
    createComponents(minCorrectedJaccard0, 100);
    knn(3);
    transitiveReduction();
    const bool colorChains = true;
    writeGraphviz("PathGraph", colorChains);
#endif
}



// Write each connected component in graphviz format.
void GlobalPathGraph1::writeGraphviz(const string& baseName, bool colorChains) const
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        const PathGraph1& component = *components[componentRank];
        ofstream out(baseName + "-" + to_string(componentRank) + ".dot");
        component.writeGraphviz(componentRank, colorChains, out);
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
            vertices.push_back(Vertex(edgeId));
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
        const MarkerGraphEdgeId edgeId = vertices[vertexId].edgeId;

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
                const MarkerGraphEdgeId edgeId = vertices[vertexId].edgeId;
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
        const MarkerGraphEdgeId edgeId0 = vertices[vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = vertices[vertexId1].edgeId;
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
        const MarkerGraphEdgeId edgeId0 = vertices[edge.vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = vertices[edge.vertexId1].edgeId;
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
        const Vertex& vertex = vertices[vertexId];
        const uint64_t componentId = disjointSets.find_set(vertexId);
        allComponents[componentId]->addVertex(vertexId,
            vertex.edgeId, vertex.chainId, vertex.positionInChain);
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
        const MarkerGraphEdgeId edgeId0 = vertices[vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = vertices[vertexId1].edgeId;
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
    MarkerGraphEdgeId edgeId,
    uint64_t chainId,
    uint64_t positionInChain)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    const vertex_descriptor v = add_vertex({vertexId, edgeId, chainId, positionInChain}, *this);
    vertexMap.insert({edgeId, v});
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
    bool colorChains,
    ostream& out) const
{
    const PathGraph1& graph = *this;
    out << "digraph PathGraphComponent" << componentId << " {\n";

    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        const PathGraph1Vertex& vertex = graph[v];
        out << vertex.edgeId;

        // Tooltip.
        out << " [tooltip=\"";
        out << vertex.edgeId;
        if(vertex.chainId != invalid<uint64_t>) {
            out << " " << vertex.chainId << ":" << vertex.positionInChain;
        }
        out << "\"";

        if(colorChains) {
            if(vertex.chainId == invalid<uint64_t>) {
                out << " color=black";
            } else {
                const uint32_t hue = MurmurHash2(&vertex.chainId, sizeof(vertex.chainId), 231) % 100;
                out << " color=\"" << 0.01 * double(hue) << ",1,1\"";
            }
        }
        out << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);


        out <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId <<
            " [tooltip=\"" <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId << " " <<
            std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << " " <<
            edge.info.offsetInBases << "\"";

        // Color.
        if(colorChains) {
            out << " style=invis";
#if 0
            const uint64_t chainId0 = graph[v0].chainId;
            const uint64_t chainId1 = graph[v1].chainId;
            if(
                (chainId0 != invalid<uint64_t>) and
                (chainId1 != invalid<uint64_t>) and
                (chainId0 == chainId1)) {
                const uint32_t hue = MurmurHash2(&chainId0, sizeof(chainId0), 231) % 100;
                out << " color=\"" << 0.01 * double(hue) << ",1,1\"";
            } else {
                // out << " color=\"0,0,0,0.1\" arrowhead=none";
                out << " style=invis";
            }
#endif
        } else {
            const double correctedJaccard = edge.info.correctedJaccard();
            const double hue = correctedJaccard / 3.;   // 0=red, 0.5=yellow, 1=green
            out << " color=\"" << hue << ",1,1\"";
        }
        out << "];\n";
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



// To create seed chains:
// - Use only very strong edges (minCorrectedJaccard >= minCorrectedJaccard1).
// - Compute connected components.
// - Keep k best outgoing/incoming edges for each vertex.
// - Transitive reduction.
// This can cause contiguity breaks, which will be recovered later using
// a more complete version of the GlobalPathGraph1.
void GlobalPathGraph1::createSeedChains(
    uint64_t minEstimatedLength)
{
    ofstream fasta("SeedChains.fasta");
    ofstream csv("SeedChains.csv");
    csv << "ChainId,Vertices,Edges,Longest path length,Estimated length\n";

    seedChains.clear();

    // Loop over connected components.
    vector<PathGraph1::vertex_descriptor> chainVertices;
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];

        // Compute the longest path in this component.
        longestPath(component, chainVertices);

        // Use this longest path to create a tentative Chain.
        Chain chain;
        for(const PathGraph1::vertex_descriptor v: chainVertices) {
            chain.vertexIds.push_back(component[v].vertexId);
        }
        for(uint64_t i=1; i<chainVertices.size(); i++) {
            const PathGraph1::vertex_descriptor v0 = chainVertices[i-1];
            const PathGraph1::vertex_descriptor v1 = chainVertices[i];
            PathGraph1::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, component);
            SHASTA_ASSERT(edgeExists);
            chain.infos.push_back(component[e].info);
        }

        // Compute total base offset for this chain.
        uint64_t totalBaseOffset = chain.totalOffset();

        // If too short, discard it.
        if(totalBaseOffset < minEstimatedLength) {
            continue;
        }


        // Store this chain as a new seed chain.
        const uint64_t chainId = seedChains.size();
        seedChains.push_back(chain);
        for(uint64_t position=0; position<chain.vertexIds.size(); position++) {
            const uint64_t vertexId = chain.vertexIds[position];
            Vertex& vertex = vertices[vertexId];
            vertex.chainId = chainId;
            vertex.positionInChain = position;
        }

        csv << chainId << ",";
        csv << num_vertices(component) << ",";
        csv << num_edges(component) << ",";
        csv << chain.vertexIds.size() << ",";
        csv << totalBaseOffset << ",";
        csv << "\n";

    }
}



// Compute total base offset for this chain.
uint64_t GlobalPathGraph1::Chain::totalOffset() const
{
    uint64_t totalBaseOffset = 0;
    for(const MarkerGraphEdgePairInfo& info: infos) {
        totalBaseOffset += info.offsetInBases;
    }
    return totalBaseOffset;
}



void GlobalPathGraph1::assembleSeedChains() const
{
    ofstream fasta("SeedChains.fasta");

    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        const Chain& chain = seedChains[chainId];

        // Generate an AssemblyPath using this chain.
        vector<MarkerGraphEdgeId> markerGraphEdgeIds;
        for(const uint64_t vertexId: chain.vertexIds) {
            markerGraphEdgeIds.push_back(vertices[vertexId].edgeId);
        }
        AssemblyPath assemblyPath(assembler, markerGraphEdgeIds, chain.infos);
        assemblyPath.writeFasta(fasta, to_string(chainId));
    }
}



void GlobalPathGraph1::connectSeedChains()
{

    // Compute the journey of each oriented read on the seed chains,
    // that is the sequence of seed chain visited by each oriented read.
    vector< vector<uint64_t> > orientedReadSeedChainsJourneys(orientedReadJourneys.size());

    // Loop over all oriented reads.
    for(uint64_t i=0; i<orientedReadJourneys.size(); i++) {
        const vector< pair<uint32_t, uint64_t> >& journey = orientedReadJourneys[i];
        vector<uint64_t>& seedChainsJourney = orientedReadSeedChainsJourneys[i];

        // Loop over the journey of this oriented read.
        for(const auto& p: journey) {

            // Find the chain.
            const uint64_t vertexId = p.second;
            SHASTA_ASSERT(vertexId < vertices.size());
            const Vertex& vertex = vertices[vertexId];
            const uint64_t chainId = vertex.chainId;

            // If no chain here, do nothing.
            if(chainId == invalid<uint64_t>) {
                continue;
            }
            SHASTA_ASSERT(chainId < seedChains.size());

            // Append the chain to the seed chain journey, if different from the last one.
            if(seedChainsJourney.empty() or chainId != seedChainsJourney.back()) {
                seedChainsJourney.push_back(chainId);
            }
        }
    }


    // Count adjacent pairs of chains in journeys.
    vector< pair<uint64_t, uint64_t> > chainPairs;
    for(uint64_t i=0; i<orientedReadJourneys.size(); i++) {
        const vector<uint64_t>& seedChainsJourney = orientedReadSeedChainsJourneys[i];
        for(uint64_t j=1; j<seedChainsJourney.size(); j++) {
            chainPairs.push_back({seedChainsJourney[j-1], seedChainsJourney[j]});
        }
    }
    vector<uint64_t> coverage;
    deduplicateAndCount(chainPairs, coverage);

    cout << "Found " << chainPairs.size() << " chain pairs for " <<
        seedChains.size() << " seed chains." << endl;
    ofstream out("SeedChains.dot");
    out << "digraph SeedChains{\n";
    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        out << chainId << " [label=\"C" << chainId << "\\n" << seedChains[chainId].totalOffset() << "\"];\n";
    }
    for(uint64_t i=0; i<chainPairs.size(); i++) {
        const auto& p = chainPairs[i];
        if(coverage[i] < 6) {
            continue;
        }
        out << p.first << "->" << p.second <<
            "[label=\"" << coverage[i] << "\"];\n";
    }
    out << "}\n";


}

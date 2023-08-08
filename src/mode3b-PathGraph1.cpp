// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>
#include <numeric>



GlobalPathGraph1::GlobalPathGraph1(const Assembler& assembler) :
    assembler(assembler)
{
    // PARAMETERS TO BE EXPOSED WHEN CODE STABILIZES
    // ARE DEFINED BEFORE EACH PHASE THAT USES THEM.

    // Create GlobalPathGraph1 vertices.
    {
        const uint64_t minPrimaryCoverage = 8;
        const uint64_t maxPrimaryCoverage = 25;
        cout << timestamp << "Creating vertices." << endl;
        createVertices(minPrimaryCoverage, maxPrimaryCoverage);
        cout << timestamp << "Found " << vertices.size() << " vertices." << endl;
        computeOrientedReadJourneys();
    }

    // Create GlobalPathGraph1 edges.
    {
        const uint64_t maxDistanceInJourney = 20;
        const uint64_t minEdgeCoverage = 3;
        const double minCorrectedJaccard = 0.8;
        cout << timestamp << "Creating edges." << endl;
        createEdges0(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard);
        // createEdges1(minEdgeCoverage, minCorrectedJaccard); // Unlimited distance (non-local).
        cout << timestamp << "Found " << edges.size() << " edges." << endl;
    }

#if 0
    // Initial display.
    {
        const double minCorrectedJaccard = 0.8;
        const uint64_t minComponentSize = 100;
        createComponents(minCorrectedJaccard, minComponentSize);

        // Clean it up a bit for easier display.
        const uint64_t k = 3;
        knn(k);
        transitiveReduction();

        writeGraphviz("PathGraph", minCorrectedJaccard, 1.);
    }
#endif

    // To create seed chains, use only the best edges and
    // do K-nn and transitive reduction, then compute the longest path
    // in each connected component.
    {

        // Compute connected components using only the best edges.
        const double minCorrectedJaccard = 0.8;
        const uint64_t minComponentSize = 3;
        createComponents(minCorrectedJaccard, minComponentSize);

        // K-nn.
        const uint64_t k = 3;
        knn(k);

        // Transitive reduction is not really needed because it does not
        // change the longest path, but it can be useful if we want to display
        // the connected components at this stage.
        // transitiveReduction();

        // Create the seed chains.
        // Only keep the ones that are long enough.
        // Minimum estimated length is in bases.
        const uint64_t minEstimatedLength = 10000;
        createSeedChains(minEstimatedLength);
        cout << "Found " << seedChains.size() << " seed chains." << endl;
        writeSeedChainsDetails();
        writeSeedChainsStatistics();
        // assembleSeedChains();
    }



    // Connect seed chains.
    vector<ChainConnector> connectors;
    {
        const uint64_t minEdgeCoverage = 4;
        const double minCorrectedJaccard = 0.6;
        connectSeedChains1(minEdgeCoverage, minCorrectedJaccard, connectors);

        // Write the chain connectors.
        ofstream csv("ChainConnectors.csv");
        csv << "ChainId0,ChainId1,Position,EdgeId0,EdgeId1,Common,Offset,J,J'\n";
        for(const ChainConnector& connector: connectors) {
            for(uint64_t position=0; position<connector.infos.size(); position++) {
                const MarkerGraphEdgePairInfo& info = connector.infos[position];
                const uint64_t vertexId0 = connector.vertexIds[position];
                const uint64_t vertexId1 = connector.vertexIds[position + 1];
                const MarkerGraphEdgeId edgeId0 = vertices[vertexId0].edgeId;
                const MarkerGraphEdgeId edgeId1 = vertices[vertexId1].edgeId;
                csv <<
                    connector.chainId0 << "," <<
                    connector.chainId1 << "," <<
                    position << "," <<
                    edgeId0 << "," <<
                    edgeId1 << "," <<
                    info.common << "," <<
                    info.offsetInBases << "," <<
                    info.jaccard() << "," <<
                    info.correctedJaccard() << "\n";
            }
        }
    }


    // Use the ChainConnectors to stitch together the seed chains.
    {
        const uint64_t minComponentSize = 3;
        stitchSeedChains(connectors, minComponentSize);
    }

}



// Write each connected component in graphviz format.
void GlobalPathGraph1::writeGraphviz(
    const string& baseName,
    double redJ,
    double greenJ) const
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        const PathGraph1& component = *components[componentRank];
        ofstream out(baseName + "-" + to_string(componentRank) + ".dot");
        component.writeGraphviz(vertices, componentRank, redJ, greenJ,out);
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
            vertices.push_back(GlobalPathGraph1Vertex(edgeId));
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



    // Store journey information in the vertices.
    for(ReadId readId=0; readId<assembler.markers.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = orientedReadJourneys[orientedReadId.getValue()];

            for(uint64_t position=0; position<journey.size(); position++) {
                const auto& p = journey[position];
                const uint64_t vertexId = p.second;
                vertices[vertexId].journeyInfoItems.push_back({orientedReadId, position});
            }
        }
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



void GlobalPathGraph1::createEdges0(
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
    cout << timestamp << "Found " << candidateEdges.size() << " candidate edges, including duplicates." << endl;

    // Deduplicate the candidate edges and count the number of times
    // each of them was found. Keep only the ones that occurred at least
    // minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateEdges, coverage, minEdgeCoverage);
    cout << timestamp << "After deduplication, there are " << candidateEdges.size() << " candidate edges." << endl;
    candidateEdges.shrink_to_fit();

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



void GlobalPathGraph1::createEdges1(
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard)
{
    // Vector to store the children of a single vertex.
    // The first element of each pair if the vertexId of the child.
    vector< pair<uint64_t, MarkerGraphEdgePairInfo> > children;

    // Loop over all vertices.
    Edge edge;
    for(edge.vertexId0=0; edge.vertexId0<vertices.size(); edge.vertexId0++) {
        if((edge.vertexId0 % 10000) == 0) {
            cout << edge.vertexId0 << "/" << vertices.size() << endl;
        }

        // Find its children.
        findChildren(edge.vertexId0, minEdgeCoverage, minCorrectedJaccard, children);

        // Store them.
        for(const auto& p: children) {
            edge.vertexId1 = p.first;
            edge.info = p.second;
            edges.push_back(edge);
        }
    }
}



// Find children edges of vertexId0.
// The first element of each pair of the children vector
// is the vertexId of the child vertex.
void GlobalPathGraph1::findChildren(
    uint64_t vertexId0,
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard,
    vector< pair<uint64_t, MarkerGraphEdgePairInfo> >& children)
{
    const GlobalPathGraph1Vertex& vertex0 = vertices[vertexId0];
    const MarkerGraphEdgeId edgeId0 = vertex0.edgeId;

    // Find vertices encountered later on journeys that go through here.
    vector<uint64_t> candidateChildren;
    for(const auto& journeyInfoItem: vertex0.journeyInfoItems) {
        const OrientedReadId orientedReadId = journeyInfoItem.orientedReadId;
        const auto& journey = orientedReadJourneys[orientedReadId.getValue()];
        for(uint64_t position = journeyInfoItem.positionInJourney+1;
            position < journey.size(); position++) {
            candidateChildren.push_back(journey[position].second);
        }
    }

    // Count the candidate children and only keep the ones
    // that occurred at least minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateChildren, coverage, minEdgeCoverage);

    // Store the ones that have sufficient correctedJaccard.
    children.clear();
    MarkerGraphEdgePairInfo info;
    for(const uint64_t vertexId1: candidateChildren) {
        const GlobalPathGraph1Vertex& vertex1 = vertices[vertexId1];
        const MarkerGraphEdgeId edgeId1 = vertex1.edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
        if(info.correctedJaccard() >= minCorrectedJaccard) {
            children.push_back({vertexId1, info});
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
        const GlobalPathGraph1Vertex& vertex = vertices[vertexId];
        const uint64_t componentId = disjointSets.find_set(vertexId);
        allComponents[componentId]->addVertex(vertexId, vertex.edgeId);
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



void GlobalPathGraph1::kClosest(uint64_t k)
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        component.kClosest(k);
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
    const vertex_descriptor v = add_vertex({vertexId, edgeId}, *this);
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



// Write a PathGraph1 in graphviz format.
void PathGraph1::writeGraphviz(
    const vector<GlobalPathGraph1Vertex>& globalVertices,
    uint64_t componentId,
    double redJ,
    double greenJ,
    ostream& out) const
{
    const PathGraph1& graph = *this;
    out << "digraph PathGraphComponent" << componentId << " {\n";

    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        const PathGraph1Vertex& vertex = graph[v];
        const GlobalPathGraph1Vertex& globalVertex = globalVertices[vertex.vertexId];
        out << vertex.edgeId;

        // Tooltip.
        out << " [tooltip=\"";
        out << vertex.edgeId;
        if(globalVertex.chainId != invalid<uint64_t>) {
            out << " " << globalVertex.chainId << ":" << globalVertex.positionInChain;
        }
        out << "\"";

        if(globalVertex.chainId != invalid<uint64_t>) {
            if(globalVertex.isFirstInChain) {
                out << " color=blue";
            } else if(globalVertex.isLastInChain) {
                out << " color=orange";
            } else {
                // const uint32_t hue = MurmurHash2(&vertex.chainId, sizeof(vertex.chainId), 231) % 100;
                // out << " color=\"" << 0.01 * double(hue) << ",0.4,1\"";
                out << " color=cyan";
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
        const double correctedJaccard = edge.info.correctedJaccard();
        if(correctedJaccard <= redJ) {
            out << " color=red";
        } else if(correctedJaccard >= greenJ) {
            out << " color=green";
        } else {
            const double hue = (correctedJaccard - redJ) / (3. * (greenJ - redJ));
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



// Keep the k closest outgoing and incoming edges for each vertex.
void PathGraph1::kClosest(uint64_t k)
{
    PathGraph1& graph = *this;

    // First mark all edges as not to be kept.
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        graph[e].keep = false;
    }



    // For each vertex, mark as to be kept the closest k outgoing
    // and incoming edges.
    vector< pair<edge_descriptor, uint64_t> > adjacentEdges;  // With estimated offset.
    BGL_FORALL_VERTICES(v, graph, PathGraph1) {

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            adjacentEdges.clear();

            if(direction == 0) {
                BGL_FORALL_OUTEDGES(v, e, graph, PathGraph1) {
                    adjacentEdges.push_back({e, graph[e].info.offsetInBases});
                }
            } else {
                BGL_FORALL_INEDGES(v, e, graph, PathGraph1) {
                    adjacentEdges.push_back({e, graph[e].info.offsetInBases});
                }
            }

            // Only keep the k best.
            if(adjacentEdges.size() > k) {
                std::nth_element(
                    adjacentEdges.begin(),
                    adjacentEdges.begin() + k,
                    adjacentEdges.end(),
                    OrderPairsBySecondOnly<edge_descriptor, uint64_t>());
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
    csv << "ChainId,Vertices,Edges,First,Last,Longest path length,Estimated length\n";

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
            GlobalPathGraph1Vertex& vertex = vertices[vertexId];
            vertex.chainId = chainId;
            vertex.positionInChain = position;
            vertex.isFirstInChain = (position == 0);
            vertex.isLastInChain = (position == chain.vertexIds.size() - 1);
        }

        csv << chainId << ",";
        csv << num_vertices(component) << ",";
        csv << num_edges(component) << ",";
        csv << vertices[chain.vertexIds.front()].edgeId << ",";
        csv << vertices[chain.vertexIds.back()].edgeId << ",";
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



void GlobalPathGraph1::writeSeedChainsDetails() const
{
    ofstream csv("SeedChainsDetails.csv");
    csv << "Chain,Position,Marker graph edge,Corrected Jaccard to next\n";

    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        const Chain& chain = seedChains[chainId];

        for(uint64_t position=0; position<chain.vertexIds.size(); position++) {
            const uint64_t vertexId = chain.vertexIds[position];

            csv <<
                chainId << "," <<
                position << "," <<
                vertices[vertexId].edgeId << ",";
            if(position != chain.vertexIds.size()-1) {
                csv << chain.infos[position].correctedJaccard() << ",";
            }
            csv << "\n";
        }

    }

}



void GlobalPathGraph1::writeSeedChainsStatistics() const
{
    const uint64_t minChainLength = 200000;
    const uint64_t binWidth = 5000;

    class Bin {
    public:
        uint64_t n;
        uint64_t commonSum = 0;
        uint64_t commonSum2 = 0;
        void add(uint64_t common)
        {
            n++;
            commonSum += common;
            commonSum2 += common * common;
        }
        double commonAverage() const
        {
            return double(commonSum) / double(n);
        }
        double commonSigma() const
        {
            return sqrt(double(commonSum2) / double(n) - commonAverage() * commonAverage());
        }
    };
    vector<Bin> bins;


    MarkerGraphEdgePairInfo info;
    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        const Chain& chain = seedChains[chainId];
        if(chain.totalOffset() < minChainLength) {
            continue;
        }

        // Choose a subset of the chain vertices to use for statistics.
        vector<uint64_t> vertexIds;
        const uint64_t desiredVertexCount = 100;
        if(desiredVertexCount < chain.vertexIds.size()) {
            const double ratio = double(desiredVertexCount) / double(chain.vertexIds.size());
            const uint32_t hashThreshold = uint32_t(ratio * double(std::numeric_limits<uint32_t>::max()));
            for(const uint64_t vertexId: chain.vertexIds) {
                if(MurmurHash2(&vertexId, sizeof(vertexId), 231) < hashThreshold) {
                    vertexIds.push_back(vertexId);
                }
            }
        } else {
            vertexIds = chain.vertexIds;
        }

        // Now loop over pairs of the vertices we selected.
        for(uint64_t i0=0; i0<vertexIds.size(); i0++) {
            const uint64_t vertexId0 = vertexIds[i0];
            const MarkerGraphEdgeId edgeId0 = vertices[vertexId0].edgeId;
            for(uint64_t i1=i0+1; i1<vertexIds.size(); i1++) {
                const uint64_t vertexId1 = vertexIds[i1];
                const MarkerGraphEdgeId edgeId1 = vertices[vertexId1].edgeId;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
                if(info.common == 0) {
                    continue;
                }

                // Access the bin for this offset.
                const uint64_t binId = info.offsetInBases / binWidth;
                if(binId >= bins.size()) {
                    bins.resize(binId + 1);
                }
                bins[binId].add(info.common);
            }
        }
    }

    ofstream csv("SeedChainsStatistics.csv");
    csv << "BinMin,BinMid,BinMax,Coverage,Coverage sigma\n";
    for(uint64_t binId=0; binId<bins.size(); binId++) {
        const Bin& bin = bins[binId];
        const uint64_t binMin = binId * binWidth;
        const uint64_t binMax = binMin + binWidth;
        const uint64_t binMid = (binMin + binMax) / 2;
        csv <<
            binMin << "," <<
            binMid << "," <<
            binMax << "," <<
            bin.commonAverage() << "," <<
            bin.commonSigma() << "\n";
    }
}



// Connect seed chains by following reads on the graph.
void GlobalPathGraph1::connectSeedChains0()
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
            const GlobalPathGraph1Vertex& vertex = vertices[vertexId];
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

        const uint64_t chainId0 = p.first;
        const uint64_t chainId1 = p.second;
        const Chain& chain0 = seedChains[chainId0];
        const Chain& chain1 = seedChains[chainId1];

        const MarkerGraphEdgeId edgeId0 = vertices[chain0.vertexIds.back()].edgeId;
        const MarkerGraphEdgeId edgeId1 = vertices[chain1.vertexIds.front()].edgeId;
        MarkerGraphEdgePairInfo info;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));

        out << std::fixed << std::setprecision(2);
        out << p.first << "->" << p.second <<
            "[label=\"" << coverage[i] << "\\n" << info.correctedJaccard() << "\"];\n";
    }
    out << "}\n";


}



// To connect a chain to the next, use a shortest path algorithm (Dijkstra algorithm)
// starting at the last vertex of the chain. The shortest path does not run
// on the GlobalGraph1 or one of its connected components.
// Instead, it runs on the implicit graph defined by findChildren,
// with the length of an edge defined by offsetInBases.
// Uses a boost::multi_index_container as the data structure in
// Dijkstra's algorithm.
// See:
// https://en.wikipedia.org/wiki/Dijkstra's_algorithm
// https://www.boost.org/doc/libs/1_82_0/libs/multi_index/doc/tutorial/basics.html#multiple_sort
void GlobalPathGraph1::connectSeedChains1(
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard,
    vector<ChainConnector>& connectors
    )
{
    using boost::multi_index_container;
    using boost::multi_index::indexed_by;
    using boost::multi_index::member;
    using boost::multi_index::ordered_unique;
    using boost::multi_index::ordered_non_unique;



    // Information about a vertex stored during the Dijkstra algorithm.
    class VertexInfo {
    public:
        uint64_t vertexId;

        // The sum of offsets starting at the start vertex.
        uint64_t distance;

        // The vertex that achieved that distance
        // and the corresponding MarkerGraphEdgePairInfo.
        uint64_t parent;
        MarkerGraphEdgePairInfo info;
    };



    // The container used to store VertexInfos,
    // with the two indices to support the required operations.
    using VertexContainer = multi_index_container<VertexInfo,
        indexed_by <
        ordered_unique< member<VertexInfo, uint64_t, &VertexInfo::vertexId> >,
        ordered_non_unique< member<VertexInfo, uint64_t, &VertexInfo::distance> >
        > >;

    // Last argument for findChildren below, defined here to reduce
    // memory allocation activity.
    vector< pair<uint64_t, MarkerGraphEdgePairInfo> > children;

    ofstream out("SeedChains.dot");
    out << "digraph SeedChains {\n";

    connectors.clear();

    // Loop over chains.
    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {

        // cout << "Working on chain " << chainId << endl;
        const Chain& chain = seedChains[chainId];
        const uint64_t chainLastVertexId = chain.vertexIds.back();
        // cout << "Last vertex of chain " <<  vertices[chainLastVertexId].edgeId << " is starting vertex for search." << endl;

        // Start the Dijkstra algorithm with this as the only seen vertex.
        VertexContainer seen;
        const auto& seenById = seen.get<0>();
        auto& seenByDistance = seen.get<1>();   // Not const so we can erase by distance.
        seen.insert({chainLastVertexId, 0, invalid<uint64_t>, MarkerGraphEdgePairInfo()});

        VertexContainer visited;
        const auto& visitedById = visited.get<0>();

        // Dijkstra loop.
        while(not seen.empty()) {

            // Visit the vertex with the lowest distance.
            auto it = seenByDistance.begin();
            VertexInfo v0 = *it;
            visited.insert(*it);
            seenByDistance.erase(it);
            // cout << "Visiting " << vertices[v0.vertexId].edgeId << endl;

            // If it belongs to a different chain, we are done.
            // Create a new connector.
            const uint64_t chainId0 = vertices[v0.vertexId].chainId;
            if(chainId0 != invalid<uint64_t> and chainId0 != chainId) {
                // cout << vertices[v0.vertexId].edgeId << " is on chain " << chainId0 <<
                //     " at distance " << v0.distance << endl;
                out << chainId << "->" << chainId0 << ";\n";

                // To construct the connector between these two chains,
                // walk back.
                connectors.resize(connectors.size() + 1);
                ChainConnector& connector = connectors.back();
                connector.chainId0 = chainId;
                connector.chainId1 = chainId0;
                uint64_t vertexId = v0.vertexId;
                while(true) {
                    connector.vertexIds.push_back(vertexId);
                    if(vertexId == chainLastVertexId) {
                        break;
                    }
                    const auto it = visitedById.find(vertexId);
                    SHASTA_ASSERT(it != visitedById.end());
                    const VertexInfo& vertexInfo = *it;
                    connector.infos.push_back(vertexInfo.info);
                    vertexId = vertexInfo.parent;
                }
                connector.reverse();
                break;
            }

            // Find its children.
            findChildren(v0.vertexId, minEdgeCoverage, minCorrectedJaccard, children);

            // Loop over the unvisited children.
            for(const auto& child: children) {
                const uint64_t vertexId1 = child.first;
                // cout << "Found child " << vertices[vertexId1].edgeId << endl;
                if(visitedById.find(vertexId1) != visitedById.end()) {
                    // cout << "Already visited." << endl;
                    continue;
                }

                const MarkerGraphEdgePairInfo& info = child.second;
                SHASTA_ASSERT(info.offsetInBases > 0);

                // Update the tentative distance of the child,
                // adding it to the seenVertices if not present.
                auto it = seenById.find(vertexId1);
                const uint64_t distance1 = v0.distance + info.offsetInBases;
                if(it == seenById.end()) {
                    seen.insert({vertexId1, distance1, v0.vertexId, info});
                } else {
                    VertexInfo seenVertex1 = *it;
                    if(distance1 < seenVertex1.distance) {
                        seenVertex1.distance = distance1;
                        seenVertex1.parent = v0.vertexId;
                        seenVertex1.info = info;
                        seen.replace(it, seenVertex1);
                    }
                }

            }
        }

    }
    out << "}\n";
}



void GlobalPathGraph1::ChainConnector::reverse()
{
    std::reverse(vertexIds.begin(), vertexIds.end());
    std::reverse(infos.begin(), infos.end());
}



// Use the ChainConnectors to stitch together the seed chains.
void GlobalPathGraph1::stitchSeedChains(
    const vector<ChainConnector>& connectors,
    uint64_t minComponentSize)
{

    // Create a PathGraph1 containing the chains and the connectors.
    PathGraph1 graph;



    // Add the vertices and edges of the chains.
    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        // cout << "Adding vertices and edges for chain " << chainId << endl;
        const Chain& chain = seedChains[chainId];

        // Add the vertices of this chain.
        for(uint64_t position=0; position<chain.vertexIds.size(); position++) {
            const uint64_t vertexId = chain.vertexIds[position];
            const MarkerGraphEdgeId edgeId = vertices[vertexId].edgeId;
            graph.addVertex(vertexId, edgeId);
            // cout << "Adding vertex for " << edgeId << endl;
        }

        // Add the edges of this chain
        for(uint64_t position=0; position<chain.vertexIds.size()-1; position++) {
            const uint64_t vertexId0 = chain.vertexIds[position];
            const uint64_t vertexId1 = chain.vertexIds[position + 1];
            const GlobalPathGraph1Vertex& vertex0 = vertices[vertexId0];
            const GlobalPathGraph1Vertex& vertex1 = vertices[vertexId1];
            graph.addEdge(vertex0.edgeId, vertex1.edgeId, chain.infos[position]);
        }
    }



    // Add the vertices and edges of the connectors.
    for(const ChainConnector& connector: connectors) {
        // cout << "Adding vertices and edges for connector between chains " <<
        //    connector.chainId0 << " " << connector.chainId1 << endl;

        // Add the vertices of this connector, except for the first and last
        // which are part of chains.
        for(uint64_t position=1; position<connector.vertexIds.size()-1; position++) {
            const uint64_t vertexId = connector.vertexIds[position];
            const MarkerGraphEdgeId edgeId = vertices[vertexId].edgeId;
            // cout << "Adding vertex for " << edgeId << endl;
            // The vertex could have already been added as part of another connector.
            if(not graph.vertexMap.contains(edgeId)) {
                graph.addVertex(vertexId, edgeId);
            }
       }

        // Add the edges of this connector.
        for(uint64_t position=0; position<connector.infos.size(); position++) {
            const uint64_t vertexId0 = connector.vertexIds[position];
            const uint64_t vertexId1 = connector.vertexIds[position + 1];
            const GlobalPathGraph1Vertex& vertex0 = vertices[vertexId0];
            const GlobalPathGraph1Vertex& vertex1 = vertices[vertexId1];
            // cout << "Adding edge for " << vertex0.edgeId << " " << vertex1.edgeId << endl;
            graph.addEdge(vertex0.edgeId, vertex1.edgeId, connector.infos[position]);
        }

    }

    ofstream out("StitchedSeedChains.dot");
    graph.writeGraphviz(vertices, 0, 0.5, 1., out);

    cout << "The stitched graph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;


    // Compute connected components of the stitched graph.
    const vector< shared_ptr<PathGraph1> > componentPointers =
        graph.createConnectedComponents(minComponentSize);
    cout << "The stitched graph has " << componentPointers.size() <<
        " connected components." << endl;
}



// Create the connected components of this PathGraph1,
// without changing the PathGraph1 itself.
vector< shared_ptr<PathGraph1> > PathGraph1::createConnectedComponents(
    uint64_t minComponentSize) const
{
    const PathGraph1& graph = *this;

    // Compute connected components.
    // We can't use boost::connected_components because it only works
    // for undirected graphs.
    const uint64_t n = num_vertices(graph);
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1::vertex_descriptor v0 = source(e, graph);
        const PathGraph1::vertex_descriptor v1 = target(e, graph);
        disjointSets.union_set(v0, v1);
    }


    // Gather the vertices in each connected component.
    vector< shared_ptr<PathGraph1> > allComponentPointers(num_vertices(graph));
    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        const PathGraph1Vertex& vertex = graph[v];
        const uint64_t componentId = disjointSets.find_set(v);
        shared_ptr<PathGraph1>& componentPointer = allComponentPointers[componentId];
        if(not componentPointer) {
            componentPointer = make_shared<PathGraph1>();
        }
        PathGraph1& component = *componentPointer;
        component.addVertex(
            vertex.vertexId,
            vertex.edgeId);
    }


    // Gather the edges in each connected component.
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1::vertex_descriptor v0 = source(e, graph);
        const PathGraph1::vertex_descriptor v1 = target(e, graph);
        const uint64_t edgeId0 = graph[v0].edgeId;
        const uint64_t edgeId1 = graph[v1].edgeId;
        const uint64_t componentId = disjointSets.find_set(v0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(v1));
        shared_ptr<PathGraph1>& componentPointer = allComponentPointers[componentId];
        SHASTA_ASSERT(componentPointer);
        PathGraph1& component = *componentPointer;
        component.addEdge(
            edgeId0,
            edgeId1,
            graph[e].info);
    }



    // Keep only the components with at least minComponentSize vertices
    // and sort them by size.
    vector< pair<shared_ptr<PathGraph1>, uint64_t> > componentPointersWithSizes;
    for(const shared_ptr<PathGraph1>& p: allComponentPointers) {
        if(p) {
            const uint64_t componentSize = num_vertices(*p);
            if(componentSize >= minComponentSize) {
                componentPointersWithSizes.push_back({p, componentSize});
            }
        }
    }
    sort(componentPointersWithSizes.begin(), componentPointersWithSizes.end(),
        OrderPairsBySecondOnlyGreater<shared_ptr<PathGraph1>, uint64_t>());


    // For now return all components, including the empty ones.
    // But we want to remove the small ones and sort them by size.
    vector< shared_ptr<PathGraph1> > componentPointers;
    for(const auto& p: componentPointersWithSizes) {
        componentPointers.push_back(p.first);
    }
    return componentPointers;
}

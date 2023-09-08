// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
#include "localTransitiveReduction.hpp"
#include "longestPath.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>
#include <numeric>



void GlobalPathGraph1::assemble(const Assembler& assembler)
{
    assemble1(assembler);
}



void GlobalPathGraph1::assemble0(const Assembler& assembler)
    {
    // PARAMETERS TO BE EXPOSED WHEN CODE STABILIZES
    // ARE DEFINED BEFORE EACH PHASE THAT USES THEM.

    // Create a first version of the GlobalPathGraph1 to find seed chains.
    vector<Chain> seedChains;
    {
        const uint64_t minPrimaryCoverage = 8;
        const uint64_t maxPrimaryCoverage = 25;
        const uint64_t maxDistanceInJourney = 20;
        const uint64_t minEdgeCoverage = 3;
        const double minCorrectedJaccard = 0.8;
        const uint64_t minComponentSize = 3;
        const uint64_t k = 3;
        const uint64_t minEstimatedLength = 10000;

        GlobalPathGraph1 graph(assembler);
        graph.createVertices(minPrimaryCoverage, maxPrimaryCoverage);
        graph.computeOrientedReadJourneys();
        graph.createEdges0(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard);
        graph.createComponents(minCorrectedJaccard, minComponentSize);
        graph.knn(k);

        // Transitive reduction is not really needed because it does not
        // change the longest path in each connected component,
        // but it can be useful if we want to display
        // the connected components at this stage.
        // graph.transitiveReduction();

        // Output each connected component in a separate file.
        // This can create a large number of files.
        // graph.writeComponentsGraphviz("PathGraphA", minCorrectedJaccard, 1.);


        graph.createChainsFromComponents(minEstimatedLength, graph.seedChains);
        cout << "Found " << graph.seedChains.size() << " seed chains." << endl;
        graph.writeSeedChains();
        graph.writeSeedChainsDetails();
        graph.writeSeedChainsStatistics();
        seedChains = graph.seedChains;

#if 0
        ofstream fasta("SeedChains.fasta");
        graph.assembleChains(graph.seedChains, fasta, "SeedChain-");
#endif
    }



    // Create a new GlobalPathGraph1 with less strict criteria
    // for vertex creation. This will be used to connect
    // the seed chains we found.
    {
        const uint64_t minPrimaryCoverage = 8;
        const uint64_t maxPrimaryCoverage = 25;
        const uint64_t minEdgeCoverage = 4;
        const double minCorrectedJaccard = 0.6;
        const uint64_t minComponentSize = 3;
        const uint64_t minEstimatedLength = 10000;
        const uint64_t connectSeedChainsMethod = 1;

        GlobalPathGraph1 graph(assembler);
        graph.createVertices(minPrimaryCoverage, maxPrimaryCoverage);
        graph.computeOrientedReadJourneys();

        // Store in this new graph the seed chains we found earlier.
        graph.storeSeedChains(seedChains);

        // If using connectSeedChains2, we also generate the edges
        // by just incrementally following the reads and accepting all edges.
        if(connectSeedChainsMethod == 2) {
            const uint64_t maxDistanceInJourney = 1;
            const uint64_t minEdgeCoverage = 1;
            const double minCorrectedJaccard = 0.;
            graph.createEdges0(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard);
        }



        // Compute connectors to join pairs of seed chains.
        vector<ChainConnector> connectors;
        if(connectSeedChainsMethod == 1) {
            graph.connectSeedChains1(minEdgeCoverage, minCorrectedJaccard, connectors);
        } else if(connectSeedChainsMethod == 2) {

            // Compute connected components.
            const double minCorrectedJaccardForComponents = 0.;
            graph.createComponents(minCorrectedJaccardForComponents, minComponentSize);
            const GlobalPathGraph1DisplayOptions options(minCorrectedJaccardForComponents, 1.);
            graph.writeComponentsGraphviz("PathGraphB", options);

            graph.connectSeedChains2(minEdgeCoverage, minCorrectedJaccard, connectors);
        }
        else {
            SHASTA_ASSERT(0);
        }
        graph.writeConnectors(connectors);



        // Use the ChainConnectors to stitch together the seed chains.
        graph.stitchSeedChains(connectors, minComponentSize);

        vector<Chain> chains;
        graph.createChainsFromComponents(minEstimatedLength, chains);
        cout << "Found " << chains.size() << " chains." << endl;

#if 0
        cout << timestamp << "Assembling the chains." << endl;
        ofstream fasta("Chains.fasta");
        graph.assembleChains(chains, fasta, "Chain-");
#endif
    }
}



void GlobalPathGraph1::assemble1(const Assembler& assembler)
{
    const uint64_t minPrimaryCoverage = 8;
    const uint64_t maxPrimaryCoverage = 45;
    const uint64_t minEdgeCoverage = 1;
    const double minCorrectedJaccard = 0.;
    const uint64_t minComponentSize = 3;
    const uint64_t transitiveReductionDistance = 20;
    const uint64_t compressedTransitiveReductionDistance = 100;
    const uint64_t minReliableLength = 200;
    const uint64_t crossEdgeCoverageThreshold1 = 3;
    const uint64_t crossEdgeCoverageThreshold2 = 6;
    const uint64_t detangleTolerance = 1;


    GlobalPathGraph1 graph(assembler);
    graph.createVertices(minPrimaryCoverage, maxPrimaryCoverage);
    graph.computeOrientedReadJourneys();
    graph.createEdges0(1, minEdgeCoverage, minCorrectedJaccard);

    graph.createComponents(minCorrectedJaccard, minComponentSize);

    // Assemble each connected component separately.
    for(uint64_t componentId=0; componentId<graph.components.size(); componentId++) {
        assemble1(graph, componentId,
            transitiveReductionDistance,
            compressedTransitiveReductionDistance,
            minReliableLength,
            crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2,
            detangleTolerance);
    }
}



void GlobalPathGraph1::assemble1(
    GlobalPathGraph1& globalGraph,
    uint64_t componentId,
    uint64_t transitiveReductionDistance,
    uint64_t compressedTransitiveReductionDistance,
    uint64_t minReliableLength,
    uint64_t crossEdgeCoverageThreshold1,
    uint64_t crossEdgeCoverageThreshold2,
    uint64_t detangleTolerance)
{
    cout << "Assembly begins for connected component " << componentId << endl;
    PathGraph1& component = *globalGraph.components[componentId];

    // Local transitive reduction.
    // This flags the non-transitive reduction edges.
    // No edges are removed.
    component.localTransitiveReduction(transitiveReductionDistance);

    // Graphviz output.
    GlobalPathGraph1DisplayOptions options;
    options.showNonTransitiveReductionEdges = false;
    component.writeGraphviz(globalGraph.verticesVector,
        "PathGraph" + to_string(componentId), options);
    options.makeCompact();
    component.writeGraphviz(globalGraph.verticesVector,
        "PathGraphCompact" + to_string(componentId), options);

    // Create a compressed representation of this connected component.
    // In this compressed representation, each linear sequence of
    // transitive reduction non-branch vertices becomes a single vertex.
    // Non-branch vertices are those with in-degree and out-degree not greater than 1.
    CompressedPathGraph1 cGraph(component, componentId, globalGraph.assembler);
    // cGraph.writeGraphviz(minReliableLength, "Initial");

    // Strict detangling, with zero detangle tolerance.
    cGraph.detangleIteration(
        compressedTransitiveReductionDistance,
        minReliableLength,
        0);

    // Remove cross-edges, then detangle again, this time with a looser detangle tolerance.
    cGraph.removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
    cGraph.detangleIteration(
        compressedTransitiveReductionDistance,
        minReliableLength,
        detangleTolerance);


    // Final output.
    cGraph.writeGraphviz(100000, "");
    cGraph.writeVerticesCsv();
    cout << "The CompressedPathGraph1 has " << num_vertices(cGraph) << " vertices and " <<
        num_edges(cGraph) << " edges." << endl;
    // cout << "Assembling sequence." << endl;
    // cGraph.assembleVertices();

}



void CompressedPathGraph1::writeVerticesCsv() const
{
    const CompressedPathGraph1& cGraph = *this;

    ofstream csv("CompressedPathGraphVertices" + to_string(componentId) + ".csv");
    csv << "Id,Begin,End,Vertex count,Estimated length\n";

    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const CompressedPathGraph1Vertex& cVertex = cGraph[cv];
        SHASTA_ASSERT(not cVertex.v.empty());
        const PathGraph1::vertex_descriptor v0 = cVertex.v.front();
        const PathGraph1::vertex_descriptor v1 = cVertex.v.back();

        uint64_t baseOffset = totalBaseOffset(cv);

        csv << componentId << "-" << cVertex.id << ",";
        csv << graph[v0].edgeId << ",";
        csv << graph[v1].edgeId << ",";
        csv << cVertex.v.size() << ",";
        csv << baseOffset << ",";
        csv << "\n";
    }
}



// This calls the lower level function twice, with and without labels.
void CompressedPathGraph1::writeGraphviz(
    uint64_t minGreenLength,
    const string& fileNamePrefix) const
{
    bool labels = true;
    writeGraphviz(labels, minGreenLength, fileNamePrefix);
    labels = false;
    writeGraphviz(labels, minGreenLength, fileNamePrefix);
}



void CompressedPathGraph1::writeGraphviz(
    bool labels,
    uint64_t minGreenLength,
    const string& fileNamePrefix) const
{
    const CompressedPathGraph1& cGraph = *this;

    const string name = "CompressedPathGraph" + to_string(componentId);
    ofstream out(name + fileNamePrefix + (labels ? ".dot" : "-NoLabels.dot"));
    out << "digraph " << name << "{\n";

    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        out << "\"" << componentId << "-" << cGraph[cv].id << "\" [";

        if(cGraph[cv].v.size() > 1) {

            const PathGraph1::vertex_descriptor v0 = cGraph[cv].v.front();
            const PathGraph1::vertex_descriptor v1 = cGraph[cv].v.back();
            const uint64_t baseOffset = totalBaseOffset(cv);

            if(baseOffset >= minGreenLength) {
                out <<
                    "style=filled color=LightGreen ";
            }

            if(labels) {
                out <<
                    "label=\""
                    << componentId << "-" << cGraph[cv].id << "\\n" <<
                    cGraph[cv].v.size() << " vertices\\n" <<
                    "First " << graph[v0].edgeId << "\\n" <<
                    "Last " << graph[v1].edgeId << "\\n" <<
                    "Length " << baseOffset <<
                    "\"";
            }

        } else {

            SHASTA_ASSERT(cGraph[cv].v.size() == 1);

            if(labels) {
                const PathGraph1::vertex_descriptor v = cGraph[cv].v.front();
                out <<
                    "label=\"" <<
                    componentId << "-" << cGraph[cv].id << "\\n" <<
                    graph[v].edgeId << "\\n" <<
                    "\"";
            }

        }

        out << "];\n";
    }

    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        const auto cv0 = source(ce, cGraph);
        const auto cv1 = target(ce, cGraph);
        out <<
            "\"" << componentId << "-" << cGraph[cv0].id << "\"->" <<
            "\"" << componentId << "-" << cGraph[cv1].id << "\"";
        if(labels) {
            const uint64_t offset = cGraph[ce].info.offsetInBases;
            out <<
                " ["
                " label=\"" <<
                offset << "\\n" << cGraph[ce].info.common <<
                "\""
                "]";
        }
        out << ";\n";
    }
    out << "}";
}



GlobalPathGraph1::GlobalPathGraph1(const Assembler& assembler) :
    assembler(assembler)
{
#if 0
    // The code below was moved to GlobalPathGraph1::assemble.

    // PARAMETERS TO BE EXPOSED WHEN CODE STABILIZES
    // ARE DEFINED BEFORE EACH PHASE THAT USES THEM.

    // Create GlobalPathGraph1 vertices.
    {
        const uint64_t minPrimaryCoverage = 8;
        const uint64_t maxPrimaryCoverage = 25;
        cout << timestamp << "Creating vertices." << endl;
        createVertices(minPrimaryCoverage, maxPrimaryCoverage);
        cout << timestamp << "Found " << verticesVector.size() << " vertices." << endl;
        computeOrientedReadJourneys();
    }

    // Create GlobalPathGraph1 edges using strict criteria suitable
    // for the generation of seed chains.
    {
        const uint64_t maxDistanceInJourney = 20;
        const uint64_t minEdgeCoverage = 3;
        const double minCorrectedJaccard = 0.8;
        cout << timestamp << "Creating edges." << endl;
        createEdges0(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard);
        cout << timestamp << "Found " << edges.size() << " edges." << endl;
    }

#if 0
    // Display.
    {
        const double minCorrectedJaccard = 0.8;
        const uint64_t minComponentSize = 100;
        createComponents(minCorrectedJaccard, minComponentSize);

        // Clean it up a bit for easier display.
        const uint64_t k = 3;
        knn(k);
        transitiveReduction();

        writeComponentsGraphviz("PathGraph", minCorrectedJaccard, 1.);
        return;
    }
#endif



    // To create seed chains, create connected components using only the best edges,
    // then compute the longest path in each connected component.
    {
        // Compute connected components using only the best edges.
        const double minCorrectedJaccard = 0.8;
        const uint64_t minComponentSize = 3;
        createComponents(minCorrectedJaccard, minComponentSize);
        writeComponentsGraphviz("PathGraphA", minCorrectedJaccard, 1.);

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
        cout << "Creating seed chains." << endl;
        createChainsFromComponents(minEstimatedLength, seedChains);
        cout << "Found " << seedChains.size() << " seed chains." << endl;
        writeSeedChains();
        writeSeedChainsDetails();
        writeSeedChainsStatistics();

#if 0
        cout << timestamp << "Assembling the seed chains." << endl;
        ofstream fasta("SeedChains.fasta");
        assembleChains(seedChains, fasta, "SeedChain-");
#endif
    }



#if 0

    // EXPERIMENT WITH connectSeedChains2.

    // Recreate GlobalPathGraph1 edges using criteria suitable
    // for the rest of the process.
    {
        const uint64_t maxDistanceInJourney = 1;
        const uint64_t minEdgeCoverage = 1;
        const double minCorrectedJaccard = 0;

        cout << timestamp << "Recreating edges." << endl;
        edges.clear();
        createEdges0(maxDistanceInJourney, minEdgeCoverage, minCorrectedJaccard);
        cout << timestamp << "Found " << edges.size() << " edges." << endl;
    }


    // Connect seed chains.
    vector<ChainConnector> connectors;
    {
        const double minCorrectedJaccardForComponents = 0.;
        const uint64_t minComponentSize = 3;
        const uint64_t minCommonForConnector = 3;
        const double minCorrectedJaccardForConnector = 0.7;

        // Compute connected components.
        createComponents(minCorrectedJaccardForComponents, minComponentSize);
        writeComponentsGraphviz("PathGraphB", minCorrectedJaccardForComponents, 1.);

        connectSeedChains2(
            minCommonForConnector,
            minCorrectedJaccardForConnector,
            connectors);
    }
#endif



    vector<ChainConnector> connectors;
    {
        const uint64_t minEdgeCoverage = 4;
        const double minCorrectedJaccard = 0.6;
        connectSeedChains1(minEdgeCoverage, minCorrectedJaccard, connectors);
    }


    writeConnectors(connectors);


    // Use the ChainConnectors to stitch together the seed chains,
    // then assemble the chains obtained in this way.
    {
        const uint64_t minComponentSize = 3;
        stitchSeedChains(connectors, minComponentSize);

        vector<Chain> chains;
        const uint64_t minEstimatedLength = 10000;
        cout << "Creating chains." << endl;
        createChainsFromComponents(minEstimatedLength, chains);
        cout << "Found " << chains.size() << " chains." << endl;

#if 1
        cout << timestamp << "Assembling the chains." << endl;
        ofstream fasta("Chains.fasta");
        assembleChains(chains, fasta, "Chain-");
#endif
    }
#endif
}



// Store in this PathGraph1 the seed chains found in another PathGraph1
void GlobalPathGraph1::storeSeedChains(const vector<Chain>& seedChainsArgument)
{
    seedChains = seedChainsArgument;

    for(uint64_t seedChainId=0; seedChainId<seedChains.size(); seedChainId++) {
        Chain& seedChain = seedChains[seedChainId];
        SHASTA_ASSERT(seedChain.vertexIds.size() == seedChain.edgeIds.size());

        for(uint64_t position=0; position<seedChain.vertexIds.size(); position++) {

            // Update the vertexId to reflect the vertex numbering
            // of this new graph.
            const MarkerGraphEdgeId edgeId = seedChain.edgeIds[position];
            const uint64_t vertexId = getVertexId(edgeId);
            SHASTA_ASSERT(vertexId != invalid<uint64_t>);
            seedChain.vertexIds[position] = vertexId;

            // Store chain information in the vertex.
            GlobalPathGraph1Vertex& vertex = verticesVector[vertexId];
            SHASTA_ASSERT(vertex.chainId == invalid<uint64_t>);
            SHASTA_ASSERT(vertex.positionInChain == invalid<uint64_t>);
            vertex.chainId = seedChainId;
            vertex.positionInChain = position;
            vertex.isFirstInChain = (position == 0);
            vertex.isLastInChain = (position == seedChain.vertexIds.size() - 1);
        }
    }
}



// Write each connected component in graphviz format.
void GlobalPathGraph1::writeComponentsGraphviz(
    const string& baseName,
    const GlobalPathGraph1DisplayOptions& options) const
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        const PathGraph1& component = *components[componentRank];
        component.writeGraphviz(
            verticesVector,
            baseName + "_Component_" + to_string(componentRank),
            options);
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
    // Requiring this is too strict and removes opportunities for
    // finding paths.
    if(not isBranchEdge(edgeId, minPrimaryCoverage)) {
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

    verticesVector.clear();

    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        if(isPrimary(edgeId, minPrimaryCoverage, maxPrimaryCoverage)) {
            verticesVector.push_back(GlobalPathGraph1Vertex(edgeId));
        }
    }
}



// Return the vertexId corresponding to a given MarkerGraphEdgeId, or
// invalid<MarkerGraphEdgeId> if no such a vertex exists.
uint64_t GlobalPathGraph1::getVertexId(MarkerGraphEdgeId edgeId) const
{
    GlobalPathGraph1Vertex targetVertex(edgeId);
    auto it = std::lower_bound(verticesVector.begin(), verticesVector.end(), targetVertex);

    if((it == verticesVector.end()) or (it->edgeId != edgeId)) {

        // Not found.
        return invalid<uint64_t>;

    } else {

        // Found it. Return its vertexId.
        return it - verticesVector.begin();

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

    for(uint64_t vertexId=0; vertexId<verticesVector.size(); vertexId++) {
        const MarkerGraphEdgeId edgeId = verticesVector[vertexId].edgeId;

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
                verticesVector[vertexId].journeyInfoItems.push_back({orientedReadId, position});
            }
        }
    }



    // Write the journeys to csv.
    ofstream csv("GlobalPathGraphJourneys.csv");
    for(ReadId readId=0; readId<assembler.markers.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = orientedReadJourneys[orientedReadId.getValue()];
            csv << orientedReadId << ",";
            for(const auto& p: journey) {
                const uint64_t vertexId = p.second;
                const MarkerGraphEdgeId edgeId = verticesVector[vertexId].edgeId;
                csv << edgeId << ",";
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
    // cout << timestamp << "Found " << candidateEdges.size() << " candidate edges, including duplicates." << endl;

    // Deduplicate the candidate edges and count the number of times
    // each of them was found. Keep only the ones that occurred at least
    // minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateEdges, coverage, minEdgeCoverage);
    // cout << timestamp << "After deduplication, there are " << candidateEdges.size() << " candidate edges." << endl;
    SHASTA_ASSERT(candidateEdges.size() == coverage.size());
    candidateEdges.shrink_to_fit();
    coverage.shrink_to_fit();

    // For each candidate edge, compute correctedJaccard, and if high enough
    // generate an edge.
    edges.clear();
    for(uint64_t i=0; i<candidateEdges.size(); i++) {
        const uint64_t c = coverage[i];
        SHASTA_ASSERT(c >= minEdgeCoverage);
        const auto& p = candidateEdges[i];
        const uint64_t vertexId0 = p.first;
        const uint64_t vertexId1 = p.second;
        GlobalPathGraph1Edge edge;
        const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1].edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, edge.info));
        if(edge.info.correctedJaccard() >= minCorrectedJaccard) {
            edge.vertexId0 = vertexId0;
            edge.vertexId1 = vertexId1;
            edge.coverage = c;
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
    GlobalPathGraph1Edge edge;
    for(edge.vertexId0=0; edge.vertexId0<verticesVector.size(); edge.vertexId0++) {
        if((edge.vertexId0 % 10000) == 0) {
            cout << edge.vertexId0 << "/" << verticesVector.size() << endl;
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
    const GlobalPathGraph1Vertex& vertex0 = verticesVector[vertexId0];
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
        const GlobalPathGraph1Vertex& vertex1 = verticesVector[vertexId1];
        const MarkerGraphEdgeId edgeId1 = vertex1.edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
        if(info.correctedJaccard() >= minCorrectedJaccard) {
            children.push_back({vertexId1, info});
        }
    }
}



void GlobalPathGraph1::findParents(
    uint64_t vertexId0,
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard,
    vector< pair<uint64_t, MarkerGraphEdgePairInfo> >& parents)
{
    const GlobalPathGraph1Vertex& vertex0 = verticesVector[vertexId0];
    const MarkerGraphEdgeId edgeId0 = vertex0.edgeId;

    // Find vertices encountered earlier on journeys that go through here.
    vector<uint64_t> candidateParents;
    for(const auto& journeyInfoItem: vertex0.journeyInfoItems) {
        const OrientedReadId orientedReadId = journeyInfoItem.orientedReadId;
        const auto& journey = orientedReadJourneys[orientedReadId.getValue()];
        for(int64_t position = int64_t(journeyInfoItem.positionInJourney)-1;
            position>=0; position--) {
            candidateParents.push_back(journey[position].second);
        }
    }

    // Count the candidate parents and only keep the ones
    // that occurred at least minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateParents, coverage, minEdgeCoverage);

    // Store the ones that have sufficient correctedJaccard.
    parents.clear();
    MarkerGraphEdgePairInfo info;
    for(const uint64_t vertexId1: candidateParents) {
        const GlobalPathGraph1Vertex& vertex1 = verticesVector[vertexId1];
        const MarkerGraphEdgeId edgeId1 = vertex1.edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId1, edgeId0, info));
        if(info.correctedJaccard() >= minCorrectedJaccard) {
            parents.push_back({vertexId1, info});
        }
    }
}



// Write the entire GlobalPathGraph in graphviz format.
void GlobalPathGraph1::writeGraphviz() const
{
    ofstream out("GlobalPathGraph.dot");
    out << "digraph GlobalPathGraph {\n";

    for(const GlobalPathGraph1Edge& edge: edges) {
        const MarkerGraphEdgeId edgeId0 = verticesVector[edge.vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[edge.vertexId1].edgeId;
        out << edgeId0 << "->";
        out << edgeId1 << ";\n";
    }

    out << "}\n";
}



// Create connected components.
// This only considers edges with corrected Jaccard at least equal to
// minCorrectedJaccard, and only stores connected components with at
// least minComponentSize vertices.
void GlobalPathGraph1::createComponents(
    double minCorrectedJaccard,
    uint64_t minComponentSize)
{
    // Compute connected components.
    const uint64_t n = verticesVector.size();
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    for(const GlobalPathGraph1Edge& edge: edges) {
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
        const GlobalPathGraph1Vertex& vertex = verticesVector[vertexId];
        const uint64_t componentId = disjointSets.find_set(vertexId);
        allComponents[componentId]->addVertex(vertexId, vertex.edgeId);
    }

    // Create edges of each connected component.
    for(const GlobalPathGraph1Edge& edge: edges) {
        if(edge.info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        const uint64_t vertexId0 = edge.vertexId0;
        const uint64_t vertexId1 = edge.vertexId1;
        const uint64_t componentId = disjointSets.find_set(vertexId0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(vertexId1));
        const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1].edgeId;
        allComponents[componentId]->addEdge(edgeId0, edgeId1, edge.info, edge.coverage);
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



// Local transitive reduction of each connected component.
void GlobalPathGraph1::localTransitiveReduction(uint64_t distance)
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        component.localTransitiveReduction(distance);
    }
}



void PathGraph1::localTransitiveReduction(uint64_t distance)
{
    PathGraph1& graph = *this;

    // We want to process edges in order of increasing coverage.
    vector<pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        edgeTable.push_back({e, graph[e].coverage});
    }
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<edge_descriptor, uint64_t>());

    // Loop over all edges v0->v1 in order of increasing coverage.
    for(const auto& p: edgeTable) {
        edge_descriptor e01 = p.first;
        const vertex_descriptor v0 = source(e01, graph);
        const vertex_descriptor v1 = target(e01, graph);

        // Do a BFS starting at v0, up to a distance maxPathLength.
        // Stop if we encounter v1.

        // The BFS queue.
        std::queue<vertex_descriptor> q;
        q.push(v0);

        // The vertices we encountered so far, with their distance from v0.
        std::map<vertex_descriptor, uint64_t> m;
        m.insert({v0, 0});

        // BFS loop.
        // cout << "BFS loop begins for " << v0 << "->" << v1 << endl;
        while(not q.empty()) {

            // Dequeue a vertex.
            const vertex_descriptor vA = q.front();
            q.pop();
            const auto itA = m.find(vA);
            SHASTA_ASSERT(itA != m.end());
            const uint64_t distanceA = itA->second;
            const uint64_t distanceB = distanceA + 1;
            // cout << "Dequeued " << vA << " at distance " << distanceA << endl;

            // Loop over the out-edges of vA.
            bool endBfs = false;
            BGL_FORALL_OUTEDGES_T(vA, eAB, graph, PathGraph1) {

                // Dont's use e01 in the BFS.
                if(eAB == e01) {
                    continue;
                }

                // If we reached v1, mark e01 as a nonTransitiveReduction edge
                // and stop the BFS.
                const vertex_descriptor vB = target(eAB, graph);
                if(vB == v1) {
                    boost::remove_edge(e01, graph);
                    endBfs = true;
                    // cout << "Reached " << v1 << endl;
                    break;
                }

                // If we already reached this vertex, do nothing.
                if(m.contains(vB)) {
                    continue;
                }

                // If not at maximum distance, enqueue vB.
                if(distanceB < distance) {
                    q.push(vB);
                    m.insert({vB, distanceB});
                    // cout << "Enqueued " << vB << " at distance " << distanceB << endl;
                }
            }
            if(endBfs) {
                break;
            }
        }
    }

}



bool CompressedPathGraph1::localTransitiveReduction(uint64_t distance)
{
    CompressedPathGraph1& cGraph = *this;

    // We want to process edges in order of increasing coverage.
    vector<pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        edgeTable.push_back({ce, cGraph[ce].info.common});
    }
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<edge_descriptor, uint64_t>());

    // Loop over all edges v0->v1 in order of increasing coverage.
    uint64_t removedEdgeCount = 0;
    for(const auto& p: edgeTable) {
        edge_descriptor e01 = p.first;
        const vertex_descriptor v0 = source(e01, cGraph);
        const vertex_descriptor v1 = target(e01, cGraph);

        // Do a BFS starting at v0, up to a distance maxPathLength.
        // Stop if we encounter v1.

        // The BFS queue.
        std::queue<vertex_descriptor> q;
        q.push(v0);

        // The vertices we encountered so far, with their distance from v0.
        std::map<vertex_descriptor, uint64_t> m;
        m.insert({v0, 0});

        // BFS loop.
        // cout << "BFS loop begins for " << v0 << "->" << v1 << endl;
        while(not q.empty()) {

            // Dequeue a vertex.
            const vertex_descriptor vA = q.front();
            q.pop();
            const auto itA = m.find(vA);
            SHASTA_ASSERT(itA != m.end());
            const uint64_t distanceA = itA->second;
            const uint64_t distanceB = distanceA + 1;
            // cout << "Dequeued " << vA << " at distance " << distanceA << endl;

            // Loop over the out-edges of vA.
            bool endBfs = false;
            BGL_FORALL_OUTEDGES_T(vA, eAB, cGraph, CompressedPathGraph1) {

                // Dont's use e01 in the BFS.
                if(eAB == e01) {
                    continue;
                }

                // If we reached v1, mark e01 as a nonTransitiveReduction edge
                // and stop the BFS.
                const vertex_descriptor vB = target(eAB, cGraph);
                if(vB == v1) {
                    ++removedEdgeCount;
                    boost::remove_edge(e01, cGraph);
                    endBfs = true;
                    // cout << "Reached " << v1 << endl;
                    break;
                }

                // If we already reached this vertex, do nothing.
                if(m.contains(vB)) {
                    continue;
                }

                // If not at maximum distance, enqueue vB.
                if(distanceB < distance) {
                    q.push(vB);
                    m.insert({vB, distanceB});
                    // cout << "Enqueued " << vB << " at distance " << distanceB << endl;
                }
            }
            if(endBfs) {
                break;
            }
        }
    }

    return removedEdgeCount > 0;
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
    const MarkerGraphEdgePairInfo& info,
    uint64_t coverage)
{
    auto it0 = vertexMap.find(edgeId0);
    auto it1 = vertexMap.find(edgeId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    add_edge(v0, v1, {info, coverage}, *this);
}



// Write a PathGraph1 in graphviz format.
void PathGraph1::writeGraphviz(
    const vector<GlobalPathGraph1Vertex>& globalVertices,
    const string& name,
    const GlobalPathGraph1DisplayOptions& options) const
{
    ofstream out(name + ".dot");

    const PathGraph1& graph = *this;
    out << "digraph " << name << " {\n";

    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        const PathGraph1Vertex& vertex = graph[v];
        const GlobalPathGraph1Vertex& globalVertex = globalVertices[vertex.vertexId];
        out << vertex.edgeId;

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "[";
        }

        if(options.labels) {
            out << "label=\"";
            out << vertex.edgeId << "\\n" << globalVertex.journeyInfoItems.size();
            if(globalVertex.chainId != invalid<uint64_t>) {
                out << "\\n" << globalVertex.chainId << ":" << globalVertex.positionInChain;
                if(globalVertex.isFirstInChain) {
                    out << "\\nFirst in chain";
                }
                if(globalVertex.isLastInChain) {
                    out << "\\nLast in chain";
                }
            }
            out << "\" ";
        }

        if(options.tooltips) {
            out << "tooltip=\"";
            out << vertex.edgeId;
            if(globalVertex.chainId != invalid<uint64_t>) {
                out << " " << globalVertex.chainId << ":" << globalVertex.positionInChain;
                if(globalVertex.isFirstInChain) {
                    out << " first in chain";
                }
                if(globalVertex.isLastInChain) {
                    out << " last in chain";
                }
            }
            out << "\" ";
        }

        // If it belongs to a chain, color it.
        if(options.colorVertices) {
            if(globalVertex.chainId != invalid<uint64_t>) {
                if(globalVertex.isFirstInChain) {
                    out << " style=filled fillcolor=blue ";
                } else if(globalVertex.isLastInChain) {
                    out << " style=filled fillcolor=orange ";
                } else {
                    // const uint32_t hue = MurmurHash2(&vertex.chainId, sizeof(vertex.chainId), 231) % 100;
                    // out << " color=\"" << 0.01 * double(hue) << ",0.4,1\"";
                    out << " style=filled fillcolor=cyan ";
                }
            }
        }

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "]";
        }
        out << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1Edge& edge = graph[e];
        if(not options.showNonTransitiveReductionEdges and edge.isNonTransitiveReductionEdge) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        out <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId;

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << " [";
        }

        if(edge.isNonTransitiveReductionEdge) {
            out << "style=dashed ";
        }

        if(options.tooltips) {
            out <<
                "tooltip=\"" <<
                graph[v0].edgeId << "->" <<
                graph[v1].edgeId << " ";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << " " <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << " " <<
                edge.info.offsetInBases << "\" ";
        }

        if(options.labels) {
            out <<
                "label=\"";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << "\\n" <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << "\\n" <<
                edge.info.offsetInBases << "\" ";

        }

        // Color.
        if(options.colorEdges) {
            const double correctedJaccard = edge.info.correctedJaccard();
            if(correctedJaccard <= options.redJ) {
                out << " color=red ";
            } else if(correctedJaccard >= options.greenJ) {
                out << " color=green ";
            } else {
                const double hue = (correctedJaccard - options.redJ) / (3. * (options.greenJ - options.redJ));
                out << " color=\"" << hue << ",1,1\" ";
            }
        }

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << "]";
        }
        out << ";\n";
    }

    out << "}\n";
}



// For each vertex, only keep the best k outgoing and k incoming edges.
// "Best" as defined by correctedJaccard of the edges.
void PathGraph1::knn(uint64_t k)
{
    PathGraph1& graph = *this;

    // Store here the edges we want to keep.
    std::set<edge_descriptor> edgesToBeKept;

    // For each vertex, mark as to be kept the best k outgoing
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
                edgesToBeKept.insert(e);
            }
        }
    }

    // Remove edges not marked as to be kept.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        if(not edgesToBeKept.contains(e)) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// For each connected component, use the longest path
// to create a Chain. Only keep the ones that are sufficiently long.
// This also stores chain information in the vertices.
void GlobalPathGraph1::createChainsFromComponents(
    uint64_t minEstimatedLength,
    vector<Chain>& chains)
{
    chains.clear();

    // Clear chain information from the vertices.
    for(GlobalPathGraph1Vertex& vertex: verticesVector) {
        vertex.chainId = invalid<uint64_t>;
        vertex.positionInChain = invalid<uint64_t>;
        vertex.isFirstInChain = false;
        vertex.isLastInChain = false;
    }

    // Loop over connected components.
    vector<PathGraph1::vertex_descriptor> chainVertices;
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];

        // Compute the longest path in this component.
        longestPath(component, chainVertices);

        // Use this longest path to create a tentative Chain.
        Chain chain;
        for(const PathGraph1::vertex_descriptor v: chainVertices) {
            const uint64_t vertexId = component[v].vertexId;
            chain.vertexIds.push_back(vertexId);
            chain.edgeIds.push_back(verticesVector[vertexId].edgeId);
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


        // Store this chain as a new chain.
        const uint64_t chainId = chains.size();
        cout << "Chain " << chainId << " has estimated length " << totalBaseOffset << endl;
        chains.push_back(chain);
        for(uint64_t position=0; position<chain.vertexIds.size(); position++) {
            const uint64_t vertexId = chain.vertexIds[position];
            GlobalPathGraph1Vertex& vertex = verticesVector[vertexId];
            vertex.chainId = chainId;
            vertex.positionInChain = position;
            vertex.isFirstInChain = (position == 0);
            vertex.isLastInChain = (position == chain.vertexIds.size() - 1);
        }
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



uint64_t GlobalPathGraph1::ChainConnector::totalOffset() const
{
    uint64_t totalBaseOffset = 0;
    for(const MarkerGraphEdgePairInfo& info: infos) {
        totalBaseOffset += info.offsetInBases;
    }
    return totalBaseOffset;
}



void GlobalPathGraph1::assembleChains(
    const vector<Chain>& chains,
    ostream& fasta,
    const string& csvPrefix) const
{
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const Chain& chain = chains[chainId];
        cout << timestamp << chainId << "/" << chains.size() << endl;

        // Generate an AssemblyPath using this chain.
        vector<MarkerGraphEdgeId> markerGraphEdgeIds;
        for(const uint64_t vertexId: chain.vertexIds) {
            markerGraphEdgeIds.push_back(verticesVector[vertexId].edgeId);
        }
        AssemblyPath assemblyPath(assembler, markerGraphEdgeIds, chain.infos);
        assemblyPath.writeFasta(fasta, to_string(chainId));

        ofstream csv(csvPrefix + to_string(chainId) + ".csv");
        assemblyPath.writeCsv(csv);
    }
}



void GlobalPathGraph1::writeSeedChains() const
{
    ofstream csv("SeedChains.csv");
    csv << "ChainId,First,Last,Length,Estimated length\n";

    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        const Chain& chain = seedChains[chainId];

        csv << chainId << ",";
        csv << verticesVector[chain.vertexIds.front()].edgeId << ",";
        csv << verticesVector[chain.vertexIds.back()].edgeId << ",";
        csv << chain.vertexIds.size() << ",";
        csv << chain.totalOffset() << ",";
        csv << "\n";

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
                verticesVector[vertexId].edgeId << ",";
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
            const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0].edgeId;
            for(uint64_t i1=i0+1; i1<vertexIds.size(); i1++) {
                const uint64_t vertexId1 = vertexIds[i1];
                const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1].edgeId;
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
            SHASTA_ASSERT(vertexId < verticesVector.size());
            const GlobalPathGraph1Vertex& vertex = verticesVector[vertexId];
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

        const MarkerGraphEdgeId edgeId0 = verticesVector[chain0.vertexIds.back()].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[chain1.vertexIds.front()].edgeId;
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

    ofstream out("SeedChains.dot");
    out << "digraph SeedChains {\n";
    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        out << chainId << " [label=\"SC" << chainId << "\"];\n";
    }

    connectors.clear();

    // Loop over chains.
    for(uint64_t chainId=0; chainId<seedChains.size(); chainId++) {
        for(uint64_t direction=0; direction<2; direction++) {
            connectSeedChain1(chainId, direction,
                minEdgeCoverage, minCorrectedJaccard, connectors, out);
        }
    }
    out << "}\n";
}



// Names here refer to direction=0 (forward).
void GlobalPathGraph1::connectSeedChain1(
    uint64_t chainId,
    uint64_t direction, // 0 = forward, 1 = backward
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard,
    vector<ChainConnector>& connectors,
    ostream& out
    )
{
    const bool debug = false; // (chainId == 16) and (direction == 0);

    if(debug) {
        cout << "Working on chain " << chainId << " direction " << direction << endl;
    }
    const Chain& chain = seedChains[chainId];
    const uint64_t chainLastVertexId = (direction == 0) ? chain.vertexIds.back() : chain.vertexIds.front();
    if(debug) {
        cout << "Last vertex of chain " <<  verticesVector[chainLastVertexId].edgeId <<
            " is starting vertex for search." << endl;
    }
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
        if(debug) {
            cout << "Visiting " << verticesVector[v0.vertexId].edgeId << endl;
        }

        // If it belongs to a different chain, we are done.
        // Create a new connector.
        const uint64_t chainId0 = verticesVector[v0.vertexId].chainId;
        if(chainId0 != invalid<uint64_t> and chainId0 != chainId) {
            if(debug) {
                cout << verticesVector[v0.vertexId].edgeId << " is on chain " << chainId0 <<
                    " at distance " << v0.distance << endl;
            }
            if(direction == 0) {
                out << chainId << "->" << chainId0 << ";\n";
            } else {
                out << chainId0 << "->" << chainId << ";\n";
            }

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
            if(direction == 0) {
                connector.reverse();
            }
            break;
        }

        // Find its children.
        if(direction == 0) {
            findChildren(v0.vertexId, minEdgeCoverage, minCorrectedJaccard, children);
        } else {
            findParents (v0.vertexId, minEdgeCoverage, minCorrectedJaccard, children);
        }


        // Loop over the unvisited children.
        for(const auto& child: children) {
            const uint64_t vertexId1 = child.first;
            if(debug) {
                cout << "Found child " << verticesVector[vertexId1].edgeId << endl;
            }
            if(visitedById.find(vertexId1) != visitedById.end()) {
                if(debug) {
                    cout << "Already visited." << endl;
                }
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
            const MarkerGraphEdgeId edgeId = verticesVector[vertexId].edgeId;
            graph.addVertex(vertexId, edgeId);
            // cout << "Adding vertex for " << edgeId << endl;
        }

        // Add the edges of this chain
        for(uint64_t position=0; position<chain.vertexIds.size()-1; position++) {
            const uint64_t vertexId0 = chain.vertexIds[position];
            const uint64_t vertexId1 = chain.vertexIds[position + 1];
            const GlobalPathGraph1Vertex& vertex0 = verticesVector[vertexId0];
            const GlobalPathGraph1Vertex& vertex1 = verticesVector[vertexId1];
            graph.addEdge(vertex0.edgeId, vertex1.edgeId, chain.infos[position], invalid<uint64_t>);
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
            const MarkerGraphEdgeId edgeId = verticesVector[vertexId].edgeId;
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
            const GlobalPathGraph1Vertex& vertex0 = verticesVector[vertexId0];
            const GlobalPathGraph1Vertex& vertex1 = verticesVector[vertexId1];
            // cout << "Adding edge for " << vertex0.edgeId << " " << vertex1.edgeId << endl;
            graph.addEdge(vertex0.edgeId, vertex1.edgeId, connector.infos[position], invalid<uint64_t>);
        }

    }

    GlobalPathGraph1DisplayOptions options(0.5, 1.);
    options.labels = false;
    graph.writeGraphviz(verticesVector, "StitchedSeedChains", options);

    cout << "The stitched graph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;


    // Replace our existing connected components with
    // the connected components of the stitched graph.
    components = graph.createConnectedComponents(minComponentSize);
    cout << "The stitched graph has " << components.size() <<
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
            graph[e].info,
            graph[e].coverage);
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



void GlobalPathGraph1::connectSeedChains2(
    uint64_t minCommonCount,
    double minCorrectedJaccard,
    vector<ChainConnector>& connectors)
{

    // Check that each seed chain appears entirely in a single connected component.
    {
        vector<uint64_t> seedChainTable(seedChains.size(), invalid<uint64_t>);
        for(uint64_t componentId=0; componentId<components.size(); componentId++) {
            const shared_ptr<const PathGraph1>& componentPointer = components[componentId];
            const PathGraph1& component = *componentPointer;
            BGL_FORALL_VERTICES(v, component, PathGraph1) {
                const uint64_t vertexId = component[v].vertexId;
                const uint64_t chainId = verticesVector[vertexId].chainId;
                if(chainId == invalid<uint64_t>) {
                    continue;
                }
                if(seedChainTable[chainId] == invalid<uint64_t>) {
                    seedChainTable[chainId] = componentId;
                } else {
                    SHASTA_ASSERT(seedChainTable[chainId] == componentId);
                }
            }
        }
    }

    ofstream out("SeedChains.dot");
    out << "digraph SeedChains {\n";
    for(uint64_t seedChainId=0; seedChainId!=seedChains.size(); seedChainId++) {
        const Chain& seedChain = seedChains[seedChainId];
        out << seedChainId << " [label=\"" << seedChainId << "\\n" << seedChain.totalOffset() << "\"];\n";
    }

    // Process each component independently.
    connectors.clear();
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const shared_ptr<const PathGraph1>& componentPointer = components[componentId];
        const PathGraph1& component = *componentPointer;
        connectSeedChains2(componentId, component, minCommonCount, minCorrectedJaccard, out, connectors);
    }

    out << "}\n";
}



void GlobalPathGraph1::connectSeedChains2(
    uint64_t componentId,
    const PathGraph1& component,
    uint64_t minCommonCount,
    double minCorrectedJaccard,
    ostream& out,
    vector<ChainConnector>& connectors)
{

    // Find the seed chains that appear in this connected component.
    std::set<uint64_t> seedChainIds;
    BGL_FORALL_VERTICES(v, component, PathGraph1) {
        const uint64_t vertexId = component[v].vertexId;
        const uint64_t chainId = verticesVector[vertexId].chainId;
        if(chainId != invalid<uint64_t>) {
            seedChainIds.insert(chainId);
        }
    }

    // Loop over the seed chains that appear in this connected component.
    for(const uint64_t seedChainId0: seedChainIds) {
        const Chain& seedChain = seedChains[seedChainId0];
        SHASTA_ASSERT(not seedChain.vertexIds.empty());

        // Extend this chain forward.
        cout << "Extending forward seed chain " << seedChainId0 <<
            " in connected component " << componentId << endl;
        Chain extendedChain = seedChain;
        const uint64_t seedChainId1 =
            extendChainForward(seedChainId0, extendedChain, component, minCommonCount, minCorrectedJaccard);

        if(seedChainId1 == invalid<uint64_t>) {
            cout << "Extension did not reach another chain." << endl;
            continue;
        }

        cout << "Reached seed chain " << seedChainId1 << endl;

        // Create a connector between these two chains.
        // Use the portion of the extendedChain that is past the end of the seedChain.
        connectors.resize(connectors.size() + 1);
        ChainConnector& connector = connectors.back();
        connector.chainId0 = seedChainId0;
        connector.chainId1 = seedChainId1;
        connector.vertexIds.push_back(seedChain.vertexIds.back());
        for(uint64_t i = seedChain.infos.size(); i < extendedChain.infos.size(); i++) {
            connector.vertexIds.push_back(extendedChain.vertexIds[i+1]);
            connector.infos.push_back(extendedChain.infos[i]);
        }
        out << seedChainId0 << "->" << seedChainId1 <<
            " [label=\"" << connector.totalOffset() << "\"];\n";
    }
}


// Extend a chain forward until we bump into another chain.
// This returns the chainId of the chain we found, or invalid<uint64_t>
// if none found.
uint64_t GlobalPathGraph1::extendChainForward(
    uint64_t chainId,
    Chain& chain,
    const PathGraph1& component,
    uint64_t minCommonCount,
    double minCorrectedJaccard) const
{
    using boost::multi_index_container;
    using boost::multi_index::indexed_by;
    using boost::multi_index::member;
    using boost::multi_index::ordered_unique;
    using boost::multi_index::ordered_non_unique;

    // Information about a vertex stored during the shortest path algorithm.
    class VertexInfo {
    public:
        PathGraph1::vertex_descriptor v;

        // The sum of offsets starting at the start vertex.
        uint64_t distance;
    };

    // The container used to store VertexInfos,
    // with the two indices to support the required operations.
    using VertexContainer = multi_index_container<VertexInfo,
        indexed_by <
        ordered_unique< member<VertexInfo, PathGraph1::vertex_descriptor, &VertexInfo::v> >,
        ordered_non_unique< member<VertexInfo, uint64_t, &VertexInfo::distance> >
        > >;

    const bool debug = false;

    // Outer iteration loop.
    // At each iteration we add one vertex to the Chain.
    while(true) {

        const uint64_t vertexIdA = chain.vertexIds.back();
        const MarkerGraphEdgeId edgeIdA = verticesVector[vertexIdA].edgeId;

        const auto itA = component.vertexMap.find(edgeIdA);
        SHASTA_ASSERT(itA != component.vertexMap.end());
        const PathGraph1::vertex_descriptor vA = itA->second;
        SHASTA_ASSERT(component[vA].vertexId == vertexIdA);

        if(debug) {
            cout << "Extending forward a chain with " << chain.vertexIds.size() <<
                " vertices, current last vertex is " << edgeIdA << endl;
        }

        // Do a shortest path search, starting at vA and moving forward.
        // Stop the search if we find a vertex with sufficient
        // commonCount and correctedJaccard relative to vA.

        // Start the Dijkstra algorithm with this as the only seen vertex.
        VertexContainer seen;
        const auto& seenByVertexDescriptor = seen.get<0>();
        auto& seenByDistance = seen.get<1>();   // Not const so we can erase by distance.
        seen.insert({vA, 0});

        VertexContainer visited;
        const auto& visitedByVertexDescriptor = visited.get<0>();

        // Dijkstra shortest path loop.
        bool wasExtended = false;
        while(not seen.empty()) {

            // Visit the vertex with the lowest distance.
            auto it0 = seenByDistance.begin();
            VertexInfo vertexInfo0 = *it0;
            const PathGraph1::vertex_descriptor v0 = vertexInfo0.v;
            const uint64_t distance0 = vertexInfo0.distance;
            const uint64_t edgeId0 = verticesVector[component[v0].vertexId].edgeId;
            visited.insert(*it0);
            seenByDistance.erase(it0);

            // Check it against vA.
            // LATER WE SHOULD ALSO CHECK AGAINST THE LAST FEW VERTICES IN THE CHAIN.
            if(v0 != vA) {
                MarkerGraphEdgePairInfo infoA0;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeIdA, edgeId0, infoA0));
                if(debug) {
                    cout << "Visiting " << edgeId0 <<
                        " at distance " << distance0 <<
                        " common " << infoA0.common;
                    if(infoA0.common > 0) {
                        cout << ", J' " << infoA0.correctedJaccard() <<
                            ", base offset " << infoA0.offsetInBases;
                    }
                    cout << endl;
                }

                // If this is good enough, use it to extend the chain and stop the shortest path search.
                if(
                    infoA0.offsetInBases>=0 and
                    infoA0.common >= minCommonCount and
                    infoA0.correctedJaccard() >= minCorrectedJaccard) {
                    wasExtended = true;
                    chain.vertexIds.push_back(component[v0].vertexId);
                    chain.infos.push_back(infoA0);
                    if(debug) {
                        cout << "Chain extended with " << edgeId0 <<
                            " at distance " << distance0 <<
                            " common " << infoA0.common << " , J' " << infoA0.correctedJaccard() <<
                            " , base offset " << infoA0.offsetInBases << endl;
                    }

                    // If we reached another chain, stop here.
                    const uint64_t chainId0 = verticesVector[component[v0].vertexId].chainId;
                    if(chainId0 != invalid<uint64_t> and chainId0 != chainId) {
                        if(debug) {
                            cout << "Reached chain " << chainId0 << endl;
                        }
                        return chainId0;
                    }

                    break;
                }
            }

            // Loop over the unvisited children.
            BGL_FORALL_OUTEDGES(v0, e, component, PathGraph1) {
                const PathGraph1::vertex_descriptor v1 = target(e, component);
                if(debug) {
                    cout << "Found child " << verticesVector[component[v1].vertexId].edgeId << endl;
                }
                if(visitedByVertexDescriptor.find(v1) != visitedByVertexDescriptor.end()) {
                    if(debug) {
                        cout << "Already visited." << endl;
                    }
                    continue;
                }

                const MarkerGraphEdgePairInfo& info = component[e].info;
                SHASTA_ASSERT(info.offsetInBases > 0);

                // Update the tentative distance of the child,
                // adding it to the seenVertices if not present.
                auto it1 = seenByVertexDescriptor.find(v1);
                const uint64_t distance1 = distance0 + info.offsetInBases;
                if(it1 == seenByVertexDescriptor.end()) {
                    seen.insert({v1, distance1});
                    if(debug) {
                        cout << "Seen: " << verticesVector[component[v1].vertexId].edgeId <<
                            " at distance " << distance1 << endl;
                    }
                } else {
                    VertexInfo seenVertex1 = *it1;
                    if(distance1 < seenVertex1.distance) {
                        seenVertex1.distance = distance1;
                        seen.replace(it1, seenVertex1);
                        if(debug) {
                            cout << "Seen: updated " << verticesVector[component[v1].vertexId].edgeId <<
                                ", new distance " << distance1 << endl;
                        }
                    }
                }
            }
        }

        if(not wasExtended) {
            if(debug) {
                cout << "Shortest path search terminated." << endl;
            }
            break;
        }
    }

    // If getting here, chain extension ended before we reached another chain.
    if(debug) {
        cout << "Forward chain extension did not reach another chain." << endl;
    }
    return invalid<uint64_t>;
}



void GlobalPathGraph1::writeConnectors(const vector<ChainConnector>& connectors) const
{
    ofstream csv("ChainConnectors.csv");
    csv << "ChainId0,ChainId1,Position,EdgeId0,EdgeId1,Common,Offset,J,J'\n";
    for(const ChainConnector& connector: connectors) {
        for(uint64_t position=0; position<connector.infos.size(); position++) {
            const MarkerGraphEdgePairInfo& info = connector.infos[position];
            const uint64_t vertexId0 = connector.vertexIds[position];
            const uint64_t vertexId1 = connector.vertexIds[position + 1];
            const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0].edgeId;
            const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1].edgeId;
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



CompressedPathGraph1::CompressedPathGraph1(
    const PathGraph1& graph,
    uint64_t componentId,
    const Assembler& assembler) :
    graph(graph),
    componentId(componentId),
    assembler(assembler)
{
    CompressedPathGraph1& cGraph = *this;
    using boost::out_degree;
    using boost::in_degree;

    // Create a filtered version of the PathGraph1, containing only the
    // transitive reduction edges.
    class EdgePredicate {
    public:
        bool operator()(const PathGraph1::edge_descriptor e) const
        {
            return not (*graph)[e].isNonTransitiveReductionEdge;
        }
        EdgePredicate(const PathGraph1& graph) : graph(&graph) {}
        EdgePredicate() : graph(0) {}
    private:
        const PathGraph1* graph;
    };
    using FilteredPathGraph1 = boost::filtered_graph<PathGraph1, EdgePredicate>;
    FilteredPathGraph1 filteredGraph(graph, EdgePredicate(graph));

    vector<uint64_t> lengthHistogram(2, 0);

    // Loop over vertices of the PathGraph1 to generate the CompressedPathGraph1 vertices.
    std::set<PathGraph1::vertex_descriptor> visited;
    vector<PathGraph1::vertex_descriptor> forwardVertices;
    vector<PathGraph1::vertex_descriptor> backwardVertices;
    BGL_FORALL_VERTICES(v, filteredGraph, FilteredPathGraph1) {

        // If already visited, skip.
        if(visited.contains(v)) {
            continue;
        }
        visited.insert(v);

        // cout << "Starting at " << filteredGraph[v].edgeId << endl;

        // If in-degree or out-degree is more than 1, generate a single CompressedPathGraph1 vertex.
        if(in_degree(v, filteredGraph) > 1 or out_degree(v, filteredGraph) > 1) {
            const CompressedPathGraph1::vertex_descriptor cv = add_vertex(cGraph);
            CompressedPathGraph1Vertex& cVertex = cGraph[cv];
            cVertex.id = nextVertexId++;
            cVertex.v.push_back(v);
            // cout << "Done: is branch vertex." << endl;
            ++lengthHistogram[1];
            continue;
        }

        // Move forward until we find a branch vertex.
        bool isCircular = false;
        forwardVertices.clear();
        PathGraph1::vertex_descriptor u = v;
        while(true) {
            if(out_degree(u, filteredGraph) == 0) {
                break;
            }
            SHASTA_ASSERT(out_degree(u, filteredGraph) == 1);

            // Move forward.
            FilteredPathGraph1::out_edge_iterator it;
            tie(it, ignore) = out_edges(u, filteredGraph);
            u = target(*it, filteredGraph);

            if(in_degree(u, filteredGraph) > 1 or out_degree(u, filteredGraph) > 1) {
                break;
            }

            if(u == v) {
                isCircular = true;
            }

            forwardVertices.push_back(u);
            SHASTA_ASSERT(not visited.contains(u));
            visited.insert(u);
            // cout << "Moving forward found " << filteredGraph[u].edgeId << endl;
        }

        // For now, don't handle the circular case.
        SHASTA_ASSERT(not isCircular);

        // Move backward until we find a branch vertex.
        backwardVertices.clear();
        u = v;
        while(true) {
            if(in_degree(u, filteredGraph) == 0) {
                break;
            }
            SHASTA_ASSERT(in_degree(u, filteredGraph) == 1);

            // Move backward.
            FilteredPathGraph1::in_edge_iterator it;
            tie(it, ignore) = in_edges(u, filteredGraph);
            u = source(*it, filteredGraph);

            if(in_degree(u, filteredGraph) > 1 or out_degree(u, filteredGraph) > 1) {
                break;
            }

            SHASTA_ASSERT(u != v);

            backwardVertices.push_back(u);
            SHASTA_ASSERT(not visited.contains(u));
            visited.insert(u);
            // cout << "Moving backward found " << filteredGraph[u].edgeId << endl;
        }

        // Create a new CompressedPathGraph1 vertex.
        const CompressedPathGraph1::vertex_descriptor cv = add_vertex(cGraph);
        CompressedPathGraph1Vertex& cVertex = cGraph[cv];
        cVertex.id = nextVertexId++;
        copy(backwardVertices.rbegin(), backwardVertices.rend(), back_inserter(cVertex.v));
        cVertex.v.push_back(v);
        copy(forwardVertices.begin(),forwardVertices.end(), back_inserter(cVertex.v));

    }
    SHASTA_ASSERT(visited.size() == num_vertices(graph));



    // To generate edges, create a map giving the compressed vertex descriptor
    // that begins with a given vertex.
    std::map<PathGraph1::vertex_descriptor, CompressedPathGraph1::vertex_descriptor> vertexMap;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        vertexMap.insert({cGraph[cv].v.front(), cv});
    }
    BGL_FORALL_VERTICES(cv0, cGraph, CompressedPathGraph1) {
        const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v.back();
        BGL_FORALL_OUTEDGES(v0, e, filteredGraph, FilteredPathGraph1) {
            const PathGraph1::vertex_descriptor v1 = target(e, filteredGraph);
            auto it1 = vertexMap.find(v1);
            SHASTA_ASSERT(it1 != vertexMap.end());
            const CompressedPathGraph1::vertex_descriptor cv1 = it1->second;
            add_edge(cv0, cv1, {graph[e].info}, cGraph);
        }
    }

    uint64_t filteredGraphEdgeCount = 0;
    BGL_FORALL_EDGES(e, filteredGraph, FilteredPathGraph1) {
        ++filteredGraphEdgeCount;
    }
    cout << "The PathGraph1 has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;
    if(false) {
        cout << "The FilteredPathGraph1 has " << num_vertices(graph) << " vertices and " <<
            filteredGraphEdgeCount << " edges." << endl;
        cout << "The CompressedPathGraph1 has " << num_vertices(cGraph) << " vertices and " <<
            num_edges(cGraph) << " edges." << endl;
    }

}



bool CompressedPathGraph1::mergeLinearChains()
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;

    vector< vector<vertex_descriptor> > linearChains;
    findLinearVertexChains(cGraph, linearChains);

    bool changesWereMade = false;
    for(const auto& linearChain: linearChains) {
        if(linearChain.size() == 1) {
            continue;
        }

        // The first and last vertex can have in/out degrees greater than one
        // and will not be merged.
        vector<vertex_descriptor> verticesToBeMerged;
        for(uint64_t i=0; i<linearChain.size(); i++) {
            const vertex_descriptor cv = linearChain[i];
            if(i==0 or i==linearChain.size()-1) {
                if(in_degree(cv, cGraph)>1 or out_degree(cv, cGraph)>1) {
                    continue;
                }
            }
            verticesToBeMerged.push_back(cv);
        }
        if(verticesToBeMerged.size() < 2) {
            continue;
        }
        changesWereMade = true;


        // Create the new vertex.
        const vertex_descriptor cvNew = add_vertex(cGraph);
        auto& cVertexNew = cGraph[cvNew];
        cVertexNew.id = nextVertexId++;
        for(const vertex_descriptor cv: verticesToBeMerged) {
            const auto& cVertex = cGraph[cv];
            copy(cVertex.v.begin(), cVertex.v.end(), back_inserter(cVertexNew.v));
        }

        if(debug) {
            cout << componentId << "-" << cVertexNew.id << " created by merging";
            for(const auto cv: verticesToBeMerged) {
                cout << " " << componentId << "-" << cGraph[cv].id;
            }
            cout << endl;
        }

        // Reroute in-edges.
        BGL_FORALL_INEDGES(verticesToBeMerged.front(), ce, cGraph, CompressedPathGraph1) {
            const auto& oldEdge = cGraph[ce];
            add_edge(source(ce, cGraph), cvNew, oldEdge, cGraph);
        }

        // Reroute out-edges.
        BGL_FORALL_OUTEDGES(verticesToBeMerged.back(), ce, cGraph, CompressedPathGraph1) {
            const auto& oldEdge = cGraph[ce];
            add_edge(cvNew, target(ce, cGraph), oldEdge, cGraph);
        }

        // Now we can remove the vertices we merged.
        for(const vertex_descriptor cv: verticesToBeMerged) {
            clear_vertex(cv, cGraph);
            remove_vertex(cv, cGraph);
        }

    }

    if(debug) {
        cout << "After merging linear chains, the CompressedPathGraph1 has " <<
            num_vertices(cGraph) << " vertices  and " <<
            num_edges(cGraph) << " edges." << endl;
    }

    return changesWereMade;
}



// An edge cv0->cv1 is a cross edge if:
// - Has coverage <= threshold1
// - cv0 has out-degree > 1 and at least one out-edge with coverage >= threshold2.
// - cv1 has in-degree  > 1 and at least one in-edge  with coverage >= threshold2.
bool CompressedPathGraph1::removeCrossEdges(
    uint64_t threshold1,
    uint64_t threshold2)
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;
    SHASTA_ASSERT(threshold2 > threshold1);



    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {

        // Check coverage.
        if(cGraph[ce].info.common > threshold1) {
            continue;
        }

        // Check degrees.
        const vertex_descriptor cv0 = source(ce, cGraph);
        const vertex_descriptor cv1 = target(ce, cGraph);
        if(out_degree(cv0, cGraph) < 2) {
            continue;
        }
        if(in_degree(cv1, cGraph) < 2) {
            continue;
        }

        // Check that v0 has at least one out-edge with coverage >= threshold2.
        bool found0 = false;
        BGL_FORALL_OUTEDGES(cv0, ce02, cGraph, CompressedPathGraph1) {
            if(cGraph[ce02].info.common >= threshold2) {
                found0 = true;
                break;
            }
        }
        if(not found0) {
            continue;
        }

        // Check that v1 has at least one in-edge with coverage >= threshold2.
        bool found1 = false;
        BGL_FORALL_INEDGES(cv1, ce21, cGraph, CompressedPathGraph1) {
            if(cGraph[ce21].info.common >= threshold2) {
                found1 = true;
                break;
            }
        }
        if(not found1) {
            continue;
        }

        edgesToBeRemoved.push_back(ce);

        if(debug) {
            cout << "Cross-edge " <<
                componentId << "-" << cGraph[cv0].id << "->" <<
                componentId << "-" << cGraph[cv1].id <<
                " will be removed." << endl;
        }

    }



    for(const edge_descriptor ce: edgesToBeRemoved) {
        boost::remove_edge(ce, cGraph);
    }

    return not edgesToBeRemoved.empty();
}



string CompressedPathGraph1::vertexIdString(vertex_descriptor cv) const
{
    const CompressedPathGraph1& cGraph = *this;

    return to_string(componentId) + "-" + to_string(cGraph[cv].id);
}



bool CompressedPathGraph1::detangleVertices(uint64_t detangleTolerance)
{
    CompressedPathGraph1& cGraph = *this;

    // Find the vertices that can potentially be detangled.
    vector<vertex_descriptor> detangleCandidates;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        if(in_degree(cv, cGraph) > 1 and out_degree(cv, cGraph) > 1) {
            detangleCandidates.push_back(cv);
        }
    }

    // Do the detangling.
    uint64_t detangledCount = 0;
    for(const auto cv: detangleCandidates) {
        bool wasDetangled = detangleVertex(cv, detangleTolerance);
        if(wasDetangled) {
            ++detangledCount;
        }
    }

    return detangledCount > 0;
}



// Attempt to detangle a compressed vertex.
// The detangle operation leaves all in-degree and out-degrees
// of all other vertices unchanged.
bool CompressedPathGraph1::detangleVertex(
    vertex_descriptor cv,
    uint64_t detangleTolerance)
{
    CompressedPathGraph1& cGraph = *this;

    const bool debug = false;
    if(debug) {
        cout << "Attempting to detangle " << componentId << "-" << cGraph[cv].id << endl;
    }

    // Gather the source vertices of incoming edges.
    vector<vertex_descriptor> incoming;
    BGL_FORALL_INEDGES(cv, e, cGraph, CompressedPathGraph1) {
        incoming.push_back(source(e, cGraph));
    }

    // Gather the target vertices of outgoing edges.
    vector<vertex_descriptor> outgoing;
    BGL_FORALL_OUTEDGES(cv, e, cGraph, CompressedPathGraph1) {
        outgoing.push_back(target(e, cGraph));
    }

    if(debug) {
        cout << "Incoming:";
        for(const auto cv1: incoming) {
            const auto& cVertex = cGraph[cv1];
            const PathGraph1::vertex_descriptor v1 = cVertex.v.back();
            cout << " " << graph[v1].edgeId;
        }
        cout << endl;

        cout << "Outgoing:";
        for(const auto cv1: outgoing) {
            const auto& cVertex = cGraph[cv1];
            const PathGraph1::vertex_descriptor v1 = cVertex.v.front();
            cout << " " << graph[v1].edgeId;
        }
        cout << endl;
    }

    // We can only detangle if the number of incoming edges
    // equals the number of outgoing edges.
    if(incoming.size() != outgoing.size()) {
        return false;
    }

    // For each incoming edge, count the number of
    // outgoing edges it has more than detangleTolerance common oriented reads with.
    vector<uint64_t> inCount(incoming.size(), 0);

    // For each outgoing edge, count the number of
    // incoming edges it has  more than detangleTolerance common oriented reads with.
    vector<uint64_t> outCount(outgoing.size(), 0);

    // Loop over pairs of incoming/outgoing edges.
    MarkerGraphEdgePairInfo info;
    class NewEdge {
    public:
        vertex_descriptor cv0;
        vertex_descriptor cv1;
        MarkerGraphEdgePairInfo info;
    };
    vector<NewEdge> newEdges;
    for(uint64_t i0=0; i0<incoming.size(); i0++) {
        const auto cv0 = incoming[i0];
        const auto& cVertex0 = cGraph[cv0];
        const PathGraph1::vertex_descriptor v0 = cVertex0.v.back();
        const MarkerGraphEdgeId edgeId0 = graph[v0].edgeId;
        for(uint64_t i1=0; i1<outgoing.size(); i1++) {
            const auto cv1 = outgoing[i1];
            const auto& cVertex1 = cGraph[cv1];
            const PathGraph1::vertex_descriptor v1 = cVertex1.v.front();
            const MarkerGraphEdgeId edgeId1 = graph[v1].edgeId;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            if(debug) {
                cout << edgeId0 << " " << edgeId1 << ": " << info.common << endl;
            }
            if(info.common > detangleTolerance) {
                ++inCount[i0];
                ++outCount[i1];
                newEdges.push_back({cv0, cv1, info});
            }
        }
    }

    if(debug) {
        cout << "inCount ";
        copy(inCount.begin(), inCount.end(), ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
        cout << "outCount ";
        copy(outCount.begin(), outCount.end(), ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
    }

    // We can only detangle if all the inCount and outCount entries are 1.
    // This means that for each incoming edge there is only one outgoing edge
    // with common oriented reads, and vice versa.
    for(const uint64_t c: inCount) {
        if(c != 1) {
            if(debug) {
                cout << "Cannot detangle." << endl;
            }
            return false;
        }
    }
    for(const uint64_t c: outCount) {
        if(c != 1) {
            if(debug) {
                cout << "Cannot detangle." << endl;
            }
            return false;
        }
    }

    // Add the edges.
    for(const auto& newEdge: newEdges) {
        edge_descriptor ce;

        // Check that the edge does not already exists.
        bool edgeExists = false;
        tie(ce, edgeExists) = boost::edge(newEdge.cv0, newEdge.cv1, cGraph);
        if(edgeExists) {
            continue;
        }

        bool edgeWasAdded = false;
        tie(ce, edgeWasAdded) = boost::add_edge(newEdge.cv0, newEdge.cv1, {newEdge.info}, cGraph);
        SHASTA_ASSERT(edgeWasAdded);
        if(debug) {
            const auto& cEdge = cGraph[ce];
            cout << "Added compressed edge " <<
                componentId << "-" << cGraph[newEdge.cv0].id << "->" <<
                componentId << "-" << cGraph[newEdge.cv1].id << ", common count " <<
                cEdge.info.common << ", offset " <<
                cEdge.info.offsetInBases << endl;
        }
    }

    // Now we can remove the vertex we detangled.
    clear_vertex(cv, cGraph);
    remove_vertex(cv, cGraph);

    return true;
}



bool CompressedPathGraph1::detangleLinearChains(uint64_t detangleTolerance)
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;

    // Find the linear chains.
    vector< vector<CompressedPathGraph1::edge_descriptor> > linearChains;
    findLinearChains(cGraph, 1, linearChains);
    if(debug) {
        cout << "Found " << linearChains.size() << " linear chains." << endl;
    }



    // Try to detangle all linear chains.
    uint64_t detangleCount = 0;
    for(const auto& linearChain: linearChains) {
        SHASTA_ASSERT(not linearChain.empty());

        // If the first vertex has in-degree less than 2, do nothing.
        const auto ce0 = linearChain.front();
        const auto cv0 = source(ce0, cGraph);
        const uint64_t inDegree = in_degree(cv0, cGraph);
        if(inDegree < 2) {
            continue;
        }

        // If the last vertex has out-degree less than 2, do nothing.
        const auto ce1 = linearChain.back();
        const auto cv1 = target(ce1, cGraph);
        const uint64_t outDegree = out_degree(cv1, cGraph);
        if(outDegree < 2) {
            continue;
        }

        // We can only detangle if the in-degree and out-degree are the same.
        if(inDegree != outDegree) {
            continue;
        }

        // Gather the vertices of this chain.
        vector<vertex_descriptor> chainVertices;
        chainVertices.push_back(cv0);
        for(const auto ce: linearChain) {
            chainVertices.push_back(target(ce, cGraph));
        }

        // Sanity check: all vertices except the first and last must have
        // in-degree and out-degree 1.
        for(uint64_t i=1; i<chainVertices.size()-1; i++) {
            const auto cv = chainVertices[i];
            SHASTA_ASSERT(in_degree(cv, cGraph) == 1);
            SHASTA_ASSERT(out_degree(cv, cGraph) == 1);
        }

        // If the first vertex has out-degree>1, we cannot detangle.
        if(out_degree(cv0, cGraph) > 1) {
            continue;
        }
        // If the last vertex has in-degree>1, we cannot detangle.
        if(in_degree(cv1, cGraph) > 1) {
            continue;
        }

        if(debug) {
            cout << "Attempting to detangle linear chain:";
            for(const auto cv: chainVertices) {
                cout << " " << componentId << "-" << cGraph[cv].id;
            }
            cout << endl;
        }

        // Gather the source vertices of incoming edges.
        vector<vertex_descriptor> incoming;
        BGL_FORALL_INEDGES(cv0, e, cGraph, CompressedPathGraph1) {
            incoming.push_back(source(e, cGraph));
        }

        // Gather the target vertices of outgoing edges.
        vector<vertex_descriptor> outgoing;
        BGL_FORALL_OUTEDGES(cv1, e, cGraph, CompressedPathGraph1) {
            outgoing.push_back(target(e, cGraph));
        }

        if(debug) {
            cout << "Incoming:";
            for(const auto cv1: incoming) {
                const auto& cVertex = cGraph[cv1];
                const PathGraph1::vertex_descriptor v1 = cVertex.v.back();
                cout << " " << graph[v1].edgeId;
            }
            cout << endl;

            cout << "Outgoing:";
            for(const auto cv1: outgoing) {
                const auto& cVertex = cGraph[cv1];
                const PathGraph1::vertex_descriptor v1 = cVertex.v.front();
                cout << " " << graph[v1].edgeId;
            }
            cout << endl;
        }
        // For each incoming edge, count the number of
        // outgoing edges it has more than detangleTolerance common oriented reads with.
        vector<uint64_t> inCount(incoming.size(), 0);

        // For each outgoing edge, count the number of
        // incoming edges it has more than detangleTolerance common oriented reads with.
        vector<uint64_t> outCount(outgoing.size(), 0);

        // Loop over pairs of incoming/outgoing edges.
        MarkerGraphEdgePairInfo info;
        class NewEdge {
        public:
            vertex_descriptor cv0;
            vertex_descriptor cv1;
            MarkerGraphEdgePairInfo info;
        };

        vector<NewEdge> newEdges;
        for(uint64_t i0=0; i0<incoming.size(); i0++) {
            const auto cv0 = incoming[i0];
            const auto& cVertex0 = cGraph[cv0];
            const PathGraph1::vertex_descriptor v0 = cVertex0.v.back();
            const MarkerGraphEdgeId edgeId0 = graph[v0].edgeId;
            for(uint64_t i1=0; i1<outgoing.size(); i1++) {
                const auto cv1 = outgoing[i1];
                const auto& cVertex1 = cGraph[cv1];
                const PathGraph1::vertex_descriptor v1 = cVertex1.v.front();
                const MarkerGraphEdgeId edgeId1 = graph[v1].edgeId;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
                if(debug) {
                    cout << edgeId0 << " " << edgeId1 << ": " << info.common << endl;
                }
                if(info.common > detangleTolerance) {
                    ++inCount[i0];
                    ++outCount[i1];
                    newEdges.push_back({cv0, cv1, info});
                }
            }
        }

        if(debug) {
            cout << "inCount ";
            copy(inCount.begin(), inCount.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
            cout << "outCount ";
            copy(outCount.begin(), outCount.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
        }

        // We can only detangle if all the inCount and outCount entries are 1.
        // This means that for each incoming edge there is only one outgoing edge
        // with common oriented reads, and vice versa.
        bool canDetangle = true;
        for(const uint64_t c: inCount) {
            if(c != 1) {
                if(debug) {
                    cout << "Cannot detangle." << endl;
                }
                canDetangle = false;
                break;
            }
        }
        for(const uint64_t c: outCount) {
            if(c != 1) {
                if(debug) {
                    cout << "Cannot detangle." << endl;
                }
                canDetangle = false;
                break;
            }
        }
        if(not canDetangle) {
            continue;
        }

        // Add the edges.
        for(const auto& newEdge: newEdges) {

            // Check that the edge does not already exists.
            edge_descriptor ce;
            bool edgeExists = false;
            tie(ce, edgeExists) = boost::edge(newEdge.cv0, newEdge.cv1, cGraph);
            if(edgeExists) {
                continue;
            }

            bool edgeWasAdded = false;
            tie(ce, edgeWasAdded) = boost::add_edge(newEdge.cv0, newEdge.cv1, {newEdge.info}, cGraph);
            SHASTA_ASSERT(edgeWasAdded);
            if(debug) {
                const auto& cEdge = cGraph[ce];
                cout << "Added compressed edge " <<
                    componentId << "-" << cGraph[newEdge.cv0].id << "->" <<
                    componentId << "-" << cGraph[newEdge.cv1].id << ", common count " <<
                    cEdge.info.common << ", offset " <<
                    cEdge.info.offsetInBases << endl;
            }
        }

        // Now we can remove the vertices of the linear chain we detangled.
        for(const auto cv: chainVertices) {
            clear_vertex(cv, cGraph);
            remove_vertex(cv, cGraph);
        }
        ++detangleCount;
    }

    return detangleCount > 0;
}




bool CompressedPathGraph1::detangleSuperbubbles(uint64_t minReliableLength)
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;

    // The edges to be added.
    class NewEdge {
    public:
        vertex_descriptor cv0;
        vertex_descriptor cv1;
        CompressedPathGraph1Edge edge;
    };
    vector<NewEdge> newEdges;

    // Loop over long vertices of the CompressedPathGraph1.
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const uint64_t baseOffset = totalBaseOffset(cv);
        if(baseOffset < minReliableLength) {
            continue;
        }
        if(debug) {
            cout << "Starting BFS at " << componentId << "-" << cGraph[cv].id << endl;
        }

        // Do a BFS starting here, stopping when we encounter a long vertices.
        std::queue<vertex_descriptor> q;
        q.push(cv);
        std::set<vertex_descriptor> visited;
        visited.insert(cv);
        while(not q.empty()) {
            const auto cv0 = q.front();
            q.pop();
            if(false) {
                cout << "Dequeued " << componentId << "-" << cGraph[cv0].id << endl;
            }

            BGL_FORALL_OUTEDGES(cv0, ce, cGraph, CompressedPathGraph1) {
                const auto cv1 = target(ce, cGraph);
                if(visited.contains(cv1)) {
                    continue;
                }
                const uint64_t baseOffset = totalBaseOffset(cv1);
                if(baseOffset >= minReliableLength) {
                    CompressedPathGraph1Edge newEdge;
                    assembler.analyzeMarkerGraphEdgePair(
                        graph[cGraph[cv].v.back()].edgeId,
                        graph[cGraph[cv1].v.front()].edgeId,
                        newEdge.info);
                    if(newEdge.info.common) {
                        newEdges.push_back({cv, cv1, newEdge});
                    }
                    if(debug) {
                        cout << "Found " << componentId << "-" << cGraph[cv1].id <<
                            ", common count " << newEdge.info.common << endl;
                    }
                }   else {
                    q.push(cv1);
                    visited.insert(cv1);
                    if(false) {
                        cout << "Enqueued " << componentId << "-" << cGraph[cv1].id << endl;
                    }
                }
            }
        }
    }



    // Now we can:
    // - Remove all short vertices.
    // - Remove all edges.
    // - Replace them with the new edges we found.

    // Remove all short vertices.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const uint64_t baseOffset = totalBaseOffset(cv);
        if(baseOffset < minReliableLength) {
            verticesToBeRemoved.push_back(cv);
        }
    }
    for(const auto cv: verticesToBeRemoved) {
        clear_vertex(cv, cGraph);
        remove_vertex(cv, cGraph);
    }

    // Remove all edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        edgesToBeRemoved.push_back(ce);
    }
    for(const auto ce: edgesToBeRemoved) {
        boost::remove_edge(ce, cGraph);
    }

    // Replace them with the new edges we found.
    for(const auto& newEdge: newEdges) {
        add_edge(newEdge.cv0, newEdge.cv1, {newEdge.edge}, cGraph);
    }
    if(debug) {
        cout << "After superbubble detangling, the CompressedPathGraph1 has " <<
            num_vertices(cGraph) << " vertices  and " <<
            num_edges(cGraph) << " edges." << endl;
    }

    return not verticesToBeRemoved.empty(); // Questionable.
}



uint64_t CompressedPathGraph1::totalBaseOffset(vertex_descriptor cv) const
{
    const CompressedPathGraph1& cGraph = *this;
    const CompressedPathGraph1Vertex& cVertex = cGraph[cv];

    SHASTA_ASSERT(not cVertex.v.empty());

    uint64_t totalBaseOffset = 0;

    for(uint64_t i=1; i<cVertex.v.size(); i++) {
        const PathGraph1::vertex_descriptor vA = cVertex.v[i-1];
        const PathGraph1::vertex_descriptor vB = cVertex.v[i];

        PathGraph1::edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(vA, vB, graph);
        if(edgeWasFound) {
            totalBaseOffset += graph[e].info.offsetInBases;
        } else {
            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(
                graph[vA].edgeId, graph[vB].edgeId, info));
            totalBaseOffset += info.offsetInBases;
        }
    }
    return totalBaseOffset;

}



void CompressedPathGraph1::detangleIteration(
    uint64_t compressedTransitiveReductionDistance,
    uint64_t minReliableLength,
    uint64_t detangleTolerance)
{

    for(uint64_t iteration=0; ; ++iteration) {

        // Try everything.
        const bool transitiveReductionChanges = localTransitiveReduction(compressedTransitiveReductionDistance);
        // writeGraphviz(minReliableLength, "A" + to_string(iteration));
        const bool detangleVerticesChanges = detangleVertices(detangleTolerance);
        // writeGraphviz(minReliableLength, "B" + to_string(iteration));
        const bool detangleLinearChainsChanges = detangleLinearChains(detangleTolerance);
        // writeGraphviz(minReliableLength, "C" + to_string(iteration));
        const bool mergeLinearChainsChanges = mergeLinearChains();
        // writeGraphviz(minReliableLength, "D" + to_string(iteration));
        const bool detangleSuperBubblesChanges = detangleSuperbubbles(minReliableLength);
        // writeGraphviz(minReliableLength, "E" + to_string(iteration));

        // If nothing changed, stop the iteration.
        if(not (
            transitiveReductionChanges or
            detangleVerticesChanges or
            detangleLinearChainsChanges or
            mergeLinearChainsChanges or
            detangleSuperBubblesChanges
            )) {
            break;
        }
    }

}



void CompressedPathGraph1::assembleVertices() const
{
    const CompressedPathGraph1& cGraph = *this;

    ofstream fasta("Component" + to_string(componentId) + ".fasta");
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        assembleVertex(cv, fasta);
    }
}



void CompressedPathGraph1::assembleVertex(
    vertex_descriptor cv,
    ostream& fasta) const
{
    const CompressedPathGraph1& cGraph = *this;
    const CompressedPathGraph1Vertex& cVertex = cGraph[cv];

    // Construct the MarkerGraphEdgeIds of the assembly path.
    vector<MarkerGraphEdgeId> edgeIds;
    for(const PathGraph1::vertex_descriptor v: cVertex.v) {
        edgeIds.push_back(graph[v].edgeId);
    }

    // Construct the MarkerGraphEdgePairInfo.
    vector<MarkerGraphEdgePairInfo> infos(edgeIds.size() - 1);
    for(uint64_t i=1; i<edgeIds.size(); i++) {
        const MarkerGraphEdgeId edgeId0 = edgeIds[i-1];
        const MarkerGraphEdgeId edgeId1 = edgeIds[i];
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, infos[i-1]));
    }

    // Create the AssemblyPath.
    AssemblyPath assemblyPath(assembler, edgeIds, infos);
    assemblyPath.writeFasta(fasta, vertexIdString(cv));
}

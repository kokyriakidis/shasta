// Shasta.
#include "mode3a-PackedAssemblyGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "removeReciprocalEdges.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include "dominatorTree.hpp"

// Standard library.
#include "fstream.hpp"
#include <map>
#include "utility.hpp"



PackedAssemblyGraph::PackedAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    uint64_t minLinkCoverage1,
    uint64_t minLinkCoverage2,
    uint64_t minLinkCoverage3,
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2,
    uint64_t minMarkerCount) :
    assemblyGraph(assemblyGraph)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    createVertices(minLinkCoverage1, minMarkerCount);
    computeJourneys();
    createEdges(minLinkCoverage2);
    removeReciprocalEdges(packedAssemblyGraph);
    writeGraphviz();
    writeJourneys();
    computePartialPaths(
        minLinkCoverage3,
        segmentCoverageThreshold1,
        segmentCoverageThreshold2);
    writePartialPaths();

    cout << "The PackedAssemblyGraph has " << num_vertices(packedAssemblyGraph) <<
        " vertices and " << num_edges(packedAssemblyGraph) << " edges." << endl;
}



void PackedAssemblyGraph::createVertices(
    uint64_t minLinkCoverage1,
    uint64_t minMarkerCount)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;

    // Gather all AssemblyGraph vertices with in-degree and out-degree
    // equal to 1 when considering only AssemblyGraph edges with coverage
    // at list equal to minLinkCoverage.
    // These are the vertices that can be referenced in PackedAssemblyGraph
    // vertices.
    // For each of them, store the previous and next vertex.
    std::map<AssemblyGraph::vertex_descriptor,
        pair<AssemblyGraph::vertex_descriptor, AssemblyGraph::vertex_descriptor> > vertexMap;
    vector<AssemblyGraph::edge_descriptor> inEdges;
    vector<AssemblyGraph::edge_descriptor> outEdges;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {

        // Gather in-edges.
        inEdges.clear();
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(assemblyGraph.edgeCoverage(e) >= minLinkCoverage1) {
                inEdges.push_back(e);
            }
        }
        if(inEdges.size() != 1) {
            // When considering only edges with sufficient coverage, v does not have in-degree 1.
            continue;
        }

        // Gather out-edges.
        outEdges.clear();
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(assemblyGraph.edgeCoverage(e) >= minLinkCoverage1) {
                outEdges.push_back(e);
            }
        }
        if(outEdges.size() != 1) {
            // When considering only edges with sufficient coverage, v does not have in-degree 1.
            continue;
        }

        // Add this vertex to our vertexMap.
        const AssemblyGraph::vertex_descriptor v0 = source(inEdges.front(), assemblyGraph);
        const AssemblyGraph::vertex_descriptor v1 = target(outEdges.front(), assemblyGraph);
        vertexMap.insert(make_pair(v, make_pair(v0, v1)));
    }


#if 0
    // Write out the vertexMap
    {
        ofstream csv("VertexMap.csv");
        for(const auto& p: vertexMap) {
            const auto v = p.first;
            const auto v0 = p.second.first;
            const auto v1 = p.second.second;
            csv << assemblyGraph.vertexStringId(v) << ",";
            csv << assemblyGraph.vertexStringId(v0) << ",";
            csv << assemblyGraph.vertexStringId(v1) << "\n";
        }
    }
#endif

    // Each sufficiently long linear chain generates a vertex of the PackedAssemblyGraph,
    // We are not interested in circular chains, so
    // each linear sequence begins at a vertex whose parent does not appear in the vertex map.
    uint64_t nextVertexId = 0;
    for(const auto& p: vertexMap) {
        const AssemblyGraph::vertex_descriptor v = p.first;
        const AssemblyGraph::vertex_descriptor v0 = p.second.first;
        if(vertexMap.contains(v0)) {
            // A linear chain does not begin at v.
            continue;
        }

        /// If getting here, a linear chain begins at v.
        vertex_descriptor pv = add_vertex(packedAssemblyGraph);
        PackedAssemblyGraphVertex& vertex = packedAssemblyGraph[pv];
        vertex.assemblyGraphVertices.push_back(v);
        AssemblyGraph::vertex_descriptor u = p.second.second;
        while(true) {
            const auto it = vertexMap.find(u);
            if(it == vertexMap.end()) {
                // u is not in our vertexMap, so it cannot be part of a linear chain.
                break;
            }
            vertex.assemblyGraphVertices.push_back(u);
            u = it->second.second;
        }
#if 0
        cout << "Chain:";
        for(const AssemblyGraph::vertex_descriptor u: vertex.path) {
            cout << " " << assemblyGraph.vertexStringId(u);
        }
        cout << "\n";
#endif

        // If the chain is too short, remove the vertex.
        // Otherwise, keep it and assign it an id.
        uint64_t totalMarkerCount = 0;
        for(const AssemblyGraph::vertex_descriptor v: vertex.assemblyGraphVertices) {
            const AssemblyGraphVertex& aVertex = assemblyGraph[v];
            const uint64_t segmentId = aVertex.segmentId;
            const uint64_t markerCount = assemblyGraph.packedMarkerGraph.segments[segmentId].size();
            totalMarkerCount += markerCount;
        }
        if(totalMarkerCount < minMarkerCount) {
            remove_vertex(pv, packedAssemblyGraph);
        } else {
            vertex.id = nextVertexId++;
        }
    }
}



void PackedAssemblyGraph::computeJourneys()
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;

    journeys.resize(assemblyGraph.journeys.size());
    for(uint64_t i=0; i<journeys.size(); i++) {
        journeys[i].resize(assemblyGraph.journeys[i].size(), null_vertex());
    }

    BGL_FORALL_VERTICES(pv, packedAssemblyGraph, PackedAssemblyGraph) {
        const auto& assemblyGraphVertices = packedAssemblyGraph[pv].assemblyGraphVertices;
        for(const auto av: assemblyGraphVertices) {
            const auto& assemblyGraphVertex = assemblyGraph[av];
            for(const JourneyEntry& journeyEntry: assemblyGraphVertex.journeyEntries) {
                const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
                const uint64_t position = journeyEntry.position;
                journeys[orientedReadId.getValue()][position] = pv;
            }
        }
    }

    // Remove from the journeys:
    // - Entries equal to null_vertex.
    // - Duplicate entries.
    for(auto& journey: journeys) {
        vector<vertex_descriptor> newJourney;
        for(const vertex_descriptor v: journey) {
            if(v == null_vertex()) {
                continue;
            }
            if(newJourney.empty() or v != newJourney.back()) {
                newJourney.push_back(v);
            }
        }
        journey.swap(newJourney);
    }


    // Store journey entries in the vertices.
    for(ReadId i=0; i<journeys.size(); i++) {
        const auto& journey = journeys[i];
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(i);
        for(uint64_t position=0; position<journey.size(); position++) {
            const vertex_descriptor v = journey[position];
            packedAssemblyGraph[v].journeyEntries.push_back({orientedReadId, position});
        }
    }

}


void PackedAssemblyGraph::createEdges(uint64_t minLinkCoverage2)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;

    vector< pair<vertex_descriptor, vertex_descriptor> > transitions;
    for(const auto& journey: journeys) {
        for(uint64_t i=1; i<journey.size(); i++) {
            transitions.push_back({journey[i-1], journey[i]});
        }
    }

    vector<uint64_t> frequency;
    deduplicateAndCount(transitions, frequency);

    for(uint64_t i=0; i<transitions.size(); i++) {
        if(frequency[i] >= minLinkCoverage2) {
            const auto& transition = transitions[i];
            edge_descriptor e;
            tie(e, ignore) = add_edge(transition.first, transition.second, {frequency[i]}, packedAssemblyGraph);
        }
    }
}



string PackedAssemblyGraph::vertexStringId(vertex_descriptor pv) const
{
    const PackedAssemblyGraph& packedAssemblyGraph = *this;
    const PackedAssemblyGraphVertex& packedAssemblyGraphVertex = packedAssemblyGraph[pv];
    const AssemblyGraph::vertex_descriptor v0 = packedAssemblyGraphVertex.assemblyGraphVertices.front();
    const AssemblyGraph::vertex_descriptor v1 = packedAssemblyGraphVertex.assemblyGraphVertices.back();

    return assemblyGraph.vertexStringId(v0) + "_" + assemblyGraph.vertexStringId(v1);
}



void PackedAssemblyGraph::writeGraphviz() const
{
    const PackedAssemblyGraph& packedAssemblyGraph = *this;

    ofstream dot("PackedAssemblyGraph.dot");
    dot << "digraph PackedAssemblyGraph {\n";

    BGL_FORALL_VERTICES(pv, packedAssemblyGraph, PackedAssemblyGraph) {
        const PackedAssemblyGraphVertex& pVertex = packedAssemblyGraph[pv];
        const AssemblyGraph::vertex_descriptor av0 = pVertex.assemblyGraphVertices.front();
        const AssemblyGraph::vertex_descriptor av1 = pVertex.assemblyGraphVertices.back();
        dot << pVertex.id;
        dot << " [label=\"P" <<
            pVertex.id << "\\n" <<
            assemblyGraph.vertexStringId(av0) << "\\n" <<
            assemblyGraph.vertexStringId(av1) << "\\n" <<
            pVertex.journeyEntries.size() <<
            "\"]";
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, packedAssemblyGraph, PackedAssemblyGraph) {
        const vertex_descriptor pv0 = source(e, packedAssemblyGraph);
        const vertex_descriptor pv1 = target(e, packedAssemblyGraph);
        dot <<
            packedAssemblyGraph[pv0].id << "->" <<
            packedAssemblyGraph[pv1].id <<
            " [label=" << packedAssemblyGraph[e].coverage << "]"
            ";\n";
    }

    dot << "}\n";
}



void PackedAssemblyGraph::writeJourneys() const
{
    const PackedAssemblyGraph& packedAssemblyGraph = *this;

    ofstream csv("PackedAssemblyGraphJourneys.csv");

    for(uint64_t i=0; i<journeys.size(); i++) {
        const vector<vertex_descriptor>& journey = journeys[i];
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        csv << orientedReadId << ",";
        for(const vertex_descriptor v: journey) {
            csv << "P" << packedAssemblyGraph[v].id << ",";
        }
        csv << "\n";
    }

}



void PackedAssemblyGraph::writePartialPaths() const
{
    const PackedAssemblyGraph& packedAssemblyGraph = *this;

    ofstream csv("PackedAssemblyGraphPartialPaths.csv");

    BGL_FORALL_VERTICES(v, packedAssemblyGraph, PackedAssemblyGraph) {
        const PackedAssemblyGraphVertex& vertex = packedAssemblyGraph[v];

        csv << "P" << vertex.id << ",Forward,";
        for(const vertex_descriptor u: vertex.forwardPartialPath) {
            csv << "P" << packedAssemblyGraph[u].id << ",";
        }
        csv << "\n";

        csv << "P" << vertex.id << ",Backward,";
        for(const vertex_descriptor u: vertex.backwardPartialPath) {
            csv << "P" << packedAssemblyGraph[u].id << ",";
        }
        csv << "\n";
    }
}



void PackedAssemblyGraph::computePartialPaths(
    uint64_t minLinkCoverage,
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;

    ofstream debugOut;
    // debugOut.open("PackedAssemblyGraph-computePartialPath.txt");
    BGL_FORALL_VERTICES(v, packedAssemblyGraph, PackedAssemblyGraph) {
        computePartialPath(
            v,
            minLinkCoverage,
            segmentCoverageThreshold1,
            segmentCoverageThreshold2,
            debugOut);
    }
}



void PackedAssemblyGraph::computePartialPath(
    vertex_descriptor vStart,
    uint64_t minLinkCoverage,
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2,
    ostream& debugOut)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    PackedAssemblyGraphVertex& startVertex = packedAssemblyGraph[vStart];

    if(debugOut) {
        debugOut << "Partial path computation for P" << startVertex.id << "\n";
    }

    // The vertices we encounter when following the reads.
    vector<vertex_descriptor> verticesEncountered;

    // The transitions we encounter when following the reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitionsEncountered;

    // Loop over JourneyEntry's of the start vertex.
    for(const JourneyEntry& journeyEntry: startVertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;

        // Store the vertices encountered in the journey of this read.
        const auto journey = journeys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            const vertex_descriptor v = journey[position];
            if(v != null_vertex()) {
                verticesEncountered.push_back(v);
            }
        }

        // Also store the transitions.
        for(uint64_t position=1; position<journey.size(); position++) {
            const vertex_descriptor v0 = journey[position-1];
            const vertex_descriptor v1 = journey[position];
            if(v0 != null_vertex() and v1 != null_vertex()) {
                transitionsEncountered.push_back(make_pair(v0, v1));
            }
        }
    }

    // Count how many times we encountered each vertex.
    vector<uint64_t> vertexFrequency;
    deduplicateAndCount(verticesEncountered, vertexFrequency);

    // Count how many times we encountered each transition.
    // Keep only the ones that appear at least minLinkCoverage times.
    vector<uint64_t> transitionFrequency;
    deduplicateAndCountWithThreshold(
        transitionsEncountered, transitionFrequency, minLinkCoverage);

    // The transitions we kept define a graph.
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
    Graph graph(verticesEncountered.size());
    for(const auto& p: transitionsEncountered) {
        array<vertex_descriptor, 2> v = {p.first, p.second};
        array<uint64_t, 2> iv;
        for(uint64_t k=0; k<2; k++) {
            const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v[k]);
            SHASTA_ASSERT(q.first != verticesEncountered.end());
            SHASTA_ASSERT(q.second - q.first == 1);
            iv[k] = q.first - verticesEncountered.begin();
        }
        add_edge(iv[0], iv[1], graph);
    }



    // Write the graph.
    if(debugOut) {
        debugOut << "digraph PartialPathGraph_P" << startVertex.id << " {\n";
        BGL_FORALL_VERTICES(iv, graph, Graph) {
            const vertex_descriptor pv = verticesEncountered[iv];
            debugOut << "P" << packedAssemblyGraph[pv].id <<
                " [label=\"P" << packedAssemblyGraph[pv].id << "\\n" << vertexFrequency[iv] << "\"]"
                ";\n";
        }
        BGL_FORALL_EDGES(e, graph, Graph) {
            const Graph::vertex_descriptor iv0 = source(e, graph);
            const Graph::vertex_descriptor iv1 = target(e, graph);
            const vertex_descriptor pv0 = verticesEncountered[iv0];
            const vertex_descriptor pv1 = verticesEncountered[iv1];
            debugOut <<
                "P" << packedAssemblyGraph[pv0].id << "->"
                "P" << packedAssemblyGraph[pv1].id << ";\n";
        }
        debugOut << "}\n";
    }



    // To compute the forward partial path, compute the dominator tree of the graph,
    // with the start vertex as the entrance.
    const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), vStart);
    SHASTA_ASSERT(q.first != verticesEncountered.end());
    SHASTA_ASSERT(q.second - q.first == 1);
    const uint64_t ivStart = q.first - verticesEncountered.begin();
    std::map<uint64_t, uint64_t> predecessorMap;
    shasta::lengauer_tarjan_dominator_tree(
        graph,
        ivStart,
        boost::make_assoc_property_map(predecessorMap));



    // Explicitly construct the forward dominator tree.
    Graph forwardTree(verticesEncountered.size());
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.second;
        const uint64_t iv1 = p.first;
        add_edge(iv0, iv1, forwardTree);
    }



    // Write the forward dominator tree.
    if(debugOut) {
        debugOut << "digraph Forward_Tree_P" << startVertex.id << " {\n";

        // Gather the vertices of the dominator tree.
        vector<vertex_descriptor> dominatorTreeVertices;
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            dominatorTreeVertices.push_back(v0);
            dominatorTreeVertices.push_back(v1);
        }
        deduplicate(dominatorTreeVertices);

        for(const vertex_descriptor v: dominatorTreeVertices) {
            const uint64_t iv = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v).first -
                verticesEncountered.begin();
            const vertex_descriptor pv = verticesEncountered[iv];
            debugOut << "P" << packedAssemblyGraph[pv].id <<
                " [label=\"P" << packedAssemblyGraph[pv].id << "\\n" << vertexFrequency[iv] << "\"]"
                ";\n";

        }
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor pv0 = verticesEncountered[iv0];
            const vertex_descriptor pv1 = verticesEncountered[iv1];
            debugOut <<
                "P" << packedAssemblyGraph[pv0].id << "->" <<
                "P" << packedAssemblyGraph[pv1].id << ";\n";
        }
        debugOut << "}\n";
    }



    // To compute the forward partial path, follow the forward dominator tree.
    startVertex.forwardPartialPath.clear();
    uint64_t iv = ivStart;
    while(true) {

        // Find the out-vertices and sort them by decreasing vertex frequency.
        vector< pair<uint64_t, uint64_t> > outVertices;
        BGL_FORALL_OUTEDGES(iv, e, forwardTree, Graph) {
            const uint64_t iv1 = target(e, forwardTree);
            outVertices.push_back(make_pair(iv1, vertexFrequency[iv1]));
        }
        sort(outVertices.begin(), outVertices.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

        // If there are no out-vertices, the forward path ends here.
        if(outVertices.empty()) {
            break;
        }

        // If the strongest out-vertex is too weak, the forward path ends here.
        if(outVertices.front().second < segmentCoverageThreshold1) {
            break;
        }

        // If the strongest in-vertex loses too much coverage compared to iv, the backward path ends here.
        const uint64_t coverageLoss =
            (outVertices.front().second >= vertexFrequency[iv]) ? 0 :
            (vertexFrequency[iv] - outVertices.front().second);
        if(coverageLoss > segmentCoverageThreshold2) {
            break;
        }

        // In all other cases, we add the strongest out-vertex to the forward path.
        iv = outVertices.front().first;
        startVertex.forwardPartialPath.push_back(verticesEncountered[iv]);
    }

    if(debugOut) {
        debugOut << "Forward partial path for P" << startVertex.id << ":\n";
        for(const vertex_descriptor v: startVertex.forwardPartialPath) {
            debugOut << "P" << packedAssemblyGraph[v].id << " ";
        }
        debugOut << "\n";
    }



    // To compute the backward partial path, compute the dominator tree of the reversed graph,
    // with the start vertex as the entrance.
    predecessorMap.clear();
    shasta::lengauer_tarjan_dominator_tree(
        boost::make_reverse_graph(graph),
        ivStart,
        boost::make_assoc_property_map(predecessorMap));



    // Explicitly construct the backward dominator tree.
    Graph backwardTree(verticesEncountered.size());
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.second;
        const uint64_t iv1 = p.first;
        add_edge(iv0, iv1, backwardTree);
    }



    // Write the backward dominator tree.
    if(debugOut) {
        debugOut << "digraph Backward_Tree_P" << startVertex.id << " {\n";

        // Gather the vertices of the dominator tree.
        vector<vertex_descriptor> dominatorTreeVertices;
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            dominatorTreeVertices.push_back(v0);
            dominatorTreeVertices.push_back(v1);
        }
        deduplicate(dominatorTreeVertices);

        for(const vertex_descriptor v: dominatorTreeVertices) {
            const uint64_t iv = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v).first -
                verticesEncountered.begin();
            const vertex_descriptor pv = verticesEncountered[iv];
            debugOut << "P" << packedAssemblyGraph[pv].id <<
                " [label=\"P" << packedAssemblyGraph[pv].id << "\\n" << vertexFrequency[iv] << "\"]"
                ";\n";

        }
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor pv0 = verticesEncountered[iv0];
            const vertex_descriptor pv1 = verticesEncountered[iv1];
            debugOut <<
                "P" << packedAssemblyGraph[pv0].id << "->" <<
                "P" << packedAssemblyGraph[pv1].id << ";\n";
        }
        debugOut << "}\n";
    }



    // To compute the backward partial path, follow the backward dominator tree.
    startVertex.backwardPartialPath.clear();
    iv = ivStart;
    while(true) {

        // Find the out-vertices and sort them by decreasing vertex frequency.
        vector< pair<uint64_t, uint64_t> > outVertices;
        BGL_FORALL_OUTEDGES(iv, e, backwardTree, Graph) {
            const uint64_t iv1 = target(e, backwardTree);
            outVertices.push_back(make_pair(iv1, vertexFrequency[iv1]));
        }
        sort(outVertices.begin(), outVertices.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

        // If there are no out-vertices, the backward path ends here.
        if(outVertices.empty()) {
            break;
        }

        // If the strongest out-vertex is too weak, the backward path ends here.
        if(outVertices.front().second < segmentCoverageThreshold1) {
            break;
        }

        // If the strongest in-vertex loses too much coverage compared to iv, the backward path ends here.
        const uint64_t coverageLoss =
            (outVertices.front().second >= vertexFrequency[iv]) ? 0 :
            (vertexFrequency[iv] - outVertices.front().second);
        if(coverageLoss > segmentCoverageThreshold2) {
            break;
        }

        // In all other cases, we add the strongest out-vertex to the forward path.
        iv = outVertices.front().first;
        startVertex.backwardPartialPath.push_back(verticesEncountered[iv]);
    }

    if(debugOut) {
        debugOut << "Backward partial path for P" << startVertex.id << ":\n";
        for(const vertex_descriptor v: startVertex.backwardPartialPath) {
            debugOut << "P" << packedAssemblyGraph[v].id << " ";
        }
        debugOut << "\n";
    }
}



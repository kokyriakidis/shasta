// Shasta.
#include "mode3a-PackedAssemblyGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "removeReciprocalEdges.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>
#include "utility.hpp"



PackedAssemblyGraph::PackedAssemblyGraph(
    AssemblyGraph& assemblyGraph,
    uint64_t minLinkCoverage,
    uint64_t minMarkerCount,
    double minJaccard,
    uint64_t threadCount) :
    assemblyGraph(assemblyGraph)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    createVertices(minLinkCoverage, minMarkerCount);
    computeJourneys();
    createEdges(minJaccard, threadCount);
    writeGraphviz(minJaccard);

    cout << "The PackedAssemblyGraph has " << num_vertices(packedAssemblyGraph) <<
        " vertices and " << num_edges(packedAssemblyGraph) << " edges." << endl;
}



void PackedAssemblyGraph::createVertices(
    uint64_t minLinkCoverage,
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
            if(assemblyGraph.edgeCoverage(e) >= minLinkCoverage) {
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
            if(assemblyGraph.edgeCoverage(e) >= minLinkCoverage) {
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
}



void PackedAssemblyGraph::createEdges(double minJaccard, uint64_t threadCount)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    cout << timestamp << "PackedAssemblyGraph::createEdges begins." << endl;

    // Compute candidate edges by following the journeys.
    vector< pair<vertex_descriptor, vertex_descriptor> > candidateEdges;
    for(const vector<vertex_descriptor>& journey: journeys) {
        for(uint64_t i=1; i<journey.size(); i++) {
            const vertex_descriptor v1 = journey[i];
            for(uint64_t j=0; j<i; j++) {
                const vertex_descriptor v0 = journey[j];
                if(v0 != v1) {
                    candidateEdges.push_back(make_pair(v0, v1));
                }
            }
        }
    }
    cout << "Before deduplication, there are " << candidateEdges.size() <<
        " candidate edges." << endl;
    deduplicate(candidateEdges);
    cout << "After deduplication, there are " << candidateEdges.size() <<
        " candidate edges." << endl;

    // To compute Jaccard similarities, we need OrientedReadIds
    // to be available in AssemblyGraph vertices.
    assemblyGraph.computeVertexOrientedReadIds(threadCount);



    // Create the edges.
    vector<OrientedReadId> commonOrientedReadIds;
    for(const auto& p: candidateEdges) {

        // Access the PackedAssemblyGraph vertices for this candidate edge.
        const vertex_descriptor pv0 = p.first;
        const vertex_descriptor pv1 = p.second;
        const PackedAssemblyGraphVertex& vertex0 = packedAssemblyGraph[pv0];
        const PackedAssemblyGraphVertex& vertex1 = packedAssemblyGraph[pv1];

        // Access the relevant AssemblyGraph vertices.
        const auto av0 = vertex0.assemblyGraphVertices.back();
        const auto av1 = vertex1.assemblyGraphVertices.front();

        // If the Jaccard similarity is sufficiently high, create an edge.
        const double jaccard = assemblyGraph.computeJaccard(av0, av1, commonOrientedReadIds);
        if(jaccard >= minJaccard) {
            add_edge(pv0, pv1, PackedAssemblyGraphEdge({jaccard}), packedAssemblyGraph);
        }
    }


    // Cleanup.
    assemblyGraph.clearVertexOrientedReadIds();

    cout << timestamp << "PackedAssemblyGraph::createEdges ends." << endl;
}



string PackedAssemblyGraph::vertexStringId(vertex_descriptor pv) const
{
    const PackedAssemblyGraph& packedAssemblyGraph = *this;
    const PackedAssemblyGraphVertex& packedAssemblyGraphVertex = packedAssemblyGraph[pv];
    const AssemblyGraph::vertex_descriptor v0 = packedAssemblyGraphVertex.assemblyGraphVertices.front();
    const AssemblyGraph::vertex_descriptor v1 = packedAssemblyGraphVertex.assemblyGraphVertices.back();

    return assemblyGraph.vertexStringId(v0) + "_" + assemblyGraph.vertexStringId(v1);
}



void PackedAssemblyGraph::writeGraphviz(double minJaccard) const
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
            assemblyGraph.vertexStringId(av1) <<
            "\"]";
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, packedAssemblyGraph, PackedAssemblyGraph) {
        const vertex_descriptor pv0 = source(e, packedAssemblyGraph);
        const vertex_descriptor pv1 = target(e, packedAssemblyGraph);

        // Color it red if jaccard=minJaccard, green if jaccard=1.
        const double jaccard = packedAssemblyGraph[e].jaccard;
        const double ratio = (jaccard - minJaccard) / (1. - minJaccard);
        const double hue = ratio / 3.;

        dot <<
            packedAssemblyGraph[pv0].id << "->" <<
            packedAssemblyGraph[pv1].id <<
            " [color=\"" << hue << ",1,1\"]"
            ";\n";
    }

    dot << "}\n";
}


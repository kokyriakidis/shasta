// Shasta.
#include "mode3a-PackedAssemblyGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>
#include "utility.hpp"



PackedAssemblyGraph::PackedAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    uint64_t minLinkCoverage) :
    assemblyGraph(assemblyGraph)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    createVertices( minLinkCoverage);
    computeJourneys();

    cout << "The PackedAssemblyGraph has " << num_vertices(packedAssemblyGraph) <<
        " vertices and " << num_edges(packedAssemblyGraph) << " edges." << endl;
}



void PackedAssemblyGraph::createVertices(uint64_t minLinkCoverage)
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

    // Each linear chain generates a vertex of the PackedAssemblyGraph.
    // We are not interested in circular chains, so
    // each linear sequence begins at a vertex whose parent does not appear in the vertex map.
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
        for(const vertex_descriptor u: vertex.path) {
            cout << " " << assemblyGraph.vertexStringId(u);
        }
        cout << "\n";
#endif
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


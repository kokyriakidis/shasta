// Shasta.
#include "mode3a-PackedAssemblyGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "deduplicate.hpp"
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
    uint64_t minLinkCoverage1,
    uint64_t minLinkCoverage2,
    uint64_t minMarkerCount) :
    assemblyGraph(assemblyGraph)
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    createVertices(minLinkCoverage1, minMarkerCount);
    computeJourneys();
    createEdges(minLinkCoverage2);
    removeRoundTripEdges();
    writeGraphviz();
    writeJourneys();

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
            assemblyGraph.vertexStringId(av1) << "\"]";
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



void PackedAssemblyGraph::removeRoundTripEdges()
{
    PackedAssemblyGraph& packedAssemblyGraph = *this;
    vector<edge_descriptor> edgesTobeRemoved;
    BGL_FORALL_EDGES(e, packedAssemblyGraph, PackedAssemblyGraph) {
        const vertex_descriptor v0 = source(e, packedAssemblyGraph);
        const vertex_descriptor v1 = target(e, packedAssemblyGraph);

        bool reverseEdgeExists = false;
        tie(ignore, reverseEdgeExists) = boost::edge(v1, v0, packedAssemblyGraph);
        if(reverseEdgeExists) {
            edgesTobeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesTobeRemoved) {
        boost::remove_edge(e, packedAssemblyGraph);
    }

}

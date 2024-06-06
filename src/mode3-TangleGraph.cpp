// Shasta.
#include "mode3-TangleGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "tuple.hpp"



TangleGraph::TangleGraph(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    TangleGraph& tangleGraph = *this;
    const uint64_t orientedReadCount = assemblyGraph.orientedReadIds.size();
    const uint64_t primaryCount = assemblyGraph.markerGraphEdgeIds.size();

    cout << "Creating a TangleGraph with " << orientedReadCount <<
        " oriented reads and " << primaryCount <<
        " primary marker graph edges." << endl;

    // Find the chains that each marker graph edge appears internally in.
    // Usually there will be 1 or 0.
    // Indexed by marker graph edge index in assemblyGraph.markerGraphEdgeIds.
    vector< vector<AssemblyGraph::edge_descriptor> >
        markerGraphEdgeIndexToChainsTable(primaryCount);
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = assemblyGraph[e];
        SHASTA_ASSERT(bubbleChain.isSimpleChain());
        const Chain& chain = bubbleChain.getOnlyChain();
        SHASTA_ASSERT(chain.size() >= 2);

        for(uint64_t position=1; position<chain.size()-1; position++) {
            const MarkerGraphEdgeId markerGraphEdgeId = chain[position];
            const uint64_t markerGraphEdgeIndex = assemblyGraph.getMarkerGraphEdgeIndex(markerGraphEdgeId);
            markerGraphEdgeIndexToChainsTable[markerGraphEdgeIndex].push_back(e);
        }
    }

    // Histogram the sizes of markerGraphEdgeIndexToChainsTable.
    {
        vector<uint64_t> histogram;
        for(const auto& v: markerGraphEdgeIndexToChainsTable) {
            const uint64_t size = v.size();
            if(size >= histogram.size()) {
                histogram.resize(size+1, 0);
            }
            ++histogram[size];
        }

        cout << "Size histogram for markerGraphEdgeIndexToChainsTable:" << endl;
        for(uint64_t size=0; size<histogram.size(); size++) {
            cout << size << " " << histogram[size] << endl;
        }
    }


    // A vector that will contain the vertex that contains each MarkerGraphEdgeId.
    // Indexed by marker graph edge index in assemblyGraph.markerGraphEdgeIds.
    vector<vertex_descriptor> markerGraphEdgeIdVertexTable(primaryCount, null_vertex());



    // Generate a vertex for each Chain.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {

        const vertex_descriptor v = boost::add_vertex(TangleGraphVertex(e), tangleGraph);

        const BubbleChain& bubbleChain = assemblyGraph[e];
        SHASTA_ASSERT(bubbleChain.isSimpleChain());
        const Chain& chain = bubbleChain.getOnlyChain();
        SHASTA_ASSERT(chain.size() >= 2);

        for(uint64_t position=1; position<chain.size()-1; position++) {
            const MarkerGraphEdgeId markerGraphEdgeId = chain[position];
            const uint64_t markerGraphEdgeIndex = assemblyGraph.getMarkerGraphEdgeIndex(markerGraphEdgeId);

            if(markerGraphEdgeIndexToChainsTable[markerGraphEdgeIndex].size() == 1) {
                SHASTA_ASSERT(markerGraphEdgeIdVertexTable[markerGraphEdgeIndex] == null_vertex());
                markerGraphEdgeIdVertexTable[markerGraphEdgeIndex] = v;
            }
        }
    }
    cout << "Created tangle graph edges for " << num_edges(assemblyGraph) << " chains." << endl;


    // Now create a vertex for each primary marker graph edge that is not internal to any Chains.
    for(uint64_t markerGraphEdgeIndex=0; markerGraphEdgeIndex<primaryCount; markerGraphEdgeIndex++) {
        if(markerGraphEdgeIndexToChainsTable[markerGraphEdgeIndex].size() == 0) {
            const vertex_descriptor v = boost::add_vertex(TangleGraphVertex(markerGraphEdgeIndex, MarkerGraphEdgeId()), tangleGraph);
            SHASTA_ASSERT(markerGraphEdgeIdVertexTable[markerGraphEdgeIndex] == null_vertex());
            markerGraphEdgeIdVertexTable[markerGraphEdgeIndex] = v;
        }
    }
    cout << "Created vertices for primary marker graph edges not internal to any chains." << endl;



    // Find the vertices encountered by each oriented read.
    vector< vector<vertex_descriptor> > orientedReadVertices(orientedReadCount);
    for(uint64_t markerGraphEdgeIndex=0; markerGraphEdgeIndex<primaryCount; markerGraphEdgeIndex++) {
        const vertex_descriptor v = markerGraphEdgeIdVertexTable[markerGraphEdgeIndex];
        if(v == null_vertex()) {
            continue;
        }

        const MarkerGraphEdgeId markerGraphEdgeId = assemblyGraph.markerGraphEdgeIds[markerGraphEdgeIndex];
        const auto markerIntervals = assemblyGraph.assembler.markerGraph.edgeMarkerIntervals[markerGraphEdgeId];

        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint64_t orientedReadIndex = assemblyGraph.getOrientedReadIndex(orientedReadId);
            orientedReadVertices[orientedReadIndex].push_back(v);
        }
    }

    // Deduplicate and count.
    vector< vector<uint64_t> > count(orientedReadCount);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        deduplicateAndCount(orientedReadVertices[i], count[i]);
    }



    // Now we can create vertices  corresponding to oriented reads and edges.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadCount; orientedReadIndex++) {
        const vertex_descriptor v0 = boost::add_vertex(TangleGraphVertex(orientedReadIndex, OrientedReadId()), tangleGraph);

        const auto& vv = orientedReadVertices[orientedReadIndex];
        const auto& cc = count[orientedReadIndex];
        SHASTA_ASSERT(vv.size() == cc.size());

        for(uint64_t i=0; i<vv.size(); i++) {
            const vertex_descriptor v1 = vv[i];
            const uint64_t coverage = cc[i];
            add_edge(v0, v1, {coverage}, tangleGraph);
        }
    }

    cout << "The tangle graph has " << num_vertices(tangleGraph) << " vertices and " <<
        num_edges(tangleGraph) << " edges." << endl;


    // Write only the edges that don't involve Chains.
    ofstream dot("TangleGraph.dot");
    dot << "graph TangleGraph {\n";
    BGL_FORALL_EDGES(e, tangleGraph, TangleGraph) {
        const vertex_descriptor v0 = source(e, tangleGraph);
        const TangleGraphVertex& vertex0 = tangleGraph[v0];
        if(vertex0.type == TangleGraphVertexType::Chain) {
            continue;
        }

        const vertex_descriptor v1 = target(e, tangleGraph);
        const TangleGraphVertex& vertex1 = tangleGraph[v1];
        if(vertex1.type == TangleGraphVertexType::Chain) {
            continue;
        }

        dot << "\"";
        if(vertex0.type == TangleGraphVertexType::OrientedRead) {
            dot << assemblyGraph.orientedReadIds[vertex0.orientedReadIndex];
        } else {
            dot << assemblyGraph.markerGraphEdgeIds[vertex0.markerGraphEdgeIndex];
        }
        dot << "\"--\"";
        if(vertex1.type == TangleGraphVertexType::OrientedRead) {
            dot << assemblyGraph.orientedReadIds[vertex1.orientedReadIndex];
        } else {
            dot << assemblyGraph.markerGraphEdgeIds[vertex1.markerGraphEdgeIndex];
        }
        dot << "\";\n";
    }
    dot << "}\n";
}


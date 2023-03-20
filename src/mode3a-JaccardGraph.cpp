// Shasta.
#include "mode3a-JaccardGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>



JaccardGraph::JaccardGraph(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{}



JaccardGraph::vertex_descriptor JaccardGraph::addVertex(
    AssemblyGraphBaseClass::vertex_descriptor av)
{
    JaccardGraph& jaccardGraph = *this;
    const vertex_descriptor v = add_vertex(JaccardGraphVertex({av}), jaccardGraph);
    vertexMap.insert({av, v});
    return v;
}



JaccardGraph::edge_descriptor JaccardGraph::addEdge(
    AssemblyGraphBaseClass::vertex_descriptor av0,
    AssemblyGraphBaseClass::vertex_descriptor av1,
    double jaccard)
{
    JaccardGraph& jaccardGraph = *this;

    auto it0 = vertexMap.find(av0);
    auto it1 = vertexMap.find(av1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    edge_descriptor e;
    bool edgeWasAdded = false;
    tie(e, edgeWasAdded) = add_edge(v0, v1, JaccardGraphEdge({jaccard}), jaccardGraph);
    SHASTA_ASSERT(edgeWasAdded);

    return e;
}



string JaccardGraph::vertexStringId(vertex_descriptor v) const
{
    const JaccardGraph& jaccardGraph = *this;
    const JaccardGraphVertex& vertex = jaccardGraph[v];
    return assemblyGraph.vertexStringId(vertex.av);
}



void JaccardGraph::makeKnn(uint64_t m)
{
    JaccardGraph& jaccardGraph = *this;

    // Mark all edges as not to be kept.
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        jaccardGraph[e].keep = false;
    }

    // Mark edges to be kept.
    vector< pair<edge_descriptor, double> > edgeTable;
    BGL_FORALL_VERTICES(v0, jaccardGraph, JaccardGraph) {

        // Children.
        edgeTable.clear();
        BGL_FORALL_OUTEDGES(v0, e, jaccardGraph, JaccardGraph) {
            const double jaccard = jaccardGraph[e].jaccard;
            edgeTable.push_back(make_pair(e, jaccard));
        }
        sort(edgeTable.begin(), edgeTable.end(),
            OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
        if(edgeTable.size() > m) {
            edgeTable.resize(m);
        }
        for(const auto& p: edgeTable) {
            jaccardGraph[p.first].keep = true;
        }

        // Parents.
        edgeTable.clear();
        BGL_FORALL_INEDGES(v0, e, jaccardGraph, JaccardGraph) {
            const double jaccard = jaccardGraph[e].jaccard;
            edgeTable.push_back(make_pair(e, jaccard));
        }
        sort(edgeTable.begin(), edgeTable.end(),
            OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
        if(edgeTable.size() > m) {
            edgeTable.resize(m);
        }
        for(const auto& p: edgeTable) {
            jaccardGraph[p.first].keep = true;
        }
    }


    // Remove all the edges not marked as to be kept.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        if(not jaccardGraph[e].keep) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, jaccardGraph);
    }

}


void JaccardGraph::writeGraphviz(const string& fileName, double minJaccard) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, minJaccard);
}



void JaccardGraph::writeGraphviz(ostream& dot, double minJaccard) const
{
    const JaccardGraph& jaccardGraph = *this;

    dot << "digraph JaccardGraph {\n";

    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        dot << "\"" << vertexStringId(v) << "\";\n";
    }


    dot << std::setprecision(2);
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const vertex_descriptor v0 = source(e, jaccardGraph);
        const vertex_descriptor v1 = target(e, jaccardGraph);

        // Color it so jaccard=1 is green, jaccard=minJaccard is red.
        const double jaccard = jaccardGraph[e].jaccard;
        const double ratio = (jaccard - minJaccard) / (1. - minJaccard);
        const double hue = ratio / 3.;

        dot << "\"" << vertexStringId(v0) << "\"->\"" << vertexStringId(v1) <<
            "\" [color=\"" << hue << ",1,1\""
            // " label=\"" << jaccard << "\""
            "]"
            ";\n";
    }
    dot << "}\n";

}



// Compute large connected components.
// The threshold is total number of bases for all vertices
// of a connected component.
void JaccardGraph::computeConnectedComponents(
    uint64_t minBaseCount,
    vector< shared_ptr<JaccardGraph> >& componentGraphs
)
{
    JaccardGraph& jaccardGraph = *this;

    // We cannot use boost::connected_components because it
    // only works for undirected graphs.

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }
    const uint64_t n = uint64_t(vertexMap.size());

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Loop over all edges to compute connected components.
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const vertex_descriptor v0 = source(e, jaccardGraph);
        const vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t iv0 = vertexIndexMap[v0];
        const uint64_t iv1 = vertexIndexMap[v1];
        disjointSets.union_set(iv0, iv1);
    }

    // Gather the vertices in each connected component.
    vector< vector<vertex_descriptor> > components(n);
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        const uint64_t iv = vertexIndexMap[v];
        const uint64_t componentId = disjointSets.find_set(iv);
        components[componentId].push_back(v);
    }



    // Get the componentId of the components we want to keep.
    vector< pair<uint64_t, uint64_t> > componentTable; // (componentId, size in bases)
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const vector<vertex_descriptor>& component = components[componentId];
        if(component.empty()) {
            continue;
        }

        // Compute the total number of bases.
        uint64_t totalBaseCount = 0;
        for(const vertex_descriptor v: component) {
            const AssemblyGraph::vertex_descriptor av = jaccardGraph[v].av;
            const uint64_t segmentId = assemblyGraph[av].segmentId;
            const uint64_t baseCount = assemblyGraph.packedMarkerGraph.segmentSequences[segmentId].size();
            totalBaseCount += baseCount;
        }

        if(totalBaseCount >= minBaseCount) {
            componentTable.push_back(make_pair(componentId, totalBaseCount));
        }
    }

    // Order the connected components by decreasing size.
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

    // Create a JaccardGraph for each of the connected components we want to keep.
    componentGraphs.clear();
    cout << "Connected components of the Jaccard graph, with their total length "
        " in bases (only primary vertices are counted)." << endl;
    uint64_t newComponentId = 0;
    for(const auto& p: componentTable) {
        const uint64_t oldComponentId = p.first;
        const vector<vertex_descriptor>& component = components[oldComponentId];
        const uint64_t totalBaseCount = p.second;
        cout << newComponentId << " " << totalBaseCount << endl;

        // Create a new JaccardGraph for this component.
        const shared_ptr<JaccardGraph> componentGraphPointer =
            make_shared<JaccardGraph>(assemblyGraph);
        componentGraphs.push_back(componentGraphPointer);
        JaccardGraph& componentGraph = *componentGraphPointer;

        // Add the vertices.
        for(const vertex_descriptor jv: component) {
            const AssemblyGraph::vertex_descriptor av = jaccardGraph[jv].av;
            componentGraph.addVertex(av);
        }

        // Add the edges.
        BGL_FORALL_VERTICES(cv0, componentGraph, JaccardGraph) {

            // The corresponding vertex in the AssemblyGraph.
            const AssemblyGraph::vertex_descriptor av0 = componentGraph[cv0].av;

            // The corresponding vertex in the global JaccardGraph.
            // cout << assemblyGraph.vertexStringId(av0) << endl;
            SHASTA_ASSERT(vertexMap.contains(av0));
            const vertex_descriptor jv0 = vertexMap[av0];

            BGL_FORALL_OUTEDGES(jv0, e, jaccardGraph, JaccardGraph) {
                const vertex_descriptor jv1 = target(e, jaccardGraph);

                // The corresponding vertex in the AssemblyGraph.
                const AssemblyGraph::vertex_descriptor av1 = jaccardGraph[jv1].av;

                // The corresponding vertex in this componentGraph.
                const vertex_descriptor cv1 = componentGraph.vertexMap[av1];
                SHASTA_ASSERT(componentGraph.vertexMap.contains(av1));

                boost::add_edge(cv0, cv1, jaccardGraph[e], componentGraph);
            }

        }

        ++newComponentId;
    }

}



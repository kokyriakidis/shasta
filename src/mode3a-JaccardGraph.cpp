// Shasta.
#include "mode3a-JaccardGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

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


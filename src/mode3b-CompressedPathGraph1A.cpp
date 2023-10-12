// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "findLinearChains.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<CompressedPathGraph1A>;



CompressedPathGraph1A::CompressedPathGraph1A(
    const PathGraph1& graph,
    uint64_t componentId,
    const Assembler& assembler) :
    MultithreadedObject<CompressedPathGraph1A>(*this),
    graph(graph),
    componentId(componentId),
    assembler(assembler)
{
    CompressedPathGraph1A& cGraph = *this;

    cout << "PathGraph1 connected component " << componentId <<
        " has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;

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

    // Find linear chains in the PathGraph1 after transitive reduction.
    vector< vector<PathGraph1::edge_descriptor> > chains;
    findLinearChains(filteredGraph, 0, chains);

    // Each chain generates an edge.
    // Vertices are added as needed.
    std::map<PathGraph1::vertex_descriptor, CompressedPathGraph1A::vertex_descriptor> vertexMap;
    for(const vector<PathGraph1::edge_descriptor>& chain: chains) {
        const PathGraph1::vertex_descriptor v0 = source(chain.front(), graph);
        const PathGraph1::vertex_descriptor v1 = target(chain.back(), graph);
        const vertex_descriptor cv0 = getCompressedVertex(v0);
        const vertex_descriptor cv1 = getCompressedVertex(v1);
        add_edge(cv0, cv1, {chain}, cGraph);
    };

    cout << "The CompressedPathGraph1A for PathGraph1 connected component " << componentId <<
        " has " << num_vertices(cGraph) << " vertices and " <<
        num_edges(cGraph) << " edges." << endl;
    writeGraphviz("Initial");


}



// Get the vertex_descriptor corresponding to a PathGraph1::vertex_descriptor,
// adding a vertex if necessary.
CompressedPathGraph1A::vertex_descriptor CompressedPathGraph1A::getCompressedVertex(PathGraph1::vertex_descriptor v)
{
    CompressedPathGraph1A& cGraph = *this;

    auto it = vertexMap.find(v);
    if(it == vertexMap.end()) {
        const vertex_descriptor cv = add_vertex({v}, cGraph);
        vertexMap.insert({v, cv});
        return cv;
    } else {
        return it->second;
    }
}



PathGraph1::vertex_descriptor CompressedPathGraph1A::firstUncompressedVertex(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;
    SHASTA_ASSERT(not chain.empty());

    const PathGraph1::edge_descriptor e = chain.front();
    const PathGraph1::vertex_descriptor v = source(e, graph);

    // Sanity check.
    auto it = vertexMap.find(v);
    SHASTA_ASSERT(it != vertexMap.end());
    const vertex_descriptor cv = it->second;
    SHASTA_ASSERT(cv == source(ce, cGraph));

    return v;
}



PathGraph1::vertex_descriptor CompressedPathGraph1A::lastUncompressedVertex(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;
    SHASTA_ASSERT(not chain.empty());

    const PathGraph1::edge_descriptor e = chain.back();
    const PathGraph1::vertex_descriptor v = target(e, graph);

    // Sanity check.
    auto it = vertexMap.find(v);
    SHASTA_ASSERT(it != vertexMap.end());
    const vertex_descriptor cv = it->second;
    SHASTA_ASSERT(cv == target(ce, cGraph));

    return v;
}



PathGraph1::vertex_descriptor CompressedPathGraph1A::firstInternalUncompressedVertex(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;
    SHASTA_ASSERT(not chain.empty());

    const PathGraph1::edge_descriptor e = chain.front();
    const PathGraph1::vertex_descriptor v = target(e, graph);

    return v;
}



PathGraph1::vertex_descriptor CompressedPathGraph1A::lastInternalUncompressedVertex(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;
    SHASTA_ASSERT(not chain.empty());

    const PathGraph1::edge_descriptor e = chain.back();
    const PathGraph1::vertex_descriptor v = source(e, graph);

    return v;
}



void CompressedPathGraph1A::writeGraphviz(const string& fileNamePrefix) const
{
    const CompressedPathGraph1A& cGraph = *this;

    ofstream dot("CompressedPathGraph1A-" + to_string(componentId) + "-" + fileNamePrefix + ".dot");
    dot << "digraph CompressedPathGraph1A_" << componentId << " {\n";

    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1A) {
        const PathGraph1::vertex_descriptor v = cGraph[cv].v;
        dot << graph[v].edgeId << ";\n";
    }

    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1A) {
        const PathGraph1::vertex_descriptor v0 = firstUncompressedVertex(ce);
        const PathGraph1::vertex_descriptor v1 = lastUncompressedVertex(ce);
        const auto& chain = cGraph[ce].chain;
        SHASTA_ASSERT(not chain.empty());
        dot << graph[v0].edgeId << "->" << graph[v1].edgeId;
        if(chain.size() == 1) {
            // No label.
        } else if(chain.size() == 2) {
            dot << " [label=\"" <<
            graph[firstInternalUncompressedVertex(ce)].edgeId <<
            "\"]";
        } else if(chain.size() == 3) {
            dot << " [label=\"" <<
            graph[firstInternalUncompressedVertex(ce)].edgeId << "\\n" <<
            graph[lastInternalUncompressedVertex(ce)].edgeId <<
            "\"]";
        } else {
            dot << " [label=\"" <<
            graph[firstInternalUncompressedVertex(ce)].edgeId << "\\n...\\n" <<
            graph[lastInternalUncompressedVertex(ce)].edgeId <<
            "\"]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}



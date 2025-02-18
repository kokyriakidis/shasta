#pragma once

// Boost libraries.
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/iteration_macros.hpp>


//Standard library.
#include <cstdint.hpp>
#include <map>
#include <vector.hpp>



namespace shasta {

    template<class Graph, class Subgraph> class InducedGraphIsomorphismsCallback;

    template<class Graph, class Subgraph>
        void inducedSubgraphIsomorphisms(
        const Graph&,
        const Subgraph&,
        vector< vector<typename Graph::vertex_descriptor> >& isomorphisms);
}



template<class Graph, class Subgraph> class shasta::InducedGraphIsomorphismsCallback {
public:

    InducedGraphIsomorphismsCallback(
        uint64_t subGraphVertexCount,
        vector< vector<typename Graph::vertex_descriptor> >& isomorphisms) :
        subGraphVertexCount(subGraphVertexCount),
        isomorphisms(isomorphisms)
        {}

    template<class Map1, class Map2> bool operator()(
        const Map1& map1,
        const Map2&) const
    {
        isomorphisms.resize(isomorphisms.size() + 1);
        vector<typename Graph::vertex_descriptor>& isomorphism = isomorphisms.back();
        isomorphism.resize(subGraphVertexCount);
        for(uint64_t i=0; i<subGraphVertexCount; i++) {
            isomorphism[i] = map1[i];
        }
        return true;
    }
private:
    uint64_t subGraphVertexCount;
    vector< vector<typename Graph::vertex_descriptor> >& isomorphisms;
};



// Find induced subgraph isomorphism.
// The Subgraph type must use vecS for its VertexList.
template<class Graph, class Subgraph>
    void shasta::inducedSubgraphIsomorphisms(
    const Graph& graph,
    const Subgraph& subgraph,
    vector< vector<typename Graph::vertex_descriptor> >& isomorphisms)
{
    using namespace boost;
    isomorphisms.clear();

    // Prepare a callback object for vf2_subgraph_iso.
    InducedGraphIsomorphismsCallback<Graph, Subgraph>
        callback(num_vertices(subgraph), isomorphisms);

    // Create a vertex index map for the Graph.
    // For the Subgraph we don't need a vertex index map:
    // we use vertex_index because Subgraph uses vecS for its VertexList.
    std::map<typename Graph::vertex_descriptor, uint64_t> vertexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        vertexMap.insert({v, vertexIndex++});
    }
    auto vertexIndexMap = make_assoc_property_map(vertexMap);

    // Use vf2_subgraph_iso to find the isomorphisms.
    // The callback will be called once for each isomorphism found.
    vf2_subgraph_iso(subgraph, graph, callback,
        get(vertex_index, subgraph),
        vertexIndexMap,
        vertex_order_by_mult(subgraph),
        always_equivalent(),
        always_equivalent());
}

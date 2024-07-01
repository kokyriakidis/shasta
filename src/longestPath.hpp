#ifndef SHASTA_LONGEST_PATH_HPP
#define SHASTA_LONGEST_PATH_HPP

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard library.
#include "algorithm.hpp"
#include <map>
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    template<class Graph> void longestPath(
        const Graph &graph,
        vector<typename Graph::vertex_descriptor>& longestPath);
    template<class Graph> void longestPath1(
        const Graph &graph,
        vector<typename Graph::vertex_descriptor>& longestPath);
    void testLongestPath();
}



// Find the longest path in a directed graph without cycles.
// Class Graph must be a boost::adjacency_list with
// the first three template arguments set to <does not matter, vecS, bidirectionalS>.
// If the graph has cycles, this throws boost::not_a_dag.
// This uses the algorithm described here:
// https://en.wikipedia.org/wiki/Longest_path_problem#Acyclic_graphs
template<class Graph> void shasta::longestPath(
    const Graph &graph,
    vector<typename Graph::vertex_descriptor>& longestPath)
{
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    // using edge_descriptor = typename Graph::edge_descriptor;
    // using edge_iterator = typename Graph::edge_iterator;

    // Check the Graph type.
    // Use C++20 concepts instead.
#if 0
    static_assert(
        std::is_same<typename Graph::out_edge_list_selector, listS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the first template argument set to boost::listS.");
#endif
    static_assert(
        std::is_same<typename Graph::vertex_list_selector, vecS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the second template argument set to boost::vecS.");
    static_assert(
        std::is_same<typename Graph::directed_selector, bidirectionalS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the third template argument set to boost::bidirectionalS.");

    // Use boost topological_sort to get a vector of vertex descriptors
    // in topological order. The output from the boost call is in
    // reverse topological order.
    vector<vertex_descriptor> sortedVertices;
    topological_sort(graph, back_inserter(sortedVertices));
    std::reverse(sortedVertices.begin(), sortedVertices.end());

    // Map to contain the length of the longest path ending at each vertex.
    std::map<vertex_descriptor, uint64_t> lengthMap;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        lengthMap.insert(make_pair(v, 0));
    }

    // Compute the maximum length of a path ending at each vertex.
    for(const vertex_descriptor v: sortedVertices) {
        uint64_t maximumLength = 0;
        BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
            maximumLength = max(maximumLength, lengthMap[source(e, graph)]);
        }
        lengthMap[v] = maximumLength + 1;
    }

    // Find the vertex with the longest length.
    // This will be the end of the longest path.
    vertex_descriptor v = Graph::null_vertex();
    uint64_t maximumLength = 0;
    for(const auto& p: lengthMap) {
        if(p.second > maximumLength) {
            v = p.first;
            maximumLength = p.second;
        }
    }

    // Constuct the path, moving backward from here.
    longestPath.clear();
    longestPath.push_back(v);
    while(true) {
        vertex_descriptor vPrevious = Graph::null_vertex();
        uint64_t maximumLength = 0;
        BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const uint64_t length = lengthMap[v0];
            if(length > maximumLength) {
                vPrevious = v0;
                maximumLength = length;
            }
        }
        if(vPrevious == Graph::null_vertex()) {
            break;
        }
        v = vPrevious;
        longestPath.push_back(v);

    }
    std::reverse(longestPath.begin(), longestPath.end());

}



// More efficient version that uses a vectorinstead of a map to record
// path lengths.
template<class Graph> void shasta::longestPath1(
    const Graph &graph,
    vector<typename Graph::vertex_descriptor>& longestPath)
{
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    // using edge_descriptor = typename Graph::edge_descriptor;
    // using edge_iterator = typename Graph::edge_iterator;

    // Check the Graph type.
    // Use C++20 concepts instead.
    static_assert(
        std::is_same<typename Graph::vertex_list_selector, vecS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the second template argument set to boost::vecS.");
    static_assert(
        std::is_same<typename Graph::directed_selector, bidirectionalS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the third template argument set to boost::bidirectionalS.");

    // Use boost topological_sort to get a vector of vertex descriptors
    // in topological order. The output from the boost call is in
    // reverse topological order.
    vector<vertex_descriptor> sortedVertices;
    topological_sort(graph, back_inserter(sortedVertices));
    std::reverse(sortedVertices.begin(), sortedVertices.end());

    // Vector to contain the length of the longest path ending at each vertex.
    vector<uint64_t> length(num_vertices(graph), 0);

    // Compute the maximum length of a path ending at each vertex.
    for(const vertex_descriptor v: sortedVertices) {
        uint64_t maximumLength = 0;
        BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
            maximumLength = max(maximumLength, length[source(e, graph)]);
        }
        length[v] = maximumLength + 1;
    }

    // Find the vertex with the longest length.
    // This will be the end of the longest path.
    vertex_descriptor v = max_element(length.begin(), length.end()) - length.begin();

    // Construct the path, moving backward from here.
    longestPath.clear();
    longestPath.push_back(v);
    while(true) {
        vertex_descriptor vPrevious = Graph::null_vertex();
        uint64_t maximumLength = 0;
        BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const uint64_t l = length[v0];
            if(l > maximumLength) {
                vPrevious = v0;
                maximumLength = l;
            }
        }
        if(vPrevious == Graph::null_vertex()) {
            break;
        }
        v = vPrevious;
        longestPath.push_back(v);

    }
    std::reverse(longestPath.begin(), longestPath.end());

}

#endif

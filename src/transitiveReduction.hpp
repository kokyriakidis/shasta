#ifndef SHASTA_TRANSITIVE_REDUCTION_HPP
#define SHASTA_TRANSITIVE_REDUCTION_HPP

// Shasta.
#include "deduplicate.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard library.
#include "iterator.hpp"
#include <map>
#include <queue>
#include <set>
#include "vector.hpp"

namespace shasta {

    // Version that requires the graph to use vecS.
    template<class Graph> void transitiveReduction(Graph&);

    // Less performant version without the above requirement.
    template<class Graph> void transitiveReductionAny(Graph&);
}



// Transitive reduction of a directed graph without cycles.
// Class Graph must be a boost::adjacency_list with
// the first three template arguments set to <listS, vecS, directedS or bidirectionalS>.
// If the graph has cycles, this throws boost::not_a_dag.
template<class Graph> void shasta::transitiveReduction(Graph &graph)
    {
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;
    using edge_iterator = typename Graph::edge_iterator;

    // Check the Graph type.
    // Use C++20 concepts instead.
    static_assert(
        std::is_same<typename Graph::out_edge_list_selector, listS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the first template argument set to boost::listS.");
    static_assert(
        std::is_same<typename Graph::vertex_list_selector, vecS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the second template argument set to boost::vecS.");
    static_assert(
        std::is_same<typename Graph::directed_selector, directedS>::value
        or
        std::is_same<typename Graph::directed_selector, bidirectionalS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the third template argument set to boost::directedS or boost::bidirectionalS.");

    // Use boost topological_sort to get a vector of vertex descriptors
    // in reverse toplogical order.
    vector<vertex_descriptor> sortedVertices;
    topological_sort(graph, back_inserter(sortedVertices));

    // Now construct a vector containing the rank of each vertex in topological order.
    vector<uint64_t> vertexRank(num_vertices(graph));
    uint64_t rank = num_vertices(graph);
    for (const vertex_descriptor v : sortedVertices) {
        vertexRank[v] = --rank;
    }

    // Find the edges that should be removed.
    vector<edge_descriptor> edgesToBeRemoved;
    vector<bool> wasVisited(num_vertices(graph), false);
    vector<vertex_descriptor> visitedVertices;
    edge_iterator it, itEnd;
    tie(it, itEnd) = edges(graph);
    while(it != itEnd) {
        edge_iterator itNext = it;
        ++itNext;
        const edge_descriptor e = *it;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Edge e should be removed if there is a path
        // from v0 to v1 that does not use edge e.

        // Do a forward BFS starting at v0 and ending at v1 but:
        // - Don't use edge e in the BFS.
        // - Don't use any vertices that have topological order
        //   greater than the topological order of v1,
        //   because there can be no paths ending at v1
        //   that use these vertices.
        // If the BFS encounters v1, edge e can be removed.

        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        q.push(v0);
        wasVisited[v0] = true;
        visitedVertices.push_back(v0);

        // BFS loop.
        while(not q.empty()) {

            // Dequeue a vertex.
            const vertex_descriptor vv0 = q.front();
            q.pop();

            // Loop over its out-edges.
            BGL_FORALL_OUTEDGES_T(vv0, ee, graph, Graph)
            {

                // Don't use edge e in the BFS.
                if (ee == e) {
                    continue;
                }

                // Get the other vertex in edge ee.
                const vertex_descriptor vv1 = target(ee, graph);

                // If vv1 was already visited in this BFS, skip it.
                if (wasVisited[vv1]) {
                    continue;
                }

                // If vv1 follows v1 in topological order, skip it.
                if (vertexRank[vv1] > vertexRank[v1]) {
                    continue;
                }

                if (vv1 == v1) {
                    // We reached v1. Edge e can be removed and we can stop the BFS.
                    boost::remove_edge(e, graph);
                    q = { };
                    break;
                } else {
                    // Continue the BFS.
                    wasVisited[vv1] = true;
                    visitedVertices.push_back(vv1);
                    q.push(vv1);
                }
            }
        }

        // Prepare for the next iteration.
        it = itNext;

        // Clean up.
        for(const vertex_descriptor v: visitedVertices) {
            wasVisited[v] = false;
        }
        visitedVertices.clear();
    }

}



// Less performant version which works on any acyclic boost directed graph.
template<class Graph> void shasta::transitiveReductionAny(Graph &graph)
    {
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;
    using edge_iterator = typename Graph::edge_iterator;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        vertexIndexMap.insert({v, vertexIndex++});
    }

    // Use boost topological_sort to get a vector of vertex descriptors
    // in reverse topological order.
    vector<vertex_descriptor> sortedVertices;
    topological_sort(
        graph,
        back_inserter(sortedVertices),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));

    // Store the rank of each vertex in topological order.
    std::map<vertex_descriptor, uint64_t> vertexRank;
    uint64_t rank = num_vertices(graph);
    for (const vertex_descriptor v : sortedVertices) {
        vertexRank.insert({v, --rank});
    }

    // Find the edges that should be removed.
    edge_iterator it, itEnd;
    tie(it, itEnd) = edges(graph);
    while(it != itEnd) {
        edge_iterator itNext = it;
        ++itNext;
        const edge_descriptor e = *it;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Edge e should be removed if there is a path
        // from v0 to v1 that does not use edge e.

        // Do a forward BFS starting at v0 and ending at v1 but:
        // - Don't use edge e in the BFS.
        // - Don't use any vertices that have topological order
        //   greater than the topological order of v1,
        //   because there can be no paths ending at v1
        //   that use these vertices.
        // If the BFS encounters v1, edge e can be removed.

        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        q.push(v0);
        std::set<vertex_descriptor> visitedVertices;
        visitedVertices.insert(v0);

        // BFS loop.
        while(not q.empty()) {

            // Dequeue a vertex.
            const vertex_descriptor vv0 = q.front();
            q.pop();

            // Loop over its out-edges.
            BGL_FORALL_OUTEDGES_T(vv0, ee, graph, Graph)
            {

                // Don't use edge e in the BFS.
                if (ee == e) {
                    continue;
                }

                // Get the other vertex in edge ee.
                const vertex_descriptor vv1 = target(ee, graph);

                // If vv1 was already visited in this BFS, skip it.
                if(visitedVertices.contains(vv1)) {
                    continue;
                }

                // If vv1 follows v1 in topological order, skip it.
                if (vertexRank[vv1] > vertexRank[v1]) {
                    continue;
                }

                if (vv1 == v1) {
                    // We reached v1. Edge e can be removed and we can stop the BFS.
                    boost::remove_edge(e, graph);
                    q = { };
                    break;
                } else {
                    // Continue the BFS.
                    visitedVertices.insert(vv1);
                    q.push(vv1);
                }
            }
        }

        // Prepare for the next iteration.
        it = itNext;

    }

}

#endif

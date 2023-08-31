#ifndef SHASTA_LOCAL_TRANSITIVE_REDUCTION_HPP
#define SHASTA_LOCAL_TRANSITIVE_REDUCTION_HPP

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>
#include "vector.hpp"

namespace shasta {
    template<class Graph> void localTransitiveReduction(
        const Graph&,
        uint64_t maxPathLength,
        vector<typename Graph::edge_descriptor>& nonTransitiveReductionEdges);
}



// For each directed edge v0->v1, look for a path that:
// - Starts at v0.
// - Ends at v1.
// - Has length at most maxPathLength;
// - Does not use edge v0>v1;
// If such a path is found, the edge descriptor is added to nonTransitiveReductionEdges.
// The edges that are not in nonTransitiveReductionEdges form a sort of "local transitive reduction"
// of the graph.
template<class Graph> void shasta::localTransitiveReduction(
    const Graph& graph,
    uint64_t maxPathLength,
    vector<typename Graph::edge_descriptor>& nonTransitiveReductionEdges)
{

    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    // using edge_descriptor = typename Graph::edge_descriptor;

    // Check the Graph type.
    static_assert(
        std::is_same<typename Graph::directed_selector, directedS>::value
        or
        std::is_same<typename Graph::directed_selector, bidirectionalS>::value,
        "shasta::transitiveReduction requires an adjacency_list "
        "with the third template argument set to boost::directedS or boost::bidirectionalS.");

    // Loop over all edges v0->v1.
    nonTransitiveReductionEdges.clear();
    BGL_FORALL_EDGES_T(e01, graph, Graph) {
        const vertex_descriptor v0 = source(e01, graph);
        const vertex_descriptor v1 = target(e01, graph);

        // Do a BFS starting at v0, up to a distance maxPathLength.
        // Stop if we encounter v1.

        // The BFS queue.
        std::queue<vertex_descriptor> q;
        q.push(v0);

        // The vertices we encountered so far, with their distance from v0.
        std::map<vertex_descriptor, uint64_t> m;
        m.insert({v0, 0});

        // BFS loop.
        // cout << "BFS loop begins for " << v0 << "->" << v1 << endl;
        while(not q.empty()) {

            // Dequeue a vertex.
            const vertex_descriptor vA = q.front();
            q.pop();
            const auto itA = m.find(vA);
            SHASTA_ASSERT(itA != m.end());
            const uint64_t distanceA = itA->second;
            const uint64_t distanceB = distanceA + 1;
            // cout << "Dequeued " << vA << " at distance " << distanceA << endl;

            // Loop over the out-edges of vA.
            bool endBfs = false;
            BGL_FORALL_OUTEDGES_T(vA, eAB, graph, Graph) {

                // Dont's use e01 in the BFS.
                if(eAB == e01) {
                    continue;
                }

                // If we reached v1, mark e01 as a nonTransitiveReduction edge
                // and stop the BFS.
                const vertex_descriptor vB = target(eAB, graph);
                if(vB == v1) {
                    nonTransitiveReductionEdges.push_back(e01);
                    endBfs = true;
                    // cout << "Reached " << v1 << endl;
                    break;
                }

                // If we already reached this vertex, do nothing.
                if(m.contains(vB)) {
                    continue;
                }

                // If not at maximum distance, enqueue vB.
                if(distanceB < maxPathLength) {
                    q.push(vB);
                    m.insert({vB, distanceB});
                    // cout << "Enqueued " << vB << " at distance " << distanceB << endl;
                }
            }
            if(endBfs) {
                break;
            }
        }
    }
}

#endif

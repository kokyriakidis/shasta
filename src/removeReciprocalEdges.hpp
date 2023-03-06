#ifndef SHASTA_REMOVE_RECIPROCAL_ERDGES_HPP
#define SHASTA_REMOVE_RECIPROCAL_ERDGES_HPP

#include <boost/graph/iteration_macros.hpp>
#include "vector.hpp"

namespace shasta {
    template<class Graph> void removeReciprocalEdges(Graph&);
}



template<class Graph> void shasta::removeReciprocalEdges(Graph& graph)
{
    vector<typename Graph::edge_descriptor> edgesTobeRemoved;

    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const typename Graph::vertex_descriptor v0 = source(e, graph);
        const typename Graph::vertex_descriptor v1 = target(e, graph);

        bool reverseEdgeExists = false;
        tie(ignore, reverseEdgeExists) = boost::edge(v1, v0, graph);
        if(reverseEdgeExists) {
            edgesTobeRemoved.push_back(e);
        }
    }

    for(const typename Graph::edge_descriptor e: edgesTobeRemoved) {
        boost::remove_edge(e, graph);
    }

}
#endif


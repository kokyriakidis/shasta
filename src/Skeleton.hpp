#pragma once


/*****************************************************************

Code to create a "Skeleton" of a directed graph g of type Graph
as follows:
- Compute dominator trees for g using as "start" for the dominator
  tree all vertices with in-degree 0. Gather the union of
  all the tree edges computed in this way and call that Es0.
- Do the same in the opposite direction and call the resulting
  set of tree edges Es1.
- Compute the intersection of Es0 and Es1, call it Es.
- Call Vs the set of G vertices that appear in Es.
- The Skeleton of g uses Vs as its vertices and Es as its edges.

*****************************************************************/

// Shasta.
#include "deduplicate.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dominator_tree.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>

// Standard library.
#include <map>
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {

    // Each vertex of the Skeleton corresponds to a vertex of
    // the original Graph.
    template<class Graph> using SkeletonBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        typename Graph::vertex_descriptor>;
    template<class Graph> class Skeleton;

}



template<class Graph> class shasta::Skeleton : public SkeletonBaseClass<Graph> {
public:

    using Vg =  typename Graph::vertex_descriptor;
    using Vs =  typename Skeleton::vertex_descriptor;


    // Map vertices of the Graph to vertices of the Skeleton.
    // Not all vertices of the Graph have a corresponding vertex in the Skeleton.
    // The inverse map is not needed because we can easily obtain the
    // Graph::vertex_descriptor corresponding to a given
    // Skeleton::vertex_descriptor (it is stored in the Skeleton vertex).
    std::map<Vg, Vs> vertexMap;


    // Construct the Skeleton given the Graph.
    Skeleton(const Graph& graph)
    {
        // Map vertices of graph to integers.
        // This is needed for the computation of dominator tree below.
        std::map<Vg, uint64_t> indexMap;
        uint64_t vertexIndex = 0;
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            indexMap.insert({v, vertexIndex++});
        }
        auto associativeIndexMap = boost::make_assoc_property_map(indexMap);
        const uint64_t vertexCount = vertexIndex;

        // Vectors used below to compute the dominator tree.
        vector<uint64_t> dfNum(vertexCount);
        vector<Vg> parent(vertexCount);
        vector<Vg> verticesByDFNum(vertexCount);

        // Tree pairs found on forward and backward dominator tree.
        vector< pair<Vg, Vg> > forwardPairs;
        vector< pair<Vg, Vg> > backwardPairs;



        // Compute dominator trees using as entrance each of the
        // vertices with zero in-degree.
        BGL_FORALL_VERTICES_T(entrance, graph, Graph) {
            if(in_degree(entrance, graph) != 0) {
                continue;
            }

            // Compute the dominator tree.
            fill(dfNum.begin(), dfNum.end(), invalid<uint64_t>);
            fill(parent.begin(), parent.end(), Graph::null_vertex());
            fill(verticesByDFNum.begin(), verticesByDFNum.end(), Graph::null_vertex());
            std::map<Vg, Vg> predecessorMap;

            boost::lengauer_tarjan_dominator_tree(
                graph,
                entrance,
                boost::make_assoc_property_map(indexMap),
                boost::make_iterator_property_map(dfNum.begin(), associativeIndexMap),
                boost::make_iterator_property_map(parent.begin(), associativeIndexMap),
                verticesByDFNum,
                boost::make_assoc_property_map(predecessorMap));

            for(const auto& p: predecessorMap) {
                const Vg cv0 = p.second;
                const Vg cv1 = p.first;
                forwardPairs.push_back(make_pair(cv0, cv1));
            }
        }



        // Compute dominator trees on the reverse graph using as entrance each of the
        // vertices with zero in-degree on the reverse graph
        // (that is, zero out-degree on the Graph).
        using ReverseGraph = boost::reverse_graph<Graph>;
        ReverseGraph reverseGraph(graph);
        BGL_FORALL_VERTICES_T(entrance, reverseGraph, ReverseGraph) {
            if(in_degree(entrance, reverseGraph) != 0) {
                continue;
            }

            // Compute the dominator tree.
            fill(dfNum.begin(), dfNum.end(), invalid<uint64_t>);
            fill(parent.begin(), parent.end(), Graph::null_vertex());
            fill(verticesByDFNum.begin(), verticesByDFNum.end(), Graph::null_vertex());
            std::map<Vg, Vg> predecessorMap;

            boost::lengauer_tarjan_dominator_tree(
                reverseGraph,
                entrance,
                boost::make_assoc_property_map(indexMap),
                boost::make_iterator_property_map(dfNum.begin(), associativeIndexMap),
                boost::make_iterator_property_map(parent.begin(), associativeIndexMap),
                verticesByDFNum,
                boost::make_assoc_property_map(predecessorMap));

            for(const auto& p: predecessorMap) {
                const Vg v0 = p.first;
                const Vg v1 = p.second;
                backwardPairs.push_back(make_pair(v0, v1));
            }
        }



        // The pairs that appear both in forwardPairs and backwardPairs are the edges of the Skeleton.
        deduplicate(forwardPairs);
        deduplicate(backwardPairs);
        vector< pair<Vg, Vg> > skeletonEdges;
        std::set_intersection(
            forwardPairs.begin(), forwardPairs.end(),
            backwardPairs.begin(), backwardPairs.end(),
            back_inserter(skeletonEdges)
            );

        // The vertices of the Skeleton are the vertices that appear in
        // the bidirectionalPairs.
        vector<Vg> skeletonVertices;
        for(const auto& p: skeletonEdges) {
            skeletonVertices.push_back(p.first);
            skeletonVertices.push_back(p.second);
        }
        deduplicate(skeletonVertices);

        // Create the Skeleton vertices.
        for(const Vg vg: skeletonVertices) {
            const Vs vs = add_vertex(vg, *this);
            vertexMap.insert(make_pair(vg, vs));
        }

        // Create the Skeleton edges.
        for(const auto& p: skeletonEdges) {
            const Vg vg0 = p.first;
            const Vg vg1 = p.second;
            const Vs vs0 = vertexMap[vg0];
            const Vs vs1 = vertexMap[vg1];
            add_edge(vs0, vs1, *this);
        }
    }
};


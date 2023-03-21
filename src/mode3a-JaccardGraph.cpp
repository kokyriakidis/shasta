// Shasta.
#include "mode3a-JaccardGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

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




// Mark the edges to be displayed.
void JaccardGraph::markDisplayEdges()
{
    JaccardGraph& jaccardGraph = *this;

    // Start with only the long path edges
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        JaccardGraphEdge& edge = jaccardGraph[e];
        edge.display = edge.isLongPathEdge;
    }

    // Mark k-NN edges.
    const uint64_t k = 1;
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
        if(edgeTable.size() > k) {
            edgeTable.resize(k);
        }
        for(const auto& p: edgeTable) {
            jaccardGraph[p.first].display = true;
        }

        // Parents.
        edgeTable.clear();
        BGL_FORALL_INEDGES(v0, e, jaccardGraph, JaccardGraph) {
            const double jaccard = jaccardGraph[e].jaccard;
            edgeTable.push_back(make_pair(e, jaccard));
        }
        sort(edgeTable.begin(), edgeTable.end(),
            OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
        if(edgeTable.size() > k) {
            edgeTable.resize(k);
        }
        for(const auto& p: edgeTable) {
            jaccardGraph[p.first].display = true;
        }
    }

}


void JaccardGraph::writeGraphviz(const string& fileName, double minJaccard)
{
    ofstream dot(fileName);
    writeGraphviz(dot, minJaccard);
}



void JaccardGraph::writeGraphviz(ostream& dot, double minJaccard)
{
    const JaccardGraph& jaccardGraph = *this;

    // To simplify the display, we only display a significant subset of the edges.
    markDisplayEdges();

    dot << "digraph JaccardGraph {\n";

    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        dot << "\"" << vertexStringId(v) << "\"";
        if(jaccardGraph[v].isLongPathVertex) {
            dot << "[style=filled fillcolor=pink]";
        }
        dot << ";\n";
    }


    dot << std::setprecision(2);
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraphEdge& edge = jaccardGraph[e];
        if(not edge.display) {
            continue;
        }
        const vertex_descriptor v0 = source(e, jaccardGraph);
        const vertex_descriptor v1 = target(e, jaccardGraph);

        // Color it so jaccard=1 is green, jaccard=minJaccard is red.
        const double jaccard = edge.jaccard;
        const double ratio = (jaccard - minJaccard) / (1. - minJaccard);
        const double H = ratio / 3.;

        double S = 1.;
        double V = 1;
        /*
        if(not edge.isLongPathEdge) {
            S = 0.4;
        }
        */

        dot << "\"" << vertexStringId(v0) << "\"->\"" << vertexStringId(v1) <<
            "\" [color=\"" << H << "," << S << "," << V << "\"";

        if(not edge.isLongPathEdge) {
            dot << " style=dashed";
        }

        // " label=\"" << jaccard << "\"";
        dot << "];\n";
    }
    dot << "}\n";

}



// Compute large connected components.
void JaccardGraph::computeConnectedComponents(
    uint64_t minComponentSize,
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
    vector< pair<uint64_t, uint64_t> > componentTable; // (componentId, size)
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const vector<vertex_descriptor>& component = components[componentId];
        const uint64_t componentSize = component.size();
        if(componentSize >= minComponentSize) {
            componentTable.push_back(make_pair(componentId, componentSize));
        }
    }

    // Order the connected components by decreasing size.
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

    // Create a JaccardGraph for each of the connected components we want to keep.
    componentGraphs.clear();
    uint64_t newComponentId = 0;
    for(const auto& p: componentTable) {
        const uint64_t oldComponentId = p.first;
        const vector<vertex_descriptor>& component = components[oldComponentId];

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



// Compute strongly connected components.
void JaccardGraph::computeStronglyConnectedComponents(
    vector< vector<vertex_descriptor> >& strongComponents)
{
    JaccardGraph& jaccardGraph = *this;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }

    // Compute strongly connected components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        jaccardGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentTable;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        const uint64_t componentId = componentMap[v];
        componentTable[componentId].push_back(v);
    }
    strongComponents.clear();
    for(const auto& p: componentTable) {
        const vector<vertex_descriptor>& component = p.second;
        if(component.size() > 1) {
            strongComponents.push_back(component);
            for(const vertex_descriptor v: component) {
                jaccardGraph[v].isCyclic = true;
            }
        }
    }

#if 0
    if(not strongComponents.empty()) {
        cout << "Found " << strongComponents.size() <<
            " strongly connected components with sizes";
        for(const auto& component: strongComponents) {
            cout << " " << component.size();
        }
        cout << endl;
    }
#endif
}

// Remove all vertices that belong to strogly connected components.
void JaccardGraph::removeStronglyConnectedComponents()
{
    JaccardGraph& jaccardGraph = *this;

    vector< vector<vertex_descriptor> > strongComponents;
    computeStronglyConnectedComponents(strongComponents);

    for(const vector<vertex_descriptor>& strongComponent: strongComponents) {
        for(const vertex_descriptor v: strongComponent) {
            vertexMap.erase(jaccardGraph[v].av);
            boost::clear_vertex(v, jaccardGraph);
            boost::remove_vertex(v, jaccardGraph);
        }
    }
}


bool JaccardGraph::markLongPathEdges(uint64_t minPathLength)
{
    JaccardGraph& jaccardGraph = *this;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        vertexIndexMap.insert({v, vertexIndex++});
    }

    // Topological sort.
    vector<vertex_descriptor> topologicallySortedVertices;
    try {
        boost::topological_sort(
            jaccardGraph,
            back_inserter(topologicallySortedVertices),
            boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));
    } catch (boost::not_a_dag&) {
        // Topological sort failed. Do nothing.
        return false;
    }
    std::reverse(topologicallySortedVertices.begin(), topologicallySortedVertices.end());



    // To compute the longest path, process vertices in topological order.
    // Set the length of the longest part ending at each vertex.
    for(const vertex_descriptor v0: topologicallySortedVertices) {
        uint64_t longestPathLength = 0;
        BGL_FORALL_INEDGES(v0, e, jaccardGraph, JaccardGraph) {
            const vertex_descriptor v1 = source(e, jaccardGraph);
            longestPathLength = max(longestPathLength,
                jaccardGraph[v1].longestPathLength + 1);
        }
        jaccardGraph[v0].longestPathLength = longestPathLength;
    }

    // Find the vertex with the longest path length.
    vertex_descriptor vLongest = null_vertex();
    uint64_t longestPathLength = 0;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        if(jaccardGraph[v].longestPathLength > longestPathLength) {
            vLongest = v;
            longestPathLength = jaccardGraph[v].longestPathLength;
        }
    }
    cout << "Longest path length is " << longestPathLength <<
        " at vertex " << vertexStringId(vLongest) << endl;




    // To compute the longest path, walk back starting here.
    vector<vertex_descriptor> longestPath;
    vertex_descriptor v0 = vLongest;
    while(true) {
        longestPath.push_back(v0);

        // Find the parent with the longest path.
        uint64_t longestParentLength = 0;
        vertex_descriptor vLongestParent = null_vertex();
        BGL_FORALL_INEDGES(v0, e, jaccardGraph, JaccardGraph) {
            const vertex_descriptor v1 = source(e, jaccardGraph);
            if(vLongestParent == null_vertex() or jaccardGraph[v1].longestPathLength > longestParentLength) {
                longestParentLength = jaccardGraph[v1].longestPathLength;
                vLongestParent = v1;
            }
        }
        if(vLongestParent == null_vertex()) {
            break;
        }
        v0 = vLongestParent;
    }
    std::reverse(longestPath.begin(), longestPath.end());
    SHASTA_ASSERT(longestPath.size() == longestPathLength + 1);

    cout << "Longest path:";
    for(const vertex_descriptor v: longestPath) {
        cout << " " << vertexStringId(v);
    }
    cout << endl;



    // Mark the edges in the longest path.
    for(uint64_t i=1; i<longestPath.size(); i++) {
        const vertex_descriptor v0 = longestPath[i-1];
        const vertex_descriptor v1 = longestPath[i];

        edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(v0, v1, jaccardGraph);
        SHASTA_ASSERT(edgeWasFound);

        jaccardGraph[e].isLongPathEdge = true;
    }

    // Also mark the vertices.
    for(const vertex_descriptor v: longestPath) {
        jaccardGraph[v].isLongPathVertex = true;
    }

    return true;
}



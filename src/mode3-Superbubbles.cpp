// Shasta.
#include "mode3-Superbubbles.hpp"
#include "mode3-AssemblyGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/dominator_tree.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include <queue>



Superbubbles::Superbubbles(
    AssemblyGraph& assemblyGraph,
    uint64_t maxOffset1     // Used to define superbubbles
    ) :
    assemblyGraph(assemblyGraph)
{
    assemblyGraph.numberVertices();
    const uint64_t vertexCount = num_vertices(assemblyGraph);

    vector<uint64_t> rank(vertexCount);
    vector<uint64_t> parent(vertexCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);

    // Compute connected components, using only edges with average offset up to maxOffset1.
    for(uint64_t i=0; i<vertexCount; i++) {
        disjointSets.make_set(i);
    }
    BGL_FORALL_EDGES(ce, assemblyGraph, AssemblyGraph) {
        uint64_t averageOffset;
        uint64_t minOffset;
        uint64_t maxOffset;
        assemblyGraph.bubbleChainOffset(assemblyGraph[ce], averageOffset, minOffset, maxOffset);
        if(averageOffset <= maxOffset1) {
            const vertex_descriptor cv0 = source(ce, assemblyGraph);
            const vertex_descriptor cv1 = target(ce, assemblyGraph);
            disjointSets.union_set(assemblyGraph[cv0].index, assemblyGraph[cv1].index);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<vertex_descriptor> > components(vertexCount);
    BGL_FORALL_VERTICES(cv, assemblyGraph, AssemblyGraph) {
        const uint64_t componentId = disjointSets.find_set(assemblyGraph[cv].index);
        components[componentId].push_back(cv);
    }

    // The superbubbles are the components with size at least 2.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<vertex_descriptor> component = components[componentId];
        if(components[componentId].size() > 1) {
            superbubbles.emplace_back(Superbubble(component));
        }
    }

    assemblyGraph.storeSuperbubblesInformation(*this);



    // Find entrances and exists of each superbubble.
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        Superbubble& superbubble = getSuperbubble(superbubbleId);

        // Find entrances. These are superbubble vertices with in-edges
        // from outside the superbubble or average offset up to maxOffset1.
        for(const vertex_descriptor cv0: superbubble) {
            BGL_FORALL_INEDGES(cv0, ce, assemblyGraph, AssemblyGraph) {
                uint64_t averageOffset;
                uint64_t minOffset;
                uint64_t maxOffset;
                assemblyGraph.bubbleChainOffset(assemblyGraph[ce], averageOffset, minOffset, maxOffset);
                const vertex_descriptor cv1 = source(ce, assemblyGraph);
                if((not isInSuperbubble(superbubbleId, cv1)) or (averageOffset > maxOffset1)) {
                    superbubble.entrances.push_back(cv0);
                    break;
                }
            }
        }

        // Find exits. These are superbubble vertices with out-edges
        // to outside the superbubble or average offset up to maxOffset1.
        vector<vertex_descriptor> exits;
        for(const vertex_descriptor cv0: superbubble) {
            BGL_FORALL_OUTEDGES(cv0, ce, assemblyGraph, AssemblyGraph) {
                uint64_t averageOffset;
                uint64_t minOffset;
                uint64_t maxOffset;
                assemblyGraph.bubbleChainOffset(assemblyGraph[ce], averageOffset, minOffset, maxOffset);
                const vertex_descriptor cv1 = target(ce, assemblyGraph);
                if((not isInSuperbubble(superbubbleId, cv1)) or (averageOffset > maxOffset1)) {
                    superbubble.exits.push_back(cv0);
                    break;
                }
            }
        }
     }

}



// This computes connected components using the set of edges
// for which the AssemblyGraphEdgePredicate returns true.
// Each connected component with more than one vertex becoens a Superbubble.
// This does not compute entrances and exits of each Superbubble.
Superbubbles::Superbubbles(
    AssemblyGraph& assemblyGraph,
    const AssemblyGraphEdgePredicate& edgePredicate) :
    assemblyGraph(assemblyGraph)
{
    assemblyGraph.numberVertices();
    const uint64_t vertexCount = num_vertices(assemblyGraph);

    vector<uint64_t> rank(vertexCount);
    vector<uint64_t> parent(vertexCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);

    // Compute connected components, using only edges
    // for which the AssemblyGraphEdgePredicate returns true.
    for(uint64_t i=0; i<vertexCount; i++) {
        disjointSets.make_set(i);
    }
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(edgePredicate(e)) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            disjointSets.union_set(assemblyGraph[v0].index, assemblyGraph[v1].index);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<vertex_descriptor> > components(vertexCount);
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const uint64_t componentId = disjointSets.find_set(assemblyGraph[v].index);
        components[componentId].push_back(v);
    }

    // The superbubbles are the components with size at least 2.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<vertex_descriptor> component = components[componentId];
        if(components[componentId].size() > 1) {
            superbubbles.emplace_back(Superbubble(component));
        }
    }

    assemblyGraph.storeSuperbubblesInformation(*this);
}



// This constructs superbubbles consisting of a single edge v0->v1
// such that:
// outDegree(v0) == 1, inDegree(v0)  > 1
// inDegree (v1) == 1, outDegree(v1) > 1
Superbubbles::Superbubbles(
    AssemblyGraph& assemblyGraph, const FromTangledEdges&) :
    assemblyGraph(assemblyGraph)
{

    // Try all edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());

        // Check outDegree of v0.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        if(out_degree(v0, assemblyGraph) != 1) {
            continue;
        }

        // Check inDegree of v1.
        const vertex_descriptor v1 = target(e, assemblyGraph);
        if(in_degree(v1, assemblyGraph) != 1) {
            continue;
        }

        // Check the number of entrances.
        const uint64_t entranceCount = in_degree(v0, assemblyGraph);
        if(entranceCount < 2) {
            continue;
        }

        // Check the number of exits.
        const uint64_t exitCount = out_degree(v1, assemblyGraph);
        if(exitCount < 2) {
            continue;
        }

        superbubbles.push_back(Superbubble());
        Superbubble& superbubble = superbubbles.back();
        superbubble.push_back(v0);
        superbubble.push_back(v1);
        superbubble.entrances.push_back(v0);
        superbubble.exits.push_back(v1);
    }

    assemblyGraph.storeSuperbubblesInformation(*this);
}



// This constructs superbubbles consisting of a single tangled vertex.
// A vertex v is tangled if:
// inDegree(v)>1, outDegree(v)>1.
Superbubbles::Superbubbles(
    AssemblyGraph& assemblyGraph, const FromTangledVertices&) :
    assemblyGraph(assemblyGraph)
{
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((in_degree(v, assemblyGraph)  > 1) and (out_degree(v, assemblyGraph) > 1)) {
            superbubbles.push_back(Superbubble());
            Superbubble& superbubble = superbubbles.back();
            superbubble.push_back(v);
            superbubble.entrances.push_back(v);
            superbubble.exits.push_back(v);
        }
    }

    assemblyGraph.storeSuperbubblesInformation(*this);
}



// This uses dominator trees.
// It only finds superbubbles with one entrance and one exit.
Superbubbles::Superbubbles(
    AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    const bool debug = false;
    if(debug) {
        cout << "Begin Superbubbles constructor using dominator trees." << endl;
    }

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> indexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        indexMap.insert({v, vertexIndex++});
    }
    auto associativeIndexMap = boost::make_assoc_property_map(indexMap);
    const uint64_t vertexCount = vertexIndex;

    // Vectors used below to compute the dominator tree.
    vector<uint64_t> dfNum(vertexCount);
    vector<vertex_descriptor> parent(vertexCount);
    vector<vertex_descriptor> verticesByDFNum(vertexCount);

    // Tree pairs found on forward and backward dominator tree.
    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;



    // Compute dominator trees using as entrance each of the
    // vertices with zero in-degree.
    BGL_FORALL_VERTICES(entrance, assemblyGraph, AssemblyGraph) {
        if(in_degree(entrance, assemblyGraph) != 0) {
            continue;
        }

        // Compute the dominator tree.
        fill(dfNum.begin(), dfNum.end(), invalid<uint64_t>);
        fill(parent.begin(), parent.end(), AssemblyGraph::null_vertex());
        fill(verticesByDFNum.begin(), verticesByDFNum.end(), AssemblyGraph::null_vertex());
        std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

        boost::lengauer_tarjan_dominator_tree(
            assemblyGraph,
            entrance,
            boost::make_assoc_property_map(indexMap),
            boost::make_iterator_property_map(dfNum.begin(), associativeIndexMap),
            boost::make_iterator_property_map(parent.begin(), associativeIndexMap),
            verticesByDFNum,
            boost::make_assoc_property_map(predecessorMap));

        if(debug) {
            cout << "Forward dominator tree with entrance at " <<
                anchorIdToString(assemblyGraph[entrance].getAnchorId()) << endl;
        }
        for(const auto& p: predecessorMap) {
            const vertex_descriptor cv0 = p.second;
            const vertex_descriptor cv1 = p.first;
            forwardPairs.push_back({cv0, cv1});
            if(debug) {
                cout << "F " << anchorIdToString(assemblyGraph[cv0].getAnchorId()) << "->" <<
                    anchorIdToString(assemblyGraph[cv1].getAnchorId()) << endl;
            }
        }
    }



    // Compute dominator trees on the reverse graph using as entrance each of the
    // vertices with zero in-degree on the reverse graph
    // (that is, zero out-degree on the AssemblyGraph).
    using ReverseAssemblyGraph = boost::reverse_graph<AssemblyGraph>;
    ReverseAssemblyGraph reverseGraph(assemblyGraph);
    BGL_FORALL_VERTICES(entrance, reverseGraph, ReverseAssemblyGraph) {
        if(in_degree(entrance, reverseGraph) != 0) {
            continue;
        }

        // Compute the dominator tree.
        fill(dfNum.begin(), dfNum.end(), invalid<uint64_t>);
        fill(parent.begin(), parent.end(), AssemblyGraph::null_vertex());
        fill(verticesByDFNum.begin(), verticesByDFNum.end(), AssemblyGraph::null_vertex());
        std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

        boost::lengauer_tarjan_dominator_tree(
            reverseGraph,
            entrance,
            boost::make_assoc_property_map(indexMap),
            boost::make_iterator_property_map(dfNum.begin(), associativeIndexMap),
            boost::make_iterator_property_map(parent.begin(), associativeIndexMap),
            verticesByDFNum,
            boost::make_assoc_property_map(predecessorMap));

        if(debug) {
            cout << "Backward dominator tree with exit at " <<
                anchorIdToString(assemblyGraph[entrance].getAnchorId()) << endl;
        }
        for(const auto& p: predecessorMap) {
            const vertex_descriptor cv0 = p.first;
            const vertex_descriptor cv1 = p.second;
            backwardPairs.push_back({cv0, cv1});
            if(debug) {
                cout << "B " << anchorIdToString(assemblyGraph[cv0].getAnchorId()) << "->" <<
                    anchorIdToString(assemblyGraph[cv1].getAnchorId()) << endl;
            }
        }
    }

    // Compute strongly connected components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        assemblyGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(indexMap)));

    // Gather the vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents(vertexCount);
    for(const auto& p: componentMap) {
        const vertex_descriptor v = p.first;
        const uint64_t componentId = p.second;
        SHASTA_ASSERT(componentId < vertexCount);
        strongComponents[componentId].push_back(v);
    }



    // The pairs that appear both in forwardPairs and backwardPairs define our superbubbles
    deduplicate(forwardPairs);
    deduplicate(backwardPairs);
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs)
        );

    if(debug) {
        cout << "Bidirectional pairs:" << endl;
        for(const auto& p: bidirectionalPairs) {
            const vertex_descriptor cv0 = p.first;
            const vertex_descriptor cv1 = p.second;
            cout << anchorIdToString(assemblyGraph[cv0].getAnchorId()) << "->" <<
                anchorIdToString(assemblyGraph[cv1].getAnchorId()) << endl;
        }
    }

    // Each bidirectional pair generates a superbubble if
    // the out-degree of the entrance and
    // the in-degree of the exit are greater than 1,
    // unless the entrance or exit or any of the
    // superbubble vertices are in a non-trivial strong component..
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor cv0 = p.first;
        const vertex_descriptor cv1 = p.second;
        if(out_degree(cv0, assemblyGraph) <= 1) {
            continue;
        }
        if(in_degree(cv1, assemblyGraph) <= 1) {
            continue;
        }
        if(strongComponents[componentMap[cv0]].size() > 1) {
            // The entrance is in a non-trivial strong component.
            continue;
        }
        if(strongComponents[componentMap[cv1]].size() > 1) {
            // The exit is in a non-trivial strong component.
            continue;
        }
        superbubbles.resize(superbubbles.size() + 1);
        Superbubble& superbubble = superbubbles.back();
        superbubble.entrances.push_back(cv0);
        superbubble.exits.push_back(cv1);
        superbubble.fillInFromEntranceAndExit(assemblyGraph);

        if(debug) {
            cout << "Tentative superbubble with entrance " << anchorIdToString(assemblyGraph[cv0].getAnchorId()) <<
                " exit " << anchorIdToString(assemblyGraph[cv1].getAnchorId()) << " and " << superbubble.size() <<
                " vertices total." << endl;
        }

        // If any vertices in the superbubble are in a non-trivial
        // strong component, remove it.
        for(const vertex_descriptor cv: superbubble) {
            if(strongComponents[componentMap[cv]].size() > 1) {
                superbubbles.pop_back();
                if(debug) {
                    cout << "This superbubble will not be stored because some vertices are in a non-trivial strong component." << endl;
                }
                break;
            }
        }
    }

    if(debug) {
        cout << "Superbubble entrance/exit pairs:" << endl;
        for(const Superbubble& superbubble: superbubbles) {
            const vertex_descriptor cv0 = superbubble.entrances.front();
            const vertex_descriptor cv1 = superbubble.exits.front();;
            cout << anchorIdToString(assemblyGraph[cv0].getAnchorId()) << "->" <<
                anchorIdToString(assemblyGraph[cv1].getAnchorId()) << endl;
        }
    }

    // Don't store superbubble informatiom in the AssemblyGraph, because the
    // superbubbles created in this way can have overlap.
    // assemblyGraph.storeSuperbubblesInformation(*this);

    if(debug) {
        cout << "End Superbubbles constructor using dominator trees." << endl;
    }

}



Superbubbles::Superbubbles(
    AssemblyGraph& assemblyGraph,
    const Empty&) :
    assemblyGraph(assemblyGraph)
{}



// Fill in the superbubble given a single entrance and exit.
void Superbubble::fillInFromEntranceAndExit(const AssemblyGraph& assemblyGraph)
{
    SHASTA_ASSERT(empty());
    SHASTA_ASSERT(entrances.size() == 1);
    SHASTA_ASSERT(exits.size() == 1);

    const vertex_descriptor entrance = entrances.front();
    const vertex_descriptor exit = exits.front();

    // Do a BFS starting at the entrance and stopping at the exit.
    std::set<vertex_descriptor> internalVertices;
    std::queue<vertex_descriptor> q;
    q.push(entrance);
    while(not q.empty()) {
        const vertex_descriptor cv0 = q.front();
        q.pop();
        BGL_FORALL_OUTEDGES(cv0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor cv1 = target(e, assemblyGraph);
            if(cv1 != exit) {
                if(not internalVertices.contains(cv1)) {
                    internalVertices.insert(cv1);
                    q.push(cv1);
                }
            }
        }
    }

    push_back(entrance);
    copy(internalVertices.begin(), internalVertices.end(), back_inserter(*this));
    push_back(exit);

}



Superbubbles::~Superbubbles()
{
    assemblyGraph.clearVertexNumbering();
}




// Figure out if a vertex is in the specified superbubble.
bool Superbubbles::isInSuperbubble(uint64_t superbubbleId, vertex_descriptor cv) const
{
    return assemblyGraph[cv].superbubbleId == superbubbleId;
}


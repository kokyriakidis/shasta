// Shasta.
#include "mode3b-PathGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "MarkerGraph.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>



GlobalPathGraph::GlobalPathGraph(const Assembler& assembler) :
    assembler(assembler)
{
}



// Write each connected component in graphviz format.
void GlobalPathGraph::writeComponentsGraphviz(
    const string& baseName,
    const GlobalPathGraphDisplayOptions& options) const
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        const PathGraph& component = *components[componentRank];
        component.writeGraphviz(
            baseName + "_Component_" + to_string(componentRank),
            options);
    }
}



void GlobalPathGraph::createVertices()
{
    const MarkerGraph& markerGraph = assembler.markerGraph;

    verticesVector.clear();

    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        if(markerGraph.edges[edgeId].isPrimary == 1) {
            verticesVector.push_back(GlobalPathGraphVertex(edgeId));
        }
    }
}



// Return the vertexId corresponding to a given MarkerGraphEdgeId, or
// invalid<MarkerGraphEdgeId> if no such a vertex exists.
uint64_t GlobalPathGraph::getVertexId(MarkerGraphEdgeId edgeId) const
{
    GlobalPathGraphVertex targetVertex(edgeId);
    auto it = std::lower_bound(verticesVector.begin(), verticesVector.end(), targetVertex);

    if((it == verticesVector.end()) or (it->edgeId != edgeId)) {

        // Not found.
        return invalid<uint64_t>;

    } else {

        // Found it. Return its vertexId.
        return it - verticesVector.begin();

    }
}



// Use the journeys stored in the MarkerGraph to create edges.
// This can be converted to multithreaded code and use mapped memory.
void GlobalPathGraph::createEdges()
{
    // Gather pairs of consecutive edges in oriented read primary journeys.
    vector< pair<MarkerGraphEdgeId, MarkerGraphEdgeId> > edgePairs;

    // Loop over all primary journeys.
    const auto& primaryJourneys = assembler.markerGraph.primaryJourneys;
    for(uint64_t i=0; i<primaryJourneys.size(); i++) {
        const auto journey = primaryJourneys[i];

        // Loop over pairs of consecutive MarkerGraphEdgeIds in this journey.
        for(uint64_t position1=1; position1<journey.size(); position1++) {
            const uint64_t position0 = position1 - 1;
            edgePairs.push_back({journey[position0].edgeId, journey[position1].edgeId});
        }
    }

    // Deduplicate the edge pairs and count the number of times each of them was found.
    vector<uint64_t> coverage;
    deduplicateAndCount(edgePairs, coverage);
    SHASTA_ASSERT(edgePairs.size() == coverage.size());
    edgePairs.shrink_to_fit();
    coverage.shrink_to_fit();

    // To generate edges, we need a table to map MarkerGraphEdgeIds
    // to vertex ids in the GlobalPathGraph.
    vector<uint64_t> primaryEdgeTable(assembler.markerGraph.edges.size(), invalid<uint64_t>);
    for(uint64_t vertexId=0; vertexId<verticesVector.size(); vertexId++) {
        const GlobalPathGraphVertex& vertex = verticesVector[vertexId];
        primaryEdgeTable[vertex.edgeId] = vertexId;
    }



    // Now we can generate the edges.
    edges.clear();
    for(uint64_t i=0; i<edgePairs.size(); i++) {
        const uint64_t c = coverage[i];
        const auto& p = edgePairs[i];

        // Get the GlobalPathGraph vertex ids corresponding to
        // the MarkerGraphedgeIds in this pair.
        const MarkerGraphEdgeId edgeId0 = p.first;
        const MarkerGraphEdgeId edgeId1 = p.second;
        const uint64_t vertexId0 = primaryEdgeTable[edgeId0];
        const uint64_t vertexId1 = primaryEdgeTable[edgeId1];
        SHASTA_ASSERT(vertexId0 != invalid<uint64_t>);
        SHASTA_ASSERT(vertexId1 != invalid<uint64_t>);

        // Create the edge.
        GlobalPathGraphEdge edge;
        edge.vertexId0 = vertexId0;
        edge.vertexId1 = vertexId1;
        edge.coverage = c;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, edge.info));
        edges.push_back(edge);
    }
}



// Write the entire GlobalPathGraph in graphviz format.
void GlobalPathGraph::writeGraphviz() const
{
    ofstream out("GlobalPathGraph.dot");
    out << "digraph GlobalPathGraph {\n";

    for(const GlobalPathGraphEdge& edge: edges) {
        const MarkerGraphEdgeId edgeId0 = verticesVector[edge.vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[edge.vertexId1].edgeId;
        out << edgeId0 << "->";
        out << edgeId1 << ";\n";
    }

    out << "}\n";
}



// Create connected components.
// This only considers edges with corrected Jaccard at least equal to
// minCorrectedJaccard, and only stores connected components with at
// least minComponentSize vertices.
void GlobalPathGraph::createComponents(
    double minCorrectedJaccard,
    uint64_t minComponentSize)
{
    // Compute connected components.
    const uint64_t n = verticesVector.size();
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    for(const GlobalPathGraphEdge& edge: edges) {
        if(edge.info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        disjointSets.union_set(edge.vertexId0, edge.vertexId1);
    }

    // Generate vertices of each connected component.
    vector< shared_ptr<PathGraph> > allComponents;
    for(uint64_t componentId=0; componentId<n; componentId++) {
        allComponents.push_back(make_shared<PathGraph>());
    }
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        const GlobalPathGraphVertex& vertex = verticesVector[vertexId];
        const uint64_t componentId = disjointSets.find_set(vertexId);
        allComponents[componentId]->addVertex(vertexId, vertex.edgeId);
    }

    // Create edges of each connected component.
    for(const GlobalPathGraphEdge& edge: edges) {
        if(edge.info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        const uint64_t vertexId0 = edge.vertexId0;
        const uint64_t vertexId1 = edge.vertexId1;
        const uint64_t componentId = disjointSets.find_set(vertexId0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(vertexId1));
        const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1].edgeId;
        allComponents[componentId]->addEdge(edgeId0, edgeId1, edge.info, edge.coverage);
    }

    // Only keep the connected components that have at least
    // minComponentSize vertices.
    components.clear();
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const shared_ptr<PathGraph> componentPointer = allComponents[componentId];
        const PathGraph& component = *componentPointer;
        if(num_vertices(component) >= minComponentSize) {
            components.push_back(componentPointer);
        }
    }



    // Sort the pointers to connected components by decreasing size.
    class ComponentSorter {
    public:
        bool operator()(
            const shared_ptr<PathGraph>& x,
            const shared_ptr<PathGraph>& y)
        {
            return num_vertices(*x) > num_vertices(*y);
        }
    };
    sort(components.begin(), components.end(), ComponentSorter());
}



// Local transitive reduction of each connected component.
void GlobalPathGraph::localTransitiveReduction(
    uint64_t distance,
    uint64_t maxCoverage)
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph& component = *components[componentRank];
        component.localTransitiveReduction(distance, maxCoverage);
    }
}



void PathGraph::localTransitiveReduction(
    uint64_t distance,
    uint64_t maxCoverage)
{
    PathGraph& graph = *this;

    // Histograms by coverage.
    // (all edges, only edges removed by transitive reduction.
    vector< pair<uint64_t, uint64_t> > histogram;

    // We want to process edges in order of increasing coverage.
    vector<pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, PathGraph) {
        edgeTable.push_back({e, graph[e].coverage});
    }
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<edge_descriptor, uint64_t>());

    // Loop over all edges v0->v1 in order of increasing coverage.
    for(const auto& p: edgeTable) {
        edge_descriptor e01 = p.first;
        const vertex_descriptor v0 = source(e01, graph);
        const vertex_descriptor v1 = target(e01, graph);

        const uint64_t coverage = graph[e01].coverage;
        if(histogram.size() <= coverage) {
            histogram.resize(coverage + 1, {0, 0});
        }
        ++histogram[coverage].first;

        if(coverage > maxCoverage) {
            continue;
        }

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
            BGL_FORALL_OUTEDGES_T(vA, eAB, graph, PathGraph) {

                // Dont's use e01 in the BFS.
                if(eAB == e01) {
                    continue;
                }

                // If eABwas already flagged as removed during transitive reduction,
                // don't use it.
                if(graph[eAB].isNonTransitiveReductionEdge) {
                    continue;
                }

                // If we reached v1, mark e01 as a nonTransitiveReduction edge
                // and stop the BFS.
                const vertex_descriptor vB = target(eAB, graph);
                if(vB == v1) {
                    graph[e01].isNonTransitiveReductionEdge = true;
                    endBfs = true;
                    // cout << "Reached " << v1 << endl;
                    ++histogram[coverage].second;
                    break;
                }

                // If we already reached this vertex, do nothing.
                if(m.contains(vB)) {
                    continue;
                }

                // If not at maximum distance, enqueue vB.
                if(distanceB < distance) {
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

    if(false) {
        cout << "PathGraph edge coverage histogram" << endl;
        cout << "Coverage,All edges,Edges removed by transitive reduction" << endl;
        for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
            const auto& p = histogram[coverage];
            if(p.first>0 or p.second>0) {
                cout << coverage << "," << p.first << "," << p.second << endl;
            }
        }
    }


}



void PathGraph::addVertex(
    uint64_t vertexId,
    MarkerGraphEdgeId edgeId)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    const vertex_descriptor v = add_vertex({vertexId, edgeId}, *this);
    vertexMap.insert({edgeId, v});
}



void PathGraph::addEdge(
    MarkerGraphEdgeId edgeId0,
    MarkerGraphEdgeId edgeId1,
    const MarkerGraphEdgePairInfo& info,
    uint64_t coverage)
{
    auto it0 = vertexMap.find(edgeId0);
    auto it1 = vertexMap.find(edgeId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    add_edge(v0, v1, {info, coverage}, *this);
}



// Write a PathGraph in graphviz format.
void PathGraph::writeGraphviz(
    const string& name,
    const GlobalPathGraphDisplayOptions& options) const
{
    ofstream out(name + ".dot");

    const PathGraph& graph = *this;
    out << "digraph " << name << " {\n";

    BGL_FORALL_VERTICES(v, graph, PathGraph) {
        const PathGraphVertex& vertex = graph[v];
        out << vertex.edgeId;

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "[";
        }

        if(options.labels) {
            out << "label=\"";
            out << vertex.edgeId /* << "\\n" << markerGraph.edgeCoverage(vertex.edgeId)*/;
            out << "\" ";
        }

        if(options.tooltips) {
            out << "tooltip=\"";
            out << vertex.edgeId;
            out << "\" ";
        }

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "]";
        }
        out << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, PathGraph) {
        const PathGraphEdge& edge = graph[e];
        if(not options.showNonTransitiveReductionEdges and edge.isNonTransitiveReductionEdge) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        out <<
            graph[v0].edgeId << "->" <<
            graph[v1].edgeId;

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << " [";
        }

        if(edge.isNonTransitiveReductionEdge) {
            out << "style=dashed ";
        }

        if(options.tooltips) {
            out <<
                "tooltip=\"" <<
                graph[v0].edgeId << "->" <<
                graph[v1].edgeId << " ";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << " " <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << " " <<
                edge.info.offsetInBases << "\" ";
        }

        if(options.labels) {
            out <<
                "label=\"";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << "\\n" <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << "\\n" <<
                edge.info.offsetInBases << "\" ";

        }

        // Color.
        if(options.colorEdges) {
            const double correctedJaccard = edge.info.correctedJaccard();
            if(correctedJaccard <= options.redJ) {
                out << " color=red ";
            } else if(correctedJaccard >= options.greenJ) {
                out << " color=green ";
            } else {
                const double hue = (correctedJaccard - options.redJ) / (3. * (options.greenJ - options.redJ));
                out << " color=\"" << hue << ",1,1\" ";
            }
        }

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << "]";
        }
        out << ";\n";
    }

    out << "}\n";
}



// Create the connected components of this PathGraph,
// without changing the PathGraph itself.
vector< shared_ptr<PathGraph> > PathGraph::createConnectedComponents(
    uint64_t minComponentSize) const
{
    const PathGraph& graph = *this;

    // Compute connected components.
    // We can't use boost::connected_components because it only works
    // for undirected graphs.
    const uint64_t n = num_vertices(graph);
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    BGL_FORALL_EDGES(e, graph, PathGraph) {
        const PathGraph::vertex_descriptor v0 = source(e, graph);
        const PathGraph::vertex_descriptor v1 = target(e, graph);
        disjointSets.union_set(v0, v1);
    }


    // Gather the vertices in each connected component.
    vector< shared_ptr<PathGraph> > allComponentPointers(num_vertices(graph));
    BGL_FORALL_VERTICES(v, graph, PathGraph) {
        const PathGraphVertex& vertex = graph[v];
        const uint64_t componentId = disjointSets.find_set(v);
        shared_ptr<PathGraph>& componentPointer = allComponentPointers[componentId];
        if(not componentPointer) {
            componentPointer = make_shared<PathGraph>();
        }
        PathGraph& component = *componentPointer;
        component.addVertex(
            vertex.vertexId,
            vertex.edgeId);
    }


    // Gather the edges in each connected component.
    BGL_FORALL_EDGES(e, graph, PathGraph) {
        const PathGraph::vertex_descriptor v0 = source(e, graph);
        const PathGraph::vertex_descriptor v1 = target(e, graph);
        const uint64_t edgeId0 = graph[v0].edgeId;
        const uint64_t edgeId1 = graph[v1].edgeId;
        const uint64_t componentId = disjointSets.find_set(v0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(v1));
        shared_ptr<PathGraph>& componentPointer = allComponentPointers[componentId];
        SHASTA_ASSERT(componentPointer);
        PathGraph& component = *componentPointer;
        component.addEdge(
            edgeId0,
            edgeId1,
            graph[e].info,
            graph[e].coverage);
    }



    // Keep only the components with at least minComponentSize vertices
    // and sort them by size.
    vector< pair<shared_ptr<PathGraph>, uint64_t> > componentPointersWithSizes;
    for(const shared_ptr<PathGraph>& p: allComponentPointers) {
        if(p) {
            const uint64_t componentSize = num_vertices(*p);
            if(componentSize >= minComponentSize) {
                componentPointersWithSizes.push_back({p, componentSize});
            }
        }
    }
    sort(componentPointersWithSizes.begin(), componentPointersWithSizes.end(),
        OrderPairsBySecondOnlyGreater<shared_ptr<PathGraph>, uint64_t>());


    // For now return all components, including the empty ones.
    // But we want to remove the small ones and sort them by size.
    vector< shared_ptr<PathGraph> > componentPointers;
    for(const auto& p: componentPointersWithSizes) {
        componentPointers.push_back(p.first);
    }
    return componentPointers;
}



// Remove cross-edges.
// This removes an edge v0->v1 if the following are all true:
// - It is not marked as removed by transitive reduction.
// - Its coverage is at most lowCoverageThreshold.
// - Its estimated offset is at least minOffset.
// - v0 has at least one out-edge with coverage at least highCoverageThreshold
//   (ignoring edges marked as removed by transitive reduction).
// - v1 has at least one in-edge with coverage at least highCoverageThreshold.
//   (ignoring edges marked as removed by transitive reduction).
void PathGraph::removeCrossEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    uint64_t minOffset)
{
    PathGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PathGraph) {
        const PathGraphEdge& edge = graph[e];

        // If it is marked as removed by transitive reduction, skip it.
        if(edge.isNonTransitiveReductionEdge) {
            continue;
        }

        // Check coverage.
        if(edge.coverage > lowCoverageThreshold) {
            continue;
        }

        // Check estimated offset.
        if(edge.info.offsetInBases < int64_t(minOffset)) {
            continue;
        }

        // Check out-edges of v0.
        const vertex_descriptor v0 = source(e, graph);
        bool v0HasStrongOutEdge = false;
        BGL_FORALL_OUTEDGES(v0, e0, graph, PathGraph) {
            // If it is marked as removed by transitive reduction, ignore it.
            if(graph[e0].isNonTransitiveReductionEdge) {
                continue;
            }
            if(graph[e0].coverage >= highCoverageThreshold) {
                v0HasStrongOutEdge = true;
                break;
            }
        }
        if(not v0HasStrongOutEdge) {
            continue;
        }

        // Check in-edges of v1.
        const vertex_descriptor v1 = target(e, graph);
        bool v1HasStrongOutEdge = false;
        BGL_FORALL_INEDGES(v1, e1, graph, PathGraph) {
            // If it is marked as removed by transitive reduction, ignore it.
            if(graph[e1].isNonTransitiveReductionEdge) {
                continue;
            }
            if(graph[e1].coverage >= highCoverageThreshold) {
                v1HasStrongOutEdge = true;
                break;
            }
        }
        if(not v1HasStrongOutEdge) {
            continue;
        }

        // If all above checks passed, this edge will be removed.
        edgesToBeRemoved.push_back(e);
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}

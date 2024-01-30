// Shasta.
#include "mode3b-PathGraph1.hpp"
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



void GlobalPathGraph1::assemble(
    const Assembler& assembler,
    uint64_t threadCount0,
    uint64_t threadCount1)
{
    assemble2(assembler, threadCount0, threadCount1);
}



GlobalPathGraph1::GlobalPathGraph1(const Assembler& assembler) :
    assembler(assembler)
{
}



// Write each connected component in graphviz format.
void GlobalPathGraph1::writeComponentsGraphviz(
    const string& baseName,
    const GlobalPathGraph1DisplayOptions& options) const
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        const PathGraph1& component = *components[componentRank];
        component.writeGraphviz(
            verticesVector,
            baseName + "_Component_" + to_string(componentRank),
            options);
    }
}



// Find out if a marker graph edge is a primary edge.
bool GlobalPathGraph1::isPrimary(
    MarkerGraphEdgeId edgeId,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage) const
{
    // Check coverage.
    const MarkerGraph& markerGraph = assembler.markerGraph;
    const uint64_t coverage = markerGraph.edgeCoverage(edgeId);
    if(coverage < minPrimaryCoverage) {
        return false;
    }
    if(coverage > maxPrimaryCoverage) {
        return false;
    }

    // Check for duplicate oriented reads on the edge.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeId)) {
        return false;
    }

    // Check for duplicate oriented reads on its vertices.
    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
    if(
        const auto& markers = assembler.markers;
        markerGraph.vertexHasDuplicateOrientedReadIds(edge.source, markers) or
        markerGraph.vertexHasDuplicateOrientedReadIds(edge.target, markers)) {
        return false;
    }

#if 0
    // Check that is also is a branch edge.
    // Requiring this is too strict and removes opportunities for
    // finding paths.
    if(not isBranchEdge(edgeId, minPrimaryCoverage)) {
        return false;
    }
#endif

    // If all above checks passed, this is a primary edge.
    return true;
}



// Find out if a marker graph edge is a branch edge.
// A marker graph edge is a branch edge if:
// - Its source vertex has more than one outgoing edge with coverage at least minPrimaryCoverage.
// OR
// - Its target vertex has more than one incoming edge with coverage at least minPrimaryCoverage.
bool GlobalPathGraph1::isBranchEdge(
    MarkerGraphEdgeId edgeId,
    uint64_t minEdgeCoverage) const
{
    // Access this marker graph edge and its vertices.
    const MarkerGraph::Edge& edge = assembler.markerGraph.edges[edgeId];
    const MarkerGraphVertexId vertexId0 = edge.source;
    const MarkerGraphVertexId vertexId1 = edge.target;

    // Check outgoing edges of vertexId0.
    const auto outgoingEdges0 = assembler.markerGraph.edgesBySource[vertexId0];
    uint64_t count0 = 0;
    for(const MarkerGraphEdgeId edgeId0: outgoingEdges0) {
        if(assembler.markerGraph.edgeCoverage(edgeId0) >= minEdgeCoverage) {
            ++count0;
        }
    }
    if(count0 > 1) {
        return true;
    }

    // Check incoming edges of vertexId1.
    const auto incomingEdges1 = assembler.markerGraph.edgesByTarget[vertexId1];
    uint64_t count1 = 0;
    for(const MarkerGraphEdgeId edgeId1: incomingEdges1) {
        if(assembler.markerGraph.edgeCoverage(edgeId1) >= minEdgeCoverage) {
            ++count1;
        }
    }
    if(count1 > 1) {
        return true;
    }

    return false;
}



void GlobalPathGraph1::createVertices(
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage)
{
    const MarkerGraph& markerGraph = assembler.markerGraph;

    verticesVector.clear();

    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        if(isPrimary(edgeId, minPrimaryCoverage, maxPrimaryCoverage)) {
            verticesVector.push_back(GlobalPathGraph1Vertex(edgeId));
        }
    }
}



// Return the vertexId corresponding to a given MarkerGraphEdgeId, or
// invalid<MarkerGraphEdgeId> if no such a vertex exists.
uint64_t GlobalPathGraph1::getVertexId(MarkerGraphEdgeId edgeId) const
{
    GlobalPathGraph1Vertex targetVertex(edgeId);
    auto it = std::lower_bound(verticesVector.begin(), verticesVector.end(), targetVertex);

    if((it == verticesVector.end()) or (it->edgeId != edgeId)) {

        // Not found.
        return invalid<uint64_t>;

    } else {

        // Found it. Return its vertexId.
        return it - verticesVector.begin();

    }
}




// The "journey" of each oriented read is the sequence of vertices it encounters.
// It stores pairs (ordinal0, vertexId) for each oriented read, sorted by ordinal0.
// The vertexId is the index in verticesVector.
// Indexed by OrientedReadId::getValue.
// Journeys are used to generate edges by "following the reads".
void GlobalPathGraph1::computeOrientedReadJourneys()
{
    orientedReadJourneys.clear();
    orientedReadJourneys.resize(assembler.markers.size());

    for(uint64_t vertexId=0; vertexId<verticesVector.size(); vertexId++) {
        const MarkerGraphEdgeId edgeId = verticesVector[vertexId].edgeId;

        // Loop over MarkerIntervals of this primary marker graph edge.
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            orientedReadJourneys[orientedReadId.getValue()].push_back(make_pair(ordinal0, vertexId));
        }
    }

    // Now sort them, for each oriented read.
    for(auto& orientedReadJourney: orientedReadJourneys) {
        sort(orientedReadJourney.begin(), orientedReadJourney.end(),
            OrderPairsByFirstOnly<uint32_t, uint64_t>());
    }



    // Store journey information in the vertices.
    for(ReadId readId=0; readId<assembler.markers.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = orientedReadJourneys[orientedReadId.getValue()];

            for(uint64_t position=0; position<journey.size(); position++) {
                const auto& p = journey[position];
                const uint64_t vertexId = p.second;
                verticesVector[vertexId].journeyInfoItems.push_back({orientedReadId, position});
            }
        }
    }



    // Write the journeys to csv.
    ofstream csv("GlobalPathGraphJourneys.csv");
    for(ReadId readId=0; readId<assembler.markers.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = orientedReadJourneys[orientedReadId.getValue()];
            csv << orientedReadId << ",";
            for(const auto& p: journey) {
                const uint64_t vertexId = p.second;
                const MarkerGraphEdgeId edgeId = verticesVector[vertexId].edgeId;
                csv << edgeId << ",";
            }
            csv << "\n";
        }
    }
}



void GlobalPathGraph1::createEdges0(
    uint64_t maxDistanceInJourney,
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard)
{
    // Candidate edges are pairs of vertices that appear near each other
    // in oriented read journeys.
    vector< pair<uint64_t, uint64_t> > candidateEdges;
    for(uint64_t i=0; i<orientedReadJourneys.size(); i++) {
        const auto& journey = orientedReadJourneys[i];

        for(uint64_t position0=0; position0 < journey.size(); position0++) {
            for(uint64_t distance = 1; distance <= maxDistanceInJourney; distance++) {
                const uint64_t position1 = position0 + distance;
                if(position1 >= journey.size()) {
                    break;
                }
                candidateEdges.push_back({journey[position0].second, journey[position1].second});
            }
        }
    }
    // cout << timestamp << "Found " << candidateEdges.size() << " candidate edges, including duplicates." << endl;

    // Deduplicate the candidate edges and count the number of times
    // each of them was found. Keep only the ones that occurred at least
    // minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateEdges, coverage, minEdgeCoverage);
    // cout << timestamp << "After deduplication, there are " << candidateEdges.size() << " candidate edges." << endl;
    SHASTA_ASSERT(candidateEdges.size() == coverage.size());
    candidateEdges.shrink_to_fit();
    coverage.shrink_to_fit();

    // For each candidate edge, compute correctedJaccard, and if high enough
    // generate an edge.
    edges.clear();
    for(uint64_t i=0; i<candidateEdges.size(); i++) {
        const uint64_t c = coverage[i];
        SHASTA_ASSERT(c >= minEdgeCoverage);
        const auto& p = candidateEdges[i];
        const uint64_t vertexId0 = p.first;
        const uint64_t vertexId1 = p.second;
        GlobalPathGraph1Edge edge;
        const MarkerGraphEdgeId edgeId0 = verticesVector[vertexId0].edgeId;
        const MarkerGraphEdgeId edgeId1 = verticesVector[vertexId1].edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, edge.info));
        if(edge.info.correctedJaccard() >= minCorrectedJaccard) {
            edge.vertexId0 = vertexId0;
            edge.vertexId1 = vertexId1;
            edge.coverage = c;
            edges.push_back(edge);
        }
    }

}



// Find children edges of vertexId0.
// The first element of each pair of the children vector
// is the vertexId of the child vertex.
void GlobalPathGraph1::findChildren(
    uint64_t vertexId0,
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard,
    vector< pair<uint64_t, MarkerGraphEdgePairInfo> >& children)
{
    const GlobalPathGraph1Vertex& vertex0 = verticesVector[vertexId0];
    const MarkerGraphEdgeId edgeId0 = vertex0.edgeId;

    // Find vertices encountered later on journeys that go through here.
    vector<uint64_t> candidateChildren;
    for(const auto& journeyInfoItem: vertex0.journeyInfoItems) {
        const OrientedReadId orientedReadId = journeyInfoItem.orientedReadId;
        const auto& journey = orientedReadJourneys[orientedReadId.getValue()];
        for(uint64_t position = journeyInfoItem.positionInJourney+1;
            position < journey.size(); position++) {
            candidateChildren.push_back(journey[position].second);
        }
    }

    // Count the candidate children and only keep the ones
    // that occurred at least minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateChildren, coverage, minEdgeCoverage);

    // Store the ones that have sufficient correctedJaccard.
    children.clear();
    MarkerGraphEdgePairInfo info;
    for(const uint64_t vertexId1: candidateChildren) {
        const GlobalPathGraph1Vertex& vertex1 = verticesVector[vertexId1];
        const MarkerGraphEdgeId edgeId1 = vertex1.edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
        if(info.correctedJaccard() >= minCorrectedJaccard) {
            children.push_back({vertexId1, info});
        }
    }
}



void GlobalPathGraph1::findParents(
    uint64_t vertexId0,
    uint64_t minEdgeCoverage,
    double minCorrectedJaccard,
    vector< pair<uint64_t, MarkerGraphEdgePairInfo> >& parents)
{
    const GlobalPathGraph1Vertex& vertex0 = verticesVector[vertexId0];
    const MarkerGraphEdgeId edgeId0 = vertex0.edgeId;

    // Find vertices encountered earlier on journeys that go through here.
    vector<uint64_t> candidateParents;
    for(const auto& journeyInfoItem: vertex0.journeyInfoItems) {
        const OrientedReadId orientedReadId = journeyInfoItem.orientedReadId;
        const auto& journey = orientedReadJourneys[orientedReadId.getValue()];
        for(int64_t position = int64_t(journeyInfoItem.positionInJourney)-1;
            position>=0; position--) {
            candidateParents.push_back(journey[position].second);
        }
    }

    // Count the candidate parents and only keep the ones
    // that occurred at least minEdgeCoverage times.
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(candidateParents, coverage, minEdgeCoverage);

    // Store the ones that have sufficient correctedJaccard.
    parents.clear();
    MarkerGraphEdgePairInfo info;
    for(const uint64_t vertexId1: candidateParents) {
        const GlobalPathGraph1Vertex& vertex1 = verticesVector[vertexId1];
        const MarkerGraphEdgeId edgeId1 = vertex1.edgeId;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId1, edgeId0, info));
        if(info.correctedJaccard() >= minCorrectedJaccard) {
            parents.push_back({vertexId1, info});
        }
    }
}



// Write the entire GlobalPathGraph in graphviz format.
void GlobalPathGraph1::writeGraphviz() const
{
    ofstream out("GlobalPathGraph.dot");
    out << "digraph GlobalPathGraph {\n";

    for(const GlobalPathGraph1Edge& edge: edges) {
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
void GlobalPathGraph1::createComponents(
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
    for(const GlobalPathGraph1Edge& edge: edges) {
        if(edge.info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        disjointSets.union_set(edge.vertexId0, edge.vertexId1);
    }

    // Generate vertices of each connected component.
    vector< shared_ptr<PathGraph1> > allComponents;
    for(uint64_t componentId=0; componentId<n; componentId++) {
        allComponents.push_back(make_shared<PathGraph1>());
    }
    for(uint64_t vertexId=0; vertexId<n; vertexId++) {
        const GlobalPathGraph1Vertex& vertex = verticesVector[vertexId];
        const uint64_t componentId = disjointSets.find_set(vertexId);
        allComponents[componentId]->addVertex(vertexId, vertex.edgeId);
    }

    // Create edges of each connected component.
    for(const GlobalPathGraph1Edge& edge: edges) {
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
        const shared_ptr<PathGraph1> componentPointer = allComponents[componentId];
        const PathGraph1& component = *componentPointer;
        if(num_vertices(component) >= minComponentSize) {
            components.push_back(componentPointer);
        }
    }



    // Sort the pointers to connected components by decreasing size.
    class ComponentSorter {
    public:
        bool operator()(
            const shared_ptr<PathGraph1>& x,
            const shared_ptr<PathGraph1>& y)
        {
            return num_vertices(*x) > num_vertices(*y);
        }
    };
    sort(components.begin(), components.end(), ComponentSorter());
}



// K-nn of each connected component:
// for each vertex, keep only the k best outgoing and incoming edges,
// as measured by correctedJaccard of each edge.
// This can break contiguity of the connected component.
void GlobalPathGraph1::knn(uint64_t k)
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        component.knn(k);
    }
}



// Local transitive reduction of each connected component.
void GlobalPathGraph1::localTransitiveReduction(
    uint64_t distance,
    uint64_t maxCoverage)
{
    for(uint64_t componentRank=0; componentRank<components.size(); componentRank++) {
        PathGraph1& component = *components[componentRank];
        component.localTransitiveReduction(distance, maxCoverage);
    }
}



void PathGraph1::localTransitiveReduction(
    uint64_t distance,
    uint64_t maxCoverage)
{
    PathGraph1& graph = *this;

    // Histograms by coverage.
    // (all edges, only edges removed by transitive reduction.
    vector< pair<uint64_t, uint64_t> > histogram;

    // We want to process edges in order of increasing coverage.
    vector<pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
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
            BGL_FORALL_OUTEDGES_T(vA, eAB, graph, PathGraph1) {

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
        cout << "PathGraph1 edge coverage histogram" << endl;
        cout << "Coverage,All edges,Edges removed by transitive reduction" << endl;
        for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
            const auto& p = histogram[coverage];
            if(p.first>0 or p.second>0) {
                cout << coverage << "," << p.first << "," << p.second << endl;
            }
        }
    }


}



void PathGraph1::addVertex(
    uint64_t vertexId,
    MarkerGraphEdgeId edgeId)
{
    SHASTA_ASSERT(not vertexMap.contains(edgeId));
    const vertex_descriptor v = add_vertex({vertexId, edgeId}, *this);
    vertexMap.insert({edgeId, v});
}



void PathGraph1::addEdge(
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



// Write a PathGraph1 in graphviz format.
void PathGraph1::writeGraphviz(
    const vector<GlobalPathGraph1Vertex>& globalVertices,
    const string& name,
    const GlobalPathGraph1DisplayOptions& options) const
{
    ofstream out(name + ".dot");

    const PathGraph1& graph = *this;
    out << "digraph " << name << " {\n";

    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        const PathGraph1Vertex& vertex = graph[v];
        const GlobalPathGraph1Vertex& globalVertex = globalVertices[vertex.vertexId];
        out << vertex.edgeId;

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "[";
        }

        if(options.labels) {
            out << "label=\"";
            out << vertex.edgeId << "\\n" << globalVertex.journeyInfoItems.size();
            if(globalVertex.chainId != invalid<uint64_t>) {
                out << "\\n" << globalVertex.chainId << ":" << globalVertex.positionInChain;
                if(globalVertex.isFirstInChain) {
                    out << "\\nFirst in chain";
                }
                if(globalVertex.isLastInChain) {
                    out << "\\nLast in chain";
                }
            }
            out << "\" ";
        }

        if(options.tooltips) {
            out << "tooltip=\"";
            out << vertex.edgeId;
            if(globalVertex.chainId != invalid<uint64_t>) {
                out << " " << globalVertex.chainId << ":" << globalVertex.positionInChain;
                if(globalVertex.isFirstInChain) {
                    out << " first in chain";
                }
                if(globalVertex.isLastInChain) {
                    out << " last in chain";
                }
            }
            out << "\" ";
        }

        // If it belongs to a chain, color it.
        if(options.colorVertices) {
            if(globalVertex.chainId != invalid<uint64_t>) {
                if(globalVertex.isFirstInChain) {
                    out << " style=filled fillcolor=blue ";
                } else if(globalVertex.isLastInChain) {
                    out << " style=filled fillcolor=orange ";
                } else {
                    // const uint32_t hue = MurmurHash2(&vertex.chainId, sizeof(vertex.chainId), 231) % 100;
                    // out << " color=\"" << 0.01 * double(hue) << ",0.4,1\"";
                    out << " style=filled fillcolor=cyan ";
                }
            }
        }

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "]";
        }
        out << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1Edge& edge = graph[e];
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



// For each vertex, only keep the best k outgoing and k incoming edges.
// "Best" as defined by correctedJaccard of the edges.
void PathGraph1::knn(uint64_t k)
{
    PathGraph1& graph = *this;

    // Store here the edges we want to keep.
    std::set<edge_descriptor> edgesToBeKept;

    // For each vertex, mark as to be kept the best k outgoing
    // and incoming edges.
    vector< pair<edge_descriptor, double> > adjacentEdges;  // With correctedJaccard.
    BGL_FORALL_VERTICES(v, graph, PathGraph1) {

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            adjacentEdges.clear();

            if(direction == 0) {
                BGL_FORALL_OUTEDGES(v, e, graph, PathGraph1) {
                    adjacentEdges.push_back({e, graph[e].info.correctedJaccard()});
                }
            } else {
                BGL_FORALL_INEDGES(v, e, graph, PathGraph1) {
                    adjacentEdges.push_back({e, graph[e].info.correctedJaccard()});
                }
            }

            // Only keep the k best.
            if(adjacentEdges.size() > k) {
                std::nth_element(
                    adjacentEdges.begin(),
                    adjacentEdges.begin() + k,
                    adjacentEdges.end(),
                    OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
                adjacentEdges.resize(k);
            }

            // Mark them as to be kept.
            for(const auto& p:adjacentEdges) {
                const edge_descriptor e = p.first;
                edgesToBeKept.insert(e);
            }
        }
    }

    // Remove edges not marked as to be kept.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        if(not edgesToBeKept.contains(e)) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// Compute total base offset for this chain.
uint64_t GlobalPathGraph1::Chain::totalOffset() const
{
    uint64_t totalBaseOffset = 0;
    for(const MarkerGraphEdgePairInfo& info: infos) {
        totalBaseOffset += info.offsetInBases;
    }
    return totalBaseOffset;
}



uint64_t GlobalPathGraph1::ChainConnector::totalOffset() const
{
    uint64_t totalBaseOffset = 0;
    for(const MarkerGraphEdgePairInfo& info: infos) {
        totalBaseOffset += info.offsetInBases;
    }
    return totalBaseOffset;
}



void GlobalPathGraph1::ChainConnector::reverse()
{
    std::reverse(vertexIds.begin(), vertexIds.end());
    std::reverse(infos.begin(), infos.end());
}



// Create the connected components of this PathGraph1,
// without changing the PathGraph1 itself.
vector< shared_ptr<PathGraph1> > PathGraph1::createConnectedComponents(
    uint64_t minComponentSize) const
{
    const PathGraph1& graph = *this;

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
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1::vertex_descriptor v0 = source(e, graph);
        const PathGraph1::vertex_descriptor v1 = target(e, graph);
        disjointSets.union_set(v0, v1);
    }


    // Gather the vertices in each connected component.
    vector< shared_ptr<PathGraph1> > allComponentPointers(num_vertices(graph));
    BGL_FORALL_VERTICES(v, graph, PathGraph1) {
        const PathGraph1Vertex& vertex = graph[v];
        const uint64_t componentId = disjointSets.find_set(v);
        shared_ptr<PathGraph1>& componentPointer = allComponentPointers[componentId];
        if(not componentPointer) {
            componentPointer = make_shared<PathGraph1>();
        }
        PathGraph1& component = *componentPointer;
        component.addVertex(
            vertex.vertexId,
            vertex.edgeId);
    }


    // Gather the edges in each connected component.
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1::vertex_descriptor v0 = source(e, graph);
        const PathGraph1::vertex_descriptor v1 = target(e, graph);
        const uint64_t edgeId0 = graph[v0].edgeId;
        const uint64_t edgeId1 = graph[v1].edgeId;
        const uint64_t componentId = disjointSets.find_set(v0);
        SHASTA_ASSERT(componentId == disjointSets.find_set(v1));
        shared_ptr<PathGraph1>& componentPointer = allComponentPointers[componentId];
        SHASTA_ASSERT(componentPointer);
        PathGraph1& component = *componentPointer;
        component.addEdge(
            edgeId0,
            edgeId1,
            graph[e].info,
            graph[e].coverage);
    }



    // Keep only the components with at least minComponentSize vertices
    // and sort them by size.
    vector< pair<shared_ptr<PathGraph1>, uint64_t> > componentPointersWithSizes;
    for(const shared_ptr<PathGraph1>& p: allComponentPointers) {
        if(p) {
            const uint64_t componentSize = num_vertices(*p);
            if(componentSize >= minComponentSize) {
                componentPointersWithSizes.push_back({p, componentSize});
            }
        }
    }
    sort(componentPointersWithSizes.begin(), componentPointersWithSizes.end(),
        OrderPairsBySecondOnlyGreater<shared_ptr<PathGraph1>, uint64_t>());


    // For now return all components, including the empty ones.
    // But we want to remove the small ones and sort them by size.
    vector< shared_ptr<PathGraph1> > componentPointers;
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
void PathGraph1::removeCrossEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    uint64_t minOffset)
{
    PathGraph1& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, PathGraph1) {
        const PathGraph1Edge& edge = graph[e];

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
        BGL_FORALL_OUTEDGES(v0, e0, graph, PathGraph1) {
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
        BGL_FORALL_INEDGES(v1, e1, graph, PathGraph1) {
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

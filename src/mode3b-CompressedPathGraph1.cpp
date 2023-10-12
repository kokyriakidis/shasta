// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "Assembler.hpp"
#include "findLinearChains.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<CompressedPathGraph1>;



void CompressedPathGraph1::writeVerticesCsv() const
{
    const CompressedPathGraph1& cGraph = *this;

    ofstream csv("CompressedPathGraphVertices" + to_string(componentId) + ".csv");
    csv << "Id,Begin,End,Vertex count,Estimated length\n";

    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const CompressedPathGraph1Vertex& cVertex = cGraph[cv];
        SHASTA_ASSERT(not cVertex.v.empty());
        const PathGraph1::vertex_descriptor v0 = cVertex.v.front();
        const PathGraph1::vertex_descriptor v1 = cVertex.v.back();

        uint64_t baseOffset = totalBaseOffset(cv);

        csv << componentId << "-" << cVertex.id << ",";
        csv << graph[v0].edgeId << ",";
        csv << graph[v1].edgeId << ",";
        csv << cVertex.v.size() << ",";
        csv << baseOffset << ",";
        csv << "\n";
    }
}



// This calls the lower level function twice, with and without labels.
void CompressedPathGraph1::writeGraphviz(const string& fileNamePrefix) const
{
    bool labels = true;
    writeGraphviz(labels, fileNamePrefix);
    labels = false;
    writeGraphviz(labels, fileNamePrefix);
}



void CompressedPathGraph1::writeGraphviz(
    bool labels,
    const string& fileNamePrefix) const
{
    const CompressedPathGraph1& cGraph = *this;

    const string name = "CompressedPathGraph" + to_string(componentId);
    ofstream out(name + "-" + fileNamePrefix + (labels ? ".dot" : "-NoLabels.dot"));
    out << "digraph " << name << "{\n";

    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        out << "\"" << componentId << "-" << cGraph[cv].id << "\" [";

        if(cGraph[cv].v.size() > 1) {

            const PathGraph1::vertex_descriptor v0 = cGraph[cv].v.front();
            const PathGraph1::vertex_descriptor v1 = cGraph[cv].v.back();
            const uint64_t baseOffset = totalBaseOffset(cv);

            if(labels) {
                out <<
                    "label=\""
                    << componentId << "-" << cGraph[cv].id << "\\n" <<
                    cGraph[cv].v.size() << " vertices\\n" <<
                    "First " << graph[v0].edgeId << "\\n" <<
                    "Last " << graph[v1].edgeId << "\\n" <<
                    "Length " << baseOffset <<
                    "\"";
            }

        } else {

            SHASTA_ASSERT(cGraph[cv].v.size() == 1);

            if(labels) {
                const PathGraph1::vertex_descriptor v = cGraph[cv].v.front();
                out <<
                    "label=\"" <<
                    componentId << "-" << cGraph[cv].id << "\\n" <<
                    graph[v].edgeId << "\\n" <<
                    "\"";
            }

        }

        out << "];\n";
    }

    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        const auto cv0 = source(ce, cGraph);
        const auto cv1 = target(ce, cGraph);
        out <<
            "\"" << componentId << "-" << cGraph[cv0].id << "\"->" <<
            "\"" << componentId << "-" << cGraph[cv1].id << "\"";
        if(labels) {
            const uint64_t offset = cGraph[ce].info.offsetInBases;
            out <<
                " ["
                " label=\"" <<
                offset << "\\n" << cGraph[ce].info.common <<
                "\""
                "]";
        }
        out << ";\n";
    }
    out << "}";
}



bool CompressedPathGraph1::localTransitiveReduction(uint64_t distance)
{
    CompressedPathGraph1& cGraph = *this;

    // We want to process edges in order of increasing coverage.
    vector<pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        edgeTable.push_back({ce, cGraph[ce].info.common});
    }
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<edge_descriptor, uint64_t>());

    // Loop over all edges v0->v1 in order of increasing coverage.
    uint64_t removedEdgeCount = 0;
    for(const auto& p: edgeTable) {
        edge_descriptor e01 = p.first;
        const vertex_descriptor v0 = source(e01, cGraph);
        const vertex_descriptor v1 = target(e01, cGraph);

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
            BGL_FORALL_OUTEDGES_T(vA, eAB, cGraph, CompressedPathGraph1) {

                // Dont's use e01 in the BFS.
                if(eAB == e01) {
                    continue;
                }

                // If we reached v1, mark e01 as a nonTransitiveReduction edge
                // and stop the BFS.
                const vertex_descriptor vB = target(eAB, cGraph);
                if(vB == v1) {
                    ++removedEdgeCount;
                    boost::remove_edge(e01, cGraph);
                    endBfs = true;
                    // cout << "Reached " << v1 << endl;
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

    return removedEdgeCount > 0;
}



CompressedPathGraph1::CompressedPathGraph1(
    const PathGraph1& graph,
    uint64_t componentId,
    const Assembler& assembler) :
    MultithreadedObject<CompressedPathGraph1>(*this),
    graph(graph),
    componentId(componentId),
    assembler(assembler)
{
    CompressedPathGraph1& cGraph = *this;
    using boost::out_degree;
    using boost::in_degree;

    // Create a filtered version of the PathGraph1, containing only the
    // transitive reduction edges.
    class EdgePredicate {
    public:
        bool operator()(const PathGraph1::edge_descriptor e) const
        {
            return not (*graph)[e].isNonTransitiveReductionEdge;
        }
        EdgePredicate(const PathGraph1& graph) : graph(&graph) {}
        EdgePredicate() : graph(0) {}
    private:
        const PathGraph1* graph;
    };
    using FilteredPathGraph1 = boost::filtered_graph<PathGraph1, EdgePredicate>;
    FilteredPathGraph1 filteredGraph(graph, EdgePredicate(graph));

    vector<uint64_t> lengthHistogram(2, 0);

    // Loop over vertices of the PathGraph1 to generate the CompressedPathGraph1 vertices.
    std::set<PathGraph1::vertex_descriptor> visited;
    vector<PathGraph1::vertex_descriptor> forwardVertices;
    vector<PathGraph1::vertex_descriptor> backwardVertices;
    BGL_FORALL_VERTICES(v, filteredGraph, FilteredPathGraph1) {

        // If already visited, skip.
        if(visited.contains(v)) {
            continue;
        }
        visited.insert(v);

        // cout << "Starting at " << filteredGraph[v].edgeId << endl;

        // If in-degree or out-degree is more than 1, generate a single CompressedPathGraph1 vertex.
        if(in_degree(v, filteredGraph) > 1 or out_degree(v, filteredGraph) > 1) {
            const CompressedPathGraph1::vertex_descriptor cv = add_vertex(cGraph);
            CompressedPathGraph1Vertex& cVertex = cGraph[cv];
            cVertex.id = nextVertexId++;
            cVertex.v.push_back(v);
            // cout << "Done: is branch vertex." << endl;
            ++lengthHistogram[1];
            continue;
        }

        // Move forward until we find a branch vertex.
        bool isCircular = false;
        forwardVertices.clear();
        PathGraph1::vertex_descriptor u = v;
        while(true) {
            if(out_degree(u, filteredGraph) == 0) {
                break;
            }
            SHASTA_ASSERT(out_degree(u, filteredGraph) == 1);

            // Move forward.
            FilteredPathGraph1::out_edge_iterator it;
            tie(it, ignore) = out_edges(u, filteredGraph);
            u = target(*it, filteredGraph);

            if(in_degree(u, filteredGraph) > 1 or out_degree(u, filteredGraph) > 1) {
                break;
            }

            if(u == v) {
                isCircular = true;
            }

            forwardVertices.push_back(u);
            SHASTA_ASSERT(not visited.contains(u));
            visited.insert(u);
            // cout << "Moving forward found " << filteredGraph[u].edgeId << endl;
        }

        // For now, don't handle the circular case.
        SHASTA_ASSERT(not isCircular);

        // Move backward until we find a branch vertex.
        backwardVertices.clear();
        u = v;
        while(true) {
            if(in_degree(u, filteredGraph) == 0) {
                break;
            }
            SHASTA_ASSERT(in_degree(u, filteredGraph) == 1);

            // Move backward.
            FilteredPathGraph1::in_edge_iterator it;
            tie(it, ignore) = in_edges(u, filteredGraph);
            u = source(*it, filteredGraph);

            if(in_degree(u, filteredGraph) > 1 or out_degree(u, filteredGraph) > 1) {
                break;
            }

            SHASTA_ASSERT(u != v);

            backwardVertices.push_back(u);
            SHASTA_ASSERT(not visited.contains(u));
            visited.insert(u);
            // cout << "Moving backward found " << filteredGraph[u].edgeId << endl;
        }

        // Create a new CompressedPathGraph1 vertex.
        const CompressedPathGraph1::vertex_descriptor cv = add_vertex(cGraph);
        CompressedPathGraph1Vertex& cVertex = cGraph[cv];
        cVertex.id = nextVertexId++;
        copy(backwardVertices.rbegin(), backwardVertices.rend(), back_inserter(cVertex.v));
        cVertex.v.push_back(v);
        copy(forwardVertices.begin(),forwardVertices.end(), back_inserter(cVertex.v));

    }
    SHASTA_ASSERT(visited.size() == num_vertices(graph));



    // To generate edges, create a map giving the compressed vertex descriptor
    // that begins with a given vertex.
    std::map<PathGraph1::vertex_descriptor, CompressedPathGraph1::vertex_descriptor> vertexMap;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        vertexMap.insert({cGraph[cv].v.front(), cv});
    }
    BGL_FORALL_VERTICES(cv0, cGraph, CompressedPathGraph1) {
        const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v.back();
        BGL_FORALL_OUTEDGES(v0, e, filteredGraph, FilteredPathGraph1) {
            const PathGraph1::vertex_descriptor v1 = target(e, filteredGraph);
            auto it1 = vertexMap.find(v1);
            SHASTA_ASSERT(it1 != vertexMap.end());
            const CompressedPathGraph1::vertex_descriptor cv1 = it1->second;
            add_edge(cv0, cv1, {graph[e].info}, cGraph);
        }
    }

    uint64_t filteredGraphEdgeCount = 0;
    BGL_FORALL_EDGES(e, filteredGraph, FilteredPathGraph1) {
        ++filteredGraphEdgeCount;
    }
    cout << "The PathGraph1 has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;
    if(false) {
        cout << "The FilteredPathGraph1 has " << num_vertices(graph) << " vertices and " <<
            filteredGraphEdgeCount << " edges." << endl;
        cout << "The CompressedPathGraph1 has " << num_vertices(cGraph) << " vertices and " <<
            num_edges(cGraph) << " edges." << endl;
    }

}



bool CompressedPathGraph1::mergeLinearChains()
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;

    vector< vector<vertex_descriptor> > linearChains;
    findLinearVertexChains(cGraph, linearChains);

    bool changesWereMade = false;
    for(const auto& linearChain: linearChains) {
        if(linearChain.size() == 1) {
            continue;
        }

        // The first and last vertex can have in/out degrees greater than one
        // and will not be merged.
        vector<vertex_descriptor> verticesToBeMerged;
        for(uint64_t i=0; i<linearChain.size(); i++) {
            const vertex_descriptor cv = linearChain[i];
            if(i==0 or i==linearChain.size()-1) {
                if(in_degree(cv, cGraph)>1 or out_degree(cv, cGraph)>1) {
                    continue;
                }
            }
            verticesToBeMerged.push_back(cv);
        }
        if(verticesToBeMerged.size() < 2) {
            continue;
        }
        changesWereMade = true;


        // Create the new vertex.
        const vertex_descriptor cvNew = add_vertex(cGraph);
        auto& cVertexNew = cGraph[cvNew];
        cVertexNew.id = nextVertexId++;
        for(const vertex_descriptor cv: verticesToBeMerged) {
            const auto& cVertex = cGraph[cv];
            copy(cVertex.v.begin(), cVertex.v.end(), back_inserter(cVertexNew.v));
        }

        if(debug) {
            cout << componentId << "-" << cVertexNew.id << " created by merging";
            for(const auto cv: verticesToBeMerged) {
                cout << " " << componentId << "-" << cGraph[cv].id;
            }
            cout << endl;
        }

        // Reroute in-edges.
        BGL_FORALL_INEDGES(verticesToBeMerged.front(), ce, cGraph, CompressedPathGraph1) {
            const auto& oldEdge = cGraph[ce];
            add_edge(source(ce, cGraph), cvNew, oldEdge, cGraph);
        }

        // Reroute out-edges.
        BGL_FORALL_OUTEDGES(verticesToBeMerged.back(), ce, cGraph, CompressedPathGraph1) {
            const auto& oldEdge = cGraph[ce];
            add_edge(cvNew, target(ce, cGraph), oldEdge, cGraph);
        }

        // Now we can remove the vertices we merged.
        for(const vertex_descriptor cv: verticesToBeMerged) {
            clear_vertex(cv, cGraph);
            remove_vertex(cv, cGraph);
        }

    }

    if(debug) {
        cout << "After merging linear chains, the CompressedPathGraph1 has " <<
            num_vertices(cGraph) << " vertices  and " <<
            num_edges(cGraph) << " edges." << endl;
    }

    return changesWereMade;
}



// An edge cv0->cv1 is a cross edge if:
// - Has coverage <= threshold1
// - cv0 has out-degree > 1 and at least one out-edge with coverage >= threshold2.
// - cv1 has in-degree  > 1 and at least one in-edge  with coverage >= threshold2.
bool CompressedPathGraph1::removeCrossEdges(
    uint64_t threshold1,
    uint64_t threshold2)
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;
    SHASTA_ASSERT(threshold2 > threshold1);



    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {

        // Check coverage.
        if(cGraph[ce].info.common > threshold1) {
            continue;
        }

        // Check degrees.
        const vertex_descriptor cv0 = source(ce, cGraph);
        const vertex_descriptor cv1 = target(ce, cGraph);
        if(out_degree(cv0, cGraph) < 2) {
            continue;
        }
        if(in_degree(cv1, cGraph) < 2) {
            continue;
        }

        // Check that v0 has at least one out-edge with coverage >= threshold2.
        bool found0 = false;
        BGL_FORALL_OUTEDGES(cv0, ce02, cGraph, CompressedPathGraph1) {
            if(cGraph[ce02].info.common >= threshold2) {
                found0 = true;
                break;
            }
        }
        if(not found0) {
            continue;
        }

        // Check that v1 has at least one in-edge with coverage >= threshold2.
        bool found1 = false;
        BGL_FORALL_INEDGES(cv1, ce21, cGraph, CompressedPathGraph1) {
            if(cGraph[ce21].info.common >= threshold2) {
                found1 = true;
                break;
            }
        }
        if(not found1) {
            continue;
        }

        edgesToBeRemoved.push_back(ce);

        if(debug) {
            cout << "Cross-edge " <<
                componentId << "-" << cGraph[cv0].id << "->" <<
                componentId << "-" << cGraph[cv1].id <<
                " will be removed." << endl;
        }

    }



    for(const edge_descriptor ce: edgesToBeRemoved) {
        boost::remove_edge(ce, cGraph);
    }

    return not edgesToBeRemoved.empty();
}



string CompressedPathGraph1::vertexIdString(vertex_descriptor cv) const
{
    const CompressedPathGraph1& cGraph = *this;

    return to_string(componentId) + "-" + to_string(cGraph[cv].id);
}



bool CompressedPathGraph1::detangleVertices(uint64_t detangleTolerance)
{
    CompressedPathGraph1& cGraph = *this;

    // Find the vertices that can potentially be detangled.
    vector<vertex_descriptor> detangleCandidates;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        if(in_degree(cv, cGraph) > 1 and out_degree(cv, cGraph) > 1) {
            detangleCandidates.push_back(cv);
        }
    }

    // Do the detangling.
    uint64_t detangledCount = 0;
    for(const auto cv: detangleCandidates) {
        bool wasDetangled = detangleVertex(cv, detangleTolerance);
        if(wasDetangled) {
            ++detangledCount;
        }
    }

    return detangledCount > 0;
}



// Attempt to detangle a compressed vertex.
// The detangle operation leaves all in-degree and out-degrees
// of all other vertices unchanged.
bool CompressedPathGraph1::detangleVertex(
    vertex_descriptor cv,
    uint64_t detangleTolerance)
{
    CompressedPathGraph1& cGraph = *this;

    const bool debug = false;
    if(debug) {
        cout << "Attempting to detangle " << componentId << "-" << cGraph[cv].id << endl;
    }

    // Gather the source vertices of incoming edges.
    vector<vertex_descriptor> incoming;
    BGL_FORALL_INEDGES(cv, e, cGraph, CompressedPathGraph1) {
        incoming.push_back(source(e, cGraph));
    }

    // Gather the target vertices of outgoing edges.
    vector<vertex_descriptor> outgoing;
    BGL_FORALL_OUTEDGES(cv, e, cGraph, CompressedPathGraph1) {
        outgoing.push_back(target(e, cGraph));
    }

    if(debug) {
        cout << "Incoming:";
        for(const auto cv1: incoming) {
            const auto& cVertex = cGraph[cv1];
            const PathGraph1::vertex_descriptor v1 = cVertex.v.back();
            cout << " " << graph[v1].edgeId;
        }
        cout << endl;

        cout << "Outgoing:";
        for(const auto cv1: outgoing) {
            const auto& cVertex = cGraph[cv1];
            const PathGraph1::vertex_descriptor v1 = cVertex.v.front();
            cout << " " << graph[v1].edgeId;
        }
        cout << endl;
    }

    // We can only detangle if the number of incoming edges
    // equals the number of outgoing edges.
    if(incoming.size() != outgoing.size()) {
        return false;
    }

    // For each incoming edge, count the number of
    // outgoing edges it has more than detangleTolerance common oriented reads with.
    vector<uint64_t> inCount(incoming.size(), 0);

    // For each outgoing edge, count the number of
    // incoming edges it has  more than detangleTolerance common oriented reads with.
    vector<uint64_t> outCount(outgoing.size(), 0);

    // Loop over pairs of incoming/outgoing edges.
    MarkerGraphEdgePairInfo info;
    class NewEdge {
    public:
        vertex_descriptor cv0;
        vertex_descriptor cv1;
        MarkerGraphEdgePairInfo info;
    };
    vector<NewEdge> newEdges;
    for(uint64_t i0=0; i0<incoming.size(); i0++) {
        const auto cv0 = incoming[i0];
        const auto& cVertex0 = cGraph[cv0];
        const PathGraph1::vertex_descriptor v0 = cVertex0.v.back();
        const MarkerGraphEdgeId edgeId0 = graph[v0].edgeId;
        for(uint64_t i1=0; i1<outgoing.size(); i1++) {
            const auto cv1 = outgoing[i1];
            const auto& cVertex1 = cGraph[cv1];
            const PathGraph1::vertex_descriptor v1 = cVertex1.v.front();
            const MarkerGraphEdgeId edgeId1 = graph[v1].edgeId;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
            if(debug) {
                cout << edgeId0 << " " << edgeId1 << ": " << info.common << endl;
            }
            if(info.common > detangleTolerance) {
                ++inCount[i0];
                ++outCount[i1];
                newEdges.push_back({cv0, cv1, info});
            }
        }
    }

    if(debug) {
        cout << "inCount ";
        copy(inCount.begin(), inCount.end(), ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
        cout << "outCount ";
        copy(outCount.begin(), outCount.end(), ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
    }

    // We can only detangle if all the inCount and outCount entries are 1.
    // This means that for each incoming edge there is only one outgoing edge
    // with common oriented reads, and vice versa.
    for(const uint64_t c: inCount) {
        if(c != 1) {
            if(debug) {
                cout << "Cannot detangle." << endl;
            }
            return false;
        }
    }
    for(const uint64_t c: outCount) {
        if(c != 1) {
            if(debug) {
                cout << "Cannot detangle." << endl;
            }
            return false;
        }
    }

    // Add the edges.
    for(const auto& newEdge: newEdges) {
        edge_descriptor ce;

        // Check that the edge does not already exists.
        bool edgeExists = false;
        tie(ce, edgeExists) = boost::edge(newEdge.cv0, newEdge.cv1, cGraph);
        if(edgeExists) {
            continue;
        }

        bool edgeWasAdded = false;
        tie(ce, edgeWasAdded) = boost::add_edge(newEdge.cv0, newEdge.cv1, {newEdge.info}, cGraph);
        SHASTA_ASSERT(edgeWasAdded);
        if(debug) {
            const auto& cEdge = cGraph[ce];
            cout << "Added compressed edge " <<
                componentId << "-" << cGraph[newEdge.cv0].id << "->" <<
                componentId << "-" << cGraph[newEdge.cv1].id << ", common count " <<
                cEdge.info.common << ", offset " <<
                cEdge.info.offsetInBases << endl;
        }
    }

    // Now we can remove the vertex we detangled.
    clear_vertex(cv, cGraph);
    remove_vertex(cv, cGraph);

    return true;
}



bool CompressedPathGraph1::detangleLinearChains(uint64_t detangleTolerance)
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;

    // Find the linear chains.
    vector< vector<CompressedPathGraph1::edge_descriptor> > linearChains;
    findLinearChains(cGraph, 1, linearChains);
    if(debug) {
        cout << "Found " << linearChains.size() << " linear chains." << endl;
    }



    // Try to detangle all linear chains.
    uint64_t detangleCount = 0;
    for(const auto& linearChain: linearChains) {
        SHASTA_ASSERT(not linearChain.empty());

        // If the first vertex has in-degree less than 2, do nothing.
        const auto ce0 = linearChain.front();
        const auto cv0 = source(ce0, cGraph);
        const uint64_t inDegree = in_degree(cv0, cGraph);
        if(inDegree < 2) {
            continue;
        }

        // If the last vertex has out-degree less than 2, do nothing.
        const auto ce1 = linearChain.back();
        const auto cv1 = target(ce1, cGraph);
        const uint64_t outDegree = out_degree(cv1, cGraph);
        if(outDegree < 2) {
            continue;
        }

        // We can only detangle if the in-degree and out-degree are the same.
        if(inDegree != outDegree) {
            continue;
        }

        // Gather the vertices of this chain.
        vector<vertex_descriptor> chainVertices;
        chainVertices.push_back(cv0);
        for(const auto ce: linearChain) {
            chainVertices.push_back(target(ce, cGraph));
        }

        // Sanity check: all vertices except the first and last must have
        // in-degree and out-degree 1.
        for(uint64_t i=1; i<chainVertices.size()-1; i++) {
            const auto cv = chainVertices[i];
            SHASTA_ASSERT(in_degree(cv, cGraph) == 1);
            SHASTA_ASSERT(out_degree(cv, cGraph) == 1);
        }

        // If the first vertex has out-degree>1, we cannot detangle.
        if(out_degree(cv0, cGraph) > 1) {
            continue;
        }
        // If the last vertex has in-degree>1, we cannot detangle.
        if(in_degree(cv1, cGraph) > 1) {
            continue;
        }

        if(debug) {
            cout << "Attempting to detangle linear chain:";
            for(const auto cv: chainVertices) {
                cout << " " << componentId << "-" << cGraph[cv].id;
            }
            cout << endl;
        }

        // Gather the source vertices of incoming edges.
        vector<vertex_descriptor> incoming;
        BGL_FORALL_INEDGES(cv0, e, cGraph, CompressedPathGraph1) {
            incoming.push_back(source(e, cGraph));
        }

        // Gather the target vertices of outgoing edges.
        vector<vertex_descriptor> outgoing;
        BGL_FORALL_OUTEDGES(cv1, e, cGraph, CompressedPathGraph1) {
            outgoing.push_back(target(e, cGraph));
        }

        if(debug) {
            cout << "Incoming:";
            for(const auto cv1: incoming) {
                const auto& cVertex = cGraph[cv1];
                const PathGraph1::vertex_descriptor v1 = cVertex.v.back();
                cout << " " << graph[v1].edgeId;
            }
            cout << endl;

            cout << "Outgoing:";
            for(const auto cv1: outgoing) {
                const auto& cVertex = cGraph[cv1];
                const PathGraph1::vertex_descriptor v1 = cVertex.v.front();
                cout << " " << graph[v1].edgeId;
            }
            cout << endl;
        }
        // For each incoming edge, count the number of
        // outgoing edges it has more than detangleTolerance common oriented reads with.
        vector<uint64_t> inCount(incoming.size(), 0);

        // For each outgoing edge, count the number of
        // incoming edges it has more than detangleTolerance common oriented reads with.
        vector<uint64_t> outCount(outgoing.size(), 0);

        // Loop over pairs of incoming/outgoing edges.
        MarkerGraphEdgePairInfo info;
        class NewEdge {
        public:
            vertex_descriptor cv0;
            vertex_descriptor cv1;
            MarkerGraphEdgePairInfo info;
        };

        vector<NewEdge> newEdges;
        for(uint64_t i0=0; i0<incoming.size(); i0++) {
            const auto cv0 = incoming[i0];
            const auto& cVertex0 = cGraph[cv0];
            const PathGraph1::vertex_descriptor v0 = cVertex0.v.back();
            const MarkerGraphEdgeId edgeId0 = graph[v0].edgeId;
            for(uint64_t i1=0; i1<outgoing.size(); i1++) {
                const auto cv1 = outgoing[i1];
                const auto& cVertex1 = cGraph[cv1];
                const PathGraph1::vertex_descriptor v1 = cVertex1.v.front();
                const MarkerGraphEdgeId edgeId1 = graph[v1].edgeId;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
                if(debug) {
                    cout << edgeId0 << " " << edgeId1 << ": " << info.common << endl;
                }
                if(info.common > detangleTolerance) {
                    ++inCount[i0];
                    ++outCount[i1];
                    newEdges.push_back({cv0, cv1, info});
                }
            }
        }

        if(debug) {
            cout << "inCount ";
            copy(inCount.begin(), inCount.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
            cout << "outCount ";
            copy(outCount.begin(), outCount.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
        }

        // We can only detangle if all the inCount and outCount entries are 1.
        // This means that for each incoming edge there is only one outgoing edge
        // with common oriented reads, and vice versa.
        bool canDetangle = true;
        for(const uint64_t c: inCount) {
            if(c != 1) {
                if(debug) {
                    cout << "Cannot detangle." << endl;
                }
                canDetangle = false;
                break;
            }
        }
        for(const uint64_t c: outCount) {
            if(c != 1) {
                if(debug) {
                    cout << "Cannot detangle." << endl;
                }
                canDetangle = false;
                break;
            }
        }
        if(not canDetangle) {
            continue;
        }

        // Add the edges.
        for(const auto& newEdge: newEdges) {

            // Check that the edge does not already exists.
            edge_descriptor ce;
            bool edgeExists = false;
            tie(ce, edgeExists) = boost::edge(newEdge.cv0, newEdge.cv1, cGraph);
            if(edgeExists) {
                continue;
            }

            bool edgeWasAdded = false;
            tie(ce, edgeWasAdded) = boost::add_edge(newEdge.cv0, newEdge.cv1, {newEdge.info}, cGraph);
            SHASTA_ASSERT(edgeWasAdded);
            if(debug) {
                const auto& cEdge = cGraph[ce];
                cout << "Added compressed edge " <<
                    componentId << "-" << cGraph[newEdge.cv0].id << "->" <<
                    componentId << "-" << cGraph[newEdge.cv1].id << ", common count " <<
                    cEdge.info.common << ", offset " <<
                    cEdge.info.offsetInBases << endl;
            }
        }

        // Now we can remove the vertices of the linear chain we detangled.
        for(const auto cv: chainVertices) {
            clear_vertex(cv, cGraph);
            remove_vertex(cv, cGraph);
        }
        ++detangleCount;
    }

    return detangleCount > 0;
}



#if 0
bool CompressedPathGraph1::detangleSuperbubbles(uint64_t minReliableLength)
{
    CompressedPathGraph1& cGraph = *this;
    const bool debug = false;

    // The edges to be added.
    class NewEdge {
    public:
        vertex_descriptor cv0;
        vertex_descriptor cv1;
        CompressedPathGraph1Edge edge;
    };
    vector<NewEdge> newEdges;

    // Loop over long vertices of the CompressedPathGraph1.
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const uint64_t baseOffset = totalBaseOffset(cv);
        if(baseOffset < minReliableLength) {
            continue;
        }
        if(debug) {
            cout << "Starting BFS at " << componentId << "-" << cGraph[cv].id << endl;
        }

        // Do a BFS starting here, stopping when we encounter a long vertices.
        std::queue<vertex_descriptor> q;
        q.push(cv);
        std::set<vertex_descriptor> visited;
        visited.insert(cv);
        while(not q.empty()) {
            const auto cv0 = q.front();
            q.pop();
            if(false) {
                cout << "Dequeued " << componentId << "-" << cGraph[cv0].id << endl;
            }

            BGL_FORALL_OUTEDGES(cv0, ce, cGraph, CompressedPathGraph1) {
                const auto cv1 = target(ce, cGraph);
                if(visited.contains(cv1)) {
                    continue;
                }
                const uint64_t baseOffset = totalBaseOffset(cv1);
                if(baseOffset >= minReliableLength) {
                    CompressedPathGraph1Edge newEdge;
                    assembler.analyzeMarkerGraphEdgePair(
                        graph[cGraph[cv].v.back()].edgeId,
                        graph[cGraph[cv1].v.front()].edgeId,
                        newEdge.info);
                    if(newEdge.info.common) {
                        newEdges.push_back({cv, cv1, newEdge});
                    }
                    if(debug) {
                        cout << "Found " << componentId << "-" << cGraph[cv1].id <<
                            ", common count " << newEdge.info.common << endl;
                    }
                }   else {
                    q.push(cv1);
                    visited.insert(cv1);
                    if(false) {
                        cout << "Enqueued " << componentId << "-" << cGraph[cv1].id << endl;
                    }
                }
            }
        }
    }



    // Now we can:
    // - Remove all short vertices.
    // - Remove all edges.
    // - Replace them with the new edges we found.

    // Remove all short vertices.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const uint64_t baseOffset = totalBaseOffset(cv);
        if(baseOffset < minReliableLength) {
            verticesToBeRemoved.push_back(cv);
        }
    }
    for(const auto cv: verticesToBeRemoved) {
        clear_vertex(cv, cGraph);
        remove_vertex(cv, cGraph);
    }

    // Remove all edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        edgesToBeRemoved.push_back(ce);
    }
    for(const auto ce: edgesToBeRemoved) {
        boost::remove_edge(ce, cGraph);
    }

    // Replace them with the new edges we found.
    for(const auto& newEdge: newEdges) {
        add_edge(newEdge.cv0, newEdge.cv1, {newEdge.edge}, cGraph);
    }
    if(debug) {
        cout << "After superbubble detangling, the CompressedPathGraph1 has " <<
            num_vertices(cGraph) << " vertices  and " <<
            num_edges(cGraph) << " edges." << endl;
    }

    return not verticesToBeRemoved.empty(); // Questionable.
}
#endif



uint64_t CompressedPathGraph1::totalBaseOffset(vertex_descriptor cv) const
{
    const CompressedPathGraph1& cGraph = *this;
    const CompressedPathGraph1Vertex& cVertex = cGraph[cv];

    SHASTA_ASSERT(not cVertex.v.empty());

    uint64_t totalBaseOffset = 0;

    for(uint64_t i=1; i<cVertex.v.size(); i++) {
        const PathGraph1::vertex_descriptor vA = cVertex.v[i-1];
        const PathGraph1::vertex_descriptor vB = cVertex.v[i];

        PathGraph1::edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(vA, vB, graph);
        if(edgeWasFound) {
            if(graph[e].info.offsetInBases > 0) {
                totalBaseOffset += graph[e].info.offsetInBases;
            }
        } else {
            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(
                graph[vA].edgeId, graph[vB].edgeId, info));
            if(info.offsetInBases > 0) {
                totalBaseOffset += info.offsetInBases;
            }
        }
    }
    return totalBaseOffset;

}



uint64_t CompressedPathGraph1::totalVertexBaseOffset() const
{
    const CompressedPathGraph1& cGraph = *this;
    uint64_t totalOffset = 0;
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        totalOffset += totalBaseOffset(cv);
    }
    return totalOffset;
}



uint64_t CompressedPathGraph1::totalEdgeBaseOffset() const
{
    const CompressedPathGraph1& cGraph = *this;
    uint64_t totalOffset = 0;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        totalOffset += cGraph[ce].info.offsetInBases;
    }
    return totalOffset;
}



void CompressedPathGraph1::detangle()
{

#if 1

    const uint64_t compressedTransitiveReductionDistance = 100;

    uint64_t detangleTolerance = 0;
    detangleIteration(
        "A",
        compressedTransitiveReductionDistance,
        detangleTolerance);

    detangleTolerance = 1;
    detangleIteration(
        "B",
        compressedTransitiveReductionDistance,
        detangleTolerance);


#endif


#if 0

    // This is the original detangle code that works well for HG002-UL-R10-Dec2022-test1
    const uint64_t compressedTransitiveReductionDistance = 100;
    const uint64_t minReliableLength = 1000;
    const uint64_t crossEdgeCoverageThreshold1 = 1;
    const uint64_t crossEdgeCoverageThreshold2 = 2;
    const uint64_t detangleTolerance = 0;

    // Strict detangling, with zero detangle tolerance.
    detangleIteration(
        "A",
        compressedTransitiveReductionDistance,
        detangleTolerance,
        minReliableLength);

    // Remove cross-edges, then detangle again, this time with a looser detangle tolerance.
    removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
    detangleIteration(
        "B",
        compressedTransitiveReductionDistance,
        detangleTolerance,
        minReliableLength);

    // EXPERIMENT
    detangleSuperbubbles(2000);
    removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
    detangleIteration(
        "C",
        compressedTransitiveReductionDistance,
        detangleTolerance,
        minReliableLength);
    detangleSuperbubbles(5000);
    removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
    detangleIteration(
        "D",
        compressedTransitiveReductionDistance,
        detangleTolerance,
        minReliableLength);
    detangleSuperbubbles(10000);
    removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
    detangleIteration(
        "E",
        compressedTransitiveReductionDistance,
        detangleTolerance,
        minReliableLength);
    detangleSuperbubbles(20000);
    removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
    detangleIteration(
        "F",
        compressedTransitiveReductionDistance,
        detangleTolerance,
        minReliableLength);
    removeCrossEdges(crossEdgeCoverageThreshold1, crossEdgeCoverageThreshold2);
#endif
}



void CompressedPathGraph1::detangleIteration(
    const string& name,     // For graphviz output
    uint64_t compressedTransitiveReductionDistance,
    uint64_t detangleTolerance)
{

    for(uint64_t iteration=0; ; ++iteration) {

        // Try everything.
        const bool transitiveReductionChanges = localTransitiveReduction(compressedTransitiveReductionDistance);
        writeGfaAndGraphviz(name + "-" + to_string(iteration) + "-0");
        const bool detangleVerticesChanges = detangleVertices(detangleTolerance);
        writeGfaAndGraphviz(name + "-" + to_string(iteration) + "-1");
        const bool detangleLinearChainsChanges = detangleLinearChains(detangleTolerance);
        writeGfaAndGraphviz(name + "-" + to_string(iteration) + "-2");
        const bool mergeLinearChainsChanges = mergeLinearChains();
        writeGfaAndGraphviz(name + "-" + to_string(iteration) + "-3");

        // If nothing changed, stop the iteration.
        if(not (
            transitiveReductionChanges or
            detangleVerticesChanges or
            detangleLinearChainsChanges or
            mergeLinearChainsChanges
            )) {
            break;
        }
    }

}



void CompressedPathGraph1::assembleVertices(
    uint64_t threadCount0,
    uint64_t threadCount1)
{
    const CompressedPathGraph1& cGraph = *this;

    // Store the information needed by the threads.
    assembleVerticesData.threadCount1 = threadCount1;
    assembleVerticesData.allVertices.clear();
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        assembleVerticesData.allVertices.push_back(cv);
    }

    // High level parallelization uses threadCount0 threads.
    // Low level parallelization uses threadCount1 threads.
    setupLoadBalancing(assembleVerticesData.allVertices.size(), 1);
    runThreads(&CompressedPathGraph1::assembleVerticesThreadFunction, threadCount0);

    ofstream fasta("Component" + to_string(componentId) + ".fasta");
    ofstream csv("Component" + to_string(componentId) + ".csv");
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {
        const CompressedPathGraph1Vertex& cVertex = cGraph[cv];
        const string name = vertexIdString(cv);
        cVertex.assemblyPath->writeFasta(fasta, name);
        cVertex.assemblyPath->writeCsv(csv, name);
    }
}



void CompressedPathGraph1::assembleVerticesThreadFunction(uint64_t threadId)
{
    // Number of threads for lower level parallelization.
    const uint64_t threadCount1 = assembleVerticesData.threadCount1;

    // High level parallel loop over all vertices.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << i << "/" << assembleVerticesData.allVertices.size() << endl;
            }
            const vertex_descriptor cv = assembleVerticesData.allVertices[i];
            assembleVertex(cv, threadCount1);
        }
    }
}



void CompressedPathGraph1::assembleVertex(
    vertex_descriptor cv,
    uint64_t threadCount1)
{
    CompressedPathGraph1& cGraph = *this;
    CompressedPathGraph1Vertex& cVertex = cGraph[cv];

    // Construct the MarkerGraphEdgeIds of the assembly path.
    vector<MarkerGraphEdgeId> edgeIds;
    for(const PathGraph1::vertex_descriptor v: cVertex.v) {
        edgeIds.push_back(graph[v].edgeId);
    }

    // Construct the MarkerGraphEdgePairInfo.
    vector<MarkerGraphEdgePairInfo> infos(edgeIds.size() - 1);
    for(uint64_t i=1; i<edgeIds.size(); i++) {
        const MarkerGraphEdgeId edgeId0 = edgeIds[i-1];
        const MarkerGraphEdgeId edgeId1 = edgeIds[i];
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, infos[i-1]));
    }

    // Create the AssemblyPath.
    cVertex.assemblyPath = make_shared<AssemblyPath>(assembler, edgeIds, infos, threadCount1);
}



void CompressedPathGraph1::writeGfa(const string& fileNamePrefix) const
{
    const CompressedPathGraph1& cGraph = *this;

    const uint64_t vertexBases = totalVertexBaseOffset();
    const uint64_t edgeBases = totalEdgeBaseOffset();
    const uint64_t totalBases = vertexBases + edgeBases;

    cout << fileNamePrefix << ": " <<
        num_vertices(cGraph) << " vertices (" << vertexBases << " bases), " <<
        num_edges(cGraph) << " edges (" << edgeBases << " bases), " << totalBases << " bases total." << endl;

    const string name = "CompressedPathGraph" + to_string(componentId);
    ofstream gfa(name + "-" + fileNamePrefix + ".gfa");

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each vertex.
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1) {

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << vertexIdString(cv) << "\t";

        // Sequence.
        gfa << "*\t";

        // Sequence length in bases.
        gfa << "LN:i:" << totalBaseOffset(cv) << "\n";

    }



    // Write a link for each edge.
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1) {
        const vertex_descriptor cv0 = source(ce, cGraph);
        const vertex_descriptor cv1 = target(ce, cGraph);

        gfa <<
            "L\t" <<
            vertexIdString(cv0) << "\t+\t" <<
            vertexIdString(cv1) << "\t+\t*\n";

    }
}



void CompressedPathGraph1::writeGfaAndGraphviz(const string& fileNamePrefix) const
{
    writeGfa(fileNamePrefix);
    writeGraphviz(fileNamePrefix);
}



bool CompressedPathGraph1::detangleKnots()
{
    vector<Knot> knots;
    findKnots(knots);
    return true;
}


// NOTE THIS DOES NOT FIND ALL POSSIBLE KNOTS.
// SEE BELOW FOR MORE INFORMATION.
void CompressedPathGraph1::findKnots(vector<Knot>& knots) const
{
    const CompressedPathGraph1& cGraph = *this;
    const bool debug = true;
    knots.clear();

    BGL_FORALL_VERTICES(cv0, cGraph, CompressedPathGraph1) {

        // If the in-degree is too small, skip this cv0.
        const uint64_t m = in_degree(cv0, cGraph);
        if(m < 2) {
            continue;
        }

        if(false) {
            cout << "Looking for knots starting at " << vertexIdString(cv0) << endl;
        }

        // Initialize a possible Knot starting at cv0.
        Knot knot;
        knot.cv0 = cv0;
        BGL_FORALL_INEDGES(cv0, ce, cGraph, CompressedPathGraph1) {
            knot.sourceVertices.push_back(source(ce, cGraph));
        }
        knot.tangleMatrix.resize(m, vector<uint64_t>(m));

        // To see if cv0 is the start of a Knot, we do a BFS starting at cv0.
        // A possible cv1 is encountered if:
        // - The queue is empty.
        // - The in-degree of cv1 is m.
        // - The tangle matrix between cv0 and cv1 is approximately diagonal.
        // THIS ALGORITHM ONLY FINDS A SUBSET OF THE KNOTS.
        // CHECKING FOR AN EMPTY QUEUE IS TOO STRINGENT.
        // THIS MIGHT BE SUFFICIENT FOR DENTAGLING, BUT IF NOT
        // USE AN ALGORITHM BASED ON LOCAL DOMINATOR TREES.

        std::queue<vertex_descriptor> q;
        q.push(cv0);
        std::set<vertex_descriptor> verticesEncountered;
        verticesEncountered.insert(cv0);
        while(not q.empty()) {
            const vertex_descriptor cv1 = q.front();
            q.pop();

            // See if cv1 is a valid cv1 for a Knot starting at cv0.
            if(q.empty() and out_degree(cv1, cGraph) == m) {
                knot.cv1 = cv1;
                BGL_FORALL_OUTEDGES(cv1, ce, cGraph, CompressedPathGraph1) {
                    knot.targetVertices.push_back(target(ce, cGraph));
                }

                // Compute the tangle matrix.
                for(uint64_t i=0; i<m; i++) {
                    const vertex_descriptor cvA = knot.sourceVertices[i];
                    const auto& cVertexA = cGraph[cvA];
                    const PathGraph1::vertex_descriptor vA = cVertexA.v.back();
                    const MarkerGraphEdgeId edgeIdA = graph[vA].edgeId;
                    for(uint64_t j=0; j<m; j++) {
                        const vertex_descriptor cvB = knot.targetVertices[j];
                        const auto& cVertexB = cGraph[cvB];
                        const PathGraph1::vertex_descriptor vB = cVertexB.v.front();
                        const MarkerGraphEdgeId edgeIdB = graph[vB].edgeId;
                        MarkerGraphEdgePairInfo info;
                        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
                        knot.tangleMatrix[i][j] = info.common;
                    }
                }


                // Count the number of non-zero tangle matrix elements in each row and column.
                // If all of these are 1, we can detangle.
                bool canDetangle = true;
                for(uint64_t i=0; i<m; i++) {
                    uint64_t nonZeroCount = 0;
                    for(uint64_t j=0; j<m; j++) {
                        if(knot.tangleMatrix[i][j] != 0) {
                            ++nonZeroCount;
                        }
                    }
                    if(nonZeroCount != 1) {
                        canDetangle = false;
                    }
                }
                if(canDetangle) {
                    for(uint64_t j=0; j<m; j++) {
                        uint64_t nonZeroCount = 0;
                        for(uint64_t i=0; i<m; i++) {
                            if(knot.tangleMatrix[i][j] != 0) {
                                ++nonZeroCount;
                            }
                        }
                        if(nonZeroCount != 1) {
                            canDetangle = false;
                        }
                    }
                }

                if(canDetangle) {
                    if(debug) {
                        cout << "Knot at " <<
                            vertexIdString(knot.cv0) << " " << vertexIdString(knot.cv1) << endl;
                        for(uint64_t i=0; i<m; i++) {
                            for(uint64_t j=0; j<m; j++) {
                                cout << knot.tangleMatrix[i][j] << " ";
                            }
                            cout << endl;
                        }
                    }
                }

                knot.targetVertices.clear();

                // If the tangle matrix is all zero, stop the BFS.
                bool foundNonZero = false;
                for(const auto& v: knot.tangleMatrix) {
                    for(const uint64_t t: v) {
                        if(t > 0) {
                            foundNonZero = true;
                            break;
                        }
                    }
                    if(foundNonZero) {
                        break;
                    }
                }
                if(not foundNonZero) {
                    break;
                }
            }

            // Enqueue the children of cv1 we did to already encounter.
            BGL_FORALL_OUTEDGES(cv1, e, cGraph, CompressedPathGraph1) {
                const vertex_descriptor cv2 = target(e, cGraph);
                if(not verticesEncountered.contains(cv2)) {
                    verticesEncountered.insert(cv2);
                    q.push(cv2);
                }
            }
        }
    }
}

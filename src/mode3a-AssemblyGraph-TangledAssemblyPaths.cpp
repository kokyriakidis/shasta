// Shasta.
#include "mode3a-AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
using namespace shasta;
using namespace mode3a;

// Boost lbraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

// Standard library.
#include "array.hpp"
#include "fstream.hpp"



void AssemblyGraph::computeTangledAssemblyPaths(uint64_t threadCount)
{
    cout << "computeTangledAssemblyPaths begins." << endl;

    // Initialize a TangledAssemblyPath for each of the
    // longest path computed by analyzePartialPaths.
    tangledAssemblyPaths.clear();
    tangledAssemblyPaths.resize(analyzePartialPathsData.longestPaths.size());

    // Compute the TangledAssemblyPaths in parallel.
    setupLoadBalancing(tangledAssemblyPaths.size(), 1);
    runThreads(&AssemblyGraph::computeTangledAssemblyPathsThreadFunction, threadCount);
}



void AssemblyGraph::computeTangledAssemblyPathsThreadFunction(uint64_t threadId)
{
    ofstream debugOut("ComputeTangledAssemblyPaths-Thread" + to_string(threadId) + ".txt");

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; ++i) {
            computeTangledAssemblyPath(
                analyzePartialPathsData.longestPaths[i],
                tangledAssemblyPaths[i],
                debugOut);
        }
    }

}



void AssemblyGraph::computeTangledAssemblyPath(
    const vector<vertex_descriptor>& primaryVertices,
    TangledAssemblyPath& tangledAssemblyPath,
    ostream& debugOut
    )
{
    // Store the primary vertices.
    tangledAssemblyPath.primaryVertices = primaryVertices;

    // Compute the secondary vertices for each pair of primary vertices.
    tangledAssemblyPath.secondaryVertices.clear();
    tangledAssemblyPath.secondaryVertices.resize(tangledAssemblyPath.primaryVertices.size() - 1);

    uint64_t linearCount = 0;
    for(uint64_t i=0; i<tangledAssemblyPath.secondaryVertices.size(); i++) {
        const bool isLinear = computeSecondaryVertices(
            tangledAssemblyPath.primaryVertices[i],
            tangledAssemblyPath.primaryVertices[i+1],
            tangledAssemblyPath.secondaryVertices[i],
            debugOut);
        if(isLinear) {
            ++linearCount;
        }
    }
    if(debugOut) {
        debugOut << "Tangled assembly path " <<
            vertexStringId(primaryVertices.front()) << "..." <<
            vertexStringId(primaryVertices.back()) << ": " <<
            linearCount << " linear out of " << primaryVertices.size() << "\n";
    }

}



// Given a pair of primary vertices in a tangled assembly path,
// compute the intervening secondary vertices.
// The code is similar to AssemblyGraph::computePartialPath2,
// but instead of using entire oriented read journeys,
// it only uses portions between v0 and v1.
bool AssemblyGraph::computeSecondaryVertices(
    vertex_descriptor v0,
    vertex_descriptor v1,
    vector<TangledAssemblyPath::SecondaryVertexInfo>& secondaryVertices,
    ostream& debugOut)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minLinkCoverage = 3;

    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];

#if 0
    const bool debug =
        vertex0.segmentId == 23361 and
        vertex1.segmentId == 7997;
#else
    const bool debug = true;
#endif

    if(debug and debugOut) {
        debugOut << "computeSecondaryVertices begins for " <<
            vertexStringId(v0) << " " << vertexStringId(v1) << endl;
    }

    secondaryVertices.clear();

    // The vertices we encounter when following the oriented read journeys.
    vector<vertex_descriptor> verticesEncountered;

    // The transitions we encounter when following the oriented read journeys.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitionsEncountered;

    // Joint loop over journey entries in v0 and v1.
    // The entries are sorted by OrientedReadId, but an
    // OrientedReadId can appear more than one.
    auto it0 = vertex0.journeyEntries.begin();
    auto it1 = vertex1.journeyEntries.begin();
    const auto end0 = vertex0.journeyEntries.end();
    const auto end1 = vertex1.journeyEntries.end();
    while(it0!=end0 and it1!=end1) {
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
            continue;
        }
        if(it0->orientedReadId > it1->orientedReadId) {
            ++it1;
            continue;
        }
        const OrientedReadId orientedReadId = it0->orientedReadId;
        SHASTA_ASSERT(orientedReadId == it1->orientedReadId);

        // Find the streaks in v0 and v1 for this oriented read.
        // In most cases these streak have length one as each
        // oriented read appears once in the journey entries
        // of each vertex.
        const auto streakBegin0 = it0;
        auto streakEnd0 = streakBegin0;
        while(streakEnd0 != end0 and streakEnd0->orientedReadId == orientedReadId) {
            ++streakEnd0;
        }
        const auto streakBegin1 = it1;
        auto streakEnd1 = streakBegin1;
        while(streakEnd1 != end1 and streakEnd1->orientedReadId == orientedReadId) {
            ++streakEnd1;
        }



        // Loop over entries in [streakBegin0, streakEnd0).
        for(auto jt0=streakBegin0; jt0!=streakEnd0; ++jt0) {
            const uint64_t position0 = jt0->position;

            // Find the best matching journey entry in [streakBegin1, streakEnd1).
            uint64_t bestPosition1 = invalid<uint64_t>;
            for(auto jt1=streakBegin1; jt1!=streakEnd1; ++jt1) {
                const uint64_t position1 = jt1->position;
                if(position1 < position0) {
                    continue;
                }
                bestPosition1 = min(bestPosition1, position1);
            }
            if(bestPosition1 == invalid<uint64_t>) {
                continue;
            }
            const uint64_t position1 = bestPosition1;

            // Consider all journey entries in [position0, position1]
            // for this oriented read.

            // Store the vertices encountered in this journey portion
            const auto journey = journeys[orientedReadId.getValue()];
            for(uint64_t position=position0; position<=position1; position++) {
                const vertex_descriptor v = journey[position];
                if(v != null_vertex()) {
                    verticesEncountered.push_back(v);
                }
            }

            // Also store the transitions.
            for(uint64_t position=position0+1; position<=position1; position++) {
                const vertex_descriptor v0 = journey[position-1];
                const vertex_descriptor v1 = journey[position];
                if(v0 != null_vertex() and v1 != null_vertex()) {
                    transitionsEncountered.push_back(make_pair(v0, v1));
                }
            }
        }

        // Position the iterators at the end of the streaks
        it0 = streakEnd0;
        it1 = streakEnd1;
    }



    // Count how many times we encountered each vertex.
    vector<uint64_t> vertexFrequency;
    deduplicateAndCount(verticesEncountered, vertexFrequency);

    // Locate v0 and v1 in verticesEncountered.
    const auto q0 = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v0);
    SHASTA_ASSERT(q0.first != verticesEncountered.end());
    SHASTA_ASSERT(q0.second - q0.first == 1);
    const uint64_t iv0 = q0.first - verticesEncountered.begin();
    const auto q1 = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v1);
    SHASTA_ASSERT(q1.first != verticesEncountered.end());
    SHASTA_ASSERT(q1.second - q1.first == 1);
    const uint64_t iv1 = q1.first - verticesEncountered.begin();

    // Count how many times we encountered each transition.
    // Keep only the ones that appear at least minLinkCoverage times.
    vector<uint64_t> transitionFrequency;
    deduplicateAndCountWithThreshold(
        transitionsEncountered, transitionFrequency, minLinkCoverage);

    // The transitions we kept define a graph.
    // The vertex stores vertexFrequency.
    // The edge stores transitionFrequency.
    using Graph = SecondaryVerticesGraph;
    Graph graph(
        assemblyGraph,
        verticesEncountered,
        vertexFrequency,
        transitionsEncountered,
        transitionFrequency);

    // Write the graph.
    if(debug and debugOut) {
        const string graphName = "ComputeSecondaryVerticesGraph" +
            vertexStringId(v0) + "_" + vertexStringId(v1);
        graph.write(debugOut, graphName);
    }

    // Figure out if the graph is linear.
    const bool isLinear = graph.isLinear(iv0, iv1);
    if(isLinear) {
        debugOut << "Graph is linear\n";
    } else {
        debugOut << "Graph is not linear\n";
    }

    if(debug and debugOut) {
        graph.handleDottedEdges1(debugOut);
    }

    return isLinear;
}



AssemblyGraph::SecondaryVerticesGraph::SecondaryVerticesGraph(
    const AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered,
    const vector<uint64_t>& vertexFrequency,
    const vector< pair<AssemblyGraph::vertex_descriptor, AssemblyGraph::vertex_descriptor> >&
        transitionsEncountered,
    const vector<uint64_t>& transitionFrequency) :
    SecondaryVerticesGraphBaseClass(verticesEncountered.size()),
    assemblyGraph(assemblyGraph),
    verticesEncountered(verticesEncountered)
{
    SHASTA_ASSERT(vertexFrequency.size() == verticesEncountered.size());
    SHASTA_ASSERT(transitionFrequency.size() == transitionsEncountered.size());

    using Graph = SecondaryVerticesGraph;
    Graph& graph = *this;

    BGL_FORALL_VERTICES(v, graph, Graph) {
        graph[v] = vertexFrequency[v];
    }

    for(uint64_t i=0; i<transitionsEncountered.size(); i++) {
        const auto& p = transitionsEncountered[i];
        array<AssemblyGraph::vertex_descriptor, 2> v = {p.first, p.second};
        array<uint64_t, 2> iv;
        for(uint64_t k=0; k<2; k++) {
            const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v[k]);
            SHASTA_ASSERT(q.first != verticesEncountered.end());
            SHASTA_ASSERT(q.second - q.first == 1);
            iv[k] = q.first - verticesEncountered.begin();
        }
        edge_descriptor e;
        tie(e, ignore) = add_edge(iv[0], iv[1], graph);
        graph[e] = transitionFrequency[i];
    }
}



AssemblyGraph::SecondaryVerticesGraph::SecondaryVerticesGraph(
    const AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered,
    uint64_t n) :
    SecondaryVerticesGraphBaseClass(n),
    assemblyGraph(assemblyGraph),
    verticesEncountered(verticesEncountered)
{
}



void AssemblyGraph::SecondaryVerticesGraph::write(
    ostream& graphOut,
    const string& graphName) const
{
    using Graph = SecondaryVerticesGraph;
    const Graph& graph = *this;

    graphOut << "digraph " << graphName << " {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const AssemblyGraph::vertex_descriptor u = verticesEncountered[v];
        const string vertexStringId = assemblyGraph.vertexStringId(u);
        graphOut << "\"" << vertexStringId << "\" [label=\"" <<
            vertexStringId << "\\n" << graph[v] <<
            "\"];\n";
    }

    // Write the edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        const AssemblyGraph::vertex_descriptor v0 = verticesEncountered[source(e, graph)];
        const AssemblyGraph::vertex_descriptor v1 = verticesEncountered[target(e, graph)];
        graphOut <<
            "\"" << assemblyGraph.vertexStringId(v0) << "\"->\"" <<
            assemblyGraph.vertexStringId(v1) << "\" [label=\"" << graph[e] <<
                "\"";
        if(not assemblyGraph.segmentsAreAdjacent(v0, v1)) {
            graphOut << " style=dashed color=red";
        }
        graphOut << "];\n";
    }
    graphOut << "}\n";

}



bool AssemblyGraph::SecondaryVerticesGraph::isLinear(
    vertex_descriptor v0,
    vertex_descriptor v1
) const
{
    using Graph = SecondaryVerticesGraph;
    const Graph& graph = *this;

    BGL_FORALL_VERTICES(v, graph, Graph) {

        // Check v0.
        if(v == v0) {
            if(in_degree(v, graph) != 0) {
                return false;
            }
            if(out_degree(v, graph) != 1) {
                return false;
            }

        // Check v1.
        } else if(v == v1) {
            if(out_degree(v, graph) != 0) {
                return false;
            }
            if(in_degree(v, graph) != 1) {
                return false;
            }

        // Check the other vertices.
        } else {
            if(out_degree(v, graph) != 1) {
                return false;
            }
            if(in_degree(v, graph) != 1) {
                return false;
            }
        }
    }

    return true;

}



template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
bool AssemblyGraph::SecondaryVerticesGraph::HandleDottedEdges1Callback::operator()(
    const CorrespondenceMap1To2& vertexMap,
    const CorrespondenceMap2To1&) const
{

    if(debugOut) {
        BGL_FORALL_VERTICES_T(v, smallGraph, SecondaryVerticesGraph) {
            const uint64_t iv = get(vertexMap, v);
            const AssemblyGraph::vertex_descriptor u = verticesEncountered[iv];
            debugOut << '(' << v << ", " <<
                assemblyGraph.vertexStringId(u) << ") ";
        }
        debugOut << "\n";

        debugOut << "Edge " <<
            assemblyGraph.vertexStringId(verticesEncountered[get(vertexMap, 0)]) << "->" <<
            assemblyGraph.vertexStringId(verticesEncountered[get(vertexMap, 2)]) << " skips vertex " <<
            assemblyGraph.vertexStringId(verticesEncountered[get(vertexMap, 1)]) << "\n";
    }

    return true;
}



// Handle dotted edges that "skip" a vertex.
void AssemblyGraph::SecondaryVerticesGraph::handleDottedEdges1(ostream& debugOut)
{
    using Graph = SecondaryVerticesGraph;
    Graph& graph = *this;

    // The small graph that we will look for.
    // Edge 0->2 "skips" vertex 1.
    Graph graphSmall(assemblyGraph, verticesEncountered, 3);
    add_edge(0, 1, graphSmall);
    add_edge(1, 2, graphSmall);
    add_edge(0, 2, graphSmall);

    HandleDottedEdges1Callback callback(assemblyGraph, verticesEncountered, graphSmall, graph, debugOut);
    boost::vf2_subgraph_iso(graphSmall, graph, callback);
}



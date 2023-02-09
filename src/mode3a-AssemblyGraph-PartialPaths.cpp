// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "longestPath.hpp"
#include "orderPairs.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
// disjoint_sets.hpp must be included firstdue to problems
// in boost include files.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include "dominatorTree.hpp"

// Standard library.
#include "array.hpp"
#include "fstream.hpp"



void AssemblyGraph::computePartialPaths(
    uint64_t threadCount,
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2,
    uint64_t minLinkCoverage)
{
    SHASTA_ASSERT(threadCount > 0);

    // Store the arguments so the thread function can see them.
    computePartialPathsData.segmentCoverageThreshold1 = segmentCoverageThreshold1;
    computePartialPathsData.segmentCoverageThreshold2 = segmentCoverageThreshold2;
    computePartialPathsData.minLinkCoverage = minLinkCoverage;

    const uint64_t segmentCount = verticesBySegment.size();
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::computePartialPathsThreadFunction, threadCount);
}



void AssemblyGraph::computePartialPathsThreadFunction(uint64_t threadId)
{
    ofstream debugOut;
    debugOut.open("ComputePartialPathsDebug-Thread-" + to_string(threadId));

    const uint64_t segmentCoverageThreshold1 = computePartialPathsData.segmentCoverageThreshold1;
    const uint64_t segmentCoverageThreshold2 = computePartialPathsData.segmentCoverageThreshold2;
    const uint64_t minLinkCoverage = computePartialPathsData.minLinkCoverage;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; segmentId++) {

            // Loop over all vertices (replicas) of this segment.
            for(const vertex_descriptor v: verticesBySegment[segmentId]) {
                if(v != null_vertex()) {
                    computePartialPath2(v,
                        segmentCoverageThreshold1, segmentCoverageThreshold2, minLinkCoverage, debugOut);
                }
            }
        }
    }
}



// Compute partial paths (forward and backward) starting from a given vertex.
// Simple version.
void AssemblyGraph::computePartialPath1(
    vertex_descriptor vStart,
    uint64_t minLinkCoverage,
    ostream& debugOut
    )
{
    AssemblyGraph& assemblyGraph = *this;

    if(debugOut) {
        debugOut << "Following reads at " << vertexStringId(vStart) << "\n";
    }

    // Access the start vertex.
    SHASTA_ASSERT(vStart != null_vertex());
    AssemblyGraphVertex& startVertex = assemblyGraph[vStart];

    // Clear the partial paths.
    vector<vertex_descriptor>& forwardPartialPath = startVertex.forwardPartialPath;
    vector<vertex_descriptor>& backwardPartialPath = startVertex.backwardPartialPath;
    forwardPartialPath.clear();
    backwardPartialPath.clear();

    // The vertices we encounter when following the reads.
    vector<vertex_descriptor> verticesEncountered;

    // The transitions we encounter when following the reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitionsEncountered;

    // Loop over JourneyEntry's of the start vertex.
    for(const JourneyEntry& journeyEntry: startVertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;

        // Store the vertices encountered in the journey of this read.
        const auto journey = journeys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            const vertex_descriptor v = journey[position];
            if(v != null_vertex()) {
                verticesEncountered.push_back(v);
            }
        }

        // Also store the transitions.
        for(uint64_t position=1; position<journey.size(); position++) {
            const vertex_descriptor v0 = journey[position-1];
            const vertex_descriptor v1 = journey[position];
            if(v0 != null_vertex() and v1 != null_vertex()) {
                transitionsEncountered.push_back(make_pair(v0, v1));
            }
        }
    }

    // Count how many times we encountered each vertex.
    vector<uint64_t> vertexFrequency;
    deduplicateAndCount(verticesEncountered, vertexFrequency);

    // Count how many times we encountered each transition.
    // Keep only the ones that appear at least minLinkCoverage times.
    vector<uint64_t> transitionFrequency;
    deduplicateAndCountWithThreshold(
        transitionsEncountered, transitionFrequency, minLinkCoverage);

    if(false) {
        debugOut << "Segments:\n";
        for(uint64_t i=0; i<verticesEncountered.size(); i++) {
            const vertex_descriptor v = verticesEncountered[i];
            debugOut << vertexStringId(v) << " " << vertexFrequency[i] << "\n";
        }

        debugOut << "Transitions:\n";
        for(uint64_t i=0; i<transitionsEncountered.size(); i++) {
            const auto& p = transitionsEncountered[i];
            const vertex_descriptor v0 = p.first;
            const vertex_descriptor v1 = p.second;
            debugOut << vertexStringId(v0) << "->" << vertexStringId(v1) <<
                " " << transitionFrequency[i] << "\n";
        }
    }



    // The transitions we kept define a graph.

    // Starting at the start vertex, follow the linear portion of the graph forward.
    // Stop when we encounter a branch or a vertex seen less than minSegmentCoverage times.
    sort(transitionsEncountered.begin(), transitionsEncountered.end(),  // Not strictly necessary.
        OrderPairsByFirstOnly<vertex_descriptor, vertex_descriptor>());
    vertex_descriptor v = vStart;
    while(true) {
        vector<pair<vertex_descriptor, vertex_descriptor> >::iterator it0, it1;
        tie(it0, it1) = std::equal_range(transitionsEncountered.begin(), transitionsEncountered.end(),
            make_pair(v, null_vertex()),
            OrderPairsByFirstOnly<vertex_descriptor, vertex_descriptor>());
        if( it0 == transitionsEncountered.end() or
            (it1 - it0) != 1) {
            break;
        }
        v = it0->second;
        if(v == vStart) {
            break;
        }
        forwardPartialPath.push_back(v);
    }

    if(debugOut) {
        debugOut << "Forward partial path:";
        for(const vertex_descriptor v: forwardPartialPath) {
            debugOut << " " << vertexStringId(v);
        }
        debugOut << "\n";
    }



    // Starting at the start vertex, follow the linear portion of the graph backward.
    // Stop when we encounter a branch or a vertex seen less than minSegmentCoverage times.
    sort(transitionsEncountered.begin(), transitionsEncountered.end(),
        OrderPairsBySecondOnly<vertex_descriptor, vertex_descriptor>());
    v = vStart;
    while(true) {
        vector<pair<vertex_descriptor, vertex_descriptor> >::iterator it0, it1;
        tie(it0, it1) = std::equal_range(transitionsEncountered.begin(), transitionsEncountered.end(),
            make_pair(null_vertex(), v),
            OrderPairsBySecondOnly<vertex_descriptor, vertex_descriptor>());
        if( it0 == transitionsEncountered.end() or
            (it1 - it0) != 1) {
            break;
        }
        v = it0->first;
        if(v == vStart) {
            break;
        }
        backwardPartialPath.push_back(v);
    }

    if(debugOut) {
        debugOut << "Backward partial path:";
        for(const vertex_descriptor v: backwardPartialPath) {
            debugOut << " " << vertexStringId(v);
        }
        debugOut << "\n";
    }

}




// Compute partial paths (forward and backward) starting from a given vertex.
// More robust version that uses dominator trees.
void AssemblyGraph::computePartialPath2(
    vertex_descriptor vStart,
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2,
    uint64_t minLinkCoverage,
    ostream& debugOut
    )
{
    AssemblyGraph& assemblyGraph = *this;

    if(debugOut) {
        debugOut << "Following reads at " << vertexStringId(vStart) << "\n";
    }

    // Access the start vertex.
    SHASTA_ASSERT(vStart != null_vertex());
    AssemblyGraphVertex& startVertex = assemblyGraph[vStart];

    // Clear the partial paths.
    vector<vertex_descriptor>& forwardPartialPath = startVertex.forwardPartialPath;
    vector<vertex_descriptor>& backwardPartialPath = startVertex.backwardPartialPath;
    forwardPartialPath.clear();
    backwardPartialPath.clear();

    // The vertices we encounter when following the reads.
    vector<vertex_descriptor> verticesEncountered;

    // The transitions we encounter when following the reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitionsEncountered;

    // Loop over JourneyEntry's of the start vertex.
    for(const JourneyEntry& journeyEntry: startVertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;

        // Store the vertices encountered in the journey of this read.
        const auto journey = journeys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            const vertex_descriptor v = journey[position];
            if(v != null_vertex()) {
                verticesEncountered.push_back(v);
            }
        }

        // Also store the transitions.
        for(uint64_t position=1; position<journey.size(); position++) {
            const vertex_descriptor v0 = journey[position-1];
            const vertex_descriptor v1 = journey[position];
            if(v0 != null_vertex() and v1 != null_vertex()) {
                transitionsEncountered.push_back(make_pair(v0, v1));
            }
        }
    }

    // Count how many times we encountered each vertex.
    vector<uint64_t> vertexFrequency;
    deduplicateAndCount(verticesEncountered, vertexFrequency);

    // Count how many times we encountered each transition.
    // Keep only the ones that appear at least minLinkCoverage times.
    vector<uint64_t> transitionFrequency;
    deduplicateAndCountWithThreshold(
        transitionsEncountered, transitionFrequency, minLinkCoverage);

    // Write the graph.
    if(debugOut) {
        debugOut << "digraph Graph_" << startVertex.segmentId << "_" << startVertex.segmentReplicaIndex << " {\n";
        for(uint64_t i=0; i<verticesEncountered.size(); i++) {
            const vertex_descriptor v = verticesEncountered[i];
            debugOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
                vertexStringId(v) << "\\n" << vertexFrequency[i] <<
                "\"];\n";
        }

        for(uint64_t i=0; i<transitionsEncountered.size(); i++) {
            const auto& p = transitionsEncountered[i];
            const vertex_descriptor v0 = p.first;
            const vertex_descriptor v1 = p.second;
            debugOut <<
                "\"" << vertexStringId(v0) << "\"->\"" <<
                vertexStringId(v1) << "\";\n";
        }
        debugOut << "}\n";
    }



    // The transitions we kept define a graph.
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
    Graph graph(verticesEncountered.size());
    for(const auto& p: transitionsEncountered) {
        array<vertex_descriptor, 2> v = {p.first, p.second};
        array<uint64_t, 2> iv;
        for(uint64_t k=0; k<2; k++) {
            const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v[k]);
            SHASTA_ASSERT(q.first != verticesEncountered.end());
            SHASTA_ASSERT(q.second - q.first == 1);
            iv[k] = q.first - verticesEncountered.begin();
        }
        add_edge(iv[0], iv[1], graph);
    }



    // To compute the forward partial path, compute the dominator tree of the graph,
    // with the start vertex as the entrance.
    const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), vStart);
    SHASTA_ASSERT(q.first != verticesEncountered.end());
    SHASTA_ASSERT(q.second - q.first == 1);
    const uint64_t ivStart = q.first - verticesEncountered.begin();
    std::map<uint64_t, uint64_t> predecessorMap;
    shasta::lengauer_tarjan_dominator_tree(
        graph,
        ivStart,
        boost::make_assoc_property_map(predecessorMap));

    // Explicitly construct the forward dominator tree.
    Graph forwardTree(verticesEncountered.size());
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.second;
        const uint64_t iv1 = p.first;
        add_edge(iv0, iv1, forwardTree);
    }



    // To compute the forward partial path, follow the forward dominator tree.
    uint64_t iv = ivStart;
    while(true) {

        // Find the out-vertices and sort them by decreasing vertex frequency.
        vector< pair<uint64_t, uint64_t> > outVertices;
        BGL_FORALL_OUTEDGES(iv, e, forwardTree, Graph) {
            const uint64_t iv1 = target(e, forwardTree);
            outVertices.push_back(make_pair(iv1, vertexFrequency[iv1]));
        }
        sort(outVertices.begin(), outVertices.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

        // If there are no out-vertices, the forward path ends here.
        if(outVertices.empty()) {
            break;
        }

        // If the strongest out-vertex is too weak, the forward path ends here.
        if(outVertices.front().second < segmentCoverageThreshold1) {
            break;
        }

        // If the strongest in-vertex loses too much coverage compared to iv, the backward path ends here.
        const uint64_t coverageLoss =
            (outVertices.front().second >= vertexFrequency[iv]) ? 0 :
            (vertexFrequency[iv] - outVertices.front().second);
        if(coverageLoss > segmentCoverageThreshold2) {
            break;
        }

        // In all other cases, we add the strongest out-vertex to the forward path.
        iv = outVertices.front().first;
        forwardPartialPath.push_back(verticesEncountered[iv]);
    }



    // Write the forward dominator tree.
    if(debugOut) {
        debugOut << "digraph Forward_Tree_" << startVertex.segmentId << "_" << startVertex.segmentReplicaIndex << " {\n";

        // Gather the vertices of the dominator tree.
        vector<vertex_descriptor> dominatorTreeVertices;
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            dominatorTreeVertices.push_back(v0);
            dominatorTreeVertices.push_back(v1);
        }
        deduplicate(dominatorTreeVertices);

        for(const vertex_descriptor v: dominatorTreeVertices) {
            const uint64_t i = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v).first -
                verticesEncountered.begin();
            debugOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
                vertexStringId(v) << "\\n" << vertexFrequency[i] <<
                "\"];\n";

        }
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            debugOut <<
                "\"" << vertexStringId(v0) << "\"->\"" <<
                vertexStringId(v1) << "\";\n";
        }
        debugOut << "}\n";
    }



    // To compute the backward partial path, compute the backward dominator tree of the graph,
    // with the start vertex as the entrance.
    predecessorMap.clear();
    shasta::lengauer_tarjan_dominator_tree(
        boost::make_reverse_graph(graph),
        ivStart,
        boost::make_assoc_property_map(predecessorMap));

    // Explicitly construct the backward dominator tree.
    Graph backwardTree(verticesEncountered.size());
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.first;
        const uint64_t iv1 = p.second;
        add_edge(iv0, iv1, backwardTree);
    }



    // To compute the backward partial path, follow the backward dominator tree.
    iv = ivStart;
    while(true) {

        // Find the in-vertices and sort them by decreasing vertex frequency.
        vector< pair<uint64_t, uint64_t> > inVertices;
        BGL_FORALL_INEDGES(iv, e, backwardTree, Graph) {
            const uint64_t iv1 = source(e, forwardTree);
            inVertices.push_back(make_pair(iv1, vertexFrequency[iv1]));
        }
        sort(inVertices.begin(), inVertices.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

        // If there are no in-vertices, the backward path ends here.
        if(inVertices.empty()) {
            break;
        }

        // If the strongest in-vertex is too weak, the backward path ends here.
        if(inVertices.front().second < segmentCoverageThreshold1) {
            break;
        }

        // If the strongest in-vertex loses too much coverage compared to iv, the backward path ends here.
        const uint64_t coverageLoss =
            (inVertices.front().second >= vertexFrequency[iv]) ? 0 :
            (vertexFrequency[iv] - inVertices.front().second);
        if(coverageLoss > segmentCoverageThreshold2) {
            break;
        }

        // In all other cases, we add the strongest in-vertex to the backward path.
        iv = inVertices.front().first;
        backwardPartialPath.push_back(verticesEncountered[iv]);
    }



    // Write the backward dominator tree.
    if(debugOut) {
        debugOut << "digraph Backward_Tree_" << startVertex.segmentId << "_" << startVertex.segmentReplicaIndex << " {\n";

        // Gather the vertices of the dominator tree.
        vector<vertex_descriptor> dominatorTreeVertices;
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the source vertex and the value is the target vertex.
            const uint64_t iv0 = p.first;
            const uint64_t iv1 = p.second;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            dominatorTreeVertices.push_back(v0);
            dominatorTreeVertices.push_back(v1);
        }
        deduplicate(dominatorTreeVertices);

        for(const vertex_descriptor v: dominatorTreeVertices) {
            const uint64_t i = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v).first -
                verticesEncountered.begin();
            debugOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
                vertexStringId(v) << "\\n" << vertexFrequency[i] <<
                "\"];\n";

        }
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.first;
            const uint64_t iv1 = p.second;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            debugOut <<
                "\"" << vertexStringId(v0) << "\"->\"" <<
                vertexStringId(v1) << "\";\n";
        }
        debugOut << "}\n";
    }

    // Compute the fraction of graph vertices that end up in the partial paths.
    if(debugOut) {
        const uint64_t totalPartialPathLength = forwardPartialPath.size() + backwardPartialPath.size() + 1;
        const double partialPathEfficiency = double(totalPartialPathLength) / double(num_vertices(graph));
        debugOut << "Partial path efficiency for " << vertexStringId(vStart) << " " <<
            partialPathEfficiency << "\n";
    }
}



void AssemblyGraph::writePartialPaths() const
{
    const AssemblyGraph& assemblyGraph = *this;

    // One line for each partial path entry.
    ofstream csv1("PartialPaths1.csv");
    csv1 << "Start,Direction,Position,Vertex\n";

    // One line for each partial path.
    ofstream csv2("PartialPaths2.csv");
    csv1 << "Start,Direction,Vertices\n";

    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];

        csv2 << vertexStringId(v0) << ",Forward,";
        for(uint64_t position=0; position<vertex0.forwardPartialPath.size(); position++) {
            const vertex_descriptor v1 = vertex0.forwardPartialPath[position];
            csv1 << vertexStringId(v0) << ",Forward,";
            csv1 << position << ",";
            csv1 << vertexStringId(v1) << "\n";
            csv2 << vertexStringId(v1) << ",";
        }
        csv2 << "\n";

        csv2 << vertexStringId(v0) << ",Backward,";
        for(uint64_t position=0; position<vertex0.backwardPartialPath.size(); position++) {
            const vertex_descriptor v1 = vertex0.backwardPartialPath[position];
            csv1 << vertexStringId(v0) << ",Backward,";
            csv1 << position << ",";
            csv1 << vertexStringId(v1) << "\n";
            csv2 << vertexStringId(v1) << ",";
        }
        csv2 << "\n";
    }
}



void AssemblyGraph::analyzePartialPaths(uint64_t threadCount)
{
    // CONSTANSTS TO EXPOSE WHEN CODE STABILIZES

    // This controls the length of the initial portion of each partial path
    // used here.
    const uint64_t m = 10;
    const uint64_t minComponentSize = 10;



    const AssemblyGraph& assemblyGraph = *this;
    const uint64_t vertexCount = num_vertices(assemblyGraph);

    // Find forward/backward pairs generated by each partial path.
    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {

        for(uint64_t i=0; i<min(m, assemblyGraph[v0].forwardPartialPath.size()); i++) {
            forwardPairs.push_back(make_pair(v0, assemblyGraph[v0].forwardPartialPath[i]));
        }

        for(uint64_t i=0; i<min(m, assemblyGraph[v0].backwardPartialPath.size()); i++) {
            backwardPairs.push_back(make_pair(assemblyGraph[v0].backwardPartialPath[i], v0));
        }

    }

    // The pairs we found in both directions are the ones we want to keep.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    sort(forwardPairs.begin(), forwardPairs.end());
    sort(backwardPairs.begin(), backwardPairs.end());
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));

    cout << "Partial paths analysis:" << endl;
    cout << "Number of vertices " << vertexCount << endl;
    cout << "Number of forward pairs " << forwardPairs.size() << endl;
    cout << "Number of backward pairs " << backwardPairs.size() << endl;
    cout << "Number of bidirectional pairs " << bidirectionalPairs.size() << endl;


#if 0
    // Write the bidirectional pairs as a Graphviz graph.
    ofstream dot("PartialPaths.dot");
    dot << "digraph PartialPaths {\n";
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        dot << "\"" << vertexStringId(v) << "\";\n";
    }
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        dot << "\"" << vertexStringId(v0) << "\"->";
        dot << "\"" << vertexStringId(v1) << "\";\n";
    }
    dot << "}\n";
#endif

    // The bidirectional pairs define a graph.
    // Each vertex of this graph corresponds to a vertex of the assembly graph.
    // To simplify processing below we want to work with a graph with vertices
    // numbered 0 to vertexCount-1.
    // The data structures below define the correspondence of these vertices with the
    // AssemblyGraph vertices.
    vector<vertex_descriptor>& vertexTable = analyzePartialPathsData.vertexTable;
    vertexTable.clear();
    std::map<vertex_descriptor, uint64_t>& vertexMap = analyzePartialPathsData.vertexMap;
    vertexMap.clear();
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexMap.insert(make_pair(v, vertexIndex++));
    }



    // Find its connected components so we can handle them one at a time.
    vector<uint64_t> rank(vertexCount);
    vector<uint64_t> parent(vertexCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<vertexCount; i++) {
        disjointSets.make_set(i);
    }
    for(const auto& p: bidirectionalPairs) {
        disjointSets.union_set(vertexMap[p.first], vertexMap[p.second]);
    }
    vector<uint64_t> componentTable(vertexCount);
    for(uint64_t i=0; i<vertexCount; i++) {
        componentTable[i] = disjointSets.find_set(i);
    }

    // Gather the vertices in each connected component.
    std::map<uint64_t, vector<uint64_t> > connectedComponentMap;
    for(uint64_t i=0; i<vertexCount; i++) {
        const uint64_t componentId = disjointSets.find_set(i);
        connectedComponentMap[componentId].push_back(i);
    }

    // Create a histogram of connected component sizes.
    {
        vector<uint64_t> histogram;
        for(const auto & p: connectedComponentMap) {
            const vector<uint64_t>& component = p.second;
            const uint64_t componentSize = component.size();
            if(histogram.size() <= componentSize) {
                histogram.resize(componentSize + 1, 0);
            }
            ++histogram[componentSize];
        }
        ofstream csv("AnalyzePartialPaths-ComponentSizeHistogram.csv");
        for(uint64_t i=0; i<histogram.size(); i++) {
            const uint64_t frequency = histogram[i];
            if(frequency) {
                csv << i << "," << frequency << "\n";
            }
        }
    }


    // Store the connected components that are sufficiently large.
    vector< vector<uint64_t> >& components = analyzePartialPathsData.components;
    components.clear();
    for(const auto & p: connectedComponentMap) {
        const vector<uint64_t>& component = p.second;
        if(component.size() >= minComponentSize) {
            components.push_back(component);
        }
    }

    // Update the componentTable to reflect the new numbering.
    fill(componentTable.begin(), componentTable.end(), invalid<uint64_t>);
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        for(const uint64_t i: component) {
            componentTable[i] = componentId;
        }
    }

    // Assign the bidirectional pairs to components.
    vector< vector< pair<vertex_descriptor, vertex_descriptor> > >& componentPairs = analyzePartialPathsData.componentPairs;
    componentPairs.clear();
    componentPairs.resize(components.size());
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        const uint64_t iv0 = vertexMap[v0];
        const uint64_t iv1 = vertexMap[v1];
        const uint64_t componentId = componentTable[iv0];
        SHASTA_ASSERT(componentId == componentTable[iv1]);
        if(componentId != invalid<uint64_t>) {
            componentPairs[componentId].push_back(p);
        }
    }


    // Process the connected components in parallel.
    analyzePartialPathsData.longestPaths.resize(components.size());
    setupLoadBalancing(components.size(), 1);
    runThreads(&AssemblyGraph::analyzePartialPathsThreadFunction, threadCount);


#if 0
    // Create the graph defined by the bidirectional pairs
    // and compute its transitive reduction.
    // This will fail if the graph has cycles.
    using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS>;
    Graph graph(vertexTable.size());
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        add_edge(vertexMap[v0], vertexMap[v1], graph);
    }
    // The transitive reduction could be computed in parallel over connected components.
    transitiveReduction(graph);
    cout << "Number of bidirectional pairs after transitive reduction " <<
        num_edges(graph) << endl;

    // Write out the transitive reduction as a Graphviz graph.
    ofstream dot("PartialPaths.dot");
    dot << "digraph PartialPaths {\n";
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = vertexTable[source(e, graph)];
        const vertex_descriptor v1 = vertexTable[target(e, graph)];
        dot << "\"" << vertexStringId(v0) << "\"->";
        dot << "\"" << vertexStringId(v1) << "\";\n";
    }
    dot << "}\n";
#endif
}



void AssemblyGraph::analyzePartialPathsThreadFunction(uint64_t threadId)
{
    ofstream graphOut("AnalyzePartialPathsThread-" + to_string(threadId) + ".dot");

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            graphOut << "digraph AnalyzePartialPathsComponent_" << i << " {\n";
            analyzePartialPathsComponent(
                analyzePartialPathsData.components[i],
                analyzePartialPathsData.componentPairs[i],
                analyzePartialPathsData.longestPaths[i],
                graphOut);
            graphOut << "}\n";

            // Write out the longest path.
            graphOut << "digraph AnalyzePartialPathsComponentLongestPath_" << i << " {\n";
            for(uint64_t j=1; j<analyzePartialPathsData.longestPaths[i].size(); j++) {
                const vertex_descriptor v0 = analyzePartialPathsData.longestPaths[i][j-1];
                const vertex_descriptor v1 = analyzePartialPathsData.longestPaths[i][j];
                graphOut << "\"" << vertexStringId(v0) << "\"->";
                graphOut << "\"" << vertexStringId(v1) << "\";\n";
            }
            graphOut << "}\n";
        }
    }
}



void AssemblyGraph::analyzePartialPathsComponent(
    const vector<uint64_t>& component,
    const vector< pair<vertex_descriptor, vertex_descriptor> >& componentPairs,
    vector<vertex_descriptor>& longestPath,
    ostream& graphOut)
{
    const vector<vertex_descriptor>& vertexTable = analyzePartialPathsData.vertexTable;

#if 0
    // Write the vertices to graphOut.
    for(const uint64_t iv: component) {
        const vertex_descriptor v = vertexTable[iv];
        graphOut << "\"" << vertexStringId(v) << "\";\n";
    }

    // Write edges to graphOut.
    for(const auto& p: componentPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        graphOut << "\"" << vertexStringId(v0) << "\"->";
        graphOut << "\"" << vertexStringId(v1) << "\";\n";
    }
#endif

    // Number the vertices in this component consecutively starting at zero.
    vector<vertex_descriptor> componentVertexTable;
    std::map<vertex_descriptor, uint64_t> componentVertexMap;
    for(uint64_t jv=0; jv<component.size(); jv++) {
        const uint64_t iv = component[jv];
        const vertex_descriptor v = vertexTable[iv];
        componentVertexTable.push_back(v);
        componentVertexMap.insert(make_pair(v, jv));
    }

    // Explicitly construct the graph for this component.
    using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS>;
    Graph graph(componentVertexTable.size());
    for(const auto& p: componentPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        add_edge(componentVertexMap[v0], componentVertexMap[v1], graph);
    }

    // Since we only cre about the longest path, we don't need to
    // compute the transitive reduction.
    // transitiveReduction(graph);

    // Write the vertices to graphOut.
    BGL_FORALL_VERTICES(iv, graph, Graph) {
        const vertex_descriptor v = componentVertexTable[iv];
        graphOut << "\"" << vertexStringId(v) << "\";\n";
    }

    // Write edges to graphOut.
    BGL_FORALL_EDGES(e, graph, Graph) {
        const uint64_t jv0 = source(e, graph);
        const uint64_t jv1 = target(e, graph);
        const vertex_descriptor v0 = componentVertexTable[jv0];
        const vertex_descriptor v1 = componentVertexTable[jv1];
        graphOut << "\"" << vertexStringId(v0) << "\"->";
        graphOut << "\"" << vertexStringId(v1) << "\";\n";
    }

    // Compute the longest path.
    vector<uint64_t> componentLongestPath;
    shasta::longestPath(graph, componentLongestPath);
    longestPath.clear();
    for(const uint64_t iv: componentLongestPath) {
        longestPath.push_back(componentVertexTable[iv]);
    }
}



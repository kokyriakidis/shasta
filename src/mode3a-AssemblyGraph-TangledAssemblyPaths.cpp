// Shasta.
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
#include "enumeratePaths.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include "dominatorTree.hpp"

// Standard library.
#include "array.hpp"
#include "fstream.hpp"



void AssemblyGraph::computeAssemblyPaths(uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;
    cout << "computeAssemblyPaths begins." << endl;

    // Initialize a TangledAssemblyPath for each of the
    // longest path computed by analyzePartialPaths.
    assemblyPaths.clear();
    assemblyPaths.resize(analyzePartialPathsData.longestPaths.size());

    // Compute the AssemblyPaths in parallel.
    setupLoadBalancing(assemblyPaths.size(), 1);
    runThreads(&AssemblyGraph::computeAssemblyPathsThreadFunction, threadCount);

    // Sort them by decreasing efficiency.
    sort(assemblyPaths.begin(), assemblyPaths.end(), OrderAssemblyPath());



    // Store path information in vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        assemblyGraph[v].pathInformation.clear();
    }
    for(uint64_t pathId=0; pathId<assemblyPaths.size(); pathId++) {
        const AssemblyPath& path = *assemblyPaths[pathId];
        SHASTA_ASSERT(path.secondaryVerticesInfos.size() == path.primaryVertices.size() - 1);

        // Information for primary vertices of this path.
        for(uint64_t positionInPath=0; positionInPath<path.primaryVertices.size(); positionInPath++) {
            const vertex_descriptor v = path.primaryVertices[positionInPath];
            auto& primaryInfo = assemblyGraph[v].pathInformation.primaryInfo;

            // Each vertex can only appear once as a primary vertex in a path.
            SHASTA_ASSERT(primaryInfo.pathId == invalid<uint64_t>);

            primaryInfo = {pathId, positionInPath};
        }

        // Information for secondary vertices of this path.
        for(uint64_t positionInPath=0; positionInPath<path.secondaryVerticesInfos.size(); positionInPath++) {
            const auto& secondaryVertexInfo = path.secondaryVerticesInfos[positionInPath];

            for(uint64_t positionInLeg=0; positionInLeg<secondaryVertexInfo.secondaryVertices.size(); positionInLeg++) {
                const vertex_descriptor v = secondaryVertexInfo.secondaryVertices[positionInLeg];
                assemblyGraph[v].pathInformation.secondaryInfos.push_back({pathId, positionInPath, positionInLeg});
            }
        }
    }



    writeAssemblyPaths1();
    writeAssemblyPaths2();
    writeAssemblyPathsVertexSummary();
    writeAssemblyPathsVertexInfo();
    writeAssemblyPathsVertexHistogram();
    writeAssemblyPathsJourneyIntervals();
    writeAssemblyPathsJourneyInfo();
}



void AssemblyGraph::writeAssemblyPaths1() const
{

    ofstream csv(debugOutputPrefix + "AssemblyPaths1.csv");
    csv << "Path,Path efficiency,Position,v0,v1,Efficiency,\n";

    for(uint64_t pathId=0; pathId<assemblyPaths.size(); pathId++) {
        const AssemblyPath& path = *assemblyPaths[pathId];
        SHASTA_ASSERT(path.secondaryVerticesInfos.size() == path.primaryVertices.size() - 1);

        for(uint64_t position=0; position<path.secondaryVerticesInfos.size(); position++) {
            const auto& secondaryVertexInfo = path.secondaryVerticesInfos[position];
            const vertex_descriptor v0 = path.primaryVertices[position];
            const vertex_descriptor v1 = path.primaryVertices[position+1];
            csv << pathId << ",";
            csv << path.efficiency << ",";
            csv << position << ",";
            csv << vertexStringId(v0) << ",";
            csv << vertexStringId(v1) << ",";
            csv << secondaryVertexInfo.efficiency << ",";

            for(const vertex_descriptor v: secondaryVertexInfo.secondaryVertices) {
                csv << vertexStringId(v) << ",";
            }
            csv << "\n";
        }
    }
}



void AssemblyGraph::writeAssemblyPaths2() const
{

    ofstream csv(debugOutputPrefix + "AssemblyPaths2.csv");
    csv << "Path,Vertices\n";

    for(uint64_t pathId=0; pathId<assemblyPaths.size(); pathId++) {
        const AssemblyPath& path = *assemblyPaths[pathId];
        SHASTA_ASSERT(path.secondaryVerticesInfos.size() == path.primaryVertices.size() - 1);

        csv << pathId << ",";
        for(uint64_t position=0; /* Check later */ ; position++) {
            const auto& secondaryVertexInfo = path.secondaryVerticesInfos[position];
            const vertex_descriptor vPrimary = path.primaryVertices[position];
            csv << vertexStringId(vPrimary) << ",";

            // If this is the last primary vertex, there are no secondary vertices.
            if(position == path.secondaryVerticesInfos.size()) {
                break;
            }

            // Write the secondary vertices between this primary vertex and the next.
            for(const vertex_descriptor vSecondary: secondaryVertexInfo.secondaryVertices) {
                csv << vertexStringId(vSecondary) << ",";
            }
        }
        csv << "\n";
    }
}



void AssemblyGraph::writeAssemblyPathsVertexSummary() const
{
    const AssemblyGraph& assemblyGraph = *this;
    ofstream csv(debugOutputPrefix + "AssemblyPathsVertexSummary.csv");
    csv << "Vertex,PrimaryCount,SecondaryCount\n";

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        csv << vertex.stringId() << ",";
        csv << vertex.pathInformation.primaryCount() << ",";
        csv << vertex.pathInformation.secondaryInfos.size() << "\n";
    }
}



void AssemblyGraph::writeAssemblyPathsVertexInfo() const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(debugOutputPrefix + "AssemblyPathsVertexInfo.csv");
    csv << "Vertex,Type,PathId,PositionInPath,PositionInLeg\n";

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const auto& primaryInfo = vertex.pathInformation.primaryInfo;

        if(primaryInfo.pathId != invalid<uint64_t>) {
            csv << vertex.stringId() << ",";
            csv << "Primary,";
            csv << primaryInfo.pathId << ",";
            csv << primaryInfo.positionInPath << "\n";
        }

        for(const auto& info: vertex.pathInformation.secondaryInfos) {
            csv << vertex.stringId() << ",";
            csv << "Secondary,";
            csv << info.pathId << ",";
            csv << info.positionInPath << ",";
            csv << info.positionInLeg << "\n";
        }
    }
}



// Histogram of the number of vertices by number of primary/secondary paths.
void AssemblyGraph::writeAssemblyPathsVertexHistogram() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< vector<uint64_t> > histogram;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const uint64_t primaryCount = vertex.pathInformation.primaryCount();
        const uint64_t secondaryCount = vertex.pathInformation.secondaryInfos.size();
        if(histogram.size() <= primaryCount) {
            histogram.resize(primaryCount + 1);
        }
        vector<uint64_t>& histogramRow = histogram[primaryCount];
        if(histogramRow.size() <= secondaryCount) {
            histogramRow.resize(secondaryCount + 1);
        }
        ++histogramRow[secondaryCount];
    }

    cout << "Histogram of number of vertices by number of appearances in paths." << endl;
    cout << "PrimaryCount,SecondaryCount,VertexCount" << endl;
    for(uint64_t primaryCount=0; primaryCount<histogram.size(); primaryCount++) {
        const vector<uint64_t>& histogramRow = histogram[primaryCount];
        for(uint64_t secondaryCount=0; secondaryCount<histogramRow.size(); secondaryCount++) {
            const uint64_t vertexCount = histogramRow[secondaryCount];
            cout << primaryCount << ",";
            cout << secondaryCount << ",";
            cout << vertexCount << endl;
        }
    }
}



void AssemblyGraph::writeAssemblyPathsJourneyInfo() const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(debugOutputPrefix + "TangledAssemblyPathsJourneyInfo.csv");

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const auto& primaryInfo = vertex.pathInformation.primaryInfo;

        // For each of the paths that this vertex is a secondary vertex of,
        // find the journey entries that belong to that path.
        vector< vector<bool> > journeyEntryFlags(vertex.pathInformation.secondaryInfos.size());
        for(uint64_t i=0; i<journeyEntryFlags.size(); i++) {
            const auto& secondaryInfo = vertex.pathInformation.secondaryInfos[i];
            const AssemblyPath& path = *assemblyPaths[secondaryInfo.pathId];
            path.secondaryVerticesInfos[secondaryInfo.positionInPath].getVertexJourneyEntries(
                vertex, journeyEntryFlags[i]);
        }


        csv << "Vertex," << vertexStringId(v) << "\n";

        // Write the header for this vertex.
        csv << "OrientedReadId,PositionInJourney,";
        if(vertex.pathInformation.primaryCount() == 1) {
            csv << "P-";
            csv << primaryInfo.pathId << "-";
            csv << primaryInfo.positionInPath << ",";

        }
        for(const auto& info: vertex.pathInformation.secondaryInfos) {
            csv << "S-";
            csv << info.pathId << "-";
            csv << info.positionInPath << "-";
            csv << info.positionInLeg << ",";
        }
        csv << "\n";

        // Write one row for each journey entry in this vertex.
        for(uint64_t i=0; i<vertex.journeyEntries.size(); i++) {
            const JourneyEntry& journeyEntry = vertex.journeyEntries[i];
            csv << journeyEntry.orientedReadId << ",";
            csv << journeyEntry.position << ",";

            // If this is a primary vertex of a path, all journey entries are in the path.
            if(vertex.pathInformation.primaryCount() == 1) {
                csv << "1,";
            }

            // Write the journey entry flags for the paths that his vertex is a secondary vertex of.
            for(uint64_t j=0; j<journeyEntryFlags.size(); j++) {
                csv << int(journeyEntryFlags[j][i]) << ",";
            }

            csv << "\n";
        }
     }

}



void AssemblyGraph::writeAssemblyPathsJourneyIntervals() const
{
    ofstream csv(debugOutputPrefix + "AssemblyPathsJourneyIntervals.csv");
    csv << "PathId,Position,OrientedReadId,Begin,End\n";

    for(uint64_t pathId=0; pathId<assemblyPaths.size(); pathId++) {
        const AssemblyPath& path = *assemblyPaths[pathId];
        SHASTA_ASSERT(path.secondaryVerticesInfos.size() == path.primaryVertices.size() - 1);

        for(uint64_t position=0; position<path.secondaryVerticesInfos.size(); position++) {
            const auto& secondaryVertexInfo = path.secondaryVerticesInfos[position];
            for(const auto& journeyInterval: secondaryVertexInfo.journeyIntervals) {
                csv << pathId << ",";
                csv << position << ",";
                csv << journeyInterval.orientedReadId << ",";
                csv << journeyInterval.begin << ",";
                csv << journeyInterval.end << "\n";
            }
        }
    }

}



// Given a vertex, find which journey entries in the vertex
// are in one of the journey intervals for this SecondaryVertexInfo.
void AssemblyGraph::AssemblyPath::SecondaryVertexInfo::getVertexJourneyEntries(
    const AssemblyGraphVertex& vertex,
    vector<bool>& flags   // True of false for each of the journey entries in the vertex.
) const
{
    const bool debug = false; // (vertex.segmentId == 14598);
    if(debug) {
        cout << " getVertexJourneys called for vertex " << vertex.stringId() << "\n";
    }

    const vector<JourneyEntry>& journeyEntries = vertex.journeyEntries;
    flags.clear();
    flags.resize(journeyEntries.size(), false);

    // Joint loop over the journey entries of the vertex and the
    // journeyIntervals of the SecondaryVertexInfo.
    // Iterators on journeyEntries are prefixed with "it".
    // Iterators on journeyIntervals are prefixed with "jt".
    const auto itBegin = journeyEntries.begin();
    const auto itEnd = journeyEntries.end();
    auto it = itBegin;
    auto jt = journeyIntervals.begin();
    const auto jtEnd = journeyIntervals.end();

    while(it!=itEnd and jt!=jtEnd) {
        if(it->orientedReadId < jt->orientedReadId) {
            ++it;
            continue;
        }
        if(jt->orientedReadId < it->orientedReadId) {
            ++jt;
            continue;
        }
        const OrientedReadId orientedReadId = it->orientedReadId;
        SHASTA_ASSERT(orientedReadId == jt->orientedReadId);
        if(debug) {
            cout << "Found common oriented read id " << orientedReadId << endl;
        }

        // journeyEntries and journeyIntervals can both contain more than one
        // entry for this OrientedReadId, although they will in most cases
        // contain only one.
        // Find the streak for this OrientedReadId in both journeyEntries and journeyIntervals.
        auto itStreakBegin = it;
        auto itStreakEnd = itStreakBegin;
        while(itStreakEnd != itEnd and itStreakEnd->orientedReadId == orientedReadId) {
            ++itStreakEnd;
        }
        auto jtStreakBegin = jt;
        auto jtStreakEnd = jtStreakBegin;
        while(jtStreakEnd != jtEnd and jtStreakEnd->orientedReadId == orientedReadId) {
            ++jtStreakEnd;
        }



        // Loop over journeyEntries in this streak.
        for(auto itStreak=itStreakBegin; itStreak!=itStreakEnd; ++itStreak) {
            const uint64_t entryPosition = itStreak->position;

            // See if we have a journey interval containing this position.
            bool found = false;
            for(auto jtStreak=jtStreakBegin; jtStreak!=jtStreakEnd; ++jtStreak) {
                if(entryPosition >= jtStreak->begin and entryPosition < jtStreak->end) {
                    found = true;
                    break;
                }

            }
            if(found) {
                flags[itStreak - itBegin] = true;
            }

        }



        // Prepare for the next loop iteration.
        it = itStreakEnd;
        jt = jtStreakEnd;
    }
}



void AssemblyGraph::computeAssemblyPathsThreadFunction(uint64_t threadId)
{
    ofstream debugOut(debugOutputPrefix + "ComputeAssemblyPaths-Thread" + to_string(threadId) + ".txt");

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; ++i) {
            assemblyPaths[i] = make_shared<AssemblyPath>();
            computeAssemblyPath(
                analyzePartialPathsData.longestPaths[i],
                *assemblyPaths[i],
                debugOut);
        }
    }

}



void AssemblyGraph::computeAssemblyPath(
    const vector<vertex_descriptor>& primaryVertices,
    AssemblyPath& assemblyPath,
    ostream& debugOut
    )
{
    // Store the primary vertices.
    assemblyPath.primaryVertices = primaryVertices;

    if(debugOut) {
        debugOut << "Computing secondary vertices for tangled assembly path " <<
            vertexStringId(primaryVertices.front()) << "..." <<
            vertexStringId(primaryVertices.back()) <<
            " with primary vertices:\n";
        for(uint64_t i=0; i<primaryVertices.size() ; i++) {
            debugOut << vertexStringId(primaryVertices[i]) << " ";
        }
        debugOut << "\n";
    }

    // Compute the secondary vertices for each pair of primary vertices.
    assemblyPath.secondaryVerticesInfos.clear();
    assemblyPath.secondaryVerticesInfos.resize(assemblyPath.primaryVertices.size() - 1);

    for(uint64_t i=0; i<assemblyPath.secondaryVerticesInfos.size(); i++) {
        computeSecondaryVertices(
            assemblyPath.primaryVertices[i],
            assemblyPath.primaryVertices[i+1],
            assemblyPath.secondaryVerticesInfos[i],
            debugOut);
    }
    assemblyPath.computeEfficiency();
    if(debugOut) {
        debugOut << "Assembly path " <<
            vertexStringId(primaryVertices.front()) << "..." <<
            vertexStringId(primaryVertices.back()) <<
            " with efficiency " << assemblyPath.efficiency <<
            "\n";
        for(uint64_t i=0; /* Check later */ ; i++) {

            debugOut << "P" << vertexStringId(primaryVertices[i]) << " ";
            if(i == primaryVertices.size() - 1) {
                break;
            }

            for(const vertex_descriptor v: assemblyPath.secondaryVerticesInfos[i].secondaryVertices) {
                debugOut << vertexStringId(v) << " ";
            }
        }
        debugOut << "\n";
    }

}



// Given a pair of primary vertices in an assembly path,
// compute the intervening secondary vertices.
// The code is similar to AssemblyGraph::computePartialPath2,
// but instead of using entire oriented read journeys,
// it only uses portions between v0 and v1.
void AssemblyGraph::computeSecondaryVertices(
    vertex_descriptor v0,
    vertex_descriptor v1,
    AssemblyPath::SecondaryVertexInfo& secondaryVerticesInfo,
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

    // The vertices we encounter when following the oriented read journeys.
    vector<vertex_descriptor> verticesEncountered;

    // The transitions we encounter when following the oriented read journeys.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitionsEncountered;

    // Joint loop over journey entries in v0 and v1.
    // The entries are sorted by OrientedReadId, but an
    // OrientedReadId can appear more than once.
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
            secondaryVerticesInfo.journeyIntervals.push_back(
                {orientedReadId, position0, position1+1});
            if(debug and debugOut) {
                debugOut << "Added " << orientedReadId << " " << position0 << " " << position1+1 << "\n";
            }

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
        transitionFrequency,
        iv0, iv1);

    // Write the graph.
    if(debug and debugOut) {
        const string graphName = "ComputeSecondaryVerticesGraph" +
            vertexStringId(v0) + "_" + vertexStringId(v1);
        graph.write(debugOut, graphName);
    }

    // Compute the path v0...v1 on the dominator tree.
    if(not graph.computeDominatorTreePath()) {
        secondaryVerticesInfo.secondaryVertices.clear();
        secondaryVerticesInfo.failed = true;
    }
    if(debug and debugOut) {
        debugOut << "Dominator tree path:";
        for(const Graph::vertex_descriptor v: graph.dominatorTreePath) {
            debugOut << " " << graph.vertexStringId(v);
        }
        debugOut << "\n";
    }
    graph.computeBestPath(debugOut);


    // Use the best path stored in the graph to fill in the secondary
    // vertices between this pair of primary vertices.
    secondaryVerticesInfo.secondaryVertices.clear();
    for(uint64_t i=1; i<graph.bestPath.size(); i++) {
        const Graph::edge_descriptor e = graph.bestPath[i];
        const Graph::vertex_descriptor v = source(e, graph);
        const vertex_descriptor u = verticesEncountered[v];
        secondaryVerticesInfo.secondaryVertices.push_back(u);
    }



    // Compute the efficiency of the secondary vertices as ratio of total vertex coverage
    // on the best path (including the primary vertices) over total vertex coverage on the entire graph.
    uint64_t sum1 = vertexFrequency[iv0] + vertexFrequency[iv1];
    for(uint64_t i=1; i<graph.bestPath.size(); i++) {
        const Graph::edge_descriptor e = graph.bestPath[i];
        const Graph::vertex_descriptor v = source(e, graph);
        sum1 += vertexFrequency[v];
    }
    uint64_t sum2 = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        sum2 += vertexFrequency[v];
    }
    secondaryVerticesInfo.efficiency = double(sum1) / double(sum2);
    if(debug and debugOut) {
        debugOut << "Secondary vertices efficiency for " <<
            vertexStringId(v0) << " " << vertexStringId(v1) << " " <<
            secondaryVerticesInfo.efficiency <<
            " " << sum1 << " " << sum2 <<
            "\n";
    }

#if 0
    // Remove edges that "skip" a vertex. These are commonly generated by errors.
    if(debug and debugOut) {
        graph.handleDottedEdges1(debugOut);
    }

    // Figure out if the graph is linear.
    const bool isLinear = graph.isLinear(iv0, iv1);
    if(isLinear) {
        debugOut << "Graph is linear\n";
    } else {
        debugOut << "Graph is not linear\n";
    }

    return isLinear;
#endif
}



AssemblyGraph::SecondaryVerticesGraph::SecondaryVerticesGraph(
    const AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraph::vertex_descriptor>& verticesEncountered,
    const vector<uint64_t>& vertexFrequency,
    const vector< pair<AssemblyGraph::vertex_descriptor, AssemblyGraph::vertex_descriptor> >&
        transitionsEncountered,
    const vector<uint64_t>& transitionFrequency,
    vertex_descriptor iv0,
    vertex_descriptor iv1) :
    SecondaryVerticesGraphBaseClass(verticesEncountered.size()),
    assemblyGraph(assemblyGraph),
    verticesEncountered(verticesEncountered),
    vertexFrequency(vertexFrequency),
    iv0(iv0),
    iv1(iv1)
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
    const vector<uint64_t>& vertexFrequency,
    uint64_t n) :
    SecondaryVerticesGraphBaseClass(n),
    assemblyGraph(assemblyGraph),
    verticesEncountered(verticesEncountered),
    vertexFrequency(vertexFrequency)
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
        graphOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
            vertexStringId(v) << "\\n" << graph[v] <<
            "\"];\n";
    }

    // Write the edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AssemblyGraph::vertex_descriptor u0 = verticesEncountered[v0];
        const AssemblyGraph::vertex_descriptor u1 = verticesEncountered[v1];
        graphOut <<
            "\"" << vertexStringId(v0) << "\"->\"" <<
            vertexStringId(v1) << "\" [label=\"" << graph[e] <<
                "\"";
        if(not assemblyGraph.segmentsAreAdjacent(u0, u1)) {
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
        const uint64_t inDegree = in_degree(v, graph);
        const uint64_t outDegree = out_degree(v, graph);

        // Check v0.
        if(v == v0) {
            if(inDegree != 0) {
                return false;
            }
            if(outDegree != 1) {
                return false;
            }

        // Check v1.
        } else if(v == v1) {
            if(outDegree != 0) {
                return false;
            }
            if(inDegree != 1) {
                return false;
            }

        // Check the other vertices.
        } else {

            // Allow isolated vertex.
            if(inDegree==0 and outDegree==0) {
                // This vertex is isolated.
                // Don't check anything else.
            } else {
                if(outDegree != 1) {
                    return false;
                }
                if(inDegree != 1) {
                    return false;
                }
            }
        }
    }

    return true;

}



// Remove edges that "skip" a vertex. These are commonly generated by errors.
template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
bool AssemblyGraph::SecondaryVerticesGraph::HandleDottedEdges1Callback::operator()(
    const CorrespondenceMap1To2& vertexMap,
    const CorrespondenceMap2To1&)
{

    if(debugOut) {
        BGL_FORALL_VERTICES_T(v, smallGraph, SecondaryVerticesGraph) {
            const uint64_t iv = get(vertexMap, v);
            const AssemblyGraph::vertex_descriptor u = verticesEncountered[iv];
            debugOut << '(' << v << ", " <<
                assemblyGraph.vertexStringId(u) << ") ";
        }
        debugOut << "\n";

        debugOut <<
            "v0  = " <<
            assemblyGraph.vertexStringId(verticesEncountered[get(vertexMap, 0)]) <<
            ", v1  = " <<
            assemblyGraph.vertexStringId(verticesEncountered[get(vertexMap, 1)]) <<
            ", v2  = " <<
            assemblyGraph.vertexStringId(verticesEncountered[get(vertexMap, 2)]) << "\n";
    }



    // The small graph consists of edges 0->1, 1->2, and 0->2,
    // that is, edge 0->2 "skips" vertex 1.
    // Find the corresponding vertices and edges in the SecondaryVerticesGraph.
    // such that the following edges exist:
    const vertex_descriptor v0 = get(vertexMap, 0);
    const vertex_descriptor v1 = get(vertexMap, 1);
    const vertex_descriptor v2 = get(vertexMap, 2);
    edge_descriptor e01;
    edge_descriptor e12;
    edge_descriptor e02;
    tie(e01, ignore) = boost::edge(v0, v1, secondaryVerticesGraph);
    tie(e12, ignore) = boost::edge(v1, v2, secondaryVerticesGraph);
    tie(e02, ignore) = boost::edge(v0, v2, secondaryVerticesGraph);

    // Find coverage in the SecondaryVerticesGraph for these vertices and edges.
    const uint64_t c0 = vertexFrequency[v0];
    const uint64_t c1 = vertexFrequency[v1];
    const uint64_t c2 = vertexFrequency[v2];
    const uint64_t c01 = secondaryVerticesGraph[e01];
    const uint64_t c12 = secondaryVerticesGraph[e12];
    const uint64_t c02 = secondaryVerticesGraph[e02];

    if(debugOut) {
        debugOut << "c0 " << c0 << "\n";
        debugOut << "c1 " << c1 << "\n";
        debugOut << "c2 " << c2 << "\n";
        debugOut << "c01 " << c01 << "\n";
        debugOut << "c12 " << c12 << "\n";
        debugOut << "c02 " << c02 << "\n";
    }

    if(
        c02<c0 and c02<c1 and c02<c2 and
        c02<c01 and c02<c12) {
        edgesToBeRemoved.push_back(e02);
    }

    // Return true to continue the processing for other triangles in the
    // SecondaryVerticesGraph.
    return true;
}



// Remove edges that "skip" a vertex. These are commonly generated by errors.
void AssemblyGraph::SecondaryVerticesGraph::handleDottedEdges1(ostream& debugOut)
{
    using Graph = SecondaryVerticesGraph;
    Graph& graph = *this;

    // The small graph that we will look for.
    // Edge 0->2 "skips" vertex 1.
    Graph graphSmall(assemblyGraph, verticesEncountered, vertexFrequency, 3);
    add_edge(0, 1, graphSmall);
    add_edge(1, 2, graphSmall);
    add_edge(0, 2, graphSmall);

    vector<edge_descriptor> edgesToBeRemoved;
    HandleDottedEdges1Callback callback(
        assemblyGraph,
        verticesEncountered, vertexFrequency,
        graphSmall, graph, edgesToBeRemoved, debugOut);
    boost::vf2_subgraph_iso(graphSmall, graph, callback);

    // Remove the edges.
    deduplicate(callback.edgesToBeRemoved);
    for(const edge_descriptor e: callback.edgesToBeRemoved) {
        if(debugOut) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            debugOut << "Removing edge " <<
                assemblyGraph.vertexStringId(verticesEncountered[v0]) << "->" <<
                assemblyGraph.vertexStringId(verticesEncountered[v1]) << "\n";
        }
        boost::remove_edge(e, graph);
    }
}



bool AssemblyGraph::SecondaryVerticesGraph::computeDominatorTreePath()
{
    SecondaryVerticesGraph& graph = *this;

    // Compute the dominator tree with v0 as the entrance.
    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
    shasta::lengauer_tarjan_dominator_tree(
        graph,
        iv0,
        boost::make_assoc_property_map(predecessorMap));

    // Construct the dominator tree path starting at v1 and moving back.
    SHASTA_ASSERT(iv0 != iv1);
    dominatorTreePath.clear();
    dominatorTreePath.push_back(iv1);
    vertex_descriptor iv = iv1;
    while(true) {
        auto it = predecessorMap.find(iv);
        if(it == predecessorMap.end()) {
            return false;   // v1 is not reachable from v0.
        }
        iv = it->second;
        dominatorTreePath.push_back(iv);
        if(iv == iv0) {
            break;
        }
    }
    std::reverse(dominatorTreePath.begin(), dominatorTreePath.end());

    return true;
}



// Compute the "best" path between iv0 and iv1.
void AssemblyGraph::SecondaryVerticesGraph::computeBestPath(ostream& debugOut)
{
    SecondaryVerticesGraph& graph = *this;
    bestPath.clear();

    // Loop over legs of the dominator tree path.
    for(uint64_t i=1; i<dominatorTreePath.size(); i++) {
        const vertex_descriptor u0 = dominatorTreePath[i-1];
        const vertex_descriptor u1 = dominatorTreePath[i];

        // We need to find the "best" path between u0 and u1.


        // Enumerate paths u0->...->u1 on the graph.
        // I was not able to get path enumeration to work on a filtered
        // graph.
        vector< vector<edge_descriptor> > allPaths;
        enumerateSelfAvoidingPaths(graph, u0, u1, allPaths);
        if(false) {
            debugOut << "Looking for paths " <<
                vertexStringId(u0) << "->...->" << vertexStringId(u1) << "\n";
            debugOut << "Found " << allPaths.size() << " paths.\n";
        }

        // Find the paths that only use solid edges.
        vector< vector<edge_descriptor> > solidPaths;
        for(const auto& path: allPaths) {
            bool hasDottedEdges = false;
            for(const edge_descriptor e: path) {
                if(not segmentsAreAdjacent(e)) {
                    hasDottedEdges = true;
                    break;
                }
            }
            if(not hasDottedEdges) {
                solidPaths.push_back(path);
            }
        }
        if(false) {
            debugOut << solidPaths.size() << " paths use only solid edges.\n";
        }

        // If paths that only uses solid edges are present, choose among those.
        const vector< vector<edge_descriptor> >& pathsToChooseFrom =
            solidPaths.empty() ? allPaths : solidPaths;

        // We hope this never happens but this is not guaranteed,
        // so we will have to behave better eventually.
        SHASTA_ASSERT(not pathsToChooseFrom.empty());

        // If there is only one path to choose from, pick that one and we are done.
        if(pathsToChooseFrom.size() == 1) {
            copy(pathsToChooseFrom.front().begin(), pathsToChooseFrom.front().end(),
                back_inserter(bestPath));
            continue;
        }



        // There is more than one path to choose from for this leg of the dominator tree.
        // Find the "best" one. Pick the one with the highest minimum edge coverage.
        const vector<edge_descriptor>* bestLegPath = 0;
        uint64_t highestMinimumEdgeCoverage = 0;
        for(const vector<edge_descriptor>& path: pathsToChooseFrom) {

            // Find minimum edge coverage for this one.
            uint64_t minimumEdgeCoverage = std::numeric_limits<uint64_t>::max();
            for(const edge_descriptor e: path) {
                minimumEdgeCoverage = min(minimumEdgeCoverage, graph[e]);
            }
            if(minimumEdgeCoverage > highestMinimumEdgeCoverage) {
                highestMinimumEdgeCoverage = minimumEdgeCoverage;
                bestLegPath = &path;
            }

        }
        SHASTA_ASSERT(bestLegPath);
        copy(bestLegPath->begin(), bestLegPath->end(),
            back_inserter(bestPath));
    }

    // Write out the best path.
    if(debugOut) {
        debugOut << "Best path:";
        for(uint64_t i=0; i<bestPath.size(); i++) {
            const edge_descriptor e = bestPath[i];
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);

            if(i == 0) {
                debugOut << " " << vertexStringId(v0);
            }
            debugOut << " " << vertexStringId(v1);
        }
        debugOut << "\n";
    }
}



AssemblyGraph::AssemblyGraph(
    DetangleUsingTangledAssemblyPaths,
    const AssemblyGraph& oldAssemblyGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    packedMarkerGraph(oldAssemblyGraph.packedMarkerGraph)
{
    createFromAssemblyPaths(oldAssemblyGraph);
}



void AssemblyGraph::createFromAssemblyPaths(
    const AssemblyGraph& oldAssemblyGraph)
{
    const bool debug = false;
    AssemblyGraph& newAssemblyGraph = *this;

    // Initialize verticesBySegment.
    verticesBySegment.clear();
    verticesBySegment.resize(packedMarkerGraph.segments.size());

    // Initialize  oriented reads journeys.
    journeys.clear();
    journeys.resize(packedMarkerGraph.journeys.size());
    for(uint64_t i=0; i<journeys.size(); i++) {
        const auto packedMarkerGraphJourney = packedMarkerGraph.journeys[i];
        journeys[i].resize(packedMarkerGraphJourney.size(), null_vertex());
    }


    // Loop over vertices of the old AssemblyGraph.
    BGL_FORALL_VERTICES(vOld, oldAssemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& oldVertex = oldAssemblyGraph[vOld];
        const uint64_t segmentId = oldVertex.segmentId;

        // For each of the paths that this vertex is a secondary vertex of,
        // find the journey entries that belong to that path.
        vector< vector<bool> > journeyEntryFlags(oldVertex.pathInformation.secondaryInfos.size());
        for(uint64_t i=0; i<journeyEntryFlags.size(); i++) {
            const auto& secondaryInfo = oldVertex.pathInformation.secondaryInfos[i];
            const AssemblyPath& path = *(oldAssemblyGraph.assemblyPaths[secondaryInfo.pathId]);
            path.secondaryVerticesInfos[secondaryInfo.positionInPath].getVertexJourneyEntries(
                oldVertex, journeyEntryFlags[i]);
        }

        // Group the journey entries with identical journeyEntryFlags.
        // For each of the journey entries we construct a binary code
        // based on its journeyEntryFlags.
        SHASTA_ASSERT(journeyEntryFlags.size() < 64);
        std::map<uint64_t, vector<uint64_t> > journeyEntryMap;
        for(uint64_t i=0; i<oldVertex.journeyEntries.size(); i++) {
            uint64_t binaryCode = 0;
            for(uint64_t j=0; j<journeyEntryFlags.size(); j++) {
                if(journeyEntryFlags[j][i]) {
                    binaryCode |= 1;
                }
                binaryCode <<= 1;
            }
            journeyEntryMap[binaryCode].push_back(i);
        }

        if(debug) {
            cout << "Journey entry groups for vertex " << oldVertex.stringId() << "\n";
            for(const auto& p: journeyEntryMap) {
                cout << "Binary code " << p.first << ":\n";
                const vector<uint64_t>& journeyEntryGroup = p.second;
                for(const uint64_t i: journeyEntryGroup) {
                    const JourneyEntry& journeyEntry = oldVertex.journeyEntries[i];
                    cout << "(" << journeyEntry.orientedReadId << "," <<
                        journeyEntry.position << ") ";
                }
                cout << "\n";
            }
        }



        // Each group with identical journey entry flags generates a new vertex
        // in the new assembly graph.
        for(const auto& p: journeyEntryMap) {
            const vertex_descriptor vNew = boost::add_vertex(
                AssemblyGraphVertex(segmentId, verticesBySegment[segmentId].size()),
                newAssemblyGraph);
            verticesBySegment[segmentId].push_back(vNew);
            AssemblyGraphVertex& newVertex = newAssemblyGraph[vNew];

            const vector<uint64_t>& journeyEntryGroup = p.second;
            for(const uint64_t i: journeyEntryGroup) {
                const JourneyEntry& journeyEntry = oldVertex.journeyEntries[i];
                newVertex.journeyEntries.push_back(journeyEntry);

                // Update the journey for this oriented read.
                const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
                const uint64_t position = journeyEntry.position;
                vector<vertex_descriptor>& orientedReadJourney = journeys[orientedReadId.getValue()];
                SHASTA_ASSERT(orientedReadJourney[position] == null_vertex());
                orientedReadJourney[position] = vNew;
            }
        }
    }
    createLinks();

    if(true) {
        cout << "The new assembly graph created by "
            "AssemblyGraph::createFromTangledAssemblyPaths has " <<
            num_vertices(newAssemblyGraph) << " segments and " <<
            num_edges(newAssemblyGraph) << " links." << endl;
    }

    // Check that the journeys have been completely filled in.
    for(const auto& journey: journeys) {
        SHASTA_ASSERT(find(journey.begin(), journey.end(), null_vertex()) == journey.end());
    }

}



void AssemblyGraph::assemble()
{
    for(uint64_t pathId=0; pathId<assemblyPaths.size(); pathId++) {
        assemble(pathId);
    }
}



void AssemblyGraph::assemble(uint64_t assemblyPathId)
{
    const AssemblyGraph& assemblyGraph = *this;
    AssemblyPath& assemblyPath = *assemblyPaths[assemblyPathId];

    const bool debug = true;
    if(debug) {
        cout << "Assembling path " << assemblyPathId << " with " <<
            assemblyPath.primaryVertices.size() << " primary vertices." << endl;
        cout << "Path begins at " << vertexStringId(assemblyPath.primaryVertices.front()) <<
            " and ends at " << vertexStringId(assemblyPath.primaryVertices.back()) << endl;
    }



    // Construct the FlattenedAssemblyPath.
    FlattenedAssemblyPath flattenedAssemblyPath;



    // Construct the FlattenedAssemblyPath segments.
    vector<bool> useForAssembly;
    for(uint64_t i=0; /* Check later */; i++) {

        // Add the primary vertex at this position.
        // All of its journey entries are marked to be used for assembly.
        const vertex_descriptor v = assemblyPath.primaryVertices[i];
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        flattenedAssemblyPath.segments.emplace_back(v, vertex.journeyEntries.size());

        // If this is the last primary vertex, there are no secondary
        // vertices and we are done.
        if(i == assemblyPath.primaryVertices.size() - 1) {
            break;
        }

        // Add the secondary vertices following this primary vertex.
        const auto& secondaryVertexInfo = assemblyPath.secondaryVerticesInfos[i];
        for(const vertex_descriptor v: secondaryVertexInfo.secondaryVertices) {
            const AssemblyGraphVertex& vertex = assemblyGraph[v];
            secondaryVertexInfo.getVertexJourneyEntries(vertex, useForAssembly);
            SHASTA_ASSERT(useForAssembly.size() == vertex.journeyEntries.size());
            flattenedAssemblyPath.segments.emplace_back(v, useForAssembly);
        }
    }



    // Construct the FlattenedAssemblyPath links.
    for(uint64_t i=1; i<flattenedAssemblyPath.segments.size(); i++) {
        const vertex_descriptor v0 = flattenedAssemblyPath.segments[i-1].v;
        const vertex_descriptor v1 = flattenedAssemblyPath.segments[i].v;

        // Find the edge.
        edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(v0, v1, assemblyGraph);
        SHASTA_ASSERT(edgeWasFound);

        flattenedAssemblyPath.links.emplace_back(e);
    }

    if(debug) {
        cout << "The FlattenedAssemblyPath has " << flattenedAssemblyPath.segments.size() <<
            " segments." << endl;
    }


    // Assemble it.
    // This fills in the fields in the segments and links that
    // were not initialized here.
    assemble(flattenedAssemblyPath, assemblyPath.assembledSequence);
}



AssemblyGraph::FlattenedAssemblyPathLink::FlattenedAssemblyPathLink(edge_descriptor e) :
    e(e) {}



// Assemble a FlattenedAssemblyPath.
// This is similar to assembly code in AssemblyGraphSnapshot,
// but it only uses a subset of the JourneyEntries of each vertex,
// as specified by the FlattenedAssemblyPath.
void AssemblyGraph::assemble(
    FlattenedAssemblyPath& flattenedAssemblyPath,
    vector<shasta::Base>& sequence)
{
    const AssemblyGraph& assemblyGraph = *this;
    const bool debug = true;

    // Fill in sequence length for all segments.
    for(auto& segment: flattenedAssemblyPath.segments) {
        const AssemblyGraphVertex& vertex = assemblyGraph[segment.v];
        segment.sequenceLength = packedMarkerGraph.segmentSequences[vertex.segmentId].size();
    }

    // Assemble the links.
    SHASTA_ASSERT(flattenedAssemblyPath.links.size() == flattenedAssemblyPath.segments.size() - 1);
    for(uint64_t i=0; i<flattenedAssemblyPath.links.size(); i++) {
        auto& link = flattenedAssemblyPath.links[i];
        const auto& previousSegment = flattenedAssemblyPath.segments[i];
        const auto& nextSegment = flattenedAssemblyPath.segments[i+1];

        if(debug) {
            cout << "Assembling link at position " << i <<
                " between vertices " << vertexStringId(previousSegment.v) <<
                " and " << vertexStringId(nextSegment.v) << endl;

        }

        assembleLink(link, previousSegment, nextSegment);
    }

}



// Assemble a link, using only journey entries permitted
// by the useForAssembly fields of the previous and next segment.
// Much of this code is similar to AssemblyGraphSnapshot::assembleLink,
// but uses only a subset of the transitions.
void AssemblyGraph::assembleLink(
    FlattenedAssemblyPathLink& link,
    const FlattenedAssemblyPathSegment& segment0,
    const FlattenedAssemblyPathSegment& segment1
) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const bool debug = true;
    using shasta::Base;

    // Gather some information we need.
    const AssemblyGraphVertex& vertex0 = assemblyGraph[segment0.v];
    const AssemblyGraphVertex& vertex1 = assemblyGraph[segment1.v];

    const uint64_t segmentId0 = vertex0.segmentId;
    const uint64_t segmentId1 = vertex1.segmentId;

    const auto leftPath = packedMarkerGraph.segments[segmentId0];
    const auto rightPath = packedMarkerGraph.segments[segmentId1];

    const auto leftSegmentSequence = packedMarkerGraph.segmentSequences[segmentId0];
    const auto rightSegmentSequence = packedMarkerGraph.segmentSequences[segmentId1];

    const auto leftVertexOffsets = packedMarkerGraph.segmentVertexOffsets[segmentId0];
    const auto rightVertexOffsets = packedMarkerGraph.segmentVertexOffsets[segmentId1];

    const uint64_t k = packedMarkerGraph.k;

    // Find the Transitions to be used to assemble this link.
    vector<Transition> transitions;
    getTransitionsForAssembly(segment0, segment1, transitions);
    SHASTA_ASSERT(not transitions.empty());

    if(debug) {
        cout << "Link assembly will use " << transitions.size() <<
            " oriented reads." << endl;
    }

    // Loop over transitions to compute the maximum left/right skip,
    // that is, the maximum number of marker graph edges skipped
    // by an oriented read at the end of the left segment or
    // at the beginning of the right segment.
    // This is needed below to compute the sequences that participate in the MSA.
    uint64_t maxLeftSkip = 0;
    uint64_t maxRightSkip = 0;
    for(const Transition& transition: transitions) {
        const OrientedReadId orientedReadId = transition.orientedReadId;
        const auto journey = packedMarkerGraph.journeys[orientedReadId.getValue()];
        SHASTA_ASSERT(transition.position0 + 1 == transition.position1);
        const auto& leftJourneyStep = journey[transition.position0];
        const auto& rightJourneyStep = journey[transition.position1];
        const uint64_t leftSkip = leftPath.size() - leftJourneyStep.positions[1];
        const uint64_t rightSkip = rightJourneyStep.positions[0];
        maxLeftSkip = max(maxLeftSkip, leftSkip);
        maxRightSkip = max(maxRightSkip, rightSkip);
    }

    // Compute the position in the left segment of the begining of the left segment
    // sequence that will be used to fill in MSA sequence.
    const uint64_t leftPositionBegin = leftVertexOffsets[leftVertexOffsets.size() -1 - maxLeftSkip];

    // Compute the position in the right segment of the end of the right segment
    // sequence that will be used to fill in MSA sequence.
    const uint64_t rightPositionEnd = rightVertexOffsets[maxRightSkip] + k;
    // A vector to contain the distinct MSA sequence we found, each with the number of times it was found.
    vector< pair<vector<Base>, uint64_t> > msaSequences;

    // A vector that, for each transition, gives the index in msaSequence.
    vector<uint64_t> msaSequenceTable(transitions.size());



    // Loop over transitions to compute the MSA sequence to be used for each oriented read.
    vector<Base> msaSequence;
    for(uint64_t iTransition=0; iTransition<transitions.size(); iTransition++) {
        const Transition& transition = transitions[iTransition];
        const OrientedReadId orientedReadId = transition.orientedReadId;
        const auto journey = packedMarkerGraph.journeys[orientedReadId.getValue()];
        SHASTA_ASSERT(transition.position0 + 1 == transition.position1);
        const auto& leftJourneyStep = journey[transition.position0];
        const auto& rightJourneyStep = journey[transition.position1];
        const uint64_t leftSkip = leftPath.size() - leftJourneyStep.positions[1];
        const uint64_t rightSkip = rightJourneyStep.positions[0];
        const uint64_t leftOrdinal = leftJourneyStep.ordinals[1];
        const uint64_t rightOrdinal = rightJourneyStep.ordinals[0];
        SHASTA_ASSERT(rightOrdinal >= leftOrdinal);
        /*
        const uint64_t ordinalSkip = rightOrdinal - leftOrdinal;
        const int64_t linkSeparation =
            int64_t(ordinalSkip) -
            int64_t(leftSkip) -
            int64_t(rightSkip);
        */
        const auto orientedReadMarkers = packedMarkerGraph.markers[orientedReadId.getValue()];
        const CompressedMarker& marker0 = orientedReadMarkers[leftOrdinal];
        const CompressedMarker& marker1 = orientedReadMarkers[rightOrdinal];

        // Now we can compute the portion of the oriented read sequence that will
        // participate in the MSA.
        // This will be extended to the left/right as necessary,
        // using the sequence of the left/right segment.
        const uint64_t positionBegin = marker0.position;
        const uint64_t positionEnd = marker1.position + k;

        // Compute the position of the left extension in the left segment.
        const uint64_t leftPositionEnd = leftVertexOffsets[leftVertexOffsets.size() - 1 - leftSkip];

        // Compute the position of the right extension in the right segment.
        const uint64_t rightPositionBegin = rightVertexOffsets[rightSkip] + k;

        // Add the left extension to the MSA sequence.
        msaSequence.clear();
        for(uint64_t position=leftPositionBegin; position!=leftPositionEnd; ++position) {
            msaSequence.push_back(leftSegmentSequence[position]);
        }

        // Add the oriented read to the MSA sequence.
        for(uint64_t position=positionBegin; position!=positionEnd; ++position) {
            msaSequence.push_back(packedMarkerGraph.reads.getOrientedReadBase(orientedReadId, uint32_t(position)));
        }

        // Add the right extension to the MSA sequence.
        for(uint64_t position=rightPositionBegin; position!=rightPositionEnd; ++position) {
            msaSequence.push_back(rightSegmentSequence[position]);
        }

        // Update the msaSequences and msaSequenceTable.
        bool done = false;
        for(uint64_t i=0; i<msaSequences.size(); i++) {
            if(msaSequence == msaSequences[i].first) {
                msaSequenceTable[iTransition] = i;
                ++msaSequences[i].second;
                done = true;
                break;
            }
        }
        if(not done) {
            msaSequenceTable[iTransition] = msaSequences.size();
            msaSequences.push_back(make_pair(msaSequence, 1));
        }
    }

    if(debug) {
        cout << "MSA sequences with coverage:" << endl;
        for(const auto& p: msaSequences) {
            const vector<Base>& sequence = p.first;
            const uint64_t coverage = p.second;
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(cout));
            cout << " " << coverage << endl;
        }
    }

    // Compute the multiple sequence alignment.
    vector<Base> consensusSequence;
    const uint64_t maxLength = 10000;
    ostream html(0);
    linkMsaUsingSpoa(msaSequences, maxLength, html, consensusSequence);

    if(debug) {
        cout << "Consensus sequence has length " << consensusSequence.size() << ":\n";
        copy(consensusSequence.begin(), consensusSequence.end(),
            ostream_iterator<Base>(cout));
        cout << "\n";
    }

    // Compute the number of bases at the beginning of the consensus sequence
    // that are identical to the corresponding bases in the left segment.
    uint64_t leftIdentical = 0;
    for(uint64_t i=leftPositionBegin; i<leftSegmentSequence.size(); i++) {
        if(i-leftPositionBegin >= consensusSequence.size()) {
            break;
        }
        if(leftSegmentSequence[i] == consensusSequence[i-leftPositionBegin]) {
            ++leftIdentical;
        } else {
            break;
        }
    }
    if(debug) {
        cout << "Number of link consensus bases identical to previous segment " << leftIdentical << "\n";
    }

    // Compute the number of bases at the end of the consensus sequence
    // that are identical to the corresponding bases in the right segment.
    uint64_t rightIdentical = 0;
    uint64_t j = rightPositionEnd - 1; // j is index in right segment sequence.
    for(uint64_t i=consensusSequence.size()-1; /* Check later */; i--, j--) { // i is index in consensus sequence.
        if(consensusSequence[i] == rightSegmentSequence[j]) {
            ++rightIdentical;
        } else {
            break;
        }
        if(i == 0) {
            break;
        }
        if(j == 0) {
            break;
        }
    }
    if(debug) {
        cout << "Number of link consensus bases identical to next segment " << rightIdentical << "\n";
    }



    // Trim consensus bases at the end of the consensusSequence that are identical
    // to the corresponding bases in the right segment sequence.
    // Here, rightOverride is the number of bases at the beginning of the right segment
    // that are overridden by bases in the consensus sequence of the link.
    link.rightOverride = rightPositionEnd;
    while(not consensusSequence.empty()) {
        if(consensusSequence.back() == rightSegmentSequence[link.rightOverride-1]) {
            consensusSequence.resize(consensusSequence.size() - 1);
            --link.rightOverride;
            if(link.rightOverride == 0) {
                break;
            }
        } else {
            break;
        }
    }

    // Same, on the left.
    link.leftOverride = leftSegmentSequence.size() - leftPositionBegin;
    uint64_t leftTrim = 0;
    for(uint64_t i=0; i<consensusSequence.size(); i++) {
        const uint64_t j = i + leftPositionBegin;
        if(j >= leftSegmentSequence.size()) {
            break;
        }
        if(consensusSequence[i] == leftSegmentSequence[j]) {
            ++leftTrim;
            --link.leftOverride;
            if(link.leftOverride == 0) {
                break;
            }
        }
    }
    copy(consensusSequence.begin() + leftTrim, consensusSequence.end(), consensusSequence.begin());
    consensusSequence.resize(consensusSequence.size() - leftTrim);

    // Store the trimmed consensus sequence.
    link.sequence = consensusSequence;


    if(debug) {
        cout << "After trimming sequence identical to adjacent segments, "
            "consensus sequence has length " << consensusSequence.size() << ":\n";
        copy(consensusSequence.begin(), consensusSequence.end(),
            ostream_iterator<Base>(cout));
        cout << "\n";
        cout << "Number of previous segment bases overridden by link consensus:" << link.leftOverride << "\n";
        cout << "Number of next segment bases overridden by link consensus:" << link.rightOverride << "\n";


        cout << "Assembly of the path consisting of this link plus the adjacent segments:\n";
        copy(leftSegmentSequence.begin(), leftSegmentSequence.end() - link.leftOverride,
            ostream_iterator<Base>(cout));
        copy(consensusSequence.begin(), consensusSequence.end(),
            ostream_iterator<Base>(cout));
        copy(rightSegmentSequence.begin()+ link.rightOverride, rightSegmentSequence.end(),
            ostream_iterator<Base>(cout));
        cout << "\n";
    }
}



// Find the Transitions to be used to assemble this link.
// These transitions are defined as follows:
// - The OrientedReadId has a JourneyEntry in the previous
//   segment that is flagged as usable for assembly.
// - The OrientedReadId has a JourneyEntry in the next
//   segment that is flagged as usable for assembly.
// - The positions of the two journey entries differ by 1.
// To find them, we use the fact that the journey entries of a vertex
// are ordered by OrientedReadId and the by position.
// Use index 0 for the previous segment/vertex
// and index 1 for the next segment/vertex.
void AssemblyGraph::getTransitionsForAssembly(
    const FlattenedAssemblyPathSegment& segment0,
    const FlattenedAssemblyPathSegment& segment1,
    vector<Transition>& transitions) const
{
    const AssemblyGraph& assemblyGraph = *this;
    transitions.clear();

    // Extract some information we need below.
    const AssemblyGraphVertex& vertex0 = assemblyGraph[segment0.v];
    const AssemblyGraphVertex& vertex1 = assemblyGraph[segment1.v];

    const vector<JourneyEntry>& entries0 = vertex0.journeyEntries;
    const vector<JourneyEntry>& entries1 = vertex1.journeyEntries;

    const uint64_t n0 = entries0.size();
    const uint64_t n1 = entries1.size();

    SHASTA_ASSERT(segment0.useForAssembly.size() == n0);
    SHASTA_ASSERT(segment1.useForAssembly.size() == n1);



    const bool debug = false;
    if(debug) {
        cout << "AssemblyGraph::getTransitionsForAssembly called at " <<
            vertexStringId(segment0.v) << " " <<
            vertexStringId(segment1.v) << endl;
        cout << "Entries0:" << endl;
        for(uint64_t i=0; i<n0; i++) {
            cout << entries0[i].orientedReadId << " " <<
                entries0[i].position  << " " <<
                int(segment0.useForAssembly[i]) << endl;
        }
        cout << "Entries1:" << endl;
        for(uint64_t i=0; i<n1; i++) {
            cout << entries1[i].orientedReadId << " " <<
                entries1[i].position  << " " <<
                int(segment1.useForAssembly[i]) << endl;
        }
    }



    // Joint iteration over the journey entries of the two vertices.
    uint64_t i0 = 0;
    uint64_t i1 = 0;
    while(i0 < n0 and i1 < n1) {
        const JourneyEntry& entry0 = entries0[i0];
        const JourneyEntry& entry1 = entries1[i1];

        // Check the OrientedReadIds.
        if(entry0.orientedReadId < entry1.orientedReadId) {
            ++i0;
            continue;
        }

        if(entry0.orientedReadId > entry1.orientedReadId) {
            ++i1;
            continue;
        }

        const OrientedReadId orientedReadId = entry0.orientedReadId;
        SHASTA_ASSERT(entry1.orientedReadId == orientedReadId);

        // We have two entries for the same OrientedReadId.
        // We are only interested if the positions differ by 1.
        if(entry0.position + 1 < entry1.position) {
            ++i0;
            continue;
        }

        if(entry0.position + 1 > entry1.position) {
            ++i1;
            continue;
        }
        SHASTA_ASSERT(entry0.position + 1 == entry1.position);

        // If these entries are usable for assembly, they generate
        // a new Transition.
        if(segment0.useForAssembly[i0] and segment1.useForAssembly[i1]) {
            Transition transition;
            transition.orientedReadId = orientedReadId;
            transition.position0 = entry0.position;
            transition.position1 = entry1.position;
            transitions.push_back(transition);
        }
        ++i0;
        ++i1;
    }
}




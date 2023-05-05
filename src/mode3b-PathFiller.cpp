// Shasta.
#include "mode3b-PathFiller.hpp"
#include "Assembler.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"



PathFiller::PathFiller(
    const Assembler& assembler,
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    ostream& html) :
    assembler(assembler),
    edgeIdA(edgeIdA),
    edgeIdB(edgeIdB),
    graph(assembler)
{
    checkAssumptions();
    if(html) {
        html << "<h1>Assembly path step between marker graph edges " <<
            edgeIdA << " and " << edgeIdB << "</h1>";
    }

    gatherOrientedReads();
    if(html) {
        html << "<p>Coverage " << orientedReadInfos.size();
    }

    createGraph();
    fillPathGreedy(html);

    if(html) {
        writeSequence(html);
        html << "<p>The local marker graph for this step assembly has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
        writeGraph(html);
    }

}



void PathFiller::checkAssumptions() const
{
    SHASTA_ASSERT(edgeIdA != edgeIdB);
    SHASTA_ASSERT(assembler.assemblerInfo->assemblyMode == 3);
    SHASTA_ASSERT(assembler.getReads().representation == 0);
    SHASTA_ASSERT(not assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA));
    SHASTA_ASSERT(not assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB));
}



// Use a joint loop over the MarkerIntervals of the two edges
// to gather the oriented reads that we will use for this assembly.
void PathFiller::gatherOrientedReads()
{
    const auto markerIntervalsA = assembler.markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = assembler.markerGraph.edgeMarkerIntervals[edgeIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();
    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common oriented read.
        OrientedReadInfo info;
        info.orientedReadId = itA->orientedReadId;
        info.ordinalA0 = itA->ordinals[0];
        info.ordinalA1 = itA->ordinals[1];
        info.ordinalB0 = itB->ordinals[0];
        info.ordinalB1 = itB->ordinals[1];
        SHASTA_ASSERT(info.ordinalA1 == info.ordinalA0 + 1);
        SHASTA_ASSERT(info.ordinalB1 == info.ordinalB0 + 1);

        // If the edges don't appear in the expected order in this oriented read,
        // skip it.
        if(info.ordinalB0 < info.ordinalA1) {
            continue;
        }

        // Store it.
        orientedReadInfos.push_back(info);

        // Continue the joint loop.
        ++itA;
        ++itB;
    }

}



PathFiller::Graph::Graph(const Assembler& assembler) :
    assembler(assembler)
{
}



PathFiller::Graph::vertex_descriptor PathFiller::Graph::addVertexIfNecessary(
    MarkerGraphVertexId vertexId)
{
    // See if we already have this vertex.
    const auto it = vertexMap.find(vertexId);

    if(it == vertexMap.end()) {

        // This vertex does not exist. Create it.
        const vertex_descriptor v = add_vertex(Vertex({vertexId}), *this);

        // Add it to the vertex map.
        vertexMap.insert(make_pair(vertexId, v));

        return v;

    } else {

        // we already have this vertex.
        return it->second;
    }
}



PathFiller::Graph::edge_descriptor PathFiller::Graph::addEdgeIfNecessary(
    MarkerGraphEdgeId edgeId)
{
    // See if we already have this edge.
    const auto it = edgeMap.find(edgeId);

    if(it == edgeMap.end()) {

        // This edge does not exist.
        // Before we can create, we need to create its vertices, if necessary.
        const MarkerGraph::Edge& markerGraphEdge = assembler.markerGraph.edges[edgeId];
        const MarkerGraphVertexId vertexId0 = markerGraphEdge.source;
        const MarkerGraphVertexId vertexId1 = markerGraphEdge.target;
        const vertex_descriptor v0 = addVertexIfNecessary(vertexId0);
        const vertex_descriptor v1 = addVertexIfNecessary(vertexId1);

        // Now we can create the edge.
        edge_descriptor e;
        bool edgeWasAdded = false;
        tie(e, edgeWasAdded) = add_edge(v0, v1, Edge({edgeId}), *this);
        SHASTA_ASSERT(edgeWasAdded);

        // Add it to the edge map.
        edgeMap.insert(make_pair(edgeId, e));

        return e;

    } else {

        // We already have this edge.
        return it->second;
    }
}



void PathFiller::createGraph()
{
    // Loop over all oriented reads and marker intervals that contribute
    // to assembling this step.
    for(const OrientedReadInfo& info: orientedReadInfos) {
        for(uint32_t ordinal0=info.ordinalA0; ordinal0!=info.ordinalB1; ordinal0++) {
            const uint32_t ordinal1 = ordinal0 + 1;
            MarkerInterval markerInterval(info.orientedReadId, ordinal0, ordinal1);

            // Find the marker graph edge that contains this MarkerInterval.
            const MarkerGraphEdgeId edgeId =
                assembler.markerGraph.locateMarkerInterval(assembler.markers, markerInterval);
            SHASTA_ASSERT(edgeId != invalid<MarkerGraphEdgeId>);

            // Create an edge corresponding to this marker graph edge,
            // if we don't already have one.
            const Graph::edge_descriptor e = graph.addEdgeIfNecessary(edgeId);

            // Add this marker interval to that edge.
            graph[e].markerIntervals.push_back(markerInterval);

        }

    }
}



void PathFiller::Graph::writeGraphviz(
    ostream& out,
    uint64_t peakCoverage,
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    const vector<MarkerGraphEdgeId>& pathSecondaryEdges
    ) const
{
    const Graph& graph = *this;
    const double S = 0.7;
    const double V = 1.;

    // Gather the path edges so we can color them differently.
    vector<MarkerGraphEdgeId> sortedPathEdges = pathSecondaryEdges;
    sortedPathEdges.push_back(edgeIdA);
    sortedPathEdges.push_back(edgeIdB);
    sort(sortedPathEdges.begin(), sortedPathEdges.end());

    out <<
        "digraph PathFillerGraph {\n"
        "ranksep=0.3;\n"
        "node [shape=point fontname=\"Courier New\"];\n"
        "edge [decorate=true penwidth=10 arrowsize=0.5];\n";

    BGL_FORALL_EDGES(e, graph, Graph) {

        const Edge& edge = graph[e];
        const uint64_t coverage = edge.markerIntervals.size();
        const bool isOnPath = find(sortedPathEdges.begin(), sortedPathEdges.end(), edge.edgeId) != sortedPathEdges.end();
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);
        const MarkerGraphVertexId vertexId0 = graph[v0].vertexId;
        const MarkerGraphVertexId vertexId1 = graph[v1].vertexId;

        // Compute the hue based on coverage.
        double H;
        if(coverage >= peakCoverage) {
            H = 1./3.;
        } else {
            H = (double(coverage-1) / (3. * double(peakCoverage - 1)));
        }
        const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";


        // Write it as an intermediate vertex, for better readability.
        out <<
            "E" << edge.edgeId <<
            "[ shape=rectangle label=\"" << edge.edgeId << "\\n" << coverage << "\\n";
        auto sequence = edgeSequence(e);
        copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(out));
        out << "\"";
        if(isOnPath) {
            out << " style=filled fillcolor=LightCyan";
        }
        out << "];\n";

        out << vertexId0 << "->E" << edge.edgeId << " [dir=none";
        out << " color=" << colorString;
        out << "];\n";
        out << "E" << edge.edgeId << "->" << vertexId1;
        out << "[color=" << colorString;
        out << "];\n";
    }

    out << "}\n";

}

#if 0
s << " label=<<font color=\"black\">";
s << "<table";
s << " color=\"black\"";
s << " bgcolor=\"" << labelColor << "\"";
s << " border=\"1\"";
s << " cellborder=\"0\"";
s << " cellspacing=\"2\"";
s << ">";
#endif


void PathFiller::writeGraph(ostream& html) const
{
    // Write out the graph in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    {
        ofstream dotFile(dotFileName);
        graph.writeGraphviz(dotFile, orientedReadInfos.size(), edgeIdA, edgeIdB, pathSecondaryEdges);
    }

    // Compute layout in svg format.
    const string command = timeoutCommand() + " 30 dot -O -T svg " + dotFileName;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    const double timeout = 30;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    if(returnCode!=0 or signalOccurred) {
        throw runtime_error("An error occurred while running the following command: " + command);
    }
    if(timeoutTriggered) {
        std::filesystem::remove(dotFileName);
        throw runtime_error("Timeout during graph layout computation.");
    }

    // Remove the .dot file.
    std::filesystem::remove(dotFileName);

    // Copy the svg file to html.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<p>" << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);
}



void PathFiller::fillPathGreedy(ostream& html)
{
    const MarkerGraph::Edge markerGraphEdgeA = assembler.markerGraph.edges[edgeIdA];
    const MarkerGraph::Edge markerGraphEdgeB = assembler.markerGraph.edges[edgeIdB];
    const Graph::vertex_descriptor vA = graph.vertexMap[markerGraphEdgeA.source];
    const Graph::vertex_descriptor vB = graph.vertexMap[markerGraphEdgeB.target];

    std::set<Graph::edge_descriptor> edgesUsed;

    // Main iteration loop.
    Graph::vertex_descriptor v = vA;
    while(v != vB) {

        // Find the edge with the most coverage.
        Graph::edge_descriptor eNext;
        uint64_t bestCoverage = 0;
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            if(edgesUsed.contains(e)) {
                continue;
            }
            const uint64_t coverage = graph[e].markerIntervals.size();
            if(coverage > bestCoverage) {
                eNext = e;
                bestCoverage = coverage;
            }
        }
        if(bestCoverage == 0) {
            throw runtime_error("Unable to fill assembly path.");
        }

        edgesUsed.insert(eNext);
        pathSecondaryEdges.push_back(graph[eNext].edgeId);
        v = target(eNext, graph);
    }

    // As constructed, pathSecondaryEdges includes edgeIdA and edgeIB.
    // Check that this is the case, then remove them.
    SHASTA_ASSERT(pathSecondaryEdges.size() >= 2);
    SHASTA_ASSERT(pathSecondaryEdges.front() == edgeIdA);
    SHASTA_ASSERT(pathSecondaryEdges.back() == edgeIdB);
    pathSecondaryEdges.resize(pathSecondaryEdges.size() - 1);
    pathSecondaryEdges.erase(pathSecondaryEdges.begin());

    if(html) {
        html << "<p>Path secondary edges:";
        for(const MarkerGraphEdgeId edgeId: pathSecondaryEdges) {
            html << " " << edgeId;
        }
    }

}



span<const shasta::Base> PathFiller::Graph::edgeSequence(
    edge_descriptor e) const
{
    return edgeSequence((*this)[e].edgeId);
}



span<const shasta::Base> PathFiller::Graph::edgeSequence(
    MarkerGraphEdgeId edgeId) const
{
    return assembler.markerGraph.edgeSequence[edgeId];
}



// Get the sequence.
// The sequences of edgeIdA and edgeIdB are only included if
// includePrimary is true.
void PathFiller::getSequence(
    vector<Base>& sequence,
    bool includePrimary) const
{
    sequence.clear();

    if(includePrimary) {
        const auto sequenceA = graph.edgeSequence(edgeIdA);
        copy(sequenceA.begin(), sequenceA.end(), back_inserter(sequence));
    }

    for(const MarkerGraphEdgeId edgeId: pathSecondaryEdges) {
        const auto edgeSequence = graph.edgeSequence(edgeId);
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
    }

    if(includePrimary) {
        const auto sequenceB = graph.edgeSequence(edgeIdB);
        copy(sequenceB.begin(), sequenceB.end(), back_inserter(sequence));
    }
}



void PathFiller::writeSequence(ostream& html) const
{
    if(not html) {
        return;
    }

    const auto sequenceA = graph.edgeSequence(edgeIdA);
    const auto sequenceB = graph.edgeSequence(edgeIdB);
    vector<Base> sequence;
    getSequence(sequence, false);

    html << "<pre style='font-family:monospace'>\n";
    html << ">" << edgeIdA << " length " << sequenceA.size() << "\n";
    copy(sequenceA.begin(), sequenceA.end(), ostream_iterator<shasta::Base>(html));
    html << "\n>Intervening length " << sequence.size() << "\n";
    copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(html));
    html << "\n>" << edgeIdB << " length " << sequenceB.size() << "\n";
    copy(sequenceB.begin(), sequenceB.end(), ostream_iterator<shasta::Base>(html));
    html << "\n>Combined length " << sequenceA.size() + sequence.size() + sequenceB.size() << "\n";
    copy(sequenceA.begin(), sequenceA.end(), ostream_iterator<shasta::Base>(html));
    copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(html));
    copy(sequenceB.begin(), sequenceB.end(), ostream_iterator<shasta::Base>(html));
    html << "\n</pre>";
}

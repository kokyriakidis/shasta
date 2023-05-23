// Shasta.
#include "mode3b-PathFiller.hpp"
#include "approximateTopologicalSort.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "platformDependent.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Seqan
#include <seqan/align.h>
#include <seqan/graph_msa.h>

// Standard library.
#include "fstream.hpp"



PathFiller::PathFiller(
    const Assembler& assembler,
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    ostream& html,
    bool showGraph,
    bool showVertices,
    bool showVertexLabels,
    bool showEdgeLabels) :
    assembler(assembler),
    edgeIdA(edgeIdA),
    edgeIdB(edgeIdB)
{
    const PathFiller& graph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxBaseSkip = 300;

    checkAssumptions();
    if(html) {
        html << "<h1>Assembly path step between marker graph edges " <<
            edgeIdA << " and " << edgeIdB << "</h1>";
    }

    gatherOrientedReads();
    if(html) {
        writeOrientedReads(html);
    }

    createGraph(maxBaseSkip);
    computeStrongComponents();
    createVirtualEdges();
    writeVerticesCsv();
    approximateTopologicalSort();
    if(html) {
        html << "<p>The local marker graph for this step assembly has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
    }

    // Assemble the path between edgeIdA and edgeIdB.
    findAssemblyPath();
    if(html) {
        writeSequence(html);
        ofstream fasta("AssemblyPath.fasta");
        writeSequenceFasta(fasta);
        ofstream csv("AssemblyPath.csv");
        writeAssemblyDetails(csv);
        if(showGraph) {
            writeGraph(
                html,
                showVertices,
                showVertexLabels,
                showEdgeLabels);
        }
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

        // Fill in the rest of the information.
        {
            OrientedReadInfo& info = orientedReadInfos.back();
            const OrientedReadId orientedReadId = info.orientedReadId;

            info.ordinalOffset = info.ordinalB1 - info.ordinalA0;

            const MarkerId markerIdA0 =assembler.getMarkerId(orientedReadId, info.ordinalA0);
            const MarkerId markerIdA1 =assembler.getMarkerId(orientedReadId, info.ordinalA1);
            const MarkerId markerIdB0 =assembler.getMarkerId(orientedReadId, info.ordinalB0);
            const MarkerId markerIdB1 =assembler.getMarkerId(orientedReadId, info.ordinalB1);

            info.positionA0 = assembler.markers.begin()[markerIdA0].position;
            info.positionA1 = assembler.markers.begin()[markerIdA1].position;
            info.positionB0 = assembler.markers.begin()[markerIdB0].position;
            info.positionB1 = assembler.markers.begin()[markerIdB1].position;

            info.baseOffset = info.positionB1 - info.positionA0;
        }

        // Continue the joint loop.
        ++itA;
        ++itB;
    }


    // Compute average ordinal and base offsets.
    uint64_t ordinalOffsetSum = 0.;
    uint64_t baseOffsetSum = 0;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        ordinalOffsetSum += info.ordinalOffset;
        baseOffsetSum += info.baseOffset;
    }
    ordinalOffset = uint32_t(std::round(double(ordinalOffsetSum) / double(orientedReadInfos.size())));
    baseOffset = uint32_t(std::round(double(baseOffsetSum) / double(orientedReadInfos.size())));

}



void PathFiller::writeOrientedReads(ostream& html) const
{
    html <<
        "<p><table>"
        "<tr><th class=left>Coverage<td class=centered>" << orientedReadInfos.size() <<
        "<tr><th class=left>Ordinal offset<td class=centered>" << ordinalOffset <<
        "<tr><th class=left>Base offset<td class=centered>" << baseOffset <<
        "</table>";

    html <<
        "<p><table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th>OrdinalA0"
        "<th>OrdinalA1"
        "<th>OrdinalB0"
        "<th>OrdinalB1"
        "<th>PositionA0"
        "<th>PositionA1"
        "<th>PositionB0"
        "<th>PositionB1"
        "<th>OrdinalOffset"
        "<th>BaseOffset";

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        html <<
            "<tr>"
            "<td class=centered>" << info.orientedReadId <<
            "<td class=centered>" << info.ordinalA0 <<
            "<td class=centered>" << info.ordinalA1 <<
            "<td class=centered>" << info.ordinalB0 <<
            "<td class=centered>" << info.ordinalB1 <<
            "<td class=centered>" << info.positionA0 <<
            "<td class=centered>" << info.positionA1 <<
            "<td class=centered>" << info.positionB0 <<
            "<td class=centered>" << info.positionB1 <<
            "<td class=centered>" << info.ordinalOffset <<
            "<td class=centered>" << info.baseOffset;
    }
    html << "</table>";
}



void PathFiller::createGraph(uint64_t maxBaseSkip)
{
    createVertices();
    splitVertices(maxBaseSkip);
    createEdges();
}



// Approximate topological sort is only use to improve the
// display of the graph and make it faster.
// It sets the isDagEdge flags in the edges.
void PathFiller::approximateTopologicalSort()
{
    PathFiller& graph = *this;

    vector< pair<edge_descriptor, uint64_t> > edgesWithCoverage;
    BGL_FORALL_EDGES(e, graph, PathFiller) {
        edgesWithCoverage.push_back({e, graph[e].coverage()});
    }
    sort(edgesWithCoverage.begin(), edgesWithCoverage.end(),
        OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());

    vector<edge_descriptor> sortedEdges;
    for(const auto& p: edgesWithCoverage) {
        sortedEdges.push_back(p.first);
    }

    shasta::approximateTopologicalSort(graph, sortedEdges);
}



void PathFiller::createVertices()
{
    PathFiller& graph = *this;

    // During this phase there is at most one vertex corresponding to
    // each marker graph vertex.
    std::map<MarkerGraphVertexId, vertex_descriptor> vertexMap;


    // Loop over our oriented reads.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        OrientedReadInfo& info = orientedReadInfos[i];
        const OrientedReadId orientedReadId = info.orientedReadId;
        SHASTA_ASSERT(info.vertices.empty());

        // For each ordinal between ordinalA0 and ordinalB1 included,
        // make sure we have a vertex and store it.
        for(uint32_t ordinal=info.ordinalA0; ordinal<=info.ordinalB1; ordinal++) {

            // Find the marker graph vertex that contains this ordinal.
            const MarkerGraphVertexId vertexId =
                assembler.getGlobalMarkerGraphVertex(orientedReadId, ordinal);
            SHASTA_ASSERT(vertexId != invalid<MarkerGraphVertexId>);

            // Get the corresponding vertex descriptor, creating the vertex
            // if necessary.
            vertex_descriptor v;
            auto it = vertexMap.find(vertexId);
            if(it == vertexMap.end()) {
                v = add_vertex(PathFillerVertex(vertexId, orientedReadInfos.size()), graph);
                vertexMap.insert(make_pair(vertexId, v));
            } else {
                v = it->second;
            }

            // Store this ordinal in the vertex.
            PathFillerVertex& vertex = graph[v];
            vertex.ordinals[i].push_back(ordinal);

            // Store this vertex in the OrientedReadInfo.
            info.vertices.push_back(v);
        }
    }
}



void PathFiller::splitVertices(uint64_t maxBaseSkip)
{
    PathFiller& graph = *this;

    // Gather the vertices we have now so we can safely iterate over them.
    vector<vertex_descriptor> initialVertices;
    BGL_FORALL_VERTICES(v, graph, PathFiller) {
        initialVertices.push_back(v);
    }

    class OrdinalInfo {
    public:
        uint64_t i;
        uint32_t ordinal;
        int32_t estimatedOffset;
        bool operator<(const OrdinalInfo& that) const
        {
            return estimatedOffset < that.estimatedOffset;
        }
    };

    // Loop over our initial vertices.
    vector<OrdinalInfo> ordinalInfos;
    vector<uint64_t> splitPositions;
    for(const vertex_descriptor v: initialVertices) {
        const PathFillerVertex& vertex = graph[v];
        SHASTA_ASSERT(vertex.ordinals.size() == orientedReadInfos.size());

        // Loop over all ordinals of all oriented reads in this vertex.
        ordinalInfos.clear();
        for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
            const OrientedReadInfo& info = orientedReadInfos[i];
            const OrientedReadId orientedReadId = info.orientedReadId;

            for(const uint32_t ordinal: vertex.ordinals[i]) {
                // Get the position of this marker.
                const MarkerId markerId =assembler.getMarkerId(orientedReadId, ordinal);
                const uint32_t position = assembler.markers.begin()[markerId].position;

                // Base offsets from A0 and to B1.
                const int32_t offsetFromA0 = int32_t(position - info.positionA0);
                const int32_t offsetToB1 = int32_t(info.positionB1 - position);

                // Estimated offset from A0.
                const int32_t estimatedOffset = (offsetFromA0 + int32_t(baseOffset) - offsetToB1) / 2;

                ordinalInfos.push_back({i, ordinal, estimatedOffset});
            }
        }

        // Sort them by estimated offset.
        sort(ordinalInfos.begin(), ordinalInfos.end());

        /*
        cout << vertex.vertexId << ":";
        for(const auto& info: ordinalInfos) {
            cout << " " << info.estimatedOffset;
        }
        cout << endl;
        */

        // Look for skips larger than maxBaseSkip.
        splitPositions.clear();
        splitPositions.push_back(0);
        for(uint64_t i=1; i<ordinalInfos.size(); i++) {
            const int32_t baseSkip = ordinalInfos[i].estimatedOffset - ordinalInfos[i-1].estimatedOffset;
            if(baseSkip > int32_t(maxBaseSkip)) {
                splitPositions.push_back(i);
            }
        }
        splitPositions.push_back(ordinalInfos.size());

        if(splitPositions.size() == 2) {
            // cout << "No skips found." << endl;
            continue;
        }
        // cout << "Found " << splitPositions.size() - 2 << " skips." << endl;
        // cout << "This vertex will be split in " << splitPositions.size() - 1 << "." << endl;

        // Create the new vertices.
        for(uint64_t i=0; i<splitPositions.size()-1; i++) {
            const vertex_descriptor vNew = add_vertex(graph);
            PathFillerVertex& vertexNew = graph[vNew];
            vertexNew.vertexId = vertex.vertexId;
            vertexNew.replicaIndex = i + 1;
            vertexNew.ordinals.resize(orientedReadInfos.size());

            // Loop over the ordinals that go in this new vertex.
            const uint64_t begin = splitPositions[i];
            const uint64_t end = splitPositions[i+1];
            for(uint64_t j=begin; j!=end; j++) {
                const OrdinalInfo& ordinalInfo = ordinalInfos[j];
                const uint64_t i = ordinalInfo.i;
                vertexNew.ordinals[i].push_back(ordinalInfo.ordinal);

                // Store this vertex in the OrientedReadInfo.
                OrientedReadInfo& info = orientedReadInfos[i];
                info.vertices[ordinalInfo.ordinal - info.ordinalA0] = vNew;
            }
        }

        // Remove the initial vertex.
        remove_vertex(v, graph);
    }
}



void PathFiller::createEdges()
{
    PathFiller& graph = *this;

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        for(uint32_t ordinal0=info.ordinalA0; ordinal0<info.ordinalB1; ordinal0++) {
            const uint32_t ordinal1 = ordinal0 + 1;

            const vertex_descriptor v0 = info.getVertex(ordinal0);
            const vertex_descriptor v1 = info.getVertex(ordinal1);

            // Find the marker graph edge that contains this MarkerInterval.
            const MarkerGraphEdgeId edgeId =
                assembler.markerGraph.locateMarkerInterval(assembler.markers,
                    MarkerInterval(info.orientedReadId, ordinal0, ordinal1));
            SHASTA_ASSERT(edgeId != invalid<MarkerGraphEdgeId>);

            // If we already have this edge, add this ordinal to it.
            bool done = false;
            BGL_FORALL_OUTEDGES(v0, e, graph, PathFiller) {
                if(target(e, graph) != v1) {
                    continue;
                }
                if(graph[e].edgeId != edgeId) {
                    continue;
                }
                graph[e].markerIntervals[i].push_back({ordinal0, ordinal1});
                done = true;
                break;
            }

            // If we did not find this edge, we have to create it.
            if(not done) {
                edge_descriptor e;
                bool edgeWasAdded = false;
                tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
                SHASTA_ASSERT(edgeWasAdded);
                PathFillerEdge& edge = graph[e];
                edge.edgeId = edgeId;
                edge.markerIntervals.resize(orientedReadInfos.size());
                edge.markerIntervals[i].push_back({ordinal0, ordinal1});
            }

        }
    }
}



PathFillerVertex::PathFillerVertex(
    MarkerGraphVertexId vertexId,
    uint64_t orientedReadCount) :
    vertexId(vertexId),
    ordinals(orientedReadCount)
{
}



PathFiller::vertex_descriptor PathFiller::OrientedReadInfo::getVertex(uint32_t ordinal) const
{
    SHASTA_ASSERT(ordinal >= ordinalA0);
    SHASTA_ASSERT(ordinal <= ordinalB1);
    const uint32_t index = ordinal - ordinalA0;
    SHASTA_ASSERT(index < vertices.size());
    return vertices[index];

}



void PathFiller::writeGraphviz(
    ostream& out,
    bool showVertices,
    bool showVertexLabels,
    bool showEdgeLabels) const
{
    const PathFiller& graph = *this;

    const double S = 0.7;
    const double V = 1.;

    out <<
        "digraph PathFillerGraph {\n"
        "mclimit=0.01;\n"       // For layout speed
        "edge [penwidth=5];\n"
        "node [fontname=\"Courier New\"];\n"
        "edge [fontname=\"Courier New\"];\n";

    if(showVertices) {
        if(showVertexLabels) {
            out << "node [shape=rectangle];\n";
        } else {
            out << "node [shape=point width=0.2];\n";
        }
    } else {
        out << "node [shape=point style=invis];\n";
    }


    // To help Graphviz compute the layout, write vertices in rank order.
    vector< pair<vertex_descriptor, uint64_t> > verticesWithRank;
    BGL_FORALL_VERTICES(v, graph, PathFiller) {
        verticesWithRank.push_back({v, graph[v].rank});
    }
    sort(verticesWithRank.begin(), verticesWithRank.end(),
        OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());
    for(const auto& p: verticesWithRank) {
        const vertex_descriptor v = p.first;
        const PathFillerVertex& vertex = graph[v];
        out << "\"" << vertex.stringId() << "\" [";
        out << "tooltip=\"" << vertex.stringId() << "\\n" << vertex.coverage() << "\"";
        if(showVertexLabels) {
            out << " label=\"" << vertex.stringId() << "\\n" << vertex.coverage() << "\"";
        }
        if(vertex.isStrongComponentVertex()) {
            out << " color=purple";
        }
        out << "];\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, graph, PathFiller) {

        const PathFillerEdge& edge = graph[e];
        const uint64_t coverage = edge.coverage();
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);
        const auto edgeSequence = getEdgeSequence(e);

        // Compute the hue based on coverage.
        double H;
        if(coverage >= orientedReadInfos.size()) {
            H = 1./3.;
        } else {
            H = (double(coverage - 1) / (3. * double(orientedReadInfos.size() - 1)));
        }
        const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";

        out << "\"" << graph[v0].stringId() << "\"->\"" << graph[v1].stringId() << "\" [";

        // Color is based on coverage.
        out << " color=" << colorString;

        // Tooltip.
        out << " tooltip=\"";
        if(edge.isVirtual) {
            out << "Virtual";
        } else {
            out << edge.edgeId;
        }
        out << "\\n" << coverage << "\\n";
        copy(edgeSequence.begin(), edgeSequence.end(), ostream_iterator<Base>(out));
        out << "\"";

        // Label.
        if(showEdgeLabels) {
            out << " label=\"";
            if(edge.isVirtual) {
                out << "Virtual";
            } else {
                out << edge.edgeId;
            }
            out << "\\n" << coverage << "\\n";
            copy(edgeSequence.begin(), edgeSequence.end(), ostream_iterator<Base>(out));
            out << "\"";
        }

        // Virtual edges and edges with source or target in a strong component
        // are shown dashed.
        if(edge.isVirtual or
            graph[v0].isStrongComponentVertex() or
            graph[v1].isStrongComponentVertex()) {
            out << " style=dashed";
        }

        if(not edge.isDagEdge) {
            out << " constraint=false";
        }
        out << "];\n";

    }

    out << "}\n";
}



string PathFillerVertex::stringId() const
{
    string s = to_string(vertexId);
    if(replicaIndex) {
        s += "." + to_string(replicaIndex);
    }
    return s;
}



void PathFiller::writeVerticesCsv() const
{
    const PathFiller& graph = *this;

    ofstream csv("PathFiller-vertices.csv");
    csv << "Vertex,i,OrientedReadId,Offsets\n";
    BGL_FORALL_VERTICES(v, graph, PathFiller) {
        const PathFillerVertex& vertex = graph[v];
        for(uint64_t i=0; i<vertex.ordinals.size(); i++) {
            const vector<uint32_t>& ordinals = vertex.ordinals[i];
            const auto& info = orientedReadInfos[i];
            csv << vertex.stringId() << ",";
            csv << i << ",";
            csv << info.orientedReadId << ",";
            for(const uint32_t ordinal: ordinals) {
                const uint64_t offset = ordinal - info.ordinalA0;
                csv << offset << ",";
            }
            csv << "\n";
        }
    }
}



// Return the total number of ordinals.
uint64_t PathFillerVertex::coverage() const
{
    uint64_t c = 0;
    for(const auto& v: ordinals) {
        c += v.size();
    }
    return c;
}



void PathFiller::writeGraph(
    ostream& html,
    bool showVertices,
    bool showVertexLabels,
    bool showEdgeLabels) const
{
    // Write out the graph in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    {
        ofstream dotFile(dotFileName);
        writeGraphviz(dotFile, showVertices, showVertexLabels, showEdgeLabels);
    }

    // Compute layout in svg format.
    const string command = "dot -O -T svg " + dotFileName;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    const double timeout = 600;
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



void PathFiller::findAssemblyPath()
{
    const PathFiller& graph = *this;
    assemblyPath.clear();

    const bool debug = false;
    if(debug) {
        cout << timestamp << "PathFiller::findAssemblyPath begins." << endl;
    }

    // Find the first and last vertex of the path we are looking for.
    const OrientedReadInfo& firstInfo = orientedReadInfos.front();
    const vertex_descriptor vA = firstInfo.vertices.front();
    const vertex_descriptor vB = firstInfo.vertices.back();



    // Main iteration loop.
    vertex_descriptor v = vA;
    while(v != vB) {
        if(debug) {
            cout << "At vertex " << graph[v].stringId() << endl;
        }

        // Find the edge with the most coverage.
        edge_descriptor eNext;
        uint64_t bestCoverage = 0;
        BGL_FORALL_OUTEDGES(v, e, graph, PathFiller) {
            const uint64_t coverage = graph[e].coverage();
            if(coverage > bestCoverage) {
                eNext = e;
                bestCoverage = coverage;
            }
        }
        SHASTA_ASSERT(bestCoverage > 0);

        // Store this edge.
        assemblyPath.push_back(eNext);
        v = target(eNext, graph);
    }

    // As constructed, the assemblyPath includes edgeIdA and edgeIB.
    SHASTA_ASSERT(assemblyPath.size() >= 2);
    SHASTA_ASSERT(graph[assemblyPath.front()].edgeId == edgeIdA);
    SHASTA_ASSERT(graph[assemblyPath.back()].edgeId == edgeIdB);

    if(debug) {
        cout << timestamp << "PathFiller::findAssemblyPath ends." << endl;
    }
}



// Get the sequence.
// The sequences of edgeIdA and edgeIdB are only included if
// includePrimary is true.
void PathFiller::getSequence(
    vector<Base>& sequence,
    bool includePrimary) const
{
    sequence.clear();

#if 0
    if(includePrimary) {
        const auto sequenceA = assembler.markerGraph.edgeSequence[edgeIdA];
        copy(sequenceA.begin(), sequenceA.end(), back_inserter(sequence));
    }

    for(const MarkerGraphEdgeId edgeId: secondaryEdges) {
        const auto thisEdgeSequence = assembler.markerGraph.edgeSequence[edgeId];
        copy(thisEdgeSequence.begin(), thisEdgeSequence.end(), back_inserter(sequence));
    }

    if(includePrimary) {
        const auto sequenceB = assembler.markerGraph.edgeSequence[edgeIdB];
        copy(sequenceB.begin(), sequenceB.end(), back_inserter(sequence));
    }
#endif


    // Use the assemblyPath to assemble sequence.
    for(uint64_t i=0; i<assemblyPath.size(); i++) {

        // Skip edgeIdA and edgeIdB if so requested.
        if(not includePrimary) {
            if(i == 0 or i==assemblyPath.size() - 1) {
                continue;
            }
        }

        // Append the sequence contributed by this edge, virtual or not.
        const auto edgeSequence = getEdgeSequence(assemblyPath[i]);
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
    }
}

// Get the edge sequence from the marker graph, for a regular edge,
// or from the edge itself, for a virtual edge.
span<const Base> PathFiller::getEdgeSequence(edge_descriptor e) const
{
    const PathFiller& graph = *this;
    const PathFillerEdge& edge = graph[e];
    if(edge.isVirtual) {
        return span<const Base>(edge.sequence.begin(), edge.sequence.end());
    } else {
        return assembler.markerGraph.edgeSequence[edge.edgeId];
    }
}


void PathFiller::writeSequence(ostream& html) const
{
    if(not html) {
        return;
    }

    const auto sequenceA = assembler.markerGraph.edgeSequence[edgeIdA];
    const auto sequenceB = assembler.markerGraph.edgeSequence[edgeIdB];
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



void PathFiller::writeSequenceFasta(ostream& fasta) const
{
    vector<Base> sequence;
    getSequence(sequence, true);

    fasta << ">Path " << sequence.size() << "\n";
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
    fasta << endl;

}



void PathFiller::writeAssemblyDetails(ostream& csv) const
{
    const PathFiller& graph = *this;

    csv << "Index,Begin,End,VertexId0,VertexId1,EdgeId,Sequence\n";

    uint64_t assembledPosition = 0;
    for(uint64_t i=0; i<assemblyPath.size(); i++) {
        const edge_descriptor e = assemblyPath[i];
        const PathFillerEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const auto edgeSequence = getEdgeSequence(e);

        csv << i << ",";
        csv << assembledPosition << ",";
        csv << assembledPosition + edgeSequence.size() << ",";
        csv << graph[v0].vertexId << ",";
        csv << graph[v1].vertexId << ",";
        if(not edge.isVirtual) {
            csv << edge.edgeId;
        }
        csv << ",";
        copy(edgeSequence.begin(), edgeSequence.end(), ostream_iterator<Base>(csv));
        csv << "\n";

        assembledPosition += edgeSequence.size();
    }
}




void PathFiller::computeStrongComponents()
{
    PathFiller& graph = *this;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, graph, PathFiller) {
        vertexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        graph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }



    // Keep the non-trivial ones.
    // A non-trivial strong component has at least one internal edge.
    // This means that it either has more than one vertex,
    // or it consists of a single vertex with a self-edge.
    strongComponents.clear();
    for(const auto& p: componentVertices) {

        // Figure out if it is non-trivial.
        bool isNonTrivial;
        if(p.second.size() > 1) {

            // More than one vertex. Certainly non-trivial.
            isNonTrivial = true;
        } else if (p.second.size() == 1) {

            // Only one vertex. Non-trivial if self-edge present.
            const vertex_descriptor v = p.second.front();
            bool selfEdgeExists = false;
            tie(ignore, selfEdgeExists) = edge(v, v, graph);
            isNonTrivial = selfEdgeExists;
        } else {

            // Empty. This should never happen.
            SHASTA_ASSERT(0);
        }

        // If non-trivial, store it.
        if(isNonTrivial) {
            strongComponents.push_back(StrongComponent(p.second));
        }
    }

    // Mark the vertices in the non-trivial strongly connected components.
    for(uint64_t strongComponentId=0; strongComponentId<strongComponents.size(); strongComponentId++) {
        const auto& strongComponent = strongComponents[strongComponentId];
        for(const vertex_descriptor v: strongComponent.vertices) {
            graph[v].strongComponentId = strongComponentId;
        }
    }

}



PathFiller::StrongComponent::StrongComponent(
    const vector<vertex_descriptor>& vertices) :
    vertices(vertices)
{}



void PathFiller::createVirtualEdges()
{
    const bool debug = false;
    PathFiller& graph = *this;
    const uint64_t coverage = orientedReadInfos.size();

    // Gather all the ordinals that will go away when we remove
    // the strong component vertices.
    vector< vector<uint32_t> > ordinalsToBeRemoved(coverage);
    BGL_FORALL_VERTICES(v, graph, PathFiller) {
        const PathFillerVertex& vertex = graph[v];
        if(not vertex.isStrongComponentVertex()) {
            continue;
        }
        for(uint64_t i=0; i<coverage; i++) {
            copy(vertex.ordinals[i].begin(), vertex.ordinals[i].end(),
                back_inserter(ordinalsToBeRemoved[i]));
        }
    }


    // For each oriented read, find consecutive streaks of ordinals
    // that will be removed. Each of these streaks will contribute to one
    // virtual edge.
    std::map<pair<vertex_descriptor, vertex_descriptor>, edge_descriptor> virtualEdges;
    for(uint64_t i=0; i<coverage; i++) {
        vector<uint32_t>& o = ordinalsToBeRemoved[i];
        sort(o.begin(), o.end());

        if(debug) {
            cout << "Ordinals to be removed for " << orientedReadInfos[i].orientedReadId << ":" << endl;
            copy(o.begin(), o.end(), ostream_iterator<uint32_t>(cout, " "));
            cout << endl;
        }

        for(uint64_t j=0; j<o.size(); /* Increment later */) {
            const uint64_t streakBegin = j;
            uint64_t streakEnd = j;
            while(streakEnd < o.size() and o[streakEnd+1] == o[streakEnd] + 1) {
                ++streakEnd;
            }
            if(debug) {
                cout << "Streak " << streakBegin << " " << streakEnd << " ordinals: " <<
                    o[streakBegin] << " to " << o[streakEnd] << endl;
            }
            j = streakEnd + 1;

            // We will remove all the vertices at these ordinals.
            // Therefore this oriented read needs a virtual edge
            // with source one position earlier than streakBegin and
            // target one position past streakEnd.
            const uint32_t ordinal0 = o[streakBegin] - 1;
            const uint32_t ordinal1 = o[streakEnd] + 1;
            if(debug) {
                cout << "ordinal0 " << ordinal0 << ", ordinal1 " << ordinal1 << endl;
            }

            // The vertices corresponding to these ordinals are
            // the source and target of the virtual edge.
            const vertex_descriptor v0 = orientedReadInfos[i].getVertex(ordinal0);
            const vertex_descriptor v1 = orientedReadInfos[i].getVertex(ordinal1);

            // Locate this virtual edge in our map, creating it if necessary.
            auto it = virtualEdges.find({v0, v1});
            edge_descriptor e;
            if(it == virtualEdges.end()) {
                bool edgeWasAdded = false;
                tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
                SHASTA_ASSERT(edgeWasAdded);
                virtualEdges.insert({{v0, v1}, e});
                PathFillerEdge& edge = graph[e];
                edge.isVirtual = true;
                edge.markerIntervals.resize(coverage);
            } else {
                e = it->second;
            }
            graph[e].markerIntervals[i].push_back({ordinal0, ordinal1});
        }
    }
    if(debug) {
        cout << "Added " << virtualEdges.size() << " virtual edges." << endl;
    }

    // Remove all the vertices in strong components ands the
    // edges that have them as source or target.
    // Those edges were replaced by the virtual edges.
    for(const StrongComponent& strongComponent: strongComponents) {
        for(const vertex_descriptor v: strongComponent.vertices) {
            clear_vertex(v, graph);
            remove_vertex(v, graph);
        }
    }
    strongComponents.clear();

    // Assemble the virtual edges.
    // This uses MSA so compute optimal MarkerGraphEdgeIds
    // for each virtual edge.
    for(const auto& p: virtualEdges) {
        assembleVirtualEdge(p.second);
    }

#if 0
    // Construct the path of each oriented read in this strong components.
    // Each oriented read can only enter and exit this strong component once.
    // (If it did exit and then reenter it, the vertices in-between
    // would also be part of the strong component).
    vector< vector<edge_descriptor> > paths(orientedReadInfos.size());
    vector< vector<vertex_descriptor> > pathVertices(orientedReadInfos.size());
    vector<uint32_t> pathOrdinals;
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        // Gather the ordinals.
        pathOrdinals.clear();
        for(const vertex_descriptor v: strongComponent.vertices) {
            const PathFillerVertex& vertex = graph[v];
            copy(vertex.ordinals[i].begin(), vertex.ordinals[i].end(),
                back_inserter(pathOrdinals));
        }

        // If this oriented read is not present in this strong component,
        // it does not contribute to this virtual edge.
        if(pathOrdinals.empty()) {
            continue;
        }

        // If this oriented read only touches this strong component
        // at a single vertex, it does not contribute to this virtual edge.
        if(pathOrdinals.size() == 1) {
            continue;
        }

        // Sort the ordinals and verify that they are all contiguous.
        sort(pathOrdinals.begin(), pathOrdinals.end());
        for(uint64_t i=1; i<pathOrdinals.size(); i++) {
            SHASTA_ASSERT(pathOrdinals[i] == pathOrdinals[i-1] + 1);
        }

        // Sanity checks and debug output.
        {
            // The first and last ordinals are the ones where
            // this oriented read enters and exits this strong component.
            const uint32_t ordinal0 = pathOrdinals.front();
            const uint32_t ordinal1 = pathOrdinals.back();

            // Find the corresponding vertices.
            const vertex_descriptor v0 = info.getVertex(ordinal0);
            const vertex_descriptor v1 = info.getVertex(ordinal1);

            if(debug) {
                cout << info.orientedReadId << ":" << endl;
                cout << "    Enters at " << ordinal0 << " " << graph[v0].stringId() << endl;
                cout << "    Exits  at " << ordinal1 << " " << graph[v1].stringId() << endl;
            }
            SHASTA_ASSERT(graph[v0].isStrongComponentEntrance);
            SHASTA_ASSERT(graph[v1].isStrongComponentExit);
        }

        // Construct the path vertices.
        for(const uint32_t ordinal: pathOrdinals) {
            pathVertices[i].push_back(info.getVertex(ordinal));
        }
        if(debug) {
            cout << "Path vertices:" << endl;
            for(const vertex_descriptor v: pathVertices[i]) {
                cout << graph[v].stringId() << " ";
            }
            cout << endl;
        }

        // Construct the path edges.
        vector<edge_descriptor>& path = paths[i];
        for(uint64_t j=1; j<pathVertices[i].size(); j++) {
            const uint32_t ordinal0 = pathOrdinals[j-1];
            const uint32_t ordinal1 = pathOrdinals[j];
            const vertex_descriptor v0 = pathVertices[i][j-1];
            const vertex_descriptor v1 = pathVertices[i][j];
            const edge_descriptor e = findEdge(v0, v1, i, ordinal0, ordinal1);
            path.push_back(e);
        }
        if(debug) {
            cout << "Path edges:" << endl;
            for(const edge_descriptor e: path) {
                cout << graph[e].edgeId << " ";
            }
            cout << endl;
        }
    }



    // Gather oriented read paths for each entrance/exit combination.
    // Each of these combinations will generate a virtual edge.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<uint64_t> > pathMap;
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        pathMap[{pathVertices[i].front(), pathVertices[i].back()}].push_back(i);
    }
    if(debug) {
        cout << "Will generate " << pathMap.size() << " virtual edges." << endl;
        for(const auto& p: pathMap) {
            const vertex_descriptor entrance = p.first.first;
            const vertex_descriptor exit = p.first.second;
            cout << graph[entrance].stringId() << " to ";
            cout << graph[exit].stringId();
            cout << " coverage " << p.second.size() << endl;
        }
    }


    // Generate the virtual edges.
    for(const auto& p: pathMap) {
        const vector<uint64_t>& orientedReadIndexes = p.second;
        VirtualEdge virtualEdge;
        virtualEdge.strongComponentId = strongComponentId;
        virtualEdge.entrance = p.first.first;
        virtualEdge.exit = p.first.second;
        virtualEdge.coverage = orientedReadIndexes.size();
        graph[virtualEdge.entrance].virtualEdgeIndexes.push_back(virtualEdges.size());
        virtualEdges.push_back(virtualEdge);

        // Gather the sequence of MargerGraphEdgeIds seen by each of the
        // oriented reads that contribute to this virtual edge.
        // SeqAN uses 0 to represent gaps, so we add 1
        // to the edge ids.
        vector< vector<MarkerGraphEdgeId> > sequences;
        for(uint64_t i: orientedReadIndexes) {
            sequences.resize(sequences.size() + 1);
            for(const edge_descriptor e: paths[i]) {
                sequences.back().push_back(graph[e].edgeId + 1);
            }
            cout << endl;
        }

        if(debug) {
            const vertex_descriptor entrance = p.first.first;
            const vertex_descriptor exit = p.first.second;
            cout << "Virtual edge ";
            cout << graph[entrance].stringId() << " to ";
            cout << graph[exit].stringId();
            cout << " coverage " << p.second.size() << endl;

            for(const auto& v: sequences) {
                for(const auto value: v) {
                    cout << " " << value - 1;
                }
                cout << endl;
            }
        }

        // Use Seqan to do an MSA of these sequences.
        {
            using namespace seqan;
            Align< String<MarkerGraphEdgeId> > align;
            const uint64_t n = sequences.size();
            resize(rows(align), n);
            for(uint64_t i=0; i<n; i++) {
                assignSource(row(align, i), sequences[i]);
            }
            globalMsaAlignment(align, Score<int64_t, Simple>(0, -1, -4));

            if(debug) {
                cout << "MSA:" << endl;
                for(uint64_t rowNumber=0; rowNumber<n; rowNumber++) {
                    for(uint64_t j=0; j<length(row(align, rowNumber)); j++) {
                        const auto value = getValue(row(align, rowNumber), j);
                        if(value != 0) {
                            cout << value - 1;
                        }
                        cout << ",";
                    }
                    cout << endl;
                }
            }

            // Compute the consensus.
            vector<uint64_t> values;
            vector<uint64_t> frequencies;
            vector<uint64_t> consensusValues;
            const uint64_t alignmentLength = length(row(align, 0));
            for(uint64_t j=0; j<alignmentLength; j++) {
                values.clear();
                frequencies.clear();
                for(uint64_t rowNumber=0; rowNumber<n; rowNumber++) {
                    const auto value = getValue(row(align, rowNumber), j);
                    values.push_back(value);
                }
                deduplicateAndCount(values, frequencies);

                // Find the value with the highest frequency.
                const uint64_t maxFrequencyIndex =
                    max_element(frequencies.begin(), frequencies.end()) - frequencies.begin();
                const uint64_t consensusValue = values[maxFrequencyIndex];

                // Store it.
                consensusValues.push_back(consensusValue);
                if(consensusValue > 0) {
                    virtualEdges.back().markerGraphEdges.push_back(consensusValue - 1);
                }
            }

            if(debug) {
                cout << "Aligned consensus:" << endl;
                for(const uint64_t value: consensusValues) {
                    if(value) {
                        cout << value - 1;
                    }
                    cout << ",";
                }
                cout << endl;
                cout << "VirtualEdge consensus:" << endl;
                for(const MarkerGraphEdgeId edgeId: virtualEdges.back().markerGraphEdges) {
                    cout << edgeId << " ";
                }
                cout << endl;
            }
        }
    }
#endif
}



// Assemble a virtual edge.
// This uses MSA so compute optimal MarkerGraphEdgeIds.
void PathFiller::assembleVirtualEdge(edge_descriptor e)
{
    const bool debug = false;
    PathFiller& graph = *this;
    PathFillerEdge& edge = graph[e];
    SHASTA_ASSERT(edge.isVirtual);


    // Gather the marker graph paths of the contributing oriented reads.
    // This contains MarkerGraphEdgeIds incremented by 1.
    // We need this because SeqAn uses 0 to represent gaps.
    vector< vector<MarkerGraphEdgeId> > sequences;

    // Loop over oriented reads.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadInfos[i].orientedReadId;
        const vector< pair<uint32_t, uint32_t> >& orientedReadMarkerIntervals = edge.markerIntervals[i];

        // Loop over MarkerIntervals of this oriented read.
        for(const pair<uint32_t, uint32_t>& markerInterval: orientedReadMarkerIntervals) {
            const uint32_t ordinal0 = markerInterval.first;
            const uint32_t ordinal1 = markerInterval.second;

            sequences.resize(sequences.size() + 1);
            vector<MarkerGraphEdgeId>& sequence = sequences.back();

            // Follow the marker graph edges in this range.
            for(uint32_t ordinal=ordinal0; ordinal<ordinal1; ordinal++) {
                const MarkerGraphEdgeId edgeId = assembler.markerGraph.locateMarkerInterval(assembler.markers,
                    MarkerInterval(orientedReadId, ordinal, ordinal+1));
                SHASTA_ASSERT(edgeId != invalid<MarkerGraphEdgeId>);
                sequence.push_back(edgeId + 1);
            }
        }
    }

    if(debug) {
        cout << "Virtual edge ";
        cout << graph[source(e, graph)].stringId() << " to ";
        cout << graph[target(e, graph)].stringId();
        cout << " coverage " << edge.coverage() << endl;

        for(const auto& v: sequences) {
            for(const auto value: v) {
                cout << " " << value - 1;   // Subtract 1 to get the actual edgeId
            }
            cout << endl;
        }
    }

    // Use Seqan to do an MSA of these sequences.
    vector<MarkerGraphEdgeId> edgeIds;
    {
        using namespace seqan;
        Align< String<MarkerGraphEdgeId> > align;
        const uint64_t n = sequences.size();
        resize(rows(align), n);
        for(uint64_t i=0; i<n; i++) {
            assignSource(row(align, i), sequences[i]);
        }
        globalMsaAlignment(align, Score<int64_t, Simple>(0, -1, -4));

        if(debug) {
            cout << "MSA:" << endl;
            for(uint64_t rowNumber=0; rowNumber<n; rowNumber++) {
                for(uint64_t j=0; j<length(row(align, rowNumber)); j++) {
                    const auto value = getValue(row(align, rowNumber), j);
                    if(value != 0) {
                        cout << value - 1; // Subtract 1 to get the actual edgeId
                    }
                    cout << ",";
                }
                cout << endl;
            }
        }

        // Compute the consensus.
        vector<uint64_t> values;
        vector<uint64_t> frequencies;
        vector<uint64_t> consensusValues;
        const uint64_t alignmentLength = length(row(align, 0));
        for(uint64_t j=0; j<alignmentLength; j++) {
            values.clear();
            frequencies.clear();
            for(uint64_t rowNumber=0; rowNumber<n; rowNumber++) {
                const auto value = getValue(row(align, rowNumber), j);
                values.push_back(value);
            }
            deduplicateAndCount(values, frequencies);

            // Find the value with the highest frequency.
            const uint64_t maxFrequencyIndex =
                max_element(frequencies.begin(), frequencies.end()) - frequencies.begin();
            const uint64_t consensusValue = values[maxFrequencyIndex];

            // Store it.
            consensusValues.push_back(consensusValue);
            if(consensusValue > 0) {
                edgeIds.push_back(consensusValue - 1);
            }
        }

        if(debug) {
            cout << "Aligned consensus:" << endl;
            for(const uint64_t value: consensusValues) {
                if(value) {
                    cout << value - 1;
                }
                cout << ",";
            }
            cout << endl;
            cout << "VirtualEdge consensus:" << endl;
            for(const MarkerGraphEdgeId edgeId: edgeIds) {
                cout << edgeId << " ";
            }
            cout << endl;
        }
    }

    // Store the sequence.
    // It is obtained from the consensus edgeIds.
    edge.sequence.clear();
    for(const MarkerGraphEdgeId edgeId: edgeIds) {
        const auto thisEdgeSequence = assembler.markerGraph.edgeSequence[edgeId];
        copy(thisEdgeSequence.begin(), thisEdgeSequence.end(),
            back_inserter(edge.sequence));
    }

}



#if 0
// Find the edge v0->v1 that contains the specified MarkerInterval
// for the i-th oriented read.
pair<PathFiller::edge_descriptor, bool> PathFiller::findEdge(
    vertex_descriptor v0,
    vertex_descriptor v1,
    uint64_t i,
    uint32_t ordinal0,
    uint32_t ordinal1) const
{
    const PathFiller& graph = *this;

    BGL_FORALL_OUTEDGES(v0, e, graph, PathFiller) {
        if(target(e, graph) != v1) {
            continue;
        }
        const PathFillerEdge& edge = graph[e];
        const auto& markerIntervals = edge.markerIntervals[i];
        if(find(markerIntervals.begin(), markerIntervals.end(),
            make_pair(ordinal0, ordinal1)) != markerIntervals.end()) {
            return make_pair(e, true);
        }
    }

    return make_pair(edge_descriptor(), false);
}
#endif

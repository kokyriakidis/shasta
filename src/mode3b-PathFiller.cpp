// Shasta.
#include "mode3b-PathFiller.hpp"
#include "approximateTopologicalSort.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "platformDependent.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
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
    if(html) {
        html << "<br>There are " << virtualEdges.size() << " virtual edges.";
    }
    writeVerticesCsv();
    approximateTopologicalSort();
    if(html) {
        html << "<p>The local marker graph for this step assembly has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
        uint64_t nonDagEdgeCount = 0;
        BGL_FORALL_EDGES(e, graph, PathFiller) {
            if(not graph[e].isDagEdge) {
                ++nonDagEdgeCount;
            }
        }
        writeStrongComponents(html);
        writeGraph(html);
    }

#if 0
    writeVerticesCsv();

    const bool success = fillPathGreedy(html);

    if(html) {
        html << "<p>The local marker graph for this step assembly has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
        writeGraph(html);

        if(success) {
            writeSequence(html);
        } else {
            cout << "Unable to fill assembly path between primary edges " <<
            edgeIdA << " " << edgeIdB << endl;
        }
    }

    if(not success) {
        throw runtime_error("Unable to fill assembly path between primary edges " +
            to_string(edgeIdA) + " " + to_string(edgeIdB));
    }
#endif
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



// Return true if any oriented reads have more than one ordinal on this vertex.
bool PathFillerVertex::hasDuplicateOrientedReads() const
{
    for(const auto& v: ordinals) {
        if(v.size() > 1) {
            return true;
        }
    }
    return false;
}



// Return true if any oriented reads have more than one
// MarkerInterval on this edge.
bool PathFillerEdge::hasDuplicateOrientedReads() const
{
    for(const auto& v: markerIntervals) {
        if(v.size() > 1) {
            return true;
        }
    }
    return false;
}



void PathFiller::writeGraphviz(ostream& out, bool showVirtualEdges) const
{
    const PathFiller& graph = *this;

    const double S = 0.7;
    const double V = 1.;

    out <<
        "digraph PathFillerGraph {\n"
        "mclimit=0.01;\n"       // For layout speed
        // "maxiter=1;\n"       // For layout speed (?)
        "node [shape=point style=invis];\n"
        "edge [penwidth=5];\n";



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
        out << "\"" << vertex.stringId() << "\"";
        if(vertex.isStrongComponentVertex()) {
            out << "[style=solid width=0.2;";
            if(vertex.isStrongComponentEntrance) {
                if(vertex.isStrongComponentExit) {
                    out << " color=Magenta"; // Entrance and exit.
                } else {
                    out << " color=Cyan";   // Entrance only.
                }
            } else if(vertex.isStrongComponentExit) {
                out << " color=Gold";       // Exit only.
            }
            out << "]"; // Override invisible style set above.
        }
        out << ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, graph, PathFiller) {

        const PathFillerEdge& edge = graph[e];
        if(showVirtualEdges and isStrongComponentEdge(e)) {
            // Skip it.
            continue;
        }
        const uint64_t coverage = edge.coverage();
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);

        // Compute the hue based on coverage.
        double H;
        if(coverage >= orientedReadInfos.size()) {
            H = 1./3.;
        } else {
            H = (double(coverage - 1) / (3. * double(orientedReadInfos.size() - 1)));
        }
        const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";

        out << "\"" << graph[v0].stringId() << "\"->\"" << graph[v1].stringId() << "\" [";
        out << " color=" << colorString;
        out << " tooltip=\"Coverage " << coverage << "\"";
        if(isStrongComponentEdge(e)) {
            out << " style=dashed";
        }
        if(not edge.isDagEdge) {
            out << " constraint=false";
        }
        out << "];\n";

    }



    // Write virtual edges, if requested.
    if(showVirtualEdges) {
        for(const VirtualEdge& virtualEdge: virtualEdges) {

            // Compute the hue based on coverage.
            const uint64_t coverage = virtualEdge.coverage;
            double H;
            if(coverage >= orientedReadInfos.size()) {
                H = 1./3.;
            } else {
                H = (double(coverage - 1) / (3. * double(orientedReadInfos.size() - 1)));
            }
            const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";

            out << "\"" << graph[virtualEdge.entrance].stringId() << "\"->\"" <<
                graph[virtualEdge.exit].stringId() << "\" [";
            out << " color=" << colorString;
            out << " tooltip=\"Coverage " << coverage << "\"";
            out << " style=dashed";
            out << "];\n";
        }
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



// Return the total number of marker intervals.
uint64_t PathFillerEdge::coverage() const
{
    uint64_t c = 0;
    for(const auto& v: markerIntervals) {
        c += v.size();
    }
    return c;
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



void PathFiller::writeGraph(ostream& html) const
{
    // Write out the graph in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    {
        ofstream dotFile(dotFileName);
        const bool showVirtualEdges = true;
        writeGraphviz(dotFile, showVirtualEdges);
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



void PathFiller::writeStrongComponents(ostream& html) const
{
    html << "<br>Found " << strongComponents.size() <<
        " non-trivial strongly connected components.";
    for(uint64_t strongComponentId=0; strongComponentId<strongComponents.size(); strongComponentId++) {
        writeStrongComponent(html, strongComponentId);
    }

}



void PathFiller::writeStrongComponent(ostream& html, uint64_t strongComponentId) const
{
    const StrongComponent& strongComponent = strongComponents[strongComponentId];
    html << "<br>Strong component " << strongComponentId << " with " <<
        strongComponent.vertices.size() << " vertices.";

    // Write out this strong component in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    {
        ofstream dotFile(dotFileName);
        writeStrongComponentGraphviz(dotFile, strongComponentId);
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



void PathFiller::writeStrongComponentGraphviz(
    ostream& out,
    uint64_t strongComponentId) const
{
    const PathFiller& graph = *this;

    const double S = 0.7;
    const double V = 1.;

    const StrongComponent& strongComponent = strongComponents[strongComponentId];
    out << "digraph StrongComponent" << strongComponentId << "{\n"
        "node [shape=rectangle style=filled fontname=\"Courier New\"];\n"
        "edge [penwidth=5 fontname=\"Courier New\"];\n";

    // For consistency with the display of the entire graph,
    // write the vertices in rank order.
    vector< pair<vertex_descriptor, uint64_t> > verticesWithRank;
    for(const vertex_descriptor v: strongComponent.vertices) {
        verticesWithRank.push_back({v, graph[v].rank});
    }
    sort(verticesWithRank.begin(), verticesWithRank.end(),
        OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());
    for(const auto& p: verticesWithRank) {
        const vertex_descriptor v = p.first;
        const PathFillerVertex& vertex = graph[v];
        SHASTA_ASSERT(vertex.strongComponentId == strongComponentId);
        out << "\"" << vertex.stringId() << "\"";
        out << " [";

        // Label.
        out << "label=\"" << vertex.stringId() << "\\n" << vertex.coverage();
        if(vertex.isStrongComponentEntrance) {
            out << "\\nEntrance";
        }
        if(vertex.isStrongComponentExit) {
            out << "\\nExit";
        }
        for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
            const auto& ordinals = vertex.ordinals[i];
            if(not ordinals.empty()) {
                out << "\\n" << i << " " << orientedReadInfos[i].orientedReadId;
                for(const uint32_t ordinal: ordinals) {
                    out << " " << ordinal;
                }
            }
        }
        out << "\"";

        // Color.
        if(vertex.isStrongComponentEntrance) {
            if(vertex.isStrongComponentExit) {
                out << " fillcolor=Magenta"; // Entrance and exit.
            } else {
                out << " fillcolor=Cyan";   // Entrance only.
            }
        } else if(vertex.isStrongComponentExit) {
            out << " fillcolor=Gold";       // Exit only.
        }

        out << "]";
        out << ";\n";
    }

    // Write the edges.
    for(const vertex_descriptor v0: strongComponent.vertices) {
        const PathFillerVertex& vertex0 = graph[v0];
        SHASTA_ASSERT(vertex0.strongComponentId == strongComponentId);

        BGL_FORALL_OUTEDGES(v0, e, graph, PathFiller) {
            const PathFillerEdge& edge = graph[e];
            const uint64_t coverage = edge.coverage();
            const vertex_descriptor v1 = target(e, graph);
            const PathFillerVertex& vertex1 = graph[v1];
            if(vertex1.strongComponentId != strongComponentId) {
                continue;
            }

            // Compute the hue based on coverage.
            double H;
            if(coverage >= orientedReadInfos.size()) {
                H = 1./3.;
            } else {
                H = (double(coverage - 1) / (3. * double(orientedReadInfos.size() - 1)));
            }
            const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";

            out <<
                "\"" << vertex0.stringId() << "\""
                "->"
                "\"" << vertex1.stringId() << "\"";
            out << "[";
            out << " color=" << colorString;
            out << " label=\"" << edge.edgeId << "\\n" << coverage << "\"";

            if(not edge.isDagEdge) {
                out << " constraint=false";
            }

            out << "]";
            out << ";\n";
        }

    }


    out << "}\n";
}



#if 0
bool PathFiller::fillPathGreedy(ostream& html)
{
    const bool debug = true;

    using vertex_descriptor = Graph::vertex_descriptor;
    using edge_descriptor = Graph::edge_descriptor;
    const MarkerGraph::Edge markerGraphEdgeA = assembler.markerGraph.edges[edgeIdA];
    const MarkerGraph::Edge markerGraphEdgeB = assembler.markerGraph.edges[edgeIdB];
    const Graph::edge_descriptor eA = graph.edgeMap[edgeIdA];
    // const Graph::edge_descriptor eB = graph.edgeMap[edgeIdB];
    const Edge& edgeA = graph[eA];
    // const Edge& edgeB = graph[eB];

    // The path we are looking for begins at vA and ends at vB.
    const vertex_descriptor vA = graph.vertexMap[markerGraphEdgeA.target];
    const vertex_descriptor vB = graph.vertexMap[markerGraphEdgeB.target];

    // The last MarkerInterval used for each of the oriented reads.
    // We will use this to make sure we only move forward.
    vector<uint32_t> lastMarkerInterval;
    SHASTA_ASSERT(not edgeA.hasCycle());    // True by construction.
    for(const auto& v: edgeA.markerIntervals) {
        SHASTA_ASSERT(v.size() == 1);
        lastMarkerInterval.push_back(v.front().first);
    }



    // Main iteration loop.
    vertex_descriptor v = vA;
    while(v != vB) {

        if(debug) {
            cout << "At vertex " << graph[v].vertexId << "\n";
            for(uint64_t i=0; i<lastMarkerInterval.size(); i++) {
                cout << i << " " << lastMarkerInterval[i] << "\n";
            }
        }

        // For each of our oriented reads,
        // find which of the possible out-edges it wants to go to,
        // if any, based on the ordinals stored in lastMarkerInterval.
        vector< pair<edge_descriptor, uint32_t> > nextEdgesAndOrdinals(orientedReadInfos.size());
        for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
            // There is no null_edge, so this can be uninitialized and it's ok,
            // but we have to suppress the warning.
            edge_descriptor eNext;
#pragma GCC diagnostic pop
            uint32_t nextOrdinal0 = std::numeric_limits<uint32_t>::max();
            BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
                const auto& markerIntervals = graph[e].markerIntervals[i];
                for(const auto& p: markerIntervals) {
                    const uint32_t ordinal0 = p.first;
                    if(ordinal0 > lastMarkerInterval[i]) {
                        if(ordinal0 < nextOrdinal0) {
                            eNext = e;
                            nextOrdinal0 = ordinal0;
                        }
                    }
                }
            }
            nextEdgesAndOrdinals[i] = {eNext, nextOrdinal0};
        }


        // Pick the next edge that most oriented reads prefer.
        vector<edge_descriptor> nextEdges;
        vector<uint64_t> nextEdgesFrequency;
        for(const auto& p: nextEdgesAndOrdinals) {
            if(p.second != std::numeric_limits<uint32_t>::max()) {
                nextEdges.push_back(p.first);
            }
        }
        SHASTA_ASSERT(not nextEdges.empty());
        deduplicateAndCount(nextEdges, nextEdgesFrequency);
        auto it = max_element(nextEdgesFrequency.begin(), nextEdgesFrequency.end());
        edge_descriptor eNext = nextEdges[it - nextEdgesFrequency.begin()];

        if(debug) {
            cout << "Next edge " << graph[eNext].edgeId << "\n";
        }


        // This edge gets added to the path.
        secondaryEdges.push_back(graph[eNext].edgeId);
        v = target(eNext, graph);

        // Update the lastMarkerInterval vector.
        for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
            const auto& p = nextEdgesAndOrdinals[i];
            if(eNext == p.first) {
                lastMarkerInterval[i] = p.second;
            }
        }
    }
    cout << flush;

    if(html) {
        html << "<p>Path secondary edges:";
        for(const MarkerGraphEdgeId edgeId: secondaryEdges) {
            html << " " << edgeId;
        }
    }

    // As constructed, pathSecondaryEdges includes edgeIB.
    // Check that this is the case, then remove it.
    SHASTA_ASSERT(secondaryEdges.size() >= 1);
    // SHASTA_ASSERT(secondaryEdges.front() == edgeIdA);
    SHASTA_ASSERT(secondaryEdges.back() == edgeIdB);
    secondaryEdges.resize(secondaryEdges.size() - 1);
    // secondaryEdges.erase(secondaryEdges.begin());

    /*
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
            return false;
        }

        edgesUsed.insert(eNext);
        secondaryEdges.push_back(graph[eNext].edgeId);
        v = target(eNext, graph);
    }

    // As constructed, pathSecondaryEdges includes edgeIdA and edgeIB.
    // Check that this is the case, then remove them.
    SHASTA_ASSERT(secondaryEdges.size() >= 2);
    SHASTA_ASSERT(secondaryEdges.front() == edgeIdA);
    SHASTA_ASSERT(secondaryEdges.back() == edgeIdB);
    secondaryEdges.resize(secondaryEdges.size() - 1);
    secondaryEdges.erase(secondaryEdges.begin());

    */

    return true;
}
#endif


span<const shasta::Base> PathFiller::edgeSequence(
    edge_descriptor e) const
{
    return edgeSequence((*this)[e].edgeId);
}



span<const shasta::Base> PathFiller::edgeSequence(
    MarkerGraphEdgeId edgeId) const
{
    return assembler.markerGraph.edgeSequence[edgeId];
}



#if 0
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

    for(const MarkerGraphEdgeId edgeId: secondaryEdges) {
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
#endif



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
    strongComponents.clear();
    for(const auto& p: componentVertices) {
        if(p.second.size() > 1) {
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

    // Mark the entrances of each strong component.
    for(StrongComponent& strongComponent: strongComponents) {
        for(const vertex_descriptor v0: strongComponent.vertices) {
            BGL_FORALL_INEDGES(v0, e, graph, PathFiller) {
                const vertex_descriptor v1 = source(e, graph);
                if(graph[v1].strongComponentId != graph[v0].strongComponentId) {
                    strongComponent.entrances.push_back(v0);
                    graph[v0].isStrongComponentEntrance = true;
                    break;
                }
            }
        }
    }

    // Mark the exits of each strong component.
    for(StrongComponent& strongComponent: strongComponents) {
        for(const vertex_descriptor v0: strongComponent.vertices) {
            BGL_FORALL_OUTEDGES(v0, e, graph, PathFiller) {
                const vertex_descriptor v1 = target(e, graph);
                if(graph[v1].strongComponentId != graph[v0].strongComponentId) {
                    strongComponent.exits.push_back(v0);
                    graph[v0].isStrongComponentExit = true;
                    break;
                }
            }
        }
    }
}



// This returns true if the specified edge is internal to a
// strongly connected component.
bool PathFiller::isStrongComponentEdge(edge_descriptor e) const
{
    const PathFiller& graph = *this;
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const uint64_t c0 = graph[v0].strongComponentId;
    const uint64_t c1 = graph[v1].strongComponentId;
    return
        (c0 != invalid<uint64_t>) and
        (c1 != invalid<uint64_t>) and
        (c0 == c1);
}



PathFiller::StrongComponent::StrongComponent(
    const vector<vertex_descriptor>& vertices) :
    vertices(vertices)
{}



void PathFiller::createVirtualEdges()
{
    for(uint64_t strongComponentId=0; strongComponentId<strongComponents.size(); strongComponentId++) {
        createVirtualEdges(strongComponentId);
    }
}



void PathFiller::createVirtualEdges(uint64_t strongComponentId)
{
    const bool debug = true;
    const PathFiller& graph = *this;
    const StrongComponent& strongComponent = strongComponents[strongComponentId];

    if(debug) {
        cout << "Creating virtual edges for strong component " << strongComponentId << endl;
    }



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
        cout << "Will generate " << pathMap.size() << " virtual edges:" << endl;
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
        VirtualEdge virtualEdge;
        virtualEdge.strongComponentId = strongComponentId;
        virtualEdge.entrance = p.first.first;
        virtualEdge.exit = p.first.second;
        virtualEdge.coverage = p.second.size();
        virtualEdges.push_back(virtualEdge);
    }

}



// Find the edge v0->v1 that contains the specified MarkerInterval
// for the i-th oriented read.
PathFiller::edge_descriptor PathFiller::findEdge(
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
            return e;
        }
    }

    SHASTA_ASSERT(0);
}


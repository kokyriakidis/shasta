// Shasta
#include "Assembler.hpp"
#include "LocalMarkerGraph1.hpp"
#include "mode3b-PathFiller.hpp"
#include "mode3b-PathFinder.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>



void Assembler::findCompleteMarkerGraphPath(
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction              // 0=forward, 1=backward
    ) const
{
    // Check our assumptions.
    SHASTA_ASSERT(assemblerInfo->assemblyMode == 3);    // Complete assembly graph.
    SHASTA_ASSERT(getReads().representation == 0);      // No RLE


    // Check that we have what we need.
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();

    // Do it.
    mode3b::PathFinder pathFinder(
        *this,
        startEdgeId,
        direction);
}



void Assembler::findCompleteMarkerGraphPaths() const
{
    mode3b::PathFinder pathFinder(*this);
}



#if 0
// Given two consecutive primary edges in the marker graph,
// find the secondary edges in between.
// This works with mode3b assembly and the complete marker graph.
// This assumes that the two edges have no duplicate oriented read ids.
void Assembler::fillMode3bAssemblyPathStep(
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    ostream& html
    ) const
{

    secondaryEdges.clear();

    // For each of the common oriented reads, gather the
    // edges in between. Each is stored with its coverage
    // from the common oriented reads.
    std::map<MarkerGraphEdgeId, uint64_t> edgeMap;

    // Joint loop over the MarkerIntervals of the two edges.
    const auto markerIntervalsA = markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = markerGraph.edgeMarkerIntervals[edgeIdB];
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
        const OrientedReadId orientedReadId = itA->orientedReadId;
        const uint32_t ordinalA0 = itA->ordinals[0];
        const uint32_t ordinalA1 = itA->ordinals[1];
        const uint32_t ordinalB0 = itB->ordinals[0];
        const uint32_t ordinalB1 = itB->ordinals[1];
        SHASTA_ASSERT(ordinalA1 == ordinalA0 + 1);
        SHASTA_ASSERT(ordinalB1 == ordinalB0 + 1);

        // If the edges don't appear in the expected order in this oriented read,
        // skip it.
        if(ordinalB0 < ordinalA1) {
            continue;
        }

        // If the edges are consecutive in this oriented read,
        // this oriented read contributes no intermediate edges.
        if(ordinalB0 == ordinalB1) {
            continue;
        }

        // This oriented read sees a gap in between these two edges.
        SHASTA_ASSERT(ordinalB0 > ordinalA1);

        // Loop over the edges in between seen by this oriented read.
        for(uint32_t ordinal0=ordinalA0; ordinal0!=ordinalB1; ordinal0++) {
            const uint32_t ordinal1 = ordinal0 + 1;
            MarkerInterval targetMarkerInterval(orientedReadId, ordinal0, ordinal1);
            const MarkerGraphEdgeId edgeId = markerGraph.locateMarkerInterval(markers, targetMarkerInterval);

            // Update the edge map.
            auto it = edgeMap.find(edgeId);
            if(it == edgeMap.end()) {
                edgeMap.insert({edgeId, 1});
            } else {
                ++it->second;
            }
        }

        // Continue the joint loop.
        ++itA;
        ++itB;
    }

    cout << "Edges with coverage from common reads:" << endl;
    for(const auto& p: edgeMap) {
        cout << p.first << " " << p.second << endl;
    }

    // Construct the LocalMarkerGraph1 with these edges.
    vector<MarkerGraphEdgeId> edgeVector;
    for(const auto& p: edgeMap) {
        edgeVector.push_back(p.first);
    }
    LocalMarkerGraph1 graph(markers, markerGraph, edgeVector);
    cout << "The local marker graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;


    // Write it out.
    {
        ofstream out("LocalMarkerGraph.dot");
        out << "digraph LocalMarkerGraph {\n";
        BGL_FORALL_EDGES(e, graph, LocalMarkerGraph1) {
            const MarkerGraphEdgeId edgeId = graph[e].edgeId;
            const uint64_t coverage = edgeMap[edgeId];
            const auto v0 = source(e, graph);
            const auto v1 = target(e, graph);
            const MarkerGraphVertexId vertexId0 = graph[v0].vertexId;
            const MarkerGraphVertexId vertexId1 = graph[v1].vertexId;
            out << vertexId0 << "->" << vertexId1 <<
                " [label=\"" << edgeId << "\\n" << coverage << "\"";
            if(coverage <= 6) {
                out << " color=red";
            }
            out << "];\n";
        }
        out << "}\n";
    }
}
#endif



// Given two consecutive primary edges in the marker graph,
// find the secondary edges in between.
// This works with mode3b assembly and the complete marker graph.
void Assembler::fillMode3bAssemblyPathStep(const vector<string>& request, ostream& html)
{
    // Check that our assumptions are satisfied.
    if(assemblerInfo->assemblyMode != 3) {
        throw runtime_error("This is only available for assembly mode 3.");
    }
    SHASTA_ASSERT(getReads().representation == 0);      // No RLE.
    SHASTA_ASSERT((assemblerInfo->k % 2) == 0);         // Marker length is even.

    // Get the parameters for the request
    uint64_t edgeIdA = invalid<uint64_t>;
    getParameterValue(request, "edgeIdA", edgeIdA);

    uint64_t edgeIdB = invalid<uint64_t>;
    getParameterValue(request, "edgeIdB", edgeIdB);

    // Write the form.
    html <<
        "<form>"
        "<table>"
        "<tr><td class=centered>Edge A<td class=centered>"
        "<input type=text required name=edgeIdA size=8 style='text-align:center' " <<
        ((edgeIdA == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeIdA) + "'")) << ">"
        "<tr><td class=centered>Edge B<td class=centered>"
        "<input type=text required name=edgeIdB size=8 style='text-align:center' " <<
        ((edgeIdB == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeIdB) + "'")) << ">"
        "</table>"
        "<br><input type=submit value='Do it'>"
        "</form>";

    // If the edge id are missing, do nothing.
    if(edgeIdA == invalid<uint64_t> or edgeIdB == invalid<uint64_t>) {
        return;
    }

    // Sanity check that the two edges are distinct.
    if(edgeIdA == edgeIdB) {
        html << "<p>Specify two distinct edges.";
        return;
    }

    // This analysis can only be done if both edges have no duplicate OrientedReadIds
    // in their MarkerIntervals.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA)) {
        html << "<p>Marker graph edge " << edgeIdA << " has duplicate oriented reads.";
        return;
    }
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB)) {
        html << "<p>Marker graph edge " << edgeIdB << " has duplicate oriented reads.";
        return;
    }

    // Check that there are common reds..
    MarkerGraphEdgePairInfo info;
    SHASTA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
    if(info.common == 0) {
        html << "<p>The two edges have no common oriented reads.";
        return;
    }

    // Fill this assembly step.
    mode3b::PathFiller filler(*this, edgeIdA, edgeIdB, html);
}

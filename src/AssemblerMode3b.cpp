// Shasta
#include "Assembler.hpp"
#include "LocalMarkerGraph1.hpp"
#include "mode3b-PathFiller1.hpp"
#include "mode3b-PathFinder.hpp"
#include "mode3b-AssemblyPath.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>



void Assembler::findCompleteMarkerGraphPath(
    MarkerGraphEdgeId startEdgeId,  // The path starts here.
    uint64_t direction              // 0=forward, 1=backward, 2=bidirectional
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
    mode3b::AssemblyPath assemblyPath(
        *this,
        startEdgeId,
        direction);
    ofstream fasta("AssemblyPath.fasta");
    assemblyPath.writeFasta(fasta);
    ofstream csv("AssemblyPath.csv");
    assemblyPath.writeCsv(csv);
}



void Assembler::findCompleteMarkerGraphPaths(uint64_t threadCount) const
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    mode3b::PathFinder pathFinder(*this, threadCount);
}



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



    // Get the parameters for the request.
    uint64_t edgeIdA = invalid<uint64_t>;
    getParameterValue(request, "edgeIdA", edgeIdA);

    uint64_t edgeIdB = invalid<uint64_t>;
    getParameterValue(request, "edgeIdB", edgeIdB);

    string showGraphString;
    const bool showGraph = getParameterValue(request, "showGraph", showGraphString);

    string showVerticesString;
    const bool showVertices = getParameterValue(request, "showVertices", showVerticesString);

    string showVertexLabelsString;
    const bool showVertexLabels = getParameterValue(request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    const bool showEdgeLabels = getParameterValue(request, "showEdgeLabels", showEdgeLabelsString);



    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr><th class=left>Edge A<td class=centered>"
        "<input type=text required name=edgeIdA size=8 style='text-align:center' " <<
        ((edgeIdA == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeIdA) + "'")) << ">"

        "<tr><th class=left>Edge B<td class=centered>"
        "<input type=text required name=edgeIdB size=8 style='text-align:center' " <<
        ((edgeIdB == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeIdB) + "'")) << ">"

        "<tr>"
        "<th class=left>Display the graph"
        "<td class=centered><input type=checkbox name=showGraph" <<
        (showGraph ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the vertices"
        "<td class=centered><input type=checkbox name=showVertices" <<
        (showVertices ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display vertex labels"
        "<td class=centered><input type=checkbox name=showVertexLabels" <<
        (showVertexLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display edge labels"
        "<td class=centered><input type=checkbox name=showEdgeLabels" <<
        (showEdgeLabels ? " checked" : "") << ">"

        "</table>"
        "<br><input type=submit value='Do it'>"
        "</form>";



    // If the edge ids are missing, do nothing.
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

    // Check that there are common reads.
    MarkerGraphEdgePairInfo info;
    SHASTA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
    if(info.common == 0) {
        html << "<p>The two edges have no common oriented reads.";
        return;
    }

    // Fill this assembly step.
    mode3b::PathFiller1 filler(*this, edgeIdA, edgeIdB, html,
        showGraph,
        showVertices,
        showVertexLabels,
        showEdgeLabels);
}

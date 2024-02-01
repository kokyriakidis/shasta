// Shasta
#include "Assembler.hpp"
#include "LocalMarkerGraph1.hpp"
#include "mode3b-PathFiller1.hpp"
#include "mode3b-PathFiller2.hpp"
#include "mode3b-PathFiller3.hpp"
#include "mode3b-PathFinder.hpp"
#include "mode3b-PathGraph.hpp"
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
    assemblyPath.writeFasta(fasta, "Path");
    ofstream csv("AssemblyPath.csv");
    assemblyPath.writeCsv(csv, "Path");
}



void Assembler::findCompleteMarkerGraphPaths(uint64_t threadCount) const
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    mode3b::PathFinder pathFinder(*this, threadCount);
}



void Assembler::flagPrimaryMarkerGraphEdges(
    uint64_t minPrimaryEdgeCoverage,
    uint64_t maxPrimaryEdgeCoverage,
    uint64_t threadCount)
{
    // Check that we have what we need.
    SHASTA_ASSERT(markers.isOpen());
    checkMarkerGraphVerticesAreAvailable();
    SHASTA_ASSERT(markerGraph.edges.isOpenWithWriteAccess);

    markerGraph.flagPrimaryEdges(
        minPrimaryEdgeCoverage,
        maxPrimaryEdgeCoverage,
        markers,
        threadCount);
}



void Assembler::createMarkerGraphPrimaryJourneys(uint64_t threadCount)
{
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();


    markerGraph.primaryJourneys.createNew(
        largeDataName("MarkerGraphPrimaryJourneys"), largeDataPageSize);
    markerGraph.createPrimaryJourneys(markers.size(), threadCount);
}



void Assembler::accessMarkerGraphPrimaryJourneys()
{
    markerGraph.primaryJourneys.accessExistingReadOnly(
        largeDataName("MarkerGraphPrimaryJourneys"));
}



void Assembler::writeMarkerGraphPrimaryJourneys()
{
    SHASTA_ASSERT(markerGraph.primaryJourneys.isOpen());
    markerGraph.writePrimaryJourneys();
}



void Assembler::findMode3bPaths(
    uint64_t threadCount0,  // High level parallelization
    uint64_t threadCount1   // Low level parallelization
    ) const
{
    mode3b::GlobalPathGraph::assemble(*this, threadCount0, threadCount1);
}



void Assembler::loadAndAssembleCompressedPathGraph(
    const string& fileName,
    uint64_t threadCount0,  // High level parallelization
    uint64_t threadCount1   // Low level parallelization
    ) const
{
    mode3b::GlobalPathGraph::loadAndAssemble(*this, fileName, threadCount0, threadCount1);
}



// Given two consecutive primary edges in the marker graph,
// find the secondary edges in between.
// This works with mode3b assembly and the complete marker graph.
void Assembler::fillMode3bAssemblyPathStep(const vector<string>& request, ostream& html)
{
    fillMode3bAssemblyPathStep3(request, html);
}




// Given two consecutive primary edges in the marker graph,
// find the secondary edges in between.
// This works with mode3b assembly and the complete marker graph.
// Version that uses class mode3b::PathFiller1.
void Assembler::fillMode3bAssemblyPathStep1(const vector<string>& request, ostream& html)
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

    string showDebugInformationString;
    const bool showDebugInformation = getParameterValue(request, "showDebugInformation", showDebugInformationString);



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

        "<tr>"
        "<th class=left>Display debug information"
        "<td class=centered><input type=checkbox name=showDebugInformation" <<
        (showDebugInformation ? " checked" : "") << ">"

        "</table>"
        "<br><input type=submit value='Do it'>"
        "</form>";



    // If the edge ids are missing, do nothing.
    if(edgeIdA == invalid<uint64_t> or edgeIdB == invalid<uint64_t>) {
        return;
    }

    // Sanity checks on the edge ids.
    if(edgeIdA >= markerGraph.edges.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdA) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edges.size()));
    }
    if(edgeIdB >= markerGraph.edges.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdB) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edges.size()));
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
        showEdgeLabels,
        showDebugInformation);
}



// Given two consecutive primary edges in the marker graph,
// find the secondary edges in between.
// This works with mode3b assembly and the complete marker graph.
// Version that uses class mode3b::PathFiller2.
void Assembler::fillMode3bAssemblyPathStep2(const vector<string>& request, ostream& html)
{
    // Check that our assumptions are satisfied.
    if(assemblerInfo->assemblyMode != 3) {
        throw runtime_error("This is only available for assembly mode 3.");
    }
    SHASTA_ASSERT(getReads().representation == 0);      // No RLE.
    SHASTA_ASSERT((assemblerInfo->k % 2) == 0);         // Marker length is even.

    mode3b::PathFiller2DisplayOptions options(html);

    // Get the parameters for the request.
    uint64_t edgeIdA = invalid<uint64_t>;
    getParameterValue(request, "edgeIdA", edgeIdA);

    uint64_t edgeIdB = invalid<uint64_t>;
    getParameterValue(request, "edgeIdB", edgeIdB);

    string showGraphString;
    options.showGraph = getParameterValue(request, "showGraph", showGraphString);

    string showVerticesString;
    options.showVertices = getParameterValue(request, "showVertices", showVerticesString);

    string showVertexLabelsString;
    options.showVertexLabels = getParameterValue(request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    options.showEdgeLabels = getParameterValue(request, "showEdgeLabels", showEdgeLabelsString);

    string showAssemblyDetailsString;
    options.showAssemblyDetails = getParameterValue(request, "showAssemblyDetails", showAssemblyDetailsString);

    string showDebugInformationString;
    options.showDebugInformation = getParameterValue(request, "showDebugInformation", showDebugInformationString);



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
        (options.showGraph ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the vertices"
        "<td class=centered><input type=checkbox name=showVertices" <<
        (options.showVertices ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display vertex labels"
        "<td class=centered><input type=checkbox name=showVertexLabels" <<
        (options.showVertexLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display edge labels"
        "<td class=centered><input type=checkbox name=showEdgeLabels" <<
        (options.showEdgeLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display assembly details"
        "<td class=centered><input type=checkbox name=showAssemblyDetails" <<
        (options.showAssemblyDetails ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display debug information"
        "<td class=centered><input type=checkbox name=showDebugInformation" <<
        (options.showDebugInformation ? " checked" : "") << ">"

        "</table>"
        "<br><input type=submit value='Do it'>"
        "</form>";



    // If the edge ids are missing, do nothing.
    if(edgeIdA == invalid<uint64_t> or edgeIdB == invalid<uint64_t>) {
        return;
    }

    // Sanity checks on the edge ids.
    if(edgeIdA >= markerGraph.edges.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdA) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edges.size()));
    }
    if(edgeIdB >= markerGraph.edges.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdB) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edges.size()));
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
    mode3b::PathFiller2 filler(*this, edgeIdA, edgeIdB, options);
}






// Given two consecutive primary edges in the marker graph,
// find the secondary edges in between.
// This works with mode3b assembly and the complete marker graph.
// Version that uses class mode3b::PathFiller3.
void Assembler::fillMode3bAssemblyPathStep3(const vector<string>& request, ostream& html)
{
    // Check that our assumptions are satisfied.
    if(assemblerInfo->assemblyMode != 3) {
        throw runtime_error("This is only available for assembly mode 3.");
    }
    SHASTA_ASSERT(getReads().representation == 0);      // No RLE.
    SHASTA_ASSERT((assemblerInfo->k % 2) == 0);         // Marker length is even.

    mode3b::PathFiller3DisplayOptions options(html);

    // Get the parameters for the request.
    uint64_t edgeIdA = invalid<uint64_t>;
    getParameterValue(request, "edgeIdA", edgeIdA);

    uint64_t edgeIdB = invalid<uint64_t>;
    getParameterValue(request, "edgeIdB", edgeIdB);

    uint64_t minVertexCoverage = 0;
    getParameterValue(request, "minVertexCoverage", minVertexCoverage);

    string showOrientedReadsString;
    options.showOrientedReads = getParameterValue(request, "showOrientedReads", showOrientedReadsString);

    string showMarkersString;
    options.showMarkers = getParameterValue(request, "showMarkers", showMarkersString);

    string showGraphString;
    options.showGraph = getParameterValue(request, "showGraph", showGraphString);

    string showVerticesString;
    options.showVertices = getParameterValue(request, "showVertices", showVerticesString);

    string showVertexLabelsString;
    options.showVertexLabels = getParameterValue(request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    options.showEdgeLabels = getParameterValue(request, "showEdgeLabels", showEdgeLabelsString);

    string showAssemblyDetailsString;
    options.showAssemblyDetails = getParameterValue(request, "showAssemblyDetails", showAssemblyDetailsString);

    string showDebugInformationString;
    options.showDebugInformation = getParameterValue(request, "showDebugInformation", showDebugInformationString);



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

        "<tr><th class=left>Minimum vertex coverage<br>(0 = automatic)<td class=centered>"
        "<input type=text required name=minVertexCoverage size=8 style='text-align:center' "
        "value='" << minVertexCoverage << "'>"

        "<tr>"
        "<th class=left>Display the oriented reads"
        "<td class=centered><input type=checkbox name=showOrientedReads" <<
        (options.showOrientedReads ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the markers"
        "<td class=centered><input type=checkbox name=showMarkers" <<
        (options.showMarkers ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the graph"
        "<td class=centered><input type=checkbox name=showGraph" <<
        (options.showGraph ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the vertices"
        "<td class=centered><input type=checkbox name=showVertices" <<
        (options.showVertices ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display vertex labels"
        "<td class=centered><input type=checkbox name=showVertexLabels" <<
        (options.showVertexLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display edge labels"
        "<td class=centered><input type=checkbox name=showEdgeLabels" <<
        (options.showEdgeLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display assembly details"
        "<td class=centered><input type=checkbox name=showAssemblyDetails" <<
        (options.showAssemblyDetails ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display debug information"
        "<td class=centered><input type=checkbox name=showDebugInformation" <<
        (options.showDebugInformation ? " checked" : "") << ">"

        "</table>"
        "<br><input type=submit value='Do it'>"
        "</form>";



    // If the edge ids are missing, do nothing.
    if(edgeIdA == invalid<uint64_t> or edgeIdB == invalid<uint64_t>) {
        return;
    }

    // Sanity checks on the edge ids.
    if(edgeIdA >= markerGraph.edges.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdA) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edges.size()));
    }
    if(edgeIdB >= markerGraph.edges.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdB) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edges.size()));
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
    mode3b::PathFiller3 filler(*this, edgeIdA, edgeIdB, minVertexCoverage, options);
}

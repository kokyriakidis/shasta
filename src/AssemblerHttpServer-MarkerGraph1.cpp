// Shasta.
#include "Assembler.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "LocalMarkerGraph1.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"



void Assembler::exploreMarkerGraph1(
    const vector<string>& request,
    ostream& html)
{
    if(assemblerInfo->assemblyMode != 3) {
        throw runtime_error("This is only available for assembly mode 3.");
    }

    // This makes the following assumptions.
    SHASTA_ASSERT(getReads().representation == 0);  // No RLE.
    SHASTA_ASSERT((assemblerInfo->k % 2) == 0);      // Marker length is even.


    // Get the request parameters.
    uint64_t vertexId = invalid<uint64_t>;
    getParameterValue(request, "vertexId", vertexId);

    uint64_t maxDistance = 2;
    getParameterValue( request, "maxDistance", maxDistance);

    uint64_t minVertexCoverage = 0;
    getParameterValue(request, "minVertexCoverage", minVertexCoverage);

    uint64_t minEdgeCoverage = 0;
    getParameterValue(request, "minEdgeCoverage", minEdgeCoverage);

    uint64_t maxPruneCoverage = 0;
    getParameterValue(request, "maxPruneCoverage", maxPruneCoverage);

    uint64_t maxLongChainCoverage = 0;
    getParameterValue(request, "maxLongChainCoverage", maxLongChainCoverage);

    uint64_t minLongChainLength = 100;
    getParameterValue(request, "minLongChainLength", minLongChainLength);

    uint64_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    double thicknessScaling = 1.;
    getParameterValue(request, "thicknessScaling", thicknessScaling);

    uint64_t layoutQuality = 2;
    getParameterValue(request, "layoutQuality", layoutQuality);

    double edgeResolution = 1.;
    getParameterValue(request, "edgeResolution", edgeResolution);

    uint64_t redCoverage = 1;
    getParameterValue(request, "redCoverage", redCoverage);

    uint64_t greenCoverage = 5;
    getParameterValue(request, "greenCoverage", greenCoverage);

    string coloring;
    getParameterValue(request, "coloring", coloring);

    uint64_t readFollowingStartEdgeId = 0;
    getParameterValue(request, "readFollowingStartEdgeId", readFollowingStartEdgeId);

    int64_t firstMarkerOffset = 0;
    getParameterValue(request, "firstMarkerOffset", firstMarkerOffset);

    int64_t lastMarkerOffset = 0;
    getParameterValue(request, "lastMarkerOffset", lastMarkerOffset);

    string showLabelsString;
    const bool showLabels = getParameterValue(request, "showLabels", showLabelsString);

    double timeout = 30;
    getParameterValue(request, "timeout", timeout);

    string outputType = "svg";
    getParameterValue(request, "outputType", outputType);


    // Write the form.
    html <<
        "<form>"

        "<h2>Local marker graph</h2>"
        "<table>"

        "<tr>"
        "<td>Start vertex id"
        "<td class=centered><input type=text required name=vertexId size=8 style='text-align:center'"
        << ((vertexId == invalid<uint64_t>) ? "" : ("value='" + to_string(vertexId) + "'")) <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td class=centered><input type=text required name=maxDistance size=8 style='text-align:center'"
        "value='"  << maxDistance << "'>"

        "<tr>"
        "<td>Minimum vertex coverage"
        "<td class=centered><input type=text required name=minVertexCoverage size=8 style='text-align:center'"
        "value='" << minVertexCoverage << "'>"

        "<tr>"
        "<td>Minimum edge coverage"
        "<td class=centered><input type=text required name=minEdgeCoverage size=8 style='text-align:center'"
        "value='" << minEdgeCoverage << "'>"

        "<tr>"
        "<td>Prune leaves with coverage up to"
        "<td class=centered><input type=text required name=maxPruneCoverage size=8 style='text-align:center'"
        "value='" << maxPruneCoverage << "'>"

        "<tr>"
        "<td>Prune long linear sections<br>with low coverage"
        "<td>"
        "<input type=text required name=maxLongChainCoverage size=8 style='text-align:center'"
        "value='" << maxLongChainCoverage << "'> Maximum coverage"
        "<br><input type=text required name=minLongChainLength size=8 style='text-align:center'"
        "value='" << minLongChainLength << "'> Minimum length (markers)"

        "<tr>"
        "<td>Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'"
        " value='" << sizePixels << "'>"

        "<tr>"
        "<td>Thickness scaling factor"
        "<td class=centered><input type=text required name=thicknessScaling size=8 style='text-align:center'"
        " value='" << thicknessScaling << "'>"

        "<tr>"
        "<td>Layout quality"
        "<td class=centered>"
        "<select name=layoutQuality style='text-align:center'>"
        "<option value=0" << (layoutQuality==0 ? " selected" : "") <<
        ">Best speed</option>"
        "<option value=1" << (layoutQuality==1 ? " selected" : "") <<
        ">Intermediate quality and speed</option>"
        "<option value=2" << (layoutQuality==2 ? " selected" : "") <<
        ">Best quality</option>"
        "</select>"

        "<tr>"
        "<td>Edge resolution ";
    writeInformationIcon(html, "Affects edge smoothness and speed of layout computation.");

    html <<
        "<td class=centered><input type=text required name=edgeResolution size=8 style='text-align:center'"
        " value='" << edgeResolution << "'>"

        "<tr>"
        "<td>Coloring"
        "<td>"
        "<select name=coloring style='text-align:center'>"
        "<option value=random" << (coloring == "random" ? " selected" : "") <<
        ">Random</option>"
        "<option value=byCoverage" << (coloring == "byCoverage" ? " selected" : "") <<
        ">By coverage</option>"
        "<option value=readFollowing" << (coloring == "readFollowing" ? " selected" : "") <<
        ">Read following</option>"
        "</select>"
        "<br><input type=text required name=redCoverage size=8 style='text-align:center'"
        " value='" << redCoverage << "'> Red coverage"
        "<br><input type=text required name=greenCoverage size=8 style='text-align:center'"
        " value='" << greenCoverage << "'> Green coverage"
        "<hr><span style='text-align:center'>Read following</span>"
        "<br><input type=text required name=readFollowingStartEdgeId size=8 style='text-align:center'"
        " value='" << readFollowingStartEdgeId << "'> Start edge for read following"
        "<br><input type=text required name=firstMarkerOffset size=8 style='text-align:center'"
        " value='" << firstMarkerOffset << "'> First marker offset"
        "<br><input type=text required name=lastMarkerOffset size=8 style='text-align:center'"
        " value='" << lastMarkerOffset << "'> Last marker offset"

        "<tr>"
        "<td>Show labels"
        "<td class=centered><input type=checkbox name=showLabels" <<
             (showLabels ? " checked" : "") <<
             ">"

        "<tr>"
        "<td>Timeout in seconds"
        "<td class=centered><input type=text required name=timeout size=8 style='text-align:center'"
        " value='" << timeout << "'>"

        "<tr>"
        "<td>Output"
        "<td>"
        "<input type=radio name=outputType value='noOutput'" <<
        (outputType == "noOutput" ? " checked=on" : "") <<
        ">Show the number of vertices and edges"
        "<br><input type=radio name=outputType value='createGfa'" <<
        (outputType == "createGfa" ? " checked=on" : "") <<
        ">Create a GFA file"
        "<br><input type=radio name=outputType value='createAndOpenGfa'" <<
       (outputType == "createAndOpenGfa" ? " checked=on" : "") <<
       ">Create a GFA file and open it in Bandage";

    html <<
       "<br><input type=radio name=outputType value='fastCanvas'" <<
       (outputType == "fastCanvas" ? " checked=on" : "") <<
       ">Display vertices only, not interactive ";
    writeInformationIcon(html, "The fastest choice. "
        "Fast display with one pixel per vertex and no edges, done using canvas. "
        "Best for large subgraphs.");

    html <<
       "<br><input type=radio name=outputType value='fastSvg'" <<
       (outputType == "fastSvg" ? " checked=on" : "") <<
       ">Display vertices only, interactive ";
    writeInformationIcon(html, "Fast display with one pixel per vertex and no edges, done using svg.");

    html <<
       "<br><input type=radio name=outputType value='svg'" <<
       (outputType == "svg" ? " checked=on" : "") <<
       ">Display vertices and edges, interactive ";


    html <<
        "</table>"

        "<br><input type=submit value='Do it'>"
        "</form>";


    // If the vertex id was not specified, stop here.
    if(vertexId == invalid<uint64_t>) {
        return;
    }

    // If the vertex id is invalid, stop here.
    if(vertexId > markerGraph.vertexCount()) {
        html << "<p>Invalid vertex id " << vertexId;
        html << ". Must be between 0 and " << markerGraph.vertexCount()-1 << " inclusive.";
        return;
    }



    // Create the local marker graph.
    LocalMarkerGraph1 graph(
        markers,
        markerGraph,
        vertexId,
        maxDistance,
        minVertexCoverage,
        minEdgeCoverage);

    // Do the requested graph cleanup.
    if(maxPruneCoverage > 0) {
        graph.pruneLowCoverageLeaves(maxPruneCoverage);
    }
    if(maxLongChainCoverage > 0) {
        graph.removeLongLowCoverageChains(maxLongChainCoverage, minLongChainLength);
    }

    html << "<p>The local marker graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges.";


    if(outputType == "noOutput") {
        return;
    }

    if(outputType == "fastCanvas") {
        graph.writeHtml0(html, sizePixels, layoutQuality, timeout, false);
    }

    else if(outputType == "fastSvg") {
        graph.writeHtml0(html, sizePixels, layoutQuality, timeout, true);
    }

    else if(outputType == "svg") {
        graph.writeHtml1(html, sizePixels, thicknessScaling, layoutQuality, edgeResolution,
            coloring, redCoverage, greenCoverage,
            readFollowingStartEdgeId, firstMarkerOffset, lastMarkerOffset,
            showLabels,
            timeout);
    }

    else {

        // Create a gfa file to represent the local marker graph.
        const string gfaFileName = tmpDirectory() + to_string(boost::uuids::random_generator()()) + ".gfa";
        graph.writeGfa(gfaFileName);
        html << "<p>The local marker graph is in "
            "<span id='SpanToBeCopied' style='color:Blue'>" << gfaFileName << "</span>"
            ". Remove it when done with it."
            "<br><button onClick='copySpanToClipboard()'>Copy GFA file name to clipboard</button>";
        html << R"###(
        <script>
        function copySpanToClipboard()
        {
              
            // Remove any previous selection.
            var selection = window.getSelection();
            selection.removeAllRanges();
            
            // Select the span.
            var element = document.getElementById("SpanToBeCopied");
            var range = document.createRange();
            range.selectNodeContents(element);
            selection.addRange(range);
            
            // Copy it to the clipboard.
            document.execCommand("copy");
    
            // Unselect it.
            selection.removeAllRanges();
            
    
        }
        </script>
        )###";


        // If requested, open it in Bandage.
        // This is done on the server side, of course. This can have unexpected
        // consequences if running remotely.
        // Also, because of this the connection with the http client is not closed
        // until Bandage terminates, so the browser thinks ore data are coming.
        if(outputType == "createAndOpenGfa") {
            ::system(("Bandage load " + gfaFileName + "&").c_str());
        }
    }
}


void Assembler::exploreMarkerGraphEdgePair(
    const vector<string>& request,
    ostream& html)
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
        html << "Specify two distinct edges.";
        return;
    }

    // This analysis can only be done if both edges have no duplicate OrientedReadIds
    // in their MarkerIntervals.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA)) {
        html << "Marker graph edge " << edgeIdA << " has duplicate oriented reads.";
        return;
    }
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB)) {
        html << "Marker graph edge " << edgeIdB << " has duplicate oriented reads.";
        return;
    }

    // Write a header.
    html << "<h1>Read composition analysis for marker graph edges " << edgeIdA <<
        " and " << edgeIdB << "</h1>";

    // Analyze read composition.
    MarkerGraphEdgePairInfo info;
    SHASTA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
    writeHtmlMarkerGraphEdgePairInfo(html, edgeIdA, edgeIdB, info);

    if(info.common == 0) {
        const uint64_t estimatedOffset = estimateBaseOffsetUnsafe(edgeIdA, edgeIdB);
        if(estimatedOffset != invalid<uint64_t>) {
            html << "<p>Estimated offset is " << estimatedOffset << " bases.";
        }
    }
}



void Assembler::writeHtmlMarkerGraphEdgePairInfo(
    ostream& html,
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    const MarkerGraphEdgePairInfo& info
    ) const
{
    // Begin the summary table.
    html <<
        "<table>"
        "<tr><th><th>On<br>edge A<th>On<br>edge B";

    // Total.
    html <<
        "<tr><th class=left>Total ";
    writeInformationIcon(html, "The total number of oriented reads on each of the two edges.");
    html << "<td class=centered>" << info.totalA << "<td class=centered>" << info.totalB;

    // Common.
    html << "<tr><th class=left>Common ";
    writeInformationIcon(html, "The number of common oriented reads between the two edges.");
    html <<
        "<td class=centered colspan = 2>" << info.common;

    // Only.
    html <<
        "<tr><th class=left>Only ";
    writeInformationIcon(html, "The number of oriented reads that appear in one edge but not the other.");
    html <<
        "<td class=centered>" << info.onlyA << "<td class=centered>" << info.onlyB;

    // The rest can only be written if there are common reads.
    if(info.common > 0) {

        // Only, short.
        html <<
            "<tr><th class=left>Only, short ";
        writeInformationIcon(html, "The number of oriented reads that appear in one edge only "
            " and are too short to appear on the other edge, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyAShort << "<td class=centered>" << info.onlyBShort;

        // Only, missing.
        html <<
            "<tr><th class=left>Only, missing ";
        writeInformationIcon(html, "The number of oriented reads that appear in one edge only "
            " and are not too short to appear on the other edge, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyA - info.onlyAShort << "<td class=centered>" << info.onlyB - info.onlyBShort;
    }

    // End the summary table.
    html << "</table>";

    // Only write out the rest if there are common reads.
    if(info.common == 0) {
        return;
    }

    // Write the table with Jaccard similarities and estimated offsets.
    using std::fixed;
    using std::setprecision;
    html <<
        "<br><table>"
        "<tr><th class=left>Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.jaccard() <<
        "<tr><th class=left>Corrected Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.correctedJaccard() <<
        "<tr><th class=left>Estimated offset in markers<td class=centered>" << info.offsetInMarkers <<
        "<tr><th class=left>Estimated offset in bases<td class=centered>" << info.offsetInBases <<
        "</table>";


    // Write the details table.
    html <<
        "<br>In the following table, positions in red are hypothetical, based on the above "
        "estimated base offset."
        "<p><table>";

    // Header row.
    html <<
        "<tr>"
        "<th class=centered rowspan=2>Oriented<br>read id"
        "<th class=centered colspan=2>Length"
        "<th colspan=4>Edge A"
        "<th colspan=4>Edge B"
        "<th rowspan=2>Ordinal offset"
        "<th rowspan=2>Base offset"
        "<th rowspan=2>Classification"
        "<tr>"
        "<th>Markers"
        "<th>Bases"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Position0"
        "<th>Position1"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Position0"
        "<th>Position1";

    // Prepare for the joint loop over OrientedReadIds of the two edges.
    const auto markerIntervalsA = markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = markerGraph.edgeMarkerIntervals[edgeIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();

    // Joint loop over the MarkerIntervals of the two edges.
    auto itA = beginA;
    auto itB = beginB;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in edge A.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(getReads().getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinals[0];
            const uint32_t ordinalA1 = itA->ordinals[1];
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Find the hypothetical positions of edge B, assuming the estimated base offset.
            const int64_t positionB0 = positionA0 + info.offsetInBases;
            const int64_t positionB1 = positionA1 + info.offsetInBases;
            const bool isShort = positionB0<0 or positionB1 >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << ordinalA0 <<
                "<td class=centered>" << ordinalA1 <<
                "<td class=centered>" << positionA0 <<
                "<td class=centered>" << positionA1 <<
                "<td><td>"
                "<td class=centered style='color:Red'>" << positionB0 <<
                "<td class=centered style='color:Red'>" << positionB1 << "<td><td>"
                "<td class=centered>OnlyA, " << (isShort ? "short" : "missing");

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in edge B.
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(getReads().getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinals[0];
            const uint32_t ordinalB1 = itB->ordinals[1];
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA0 = positionB0 - info.offsetInBases;
            const int64_t positionA1 = positionB1 - info.offsetInBases;
            const bool isShort = positionA0<0 or positionA1 >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td><td>"
                "<td class=centered style='color:Red'>" << positionA0 <<
                "<td class=centered style='color:Red'>" << positionA1 <<
                "<td class=centered>" << ordinalB0 <<
                "<td class=centered>" << ordinalB1 <<
                "<td class=centered>" << positionB0 <<
                "<td class=centered>" << positionB1 << "<td><td>"
                "<td class=centered>OnlyB, " << (isShort ? "short" : "missing");

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both edges.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(getReads().getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinals[0];
            const uint32_t ordinalA1 = itA->ordinals[1];
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinals[0];
            const uint32_t ordinalB1 = itB->ordinals[1];
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Compute estimated offsets.
            const int64_t ordinalOffset = uint64_t(ordinalB1) - uint64_t(ordinalA0);
            const int64_t baseOffset = positionB1 - positionA0;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << ordinalA0 <<
                "<td class=centered>" << ordinalA1 <<
                "<td class=centered>" << positionA0 <<
                "<td class=centered>" << positionA1 <<
                "<td class=centered>" << ordinalB0 <<
                "<td class=centered>" << ordinalB1 <<
                "<td class=centered>" << positionB0 <<
                "<td class=centered>" << positionB1 <<
                "<td class=centered>" << ordinalOffset <<
                "<td class=centered>" << baseOffset <<
                "<td class=centered>Common";

            ++itA;
            ++itB;
        }
    }

    // Finish the details table.
    html << "</table>";


}


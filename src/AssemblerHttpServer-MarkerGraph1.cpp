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

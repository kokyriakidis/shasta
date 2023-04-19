// Shasta.
#include "Assembler.hpp"
#include "invalid.hpp"
#include "LocalMarkerGraph1.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>



void Assembler::exploreMarkerGraph1(
    const vector<string>& request,
    ostream& html)
{
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

    uint64_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    string outputType = "createGfa";
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
        "<td>Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'"
        " value='" << sizePixels << "'>"

        "<tr>"
        "<td>Output"
        "<td>"
        "<input type=radio name=outputType value='createGfa'" <<
        (outputType == "createGfa" ? " checked=on" : "") <<
        ">Create a GFA file"
         "<br><input type=radio required name=outputType value='createAndOpenGfa'" <<
        (outputType == "createAndOpenGfa" ? " checked=on" : "") <<
        ">Create a GFA file and open it in Bandage"

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
        markerGraph,
        vertexId,
        maxDistance,
        minVertexCoverage,
        minEdgeCoverage);
    html << "<p>The local marker graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges.";

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
    if(outputType == "createAndOpenGfa") {
        ::system(("Bandage load " + gfaFileName + "&").c_str());
    }
}

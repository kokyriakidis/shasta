#include "Assembler.hpp"
#include "invalid.hpp"
using namespace shasta;



void Assembler::exploreMarkerGraph1(
    const vector<string>& request,
    ostream& html)
{
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

        "</table>"

        "<br><input type=submit value='Display'>"
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
}
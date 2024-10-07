// Shasta.
#include "Mode3Assembler.hpp"
#include "HttpServer.hpp"
#include "Marker.hpp"
using namespace shasta;
using namespace mode3;



void Mode3Assembler::exploreAnchor(const vector<string>& request, ostream& html)
{
    // Get the AnchorId.
    string anchorIdString;
    const bool anchorIdStringIsPresent = HttpServer::getParameterValue(request, "anchorIdString", anchorIdString);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for anchor'> "
        "<input type=text name=anchorIdString required";
    if(anchorIdStringIsPresent) {
        html << " value='" << anchorIdString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>";
    html << "</form>";

    // If the anchor id missing or invalid, stop here.
    if(not anchorIdStringIsPresent) {
        return;
    }
    const AnchorId anchorId = anchorIdFromString(anchorIdString);

    if((anchorId == invalid<AnchorId>) or (anchorId >= anchors().size())) {
        html << "<p>Invalid anchor id. Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }


    html << "<h1>Anchor " << anchorIdString << "</h1>";

    const auto markerIntervals = anchors()[anchorId];
    const uint64_t coverage = markerIntervals.size();
    const auto sequence = anchors().anchorSequences[anchorId];

    // Write a summary table.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "<tr><th class=left>Sequence length<td class=centered>" << sequence.size() <<
        "<tr><th class=left>Sequence<td class=centered style='font-family:courier'>";
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));
    html << "</table>";



    // Write the marker intervals of this Anchor.
    html <<
        "<h2>Marker intervals</h2>"
        "<table>"
        "<tr>"
        "<th>Index"
        "<th>Oriented<br>read<br>id"
        "<th>Position<br>in<br>journey"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Position0"
        "<th>Position1";

    // Loop over the marker intervals.
    for(uint64_t i=0; i<coverage; i++) {
        const AnchorMarkerInterval& markerInterval = markerIntervals[i];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;

        const uint32_t ordinal0 = markerInterval.ordinal0;
        const uint32_t ordinal1 = ordinal0 + anchors().ordinalOffset(anchorId);

        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t position0 = orientedReadMarkers[ordinal0].position;
        const uint32_t position1 = orientedReadMarkers[ordinal1].position;

        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << markerInterval.positionInJourney <<
            "<td class=centered>" << ordinal0 <<
            "<td class=centered>" << ordinal1 <<
            "<td class=centered>" << position0 <<
            "<td class=centered>" << position1;
    }
}



void Mode3Assembler::exploreAnchorPair(const vector<string>& /* request */, ostream& html)
{
    html << "Not implemented.";
}



void Mode3Assembler::exploreLocalAssembly(const vector<string>& /* request */, ostream& html)
{
    html << "Not implemented.";
}



void Mode3Assembler::exploreLocalAnchorGraph(const vector<string>& /* request */, ostream& html)
{
    html << "Not implemented.";
}

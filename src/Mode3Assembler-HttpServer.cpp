// Shasta.
#include "Mode3Assembler.hpp"
#include "deduplicate.hpp"
#include "HttpServer.hpp"
#include "Marker.hpp"
#include "mode3-AssemblyGraphPostprocessor.hpp"
#include "mode3-LocalAssembly.hpp"
#include "mode3-LocalAnchorGraph.hpp"
#include "mode3-LocalAssemblyGraph.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

namespace shasta {
    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.
}

// Boost libraries.
#include <boost/tokenizer.hpp>

// Standard library.
#include "fstream.hpp"



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

    vector<AnchorId> parents;
    vector<uint64_t> parentsCoverage;
    anchors().findParents(anchorId, parents, parentsCoverage);

    vector<AnchorId> children;
    vector<uint64_t> childrenCoverage;
    anchors().findChildren(anchorId, children, childrenCoverage);

    // Write a summary table.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "<tr><th class=left>Sequence length<td class=centered>" << sequence.size() <<
        "<tr><th class=left>Sequence<td class=centered style='font-family:courier'>";
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));

    html <<
        "<tr><th class=left>Parent anchors"
        "<td>";
    for(uint64_t i=0; i<parents.size(); i++) {
        const string parentAnchorIdString = anchorIdToString(parents[i]);
        if(i != 0) {
            html << "<br>";
        }
        html <<
            "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(parentAnchorIdString) << "'>" <<
            parentAnchorIdString << "</a> coverage " << parentsCoverage[i];

    }

    html <<
        "<tr><th class=left>Children anchors"
        "<td>";
    for(uint64_t i=0; i<children.size(); i++) {
        const string childAnchorIdString = anchorIdToString(children[i]);
        if(i != 0) {
            html << "<br>";
        }
        html <<
            "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(childAnchorIdString) << "'>" <<
            childAnchorIdString << "</a> coverage " << childrenCoverage[i];

    }


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
            "<td class=centered>" << i;

        // The OrientedReadId is written with an hyperlink that will
        // display the portion of the oriented read around this Anchor.
        const string url =
            "exploreReadSequence?"
            "readId=" + to_string(orientedReadId.getReadId()) +
            "&strand=" + to_string(orientedReadId.getStrand()) +
            "&beginPosition=" + to_string((position0 > 2 * k) ? (position0 - 2 * k) : 0) +
            "&endPosition=" + to_string(position1 + 3 * k - 1);
        html <<
            "<td class=centered>" <<
            "<a href='" << url << "'>" <<
            orientedReadId << "</a>";

       html <<
            "<td class=centered>" << markerInterval.positionInJourney <<
            "<td class=centered>" << ordinal0 <<
            "<td class=centered>" << ordinal1 <<
            "<td class=centered>" << position0 <<
            "<td class=centered>" << position1;
    }
}



void Mode3Assembler::exploreAnchorPair(const vector<string>& request, ostream& html)
{

    // Get the parameters for the request
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);
    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);

    // Write the form.
    html << "<form><table>";

    html <<
        "<tr><th class=left>Anchor A"
        "<td class=centered><input type=text name=anchorIdAString required";
    if(anchorIdAStringIsPresent) {
        html << " value='" << anchorIdAString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "<tr><th class=left>Anchor B"
        "<td class=centered><input type=text name=anchorIdBString required";
    if(anchorIdBStringIsPresent) {
        html << " value='" << anchorIdBString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "</table>"
        "<input type=submit value='Analyze anchor pair'>"
        "</form>";



    // Check the AnchorIds
    if(not (anchorIdAStringIsPresent and anchorIdBStringIsPresent)) {
        return;
    }
    const AnchorId anchorIdA = anchorIdFromString(anchorIdAString);
    const AnchorId anchorIdB = anchorIdFromString(anchorIdBString);

    if((anchorIdA == invalid<AnchorId>) or (anchorIdA >= anchors().size())) {
        html << "<p>Invalid anchor id " << anchorIdAString << ". Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }

    if((anchorIdB == invalid<AnchorId>) or (anchorIdB >= anchors().size())) {
        html << "<p>Invalid anchor id " << anchorIdBString << " .Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }

    if(anchorIdA == anchorIdB) {
        html << "Specify two distinct anchors.";
        return;
    }



    // Write a header.
    html << "<h1>Read composition analysis for anchors " << anchorIdToString(anchorIdA) <<
        " and " << anchorIdToString(anchorIdB) << "</h1>";

    // Analyze this anchor pair and write the result to html.
    AnchorPairInfo info;
    anchors().analyzeAnchorPair(anchorIdA, anchorIdB, info);
    anchors().writeHtml(anchorIdA, anchorIdB, info, html);
}



void Mode3Assembler::exploreJourney(const vector<string>& request, ostream& html)
{
    html <<
        "<h2>Oriented read journey</h2>"
        "<p>The journey of an oriented read is the sequence of anchors it visits.";


    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = HttpServer::getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = HttpServer::getParameterValue(request, "strand", strand);
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = HttpServer::getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = HttpServer::getParameterValue(request, "endPosition", endPosition);

    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads.readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "<tr>"
        "<th class=left>Begin position"
        "<td><input type=text name=beginPosition"
        " title='Leave blank to begin display at beginning of read.'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>End position"
        "<td><input type=text name=endPosition"
        " title='Leave blank to end display at end of read.'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << ">";


    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(not readIdIsPresent) {
        html << "Specify a numeric read id.";
        return;
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= reads.readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }

    // Adjust the position range, if necessary.
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(reads.getReadRawSequenceLength(readId));
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }

    // Access the information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const span<const CompressedMarker> orientedReadMarkers = markers[orientedReadId.getValue()];
    const auto& journeys = anchors().journeys;
    SHASTA_ASSERT(journeys.isOpen());
    SHASTA_ASSERT(journeys.size() == 2 * reads.readCount());
    const span<const AnchorId> journey = journeys[orientedReadId.getValue()];

    // Page title.
    html << "<h2>Journey of oriented read " << orientedReadId << "</h2>"
        "<p>Each anchor begins at the midpoint of the first marker (Marker0) "
        "and ends at the midpoint of the second marker (Marker1).";

    // Begin the main table.
    html <<
        "<table><tr>"
        "<th>Position<br>in journey"
        "<th>Anchor"
        "<th>Anchor<br>coverage"
        "<th>Marker0<br>ordinal"
        "<th>Marker1<br>ordinal"
        "<th>Marker0<br>position"
        "<th>Marker1<br>position"
        "<th>Anchor<br>begin"
        "<th>Anchor<br>end"
        "<th>Anchor<br>length"
        "<th>Anchor<br>sequence";

    // Loop over the anchors in the journey of this oriented read.
    for(uint64_t positionInJourney=0; positionInJourney<journey.size(); positionInJourney++) {
        const AnchorId anchorId = journey[positionInJourney];
        const uint64_t anchorCoverage = anchors()[anchorId].coverage();
        const string anchorIdString = anchorIdToString(anchorId);

        const uint64_t ordinal0 = anchors().getFirstOrdinal(anchorId, orientedReadId);
        const uint64_t ordinal1 = ordinal0 + anchors().ordinalOffset(anchorId);

        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t position0 = orientedReadMarkers[ordinal0].position;
        const uint32_t position1 = orientedReadMarkers[ordinal1].position;

        if(position0 < beginPosition) {
            continue;
        }
        if(position1 >= endPosition) {
            continue;
        }

        const span<const Base> anchorSequence = anchors().anchorSequence(anchorId);
        SHASTA_ASSERT(anchorSequence.size() == position1 - position0);

        html <<
            "<tr>"
            "<td class=centered>" << positionInJourney <<
            "<td class=centered>" <<
            "<a href='exploreAnchor?anchorIdString=" << anchorIdString << "'>" <<
            anchorIdString << "</a>"
            "<td class=centered>" << anchorCoverage <<
            "<td class=centered>" << ordinal0 <<
            "<td class=centered>" << ordinal1 <<
            "<td class=centered>" << position0 <<
            "<td class=centered>" << position1 <<
            "<td class=centered>" << position0 + k/2 <<
            "<td class=centered>" << position1 + k/2 <<
            "<td class=centered>" << anchorSequence.size() <<
            "<td class=centered style='font-family:courier'>";
        copy(anchorSequence.begin(), anchorSequence.end(), ostream_iterator<Base>(html));
    }

    html << "</table>";
}



void Mode3Assembler::exploreLocalAssembly(
    const vector<string>& request,
    ostream& html)
{
    const Mode3AssemblyOptions::LocalAssemblyOptions& localAssemblyOptions = options.localAssemblyOptions;
    LocalAssemblyDisplayOptions options(html);

    // Get the parameters for the request.
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);

    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);

    string useAString;
    const bool useA = HttpServer::getParameterValue(request, "useA", useAString);

    string useBString;
    const bool useB = HttpServer::getParameterValue(request, "useB", useBString);

    uint64_t minVertexCoverage = 0;
    HttpServer::getParameterValue(request, "minVertexCoverage", minVertexCoverage);

    string showOrientedReadsString;
    options.showOrientedReads = HttpServer::getParameterValue(request, "showOrientedReads", showOrientedReadsString);

    string showMarkersString;
    options.showMarkers = HttpServer::getParameterValue(request, "showMarkers", showMarkersString);

    string showGraphString;
    options.showGraph = HttpServer::getParameterValue(request, "showGraph", showGraphString);

    string showVerticesString;
    options.showVertices = HttpServer::getParameterValue(request, "showVertices", showVerticesString);

    string showVertexLabelsString;
    options.showVertexLabels = HttpServer::getParameterValue(request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    options.showEdgeLabels = HttpServer::getParameterValue(request, "showEdgeLabels", showEdgeLabelsString);

    string showAssemblyDetailsString;
    options.showAssemblyDetails = HttpServer::getParameterValue(request, "showAssemblyDetails", showAssemblyDetailsString);

    string showDebugInformationString;
    options.showDebugInformation = HttpServer::getParameterValue(request, "showDebugInformation", showDebugInformationString);



    // Write the form.
    html <<
        "<form>"
        "<table>";

    html <<
        "<tr><th class=left>Anchor A"
        "<td class=centered><input type=text name=anchorIdAString required";
    if(anchorIdAStringIsPresent) {
        html << " value='" << anchorIdAString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "<tr><th class=left>Anchor B"
        "<td class=centered><input type=text name=anchorIdBString required";
    if(anchorIdBStringIsPresent) {
        html << " value='" << anchorIdBString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "<tr>"
        "<th class=left>Use for assembly oriented reads that appear only on anchor A"
        "<td class=centered><input type=checkbox name=useA" <<
        (useA ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Use for assembly oriented reads that appear only on anchor B"
        "<td class=centered><input type=checkbox name=useB" <<
        (useB ? " checked" : "") << ">"

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



    // Check the AnchorIds
    if(not (anchorIdAStringIsPresent and anchorIdBStringIsPresent)) {
        return;
    }
    const AnchorId anchorIdA = anchorIdFromString(anchorIdAString);
    const AnchorId anchorIdB = anchorIdFromString(anchorIdBString);

    if((anchorIdA == invalid<AnchorId>) or (anchorIdA >= anchors().size())) {
        html << "<p>Invalid anchor id " << anchorIdAString << ". Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }

    if((anchorIdB == invalid<AnchorId>) or (anchorIdB >= anchors().size())) {
        html << "<p>Invalid anchor id " << anchorIdBString << " .Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }

    if(anchorIdA == anchorIdB) {
        html << "Specify two distinct anchors.";
        return;
    }

    // Local assembly for this assembly step.
    LocalAssembly localAssembly(
        k, reads, markers, anchors(),
        anchorIdA, anchorIdB, minVertexCoverage,
        options,
        localAssemblyOptions,
        useA, useB);
}



void Mode3Assembler::exploreLocalAnchorGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the options that control graph creation.
    string anchorIdsString;
    HttpServer::getParameterValue(request, "anchorIdsString", anchorIdsString);

    uint64_t distance = 10;
    HttpServer::getParameterValue(request, "distance", distance);

    string filterEdgesByCoverageLossString;
    const bool filterEdgesByCoverageLoss = HttpServer::getParameterValue(request,
        "filterEdgesByCoverageLoss", filterEdgesByCoverageLossString);

    double maxCoverageLoss =  options.primaryGraphOptions.maxLoss;
    HttpServer::getParameterValue(request, "maxCoverageLoss", maxCoverageLoss);

    // Get the options that control graph display.
    const LocalAnchorGraphDisplayOptions displayOptions(request);



    // Start the form.
    html << "<form><table>";

    // Form items for options that control graph creation.
    html <<
        "<tr title='Enter comma separated anchor ids, each a number between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>"
        "<th class=left>Starting anchor ids"
        "<td class=centered><input type=text name=anchorIdsString style='text-align:center' required";
    if(not anchorIdsString.empty()) {
        html << " value='" << anchorIdsString + "'";
    }
    html <<
        " size=8 title='Enter comma separated anchor ids, each between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>";

    html << "<tr>"
        "<th class=left>Distance"
        "<td class=centered>"
        "<input type=text name=distance style='text-align:center' required size=8 value=" <<
        distance << ">";

    html <<
        "<tr><th>Edge filtering"
        "<td><input type=checkbox name=filterEdgesByCoverageLoss" <<
        (filterEdgesByCoverageLoss ? " checked" : "") <<
        ">Filter edges by coverage loss"
        "<br><input type=text name=maxCoverageLoss style='text-align:center' required size=6 value=" <<
        maxCoverageLoss << "> Maximum coverage loss"
        "<hr>";

    // Form items for options that control graph display.
    displayOptions.writeForm(html);

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Create local anchor graph'>"
        "</form>";



    // If the anchor id are missing, stop here.
    if(anchorIdsString.empty()) {
        return;
    }


    // Extract the AnchorIds.
    vector<AnchorId> anchorIds;
    boost::tokenizer< boost::char_separator<char> > tokenizer(anchorIdsString, boost::char_separator<char>(","));

    for(const string& anchorIdString: tokenizer) {
        const AnchorId anchorId = anchorIdFromString(anchorIdString);

        if((anchorId == invalid<AnchorId>) or (anchorId >= anchors().size())) {
            html << "<p>Invalid anchor id " << anchorIdString << ". Must be a number between 0 and " <<
                anchors().size() / 2 - 1 << " followed by + or -.";
            return;
        }

        anchorIds.push_back(anchorId);
    }
    deduplicate(anchorIds);



    // Create the LocalAnchorGraph starting from these AnchorIds and moving
    // away up to the specified distance.
    LocalAnchorGraph graph(
        anchors(),
        anchorIds,
        distance,
        filterEdgesByCoverageLoss,
        maxCoverageLoss
        );

    html << "<h1>Local anchor graph</h1>";
    html << "The local anchor graph has " << num_vertices(graph) <<
         " vertices and " << num_edges(graph) << " edges.";

    // Write it to html.
    graph.writeHtml(html, displayOptions);

}



void Mode3Assembler::exploreAssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the options that control graph creation.
    string assemblyStage = "Final";
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    uint64_t componentId = 0;
    HttpServer::getParameterValue(request, "componentId", componentId);

    string chainStringIds;
    HttpServer::getParameterValue(request, "chainStringIds", chainStringIds);

    uint64_t distance = 10;
    HttpServer::getParameterValue(request, "distance", distance);

    // Get the options that control graph display.
    const LocalAssemblyGraphDisplayOptions displayOptions(request);



    // Start the form.
    html << "<h2>Local assembly graph</h2><form><table>";

    // Form items for options to choose the assembly graph to be used.
    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=8>";

    html <<
        "<tr>"
        "<th class=left>Component"
        "<td class=centered><input type=text name=componentId style='text-align:center' required"
        " value='" << componentId << "' size=8>";

    // Form items for options that control graph creation.
    html <<
        "<tr title='Enter comma separated Chain (Segment) ids, each of the form a-b-c-d-Pn'>"
        "<th class=left>Starting segments"
        "<td class=centered><input type=text name=chainStringIds style='text-align:center' required";
    if(not chainStringIds.empty()) {
        html << " value='" << chainStringIds + "'";
    }
    html << " size=40>";

    html << "<tr>"
        "<th class=left>Distance"
        "<td class=centered>"
        "<input type=text name=distance style='text-align:center' required size=8 value=" <<
        distance << ">";

    // Form items for options that control graph display.
    displayOptions.writeForm(html);

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Show local assembly graph'>"
        "</form>";


    if(assemblyStage.empty()) {
        return;
    }

    if(chainStringIds.empty()) {
        return;
    }

    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage, componentId);


    // Extract the ChainIdentifiers for the starting chains.
    vector<ChainIdentifier> startingChains;
    boost::tokenizer< boost::char_separator<char> > tokenizer(chainStringIds, boost::char_separator<char>(","));
    for(const string& chainStringId: tokenizer) {
        startingChains.push_back(assemblyGraph.getChainIdentifier(chainStringId));
    }


    // Create the LocalAssemblyGraph.
    LocalAssemblyGraph graph(assemblyGraph, startingChains, distance, assemblyStage);

    html <<
        "<h1>Local assembly graph</h1>"
        "<p>The local assembly graph has " << num_vertices(graph) <<
         " vertices and " << num_edges(graph) << " edges.";

    // Write it to html.
    graph.writeHtml(html, displayOptions, assemblyStage);

}



void Mode3Assembler::exploreSegment(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage = "Final";
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    string displayAnchors = "none";
    HttpServer::getParameterValue(request, "displayAnchors", displayAnchors);

    string beginString;
    HttpServer::getParameterValue(request, "begin", beginString);

    string endString;
    HttpServer::getParameterValue(request, "end", endString);

    string firstAnchorsCountString;
    HttpServer::getParameterValue(request, "firstAnchorsCount", firstAnchorsCountString);

    string lastAnchorsCountString;
    HttpServer::getParameterValue(request, "lastAnchorsCount", lastAnchorsCountString);


    // Start the form.
    html << "<h2>Assembly graph segment</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=30>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << " title='Enter a segment name of the form a-b-c-d-Pn' size=30>";



    // Options to control which anchors are shown.
    html <<
        "<tr>"
        "<th class=left>Show anchors"
        "<td class=left>"

        "<input type=radio required name=displayAnchors value='none'" <<
        (displayAnchors == "none" ? " checked=on" : "") << "> None"

        "<br><input type=radio required name=displayAnchors value='all'" <<
        (displayAnchors == "all" ? " checked=on" : "") << "> All"

        "<br><input type=radio required name=displayAnchors value='range'" <<
        (displayAnchors == "range" ? " checked=on" : "") << "> Anchors in position range "
        "<input type=text name=begin size=8 style='text-align:center' value='" << beginString << "'> to "
        "<input type=text name=end size=8 style='text-align:center' value='" << endString << "'>"

        "<br><input type=radio required name=displayAnchors value='first'" <<
        (displayAnchors == "first" ? " checked=on" : "") << "> First "
        "<input type=text name=firstAnchorsCount size=8 style='text-align:center' value='" << firstAnchorsCountString << "'>"
        " anchors"

        "<br><input type=radio required name=displayAnchors value='last'" <<
        (displayAnchors == "last" ? " checked=on" : "") << "> Last "
        "<input type=text name=lastAnchorsCount size=8 style='text-align:center' value='" << lastAnchorsCountString << "'>"
        " anchors"
        ;


    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Get segment information'>"
        "</form>";


    if(segmentName.empty()) {
        return;
    }

    // Parse the segment name.
    uint64_t componentId;
    uint64_t bubbleChainId;
    uint64_t positionInBubbleChain;
    uint64_t indexInBubble;
    uint64_t bubblePloidy;
    AssemblyGraphPostprocessor::parseChainStringId(
        segmentName,
        componentId,
        bubbleChainId,
        positionInBubbleChain,
        indexInBubble,
        bubblePloidy);

    // Get the AssemblyGraph for this component.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage, componentId);

    // Extract the Chain (Segment) we want.
    const Chain& chain = assemblyGraph.getChain(segmentName);

    html << "<h2>Assembly graph segment (chain) " << segmentName <<
        " at assembly stage " << assemblyStage <<
        "</h2>";



    // Write a summary table for this chain.
    html <<
        "<table>"
        "<tr><th class=left>Name<td class=centered>" << segmentName <<
        "<tr><th class=left>Component<td class=centered>" << componentId <<
        "<tr><th class=left>Bubble chain<td class=centered>" << bubbleChainId <<
        "<tr><th class=left>Bubble position in bubble chain<td class=centered>" << positionInBubbleChain <<
        "<tr><th class=left>Index in bubble<td class=centered>" << indexInBubble <<
        "<tr><th class=left>Bubble ploidy<td class=centered>" << bubblePloidy <<
        "<tr><th class=left>Number of anchors<td class=centered>" << chain.size();
    if(assemblyGraph.sequenceWasAssembled) {
        html << "<tr><th class=left>Sequence length<td class=centered>" << chain.sequence.size();

    }
    html << "</table>";



    // Write the anchors.
    if(displayAnchors != "none") {

        // Figure out the anchor position range to use.
        uint64_t begin = invalid<uint64_t>;
        uint64_t end = invalid<uint64_t>;
        if(displayAnchors == "all") {
            begin = 0;
            end = chain.size();
        } else if(displayAnchors == "range") {
            try {
                begin = atoul(beginString);
            } catch(std::exception& e) {
                throw runtime_error("Begin string " + beginString + " is not valid. Must be a number.");
            }
            try {
                end = atoul(endString);
            } catch(std::exception& e) {
                throw runtime_error("End string " + endString + " is notvalid. Must be a number.");
            }
            if(begin > chain.size()) {
                begin = chain.size() - 1;
            }
            if(end > chain.size()) {
                end = chain.size();
            }
            if(end < begin) {
                end = begin + 1;
            }
        } else if(displayAnchors == "first") {
            begin = 0;
            try {
                end = atoul(firstAnchorsCountString);
            } catch(std::exception& e) {
                throw runtime_error("First anchors count string " + firstAnchorsCountString + " is not valid. Must be a number.");
            }
            if(end > chain.size()) {
                end = chain.size();
            }
        } else if(displayAnchors == "last") {
            end = chain.size();
            uint64_t count = invalid<uint64_t>;
            try {
                count = atoul(lastAnchorsCountString);
            } catch(std::exception& e) {
                throw runtime_error("Last anchors count string " + lastAnchorsCountString + " is not valid. Must be a number.");
            }
            if(count > chain.size()) {
                begin = 0;
            } else {
                begin = end - count;
            }
        } else {
            SHASTA_ASSERT(0);
        }
        SHASTA_ASSERT(end > begin);

        html << "<h3>Anchors</h3>";

        // Link to a local anchor graph with these anchors.
        {
            string anchorIds;
            for(uint64_t positionInChain=begin; positionInChain!=end; positionInChain++) {
                const AnchorId anchorId = chain[positionInChain];
                const string anchorIdString = anchorIdToString(anchorId);
                anchorIds += anchorIdString;
                anchorIds += ",";
            }
            anchorIds.pop_back();

            html << "<p><a href='exploreLocalAnchorGraph?anchorIdsString=" << HttpServer::urlEncode(anchorIds) <<
                "&filterEdgesByCoverageLoss=on'>See these anchors in a local anchor graph</a>";
        }

        // Write a table with the requested anchors.
        html <<
            "<p><table>"
            "<tr>"
            "<th>Position<br>in segment"
            "<th>Anchor<br>id"
            "<th>Link<br>to<br>anchor<br>graph"
            "<th>Coverage";
        for(uint64_t positionInChain=begin; positionInChain!=end; positionInChain++) {
            const AnchorId anchorId = chain[positionInChain];
            const auto markerIntervals = anchors()[anchorId];
            const uint64_t coverage = markerIntervals.size();
            const string anchorIdString = anchorIdToString(anchorId);

            html <<
                "<tr><td class=centered>" << positionInChain <<
                "<td class=centered>" <<
                "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "'>" <<
                anchorIdString << "</a>"
                "<td class=centered>" <<
                "<a href='exploreLocalAnchorGraph?anchorIdsString=" << HttpServer::urlEncode(anchorIdString) <<
                "&filterEdgesByCoverageLoss=on'>" <<
                "&#x22B6;</a>"
                "<td class=centered>" << coverage;
        }

        html << "</table>";

    }


}



const AssemblyGraphPostprocessor& Mode3Assembler::getAssemblyGraph(
    const string& assemblyStage,
    uint64_t componentId
    )
{
    auto it = assemblyGraphsMap.find({assemblyStage, componentId});

    if(it == assemblyGraphsMap.end()) {

        // This AssemblyGraph is not among the ones we already loaded. Load it now.
        cout << timestamp << "Loading assembly graph for stage " << assemblyStage <<
            " component " << componentId << endl;
        shared_ptr<const AssemblyGraphPostprocessor> assemblyGraphPointer =
            make_shared<const AssemblyGraphPostprocessor>(
            assemblyStage,
            componentId,
            anchors(),
            options);
        cout << timestamp << "Done loading assembly graph for stage " << assemblyStage <<
            " component " << componentId << endl;
        assemblyGraphsMap.insert({{assemblyStage, componentId}, assemblyGraphPointer});
        SHASTA_ASSERT(assemblyGraphPointer->componentId == componentId);
        return *assemblyGraphPointer;

    } else {

        // This AssemblyGraph is among the ones we already loaded. Return a reference to it.
        const auto assemblyGraphPointer = it->second;
        SHASTA_ASSERT(assemblyGraphPointer->componentId == componentId);
        return *assemblyGraphPointer;
    }
}

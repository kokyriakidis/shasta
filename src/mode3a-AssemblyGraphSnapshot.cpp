#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "Base.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "algorithm.hpp"
#include "fstream.hpp"



// This creates a snapshot of the AssemblyGraph in the current state.
AssemblyGraphSnapshot::AssemblyGraphSnapshot(
    const AssemblyGraph& assemblyGraph,
    const string& name,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    name(name),
    packedMarkerGraph(assemblyGraph.packedMarkerGraph)
{
    // Store the segments.
    createNew(vertexVector, name + "-vertices");
    createNew(vertexJourneyEntries, name + "-vertexJourneyEntries");
    std::map<AssemblyGraph::vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& assemblyGraphVertex = assemblyGraph[v];
        vertexMap.insert(make_pair(v, vertexVector.size()));
        vertexVector.push_back(Vertex(assemblyGraphVertex));
        vertexJourneyEntries.appendVector(assemblyGraphVertex.journeyEntries);
    }

    // Store the edges.
    createNew(edgeVector, name + "-edges");
    std::map<AssemblyGraph::edge_descriptor, uint64_t> edgeMap;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
        const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
        edgeMap.insert(make_pair(e, edgeVector.size()));
        edgeVector.push_back({vertexMap[v0], vertexMap[v1]});
    }



    // Copy the journeys to the snapshot.
    createNew(journeys, name + "-journeys");
    for(uint64_t i=0; i<assemblyGraph.journeys.size(); i++) {
        const vector<AssemblyGraph::vertex_descriptor>& assemblyGraphJourney = assemblyGraph.journeys[i];
        journeys.appendVector();
        for(const AssemblyGraph::vertex_descriptor v: assemblyGraphJourney) {
            uint64_t vertexId = invalid<uint64_t>;
            auto it = vertexMap.find(v);
            if(it != vertexMap.end()) {
                vertexId = it->second;
            }
            journeys.append(vertexId);
        }
    }



    // Compute connectivity.
    createNew(edgesBySource, name + "-edgesBySource");
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        edgesBySource.appendVector();
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            edgesBySource.append(edgeMap[e]);
        }
    }
    createNew(edgesByTarget, name + "-edgesByTarget");
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        edgesByTarget.appendVector();
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            edgesByTarget.append(edgeMap[e]);
        }
    }

    createVertexTable(assemblyGraph.packedMarkerGraph);


    // Store AssemblyPaths.
    createNew(assemblyPaths, name + "-assemblyPaths");
    for(const auto& assemblyPathPointer: assemblyGraph.assemblyPaths) {
        const auto& assemblyPath = *assemblyPathPointer;
        assemblyPaths.appendVector();
        SHASTA_ASSERT(
            assemblyPath.primaryVertices.size() ==
            assemblyPath.secondaryVerticesInfos.size() + 1);
        for(uint64_t i=0; /* Check later */ ; i++) {

            // Add the primary vertex.
            const AssemblyGraph::vertex_descriptor v = assemblyPath.primaryVertices[i];
            const AssemblyGraphVertex& vertex = assemblyGraph[v];
            const uint64_t vertexId = getVertexId(vertex.segmentId, vertex.segmentReplicaIndex);
            assemblyPaths.append({vertexId, true});

            // If this is the last primary vertex, we are done.
            // There are no secondary vertices following it.
            if(i == assemblyPath.primaryVertices.size() - 1) {
                break;
            }

            // Add the secondary vertices.
            const auto& secondaryVerticesInfos = assemblyPath.secondaryVerticesInfos[i];
            for(const AssemblyGraph::vertex_descriptor v: secondaryVerticesInfos.secondaryVertices) {
                const AssemblyGraphVertex& vertex = assemblyGraph[v];
                const uint64_t vertexId = getVertexId(vertex.segmentId, vertex.segmentReplicaIndex);
                assemblyPaths.append({vertexId, false});
            }
        }
    }
}



// The vertex table is a data structure that allows to get a vertexId (index in vertexVector)
// given a segmentId and segmentCopyIndex.
// Indexed by segmentId.
// vertexTable[segmentId][segmentCopyIndex] contains the vertexId
// with the given segmentId and segmentCopyIndex, or invalid<uint64_t>
// if no such vertex.
void AssemblyGraphSnapshot::createVertexTable(
    const PackedMarkerGraph& packedMarkerGraph)
{
    const uint64_t segmentCount = packedMarkerGraph.segments.size();

    // Create a temporary vertex table stored in a way that is easy to manipulate.
    vector< vector<uint64_t> > tmpVertexTable(segmentCount);
    for(uint64_t vertexId=0; vertexId<vertexVector.size(); vertexId++) {
        const Vertex& vertex = vertexVector[vertexId];
        auto& v = tmpVertexTable[vertex.segmentId];
        if(v.size() <= vertex.segmentReplicaIndex) {
            v.resize(vertex.segmentReplicaIndex + 1, invalid<uint64_t>);
        }
        v[vertex.segmentReplicaIndex] = vertexId;
    }

    // Now copy it to its permanent location.
    createNew(vertexTable, name + "-vertexTable");
    for(const auto& v: tmpVertexTable) {
        vertexTable.appendVector(v);
    }
}



// Use the vertex table to get a vertex id given segment id
// and replica index. If not found, returns false
// and also fills in a message string.
uint64_t AssemblyGraphSnapshot::getVertexId(
    uint64_t segmentId,
    uint64_t segmentReplicaIndex,
    string& message
) const
{
    uint64_t vertexId = invalid<uint64_t>;

    // Check that we have a valid segmentId.
    if(segmentId >= vertexTable.size()) {
        message = "Invalid segment id. Valid segment ids are 0 through "  +
            to_string(vertexTable.size()) + ".";
        return vertexId;
    }

    // Get the vertexIds for this segmentId.
    auto v = vertexTable[segmentId];

    // Extract the vertexId for this segmentReplicaIndex.
    if(segmentReplicaIndex >= v.size() or v[segmentReplicaIndex] == invalid<uint64_t> ) {
        message = "Invalid segment replica index. Valid segment replica indexes for this segment are:";
        for(uint64_t segmentReplicaIndex=0; segmentReplicaIndex<v.size(); segmentReplicaIndex++) {
            if(v[segmentReplicaIndex] != invalid<uint64_t>) {
                message = message + " " + to_string(segmentReplicaIndex);
            }
        }
        return vertexId;
    }
    vertexId = v[segmentReplicaIndex];
    SHASTA_ASSERT(vertexId < vertexVector.size());
    return vertexId;
}



uint64_t AssemblyGraphSnapshot::getVertexId(
    uint64_t segmentId,
    uint64_t segmentReplicaIndex
) const
{
    string message;
    return getVertexId(segmentId, segmentReplicaIndex, message);
}



// This accesses an existing snapshot.
AssemblyGraphSnapshot::AssemblyGraphSnapshot(
    const string& name,
    const MappedMemoryOwner& mappedMemoryOwner,
    const PackedMarkerGraph& packedMarkerGraph) :
    MappedMemoryOwner(mappedMemoryOwner),
    name(name),
    packedMarkerGraph(packedMarkerGraph)
{
    accessExistingReadOnly(vertexVector, name + "-vertices");
    accessExistingReadOnly(edgeVector, name + "-edges");
    accessExistingReadOnly(journeys, name + "-journeys");
    accessExistingReadOnly(vertexJourneyEntries, name + "-vertexJourneyEntries");
    accessExistingReadOnly(edgesBySource, name + "-edgesBySource");
    accessExistingReadOnly(edgesByTarget, name + "-edgesByTarget");
    accessExistingReadOnly(vertexTable, name + "-vertexTable");
    accessExistingReadOnly(assemblyPaths, name + "-assemblyPaths");
}



AssemblyGraphSnapshot::Vertex::Vertex(const AssemblyGraphVertex& vertex) :
    segmentId(vertex.segmentId),
    segmentReplicaIndex(vertex.segmentReplicaIndex),
    pathId(vertex.pathInformation.primaryInfo.pathId),
    positionInPath(vertex.pathInformation.primaryInfo.positionInPath),
    packedAssemblyGraphVertexId(vertex.packedAssemblyGraphVertexId),
    positionInPackedAssemblyGraph(vertex.positionInPackedAssemblyGraph)
{}



// Get the stringId for a given vertexId, or "None" if vertexId is invalid<uint64_t>.
string AssemblyGraphSnapshot::vertexStringId(uint64_t vertexId) const
{
    if(vertexId == invalid<uint64_t>) {
        return "None";
    } else {
        SHASTA_ASSERT(vertexId < vertexVector.size());
        const Vertex& vertex = vertexVector[vertexId];
        return vertex.stringId();
    }

}



void AssemblyGraphSnapshot::getEdgeTransitions(
    uint64_t edgeId,
    vector<Transition>& transitions) const
{
    transitions.clear();
    const Edge& edge = edgeVector[edgeId];
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;

    // Loop over journey entries of vertex0.
    for(const JourneyEntry& journeyEntry: vertexJourneyEntries[vertexId0]) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const uint64_t position0 = journeyEntry.position;
        const uint64_t position1 = position0 + 1;
        const auto journey = journeys[orientedReadId.getValue()];
        if(position1 < journey.size()) {
            if(journey[position1] == vertexId1) {
                transitions.push_back({orientedReadId, position0});
            }
        }
    }
}



uint64_t AssemblyGraphSnapshot::getEdgeCoverage(uint64_t edgeId) const
{
    uint64_t coverage = 0;

    const Edge& edge = edgeVector[edgeId];
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;

    // Loop over journey entries of vertex0.
    for(const JourneyEntry& journeyEntry: vertexJourneyEntries[vertexId0]) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const uint64_t position0 = journeyEntry.position;
        const uint64_t position1 = position0 + 1;
        const auto journey = journeys[orientedReadId.getValue()];
        if(position1 < journey.size()) {
            if(journey[position1] == vertexId1) {
                ++coverage;
            }
        }
    }

    return coverage;
}



// Find out if the segments of an edge are adjacent in the marker graph.
bool AssemblyGraphSnapshot::segmentsAreAdjacent(uint64_t edgeId) const
{
    // Get this edge.
    const Edge& edge = edgeVector[edgeId];

    // Get the corresponding vertices.
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;
    const Vertex& vertex0 = vertexVector[vertexId0];
    const Vertex& vertex1 = vertexVector[vertexId1];

    // Get the corresponding segments.
    const uint64_t segmentId0 = vertex0.segmentId;
    const uint64_t segmentId1 = vertex1.segmentId;

    // Get the marker graph paths of these segments.
    const auto path0 = packedMarkerGraph.segments[segmentId0];
    const auto path1 = packedMarkerGraph.segments[segmentId1];

    // The the last marker graph edge of path0
    // and the first marker graph edge of path1.
    const uint64_t markerGraphEdgeId0 = path0.back();
    const uint64_t markerGraphEdgeId1 = path1.front();
    const MarkerGraph::Edge& markerGraphEdge0 = packedMarkerGraph.markerGraph.edges[markerGraphEdgeId0];
    const MarkerGraph::Edge& markerGraphEdge1 = packedMarkerGraph.markerGraph.edges[markerGraphEdgeId1];

    return markerGraphEdge0.target == markerGraphEdge1.source;
}



void AssemblyGraphSnapshot::write() const
{
    writeJourneys();
    writeJourneyEntries();
    writeTransitions();

#if 0
    for(uint64_t minLinkCoverage=2; minLinkCoverage<=10; minLinkCoverage++) {
        writeGfa(minLinkCoverage);
    }
#endif
}



void AssemblyGraphSnapshot::writeGfa(uint64_t minLinkCoverage) const
{

    ofstream gfa(name + "-minLinkCoverage-" + to_string(minLinkCoverage) + ".gfa");

    // Write the header.
    gfa << "H\tVN:Z:1.0\n";

    // Write the vertices.
    for(const Vertex& vertex: vertexVector) {
        const span<const Base> sequence = packedMarkerGraph.segmentClippedSequence(vertex.segmentId);
        gfa <<"S\t" << vertex.stringId() << "\t";
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(gfa));
        gfa << "\n";
    }

    // Write the edges.
    for(uint64_t edgeId=0; edgeId<edgeVector.size(); edgeId++) {
        if(getEdgeCoverage(edgeId) < minLinkCoverage) {
            continue;
        }
        const Edge& edge = edgeVector[edgeId];
        gfa << "L\t" <<
            vertexVector[edge.vertexId0].stringId() << "\t+\t" <<
            vertexVector[edge.vertexId1].stringId() << "\t+\t0M\n";
    }

}



void AssemblyGraphSnapshot::writeJourneys() const
{
    ofstream csv(name + "-journeys.csv");

    for(uint64_t i=0; i<journeys.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto journey = journeys[orientedReadId.getValue()];

        csv << orientedReadId << ",";
        for(const uint64_t vertexId: journey) {
            csv << vertexVector[vertexId].stringId() << ",";
        }
        csv << "\n";
    }

}



void AssemblyGraphSnapshot::writeJourneyEntries() const
{
    ofstream csv(name + "-journeyEntries.csv");
    csv << "VertexId,SegmentId,Replica,OrientedReadId,Position\n";

    for(uint64_t vertexId=0; vertexId<vertexJourneyEntries.size(); vertexId++)
    {
        const Vertex& vertex = vertexVector[vertexId];
        const auto journeyEntries = vertexJourneyEntries[vertexId];

        for(const JourneyEntry& journeyEntry: journeyEntries) {
            csv << vertexId << ",";
            csv << vertex.segmentId << ",";
            csv << vertex.segmentReplicaIndex << ",";
            csv << journeyEntry.orientedReadId << ",";
            csv << journeyEntry.position << "\n";
        }
    }

}



void AssemblyGraphSnapshot::writeTransitions() const
{
    ofstream csv(name + "-transitions.csv");
    csv << "LinkId,SegmentReplica0,SegmentReplica1,OrientedReadId,Position0,Position1\n";

    vector<Transition> transitions;
    for(uint64_t edgeId=0; edgeId<edgeVector.size(); edgeId++) {
        const Edge& edge = edgeVector[edgeId];
        getEdgeTransitions(edgeId, transitions);
        for(const Transition& transition: transitions) {
            csv << edgeId << ",";
            csv << vertexVector[edge.vertexId0].stringId() << ",";
            csv << vertexVector[edge.vertexId1].stringId() << ",";
            csv << transition.orientedReadId << ",";
            csv << transition.position << ",";
            csv << transition.position + 1 << "\n";
        }

    }
}



void AssemblyGraphSnapshot::writeLinkTransitionsHtml(uint64_t linkId, ostream& html) const
{
    // Get some information about this link.
    const Edge& edge = edgeVector[linkId];
    const Vertex& vertex0 = vertexVector[edge.vertexId0];
    const Vertex& vertex1 = vertexVector[edge.vertexId1];
    vector<Transition> transitions;
    getEdgeTransitions(linkId, transitions);
    const auto leftPath = packedMarkerGraph.segments[vertex0.segmentId];
    const auto rightPath = packedMarkerGraph.segments[vertex1.segmentId];

    html << "<h2>Oriented read transitions from " <<
        vertexStringId(edge.vertexId0) << " to " <<
        vertexStringId(edge.vertexId1) << "</h2>";

    html<< "<table>"
        "<tr>"
        "<th>Oriented<br>read<br>id"
        "<th>Journey<br>Length"
        "<th>Position<br>of " << vertexStringId(edge.vertexId0) << "<br>in<br>journey"
        "<th>Position<br>of " << vertexStringId(edge.vertexId1) << "<br>in<br>journey";

    html << "<th>Left<br>skip<br>";
    writeInformationIcon(html,
        "Number of marker graph edges skipped at the end of " + vertexStringId(edge.vertexId0));

    html << "<th>Right<br>skip<br>";
    writeInformationIcon(html,
        "Number of marker graph edges skipped at the beginning of " + vertexStringId(edge.vertexId1));

    html << "<th>Left<br>last<br>ordinal<br>";
    writeInformationIcon(html,
        "Last ordinal on " + vertexStringId(edge.vertexId0));

    html << "<th>Right<br>first<br>ordinal<br>";
    writeInformationIcon(html,
        "First ordinal on " + vertexStringId(edge.vertexId1));

    html << "<th>Ordinal<br>skip"
        "<th>Estimated<br>link<br>separation";

    html << "<th>Sequence<br>";
    writeInformationIcon(html,
        "Sequence between left last ordinal and right first ordinal");



    // Write one row for each transition.
    for(const AssemblyGraphSnapshot::Transition& transition: transitions) {
        const OrientedReadId orientedReadId = transition.orientedReadId;
        const auto journey = packedMarkerGraph.journeys[orientedReadId.getValue()];
        const auto& leftJourneyStep = journey[transition.position];
        const auto& rightJourneyStep = journey[transition.position+1];
        const uint64_t leftSkip = leftPath.size() - leftJourneyStep.positions[1];
        const uint64_t rightSkip = rightJourneyStep.positions[0];
        const uint64_t leftOrdinal = leftJourneyStep.ordinals[1];
        const uint64_t rightOrdinal = rightJourneyStep.ordinals[0];
        SHASTA_ASSERT(rightOrdinal >= leftOrdinal);
        const uint64_t ordinalSkip = rightOrdinal - leftOrdinal;
        const int64_t linkSeparation =
            int64_t(ordinalSkip) -
            int64_t(leftSkip) -
            int64_t(rightSkip);
        html <<
            "<tr>"
            "<td class=centered>"
            "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
            "&strand=" << orientedReadId.getStrand() <<
            "'>"
            << transition.orientedReadId << "</a>"
            "<td class=centered>" << journey.size() <<
            "<td class=centered>" << transition.position <<
            "<td class=centered>" << transition.position+1 <<
            "<td class=centered>" << leftSkip <<
            "<td class=centered>" << rightSkip <<
            "<td class=centered>" << leftOrdinal <<
            "<td class=centered>" << rightOrdinal <<
            "<td class=centered>" << ordinalSkip <<
            "<td class=centered>" << linkSeparation;

        // Write the sequence.
        const auto orientedReadMarkers = packedMarkerGraph.markers[orientedReadId.getValue()];
        const CompressedMarker& marker0 = orientedReadMarkers[leftOrdinal];
        const CompressedMarker& marker1 = orientedReadMarkers[rightOrdinal];
        const uint64_t positionBegin = marker0.position;
        const uint64_t positionEnd = marker1.position + packedMarkerGraph.k;
        html << "<td class=centered style='font-family:courier'>";
        for(uint64_t position=positionBegin; position!=positionEnd; ++position) {
            html << packedMarkerGraph.reads.getOrientedReadBase(orientedReadId, uint32_t(position));
        }

    }
    html<< "</table>";
}


// Compute the consensus sequence for a link.
// Also compute the number of bases at the end
// of the left segment and at the beginning of the
// right segment that are overridden by the link
// and should be ignored during assembly.
void AssemblyGraphSnapshot::assembleLink(
    uint64_t linkId,
    ostream& html,
    vector<Base>& consensusSequence,
    uint64_t& leftOverride,
    uint64_t& rightOverride
    ) const
{
    // Get some information about this link.
    const Edge& edge = edgeVector[linkId];
    const Vertex& vertex0 = vertexVector[edge.vertexId0];
    const Vertex& vertex1 = vertexVector[edge.vertexId1];

    vector<Transition> transitions;
    getEdgeTransitions(linkId, transitions);

    const uint64_t segmentId0 = vertex0.segmentId;
    const uint64_t segmentId1 = vertex1.segmentId;

    const auto leftPath = packedMarkerGraph.segments[segmentId0];
    const auto rightPath = packedMarkerGraph.segments[segmentId1];

    const auto leftSegmentSequence = packedMarkerGraph.segmentSequences[segmentId0];
    const auto rightSegmentSequence = packedMarkerGraph.segmentSequences[segmentId1];

    const auto leftVertexOffsets = packedMarkerGraph.segmentVertexOffsets[segmentId0];
    const auto rightVertexOffsets = packedMarkerGraph.segmentVertexOffsets[segmentId1];

    const uint64_t k = packedMarkerGraph.k;


    // Page header.
    if(html) {
        html << "<h2>Assembly of link " << linkId << " from " <<
            vertexStringId(segmentId0) << " to " <<
            vertexStringId(segmentId1) << "</h2>";
    }

    // Loop over transitions to compute the maximum left/right skip,
    // that is, the maximum number of marker graph edges skipped
    // by an oriented read at the end of the left segment or
    // at the beginning of the right segment.
    // This is needed below to compute the sequences that participate in the MSA.
    uint64_t maxLeftSkip = 0;
    uint64_t maxRightSkip = 0;
    for(const AssemblyGraphSnapshot::Transition& transition: transitions) {
        const OrientedReadId orientedReadId = transition.orientedReadId;
        const auto journey = packedMarkerGraph.journeys[orientedReadId.getValue()];
        const auto& leftJourneyStep = journey[transition.position];
        const auto& rightJourneyStep = journey[transition.position+1];
        const uint64_t leftSkip = leftPath.size() - leftJourneyStep.positions[1];
        const uint64_t rightSkip = rightJourneyStep.positions[0];
        maxLeftSkip = max(maxLeftSkip, leftSkip);
        maxRightSkip = max(maxRightSkip, rightSkip);
    }

    // Compute the position in the left segment of the begining of the left segment
    // sequence that will be used to fill in MSA sequence.
    const uint64_t leftPositionBegin = leftVertexOffsets[leftVertexOffsets.size() -1 - maxLeftSkip];

    // Compute the position in the right segment of the end of the right segment
    // sequence that will be used to fill in MSA sequence.
    const uint64_t rightPositionEnd = rightVertexOffsets[maxRightSkip] + k;

    if(html) {
        html << "<h3>Contribution of each oriented read to the MSA</h3>";
        html << "Maximum left skip is " << maxLeftSkip << " markers, " <<
            leftSegmentSequence.size() - leftPositionBegin << " bases.";
        html << "<br>Maximum right skip is " << maxRightSkip << " markers, " <<
            rightPositionEnd << " bases.";
    }



    // Table header.
    if(html) {
        html<< "<p><table>"
            "<tr>"
            "<th>Oriented<br>read<br>id"
            "<th>Journey<br>Length"
            "<th>Position<br>of " << vertexStringId(segmentId0) << "<br>in<br>journey"
            "<th>Position<br>of " << vertexStringId(segmentId1) << "<br>in<br>journey";

        html << "<th>Left<br>skip<br>";
        writeInformationIcon(html,
            "Number of marker graph edges skipped at the end of " + vertexStringId(segmentId0));

        html << "<th>Right<br>skip<br>";
        writeInformationIcon(html,
            "Number of marker graph edges skipped at the beginning of " + vertexStringId(segmentId1));

        html << "<th>Left<br>last<br>ordinal<br>";
        writeInformationIcon(html,
            "Last ordinal on " + vertexStringId(edge.vertexId0));

        html << "<th>Right<br>first<br>ordinal<br>";
        writeInformationIcon(html,
            "First ordinal on " + vertexStringId(edge.vertexId1));

        html << "<th>Ordinal<br>skip"
            "<th>Estimated<br>link<br>separation";

        html << "<th>MSA sequence<br>";
        writeInformationIcon(html,
            "Sequence between left last ordinal and right first ordinal, "
            "extended out to the maximum left/right skip using the sequence of the "
            "left/right segments.");

        html << "<th>MSA<br>sequence<br>index";

        html << "<th>MSA<br>sequence<br>length";
    }

    // A vector to contain the distinct MSA sequence we found, each with the number of times it was found.
    vector< pair<vector<Base>, uint64_t> > msaSequences;

    // A vector that, for each transition, gives the index in msaSequence.
    vector<uint64_t> msaSequenceTable(transitions.size());



    // Loop over transitions to compute the MSA sequence to be used for each oriented read.
    vector<Base> msaSequence;
    for(uint64_t iTransition=0; iTransition<transitions.size(); iTransition++) {
        const AssemblyGraphSnapshot::Transition& transition = transitions[iTransition];
        const OrientedReadId orientedReadId = transition.orientedReadId;
        const auto journey = packedMarkerGraph.journeys[orientedReadId.getValue()];
        const auto& leftJourneyStep = journey[transition.position];
        const auto& rightJourneyStep = journey[transition.position+1];
        const uint64_t leftSkip = leftPath.size() - leftJourneyStep.positions[1];
        const uint64_t rightSkip = rightJourneyStep.positions[0];
        const uint64_t leftOrdinal = leftJourneyStep.ordinals[1];
        const uint64_t rightOrdinal = rightJourneyStep.ordinals[0];
        SHASTA_ASSERT(rightOrdinal >= leftOrdinal);
        const uint64_t ordinalSkip = rightOrdinal - leftOrdinal;
        const int64_t linkSeparation =
            int64_t(ordinalSkip) -
            int64_t(leftSkip) -
            int64_t(rightSkip);
        const auto orientedReadMarkers = packedMarkerGraph.markers[orientedReadId.getValue()];
        const CompressedMarker& marker0 = orientedReadMarkers[leftOrdinal];
        const CompressedMarker& marker1 = orientedReadMarkers[rightOrdinal];

        // Now we can compute the portion of the oriented read sequence that will
        // participate in the MSA.
        // This will be extended to the left/right as necessary,
        // using the sequence of the left/right segment.
        const uint64_t positionBegin = marker0.position;
        const uint64_t positionEnd = marker1.position + k;

        // Compute the position of the left extension in the left segment.
        const uint64_t leftPositionEnd = leftVertexOffsets[leftVertexOffsets.size() - 1 - leftSkip];

        // Compute the position of the right extension in the right segment.
        const uint64_t rightPositionBegin = rightVertexOffsets[rightSkip] + k;

        // Add the left extension to the MSA sequence.
        msaSequence.clear();
        for(uint64_t position=leftPositionBegin; position!=leftPositionEnd; ++position) {
            msaSequence.push_back(leftSegmentSequence[position]);
        }

        // Add the oriented read to the MSA sequence.
        for(uint64_t position=positionBegin; position!=positionEnd; ++position) {
            msaSequence.push_back(packedMarkerGraph.reads.getOrientedReadBase(orientedReadId, uint32_t(position)));
        }

        // Add the right extension to the MSA sequence.
        for(uint64_t position=rightPositionBegin; position!=rightPositionEnd; ++position) {
            msaSequence.push_back(rightSegmentSequence[position]);
        }

        // Update the msaSequences and msaSequenceTable.
        bool done = false;
        for(uint64_t i=0; i<msaSequences.size(); i++) {
            if(msaSequence == msaSequences[i].first) {
                msaSequenceTable[iTransition] = i;
                ++msaSequences[i].second;
                done = true;
                break;
            }
        }
        if(not done) {
            msaSequenceTable[iTransition] = msaSequences.size();
            msaSequences.push_back(make_pair(msaSequence, 1));
        }


        // Write a row to the table for this oriented read.
        if(html) {
            html <<
                "<tr>"
                "<td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() <<
                "'>"
                << transition.orientedReadId << "</a>"
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << transition.position <<
                "<td class=centered>" << transition.position+1 <<
                "<td class=centered>" << leftSkip <<
                "<td class=centered>" << rightSkip <<
                "<td class=centered>" << leftOrdinal <<
                "<td class=centered>" << rightOrdinal <<
                "<td class=centered>" << ordinalSkip <<
                "<td class=centered>" << linkSeparation;

            // Write the sequence.
            html << "<td class=centered style='font-family:courier'>";
            html << "<span style='background-color:LightGrey'>";
            for(uint64_t position=leftPositionBegin; position!=leftPositionEnd; ++position) {
                html << leftSegmentSequence[position];
            }
            html << "</span>";
            for(uint64_t position=positionBegin; position!=positionEnd; ++position) {
                html << packedMarkerGraph.reads.getOrientedReadBase(orientedReadId, uint32_t(position));
            }
            html << "<span style='background-color:LightGrey'>";
            for(uint64_t position=rightPositionBegin; position!=rightPositionEnd; ++position) {
                html << rightSegmentSequence[position];
            }
            html << "</span>";
            html << "<td class=centered>" << msaSequenceTable[iTransition];
            html << "<td class=centered>" << msaSequence.size();
        }

    }

    if(html) {
        html<< "</table>";
    }

    html <<
        "<h3>Distinct MSA sequences</h3>"
        "Found " << msaSequences.size() << " distinct MSA sequences."
        "<table>"
        "<tr><th>MSA<br>sequence<br>index<th>MSA Sequence<th>Length<th>Frequency";
    for(uint64_t msaSequenceIndex=0; msaSequenceIndex<msaSequences.size(); msaSequenceIndex++) {
        const auto& p = msaSequences[msaSequenceIndex];
        const vector<Base>& msaSequence = p.first;
        const uint64_t frequency = p.second;
        html << "<tr><td class=centered>" << msaSequenceIndex;
        html << "<td class=centered style='font-family:courier'>";
        copy(msaSequence.begin(), msaSequence.end(), ostream_iterator<Base>(html));
        html << "<td class=centered>" << msaSequence.size();
        html << "<td class=centered>" << frequency;
    }
    html << "</table>";

    // Compute the multiple sequence alignment.
    const uint64_t maxLength = 10000;
    linkMsaUsingSpoa(msaSequences, maxLength, html, consensusSequence);

    // Compute the number of bases at the beginning of the consensus sequence
    // that are identical to the corresponding bases in the left segment.
    uint64_t leftIdentical = 0;
    for(uint64_t i=leftPositionBegin; i<leftSegmentSequence.size(); i++) {
        if(i-leftPositionBegin >= consensusSequence.size()) {
            break;
        }
        if(leftSegmentSequence[i] == consensusSequence[i-leftPositionBegin]) {
            ++leftIdentical;
        } else {
            break;
        }
    }

    // Compute the number of bases at the end of the consensus sequence
    // that are identical to the corresponding bases in the right segment.
    uint64_t rightIdentical = 0;
    uint64_t j = rightPositionEnd - 1; // j is index in right segment sequence.
    for(uint64_t i=consensusSequence.size()-1; /* Check later */; i--, j--) { // i is index in consensus sequence.
        if(consensusSequence[i] == rightSegmentSequence[j]) {
            ++rightIdentical;
        } else {
            break;
        }
        if(i == 0) {
            break;
        }
        if(j == 0) {
            break;
        }
    }

    if(html) {
        html << "<br>" << leftIdentical << " bases at the beginning of the consensus sequence "
            "are identical to the corresponding bases in the left segment.";
        html << "<br>" << rightIdentical << " bases at the end of the consensus sequence "
            "are identical to the corresponding bases in the right segment.";
    }


    // Trim consensus bases at the end of the consensusSequence that are identical
    // to the corresponding bases in the right segment sequence.
    // Here, rightOverride is the number of bases at the beginning of the right segment
    // that are overridden by bases in the consensus sequence of the link.
    rightOverride = rightPositionEnd;
    while(not consensusSequence.empty()) {
        if(consensusSequence.back() == rightSegmentSequence[rightOverride-1]) {
            consensusSequence.resize(consensusSequence.size() - 1);
            --rightOverride;
            if(rightOverride == 0) {
                break;
            }
        } else {
            break;
        }
    }

    // Same, on the left.
    leftOverride = leftSegmentSequence.size() - leftPositionBegin;
    uint64_t leftTrim = 0;
    for(uint64_t i=0; i<consensusSequence.size(); i++) {
        const uint64_t j = i + leftPositionBegin;
        if(j >= leftSegmentSequence.size()) {
            break;
        }
        if(consensusSequence[i] == leftSegmentSequence[j]) {
            ++leftTrim;
            --leftOverride;
            if(leftOverride == 0) {
                break;
            }
        }
    }
    copy(consensusSequence.begin() + leftTrim, consensusSequence.end(), consensusSequence.begin());
    consensusSequence.resize(consensusSequence.size() - leftTrim);

    if(html) {
        html << "<h3>Link assembly results</h3>"
            "After trimming portions of consensus sequence identical to adjacent segment sequences, "
            "the consensus sequence has length " << consensusSequence.size() << ":"
            "<div style='font-family:courier'>";
        copy(consensusSequence.begin(), consensusSequence.end(), ostream_iterator<Base>(html));
        html << "</div>"
            "<br>Number of bases at the end of the left segment sequence overridden by link sequence: " <<
            leftOverride <<
            "<br>Number of bases at the beginning of the right segment sequence overridden by link sequence: " <<
            rightOverride;

        html << "<p>Assembly of the path consisting of this link plus the adjacent segments:"
            "<div style='font-family:courier'>";
        copy(leftSegmentSequence.begin(), leftSegmentSequence.end() - leftOverride,
            ostream_iterator<Base>(html));
        copy(consensusSequence.begin(), consensusSequence.end(),
            ostream_iterator<Base>(html));
        copy(rightSegmentSequence.begin()+ rightOverride, rightSegmentSequence.end(),
            ostream_iterator<Base>(html));
        html << "</div>";


    }
}


#

// Analyze the simple tangle at a given vertex vertexId1 by
// "following the reads" one step forward or backward.
// On return:
// - The previousVertices and nextVertices vectors
//   have size equal to the number of journey entries for vertexId1,
//   and contain the vertices vertexId0 and vertexId2 each of
//   the oriented reads visits immediately before/after visiting vertexId1.
//   This can be invalid<uint64_t> if vertexId1 is the
//   first or last vertex of a journey.
// - inCoverage[vertexId0] gives the number of oriented reads that
//   visit vertexId0 immediately before visiting vertexId1.
// - outCoverage[vertexId2] gives the number of oriented reads that
//   visit vertexId0 immediately after visiting vertexId1.
// - tangle_matrix[make_pair(vertexId0, vertexId2)] is
//   the number of oriented reads that visit vertexId0 immediately before vertexId1, and
//   visit vertexId2 immediately after  vertexId1.
// The three maps can include entries for which vertexId0 and/or vertexId2
// are invalid<uint64_t>, if vertexId1 is the first or last
// vertex of one or more oriented read journeys.
void AssemblyGraphSnapshot::analyzeSimpleTangleAtVertex(
    uint64_t vertexId1,
    vector<uint64_t>& inVertices,
    vector<uint64_t>& outVertices,
    std::map<uint64_t, uint64_t>& inCoverage,
    std::map<uint64_t, uint64_t>& outCoverage,
    std::map< pair<uint64_t, uint64_t>, uint64_t>& tangleMatrix) const
{
    // Access the journey entries for vertexId1.
    auto journeyEntries = vertexJourneyEntries[vertexId1];

    // Initialize output arguments.
    inVertices.clear();
    inVertices.reserve(journeyEntries.size());
    outVertices.clear();
    outVertices.reserve(journeyEntries.size());
    inCoverage.clear();
    outCoverage.clear();
    tangleMatrix.clear();

    // Loop over the journey entries for this vertex.
    for(const JourneyEntry& journeyEntry: journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const uint64_t position1 = journeyEntry.position;
        const auto journey = journeys[orientedReadId.getValue()];

        // Find the vertex this oriented read visits before vertexId1.
        uint64_t vertexId0 = invalid<uint64_t>;
        if(position1 > 0) {
            const uint64_t position0 = position1 - 1;
            vertexId0 = journey[position0];
        }

        // Find the vertex this oriented read visits after vertexId1.
        uint64_t vertexId2 = invalid<uint64_t>;
        if(position1 < journey.size() - 1) {
            const uint64_t position2 = position1 + 1;
            vertexId2 = journey[position2];
        }

        // Store vertexId0 and vertexId2.
        inVertices.push_back(vertexId0);
        outVertices.push_back(vertexId2);

        // Update the output maps.
        ++inCoverage[vertexId0];
        ++outCoverage[vertexId2];
        ++tangleMatrix[make_pair(vertexId0, vertexId2)];
    }
}



// Get deduplicated oriented read ids for a vertex.
// The journey entries are sorted by OrientedReadId,
// so this can be done quickly.
void AssemblyGraphSnapshot::getDeduplicatedOrientedReads(
    uint64_t vertexId,
    vector<OrientedReadId>& orientedReadIds
    ) const
{
    orientedReadIds.clear();
    for(const JourneyEntry& journeyEntry: vertexJourneyEntries[vertexId]) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        if(orientedReadIds.empty() or orientedReadId != orientedReadIds.back()) {
            orientedReadIds.push_back(orientedReadId);
        }
    }
}



// Compute the Jaccard similarity of the oriented read composition of two vertices.
// Duplicate oriented reads in the path entries for the vertices are ignored.
// This also computes the number of oriented reads in the union and intersection
// of the two read compositions, as well as vectors containing the deduplicated
// oriented reads for each of the two vertices.
double AssemblyGraphSnapshot::jaccard(
    uint64_t vertexId0,
    uint64_t vertexId1,
    vector<OrientedReadId>& orientedReadIds0,
    vector<OrientedReadId>& orientedReadIds1,
    vector<OrientedReadId>& unionOrientedReads,
    vector<OrientedReadId>& intersectionOrientedReads
) const
{
    // Get deduplicated OrientedReadId's.
    getDeduplicatedOrientedReads(vertexId0, orientedReadIds0);
    getDeduplicatedOrientedReads(vertexId1, orientedReadIds1);

    // Compute the union.
    unionOrientedReads.clear();
    std::set_union(
        orientedReadIds0.begin(), orientedReadIds0.end(),
        orientedReadIds1.begin(), orientedReadIds1.end(),
        back_inserter(unionOrientedReads));

    // Compute the intersection.
    intersectionOrientedReads.clear();
    std::set_intersection(
        orientedReadIds0.begin(), orientedReadIds0.end(),
        orientedReadIds1.begin(), orientedReadIds1.end(),
        back_inserter(intersectionOrientedReads));

    // Return the Jaccard similarity.
    return
        double(intersectionOrientedReads.size()) /
        double(unionOrientedReads.size());
}



void AssemblyGraphSnapshot::createSimpleAssemblyPath(
    const vector<uint64_t>& vertexIds,
    SimpleAssemblyPath& path) const
{
    SHASTA_ASSERT(vertexIds.size() >= 2);
    path.clear();
    ofstream html;  // Not open so no output takes place.

    // Store the vertex ids.
    // Vertices don't need to be assembled because their sequence
    // is stored in the PackedMarkerGraph.
    path.vertices.resize(vertexIds.size());
    for(uint64_t i=0; i<vertexIds.size(); i++) {
        SimpleAssemblyPath::Vertex& vertex = path.vertices[i];
        vertex.id = vertexIds[i];
        vertex.sequenceLength =
            packedMarkerGraph.segmentSequences[vertexVector[vertex.id].segmentId].size();
    }

    // Find the edge ids and assemble the edges (links).
    path.edges.resize(vertexIds.size() - 1);
    for(uint64_t i=1; i<vertexIds.size(); i++) {
        const uint64_t vertexId0 = path.vertices[i-1].id;
        const uint64_t vertexId1 = path.vertices[i].id;
        const auto outEdges0 = edgesBySource[vertexId0];
        uint64_t edgeId01 = invalid<uint64_t>;
        for(const uint64_t edgeId: outEdges0) {
            const Edge& edge = edgeVector[edgeId];
            if(edge.vertexId1 == vertexId1) {
                edgeId01 = edgeId;
                break;
            }
        }
        if(edgeId01 == invalid<uint64_t>) {
            throw runtime_error("There is no edge between vertex " + vertexStringId(vertexId0) +
                " and vertex " + vertexStringId(vertexId1));
        }

        // Fill in the information for this edge (link).
        SimpleAssemblyPath::Edge& edge = path.edges[i-1];
        edge.edgeId = edgeId01;
        assembleLink(edgeId01, html,
            edge.sequence, edge.leftOverride, edge.rightOverride);

        // Except for special cases (see below), the entire link sequence is
        // used to assemble the path.
        edge.sequenceBegin = 0;
        edge.sequenceEnd = edge.sequence.size();
    }


    // Use the edges leftOverride/rightOverride to choose the portion of
    // sequence of each vertex (segment) that will be used to assemble the path.
    path.vertices.front().sequenceBegin = 0;
    path.vertices.back().sequenceEnd = path.vertices.back().sequenceLength;
    for(uint64_t i=0; i<path.edges.size(); i++) {
        const auto& edge = path.edges[i];
        auto& vertex0 = path.vertices[i];
        auto& vertex1 = path.vertices[i+1];
        vertex0.sequenceEnd = vertex0.sequenceLength - edge.leftOverride;
        vertex1.sequenceBegin = edge.rightOverride;
    }

    for(const auto& vertex:path.vertices) {
        // This can happen and has to be handled.
        SHASTA_ASSERT(vertex.sequenceBegin <= vertex.sequenceEnd);
    }


    // Compute the begin position in the path sequence of the sequence portion
    // contributed by each segment or link.
    for(uint64_t i=0; /* Check later */ ; i++) {
        auto& segment = path.vertices[i];

        // Do it for this segment.
        if(i == 0) {
            segment.pathPosition = 0;
        } else {
            const auto& previousLink = path.edges[i-1];
            segment.pathPosition = previousLink.pathPosition +
                previousLink.sequenceEnd - previousLink.sequenceBegin;
        }

        // If this is the last segment, we are done.
        if(i == path.vertices.size() - 1) {
            break;
        }

        // Do it for this link.
        auto& link = path.edges[i];
        link.pathPosition = segment.pathPosition +
            segment.sequenceEnd - segment.sequenceBegin;

    }
}



void AssemblyGraphSnapshot::writeSimpleAssemblyPath(
    const SimpleAssemblyPath& path,
    ostream& html
) const
{
    html << "<h3>Assembly path with " << path.vertices.size() << " segments and " <<
        path.edges.size() << " links</h3>";

    vector<Base> assembledSequence;
    path.getAssembledSequence(packedMarkerGraph, *this, assembledSequence);
    html << "<h4>Assembled sequence </h4>" << assembledSequence.size() << " bases.<p>"
        "<div style='max-width:600px;overflow-wrap:break-word;font-family:courier'>";
    copy(assembledSequence.begin(), assembledSequence.end(),
        ostream_iterator<Base>(html));
    html << "</div>";

    html << "<h4>Assembled sequence breakdown</h4>"
        "This table uses red for segments and green for links.<p>"
        "<table><tr>"
        "<th>Segment<br>or<br>link"
        "<th>Position"
        "<th>Sequence<br><tr>";

    // Write the ids in a single cell.
    html << "<td class=centered>";
    for(uint64_t i=0; /* Check later */; i++) {
        const SimpleAssemblyPath::Vertex& vertex = path.vertices[i];
        html << "<span style='font-family:courier;color:Green'>" <<
            vertexStringId(vertex.id) << "</span><br>";
        if(i == path.vertices.size() -1) {
            break;
        }
        const SimpleAssemblyPath::Edge& edge = path.edges[i];
        html << "<span style='font-family:courier;color:Red'>" <<
            edge.edgeId << "</span><br>";
    }

    // Write the positions in a single cell.
    html << "<td class=centered>";
    for(uint64_t i=0; /* Check later */; i++) {
        const SimpleAssemblyPath::Vertex& vertex = path.vertices[i];
        html << "<span style='font-family:courier;color:Green'>" <<
            vertex.pathPosition << "</span><br>";
        if(i == path.vertices.size() -1) {
            break;
        }
        const SimpleAssemblyPath::Edge& edge = path.edges[i];
        html << "<span style='font-family:courier;color:Red'>" <<
            edge.pathPosition << "</span><br>";
    }

    // Write the sequence in a single cell.
    html << "<td>";
    for(uint64_t i=0; /* Check later */; i++) {
        const SimpleAssemblyPath::Vertex& vertex = path.vertices[i];
        const uint64_t segmentId = vertexVector[vertex.id].segmentId;
        const auto vertexSequence = packedMarkerGraph.segmentSequences[segmentId];
        html << "<span style='font-family:courier;color:Green'>";
        copy(vertexSequence.begin() + vertex.sequenceBegin, vertexSequence.begin() + vertex.sequenceEnd,
            ostream_iterator<Base>(html));
        html << "</span><br>";
        if(i == path.vertices.size() -1) {
            break;
        }
        const SimpleAssemblyPath::Edge& edge = path.edges[i];
        html << "<span style='font-family:courier;color:Red'>";
        copy(edge.sequence.begin() + edge.sequenceBegin, edge.sequence.begin() + edge.sequenceEnd,
            ostream_iterator<Base>(html));
        html << "</span><br>";
    }


    html << "</table>";



    html <<
        "<h4>Assembly details</h4>"
        "<table><tr>"
        "<th>Segment<br>id"
        "<th>Link<br>id"
        "<th>Total<br>sequence<br>length<br>";
    writeInformationIcon(html,
        "Total length of the sequence of this segment or link, including "
        "the portion not used to assemble this path.");
    html << "<th>Path<br>sequence<br>length<br>";
    writeInformationIcon(html,
        "Length of the sequence portion of this segment or link used to assemble this path.");
    html << "<th>Link<br>or<br>segment<br>begin<br>";
    writeInformationIcon(html,
        "Begin position in the sequence of this segment or link of "
        "the sequence portion used to assemble the path.");
    html << "<th>Link<br>or<br>segment<br>end<br>";
    writeInformationIcon(html,
        "End position in the sequence of this segment or link of "
        "the sequence portion used to assemble the path. "
        "End position is one past the end.");
    html << "<th>Path<br>begin<br>";
    writeInformationIcon(html,
        "Begin position in the assembled path sequence of "
        "the sequence portion contributed by this segment or link.");
    html << "<th>Path<br>end<br>";
    writeInformationIcon(html,
        "End position in the assembled path sequence of "
        "the sequence portion contributed by this segment or link. "
        "End position is one past the end.");
    html << "<th style='max-width:600px'>Sequence<br>";
    writeInformationIcon(html,
        "The portion greyed out is not used for assembly.");

    // Loop over all segments and links.
    for(uint64_t i=0; /* Check later */; i++) {

        // Write the vertex (segment).
        const SimpleAssemblyPath::Vertex& vertex = path.vertices[i];
        const uint64_t segmentId = vertexVector[vertex.id].segmentId;

        html << "<tr>"
            "<td class=centered>" << vertexStringId(vertex.id) <<
            "<td>"
            "<td class=centered>" << vertex.sequenceLength <<
            "<td class=centered>" << vertex.sequenceEnd - vertex.sequenceBegin <<
            "<td class=centered>" << vertex.sequenceBegin <<
            "<td class=centered>" << vertex.sequenceEnd <<
            "<td class=centered>" << vertex.pathPosition <<
            "<td class=centered>" << vertex.pathPosition + vertex.sequenceEnd - vertex.sequenceBegin <<
            "<td class=centered style='max-width:400px;overflow-wrap:break-word;font-family:courier'>";
        const auto vertexSequence = packedMarkerGraph.segmentSequences[segmentId];
        html << "<span style='background-color:LightGrey'>";
        copy(vertexSequence.begin(), vertexSequence.begin() + vertex.sequenceBegin,
            ostream_iterator<Base>(html));
        html << "</span>";
        copy(vertexSequence.begin() + vertex.sequenceBegin, vertexSequence.begin() + vertex.sequenceEnd,
            ostream_iterator<Base>(html));
        html << "<span style='background-color:LightGrey'>";
        copy(vertexSequence.begin() + vertex.sequenceEnd, vertexSequence.end(),
            ostream_iterator<Base>(html));
        html << "</span>";



        if(i == path.edges.size()) {
            break;
        }

        // Write the edge (link).
        const SimpleAssemblyPath::Edge& edge = path.edges[i];
        html << "<tr>"
            "<td>" <<
            "<td class=centered>" << edge.edgeId <<
            "<td class=centered>" << edge.sequence.size() <<
            "<td class=centered>" << edge.sequenceEnd - edge.sequenceBegin <<
            "<td class=centered>" << edge.sequenceBegin <<
            "<td class=centered>" << edge.sequenceEnd <<
            "<td class=centered>" << edge.pathPosition <<
            "<td class=centered>" << edge.pathPosition + edge.sequenceEnd - edge.sequenceBegin <<
            "<td class=centered style='max-width:600px;overflow-wrap:break-word;font-family:courier'>";
        html << "<span style='background-color:LightGrey'>";
        copy(edge.sequence.begin(), edge.sequence.begin() + edge.sequenceBegin,
            ostream_iterator<Base>(html));
        html << "</span>";
        copy(edge.sequence.begin() + edge.sequenceBegin, edge.sequence.begin() + edge.sequenceEnd,
            ostream_iterator<Base>(html));
        html << "<span style='background-color:LightGrey'>";
        copy(edge.sequence.begin() + edge.sequenceEnd, edge.sequence.end(),
            ostream_iterator<Base>(html));
        html << "</span>";
    }

    html << "</table>";
}



void AssemblyGraphSnapshot::SimpleAssemblyPath::clear()
{
    vertices.clear();
    edges.clear();
}



void AssemblyGraphSnapshot::SimpleAssemblyPath::getAssembledSequence(
    const PackedMarkerGraph& packedMarkerGraph,
    const AssemblyGraphSnapshot& snapshot,
    vector<Base>& sequence) const
{
    sequence.clear();
    for(uint64_t i=0; /* Check later */; i++) {
        const SimpleAssemblyPath::Vertex& vertex = vertices[i];
        const uint64_t segmentId = snapshot.vertexVector[vertex.id].segmentId;
        const auto vertexSequence = packedMarkerGraph.segmentSequences[segmentId];
        copy(vertexSequence.begin() + vertex.sequenceBegin, vertexSequence.begin() + vertex.sequenceEnd,
            back_inserter(sequence));
        if(i == vertices.size() -1) {
            break;
        }
        const SimpleAssemblyPath::Edge& edge = edges[i];
        copy(edge.sequence.begin() + edge.sequenceBegin, edge.sequence.begin() + edge.sequenceEnd,
            back_inserter(sequence));
    }
}


// Return true if the specified edge edgeId = v0->v1
// satisfies both of the following two conditions:
// - It has maximum coverage among the out-edges of v0.
// - It has maximum coverage among the in-edges of v1.
// This could be slow because of repeated calls to getEdgeCoverage.
// If this becomes a problem, store edge coverage in the AssemblyGraphSnapshot::Edge.
bool AssemblyGraphSnapshot::isMaximalCoverageEdge(uint64_t edgeId) const
{
    const Edge& edge = edgeVector[edgeId];
    const uint64_t coverage = getEdgeCoverage(edgeId);
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;
    const auto outEdges0 = edgesBySource[vertexId0];
    const auto inEdges1  = edgesByTarget[vertexId1];

    // Check the out-edges of vertexId0.
    uint64_t maxCoverage0 = 0;
    for(const uint64_t edgeId0: outEdges0) {
        const uint64_t coverage0 = getEdgeCoverage(edgeId0);
        maxCoverage0 = max(maxCoverage0, coverage0);
    }
    SHASTA_ASSERT(maxCoverage0 > 0);
    if(coverage < maxCoverage0) {
        // edgeId is not one of the outEdges0 with maximum coverage.
        return false;
    }

    // Check the in-edges of vertexId10.
    uint64_t maxCoverage1 = 0;
    for(const uint64_t edgeId1: inEdges1) {
        const uint64_t coverage1 = getEdgeCoverage(edgeId1);
        maxCoverage1 = max(maxCoverage1, coverage1);
    }
    SHASTA_ASSERT(maxCoverage1 > 0);
    if(coverage < maxCoverage1) {
        // edgeId is not one of the inEdges1 with maximum coverage.
        return false;
    }

    return true;
}

// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<shasta::mode3a::AssemblyGraph>;



AssemblyGraph::AssemblyGraph(
    const PackedMarkerGraph& packedMarkerGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    packedMarkerGraph(packedMarkerGraph)
{
    createSegmentsAndJourneys();
    createLinks();
}



void AssemblyGraph::createSegmentsAndJourneys()
{
    AssemblyGraph& assemblyGraph = *this;

    // Initially, create a vertex for each segment in the packedMarkerGraph.
    verticesBySegment.clear();
    verticesBySegment.resize(packedMarkerGraph.segments.size());
    for(uint64_t segmentId=0; segmentId<packedMarkerGraph.segments.size(); segmentId++) {
        const vertex_descriptor v = add_vertex(AssemblyGraphVertex(segmentId), assemblyGraph);
        verticesBySegment[segmentId].push_back(v);
    }

    // The journey of an oriented read is the sequence of segments
    // visited by the oriented read.
    // Initially, we construct it from the corresponding journey
    // in the PackedMarkerGraph.
    // While constructing the journeys, we also store journey entries in the vertices.
    journeys.resize(packedMarkerGraph.journeys.size());
    for(uint64_t i=0; i<journeys.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto packedMarkerGraphJourney = packedMarkerGraph.journeys[i];
        auto& journey = journeys[i];
        for(uint64_t position=0; position<packedMarkerGraphJourney.size(); position++) {
            const uint64_t segmentId = packedMarkerGraphJourney[position].segmentId;
            const vertex_descriptor v = verticesBySegment[segmentId].front();
            journey.push_back(v);
            assemblyGraph[v].journeyEntries.push_back({orientedReadId, position});
        }
    }
}



// Get the stringId for a given vertex_descriptor, or "None" if v is null_vertex().
string AssemblyGraph::vertexStringId(vertex_descriptor v) const
{
    if(v == null_vertex()) {
        return "None";
    } else {
        return (*this)[v].stringId();
    }
}




void AssemblyGraph::createLinks()
{
    AssemblyGraph& assemblyGraph = *this;

    // Gather transitions for all oriented reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitions;
    for(uint64_t i=0; i<journeys.size(); i++) {

        // Loop over the journey for this oriented read.
        const auto& journey = journeys[i];
        for(uint64_t position1=1; position1<journey.size(); position1++) {
            const uint64_t position0 = position1 - 1;
            const vertex_descriptor v0 = journey[position0];
            const vertex_descriptor v1 = journey[position1];
            // This should be done at the beginning, when the journeys
            // don't contain any null vertices.
            SHASTA_ASSERT(v0 != null_vertex());
            SHASTA_ASSERT(v1 != null_vertex());
            transitions.push_back({v0, v1});
        }
    }
    deduplicate(transitions);

    for(const pair<vertex_descriptor, vertex_descriptor>& transition: transitions) {
        const vertex_descriptor v0 = transition.first;
        const vertex_descriptor v1 = transition.second;
        add_edge(v0, v1, assemblyGraph);
    }
}



// Find out if two segments are adjacent in the marker graph.
bool AssemblyGraph::segmentsAreAdjacent(edge_descriptor e) const
{
    const AssemblyGraph& assemblyGraph = *this;
    return segmentsAreAdjacent(
        source(e, assemblyGraph),
        target(e, assemblyGraph));
}



bool AssemblyGraph::segmentsAreAdjacent(
    vertex_descriptor v0,
    vertex_descriptor v1) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const uint64_t segmentId0 = assemblyGraph[v0].segmentId;
    const uint64_t segmentId1 = assemblyGraph[v1].segmentId;

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



// Get the transitions for an edge.
void AssemblyGraph::getEdgeTransitions(
    edge_descriptor e,
    vector<Transition>& transitions) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Access the vertices of this edge.
    const vertex_descriptor v0 = source(e, assemblyGraph);
    const vertex_descriptor v1 = target(e, assemblyGraph);
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];



    // Loop over journey entries of vertex1.
    transitions.clear();
    for(const JourneyEntry& journeyEntry: vertex1.journeyEntries) {
        const uint64_t position1 = journeyEntry.position;
        if(position1 == 0) {
            // v1 is at the beginning of the journey.
            // There is no previous vertex in the journey.
            continue;
        }
        const uint64_t position0 = position1 - 1;

        // Access the journey for this oriented read.
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const vector<vertex_descriptor>& journey = journeys[orientedReadId.getValue()];

        // If the previous entry is not on v0, this does not
        // correspond to a transiton for the edge we are working on.
        if(journey[position0] != v0) {
            continue;
        }

        // Store this transition.
        transitions.push_back({position0, position1, orientedReadId});
    }
}



// This is similar to getTransitions above, but it
// just counts the transitions instead of storing them.
uint64_t AssemblyGraph::edgeCoverage(edge_descriptor e) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Access the vertices of this edge.
    const vertex_descriptor v0 = source(e, assemblyGraph);
    const vertex_descriptor v1 = target(e, assemblyGraph);
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];



    // Loop over journey entries of vertex1.
    uint64_t coverage = 0;
    for(const JourneyEntry& journeyEntry: vertex1.journeyEntries) {
        const uint64_t position1 = journeyEntry.position;
        if(position1 == 0) {
            // v1 is at the beginning of the journey.
            // There is no previous vertex in the journey.
            continue;
        }
        const uint64_t position0 = position1 - 1;

        // Access the journey for this oriented read.
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const vector<vertex_descriptor>& journey = journeys[orientedReadId.getValue()];

        // If the previous entry is not on v0, this does not
        // correspond to a transiton for the edge we are working on.
        if(journey[position0] != v0) {
            continue;
        }

        ++coverage;
    }

    return coverage;
}



void AssemblyGraph::write(const string& name) const
{
    for(uint64_t minLinkCoverage=2; minLinkCoverage<=6; minLinkCoverage++) {
        writeGfa(name + "-minLinkCoverage-" + to_string(minLinkCoverage) + ".gfa", minLinkCoverage);
    }
    writeLinkCoverageHistogram(name + "-LinkCoverageHistogram.csv");
    writeJourneys(name + "-journeys.csv");
}



void AssemblyGraph::writeLinkCoverageHistogram(const string& name) const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector<uint64_t> histogram;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const uint64_t coverage = edgeCoverage(e);
        if(histogram.size() <= coverage) {
            histogram.resize(coverage+1, 0);
        }
        ++histogram[coverage];
    }

    ofstream csv(name);
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



void AssemblyGraph::writeJourneys(const string& name) const
{
    const AssemblyGraph& assemblyGraph = *this;
    ofstream csv(name);

    for(uint64_t i=0; i<journeys.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto journey = journeys[orientedReadId.getValue()];

        csv << orientedReadId << ",";
        for(const vertex_descriptor v: journey) {
            csv << assemblyGraph.vertexStringId(v) << ",";
        }
        csv << "\n";
    }

}



void AssemblyGraph::writeGfa(const string& name, uint64_t minLinkCoverage) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream gfa(name);

    // Write the headers.
    gfa << "H\tVN:Z:1.0\n";

    // Write the segments.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        gfa <<"S\t" << assemblyGraph[v].stringId() << "\t*\n";
    }

    // Write the links.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(edgeCoverage(e) < minLinkCoverage) {
            continue;
        }
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        gfa << "L\t" <<
            assemblyGraph[v0].stringId() << "\t+\t" <<
            assemblyGraph[v1].stringId() << "\t+\t0M\n";
    }
}



void AssemblyGraph::simpleDetangle(
    uint64_t minLinkCoverage,
    uint64_t minTangleCoverage)
{
    AssemblyGraph& assemblyGraph = *this;

    vector<vertex_descriptor> allVertices;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        allVertices.push_back(v);
    }
    for(const vertex_descriptor v: allVertices) {
        simpleDetangle(v, minLinkCoverage, minTangleCoverage);
    }
}



void AssemblyGraph::simpleDetangle(
    vertex_descriptor v1,
    uint64_t minLinkCoverage,
    uint64_t minTangleCoverage)
{
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];

    // Find adjacent vertices by following the reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > adjacentVertices;
    findAdjacentVertices(v1, adjacentVertices);

    // Group them.
    // Store journey entry indexes, that is indexes into the journeyEntries vector
    // of v1 and into the adjacentVertices vector.
    std::map<vertex_descriptor, vector<uint64_t> > map0;  // By previous vertex
    std::map<vertex_descriptor, vector<uint64_t> > map2;  // By next vertex.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<uint64_t> > map02;    // By previous and next vertex.
    for(uint64_t i=0; i<adjacentVertices.size(); i++) {
        const pair<vertex_descriptor, vertex_descriptor>& v02 = adjacentVertices[i];
        const vertex_descriptor v0 = v02.first;
        const vertex_descriptor v2 = v02.second;
        map0[v0].push_back(i);
        map2[v2].push_back(i);
        map02[v02].push_back(i);
    }

    // For detangling, we only consider parent/children
    // that are not null_vertex() and with at least minLinkCoverage
    // oriented reads.
    vector<vertex_descriptor> parents;
    for(const auto& p: map0) {
        const vertex_descriptor v0 = p.first;
        if(v0 != null_vertex() and p.second.size() >= minLinkCoverage) {
            parents.push_back(v0);
        }
    }
    vector<vertex_descriptor> children;
    for(const auto& p: map2) {
        const vertex_descriptor v2 = p.first;
        if(v2 != null_vertex() and p.second.size() >= minLinkCoverage) {
            children.push_back(v2);
        }
    }

    // For now we only attempt detangling if there are at least two
    // parents and two children.
    if(parents.size() < 2 or children.size() < 2){
        return;
    }

    if(debug) {
        cout << "Detangling " << vertex1.stringId() << " with " <<
            parents.size() << " parents and " <<
            children.size() << " children.\n";
        cout << "Parents: ";
        for(const vertex_descriptor parent: parents) {
            cout << " " << assemblyGraph[parent].stringId();
        }
        cout << "\n";
        cout << "Children: ";
        for(const vertex_descriptor child: children) {
            cout << " " << assemblyGraph[child].stringId();
        }
        cout << "\n";

        if(false) {
            cout << "Details of journey entries:\n";
            for(uint64_t i=0; i<vertex1.journeyEntries.size(); i++) {
                cout << i << " " << vertex1.journeyEntries[i].orientedReadId << " " <<
                    assemblyGraph[adjacentVertices[i].first].stringId() << " " <<
                    assemblyGraph[adjacentVertices[i].second].stringId() << "\n";
            }
        }
    }


    // Find the pairs that will generate a new vertex.
    // These are called the "active pairs" here.
    // They are the ones for which map02 contains at least
    // minTangleCoverage entries.
    vector< pair<vertex_descriptor, vertex_descriptor> > activePairs;
    for(const auto& p: map02) {
        const auto& v02 = p.first;
        const vertex_descriptor v0 = v02.first;
        const vertex_descriptor v2 = v02.second;
        if(v0 != null_vertex() and v2 != null_vertex() and p.second.size() >= minTangleCoverage) {
            activePairs.push_back(v02);
        }
    }


    // Each active pair generates a new vertex with the same segmentId as v1.
    vector<vertex_descriptor> newVertices;
    const uint64_t segmentId1 = vertex1.segmentId;
    for(uint64_t i=0; i<activePairs.size(); i++) {
        const vertex_descriptor vNew = add_vertex(
            AssemblyGraphVertex(segmentId1, verticesBySegment[segmentId1].size()),
            assemblyGraph);
        verticesBySegment[segmentId1].push_back(vNew);
        newVertices.push_back(vNew);
    }
    if(debug) {
        cout << "Active pairs:\n";
        for(uint64_t i=0; i<activePairs.size(); i++) {
            const auto& v02 = activePairs[i];
            const vertex_descriptor v0 = v02.first;
            const vertex_descriptor v2 = v02.second;
            cout << assemblyGraph[v0].stringId() << " " <<
                assemblyGraph[v2].stringId() << " new vertex " <<
                assemblyGraph[newVertices[i]].stringId() << "\n";
        }
    }



    // In addition, we create a vertex that will receive journey entries
    // that are not in active pairs.
    const vertex_descriptor vNew = add_vertex(
        AssemblyGraphVertex(segmentId1, verticesBySegment[segmentId1].size()),
        assemblyGraph);
    verticesBySegment[segmentId1].push_back(vNew);
    if(false) {
        cout << "New vertex not associated with an active pair: " << assemblyGraph[vNew].stringId() << endl;
    }



    // Assign each journey entry of v1 to one of the new vertices we just created.
    for(uint64_t i=0; i<vertex1.journeyEntries.size(); i++) {
        const JourneyEntry& journeyEntry = vertex1.journeyEntries[i];
        const auto& v02 = adjacentVertices[i];
        const vertex_descriptor v0 = v02.first;
        const vertex_descriptor v2 = v02.second;
        if(false) {
            cout << "Assigning journey entry " << i << " " << flush <<
            assemblyGraph[v0].stringId() << " " << flush <<
            assemblyGraph[v2].stringId() << " to a new vertex." << endl;
        }

        // Look it up in the active pairs.
        auto it = find(activePairs.begin(), activePairs.end(), v02);

        // Find the vertex we are going to add this JourneyEntry to.
        const vertex_descriptor v = (it == activePairs.end()) ? vNew : newVertices[it - activePairs.begin()];
        if(false) {
            cout << "This journey entry will be assigned to vertex " << assemblyGraph[v].stringId() << endl;
        }

        // Add the journey entry to this vertex.
        assemblyGraph[v].journeyEntries.push_back(journeyEntry);

        // Update the journey of this oriented read to reflect this change.
        journeys[journeyEntry.orientedReadId.getValue()][journeyEntry.position] = v;

        // Make sure the edges v0->v and v->v2 exist.
        if(v0 != null_vertex()) {
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v, assemblyGraph);
            if(not edgeExists) {
                add_edge(v0, v, assemblyGraph);
            }
        }
        if(v2 != null_vertex()) {
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v, v2, assemblyGraph);
            if(not edgeExists) {
                add_edge(v, v2, assemblyGraph);
            }
        }
    }

    // Now we can remove v1.
    verticesBySegment[segmentId1][vertex1.segmentReplicaIndex] = null_vertex();
    clear_vertex(v1, assemblyGraph);
    remove_vertex(v1, assemblyGraph);

}



// Find the previous and next vertex for each JourneyEntry in a given vertex.
// On return, adjacentVertices contains a pair of vertex descriptors for
// each JourneyEntry in vertex v, in the same order.
// Those vertex descriptors are the previous and next vertex visited
// by the oriented read for that JourneyEntry, and can be null_vertex()
// if v is at the beginning or end of the journey of an oriented read.
void AssemblyGraph::findAdjacentVertices(
    vertex_descriptor v,
    vector< pair<vertex_descriptor, vertex_descriptor> >& adjacentVertices
) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphVertex& vertex = assemblyGraph[v];

    adjacentVertices.clear();
    for(const JourneyEntry& journeyEntry: vertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const vector<vertex_descriptor>& journey = journeys[orientedReadId.getValue()];

        const uint64_t position1 = journeyEntry.position;

        // Find the previous vertex visited by this journey, if any.
        vertex_descriptor v0 = null_vertex();
        if(position1 > 0) {
            const uint64_t position0 = position1 - 1;
            v0 = journey[position0];
        }

        // Find the next vertex visited by this journey, if any.
        vertex_descriptor v2 = null_vertex();
        const uint64_t position2 = position1 + 1;
        if(position2 < journey.size()) {
            v2 = journey[position2];
        }

        // Store this pair of adjacent vertices.
        adjacentVertices.push_back({v0, v2});
    }
}




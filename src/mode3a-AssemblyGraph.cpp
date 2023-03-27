// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-JaccardGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "shastaLapack.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Spoa.
#include "spoa/spoa.hpp"

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
    // Because of the way they are constructed, the journey entries of a vertex
    // are ordered by OrientedReadId and then by position.
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
            if(v0 == null_vertex()) {
                continue;
            }
            const vertex_descriptor v1 = journey[position1];
            if(v1 == null_vertex()) {
                continue;
            }
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



// Find descendants of a given vertex up to a specified distance.
// This is done by following the journeys.
void AssemblyGraph::findDescendants(
    vertex_descriptor v,
    uint64_t distance,
    vector<vertex_descriptor>& descendants) const
{
    const AssemblyGraph& assemblyGraph = *this;

    descendants.clear();
    const AssemblyGraphVertex& vertex = assemblyGraph[v];

    for(const JourneyEntry& journeyEntry: vertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const vector<vertex_descriptor>& journey = journeys[orientedReadId.getValue()];

        const uint64_t firstPosition = journeyEntry.position + 1;
        const uint64_t lastPosition  = min(journeyEntry.position + distance, journey.size() - 1);

        for(uint64_t position=firstPosition; position<=lastPosition; position++) {
            const vertex_descriptor u = journey[position];
            if(u != null_vertex()) {
                descendants.push_back(journey[position]);
            }
        }
    }
    deduplicate(descendants);
}



// Find ancestors of a given vertex up to a specified distance.
// This is done by following the journeys.
void AssemblyGraph::findAncestors(
    vertex_descriptor v,
    uint64_t distance,
    vector<vertex_descriptor>& ancestors) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ancestors.clear();
    const AssemblyGraphVertex& vertex = assemblyGraph[v];

    for(const JourneyEntry& journeyEntry: vertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const vector<vertex_descriptor>& journey = journeys[orientedReadId.getValue()];
        const int64_t vPosition = int64_t(journeyEntry.position);

        const int64_t firstPosition = max(vPosition - int64_t(distance), 0L);
        const int64_t lastPosition  = vPosition - 1;

        for(int64_t position=firstPosition; position<=lastPosition; position++) {
            const vertex_descriptor u = journey[position];
            if(u != null_vertex()) {
                ancestors.push_back(journey[position]);
            }
        }
    }
    deduplicate(ancestors);
}



// Create a detangled AssemblyGraph using tangle matrices to split vertices
// of another AssemblyGraph.
// Only tangle matrix entries that are at least equal to minCoverage as used.
// As a result, the new AssemblyGraph can have missing journey entries.
// That is, some journey entries will remain set to null_vertex().
AssemblyGraph::AssemblyGraph(
    DetangleUsingTangleMatrices,
    const AssemblyGraph& oldAssemblyGraph,
    uint64_t minCoverage) :
    MultithreadedObject<AssemblyGraph>(*this),
    packedMarkerGraph(oldAssemblyGraph.packedMarkerGraph)
{
    createFromTangledMatrices(oldAssemblyGraph, minCoverage);

}
void AssemblyGraph::createFromTangledMatrices(
    const AssemblyGraph& oldAssemblyGraph,
    uint64_t minCoverage)
{
    AssemblyGraph& newAssemblyGraph = *this;

    // Initialize verticesBySegment.
    verticesBySegment.clear();
    verticesBySegment.resize(packedMarkerGraph.segments.size());

    // Initialize  oriented reads journeys.
    journeys.clear();
    journeys.resize(packedMarkerGraph.journeys.size());
    for(uint64_t i=0; i<journeys.size(); i++) {
        const auto packedMarkerGraphJourney = packedMarkerGraph.journeys[i];
        journeys[i].resize(packedMarkerGraphJourney.size(), null_vertex());
    }



    // Loop over vertices of the old assembly graph.
    // Each of them can create one or more vertices in the new assembly graph.
    BGL_FORALL_VERTICES(vOld, oldAssemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& oldVertex = oldAssemblyGraph[vOld];
        const uint64_t segmentId = oldVertex.segmentId;

        // The tangle matrix contains journey entries index for each pair
        // (previous vertex, next vertex).
        std::map< pair<vertex_descriptor, vertex_descriptor>, vector<uint64_t> > tangleMatrix;

        // To construct the tangle matrix,
        // loop over journey entries of this vertex of the old assembly graph.
        for(uint64_t journeyEntryIndex=0; journeyEntryIndex<oldVertex.journeyEntries.size(); journeyEntryIndex++) {
            const JourneyEntry& journeyEntry = oldVertex.journeyEntries[journeyEntryIndex];
            const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
            const uint64_t position = journeyEntry.position;

            // Access the journey of this oriented read.
            const vector<vertex_descriptor>& journey = oldAssemblyGraph.journeys[orientedReadId.getValue()];

            // Find the previous vertex in the journey.
            vertex_descriptor vPrevious = null_vertex();
            if(position > 0) {
                vPrevious = journey[position - 1];
            };

            // Find the next vertex in the journey.
            vertex_descriptor vNext = null_vertex();
            if(position < journey.size() - 1) {
                vNext = journey[position + 1];
            };

            // Update the tangle matrix.
            tangleMatrix[make_pair(vPrevious, vNext)].push_back(journeyEntryIndex);
        }



        // Each tangle matrix entry at least equal to minCoverage generates a vertex
        // in the new AssemblyGraph.
        for(const auto& p: tangleMatrix) {
            const vector<uint64_t>& journeyEntryIndexes = p.second;
            if(journeyEntryIndexes.size() < minCoverage) {
                continue;
            }

            // This tangle matrix entry is large enough.
            // Create a new vertex.
            const vertex_descriptor vNew = boost::add_vertex(
                AssemblyGraphVertex(segmentId, verticesBySegment[segmentId].size()),
                newAssemblyGraph);
            verticesBySegment[segmentId].push_back(vNew);
            AssemblyGraphVertex& newVertex = newAssemblyGraph[vNew];

            // Add these journey entries to the vertex.
            for(const uint64_t journeyEntryIndex: journeyEntryIndexes) {
                const JourneyEntry& journeyEntry = oldVertex.journeyEntries[journeyEntryIndex];
                newVertex.journeyEntries.push_back(journeyEntry);

                // Update the journey for this oriented read.
                const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
                const uint64_t position = journeyEntry.position;
                vector<vertex_descriptor>& orientedReadJourney = journeys[orientedReadId.getValue()];
                SHASTA_ASSERT(orientedReadJourney[position] == null_vertex());
                orientedReadJourney[position] = vNew;
            }
        }

    }

    createLinks();
}



// Detangle by local clustering.
// Create a detangled AssemblyGraph using clustering on another AssemblyGraph.
AssemblyGraph::AssemblyGraph(
    DetangleUsingLocalClustering,
    const AssemblyGraph& oldAssemblyGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    packedMarkerGraph(oldAssemblyGraph.packedMarkerGraph)
{
    createByLocalClustering(oldAssemblyGraph);
}
void AssemblyGraph::createByLocalClustering(
    const AssemblyGraph& oldAssemblyGraph)
{
    // EXPOSE WHEN CODE STABILIZES
    const uint64_t distance = 2;
    const double singularValueThreshold = 5.;
    const double scaledFiedlerComponentThreshold = 0.2;

    AssemblyGraph& newAssemblyGraph = *this;
    const bool debug = true;

    // Initialize verticesBySegment.
    verticesBySegment.clear();
    verticesBySegment.resize(packedMarkerGraph.segments.size());

    // Work vectors used below.
    vector<vertex_descriptor> descendants;
    vector<vertex_descriptor> ancestors;
    vector<vertex_descriptor> neighbors;

    // Loop over vertices of the old AssemblyGraph.
    uint64_t splitCount = 0;
    BGL_FORALL_VERTICES(vOld, oldAssemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& oldVertex = oldAssemblyGraph[vOld];
        const uint64_t segmentId = oldVertex.segmentId;
        if(debug) {
            cout << "Local clustering on " << oldVertex.stringId() << endl;
        }

        // To create the neighborhood of this vertex on which we will do
        // clustering, we need the descendants and ancestors.
        oldAssemblyGraph.findDescendants(vOld, distance, descendants);
        oldAssemblyGraph.findAncestors(vOld, distance, ancestors);

        // Combine descendants and ancestors to create the neighborhood.
        neighbors.clear();
        set_union(
            descendants.begin(), descendants.end(),
            ancestors.begin(), ancestors.end(),
            back_inserter(neighbors));

        // For each of the journey entries in oldVertex, follow the journey
        // until we exit it. Keep track of the neighborhood vertices we encounter.
        // Matrix m is indexed by [journeyEntryIndex][neighborIndex].
        vector< vector<bool> > m(
            oldVertex.journeyEntries.size(),
            vector<bool>(neighbors.size(), false));
        for(uint64_t journeyEntryIndex=0; journeyEntryIndex<oldVertex.journeyEntries.size();
            journeyEntryIndex++) {
            const JourneyEntry& journeyEntry = oldVertex.journeyEntries[journeyEntryIndex];
            const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
            const vector<vertex_descriptor>& journey = oldAssemblyGraph.journeys[orientedReadId.getValue()];

            // Follow the journey forward until we leave the neighborhood.
            for(uint64_t position=journeyEntry.position+1; position<journey.size(); position++) {
                const vertex_descriptor u = journey[position];

                // Look it up in the neighborhood.
                auto it = lower_bound(neighbors.begin(), neighbors.end(), u);
                if(it == neighbors.end() or *it != u) {
                    // Not found.
                    break;
                }

                // We found it. Update the matrix.
                const uint64_t neighborIndex = it - neighbors.begin();
                m[journeyEntryIndex][neighborIndex] = true;
            }

            // Follow the journey backward until we leave the neighborhood.
            for(int64_t position=int64_t(journeyEntry.position)-1L; position>=0; position--) {
                const vertex_descriptor u = journey[position];

                // Look it up in the neighborhood.
                auto it = lower_bound(neighbors.begin(), neighbors.end(), u);
                if(it == neighbors.end() or *it != u) {
                    // Not found.
                    break;
                }

                // We found it. Update the matrix.
                const uint64_t neighborIndex = it - neighbors.begin();
                m[journeyEntryIndex][neighborIndex] = true;
            }
        }



        // Compute a singular value decomposition of this matrix.
        const int M = int(m.size());
        const int N = int(m.front().size());
        vector<double> A(M*N);
        for(uint64_t journeyEntryIndex=0; journeyEntryIndex<m.size(); journeyEntryIndex++) {
            for(uint64_t neighborIndex=0; neighborIndex<neighbors.size(); neighborIndex++) {
                A[journeyEntryIndex + M * neighborIndex] = double(int(m[journeyEntryIndex][neighborIndex]));
            }
        }
        const string JOBU = "A";
        const string JOBVT = "A";
        const int LDA = M;
        vector<double> S(min(M, N));
        vector<double> U(M*M);
        const int LDU = M;
        vector<double> VT(N*N);
        const int LDVT = N;
        const int LWORK = 10 * max(M, N);
        vector<double> WORK(LWORK);
        int INFO = 0;
        dgesvd_(
            JOBU.data(), JOBVT.data(),
            M, N,
            &A[0], LDA, &S[0], &U[0], LDU, &VT[0], LDVT, &WORK[0], LWORK, INFO);
        if(INFO != 0) {
            throw runtime_error("Error " + to_string(INFO) +
                " computing SVD decomposition in createByLocalClustering.");
        }

        // See how many large singular values we have.
        uint64_t nLarge = 0;
        for(; nLarge<S.size(); nLarge++) {
            SHASTA_ASSERT(S[nLarge] >= 0.);
            if(S[nLarge] < singularValueThreshold) {
                break;
            }
        }

        if(debug) {
            cout << "Singular values:";
            for(const double s: S) {
                cout << " " << s;
            }
            cout << "\n";
            cout << "Number of large singular values is " << nLarge << "\n";
            cout << "OrientedReadId,Position,";
            for(uint64_t i=0; i<nLarge; i++) {
                cout << "U" << i << ",";
            }
            for(const vertex_descriptor u: neighbors) {
                cout << oldAssemblyGraph.vertexStringId(u) << ",";
            }
            cout << "\n";
            for(uint64_t journeyEntryIndex=0; journeyEntryIndex<m.size(); journeyEntryIndex++) {
                const JourneyEntry& journeyEntry = oldVertex.journeyEntries[journeyEntryIndex];
                const vector<bool>& row = m[journeyEntryIndex];
                cout << journeyEntry.orientedReadId << ",";
                cout << journeyEntry.position << ",";
                for(uint64_t i=0; i<nLarge; i++) {
                    cout << U[journeyEntryIndex + i * m.size()] << ",";
                }
                for(uint64_t neighborIndex=0; neighborIndex<neighbors.size(); neighborIndex++) {
                    cout << int(row[neighborIndex]) << ",";
                }
                cout << endl;
            }

        }

        // If the number of large singular values is 1, don't split this vertex.
        // Create a vertex in the new assembly graph, with the same journey entries.
        if(nLarge == 1) {
            if(debug) {
                cout << "This vertex will not be split.\n";
            }
            const vertex_descriptor vNew = boost::add_vertex(
                AssemblyGraphVertex(segmentId, verticesBySegment[segmentId].size()),
                newAssemblyGraph);
            verticesBySegment[segmentId].push_back(vNew);
            AssemblyGraphVertex& newVertex = newAssemblyGraph[vNew];
            newVertex.journeyEntries = oldVertex.journeyEntries;
            continue;
        }

        // If getting here, the number of large singular values is greater than one,
        // and we will split this old vertex into two.
        // Each journey entry gets assigned to a new vertex based on the sign
        // of the Fiedler vector component for that entry.
        // The Fiedler vector is the left eigenvector corresponding to the second
        // highest singular value.
        // For the case of two clusters, the sign of its components
        // determines the cluster each journey entry belongs to.
        ++splitCount;
        const double* fiedlerVector = U.data() + m.size();


        // Create the two new vertices.
        const vertex_descriptor v0New = boost::add_vertex(
            AssemblyGraphVertex(segmentId, verticesBySegment[segmentId].size()),
            newAssemblyGraph);
        verticesBySegment[segmentId].push_back(v0New);
        AssemblyGraphVertex& newVertex0 = newAssemblyGraph[v0New];

        const vertex_descriptor v1New = boost::add_vertex(
            AssemblyGraphVertex(segmentId, verticesBySegment[segmentId].size()),
            newAssemblyGraph);
        verticesBySegment[segmentId].push_back(v1New);
        AssemblyGraphVertex& newVertex1 = newAssemblyGraph[v1New];

        // Assign journey entries to the new vertices.
        const double fiedlerComponentThreshold = scaledFiedlerComponentThreshold / sqrt(double(m.size()));
        uint64_t n0 = 0;
        uint64_t n1 = 0;
        uint64_t nOther = 0;
        for(uint64_t journeyEntryIndex=0; journeyEntryIndex<m.size(); journeyEntryIndex++) {
            const JourneyEntry& journeyEntry = oldVertex.journeyEntries[journeyEntryIndex];
            const double fiedlerComponent = fiedlerVector[journeyEntryIndex];

            if(fiedlerComponent > fiedlerComponentThreshold) {
                newVertex0.journeyEntries.push_back(journeyEntry);
                ++n0;
            } else if(fiedlerComponent < -fiedlerComponentThreshold) {
                newVertex1.journeyEntries.push_back(journeyEntry);
                ++n1;
            } else {
                ++nOther;
            }
        }
        if(debug) {
            cout << "Vertex split: n0 " << n0 << ", n1 " << n1 <<
                ", other " << nOther <<
                ", total " << m.size() << "\n";
        }

    }



    // Initialize  oriented reads journeys.
    journeys.clear();
    journeys.resize(packedMarkerGraph.journeys.size());
    for(uint64_t i=0; i<journeys.size(); i++) {
        const auto packedMarkerGraphJourney = packedMarkerGraph.journeys[i];
        journeys[i].resize(packedMarkerGraphJourney.size(), null_vertex());
    }


    // Construct oriented read journeys from the journey entries
    // stored in the vertices we created.
    BGL_FORALL_VERTICES(v, newAssemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = newAssemblyGraph[v];

        for(const JourneyEntry& journeyEntry: vertex.journeyEntries) {
            const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
            const uint64_t position = journeyEntry.position;

            vector<vertex_descriptor>& journey = journeys[orientedReadId.getValue()];
            journey[position] = v;
        }
    }

    if(debug) {
        cout << splitCount << " vertices of the old assembly graph were split, "
            "out of " << num_vertices(oldAssemblyGraph) << " total." << endl;
        cout << "The new assembly graph has " << num_vertices(newAssemblyGraph) <<
            " vertices." << endl;
    }

    // Now that we have the vertices and the journeys we can create the links.
    createLinks();
}



void AssemblyGraph::computeJaccardGraph(
    uint64_t threadCount,
    double minJaccard,
    uint64_t minComponentSize,
    uint64_t minPathLength)
{
    AssemblyGraph& assemblyGraph = *this;

    // Find candidate pairs.
    findJaccardGraphCandidatePairs(threadCount);

    // To compute Jaccard similarity, we need to get distinct oriented reads
    // for each vertex.
    computeVertexOrientedReadIds(threadCount);

    // Now we can compute Jaccard similarity for all the candidate pairs.
    computeJaccardPairs(threadCount, minJaccard);
    computeJaccardGraphData.candidatePairs.clear();
    computeJaccardGraphData.candidatePairs.shrink_to_fit();
    clearVertexOrientedReadIds();

    // Create the JaccardGraph.
    JaccardGraph jaccardGraph(*this);
    BGL_FORALL_VERTICES(av, assemblyGraph, AssemblyGraph) {
        jaccardGraph.addVertex(av);
    }
    for(const auto& p: computeJaccardGraphData.goodPairs) {
        const double jaccard = p.second;
        jaccardGraph.addEdge(p.first.first,  p.first.second, jaccard);
    }
    cout << "The initial Jaccard graph has " << num_vertices(jaccardGraph) <<
        " vertices and " << num_edges(jaccardGraph) << " edges." << endl;

    // Cleanup.
    computeJaccardGraphData.goodPairs.clear();
    computeJaccardGraphData.goodPairs.shrink_to_fit();

    // Remove Jaccard graph vertices belonging to strongly connected
    // component. After this, the Jaccard graph is acyclic.
    jaccardGraph.removeStronglyConnectedComponents();

    // Compute connected components so we can process them one at a time.
    vector< shared_ptr<JaccardGraph> > components;
    jaccardGraph.computeConnectedComponents(minComponentSize, components);



    // Process each connected component separately.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        JaccardGraph& componentGraph = *components[componentId];
        cout << "Processing Jaccard graph connected component " << componentId <<
            " with " << num_vertices(*components[componentId]) << " vertices." << endl;

        if(not componentGraph.findLongPaths(minPathLength)) {
            cout << "Jaccard graph component " << componentId << " is cyclic." << endl;
        }

        // Write it out.
        componentGraph.writeGraphviz(
            "JaccardGraphComponent-" + to_string(componentId) + ".dot", minJaccard);
    }


}



void AssemblyGraph::findJaccardGraphCandidatePairs(uint64_t threadCount)
{
    computeJaccardGraphData.threadCandidatePairs.clear();
    computeJaccardGraphData.threadCandidatePairs.resize(threadCount);
    setupLoadBalancing(journeys.size(), 1);
    runThreads(&AssemblyGraph::findJaccardGraphCandidatePairsThreadFunction, threadCount);

    // Consolidate and deduplicate the candidate pairs found by all threads.
    computeJaccardGraphData.candidatePairs.clear();
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        const auto& v = computeJaccardGraphData.threadCandidatePairs[threadId];
        copy(v.begin(), v.end(), back_inserter(computeJaccardGraphData.candidatePairs));
    }
    computeJaccardGraphData.threadCandidatePairs.clear();
    deduplicate(computeJaccardGraphData.candidatePairs);
    cout << "Found " << computeJaccardGraphData.candidatePairs.size() <<
        " candidate pairs for the Jaccard graph." << endl;
}



void AssemblyGraph::findJaccardGraphCandidatePairsThreadFunction(uint64_t threadId)
{

    auto& threadCandidatePairs = computeJaccardGraphData.threadCandidatePairs[threadId];

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(uint64_t orientedReadId=begin; orientedReadId!=end; orientedReadId++) {
            const vector<vertex_descriptor>& journey = journeys[orientedReadId];

            for(uint64_t i=1; i<journey.size(); i++) {
                if(journey[i] == null_vertex()) {
                    continue;
                }
                for(uint64_t j=0; j<i; j++) {
                    if(journey[j] == null_vertex()) {
                        continue;
                    }
                    auto p = make_pair(journey[j], journey[i]);
                    if(p.first != p.second) {
                        threadCandidatePairs.push_back(p);
                    }
                }
            }
        }
    }
}



void AssemblyGraph::computeVertexOrientedReadIds(uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;

    // Find distinct oriented read ids in each vertex.
    computeVertexOrientedReadIdsData.allVertices.clear();
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        computeVertexOrientedReadIdsData.allVertices.push_back(v);
    }
    setupLoadBalancing(computeVertexOrientedReadIdsData.allVertices.size(), 100);
    runThreads(&AssemblyGraph::computeVertexOrientedReadIdsThreadFunction, threadCount);
    computeVertexOrientedReadIdsData.allVertices.clear();
    computeVertexOrientedReadIdsData.allVertices.shrink_to_fit();
}



void AssemblyGraph::computeVertexOrientedReadIdsThreadFunction(uint64_t threadId)
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const vertex_descriptor v = computeVertexOrientedReadIdsData.allVertices[i];
            assemblyGraph[v].computeOrientedReadIds();
        }
    }

}



// Compute Jaccard similarity between two vertices.
// This requires vertex oriented read ids to be available.
double AssemblyGraph::computeJaccard(
    vertex_descriptor v0,
    vertex_descriptor v1,
    vector<OrientedReadId>& commonOrientedReadIds) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];;

    commonOrientedReadIds.clear();
    set_intersection(
        vertex0.orientedReadIds.begin(), vertex0.orientedReadIds.end(),
        vertex1.orientedReadIds.begin(), vertex1.orientedReadIds.end(),
        back_inserter(commonOrientedReadIds));

    const uint64_t intersectionSize = commonOrientedReadIds.size();
    const uint64_t unionSize =
        vertex0.orientedReadIds.size() + vertex1.orientedReadIds.size() - intersectionSize;

    return double(intersectionSize) / double(unionSize);
}



void AssemblyGraphVertex::computeOrientedReadIds()
{
    orientedReadIds.clear();
    for(const JourneyEntry& journeyEntry: journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        if(orientedReadIds.empty() or orientedReadId!=orientedReadIds.back()) {
            orientedReadIds.push_back(orientedReadId);
        }
    }
}



void AssemblyGraph::clearVertexOrientedReadIds()
{
    AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        assemblyGraph[v].orientedReadIds.clear();
    }
}



void AssemblyGraph::computeJaccardPairs(uint64_t threadCount, double minJaccard)
{
    computeJaccardGraphData.minJaccard = minJaccard;
    computeJaccardGraphData.threadGoodPairs.clear();
    computeJaccardGraphData.threadGoodPairs.resize(threadCount);
    setupLoadBalancing(computeJaccardGraphData.candidatePairs.size(), 1000);
    runThreads(&AssemblyGraph::computeJaccardPairsThreadFunction, threadCount);

    // Consolidate the good pairs found by all threads.
    computeJaccardGraphData.goodPairs.clear();
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        const auto& v = computeJaccardGraphData.threadGoodPairs[threadId];
        copy(v.begin(), v.end(), back_inserter(computeJaccardGraphData.goodPairs));
    }
    computeJaccardGraphData.threadGoodPairs.clear();
    cout << "Found " << computeJaccardGraphData.goodPairs.size() <<
        " good pairs for the Jaccard graph." << endl;
}



void AssemblyGraph::computeJaccardPairsThreadFunction(uint64_t threadId)
{
    const double minJaccard = computeJaccardGraphData.minJaccard;
    vector<OrientedReadId> commonOrientedReadIds;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate pairs assigned to this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const auto& p = computeJaccardGraphData.candidatePairs[i];
            const double jaccard = computeJaccard(p.first, p.second, commonOrientedReadIds);
            if(jaccard >= minJaccard) {
                computeJaccardGraphData.threadGoodPairs[threadId].push_back(make_pair(p, jaccard));
            }
        }
    }
}



// Compute the MSA for a link using spoa.
// Takes as input a vector of (sequence, frequency) containing
// the MSA sequences for the oriented reads of the link
// and the number of times each was found.
void shasta::mode3a::linkMsaUsingSpoa(
    const vector< pair<vector<Base>, uint64_t> >& msaSequences,
    uint64_t maxLength,
    ostream& html,
    vector<Base>& consensusSequence
    )
{
    // We want to enter the msaSequences in order of decreasing frequency.
    // Create a table of pairs (msaSequenceIndex, frequency)
    // where msaSequenceIndex is the index in the msaSequences vector.
    // Then sort by decreasing frequency.
    vector< pair<uint64_t, uint64_t> > msaSequencesTable;
    for(uint64_t msaSequenceIndex=0; msaSequenceIndex<msaSequences.size(); msaSequenceIndex++) {
        const auto& p = msaSequences[msaSequenceIndex];
        const uint64_t frequency = p.second;
        msaSequencesTable.push_back(make_pair(msaSequenceIndex, frequency));
    }
    sort(msaSequencesTable.begin(), msaSequencesTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());


    // If any of the sequences are too long, just return
    // the most frequent sequence as the consensus.
    bool isTooLong = false;
    for(const auto& p: msaSequences) {
        const vector<Base>& msaSequence = p.first;
        if(msaSequence.size() > maxLength) {
            isTooLong = true;
            break;
        }
    }
    if(isTooLong) {
        const uint64_t mostFrequentSequenceIndex = msaSequencesTable.front().first;
        consensusSequence = msaSequences[mostFrequentSequenceIndex].first;
        return;
    }


    // Create the spoa alignment engine and alignment graph.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto spoaAlignmentEngine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);
    spoa::Graph spoaAlignmentGraph;



    // Add the sequences to the MSA in order of decreasing frequency.
    // Give each sequence a weight equal to its frequency.
    if(html) {
        html <<
            "<h3>MSA input</h3>"
            "Oriented read MSA sequences are used in the following order:"
            "<table><tr>"
            "<th>Index<br>by<br>frequency"
            "<th>MSA<br>sequence<br>index"
            "<th>MSA<br>sequence<br>frequency"
            "<th>MSA<br>sequence"
            "<th>MSA<br>sequence<br>length";
    }
    string sequenceString;
    for(uint64_t indexByFrequency=0; indexByFrequency<msaSequencesTable.size(); indexByFrequency++) {
        const auto& p = msaSequencesTable[indexByFrequency];
        const uint64_t msaSequenceIndex = p.first;
        const uint64_t frequency = p.second;
        const auto& q = msaSequences[msaSequenceIndex];
        SHASTA_ASSERT(q.second == frequency);
        const vector<Base>& msaSequence = q.first;

        if(html) {
            html << "<tr>"
                "<td class=centered>" << indexByFrequency <<
                "<td class=centered>" << msaSequenceIndex <<
                "<td class=centered>" << frequency <<
                "<td class=centered style='font-family:courier'>";
            copy(msaSequence.begin(), msaSequence.end(), ostream_iterator<Base>(html));
            html << "<td class=centered>" << msaSequence.size();
        }

        sequenceString.clear();
        for(const Base base: msaSequence) {
            sequenceString += base.character();
        }
        auto alignment = spoaAlignmentEngine->Align(sequenceString, spoaAlignmentGraph);
        spoaAlignmentGraph.AddAlignment(alignment, sequenceString, uint32_t(frequency));
    }
    if(html) {
        html << "</table>";
    }



#if 0
    // Get the consensus as computed by spoa.
    // This is generally not as accurate as coverage-based consensus computed below.
    // See here for some discussion:
    // https://github.com/rvaser/spoa/issues/63
    const string spoaConsensusString = spoaAlignmentGraph.GenerateConsensus();
    vector<Base> spoaConsensus;
    for(const char c: spoaConsensusString) {
        spoaConsensus.push_back(Base::fromCharacter(c));
    }
#endif



    // Get the MSA alignment.
    // The true argument causes a final alignment entry equal to the consensus.
    vector<string> alignment = spoaAlignmentGraph.GenerateMultipleSequenceAlignment(false);
    SHASTA_ASSERT(alignment.size() == msaSequencesTable.size());

    // Compute coverage at each alignment position for each of the 5 AlignedBases.
    const uint64_t alignmentLength = alignment.front().size();
    vector< array<uint64_t, 5> > coverage(alignmentLength, {0, 0, 0, 0, 0});
    for(uint64_t indexByFrequency=0; indexByFrequency<msaSequencesTable.size(); indexByFrequency++) {
        const string& alignmentRow = alignment[indexByFrequency];
        SHASTA_ASSERT(alignmentRow.size() == alignmentLength);
        for(uint64_t position=0; position<alignmentLength; position++) {
            const AlignedBase b = AlignedBase::fromCharacter(alignmentRow[position]);
            coverage[position][b.value] += msaSequencesTable[indexByFrequency].second;
        }
    }

    // Compute coverage-based consensus at each alignment position.
    vector<AlignedBase> coverageBasedAlignedConsensus;
    for(const auto& c: coverage) {
        const uint64_t iBase = std::max_element(c.begin(), c.end()) - c.begin();
        coverageBasedAlignedConsensus.push_back(AlignedBase::fromInteger(iBase));
    }
    SHASTA_ASSERT(coverageBasedAlignedConsensus.size() == alignmentLength);

    // Take out the gaps.
    vector<Base> coverageBasedConsensus;
    for(const AlignedBase b: coverageBasedAlignedConsensus) {
        if(not b.isGap()) {
            coverageBasedConsensus.push_back(Base(b));
        }
    }
    consensusSequence = coverageBasedConsensus;

    // Compute concordant and discordant count at each alignment position.
    vector<uint64_t> concordantCount(alignmentLength, 0);
    vector<uint64_t> discordantCount(alignmentLength, 0);
    for(uint64_t position=0; position<alignmentLength; position++) {
        for(uint64_t iBase=0; iBase<5; iBase++) {
            if(AlignedBase::fromInteger(iBase) == coverageBasedAlignedConsensus[position]) {
                concordantCount[position] += coverage[position][iBase];
            } else {
                discordantCount[position] += coverage[position][iBase];
            }
        }
    }



    if(html) {
        html <<
            "<h2>Consensus sequence</h2>"
            "The consensus sequence has length " << coverageBasedConsensus.size() <<
            ":<div style='font-family:courier'>";
        for(const Base base: coverageBasedConsensus) {
            html << base;
        }
        html << "</div>";

        // See if the consensus is equal to one of our MSA sequences.
        // This is often the case but does not have to be.
        bool found = false;
        for(uint64_t msaSequenceIndex=0; msaSequenceIndex<msaSequences.size(); msaSequenceIndex++) {
            const auto& p = msaSequences[msaSequenceIndex];
            if(coverageBasedConsensus == p.first) {
                html << "The consensus sequences is the same as the MSA sequence with index " <<
                    msaSequenceIndex << ".";
                found = true;
                break;
            }
        }
        if(not found) {
            html << "<br>The consensus sequence is not equal to any of the MSA sequences.";
        }

        // Write MSA details.
        html << "<h3>MSA details</h3>"
            "<table><tr>"
            "<th>Index<br>by<br>frequency"
            "<th>MSA<br>sequence<br>index"
            "<th>MSA<br>sequence<br>frequency"
            "<th>MSA<br>sequence"
            "<th>MSA<br>sequence<br>length";
        for(uint64_t indexByFrequency=0; indexByFrequency<msaSequencesTable.size(); indexByFrequency++) {
            const auto& p = msaSequencesTable[indexByFrequency];
            const uint64_t msaSequenceIndex = p.first;
            const uint64_t frequency = p.second;
            const auto& q = msaSequences[msaSequenceIndex];
            SHASTA_ASSERT(q.second == frequency);
            const vector<Base>& msaSequence = q.first;

            html << "<tr>"
                "<td class=centered>" << indexByFrequency <<
                "<td class=centered>" << msaSequenceIndex <<
                "<td class=centered>" << frequency <<
                "<td class=centered style='font-family:courier'>";
            for(uint64_t position=0; position<alignment[indexByFrequency].size(); position++) {
                const char c = alignment[indexByFrequency][position];
                const bool isConcordant = (c == coverageBasedAlignedConsensus[position].character());
                if(not isConcordant) {
                    html << "<span style='color:Red;font-weight:bold'>";
                }
                html << c;
                if(not isConcordant) {
                    html << "</span>";
                }
            }
            html << "<td class=centered>" << msaSequence.size();
        }

        // Add rows with coverage for each of the 5 AlignedBases.
        for(uint64_t iBase=0; iBase<5; iBase++) {
            if(iBase==0) {
                html << "<tr><td colspan=2 rowspan=5 class=centered>Coverage";
            } else {
                html << "<tr>";
            }
            const AlignedBase base = AlignedBase::fromInteger(iBase);
            html << "<td class=centered style='font-family:courier'>" << base <<
                "<td class=centered style='font-family:courier'>";
            for(uint64_t position=0; position<alignmentLength; position++) {
                const uint64_t n = coverage[position][iBase];
                if(n == 0) {
                    html << "&nbsp;";
                } else if(n < 10) {
                    html << n;
                } else {
                    html << "*";
                }
            }
        }

        // Add a row with consensus.
        html << "<tr><td class=centered colspan=3>Consensus"
            "<td class=centered style='font-family:courier'>";
        for(uint64_t position=0; position<coverageBasedAlignedConsensus.size(); position++) {
            const AlignedBase b = coverageBasedAlignedConsensus[position];
            html << b;
        }
        html <<
            "<td class=centered>" << coverageBasedConsensus.size();


        // Add a row with concordant count.
        html <<
            "<tr><td class=centered colspan=3>Concordant"
            "<td class=centered style='font-family:courier'>";
        for(uint64_t i=0; i<concordantCount.size(); i++) {
            const uint64_t n = concordantCount[i];
            if(n == 0) {
                html << "&nbsp;";
            } else if(n < 10) {
                html << n;
            } else {
                html << "*";
            }
        }
        html << "<td>";

        // Add a row with discordant count.
        html <<
            "<tr><td class=centered colspan=3>Discordant"
            "<td class=centered style='font-family:courier'>";
        for(uint64_t i=0; i<discordantCount.size(); i++) {
            const uint64_t n = discordantCount[i];
            if(n == 0) {
                html << "&nbsp;";
            } else if(n < 10) {
                html << n;
            } else {
                html << "*";
            }
        }
        html << "<td>";

        html << "</table>";

    }
}


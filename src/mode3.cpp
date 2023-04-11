
// Shasta
#include "mode3.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "deduplicate.hpp"
#include "findMarkerId.hpp"
#include "html.hpp"
#include "MarkerGraph.hpp"
#include "mode3-AssemblyPath.hpp"
#include "mode3-JaccardGraph.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "ReadFlags.hpp"
#include "mode3-SegmentPairInformation.hpp"
#include "SubsetGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
// Include disjoint_sets.hpp first to avoid Boost problems.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include <bitset>
#include "fstream.hpp"
#include <map>
#include <queue>
#include <set>
#include <unordered_set>

#include "MultithreadedObject.tpp"
template class MultithreadedObject<mode3::AssemblyGraph>;


// Each  linear chain of marker graph edges generates a segment.
void AssemblyGraph::createSegmentPaths()
{
    const bool debug = false;

    createNew(markerGraphPaths, "Mode3-MarkerGraphPaths");
    const MarkerGraph::EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    using MarkerGraphPath = vector<MarkerGraph::EdgeId>;
    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;

    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(MarkerGraph::EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        if(debug) {
            cout << "Starting a new path at edge " << startEdgeId << endl;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v1 = edge.target;
            const auto outEdges = markerGraph.edgesBySource[v1];
            if(outEdges.size() != 1) {
                break;
            }
            const auto inEdges = markerGraph.edgesByTarget[v1];
            if(inEdges.size() != 1) {
                break;
            }
            edgeId = outEdges[0];
            if(edgeId == startEdgeId) {
                isCircular = true;
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
            if(debug) {
                cout << "Moving forward: added " << edgeId << endl;
            }
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v0 = edge.source;
                const auto outEdges = markerGraph.edgesBySource[v0];
                if(outEdges.size() != 1) {
                    break;
                }
                const auto inEdges = markerGraph.edgesByTarget[v0];
                if(inEdges.size() != 1) {
                    break;
                }
                edgeId = inEdges[0];
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
                if(debug) {
                    cout << "Moving backward: added " << edgeId << endl;
                }
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }

        // Store this path as a new segment.
        markerGraphPaths.appendVector();
        for(const MarkerGraphEdgeId edgeId: path) {
            markerGraphPaths.append(edgeId);
        }
    }



    // Check that all edges of the marker graph were found.
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());


    // Debug output: write the paths.
    if(debug) {
        ofstream csv("Paths.csv");
        for(uint64_t segmentId=0; segmentId<markerGraphPaths.size(); segmentId++) {
            const auto path = markerGraphPaths[segmentId];
            for(const MarkerGraphEdgeId edgeId: path) {
                csv << segmentId << ",";
                csv << edgeId << "\n";
            }
        }
    }

}



// Compute coverage for all segments.
// It is computed as average marker graph edge coverage
// over the marker graph edges in the path of each segment.
void AssemblyGraph::computeSegmentCoverage()
{
    // Initialize segmentCoverage.
    createNew(segmentCoverage, "Mode3-SegmentCoverage");
    const uint64_t segmentCount = markerGraphPaths.size();
    segmentCoverage.resize(segmentCount);

    // Loop over all segments.
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {

        // Access the marker graph path for this segment.
        const span<MarkerGraphEdgeId> path = markerGraphPaths[segmentId];


        // Loop over this path.
        uint64_t coverageSum = 0.;
        for(uint64_t position=0; position<path.size(); position++) {
            MarkerGraphEdgeId& edgeId = path[position];

            // Add the marker intervals on this marker graph edge.
            const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            coverageSum += markerIntervals.size();
        }

        segmentCoverage[segmentId] = float(coverageSum) / float(path.size());

    }


    // Write a histogram of segment coverage.
    vector<uint64_t> histogram;
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const uint64_t coverage = uint64_t(std::round(segmentCoverage[segmentId]));
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }
    ofstream csv("SegmentCoverageHistogram.csv");
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



void AssemblyGraph::computeMarkerGraphEdgeTable(size_t threadCount)
{

    // Initialize the marker graph edge table.
    createNew(markerGraphEdgeTable, "Mode3-MarkerGraphEdgeTable");
    markerGraphEdgeTable.resize(markerGraph.edges.size());
    fill(markerGraphEdgeTable.begin(), markerGraphEdgeTable.end(), make_pair(
        std::numeric_limits<uint64_t>::max(),
        std::numeric_limits<uint32_t>::max()
        ));

    // Fill in the marker graph edge table.
    const uint64_t batchSize = 100;
    setupLoadBalancing(markerGraphPaths.size(), batchSize);
    runThreads(&AssemblyGraph::computeMarkerGraphEdgeTableThreadFunction, threadCount);
}



void AssemblyGraph::computeMarkerGraphEdgeTableThreadFunction(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {
            const span<MarkerGraphEdgeId> path = markerGraphPaths[segmentId];

            // Loop over the path of this segment.
            for(uint64_t position=0; position<path.size(); position++) {
                const MarkerGraphEdgeId edgeId = path[position];

                // Store the marker graph edge table entry for this edge.
                SHASTA_ASSERT(edgeId < markerGraphEdgeTable.size());
                markerGraphEdgeTable[edgeId] = make_pair(segmentId, position);
            }
        }

    }
}



void AssemblyGraph::computeMarkerGraphJourneys(size_t threadCount)
{
    const bool debug = true;

    createNew(markerGraphJourneys, "tmp-mode3-MarkerGraphJourneys");

    uint64_t batchSize = 1000;
    markerGraphJourneys.beginPass1(markers.size());
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&AssemblyGraph::computeMarkerGraphJourneysPass1, threadCount);
    markerGraphJourneys.beginPass2();
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&AssemblyGraph::computeMarkerGraphJourneysPass2, threadCount);
    markerGraphJourneys.endPass2();

    batchSize = 100;
    setupLoadBalancing(markerGraphJourneys.size(), batchSize);
    runThreads(&AssemblyGraph::sortMarkerGraphJourneys, threadCount);

    if(debug) {
        ofstream csv("MarkerGraphJourneys.csv");
        csv << "OrientedReadId,SegmentId,Position,ordinal0,Ordinal1\n";
        for(uint64_t i=0; i<markers.size(); i++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
            const auto markerGraphJourney = markerGraphJourneys[i];
            for(uint64_t position=0; position<markerGraphJourney.size(); position++) {
                const MarkerGraphJourneyEntry& entry = markerGraphJourney[position];
                csv << orientedReadId << ",";
                csv << entry.segmentId << ",";
                csv << entry.position << ",";
                csv << entry.ordinals[0] << ",";
                csv << entry.ordinals[1] << "\n";
            }
        }

    }
}



void AssemblyGraph::computeMarkerGraphJourneysPass1(size_t threadId)
{
    computeMarkerGraphJourneysPass12(1);
}



void AssemblyGraph::computeMarkerGraphJourneysPass2(size_t threadId)
{
    computeMarkerGraphJourneysPass12(2);
}



void AssemblyGraph::computeMarkerGraphJourneysPass12(uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(MarkerGraph::EdgeId edgeId=begin; edgeId!=end; ++edgeId) {
            const auto& p = markerGraphEdgeTable[edgeId];
            const uint64_t segmentId = p.first;
            const uint32_t position = p.second;
            SHASTA_ASSERT(segmentId != std::numeric_limits<uint64_t>::max());
            SHASTA_ASSERT(position != std::numeric_limits<uint32_t>::max());

            // Loop over the marker intervals of this marker graph edge..
            const auto markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                const OrientedReadId orientedReadId = markerInterval.orientedReadId;

                if(pass == 1) {
                    markerGraphJourneys.incrementCountMultithreaded(orientedReadId.getValue());
                } else {
                    MarkerGraphJourneyEntry markerGraphJourneyEntry;
                    markerGraphJourneyEntry.segmentId = segmentId;
                    markerGraphJourneyEntry.position = position;
                    markerGraphJourneyEntry.ordinals = markerInterval.ordinals;
                    markerGraphJourneys.storeMultithreaded(orientedReadId.getValue(), markerGraphJourneyEntry);
                }
            }
        }
    }
}



void AssemblyGraph::sortMarkerGraphJourneys(size_t threadId)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            auto markerGraphJourney = markerGraphJourneys[i];
            sort(markerGraphJourney.begin(), markerGraphJourney.end());
        }
    }
}


// The assembly graph journey of an oriented read
// is the sequence of segmentIds it encounters.
void AssemblyGraph::computeAssemblyGraphJourneys()
{
    const bool debug = true;

    // Initialize the assembly graph journeys.
    createNew(assemblyGraphJourneys, "Mode3-AssemblyGraphJourneys");

    // Work vector defined outside the loop to reduce memory allocation overhead.
    vector<AssemblyGraphJourneyEntry> assemblyGraphJourney;

    // Loop over all oriented reads.
    for(uint64_t i=0; i<markerGraphJourneys.size(); i++) {

        // Access the marker graph journey for this oriented read.
        const span<MarkerGraphJourneyEntry> markerGraphJourney = markerGraphJourneys[i];

        // Compute the assembly graph journey.
        computeAssemblyGraphJourney(markerGraphJourney, assemblyGraphJourney);

        // Store it.
        assemblyGraphJourneys.appendVector(assemblyGraphJourney);
    }



    // Write them out.
    if(debug) {
        ofstream csv("AssemblyGraphJourneys.csv");
        for(uint64_t i=0; i<assemblyGraphJourneys.size(); i++) {
            const ReadId readId = ReadId(i >> 1);
            const Strand strand = i & 1;
            const OrientedReadId orientedReadId(readId, strand);
            const span<AssemblyGraphJourneyEntry> assemblyGraphJourney = assemblyGraphJourneys[i];

            csv << orientedReadId << ",";
            for(const AssemblyGraphJourneyEntry entry: assemblyGraphJourney) {
                csv << entry.segmentId << ",";
            }
            csv << endl;
        }
    }



    // Write them out again, with more details.
    if(debug) {
        ofstream csv("AssemblyGraphJourneysDetails.csv");
        csv << "OrientedReadId,Position,SegmentId,"
            "First position,First ordinal0,First ordinal1,"
            "Last position,Last ordinal0,Last ordinal1\n";
        for(uint64_t i=0; i<assemblyGraphJourneys.size(); i++) {
            const ReadId readId = ReadId(i >> 1);
            const Strand strand = i & 1;
            const OrientedReadId orientedReadId(readId, strand);
            const span<AssemblyGraphJourneyEntry> assemblyGraphJourney = assemblyGraphJourneys[i];

            for(uint64_t position=0; position<assemblyGraphJourney.size(); position++) {
                const AssemblyGraphJourneyEntry& entry = assemblyGraphJourney[position];
                const MarkerGraphJourneyEntry& first = entry.markerGraphJourneyEntries[0];
                const MarkerGraphJourneyEntry& last = entry.markerGraphJourneyEntries[1];
                csv << orientedReadId << ",";
                csv << position << ",";
                csv << entry.segmentId << ",";
                csv << first.position << ",";
                csv << first.ordinals[0] << ",";
                csv << first.ordinals[1] << ",";
                csv << last.position << ",";
                csv << last.ordinals[0] << ",";
                csv << last.ordinals[1] << "\n";
            }
        }
    }


}



// Given the marker graph journey of an oriented read,
// find the corresponding assembly graph journey.
void AssemblyGraph::computeAssemblyGraphJourney(
    const span<MarkerGraphJourneyEntry> markerGraphJourney,
    vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney)
{
    // Start with an empty journey.
    assemblyGraphJourney.clear();

    // Loop over the marker graph journey, looking for places
    // where the segmentId changes.
    for(uint32_t i=0; i<markerGraphJourney.size(); /* Increment later */) {
        const MarkerGraphJourneyEntry& markerGraphJourneyEntry = markerGraphJourney[i];
        const uint64_t segmentId = markerGraphJourneyEntry.segmentId;

        // Move to the end of the streak with the same segmentId.
        const uint32_t streakBegin = i;
        uint32_t streakEnd = streakBegin + 1;
        for(;
            streakEnd<markerGraphJourney.size() and
            (markerGraphJourney[streakEnd].segmentId == segmentId);
            streakEnd++) {
        }

        // Store this segmentId in the assembly graph journey.
        AssemblyGraphJourneyEntry assemblyGraphJourneyEntry;
        assemblyGraphJourneyEntry.segmentId = segmentId;
        assemblyGraphJourneyEntry.markerGraphJourneyEntries[0] = markerGraphJourney[streakBegin];
        assemblyGraphJourneyEntry.markerGraphJourneyEntries[1] = markerGraphJourney[streakEnd - 1];
        assemblyGraphJourney.push_back(assemblyGraphJourneyEntry);

        // Prepare to handle the next segment.
        i = streakEnd;
    }
}



void AssemblyGraph::computeAssemblyGraphJourneyInfos()
{
    const bool debug = true;

    const uint64_t segmentCount = markerGraphPaths.size();
    const uint64_t readCount = assemblyGraphJourneys.size()/2;

    createNew(assemblyGraphJourneyInfos, "Mode3-AssemblyGraphJourneyInfos");

    // Pass 1.
    assemblyGraphJourneyInfos.beginPass1(segmentCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto assemblyGraphJourney = assemblyGraphJourneys[orientedReadId.getValue()];

            for(uint64_t position=0; position<assemblyGraphJourney.size(); position++) {
                const AssemblyGraphJourneyEntry& entry = assemblyGraphJourney[position];
                assemblyGraphJourneyInfos.incrementCount(entry.segmentId);
            }
        }
    }

    // Pass 2.
    assemblyGraphJourneyInfos.beginPass2();
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto assemblyGraphJourney = assemblyGraphJourneys[orientedReadId.getValue()];

            for(uint64_t position=0; position<assemblyGraphJourney.size(); position++) {
                const AssemblyGraphJourneyEntry& entry = assemblyGraphJourney[position];
                assemblyGraphJourneyInfos.store(entry.segmentId, make_pair(orientedReadId, position));
            }
        }
    }
    assemblyGraphJourneyInfos.endPass2();

    // Sort.
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const auto v = assemblyGraphJourneyInfos[segmentId];
        sort(v.begin(), v.end());
    }


    if(debug) {
        ofstream csv("SegmentJourneyInfo.csv");
        csv << "SegmentId,OrientedReadId,Position in assembly graph journey\n";
        for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
            const auto v = assemblyGraphJourneyInfos[segmentId];
            for(const auto& p: v) {
                csv << segmentId << ",";
                csv << p.first << ",";
                csv << p.second << "\n";
            }
        }
    }
}



// Find out if a segment contains a given OrientedReadId.
// This returns true if assemblyGraphJourneyInfos[segmentId]
// contains an entry with the given OrientedReadId.
bool AssemblyGraph::segmentContainsOrientedRead(
    uint64_t segmentId,
    OrientedReadId orientedReadId) const
{
    for(const auto& p: assemblyGraphJourneyInfos[segmentId]) {
        if(p.first == orientedReadId) {
            return true;
        }
    }
    return false;
}



void AssemblyGraph::findTransitions(std::map<SegmentPair, Transitions>& transitionMap)
{
    transitionMap.clear();

    for(ReadId readId=0; readId<assemblyGraphJourneys.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = assemblyGraphJourneys[orientedReadId.getValue()];

            for(uint64_t i=1; i<journey.size(); i++) {
                const auto& previous = journey[i-1];
                const auto& current = journey[i];
                SHASTA_ASSERT(previous.segmentId != current.segmentId);

                const SegmentPair segmentPair = make_pair(previous.segmentId, current.segmentId);
                transitionMap[segmentPair].push_back(
                    make_pair(orientedReadId, Transition({
                    previous.markerGraphJourneyEntries[1],
                    current.markerGraphJourneyEntries[0]})));

            }
        }
    }
}



void AssemblyGraph::createLinks(
    const std::map<SegmentPair, Transitions>& transitionMap,
    uint64_t minCoverage)
{
    createNew(links, "Mode3-Links");
    createNew(transitions, "Mode3-Transitions");
    for(const auto& p: transitionMap) {
        const auto& transitionVector = p.second;
        const uint64_t coverage = transitionVector.size();
        if(coverage >= minCoverage) {
            const uint64_t segmentId0 = p.first.first;
            const uint64_t segmentId1 = p.first.second;
            links.push_back(Link(segmentId0, segmentId1));
            transitions.appendVector(transitionVector);
        }
    }

    // Store link separation.
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        Link& link = links[linkId];

        // Check if these two segments are adjacent in the marker graph.
        const uint64_t segmentId0 = link.segmentId0;
        const uint64_t segmentId1 = link.segmentId1;
        const auto path0 = markerGraphPaths[segmentId0];
        const auto path1 = markerGraphPaths[segmentId1];
        const MarkerGraph::Edge lastEdge0 = markerGraph.edges[path0.back()];
        const MarkerGraph::Edge firstEdge1 = markerGraph.edges[path1.front()];
        if(lastEdge0.target == firstEdge1.source) {
            // The segments are adjacent. Set the link separation to 0.
            link.segmentsAreAdjacent = true;
            link.separation = 0;
        } else {
            // The segments are not adjacent.
            // Use the transitions to estimate the separation.
            const auto linkTransitions = transitions[linkId];
            const double separation = linkSeparation(linkTransitions, path0.size());

            link.segmentsAreAdjacent = false;
            link.separation = int32_t(std::round(separation));
        }
    }



    ofstream csv("Links.csv");
    csv << "LinkId,SegmentId0,SegmentId1,Coverage,Adjacent,Separation\n";
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        Link& link = links[linkId];

        csv << linkId << ",";
        csv << link.segmentId0 << ",";
        csv << link.segmentId1 << ",";
        csv << transitions[linkId].size() << ",";
        csv << (link.segmentsAreAdjacent ? "Yes" : "No") << ",";
        csv << link.separation << "\n";
    }

}



// Initial construction of the AssemblyGraph.
AssemblyGraph::AssemblyGraph(
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize,
    size_t threadCount,
    uint64_t readRepresentation,
    uint64_t k, // Marker length
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    const ConsensusCaller& consensusCaller) :
    MultithreadedObject<AssemblyGraph>(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    readRepresentation(readRepresentation),
    k(k),
    reads(reads),
    markers(markers),
    markerGraph(markerGraph),
    consensusCaller(consensusCaller)
{
    // K must be even.
    SHASTA_ASSERT((k % 2) == 0);

    // This assumes RLE is not used.
    SHASTA_ASSERT(reads.representation == 0);

    // Minimum number of transitions (oriented reads) to create a link.
    // If this equals 1, then the sequence of segments visited by every
    // oriented read is a path in the graph.
    // But that is not desirable because of the extra edges it causes.
    const uint64_t minLinkCoverage = 4; // EXPOSE WHEN CODE STABILIZES

    // Create a segment for each linear chain of marker graph edges.
    createSegmentPaths();
    computeSegmentCoverage();

    // Assembled sequence for each segment.
    assembleSegments();

    // Keep track of the segment and position each marker graph edge corresponds to.
    computeMarkerGraphEdgeTable(threadCount);

    // Compute marker graph and assembly graph journeys of all oriented reads.
    // We permanently store only the assembly graph journeys.
    computeMarkerGraphJourneys(threadCount);
    computeAssemblyGraphJourneys();
    markerGraphJourneys.remove();
    computeAssemblyGraphJourneyInfos();

    // Find transitions from segment to segment in the marker graph
    // journeys of all oriented reads, and store them keyed by the pair of segments.
    std::map<SegmentPair, Transitions> transitionMap;
    findTransitions(transitionMap);

    // Create a links between pairs of segments with a sufficient number of transitions.
    createLinks(transitionMap, minLinkCoverage);
    createConnectivity();
    flagBackSegments();

    cout << "The mode 3 assembly graph has " << markerGraphPaths.size() << " segments and " <<
        links.size() << " links." << endl;
}



string AssemblyGraph::largeDataName(const string& name) const
{
    if(largeDataFileNamePrefix.empty()) {
        return "";  // Anonymous;
    } else {
        return largeDataFileNamePrefix + name;
    }
}



// Constructor from binary data.
AssemblyGraph::AssemblyGraph(
    const string& largeDataFileNamePrefix,
    uint64_t readRepresentation,
    uint64_t k, // Marker length
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    const ConsensusCaller& consensusCaller) :
    MultithreadedObject<AssemblyGraph>(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    readRepresentation(readRepresentation),
    k(k),
    reads(reads),
    markers(markers),
    markerGraph(markerGraph),
    consensusCaller(consensusCaller)
{
    accessExistingReadOnly(markerGraphPaths, "Mode3-MarkerGraphPaths");
    accessExistingReadOnly(segmentCoverage, "Mode3-SegmentCoverage");
    accessExistingReadOnly(segmentSequences, "Mode3-SegmentSequences");
    accessExistingReadOnly(segmentVertexOffsets, "Mode3-SegmentVertexOffsets");
    accessExistingReadOnly(markerGraphEdgeTable, "Mode3-MarkerGraphEdgeTable");
    accessExistingReadOnly(assemblyGraphJourneys, "Mode3-AssemblyGraphJourneys");
    accessExistingReadOnly(assemblyGraphJourneyInfos, "Mode3-AssemblyGraphJourneyInfos");
    accessExistingReadOnly(links, "Mode3-Links");
    accessExistingReadOnly(transitions, "Mode3-Transitions");
    accessExistingReadOnly(linksBySource, "Mode3-LinksBySource");
    accessExistingReadOnly(linksByTarget, "Mode3-LinksByTarget");
    accessExistingReadOnly(isBackSegment, "Mode3-IsBackSegment");
    accessExistingReadOnly(clusterIds, "Mode3-ClusterIds");
}



void AssemblyGraph::createConnectivity()
{
    createNew(linksBySource, "Mode3-LinksBySource");
    createNew(linksByTarget, "Mode3-LinksByTarget");

    linksBySource.beginPass1(markerGraphPaths.size());
    linksByTarget.beginPass1(markerGraphPaths.size());
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        const Link& link = links[linkId];
        linksBySource.incrementCount(link.segmentId0);
        linksByTarget.incrementCount(link.segmentId1);
    }
    linksBySource.beginPass2();
    linksByTarget.beginPass2();
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        const Link& link = links[linkId];
        linksBySource.store(link.segmentId0, linkId);
        linksByTarget.store(link.segmentId1, linkId);
    }
    linksBySource.endPass2();
    linksByTarget.endPass2();
}



uint64_t AssemblyGraph::findLink(uint64_t segmentId0, uint64_t segmentId1) const
{
    for(const uint64_t linkId: linksBySource[segmentId0]) {
        if(links[linkId].segmentId1 == segmentId1) {
            return linkId;
        }
    }
    SHASTA_ASSERT(0);
}



// Flag back-segments.
// This does not do a full blown search for locally strongly connected components.
// A segment is marked as a back-segment if:
// - It has only a single incoming link.
// - It has a single outgoing link.
// - The incoming and outgoing links both connect to/from the same segment.
void AssemblyGraph::flagBackSegments()
{
    const uint64_t segmentCount = markerGraphPaths.size();
    createNew(isBackSegment, "Mode3-IsBackSegment");
    isBackSegment.resize(segmentCount);

    uint64_t backSegmentCount = 0;
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {

        // Initially flag it as not a back-segment.
        isBackSegment[segmentId] = false;

        // For a back-segment, there must be a single incoming link.
        const auto incomingLinks = linksByTarget[segmentId];
        if(incomingLinks.size() != 1) {
            continue;
        }

        // For a back-segment, there must be a single outgoing link.
        const auto outgoingLinks = linksBySource[segmentId];
        if(outgoingLinks.size() != 1) {
            continue;
        }

        // For a back-segment, the incoming and outgoing links
        // both connect to/from the same segment.
        const uint64_t incomingLinkId = incomingLinks[0];
        const uint64_t outgoingLinkId = outgoingLinks[0];
        const Link& incomingLink = links[incomingLinkId];
        const Link& outgoingLink = links[outgoingLinkId];
        if(incomingLink.segmentId0 != outgoingLink.segmentId1) {
            continue;
        }

        // Flag it as a back-segment.
        isBackSegment[segmentId] = true;
        ++backSegmentCount;
    }

    cout << "Found " << backSegmentCount << " back-segments." << endl;
}



// Get the children or parents of a given segment.
// Only use links with at least a specified coverage.
void AssemblyGraph::getChildrenOrParents(
    uint64_t segmentId,
    uint64_t direction, // 0=forward (children), 1=backward (parents).
    uint64_t minimumLinkCoverage,
    vector<uint64_t>& childrenOrParents) const
{
    switch(direction) {
    case 0:
        getChildren(segmentId, minimumLinkCoverage, childrenOrParents);
        break;
    case 1:
        getParents(segmentId, minimumLinkCoverage, childrenOrParents);
        break;
    default:
        SHASTA_ASSERT(0);
    }
}



void AssemblyGraph::getChildren(
    uint64_t segmentId,
    uint64_t minimumLinkCoverage,
    vector<uint64_t>& children) const
{
    children.clear();
    for(const auto linkId: linksBySource[segmentId]) {
        if(transitions.size(linkId) >= minimumLinkCoverage) {
            const Link& link = links[linkId];
            children.push_back(link.segmentId1);
        }
    }
}



void AssemblyGraph::getParents(
    uint64_t segmentId,
    uint64_t minimumLinkCoverage,
    vector<uint64_t>& parents) const
{
    parents.clear();
    for(const auto linkId: linksByTarget[segmentId]) {
        if(transitions.size(linkId) >= minimumLinkCoverage) {
            const Link& link = links[linkId];
            parents.push_back(link.segmentId0);
        }
    }
}



void AssemblyGraph::writeGfa(const string& baseName) const
{
    ofstream gfa(baseName + ".gfa");
    ofstream csv(baseName + ".csv");

    // Write the headers.
    gfa << "H\tVN:Z:1.0\n";
    csv << "Segment,Path Length,Sequence Length,Average coverage,Read count\n";

    // Write the segments.
    for(uint64_t segmentId=0; segmentId<markerGraphPaths.size(); segmentId++) {

        const auto sequence = segmentSequences[segmentId];
        gfa <<"S\t" << segmentId << "\t";
        copy(sequence.begin()+k/2, sequence.end()-k/2, ostream_iterator<Base>(gfa));
        gfa << "\n";

        const auto path = markerGraphPaths[segmentId];
        csv << segmentId << ",";
        csv << path.size() << ",";
        csv << sequence.size() << ",";
        csv << segmentCoverage[segmentId] << ",";
        csv << assemblyGraphJourneyInfos[segmentId].size() << "\n";
    }

    // Write the links.
    for(const Link& link: links) {
        if(true /*link.segmentsAreAdjacent*/) {
            gfa << "L\t" <<
                link.segmentId0 << "\t+\t" <<
                link.segmentId1 << "\t+\t0M\n";
        } else {
            // This writes non-adjacent links as Jumps (GFA 1.2).
            // The original Bandage does not display them.
            // BandageNG does, but they are not taken into account during graph
            // creation, so it i not useful to write them like this.
            // For this reason, the if condition above was set to true,
            // so this branch is never reached.
            // Leaveing the code in place for possible future use.
            gfa << "J\t" <<
                link.segmentId0 << "\t+\t" <<
                link.segmentId1 << "\t+\t" << k * link.separation << "\n";
        }
    }

}



// Find the distinct oriented reads that appear on the path
// of a segment. Also return the average edge coverage for the path.
double AssemblyGraph::findOrientedReadsOnSegment(
    uint64_t segmentId,
    vector<OrientedReadId>& orientedReadIdsArgument) const
{
    // Loop over the marker graph path corresponding to this segment.
    const span<const MarkerGraphEdgeId> path = markerGraphPaths[segmentId];
    double coverage = 0.;
    std::set<OrientedReadId> orientedReadIds;
    for(const MarkerGraphEdgeId& edgeId: path) {

        // Loop over the marker intervals for this marker graph edge.
        const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
        coverage += double(markerIntervals.size());
        for(const MarkerInterval& markerInterval: markerIntervals) {
            orientedReadIds.insert(markerInterval.orientedReadId);
        }
    }

    // Copy the oriented reads to the vector passed as an argument.
    orientedReadIdsArgument.clear();
    orientedReadIdsArgument.insert(orientedReadIdsArgument.end(),
        orientedReadIds.begin(), orientedReadIds.end());

    return coverage / double(path.size());
}



// Get information about the oriented reads that appear on the
// marker graph path of a segment.
void AssemblyGraph::getOrientedReadsOnSegment(
    uint64_t segmentId,
    SegmentOrientedReadInformation& information) const
{
    // A data structure that, for each oriented read we find,
    // contains a sum of offsets and the number of marker graph vertices
    // that contributed to the sum.
    std::map<OrientedReadId, pair<uint64_t, int64_t>  > table;

    // Loop over the marker graph path corresponding to this segment.
    const span<const MarkerGraphEdgeId> path = markerGraphPaths[segmentId];
    std::set<OrientedReadId> orientedReadIds;
    for(uint64_t position=0; position<path.size(); position++) {
        const MarkerGraphEdgeId& edgeId = path[position];

        // Loop over the marker intervals for this marker graph edge.
        const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;

            // Update our table for this oriented read.
            auto it = table.find(orientedReadId);
            if(it == table.end()) {
                tie(it, ignore) = table.insert(make_pair(orientedReadId, make_pair(0ULL, 0LL)));
            }
            auto& p = it->second;
            p.first += 2;
            p.second += int32_t(position) - int32_t(markerInterval.ordinals[0]);
            p.second += int32_t(position + 1) -int32_t(markerInterval.ordinals[1]);
        }
    }



    // Store what we found.
    information.infos.clear();
    for(const auto& p: table) {
        SegmentOrientedReadInformation::Info info;
        info.orientedReadId = p.first;
        const uint64_t n = p.second.first;
        const int64_t sum = p.second.second;
        info.averageOffset = int32_t(std::round(double(sum) / double(n)));
        information.infos.push_back(info);
    }
 }



// Estimate the offset between two segments.
// Takes as input SegmentOrientedReadInformation objects
// for the two segments.
// Common oriented reads between the two segments are used
// to estimate the average offset, in markers,
// between the beginning of the segments.
// The number of common oriented reads
// is computed and stored in the last argument.
// If that is zero, the computed offset is not valid.
void AssemblyGraph::estimateOffset(
    const SegmentOrientedReadInformation& info0,
    const SegmentOrientedReadInformation& info1,
    int64_t& offset,
    uint64_t& commonOrientedReadCount
    ) const
{
    offset = 0;
    commonOrientedReadCount = 0;

    // Joint loop over common oriented reads in the two segments.
    const auto begin0 = info0.infos.begin();
    const auto begin1 = info1.infos.begin();
    const auto end0 = info0.infos.end();
    const auto end1 = info1.infos.end();
    auto it0 = begin0;
    auto it1 = begin1;
    while((it0 != end0) and (it1 != end1)) {

        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
        } else if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
        } else {
            SHASTA_ASSERT(it0->orientedReadId == it1->orientedReadId);

            commonOrientedReadCount++;
            offset += (int64_t(it0->averageOffset) - int64_t(it1->averageOffset));

            ++it0;
            ++it1;
        }
    }

    if(commonOrientedReadCount) {
        offset = int64_t(std::round(double(offset) / double(commonOrientedReadCount)));
    } else {
        offset = std::numeric_limits<uint64_t>::max();
    }

}



// Analyze a pair of segments for common oriented reads,
// offsets, missing reads, etc.
void AssemblyGraph::analyzeSegmentPair(
    uint64_t segmentId0,
    uint64_t segmentId1,
    const SegmentOrientedReadInformation& info0,
    const SegmentOrientedReadInformation& info1,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    SegmentPairInformation& info01
    ) const
{
    using boost::icl::discrete_interval;
    using boost::icl::intersects;

    // Store the number of oriented reads in each segment.
    info01.totalCount[0] = info0.infos.size();
    info01.totalCount[1] = info1.infos.size();

    // Use common oriented reads to estimate the offset between the two segments.
    // If there are no common oriented reads, stop here.
    estimateOffset(info0, info1, info01.offset, info01.commonCount);
    if(info01.commonCount == 0) {
        return;
    }


    // Count the oriented reads missing from each segment,
    // and which should have been present based on
    // the known relative offsets.
    info01.unexplainedCount = {0, 0};
    info01.shortCount = {0, 0};

    // Set up a joint loop over oriented reads in the two segments.
    const auto begin0 = info0.infos.begin();
    const auto begin1 = info1.infos.begin();
    const auto end0 = info0.infos.end();
    const auto end1 = info1.infos.end();
    auto it0 = begin0;
    auto it1 = begin1;

    const uint64_t length0 = markerGraphPaths.size(segmentId0);
    const uint64_t length1 = markerGraphPaths.size(segmentId1);
    while(true) {

        // At end of both segments.
        if((it0 == end0) and (it1 == end1)) {
            break;
        }



        // This read only appears on segment 0.
        if((it1 == end1) or ((it0 != end0) and (it0->orientedReadId < it1->orientedReadId))) {
            const int64_t orientedReadLength = markers.size(it0->orientedReadId.getValue());

            // Compute the hypothetical range of the oriented read relative
            // to the beginning of segment 1.
            const discrete_interval<int64_t> orientedReadRange1(
                it0->averageOffset - info01.offset,
                it0->averageOffset - info01.offset + orientedReadLength);
            const discrete_interval<int64_t> segment1Range(0, length1);

            // Figure out if it the oriented read would overlap segment 1.
            const bool wouldOverlap = intersects(orientedReadRange1, segment1Range);

            if(wouldOverlap) {
                ++info01.unexplainedCount[0];
            } else {
                ++info01.shortCount[0];
            }

            SHASTA_ASSERT(it0 != end0);
            ++it0;
        }



        // Only on segment 1
        else if((it0 == end0) or ((it1 != end1) and (it1->orientedReadId < it0->orientedReadId))) {
            const int64_t orientedReadLength = markers.size(it1->orientedReadId.getValue());

            // Compute the hypothetical range of the oriented read relative
            // to the beginning of segment 0.
            const discrete_interval<int64_t> orientedReadRange0(
                it1->averageOffset + info01.offset,
                it1->averageOffset + info01.offset + orientedReadLength);
            const discrete_interval<int64_t> segment0Range(0, length0);

            // Figure out if it the oriented read would overlap segment 0.
            const bool wouldOverlap = intersects(orientedReadRange0, segment0Range);

            if(wouldOverlap) {
                ++info01.unexplainedCount[1];
            } else {
                ++info01.shortCount[1];
            }

            SHASTA_ASSERT(it1 != end1);
            ++it1;
        }

        // On both segments.
        else {
            SHASTA_ASSERT(it0 != end0);
            SHASTA_ASSERT(it1 != end1);
            ++it0;
            ++it1;
        }
    }

    info01.check();

}



// Gather oriented read information for each segment.
void AssemblyGraph::storeSegmentOrientedReadInformation(size_t threadCount)
{
    const uint64_t segmentCount = markerGraphPaths.size();
    segmentOrientedReadInformation.resize(segmentCount);
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::storeSegmentOrientedReadInformationThreadFunction, threadCount);
}




// Gather oriented read information for each segment.
void AssemblyGraph::storeSegmentOrientedReadInformationThreadFunction(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {

            // Get oriented read information for this segment.
            getOrientedReadsOnSegment(segmentId, segmentOrientedReadInformation[segmentId]);
        }
    }
}



#if 0
void AssemblyGraph::clusterSegments(size_t threadCount, uint64_t minClusterSize)
{
    // Gather oriented read information for all segments.
    storeSegmentOrientedReadInformation(threadCount);

    // Find the segment pairs.
    const uint64_t segmentCount = markerGraphPaths.size();
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    clusterSegmentsData.threadPairs.resize(threadCount);
    runThreads(&AssemblyGraph::clusterSegmentsThreadFunction1, threadCount);

    // For now, write a dot file with the pairs.
    ofstream dot("SegmentGraph.dot");
    dot << "graph segmentGraph {\n";
    for(const auto& threadPairs: clusterSegmentsData.threadPairs) {
        for(const auto& p: threadPairs) {
            dot << p.first << "--" << p.second << ";\n";
        }
    }
    dot << "}\n";



    // The segment pairs we found define a subgraph of the assembly graph.
    // Compute connected components of this subgraph.
    // The connected components of sufficient size become clusters.
    vector<uint64_t> rank(segmentCount);
    vector<uint64_t> parent(segmentCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        disjointSets.make_set(segmentId);
    }
    for(const auto& threadPairs: clusterSegmentsData.threadPairs) {
        for(const auto& p: threadPairs) {
            disjointSets.union_set(p.first, p.second);
        }
    }

    // Gather the segments in each connected component.
    vector< vector<uint64_t> > components(segmentCount);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const uint64_t componentId = disjointSets.find_set(segmentId);
        components[componentId].push_back(segmentId);
    }

    // Each connected components of size at least minClusterSize
    // becomes a cluster.
    vector< pair<uint64_t, uint64_t> > clusterTable;
    for(uint64_t componentId=0; componentId<segmentCount; componentId++) {
        const vector<uint64_t>& component = components[componentId];
        const uint64_t componentSize = component.size();
        if(component.size() >= minClusterSize) {
            clusterTable.push_back(make_pair(componentId, componentSize));
        }
    }

    // Sort the clusters by decreasing size.
    sort(clusterTable.begin(), clusterTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    cout << "Found " << clusterTable.size() << " segment clusters with the following sizes:" << endl;
    uint64_t clusteredSegmentCount = 0;
    for(uint64_t clusterId=0; clusterId<clusterTable.size(); clusterId++) {
        const auto& p = clusterTable[clusterId];
        const uint64_t componentSize = p.second;
        cout << " " << componentSize;
        clusteredSegmentCount += componentSize;
    }
    cout << endl;
    cout << "Out of " << segmentCount << " segments, " <<
        clusteredSegmentCount << " were assigned to a cluster." << endl;



    // Store the cluster id of each segment.
    createNew(clusterIds, "Mode3-ClusterIds");
    clusterIds.resize(segmentCount);
    fill(clusterIds.begin(), clusterIds.end(), std::numeric_limits<uint64_t>::max());
    for(uint64_t clusterId=0; clusterId<clusterTable.size(); clusterId++) {
        const auto& p = clusterTable[clusterId];
        const uint64_t componentId = p.first;
        const vector<uint64_t>& cluster = components[componentId];
        for(const uint64_t segmentId: cluster) {
            clusterIds[segmentId] = clusterId;
        }
    }



    // Clean up.
    clusterSegmentsData.threadPairs.clear();
    clusterSegmentsData.threadPairs.shrink_to_fit();
    segmentOrientedReadInformation.clear();
    segmentOrientedReadInformation.shrink_to_fit();
}



void AssemblyGraph::clusterSegmentsThreadFunction1(size_t threadId)
{

    auto& threadPairs = clusterSegmentsData.threadPairs[threadId];
    threadPairs.clear();
    vector<uint64_t> descendants;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segments assigned to this batch.
        for(uint64_t segmentId0=begin; segmentId0!=end; ++segmentId0) {

            // Add pairs for which the lowest numbered segment is segmentId0.
            addClusterPairs(threadId, segmentId0);
        }
    }
}



void AssemblyGraph::addClusterPairs(size_t threadId, uint64_t startSegmentId)
{
    // EXPOSE THESE CONSTANTS WHEN CODE STABILIZES.
    const uint64_t minCommonReadCount = 10;
    const double maxUnexplainedFraction = 0.25;
    const double minJaccard = 0.7;
    const uint64_t pairCountPerSegment = 1;
    const uint64_t maxDistance = 200;

    // std::lock_guard<std::mutex> lock(mutex);    // *********** TAKE OUT

    // Do a BFS and check each pair as we encounter it.
    // The BFS terminates when we found enough pairs.

    // Do the BFS in both directions.
    for(uint64_t direction=0; direction<1; direction++) { // ********* ONE DIRECTION ONLY
        // cout << startSegmentId << " direction " << direction << endl;

        // Initialize the BFS.
        std::queue<uint64_t> q;
        q.push(startSegmentId);
        std::map<uint64_t, uint64_t> distanceMap;
        distanceMap.insert(make_pair(startSegmentId, 0));
        uint64_t foundCount = 0;

        // BFS loop.
        while(not q.empty()) {
            const uint64_t segmentId0 = q.front();
            // cout << "Dequeued " << segmentId0 << endl;
            q.pop();

            const uint64_t distance0 = distanceMap[segmentId0];
            const uint64_t distance1 = distance0 + 1;

            // Loop over children or parents of segmentId0.
            const auto neighbors = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
            for(const uint64_t linkId01: neighbors) {
                const Link& link01 = links[linkId01];
                const uint64_t segmentId1 = (direction==0) ? link01.segmentId1 : link01.segmentId0;

                // If we already encountered segmentId1, skip it.
                if(distanceMap.find(segmentId1) != distanceMap.end()) {
                    continue;
                }

                // Enqueue it.
                if(distance1 < maxDistance) {
                    q.push(segmentId1);
                }
                distanceMap.insert(make_pair(segmentId1, distance1));

                // cout << "Found " << segmentId1 << endl;

                // Check the pair (startSegmentId, segmentId1).
                SegmentPairInformation info;
                analyzeSegmentPair(startSegmentId, segmentId1,
                    segmentOrientedReadInformation[startSegmentId],
                    segmentOrientedReadInformation[segmentId1],
                    markers, info);
                if(info.commonCount < minCommonReadCount) {
                    continue;
                }
                if(info.maximumUnexplainedFraction() > maxUnexplainedFraction) {
                    continue;
                }
                if(info.jaccard() < minJaccard) {
                    continue;
                }

                // Store it.
                // cout << "Stored " << segmentId1 << endl;
                clusterSegmentsData.threadPairs[threadId].push_back(make_pair(startSegmentId, segmentId1));
                ++foundCount;
                if(foundCount >= pairCountPerSegment) {
                    break;
                }
            }

            if(foundCount >= pairCountPerSegment) {
                break;
            }
        }
    }
}
#endif



// Find descendants of a given segment, up to a given distance in the graph.
void AssemblyGraph::findDescendants(
    uint64_t startSegmentId,
    uint64_t maxDistance,
    vector<uint64_t>& descendants
    ) const
{
    // Initialize the BFS.
    descendants.clear();
    std::queue<uint64_t> q;
    q.push(startSegmentId);
    std::map<uint64_t, uint64_t> distanceMap;
    distanceMap.insert(make_pair(startSegmentId, 0));

    // BFS loop.
    while(not q.empty()) {
        const uint64_t segmentId0 = q.front();
        q.pop();

        const uint64_t distance0 = distanceMap[segmentId0];
        const uint64_t distance1 = distance0 + 1;

        // Loop over children of segmentId0.
        for(const uint64_t linkId01: linksBySource[segmentId0]) {
            const Link& link01 = links[linkId01];
            const uint64_t segmentId1 = link01.segmentId1;

            // If we already encountered segmentId1, skip it.
            if(distanceMap.find(segmentId1) != distanceMap.end()) {
                continue;
            }

            descendants.push_back(segmentId1);
            distanceMap.insert(make_pair(segmentId1, distance1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }
    }
}



void AssemblyGraph::analyzeSubgraph(
    const vector<uint64_t>& segmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    if(segmentIds.size() <= 64) {
        analyzeSubgraphTemplate<64>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 128) {
        analyzeSubgraphTemplate<128>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 192) {
        analyzeSubgraphTemplate<192>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 256) {
        analyzeSubgraphTemplate<256>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 320) {
        analyzeSubgraphTemplate<320>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 384) {
        analyzeSubgraphTemplate<384>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 448) {
        analyzeSubgraphTemplate<448>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 512) {
        analyzeSubgraphTemplate<512>(segmentIds, clusters, debug);
    } else {
        SHASTA_ASSERT(0);
    }
}



template<uint64_t N> void AssemblyGraph::analyzeSubgraphTemplate(
    const vector<uint64_t>& unsortedSegmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const double fractionThreshold = 0.05;
    const uint64_t minClusterCoverage = 6;
    const uint64_t minSegmentCoverage = 6;

    using BitVector = std::bitset<N>;
    using JourneySnippet = AnalyzeSubgraphClasses::JourneySnippet;
    using Cluster = AnalyzeSubgraphClasses::Cluster;
    using SnippetGraphVertex = AnalyzeSubgraphClasses::SnippetGraphVertex;
    using SnippetGraph = AnalyzeSubgraphClasses::SnippetGraph;
    using vertex_descriptor = SnippetGraph::vertex_descriptor;

    // Create a sorted version of the segmentIds. We will need it later.
    vector<uint64_t> segmentIds = unsortedSegmentIds;
    sort(segmentIds.begin(), segmentIds.end());

    // Gather triplets (orientedReadId, position in assembly graph journey, segmentId).
    using Triplet = tuple<OrientedReadId, uint64_t, uint64_t>;
    vector<Triplet> triplets;
    for(const uint64_t segmentId: segmentIds) {
        const auto v = assemblyGraphJourneyInfos[segmentId];
        for(const auto& p: v) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t position = p.second;
            triplets.push_back(Triplet(orientedReadId, position, segmentId));
        }
    }
    sort(triplets.begin(), triplets.end());

    // Write the triplets.
    if(debug) {
        ofstream csv("Triplets.csv");
        for(const Triplet& triplet: triplets) {
            csv << get<0>(triplet) << ",";
            csv << get<1>(triplet) << ",";
            csv << get<2>(triplet) << "\n";
        }
    }



    // Find streaks for the same OrientedReadId where the position
    // increases by 1 each time.
    // Each streak generates a JourneySnippet.
    vector<JourneySnippet> snippets;
    for(uint64_t i=0; i<triplets.size(); /* Increment later */) {
        const OrientedReadId orientedReadId = get<0>(triplets[i]);

        // Find this streak.
        uint64_t streakBegin = i;
        uint64_t streakEnd = streakBegin + 1;
        for(; streakEnd<triplets.size(); streakEnd++) {
            if(get<0>(triplets[streakEnd]) != orientedReadId) {
                break;
            }
            if(get<1>(triplets[streakEnd]) != get<1>(triplets[streakEnd-1]) + 1) {
                break;
            }
        }

        // Add a snippet.
        JourneySnippet snippet;
        snippet.orientedReadId = orientedReadId;
        snippet.firstPosition = get<1>(triplets[streakBegin]);
        for(uint64_t j=streakBegin; j!=streakEnd; ++j) {
            snippet.segmentIds.push_back(get<2>(triplets[j]));
        }
        snippets.push_back(snippet);

        // Prepare to process the next streak.
        i = streakEnd;
    }



    // Write the snippets.
    if(debug) {
        ofstream csv("JourneySnippets.csv");
        csv << "SnippetIndex,OrientedReadId,First position,LastPosition,SegmentIds\n";
        for(uint64_t snippetIndex=0; snippetIndex<snippets.size(); snippetIndex++) {
            const JourneySnippet& snippet = snippets[snippetIndex];
            csv << snippetIndex << ",";
            csv << snippet.orientedReadId << ",";
            csv << snippet.firstPosition << ",";
            csv << snippet.lastPosition() << ",";
            for(const uint64_t segmentId: snippet.segmentIds) {
                csv << segmentId << ",";
            }
            csv << "\n";
        }
    }



    // For each snippet, create a BitVector that describes the segments
    // the snippet visits.
    const uint64_t snippetCount = snippets.size();
    vector<BitVector> bitVectors(snippetCount);
    vector<uint64_t> bitVectorsPopCount(snippetCount);
    for(uint64_t snippetIndex=0; snippetIndex<snippetCount; snippetIndex++) {
        const JourneySnippet& snippet = snippets[snippetIndex];
        BitVector& bitVector = bitVectors[snippetIndex];

        for(const uint64_t segmentId: snippet.segmentIds) {
            auto it = lower_bound(segmentIds.begin(), segmentIds.end(), segmentId);
            SHASTA_ASSERT(it != segmentIds.end());
            SHASTA_ASSERT(*it == segmentId);
            const uint64_t bitIndex = it - segmentIds.begin();
            bitVector.set(bitIndex);
        }
        bitVectorsPopCount[snippetIndex] = bitVector.count();
    }



    // Create the SnippetGraph.
    // A vertex represents a set of snippets and stores
    // the corresponding snippet indexes.
    // An edge x->y is created if there is at least one snippet in y
    // that is an approximate subset of a snippet in x.
    // We express this condition as |y-x| < fractionThreshold * |y|
    // We start with one snippet per vertex.
    SnippetGraph graph;
    vector<vertex_descriptor> vertexTable;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    for(uint64_t snippetIndex=0; snippetIndex<snippetCount; snippetIndex++) {
        const auto v = add_vertex(SnippetGraphVertex(snippetIndex), graph);
        vertexTable.push_back(v);
        vertexMap.insert(make_pair(v, snippetIndex));
    }
    for(uint64_t iy=0; iy<snippetCount; iy++) {
        const BitVector& y = bitVectors[iy];
        const uint64_t threshold = uint64_t(std::round(fractionThreshold * double(bitVectorsPopCount[iy])));
        const vertex_descriptor vy = vertexTable[iy];
        for(uint64_t ix=0; ix<snippetCount; ix++) {
            if(ix == iy) {
                continue;
            }
            const BitVector& x = bitVectors[ix];

            // Compute z = y-x.
            BitVector z = y;
            z &= (~x);

            if(z.count() <= threshold) {
                const vertex_descriptor vx = vertexTable[ix];
                add_edge(vx, vy, graph);
            }
        }
    }



    // Compute strongly connected components of the SnippetGraph.
    std::map<vertex_descriptor, uint64_t> componentMap;
    const uint64_t componentCount = boost::strong_components(
        graph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));
    // cout << "Found " << componentCount << " strongly connected components." << endl;

    // Gather the vertices of each strongly connected component.
    vector< vector<vertex_descriptor> > components(componentCount);
    BGL_FORALL_VERTICES_T(v, graph, SnippetGraph) {
        const uint64_t componentId = componentMap[v];
        SHASTA_ASSERT(componentId < componentCount);
        components[componentId].push_back(v);
    }
    if(false) {
        cout << "Strongly connected components:\n";
        for(uint64_t componentId=0; componentId<componentCount; componentId++) {
            cout << componentId << ": ";
            for(const vertex_descriptor v: components[componentId]) {
                cout << vertexMap[v] << " ";
            }
            cout << "\n";
        }
    }



    // Condense the strongly connected components.
    // After this, the SnippetGraph is guaranteed to be acyclic.
    for(const vector<vertex_descriptor>& component: components) {
        if(component.size() == 1) {
            continue;
        }

        // Create a new vertex to represent this component.
        const auto vNew = add_vertex(graph);
        vector<uint64_t>& snippetsNew = graph[vNew].snippetIndexes;
        for(const vertex_descriptor v: component) {
            const vector<uint64_t>& snippets = graph[v].snippetIndexes;
            SHASTA_ASSERT(snippets.size() == 1);
            snippetsNew.push_back(snippets.front());
        }

        // Create the new edges.
        for(const vertex_descriptor v0: component) {

            // Out-edges.
            BGL_FORALL_OUTEDGES_T(v0, e01, graph, SnippetGraph) {
                const vertex_descriptor v1 = target(e01, graph);
                if(v1 != vNew) {
                    add_edge(vNew, v1,graph);
                }
            }

            // In-edges.
            BGL_FORALL_INEDGES_T(v0, e10, graph, SnippetGraph) {
                const vertex_descriptor v1 = source(e10, graph);
                if(v1 != vNew) {
                    add_edge(v1, vNew, graph);
                }
            }
        }

        // Remove the old vertices and their edges.
        for(const vertex_descriptor v: component) {
            clear_vertex(v, graph);
            remove_vertex(v, graph);
        }
    }


    // Compute which maximal vertices each vertex is a descendant of.
    std::map<vertex_descriptor, vector<vertex_descriptor> > ancestorMap;
    BGL_FORALL_VERTICES_T(v0, graph, SnippetGraph) {
        if(in_degree(v0, graph) != 0) {
            continue;   // Not a maximal vertex.
        }

        // Find the descendants of this maximal vertex.
        vector<vertex_descriptor> descendants;
        graph.findDescendants(v0, descendants);

        // Update the ancestor map.
        for(const vertex_descriptor v1: descendants) {
            ancestorMap[v1].push_back(v0);
        }
    }



    // Each maximal vertex generates a cluster consisting of the vertices
    // that descend from it and from no other maximal vertex.
    // Gather the vertices in each cluster.
    std::map<vertex_descriptor, vector<vertex_descriptor> > clusterMap;
    uint64_t unclusterVertexCount = 0;
    BGL_FORALL_VERTICES_T(v1, graph, SnippetGraph) {
        const vector<vertex_descriptor>& ancestors = ancestorMap[v1];
        if(ancestors.size() == 1) {
            const vertex_descriptor v0 = ancestors.front();
            clusterMap[v0].push_back(v1);
        } else {
            ++unclusterVertexCount;
        }
    }
    cout << "Found " << unclusterVertexCount << " unclustered vertices." << endl;



    // Gather the snippets in each cluster.
    clusters.clear();
    for(const auto& p: clusterMap) {
        const vector<vertex_descriptor>& clusterVertices = p.second;
        clusters.resize(clusters.size() + 1);
        Cluster& cluster = clusters.back();

        vector<uint64_t> clusterSnippetIndexes; // Only used for debug output.
        for(const vertex_descriptor v: clusterVertices) {
            const vector<uint64_t>& snippetIndexes = graph[v].snippetIndexes;
            for(const uint64_t snippetIndex: snippetIndexes) {
                cluster.snippets.push_back(snippets[snippetIndex]);
                clusterSnippetIndexes.push_back(snippetIndex);
            }
        }
        cluster.constructSegments();
        cluster.cleanupSegments(minSegmentCoverage);
        cout << "Found a cluster candidate with " <<
            clusterVertices.size() << " vertices and " <<
            cluster.snippets.size() << " snippets:" << endl;
        for(const uint64_t snippetIndex: clusterSnippetIndexes) {
            cout << snippetIndex << " ";
        }
        cout << endl;

        // If coverage on this cluster is too low, discard it.
        if(cluster.coverage() < minClusterCoverage) {
            clusters.resize(clusters.size() - 1);
            cout << "This cluster candidate was discarded because of low coverage." << endl;
            continue;
        }

        // This cluster will be stored and is assigned this clusterId;
        const uint64_t clusterId = clusters.size() - 1;

        if(debug) {

            cout << "This cluster was stored as cluster " << clusterId << endl;
            cout << "Segment(coverage) for this cluster:\n";
            for(const auto& p: cluster.segments) {
                cout << p.first << "(" << p.second << ") ";
            }
            cout << endl;
        }

        // Mark the vertices of this cluster.
        for(const vertex_descriptor v: clusterVertices) {
            graph[v].clusterId = clusterId;
        }
    }
    graph.clusterCount = clusters.size();



    // Write out the SnippetGraph.
    if(debug) {
        graph.writeGraphviz("SnippetGraph.dot");
    }
}



void AssemblyGraph::AnalyzeSubgraphClasses::SnippetGraph::findDescendants(
    const vertex_descriptor vStart,
    vector<vertex_descriptor>& descendants) const
{
    const SnippetGraph& graph = *this;

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    std::set<vertex_descriptor> descendantsSet;
    descendantsSet.insert(vStart);

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(v0, e01, graph, SnippetGraph) {
            const vertex_descriptor v1 = target(e01, graph);
            if(descendantsSet.find(v1) == descendantsSet.end()) {
                q.push(v1);
                descendantsSet.insert(v1);
            }
        }
    }

    descendants.clear();
    copy(descendantsSet.begin(), descendantsSet.end(), back_inserter(descendants));
}



void AssemblyGraph::AnalyzeSubgraphClasses::SnippetGraph::writeGraphviz(
    const string& fileName) const
{
    const SnippetGraph& graph = *this;

    ofstream dot(fileName);
    dot << "digraph SnippetGraph{\n"
        "node [shape=rectangle];\n";
    BGL_FORALL_VERTICES(v, graph, SnippetGraph) {
        dot << "\"" << v << "\" [label=\"";
        const vector<uint64_t>& snippetIndexes = graph[v].snippetIndexes;
        for(const uint64_t snippetIndex: snippetIndexes) {
            dot << snippetIndex;
            dot << "\\n";
        }
        dot << "\"";
        const uint64_t clusterId = graph[v].clusterId;
        if(clusterId != std::numeric_limits<uint64_t>::max()) {
            dot << " style=filled fillcolor=\"" <<
                float(clusterId)/float(clusterCount) <<
                ",0.3,1\"";
        }
        dot << "];\n";
    }
    BGL_FORALL_EDGES(e, graph, SnippetGraph) {
        const vertex_descriptor vx = source(e, graph);
        const vertex_descriptor vy = target(e, graph);
        dot << "\"" << vx << "\"->\"" << vy << "\";\n";
    }
    dot << "}\n";

}



void AssemblyGraph::AnalyzeSubgraphClasses::Cluster::constructSegments()
{
    // A map with Key=segmentId, value = coverage.
    std::map<uint64_t, uint64_t> segmentMap;

    for(const JourneySnippet& snippet: snippets) {
        for(const uint64_t segmentId: snippet.segmentIds) {
            auto it = segmentMap.find(segmentId);
            if(it == segmentMap.end()) {
                segmentMap.insert(make_pair(segmentId, 1));
            } else {
                ++(it->second);
            }
        }
    }

    segments.clear();
    copy(segmentMap.begin(), segmentMap.end(), back_inserter(segments));
}



void AssemblyGraph::AnalyzeSubgraphClasses::Cluster::cleanupSegments(uint64_t minSegmentCoverage)
{
    vector< pair<uint64_t, uint64_t > > newSegments;
    for(const auto& p: segments) {
        if(p.second >= minSegmentCoverage) {
            newSegments.push_back(p);
        }
    }
    segments.swap(newSegments);
}



vector<uint64_t> AssemblyGraph::AnalyzeSubgraphClasses::Cluster::getSegments() const
{
    vector<uint64_t> v;
    for(const auto& p: segments) {
        v.push_back(p.first);
    }
    return v;
}



// Create an assembly path starting at a given segment.
void AssemblyGraph::createAssemblyPath(
    uint64_t startSegmentId,
    uint64_t direction,    // 0 = forward, 1 = backward
    AssemblyPath& path
    ) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonForLink = 3;
    const uint64_t minCommonForReference = 3;
    const double minJaccard = 0.75;
    const int32_t minLinkSeparation = -20;

    const bool debug = false;
    if(true) {
        cout << timestamp << "createAssemblyPath begins at segment " << startSegmentId <<
            ", direction " << direction << endl;
    }



    // At each iteration, we start from segmentIdA (the current "primary segment")
    // and move in the specified direction until we find segmentIdB with
    // sufficiently high Jaccard similarity and number of
    // common oriented reads with segmentIdA.
    // At each step, we choose the links that has the most common oriented
    // reads with the current primary segment.
    uint64_t referenceSegmentId = startSegmentId;
    SegmentOrientedReadInformation infoReference;
    getOrientedReadsOnSegment(referenceSegmentId, infoReference);
    uint64_t segmentId0 = startSegmentId;
    path.clear();
    path.segments.push_back(AssemblyPathSegment(startSegmentId, true));
    vector<uint64_t> lastIterationSegments;
    std::set< pair<uint64_t, uint64_t> > previousPairs;  // (reference segment, current segment).
    while(true) {

        if(debug) {
            cout << "Reference segment " << referenceSegmentId <<
                ", segmentId0 " << segmentId0 << endl;
        }

        // Loop over outgoing or incoming links of the current segment.
        // Find the link with the most common reads with the reference segment.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        if(linkIds.empty()) {
            if(debug) {
                cout << "No links in this direction." << endl;
            }
            break;
        }
        uint64_t linkIdBest = invalid<uint64_t>;
        uint64_t commonOrientedReadCountBest = 0;
        for(const uint64_t linkId: linkIds) {

            // If link separation is too negative, skip it.
            // The goal here is to avoid cycles in paths.
            const Link& link = links[linkId];
            if(link.separation < minLinkSeparation) {
                continue;
            }

            // Count the number of common oriented reads between the reference segment and this link.
            uint64_t commonOrientedReadCount;
            analyzeSegmentLinkPair(referenceSegmentId, linkId, commonOrientedReadCount);

            // If better than the one we have it, record it.
            if(commonOrientedReadCount > commonOrientedReadCountBest) {
                linkIdBest = linkId;
                commonOrientedReadCountBest = commonOrientedReadCount;
            }
        }
        if(commonOrientedReadCountBest < minCommonForLink) {
            if(debug) {
                cout << "No good links found." << endl;
            }
            break;
        }
        const uint64_t linkId = linkIdBest;
        if(debug) {
            cout << "Best link " << linkId <<
                ", " << commonOrientedReadCountBest <<
                " common oriented reads with the reference segment." << endl;
        }

        // Get the segment at the other side of this link.
        const Link& link = links[linkId];
        const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;
        if(debug) {
            cout << "segmentId1 " << segmentId1 << endl;
        }
        lastIterationSegments.push_back(segmentId1);

        // Check that we haven't been here before.
        if(previousPairs.contains(make_pair(referenceSegmentId, segmentId1))) {
            break;
        }
        previousPairs.insert(make_pair(referenceSegmentId, segmentId1));

        // Check segmentId1 against the reference segment.
        SegmentOrientedReadInformation info1;
        getOrientedReadsOnSegment(segmentId1, info1);
        SegmentPairInformation info;
        analyzeSegmentPair(
            referenceSegmentId, segmentId1,
            infoReference, info1,
            markers, info);
        if(debug) {
            cout << "Jaccard " << info.jaccard() << endl;
        }

        // If the Jaccard similarity is high, this becomes the new reference segment.
        if(info.commonCount >= minCommonForReference and info.jaccard() >= minJaccard) {
            referenceSegmentId = segmentId1;
            getOrientedReadsOnSegment(referenceSegmentId, infoReference);
            const uint64_t lastPrimarySegmentId = path.segments.back().id;
            if(debug) {
                cout << "New reference segment is " << segmentId1 << endl;
                cout << "Previous reference segment is " << lastPrimarySegmentId << endl;
            }
            for(const uint64_t segmentId: lastIterationSegments) {
                path.segments.push_back(AssemblyPathSegment(segmentId, false));
                if(debug) {
                    cout << "Added segment " << segmentId << " to path." << endl;
                }
                if(segmentId != segmentId1) {
                    if(direction == 0) {
                        path.segments.back().previousPrimarySegmentId = lastPrimarySegmentId;
                        path.segments.back().nextPrimarySegmentId = segmentId1;
                    } else {
                        path.segments.back().previousPrimarySegmentId = segmentId1;
                        path.segments.back().nextPrimarySegmentId = lastPrimarySegmentId;
                    }
                }
            }
            path.segments.back().isPrimary = true;
            lastIterationSegments.clear();
        }

        segmentId0 = segmentId1;
    }



    if(true) {
        cout << timestamp << "createAssemblyPath3 ends." << endl;
    }
}



// Count the number of common oriented reads between a segment and a link,
// without counting oriented reads that appear more than once on the
// segment or on the link.
void AssemblyGraph::analyzeSegmentLinkPair(
    uint64_t segmentId,
    uint64_t linkId,
    uint64_t& commonOrientedReadCount
) const
{
    // The oriented reads in this segment,
    // with some extra information that we don't care about here.
    const auto segmentOrientedReads = assemblyGraphJourneyInfos[segmentId];

    // The oriented reads in this link,
    // with some extra information that we don't care about here.
    const auto linkOrientedReads = transitions[linkId];

    // Joint loop over oriented reads.
    commonOrientedReadCount = 0;
    const auto segmentBegin = segmentOrientedReads.begin();
    const auto segmentEnd = segmentOrientedReads.end();
    const auto linkBegin = linkOrientedReads.begin();
    const auto linkEnd = linkOrientedReads.end();
    auto itSegment = segmentBegin;
    auto itLink = linkBegin;
    while(itSegment != segmentEnd and itLink != linkEnd) {

        if(itSegment->first < itLink->first) {
            ++itSegment;
            continue;
        }
        if(itLink->first < itSegment->first) {
            ++itLink;
            continue;
        }
        SHASTA_ASSERT(itSegment->first == itLink->first);

        // If it appears more than once in the segment, skip it.
        auto itSegmentNext = itSegment + 1;
        if(itSegmentNext != segmentEnd and  itSegmentNext->first == itSegment->first) {
            ++itSegment;
            ++itLink;
            continue;
        }
        if(itSegment != segmentBegin) {
            auto itSegmentPrevious = itSegment - 1;
            if(itSegmentPrevious->first == itSegment->first) {
                ++itSegment;
                ++itLink;
                continue;
            }
        }

        // If it appears more than once in the link, skip it.
        auto itLinkNext = itLink + 1;
        if(itLinkNext != linkEnd and  itLinkNext->first == itLink->first) {
            ++itSegment;
            ++itLink;
            continue;
        }
        if(itLink != linkBegin) {
             auto itLinkPrevious = itLink - 1;
            if(itLinkPrevious->first == itLink->first) {
                ++itSegment;
                ++itLink;
                continue;
            }
        }

        // Ok, this is a common oriented read that appears only once
        // in both the segment and the link.
        ++commonOrientedReadCount;
        ++itSegment;
        ++itLink;
    }

}



// Given a segment, use a BFS to move in the specified direction until
// we find a segment with sufficiently high Jaccard similarity
// and number of common reads.
// This returns invalid<uint64_t> if no such segment is found
// within the specified distance.
uint64_t AssemblyGraph::findSimilarSegmentBfs(
    uint64_t segmentIdA,
    uint64_t direction, // 0 = forward, 1 = backward
    uint64_t maxDistance,
    uint64_t minCommon,
    double minJaccard) const
{
    const bool debug = true;
    if(debug) {
        cout << "findSimilarSegmentBfs starts " << segmentIdA << " " << direction << endl;
    }

    // Sanity check.
    SHASTA_ASSERT(maxDistance > 0);

    // Get the oriented reads on segmentIdA.
    SegmentOrientedReadInformation infoA;
    getOrientedReadsOnSegment(segmentIdA, infoA);

    // Initialize a BFS.
    std::queue<uint64_t> q;
    q.push(segmentIdA);

    // Keep track of segments we already encountered and their distance.
    // Key = segmentId;
    // Value = distance.
    std::map<uint64_t, uint64_t> distanceMap;
    distanceMap.insert(make_pair(segmentIdA, 0));



    // BFS loop.
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();
        const uint64_t distance0 = distanceMap[segmentId0];
        SHASTA_ASSERT(distance0 < maxDistance);
        const uint64_t distance1 = distance0 + 1;
        if(debug) {
            cout << "dequeued " << segmentId0 << " " << distance0 << endl;
        }

        // Loop over outgoing or incoming links.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        for(const uint64_t linkId: linkIds) {
            const Link& link = links[linkId];

            // Get the segment at the other side of this link.
            const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;

            // If we already found it, skip it.
            if(distanceMap.contains(segmentId1)) {
                continue;
            }

            if(debug) {
                cout << "found " << segmentId1 << " " << distance1 << endl;
            }

            // Get the oriented reads on segmentId1.
            SegmentOrientedReadInformation info1;
            getOrientedReadsOnSegment(segmentId1, info1);

            // See how similar this is to segmentIdA.
            SegmentPairInformation infoA1;
            analyzeSegmentPair(
                segmentIdA, segmentId1,
                infoA, info1,
                markers, infoA1);

            // If this satisfies our criteria, we are done.
            if(infoA1.commonCount >= minCommon and
                infoA1.jaccard() >= minJaccard) {
                if(debug) {
                    cout << "findSimilarSegmentBFS returns " << segmentId1 << " " << direction << endl;
                }
                return segmentId1;
            }

            // This segment did not satisfy our criteria, so we
            // have to continue the BFS.
            if(distance1 < maxDistance) {
                q.push(segmentId1);
                distanceMap.insert(make_pair(segmentId1, distance1));
                if(debug) {
                    cout << "enqueued " << segmentId1 << " " << distance1 << endl;
                }
            }
        }

    }



    // If getting here, we did not find a segment that satisfies
    // the requested criteria.
    if(debug) {
        cout << "findSimilarSegmentBfs returns invalid" << endl;
    }
    return invalid<uint64_t>;
}



// Given a segment, move in the specified direction,
// in order of increasing distance in markers, until
// we find a segment with sufficiently high Jaccard similarity
// and number of common reads.
// This returns invalid<uint64_t> if no such segment is found
// within the specified distance.
uint64_t AssemblyGraph::findSimilarSegment(
    uint64_t segmentIdA,
    uint64_t direction,     // 0 = forward, 1 = backward
    uint64_t maxDistance,   // In markers
    uint64_t minLinkCoverage,
    int32_t minLinkSeparation,
    uint64_t minCommon,
    double maxUnexplainedFraction,
    double minJaccard,
    vector<uint64_t>& segments) const
{
    const bool debug = false;
    if(debug) {
        cout << "findSimilarSegment begins, segmentIdA " << segmentIdA << endl;
    }
    // Sanity check.
    SHASTA_ASSERT(maxDistance > 0);

    segments.clear();

    // Get the oriented reads on segmentIdA.
    SegmentOrientedReadInformation infoA;
    getOrientedReadsOnSegment(segmentIdA, infoA);

    // (Offset, segmentId) for queued segments.
    std::multimap<uint64_t, uint64_t> q;
    q.insert(make_pair(0, segmentIdA));

    // The segments that we already encountered.
    std::set<uint64_t> visitedSegmentSet;
    visitedSegmentSet.insert(segmentIdA);



    // Search loop.
    while(not q.empty()) {

        // Dequeue the segment with the smallest offset.
        const auto it0 = q.begin();
        const uint64_t segmentId0 = it0->second;
        q.erase(it0);

        // Analyze against segmentIdA.
        SegmentOrientedReadInformation info0;
        getOrientedReadsOnSegment(segmentId0, info0);
        SegmentPairInformation infoA0;
        analyzeSegmentPair(
            segmentIdA, segmentId0,
            infoA, info0,
            markers, infoA0);

        // Add it to our list of segments, if possible.
        const double unexplainedFraction = infoA0.unexplainedFraction(0);
        if(unexplainedFraction < maxUnexplainedFraction) {
            segments.push_back(segmentId0);
        }

        // If unexplained fraction and Jaccard similarity are low, we are done.
        if(segmentId0 != segmentIdA) {
            if((unexplainedFraction < maxUnexplainedFraction) and (infoA0.jaccard() >= minJaccard)) {
                SHASTA_ASSERT(segments.back() == segmentId0);
                return segmentId0;
            }
        }

        if(debug) {
            cout << "Dequeued " << segmentId0 << endl;
        }

        // Loop over outgoing or incoming links.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        for(const uint64_t linkId: linkIds) {
            // If link coverage is too low, skip.
            if(transitions.size(linkId) < minLinkCoverage) {
                continue;
            }

            // If link separation is too negative, skip it.
            // The goal here is to avoid cycles in paths.
            const Link& link = links[linkId];
            if(link.separation < minLinkSeparation) {
                continue;
            }

            // Get the segment at the other side of this link.
            const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;
            if(debug) {
                cout << "Found " << segmentId1 << endl;
            }

            // If we already found it, skip it.
            if(visitedSegmentSet.contains(segmentId1)) {
                if(debug) {
                    cout << "Already found, skipping." << endl;
                }
                continue;
            }
            visitedSegmentSet.insert(segmentId1);

            // Get the oriented reads on segmentId1.
            SegmentOrientedReadInformation info1;
            getOrientedReadsOnSegment(segmentId1, info1);

            // Analyze similarity to segmentIdA.
            SegmentPairInformation infoA1;
            analyzeSegmentPair(
                segmentIdA, segmentId1,
                infoA, info1,
                markers, infoA1);

            // If not enough common segments, skip it.
            if(infoA1.commonCount < minCommon) {
                if(debug) {
                    cout << "Not enough common reads." << endl;
                }
                continue;
            }

            // Offset estimates are not reliable.
            // Don't use them to rule out segments.
#if 0
            // If not in the expected direction, skip it.
            uint64_t offset;
            if(direction == 0) {
                if(infoA1.offset < 0) {
                    if(debug) {
                        cout << "Not in the forward direction." << endl;
                    }
                    continue;
                } else {
                    offset = uint64_t(infoA1.offset);
                }
            } else {
                if(infoA1.offset > 0) {
                    if(debug) {
                        cout << "Not in the backward direction." << endl;
                    }
                    continue;
                } else {
                    offset = uint64_t(-infoA1.offset);
                }
            }
#endif

            // If we went too far, skip it.
            if(labs(infoA1.offset) > maxDistance) {
                if(debug) {
                    cout << "Too far." << endl;
                }
                continue;
            }

            // Queue it.
            q.insert(make_pair(labs(infoA1.offset), segmentId1));
            if(debug) {
                cout << "Queued " << segmentId1 << endl;
            }

        }

    }



    // If getting here, we did not find a segment that satisfies
    // the requested criteria.
    return invalid<uint64_t>;
}



// BFS with given begin/end.
// Does a BFS which starts at segmentIdA.
// and ends when segmentIdB is encountered.
// The BFS if forward if direction is 0
// and backward if direction is 1.
// Computes a vector of all the segments encountered,
// excluding segmentIdA and segmentIdB,
// in the order in which they are encountered in the BFS.
void AssemblyGraph::targetedBfs(
    uint64_t segmentIdA,
    uint64_t segmentIdB,
    uint64_t direction,
    vector<uint64_t>& segments
    ) const
{

    // Initialize the BFS.
    std::queue<uint64_t> q;
    q.push(segmentIdA);

    // Keep track of segments we already encountered.
    std::set<uint64_t> segmentSet;
    segmentSet.insert(segmentIdA);



    // BFS loop.
    segments.clear();
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();

        // Loop over outgoing or incoming links.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        for(const uint64_t linkId: linkIds) {
            const Link& link = links[linkId];

            // Get the segment at the other side of this link.
            const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;

            // If we found segmentIdB, we are done.
            if(segmentId1 == segmentIdB) {
                break;
            }

            // If we already found it, skip it.
            if(segmentSet.contains(segmentId1)) {
                continue;
            }

            // Queue and store this segment.
            q.push(segmentId1);
            segments.push_back(segmentId1);
            segmentSet.insert(segmentId1);
        }
    }

}


// Assemble the assembly paths stored in the JaccardGraph.
void AssemblyGraph::assembleJaccardGraphPaths()
{
    const JaccardGraph& jaccardGraph = *jaccardGraphPointer;
    ofstream fasta("JaccardGraphPaths.fasta");

    uint64_t totalSequenceAssembled = 0;
    for(uint64_t clusterId=0; clusterId<jaccardGraph.assemblyPaths.size(); clusterId++) {
        const vector<uint64_t>& primarySegments = jaccardGraph.assemblyPaths[clusterId];
        AssemblyPath assemblyPath;
        assembleJaccardGraphPath(primarySegments, assemblyPath);

        const auto& sequence = assemblyPath.rawSequence;
        totalSequenceAssembled += sequence.size();
        fasta << ">" << clusterId << " " << sequence.size() << "\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
        fasta << "\n";
    }
    cout << "Assembled a total " << totalSequenceAssembled << " bases." << endl;
}



void AssemblyGraph::assembleJaccardGraphPath(
    const vector<uint64_t>& primarySegments,
    AssemblyPath& assemblyPath)
{
    SHASTA_ASSERT(primarySegments.size() >= 2);

    const JaccardGraph& jaccardGraph = *jaccardGraphPointer;

    // Initialize the path and add the first primary segment.
    assemblyPath.segments.clear();
    assemblyPath.links.clear();
    assemblyPath.segments.push_back(AssemblyPathSegment(primarySegments.front(), true));

    // Add the remaining primary and secondary segments.
    for(uint64_t i=1; i<primarySegments.size(); i++) {
        const uint64_t primarySegment0 = primarySegments[i-1];
        const uint64_t primarySegment1 = primarySegments[i];

        // Get the JaccardGraph vertices corresponding to these primary segments.
        const JaccardGraph::vertex_descriptor v0 = jaccardGraph.vertexTable[primarySegment0];
        const JaccardGraph::vertex_descriptor v1 = jaccardGraph.vertexTable[primarySegment1];

        // Get the JaccardGraph edge between these two vertices.
        JaccardGraph::edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(v0, v1, jaccardGraph);
        SHASTA_ASSERT(edgeWasFound);
        const JaccardGraphEdge& edge = jaccardGraph[e];

        // Access the secondary segments on this edge.
        const vector<uint64_t>& secondarySegmentIds = edge.segmentIds;

        // Add the secondary segments between these two primary segments.
        for(const uint64_t segmentId: secondarySegmentIds) {
            AssemblyPathSegment assemblyPathSegment(segmentId, false);
            assemblyPathSegment.previousPrimarySegmentId = primarySegment0;
            assemblyPathSegment.nextPrimarySegmentId = primarySegment1;
            assemblyPath.segments.push_back(assemblyPathSegment);

        }

        // Add the next primary segment.
        assemblyPath.segments.push_back(AssemblyPathSegment(primarySegment1, true));
    }

    // Assemble sequence for this path.
    assemblyPath.assemble(*this, reads);

}



// De Bruijn graph of the assembly graph journeys of all oriented reads.
// Each assembly graph journey is interpreted as a sequence of segment ids.
// Each vertex represents a sequence of K segment ids.
template<uint64_t K> void AssemblyGraph::createDeBruijnGraphTemplated() const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 8;

    // Type used to store the sequence of a vertex (K segment ids).
    using VertexSequence = array<uint64_t, K>;

    // Loop over all oriented reads to gather vertex sequences.
    vector<VertexSequence> vertexSequences;
    for(uint64_t i=0; i<assemblyGraphJourneys.size(); i++) {

        // Get the assembly graph journey for this oriented read.
        const span<const AssemblyGraphJourneyEntry>& journey = assemblyGraphJourneys[i];
        const uint64_t journeyLength = journey.size();

        // If too short, skip.
        if(journeyLength < K) {
            continue;
        }

        // Extract sequences of length K from the journey.
        // Loop over starting positions.
        for(uint64_t j=0; j<=journeyLength-K; j++) {

            // Fill in the seqyence of length K starting here.
            VertexSequence vertexSequence;
            for(uint64_t k=0; k<K; k++) {
                vertexSequence[k] = journey[j+k].segmentId;
            }
            // Store it.
            vertexSequences.push_back(vertexSequence);
        }
    }

    // Count how many times each sequence was found.
    vector<uint64_t> coverage;
    deduplicateAndCount(vertexSequences, coverage);
    SHASTA_ASSERT(vertexSequences.size() == coverage.size());

    // Each sequence with sufficient coverage generates a vertex.
    using Vertex = pair<VertexSequence, uint64_t>;  // Stores sequence and coverage.
    using Graph = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, Vertex>;
    using vertex_descriptor = Graph::vertex_descriptor;
    Graph graph;

    for(uint64_t i=0; i<vertexSequences.size(); i++) {
        const uint64_t c = coverage[i];
        if(c >= minCoverage) {
            add_vertex(Vertex(vertexSequences[i], c), graph);
        }
    }

    // To generate edges, index the vertices by their (K-1)-prefix.
    using Prefix = array<uint64_t, K-1>;
    std::map<Prefix, vector<vertex_descriptor> > vertexMap;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const VertexSequence& sequence = graph[v].first;
        Prefix prefix;
        copy(sequence.begin(), sequence.begin()+K-1, prefix.begin());
        vertexMap[prefix].push_back(v);
    }

    // Now we can generate the edges.
    using Suffix = Prefix;
    BGL_FORALL_VERTICES_T(v0, graph, Graph) {
        const VertexSequence& sequence = graph[v0].first;
        Suffix suffix;
        copy(sequence.begin()+1, sequence.end(), suffix.begin());
        auto it = vertexMap.find(suffix);
        if(it == vertexMap.end()) {
            continue;
        }
        for(const vertex_descriptor v1: it->second) {
            add_edge(v0, v1, graph);
        }
    }
    cout << "The DeBruijn graph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;

    // Write it out.
    ofstream dot("DeBruijnGraph.dot");
    dot << "digraph DeBruijnGraph {\n";
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const uint64_t coverage = graph[v].second;
        dot << "\"" << v << "\" ["
            "tooltip=\"" << coverage << "\" "
            "width=" << 0.001*double(coverage) <<
            "];\n";
    }
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dot << "\"" << v0 << "\"->\"" << v1 << "\";\n";
    }
    dot << "}\n";
}



void AssemblyGraph::createDeBruijnGraph() const
{
    createDeBruijnGraphTemplated<3>();
}



void AssemblyGraph::assembleSegments()
{
    createNew(segmentSequences, "Mode3-SegmentSequences");
    createNew(segmentVertexOffsets, "Mode3-SegmentVertexOffsets");
    for(uint64_t segmentId=0; segmentId<markerGraphPaths.size(); segmentId++) {
        assembleSegment(segmentId);
    }
}
void AssemblyGraph::assembleSegment(uint64_t segmentId)
{
    // Assemble it.
    AssembledSegment assembledSegment;
    assembleMarkerGraphPath(
        readRepresentation,
        k,
        reads,
        markers,
        markerGraph,
        markerGraphPaths[segmentId],
        false,
        assembledSegment);

    // Store assembled sequence and vertex offsets.
    segmentSequences.appendVector(assembledSegment.rawSequence);
    segmentVertexOffsets.appendVector(assembledSegment.vertexOffsets);

}



// Assemble a link, given a set of allowed OrientedReadId(s).
// The returned sequence overrides:
// - The trim0 last bases of the preceding segment.
// - The trim1 first bases of the following segment.
void AssemblyGraph::assembleLink(
    uint64_t linkId,
    const vector<OrientedReadId>& allowedOrientedReadIds,
    vector<Base>& sequence, // The entire MSA sequence
    uint64_t& leftTrim,     // The number of MSA sequence to be trimmed on the left for assembly
    uint64_t& rightTrim,    // The number of MSA sequence to be trimmed on the left for assembly
    uint64_t& trim0,        // The number of bases at the end of segment0 to be trimmed for assembly
    uint64_t& trim1,        // The number of bases at the beginning of segment1 to be trimmed for assembly
    ostream& html
) const
{
    const bool debug = false;
    const Link& link = links[linkId];

    SHASTA_ASSERT(std::is_sorted(allowedOrientedReadIds.begin(), allowedOrientedReadIds.end()));

    // If the preceding and last segment are adjacent,
    // assembling the link is trivial.
    // We return an empty sequence which overrides
    // the last k/2 bases of the preceding segment
    // and the first k/2 bases of the following segment.
    // This was the resulting segment sequences are exactly adjacent.
    if(link.segmentsAreAdjacent) {
        sequence.clear();
        leftTrim = 0;
        rightTrim = 0;
        trim0 = k / 2;
        trim1 = k / 2;
        return;
    }

    // If getting here, the two segments of this link are not adjacent.
    // Get some infomation we are going to need below.
    const uint64_t segmentId0 = link.segmentId0;
    const uint64_t segmentId1 = link.segmentId1;
    const auto sequence0 = segmentSequences[segmentId0];
    const auto sequence1 = segmentSequences[segmentId1];
    SHASTA_ASSERT(not sequence0.empty());
    SHASTA_ASSERT(not sequence1.empty());
    const auto vertexOffsets0 = segmentVertexOffsets[segmentId0];
    const auto vertexOffsets1 = segmentVertexOffsets[segmentId1];


    // First, find:
    // - The position in segmentId0 of the leftmost transition.
    // - The position in segmentId1 of the rightmost transition.
    uint64_t minEdgePosition0 = markerGraphPaths[segmentId0].size();
    uint64_t maxEdgePosition1 = 0;
    for(const auto &p : transitions[linkId]) {
        const OrientedReadId orientedReadId = p.first;

        // If not one of the allowed OrientedReadId(s), skip it.
        if(not binary_search(allowedOrientedReadIds.begin(), allowedOrientedReadIds.end(), orientedReadId)) {
            continue;
        }

        // Access the transition from segmentId0 to segmentId1 for this oriented read.
        const Transition &transition = p.second;

        minEdgePosition0 = min(minEdgePosition0,
            uint64_t(transition[0].position));
        maxEdgePosition1 = max(maxEdgePosition1,
            uint64_t(transition[1].position));
    }

    // When getting here:
    // - minEdgePosition0 is the leftmost position of the transitions in path0.
    // - maxEdgePosition1 is the rightmost position of the transitions in path1.
    // These positions are edge positions in markerGraphPath0 and markerGraphPath1.
    // We will do a multiple sequence alignment of the oriented reads,
    // using the sequence of segmentId0 to extend to the left all reads to minEdgePosition0,
    // and using the sequence of segmentId1 to extend to the right all reads to maxEdgePosition1,

    // Get the corresponding vertex positions in segmentId0 and segmentId1.
    const uint64_t minVertexPosition0 = minEdgePosition0 + 1;
    const uint64_t maxVertexPosition1 = maxEdgePosition1;

    // To compute an MSA anchored at both sides,we will extend the
    // sequence of each read to the left/right using the sequence of
    // adjacent segments.


   // Now extract the portion of each oriented read sequence that
    // will be used to assemble this link.
    vector<OrientedReadId> orientedReadIdsForAssembly;
    vector<vector<Base> > orientedReadsSequencesForAssembly;
    for(const auto &p : transitions[linkId]) {
        const OrientedReadId orientedReadId = p.first;

        // If not one of the allowed OrientedReadId(s), skip it.
        if(not binary_search(allowedOrientedReadIds.begin(), allowedOrientedReadIds.end(), orientedReadId)) {
            continue;
        }

        // Access the transition from segmentId0 to segmentId1 for this oriented read.
        const Transition &transition = p.second;

        // Get the ordinals of the last appearance of this oriented
        // read on segmentId0 and the first on segmentId1,
        // and the corresponding markers.
        const uint32_t ordinal0 = transition[0].ordinals[1];
        const uint32_t ordinal1 = transition[1].ordinals[0];
        const CompressedMarker &marker0 = markers[orientedReadId.getValue()][ordinal0];
        const CompressedMarker &marker1 = markers[orientedReadId.getValue()][ordinal1];

        // Get the positions of these markers on the oriented read.
        // If using RLE, these are RLE positions.
        const uint32_t position0 = marker0.position;
        const uint32_t position1 = marker1.position;

        // Extract the sequence between these markers (including the markers).
        vector<Base> orientedReadSequence;
        for(uint64_t position = position0;
            position < position1 + k; position++) {
            const Base b = reads.getOrientedReadBase(orientedReadId, uint32_t(position));
            orientedReadSequence.push_back(b);
        }

        // We need to extend the sequence of this read to the left,
        // using segmentId0 sequence, up to minVertexPosition0,
        // so the portions of all reads we will be using for the MSA
        // all begin in the same place.
        vector<Base> leftSequence;
        vector<uint32_t> leftRepeatCounts;
        const uint64_t vertexPosition0 = transition[0].position + 1; // Add 1 to get vertex position.
        const uint64_t begin0 =
            vertexOffsets0[minVertexPosition0];
        const uint64_t end0 = vertexOffsets0[vertexPosition0];
        for(uint64_t position = begin0; position != end0; position++) {
            leftSequence.push_back(sequence0[position]);
        }

        vector<Base> rightSequence;
        const uint64_t vertexPosition1 = transition[1].position;
        const uint64_t begin1 = vertexOffsets1[vertexPosition1] + k;
        const uint64_t end1 = vertexOffsets1[maxVertexPosition1] + k;
        for(uint64_t position = begin1; position != end1; position++) {
            rightSequence.push_back(sequence1[position]);
       }

        // Construct the extended sequence for this oriented read,
        // to be used in the MSA.
        vector<Base> orientedReadExtendedSequence;
        const auto addToExtendedSequence = back_inserter(orientedReadExtendedSequence);
        copy(leftSequence, addToExtendedSequence);
        copy(orientedReadSequence, addToExtendedSequence);
        copy(rightSequence, addToExtendedSequence);

        orientedReadIdsForAssembly.push_back(orientedReadId);
        orientedReadsSequencesForAssembly.push_back(orientedReadExtendedSequence);

        if(debug) {
            copy(orientedReadExtendedSequence, ostream_iterator<Base>(cout));
            cout << " " << orientedReadId << endl;
        }
    }



   // Compute the consensus sequence for the link.
   if(html) {
       html << "<h2>Link " << linkId << "</h2>\n";
   }
   vector<Base> msaRleSequence;
   computeLinkConsensusUsingSpoa(
       orientedReadIdsForAssembly,
       orientedReadsSequencesForAssembly,
       consensusCaller,
       debug,
       html,
       msaRleSequence);

   if(debug) {
       cout << "Consensus RLE sequence length before trimming " << msaRleSequence.size() << endl;
       cout << "Portion of segment on left involved in the MSA begins at position " <<
           vertexOffsets0[minVertexPosition0] << endl;
       cout << "Portion of segment on right involved in the MSA ends at position " <<
           vertexOffsets1[maxVertexPosition1] + k << endl;
   }

   // Count the number of identical (RLE) bases at the beginning of the
   // link consensus sequence and of the segmentId0 sequence portion
   // involved in assembling this link.
   uint64_t identicalOnLeft = 0;
   const uint64_t begin0 = vertexOffsets0[minVertexPosition0];
   const uint64_t end0 = sequence0.size();
   for(uint64_t i=begin0; (i!=end0 and (i-begin0)<msaRleSequence.size()); i++) {
       if(msaRleSequence[i-begin0] == sequence0[i]) {
           // cout << "*** " << begin0 << " " << end0 << " " << i << endl;
           ++identicalOnLeft;
       } else {
           break;
       }
   }
   if(debug) {
       cout << "Identical on left: " << identicalOnLeft << endl;
   }

   // Count the number of identical (RLE) bases at the end of the
   // link consensus sequence and the beginning of segmentId1 .
   uint64_t identicalOnRight = 0;
   const uint64_t end1 = vertexOffsets1[maxVertexPosition1] + k;
   for(uint64_t i=end1-1; ; i--) {
       const uint64_t j = msaRleSequence.size() - (end1 - i);
       if(msaRleSequence[j] == sequence1[i]) {
           // cout << "*** " << i << " " << assembledSegment1.runLengthSequence[i] << " " <<
           //     j << " " << consensusRleSequence[j] << endl;
           ++identicalOnRight;
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
   identicalOnRight = min(identicalOnRight, msaRleSequence.size()-identicalOnLeft);
   if(debug) {
       cout << "Identical on right: " << identicalOnRight << endl;
   }

   // Trim these identical bases from the link consensus sequence.
   leftTrim = identicalOnLeft;
   rightTrim = identicalOnRight;

   // Compute and store the number of bases to be trimmed at the end of segmentId0
   // and at the beginning of segmentId1.
   trim0 =
       sequence0.size() -
       vertexOffsets0[minVertexPosition0] -
       identicalOnLeft;
   trim1 =
       vertexOffsets1[maxVertexPosition1] + k
       - identicalOnRight;
}



void AssemblyGraph::computeLinkConsensusUsingSpoa(
    const vector<OrientedReadId> orientedReadIds,
    const vector< vector<Base> > rleSequences,
    const ConsensusCaller&,
    bool debug,
    ostream& html,
    vector<Base>& consensusRleSequence
    )
{
    SHASTA_ASSERT(0);
}

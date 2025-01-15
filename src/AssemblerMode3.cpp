// Shasta
#include "Assembler.hpp"
#include "LocalMarkerGraph1.hpp"
#include "mode3-LocalAssembly.hpp"
#include "mode3-AssemblyGraph.hpp"
#include "mode3-AnchorGraph.hpp"
#include "Mode3Assembler.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>



// If the coverage range for primary marker graph edges is not
// specified, this uses the disjoint sets histogram to compute reasonable values.
pair<uint64_t, uint64_t> Assembler::getPrimaryCoverageRange()
{
    const auto& disjointSetsHistogram = markerGraph.disjointSetsHistogram;

    // Set minPrimaryCoverage to the first value where the
    // disjointSetsHistogram starts increasing.
    uint64_t minPrimaryCoverage = invalid<uint64_t>;
    uint64_t frequencyAtMinPrimaryCoverage = 0;
    for(uint64_t i=1; i<disjointSetsHistogram.size(); i++) {
        const uint64_t coverage = disjointSetsHistogram[i].first;
        const uint64_t frequency = disjointSetsHistogram[i].second;
        const uint64_t previousCoverage = disjointSetsHistogram[i-1].first;
        const uint64_t previousFrequency = disjointSetsHistogram[i-1].second;
        if(
            (coverage != previousCoverage+1) // Frequency at coverage-1 is zero, so the histogram went up.
            or
            frequency > previousFrequency    // The histogram went up.
            ) {
            minPrimaryCoverage = coverage;
            frequencyAtMinPrimaryCoverage = frequency;
            break;
        }
    }
    SHASTA_ASSERT(minPrimaryCoverage != invalid<uint64_t>);

    // Set maxPrimaryCoverage to the last coverage with frequency
    // at least equal to frequencyAtMinPrimaryCoverage.
    uint64_t maxPrimaryCoverage = invalid<uint64_t>;
    for(uint64_t i=disjointSetsHistogram.size()-1; i>0; i--) {
        const uint64_t coverage = disjointSetsHistogram[i].first;
        const uint64_t frequency = disjointSetsHistogram[i].second;
        if(frequency >= frequencyAtMinPrimaryCoverage) {
            maxPrimaryCoverage = coverage;
            break;
        }
    }
    SHASTA_ASSERT(maxPrimaryCoverage != invalid<uint64_t>);

    return {minPrimaryCoverage, maxPrimaryCoverage};
}



void Assembler::mode3Assembly(
    uint64_t threadCount,
    shared_ptr<mode3::Anchors> anchorsPointer,
    const Mode3AssemblyOptions& options,
    bool debug
    )
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    mode3Assembler = make_shared<Mode3Assembler>(
        mappedMemoryOwner,
        assemblerInfo->k, getReads(), markers,
        anchorsPointer, threadCount, options, debug);
}



// Same, but use existing Anchors. Python callable.
void Assembler::mode3Reassembly(
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool debug
    )
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    // Create the Anchors from binary data.
    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(mappedMemoryOwner, getReads(), assemblerInfo->k, markers);

    // Run the Mode 3 assembly.
    mode3Assembler = make_shared<Mode3Assembler>(
        mappedMemoryOwner,
        assemblerInfo->k, getReads(), markers,
        anchorsPointer, threadCount, options, debug);
}

void Assembler::accessMode3Assembler()
{
    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(MappedMemoryOwner(*this), getReads(), assemblerInfo->k, markers);
    mode3Assembler = make_shared<Mode3Assembler>(*this,
        assemblerInfo->k, getReads(), markers,
        anchorsPointer, httpServerData.assemblerOptions->assemblyOptions.mode3Options);
}



// Assemble sequence between two Anchors.
void Assembler::fillMode3AssemblyPathStep(const vector<string>& request, ostream& html)
{
    // Check that our assumptions are satisfied.
    if(assemblerInfo->assemblyMode != 3) {
        throw runtime_error("This is only available for assembly mode 3.");
    }
    SHASTA_ASSERT(getReads().representation == 0);      // No RLE.
    SHASTA_ASSERT((assemblerInfo->k % 2) == 0);         // Marker length is even.

    mode3::LocalAssemblyDisplayOptions options(html);

    // Get the parameters for the request.
    uint64_t edgeIdA = invalid<uint64_t>;
    getParameterValue(request, "edgeIdA", edgeIdA);

    uint64_t edgeIdB = invalid<uint64_t>;
    getParameterValue(request, "edgeIdB", edgeIdB);

    string useAString;
    const bool useA = getParameterValue(request, "useA", useAString);

    string useBString;
    const bool useB = getParameterValue(request, "useB", useBString);

    uint64_t minVertexCoverage = 0;
    getParameterValue(request, "minVertexCoverage", minVertexCoverage);

    string showOrientedReadsString;
    options.showOrientedReads = getParameterValue(request, "showOrientedReads", showOrientedReadsString);

    string showMarkersString;
    options.showMarkers = getParameterValue(request, "showMarkers", showMarkersString);

    string showGraphString;
    options.showGraph = getParameterValue(request, "showGraph", showGraphString);

    string showVerticesString;
    options.showVertices = getParameterValue(request, "showVertices", showVerticesString);

    string showVertexLabelsString;
    options.showVertexLabels = getParameterValue(request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    options.showEdgeLabels = getParameterValue(request, "showEdgeLabels", showEdgeLabelsString);

    string showAssemblyDetailsString;
    options.showAssemblyDetails = getParameterValue(request, "showAssemblyDetails", showAssemblyDetailsString);

    string showDebugInformationString;
    options.showDebugInformation = getParameterValue(request, "showDebugInformation", showDebugInformationString);



    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr><th class=left>Edge A<td class=centered>"
        "<input type=text required name=edgeIdA size=8 style='text-align:center' " <<
        ((edgeIdA == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeIdA) + "'")) << ">"

        "<tr><th class=left>Edge B<td class=centered>"
        "<input type=text required name=edgeIdB size=8 style='text-align:center' " <<
        ((edgeIdB == invalid<uint64_t>) ? "" : ("value='" + to_string(edgeIdB) + "'")) << ">"

        "<tr>"
        "<th class=left>Use for assembly oriented reads that appear only on edge A"
        "<td class=centered><input type=checkbox name=useA" <<
        (useA ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Use for assembly oriented reads that appear only on edge B"
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



    // If the edge ids are missing, do nothing.
    if(edgeIdA == invalid<uint64_t> or edgeIdB == invalid<uint64_t>) {
        return;
    }

    // Sanity checks on the edge ids.
    if(edgeIdA >= markerGraph.edgeMarkerIntervals.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdA) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edgeMarkerIntervals.size()));
    }
    if(edgeIdB >= markerGraph.edgeMarkerIntervals.size()) {
        throw runtime_error("Marker graph edge " + to_string(edgeIdB) +
            " is not valid. Maximum valid edge id is " + to_string(markerGraph.edgeMarkerIntervals.size()));
    }

    // Sanity check that the two edges are distinct.
    if(edgeIdA == edgeIdB) {
        html << "<p>Specify two distinct edges.";
        return;
    }

    // This analysis can only be done if both edges have no duplicate OrientedReadIds
    // in their MarkerIntervals.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA)) {
        html << "<p>Marker graph edge " << edgeIdA << " has duplicate oriented reads.";
        return;
    }
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB)) {
        html << "<p>Marker graph edge " << edgeIdB << " has duplicate oriented reads.";
        return;
    }

    // Check that there are common reads.
    MarkerGraphEdgePairInfo info;
    SHASTA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
    if(info.common == 0) {
        html << "<p>The two edges have no common oriented reads.";
        return;
    }

    // Local assembly for this assembly step.
    mode3::LocalAssembly localAssembly(
        assemblerInfo->k, getReads(), markers, mode3Assembler->anchors(),
        edgeIdA, edgeIdB, minVertexCoverage,
        options,
        httpServerData.assemblerOptions->assemblyOptions.mode3Options.localAssemblyOptions,
        useA, useB);
}



// Creation of marker graph edges (anchors) for Mode 3 assembly.
// - Like createMarkerGraphEdgesStrict, this
//   will only create edges in which all contributing oriented reads have
//   exactly the same sequence. If more than one distinct sequence
//   is present, the edge is split into two parallel edges.
// - This only generates primary marker graph edges, defined as follows:
//   * Edge coverage is >= minPrimaryCoverage and <= maxPrimaryCoverage.
//   * Both its vertices have no duplicate ReadId.
//   * The edge marker interval has no duplicate ReadId.
// - This will only create the MarkerGraph::edgeMarkerIntervals
//   and nothing else.
void Assembler::createPrimaryMarkerGraphEdges(
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t threadCount)
{
    performanceLog << timestamp << "createPrimaryMarkerGraphEdges begins." << endl;

    // Check that we have what we need.
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the arguments so the threads can see them.
    createPrimaryMarkerGraphEdgesData.minPrimaryCoverage = minPrimaryCoverage;
    createPrimaryMarkerGraphEdgesData.maxPrimaryCoverage = maxPrimaryCoverage;

    // Make space for the edges found by each thread.
    createPrimaryMarkerGraphEdgesData.threadMarkerIntervals.resize(threadCount);
    createPrimaryMarkerGraphEdgesData.threadSequences.resize(threadCount);

    // Parallelize over the source vertex of primary marker graph edges.
    // Each thread stores in a separate vector the marker intervals for the primary edges it finds.
    uint64_t batchSize = 1000;
    const MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::CompressedVertexId>& markerGraphVertices =
        *markerGraph.verticesPointer;
    setupLoadBalancing(markerGraphVertices.size(), batchSize);
    runThreads(&Assembler::createPrimaryMarkerGraphEdgesThreadFunction, threadCount);

    // Gather the primary edges found by all threads.
    markerGraph.edgeMarkerIntervals.createNew(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"),
            largeDataPageSize);
    markerGraph.edgeSequence.createNew(
        largeDataName("MarkerGraphEdgesSequence"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadMarkerIntervals = *createPrimaryMarkerGraphEdgesData.threadMarkerIntervals[threadId];
        auto& threadSequences = *createPrimaryMarkerGraphEdgesData.threadSequences[threadId];
        SHASTA_ASSERT(threadMarkerIntervals.size() == threadSequences.size());
        for(uint64_t i=0; i<threadMarkerIntervals.size(); i++) {
            const span<MarkerInterval> markerIntervals = threadMarkerIntervals[i];
            markerGraph.edgeMarkerIntervals.appendVector(markerIntervals.begin(), markerIntervals.end());
            const span<Base> sequence = threadSequences[i];
            markerGraph.edgeSequence.appendVector(sequence.begin(), sequence.end());
        }
        threadMarkerIntervals.remove();
        threadSequences.remove();
        createPrimaryMarkerGraphEdgesData.threadMarkerIntervals[threadId] = 0;
        createPrimaryMarkerGraphEdgesData.threadSequences[threadId] = 0;
    }
    createPrimaryMarkerGraphEdgesData.threadMarkerIntervals.clear();
    createPrimaryMarkerGraphEdgesData.threadSequences.clear();
    cout << "Found " << markerGraph.edgeMarkerIntervals.size() << " primary marker graph edges." << endl;

    performanceLog << timestamp << "createPrimaryMarkerGraphEdges ends." << endl;
}



void Assembler::createPrimaryMarkerGraphEdgesThreadFunction(uint64_t threadId)
{
    // Sanity check and get k/2.
    const uint64_t k = assemblerInfo->k;
    SHASTA_ASSERT((k % 2) == 0);
    const uint32_t kHalf = uint32_t(k / 2);
    SHASTA_ASSERT(getReads().representation == 0);

    // Access the data set up by createPrimaryMarkerGraphEdges.
    auto& data = createPrimaryMarkerGraphEdgesData;

    // Get the primary coverage range.
    const uint64_t minPrimaryCoverage = data.minPrimaryCoverage;
    const uint64_t maxPrimaryCoverage = data.maxPrimaryCoverage;

    // Create the vector to contain the marker intervals for the edges found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > markerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> >();
    data.threadMarkerIntervals[threadId] = markerIntervalsPointer;
    MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t>&
        markerIntervals = *markerIntervalsPointer;
    markerIntervals.createNew(
            largeDataName("tmp-ThreadMarkerGraphEdgeMarkerIntervals-" + to_string(threadId)),
            largeDataPageSize);

    // Create the vector to contain the sequences for the edges found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > sequencesPointer =
        make_shared< MemoryMapped::VectorOfVectors<Base, uint64_t> >();
    data.threadSequences[threadId] = sequencesPointer;
    MemoryMapped::VectorOfVectors<Base, uint64_t>&
        sequences = *sequencesPointer;
    sequences.createNew(
            largeDataName("tmp-ThreadMarkerGraphSequences-" + to_string(threadId)),
            largeDataPageSize);


    // Access data structures we need.
    const Reads& reads = getReads();
    const MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::CompressedVertexId>& markerGraphVertices =
        markerGraph.vertices();

    // For each sequence and next vertex we encounter, we store the corresponding MarkerIntervals.
    class Info {
    public:
        vector<Base> sequence;
        MarkerGraphVertexId vertexId1;
        vector<MarkerInterval> markerIntervals;
        Info(
            const vector<Base>& sequence,
            MarkerGraphVertexId vertexId1,
            const MarkerInterval& markerInterval) :
            sequence(sequence),
            vertexId1(vertexId1),
            markerIntervals(1, markerInterval) {}
    };
    vector<Info> infos;

    // Work vector used in the loop below.
    vector<Base> sequence;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over source vertex ids assigned to this batch.
        for(MarkerGraphVertexId vertexId0=begin; vertexId0!=end; vertexId0++) {

            infos.clear();

            // Access the markers of this vertex.
            const auto vertexMarkerIds = markerGraphVertices[vertexId0];

            // If vertex coverage of this vertex is less than minPrimaryCoverage,
            // no primary edges start here for sure.
            if(vertexMarkerIds.size() < minPrimaryCoverage) {
                continue;
            }

            // If this vertex has duplicate ReadIds, no primary edge can begin here.
            if(markerGraph.vertexHasDuplicateReadIds(vertexId0, markers)) {
                continue;
            }

            // Loop over the markers of this vertex.
            for(const MarkerId markerId0: vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal0;
                tie(orientedReadId, ordinal0) = findMarkerId(markerId0);

                // Find the next marker graph vertex visited by this OrientedReadId.
                const uint32_t ordinal1 = ordinal0 + 1;
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                if(ordinal1 == orientedReadMarkers.size()) {
                    // There is no next vertex.
                    continue;
                }
                SHASTA_ASSERT(ordinal1 < orientedReadMarkers.size());
                const MarkerId markerId1 = markerId0 + 1;
                const MarkerGraphVertexId vertexId1 = markerGraph.vertexTable[markerId1];

                // If vertexId1 has coverage less than minPrimaryCoverage, there
                // cannot be a primary edge vertexId0->vertexId1.
                if(markerGraphVertices.size(vertexId1) < minPrimaryCoverage) {
                    continue;
                }

                // Get the sequence between these two markers.
                const uint32_t position0 = uint32_t(markers.begin()[markerId0].position) + kHalf;
                const uint32_t position1 = uint32_t(markers.begin()[markerId1].position) + kHalf;
                sequence.clear();
                for(uint32_t position=position0; position<position1; position++) {
                    const Base base = reads.getOrientedReadBase(orientedReadId, position);
                    sequence.push_back(base);
                }

                // Construct this MarkerInterval.
                MarkerInterval markerInterval;
                markerInterval.orientedReadId = orientedReadId;
                markerInterval.ordinals[0] = ordinal0;
                markerInterval.ordinals[1] = ordinal1;

                // If we already have an Info with this sequence and vertexId1, add id there.
                // Otherwise create a new Info.
                bool found = false;
                for(Info& info: infos) {
                    if(info.sequence == sequence and info.vertexId1 == vertexId1) {
                        info.markerIntervals.push_back(markerInterval);
                        found = true;
                        break;
                    }
                }
                if(not found) {
                    infos.push_back(Info(sequence, vertexId1, markerInterval));
                }
            }

            // Each of our Infos object can generate a primary edge, but we have to check that:
            // - Coverage is in the primary coverage range.
            // - vertexId1 does not have duplicate ReadIds.
            for(const Info& info: infos) {
                if(info.markerIntervals.size() < minPrimaryCoverage) {
                    continue;
                }
                if(info.markerIntervals.size() > maxPrimaryCoverage) {
                    continue;
                }
                if(markerGraph.vertexHasDuplicateReadIds(info.vertexId1, markers)) {
                    continue;
                }

                // Generate a new primary marker graph edge.
                markerIntervals.appendVector(info.markerIntervals);
                sequences.appendVector(info.sequence);
            }

        }
    }

}



void Assembler::exploreAnchor(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreAnchor(request, html);
}



void Assembler::exploreAnchorPair(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreAnchorPair(request, html);
}



void Assembler::exploreJourney(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 0);
    SHASTA_ASSERT(assemblerInfo->assemblyMode == 3);
    mode3Assembler->exploreJourney(request, html);
}



void Assembler::exploreReadFollowing(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 0);
    SHASTA_ASSERT(assemblerInfo->assemblyMode == 3);
    mode3Assembler->exploreReadFollowing(request, html);
}



void Assembler::exploreLocalAssembly(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreLocalAssembly(request, html);
}



void Assembler::exploreLocalAnchorGraph(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreLocalAnchorGraph(request, html);
}



void Assembler::exploreMode3AssemblyGraph(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreAssemblyGraph(request, html);
}



void Assembler::exploreSegment(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreSegment(request, html);
}



void Assembler::exploreReadFollowingAssemblyGraph(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 0);
    SHASTA_ASSERT(assemblerInfo->assemblyMode == 3);
    mode3Assembler->exploreReadFollowingAssemblyGraph(request, html);
}



// Alignment-free version of mode 3 assembly.
void Assembler::alignmentFreeAssembly(
    const Mode3AssemblyOptions& mode3Options,
    const vector<string>& anchorFileAbsolutePaths,
    uint64_t threadCount)
{
    cout << timestamp << "Alignment free mode 3 assembly begins." << endl;

    // Create the Anchors.
    shared_ptr<mode3::Anchors> anchorsPointer;
    if(mode3Options.anchorCreationMethod == "FromMarkerKmers") {
        anchorsPointer =
            make_shared<mode3::Anchors>(
                MappedMemoryOwner(*this),
                getReads(),
                assemblerInfo->k,
                markers,
                mode3Options.minAnchorCoverage,
                mode3Options.maxAnchorCoverage,
                threadCount);
    } else if(mode3Options.anchorCreationMethod == "FromJson") {
        anchorsPointer =
            make_shared<mode3::Anchors>(
                MappedMemoryOwner(*this),
                getReads(),
                assemblerInfo->k,
                markers,
                anchorFileAbsolutePaths,
                mode3Options.minAnchorCoverage,
                mode3Options.maxAnchorCoverage,
                threadCount);
    } else {
        throw runtime_error("Invalid value for --Assembly.mode3.anchorCreationMethod.");
    }

    // Compute oriented read journeys.
    anchorsPointer->computeJourneys(threadCount);

    // Run Mode 3 assembly.
    mode3Assembly(threadCount, anchorsPointer, mode3Options, false);

    cout << timestamp << "Alignment free mode 3 assembly ends." << endl;
}

/********************************************************************************

Creation of Anchors from the MarkerGraph.

Each Anchor corresponds to a "primary marker graph edge"
defined as follows:

- All contributing oriented reads have exactly the same sequence.
  If more than one distinct sequence is present, the edge (anchor)
  is split into two edges (anchors).
- Edge coverage is >= minPrimaryCoverage and <= maxPrimaryCoverage.
- Both vertices have no duplicate ReadIds, and as a result the
  resulting anchor has no duplicate ReadIds.

This uses as input the MarkerGraph vertices. MarkerGraph edges
are created implicitly and never stored. Instead, the corresponding
information is stored directly in the Anchors.

This code is equivalent to Assembler::createPrimaryMarkerGraphEdges.

********************************************************************************/

#include "mode3-Anchor.hpp"
#include "findMarkerId.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;


Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
#if 0
    // For now copy the marker intervals from the marker graph.
    anchorMarkerIntervals.createNew(largeDataName("AnchorMarkerIntervals"), largeDataPageSize);
    for(uint64_t anchorId=0; anchorId<markerGraph.edgeMarkerIntervals.size(); anchorId++) {
        const auto v = markerGraph.edgeMarkerIntervals[anchorId];
        anchorMarkerIntervals.appendVector(v.begin(), v.end());
    }

    // Also copy the Anchor sequences from the marker graph.
    anchorSequences.createNew(largeDataName("AnchorSequences"), largeDataPageSize);
    for(uint64_t anchorId=0; anchorId<markerGraph.edgeSequence.size(); anchorId++) {
        const auto v = markerGraph.edgeSequence[anchorId];
        anchorSequences.appendVector(v.begin(), v.end());
    }

    SHASTA_ASSERT(anchorSequences.size() == anchorMarkerIntervals.size());

    check();
#endif


    performanceLog << timestamp << "Anchor creation from the marker graph begins." << endl;

    // Sanity checks and get kHalf.
    SHASTA_ASSERT((k %2) == 0);
    kHalf = k / 2;
    SHASTA_ASSERT(reads.representation == 0);
    SHASTA_ASSERT(markers.isOpen());
    SHASTA_ASSERT(markerGraph.vertices().isOpen());
    SHASTA_ASSERT(markerGraph.vertexTable.isOpen);

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the arguments so the threads can see them.
    auto& data = constructFromMarkerGraphData;
    data.minPrimaryCoverage = minPrimaryCoverage;
    data.maxPrimaryCoverage = maxPrimaryCoverage;

    data.markerGraphPointer = &markerGraph;

    // Make space for the edges found by each thread.
    constructFromMarkerGraphData.threadMarkerIntervals.resize(threadCount);
    constructFromMarkerGraphData.threadSequences.resize(threadCount);

    // Parallelize over the MarkerGraph source vertex.
    // Each thread stores in a separate vector the Anchors it finds.
    uint64_t batchSize = 1000;
    setupLoadBalancing(markerGraph.vertices().size(), batchSize);
    runThreads(&Anchors::constructFromMarkerGraphThreadFunction, threadCount);

    // Gather the Anchors found by all threads.
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadMarkerIntervals = *data.threadMarkerIntervals[threadId];
        auto& threadSequences = *data.threadSequences[threadId];
        SHASTA_ASSERT(threadMarkerIntervals.size() == threadSequences.size());
        for(uint64_t i=0; i<threadMarkerIntervals.size(); i++) {
            const auto thisAnchorMarkerIntervals = threadMarkerIntervals[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& threadMarkerInterval: thisAnchorMarkerIntervals) {
                anchorMarkerIntervals.append(
                    AnchorMarkerInterval(threadMarkerInterval.orientedReadId, threadMarkerInterval.ordinal0));
            }
            const span<Base> sequence = threadSequences[i];
            anchorSequences.appendVector(sequence.begin(), sequence.end());
        }
        threadMarkerIntervals.remove();
        threadSequences.remove();
        data.threadMarkerIntervals[threadId] = 0;
        data.threadSequences[threadId] = 0;
    }
    data.threadMarkerIntervals.clear();
    data.threadSequences.clear();
    cout << "Found " << anchorMarkerIntervals.size() << " anchors." << endl;

    performanceLog << timestamp << "Anchor creation from the marker graph ends." << endl;
}



void Anchors::constructFromMarkerGraphThreadFunction(uint64_t threadId)
{
    // Access the data set up by createPrimaryMarkerGraphEdges.
    using ThreadMarkerInterval = ConstructFromMarkerGraphData::ThreadMarkerInterval;
    auto& data = constructFromMarkerGraphData;

    // Get the primary coverage range.
    const uint64_t minPrimaryCoverage = data.minPrimaryCoverage;
    const uint64_t maxPrimaryCoverage = data.maxPrimaryCoverage;

    const MarkerGraph& markerGraph = *data.markerGraphPointer;
    const MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::CompressedVertexId>& markerGraphVertices =
        markerGraph.vertices();

    // Create the vector to contain the marker intervals for the Anchors found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> > markerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> >();
    data.threadMarkerIntervals[threadId] = markerIntervalsPointer;
    MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t>&
        markerIntervals = *markerIntervalsPointer;
    markerIntervals.createNew(
            largeDataName("tmp-ThreadAnchorMarkerIntervals-" + to_string(threadId)),
            largeDataPageSize);

    // Create the vector to contain the sequences for the Anchors found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > sequencesPointer =
        make_shared< MemoryMapped::VectorOfVectors<Base, uint64_t> >();
    data.threadSequences[threadId] = sequencesPointer;
    MemoryMapped::VectorOfVectors<Base, uint64_t>& sequences = *sequencesPointer;
    sequences.createNew(
            largeDataName("tmp-ThreadAnchorSequences-" + to_string(threadId)),
            largeDataPageSize);

    // For each sequence and next vertex we encounter, we store the corresponding MarkerIntervals.
    class Info {
    public:
        vector<Base> sequence;
        MarkerGraphVertexId vertexId1;
        vector<ThreadMarkerInterval> markerIntervals;
        Info(
            const vector<Base>& sequence,
            MarkerGraphVertexId vertexId1,
            const ThreadMarkerInterval& markerInterval) :
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
                tie(orientedReadId, ordinal0) = shasta::findMarkerId(markerId0, markers);

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
                const uint32_t position0 = uint32_t(markers.begin()[markerId0].position + kHalf);
                const uint32_t position1 = uint32_t(markers.begin()[markerId1].position + kHalf);
                sequence.clear();
                for(uint32_t position=position0; position<position1; position++) {
                    const Base base = reads.getOrientedReadBase(orientedReadId, position);
                    sequence.push_back(base);
                }

                // Construct this ThreadMarkerInterval.
                ThreadMarkerInterval markerInterval;
                markerInterval.orientedReadId = orientedReadId;
                markerInterval.ordinal0 = ordinal0;

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



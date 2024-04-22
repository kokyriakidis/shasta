// Shasta.
#include "MarkerGraph.hpp"
#include "Coverage.hpp"
#include "deduplicate.hpp"
#include "findMarkerId.hpp"
#include "invalid.hpp"
#include "markerAccessFunctions.hpp"
using namespace shasta;

// Standard library.
#include "fstream.hpp"

#include "MultithreadedObject.tpp"
template class MultithreadedObject<MarkerGraph>;


const MarkerGraph::VertexId MarkerGraph::invalidVertexId = std::numeric_limits<VertexId>::max();
const MarkerGraph::EdgeId MarkerGraph::invalidEdgeId = std::numeric_limits<EdgeId>::max();
const MarkerGraph::CompressedVertexId
    MarkerGraph::invalidCompressedVertexId = std::numeric_limits<uint64_t>::max();



MarkerGraph::MarkerGraph() :
    MultithreadedObject<MarkerGraph>(*this)
    {}



void MarkerGraph::remove()
{
    destructVertices();
    if(vertexTable.isOpen) {
        vertexTable.remove();
    }
    if(reverseComplementVertex.isOpen) {
        reverseComplementVertex.remove();
    }
    if(edges.isOpen) {
        edges.remove();
    }
    if(reverseComplementEdge.isOpen) {
        reverseComplementEdge.remove();
    }
    if(edgeMarkerIntervals.isOpen()) {
        edgeMarkerIntervals.remove();
    }
    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }
    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }
    if(vertexRepeatCounts.isOpen) {
        vertexRepeatCounts.remove();
    }
    if(edgeConsensus.isOpen()) {
        edgeConsensus.remove();
    }
    if(edgeConsensusOverlappingBaseCount.isOpen) {
        edgeConsensusOverlappingBaseCount.remove();
    }
    if(vertexCoverageData.isOpen()) {
        vertexCoverageData.remove();
    }
    if(edgeCoverageData.isOpen()) {
        edgeCoverageData.remove();
    }
}



// Locate the edge given the vertices.
const MarkerGraph::Edge*
    MarkerGraph::findEdge(Uint40 source, Uint40 target) const
{
    const auto edgesWithThisSource = edgesBySource[source];
    for(const uint64_t i: edgesWithThisSource) {
        const Edge& edge = edges[i];
        if(edge.target == target) {
            return &edge;
        }
    }
    return 0;
}


MarkerGraph::EdgeId MarkerGraph::findEdgeId(Uint40 source, Uint40 target) const
{
    const Edge *edgePointer = findEdge(source, target);
    SHASTA_ASSERT(edgePointer);
    return edgePointer - edges.begin();
}

// Compute in-degree or out-degree of a vertex,
// counting only edges that were not removed.
uint64_t MarkerGraph::inDegree(VertexId vertexId) const
{
    uint64_t degree = 0;
    for(const EdgeId edgeId: edgesByTarget[vertexId]) {
        if(not edges[edgeId].wasRemoved()) {
            ++degree;
        }
    }
    return degree;
}
uint64_t MarkerGraph::outDegree(VertexId vertexId) const
{
    uint64_t degree = 0;
    for(const EdgeId edgeId: edgesBySource[vertexId]) {
        if(not edges[edgeId].wasRemoved()) {
            ++degree;
        }
    }
    return degree;
}



void MarkerGraph::Edge::writeFlags(ostream& s) const
{
    s << "wasRemovedByTransitiveReduction " << int(wasRemovedByTransitiveReduction) << "\n";
    s << "wasPruned " << int(wasPruned) << "\n";
    s << "isSuperBubbleEdge " << int(isSuperBubbleEdge) << "\n";
    s << "isLowCoverageCrossEdge " << int(isLowCoverageCrossEdge) << "\n";
    s << "wasAssembled " << int(wasAssembled) << "\n";
    s << "isSecondary " << int(isSecondary) << "\n";
    s << "wasRemovedWhileSplittingSecondaryEdges " << int(wasRemovedWhileSplittingSecondaryEdges) << "\n";
    s << flush;
}



MarkerGraph::EdgeId MarkerGraph::getFirstNonRemovedOutEdge(
    MarkerGraph::VertexId vertexId) const
{
    for(const EdgeId edgeId: edgesBySource[vertexId]) {
        if(not edges[edgeId].wasRemoved()) {
            return edgeId;
        }
    }
    return invalidEdgeId;
}
MarkerGraph::EdgeId MarkerGraph::getFirstNonRemovedInEdge(
    MarkerGraph::VertexId vertexId) const
{
    for(const EdgeId edgeId: edgesByTarget[vertexId]) {
        if(not edges[edgeId].wasRemoved()) {
            return edgeId;
        }
    }
    return invalidEdgeId;
}



// Remove marker graph vertices and update vertices and vertexTable.
void MarkerGraph::removeVertices(
    const MemoryMapped::Vector<VertexId>& verticesToBeKept,
    uint64_t pageSize,
    uint64_t threadCount)
{
    // Get the names of the data structures we will work with.
    const string verticesName = vertices().getName();
    const string vertexTableName = vertexTable.fileName;
    const string newVerticesName =
        verticesName.empty() ? "" : (verticesName + "-tmp");
    const string newVertexTableName =
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp");

    // Create the new vertices.
    removeVerticesData.verticesToBeKept = &verticesToBeKept;
    removeVerticesData.newVerticesPointer =
        make_shared<MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> >();
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& newVertices
        = *removeVerticesData.newVerticesPointer;
    newVertices.createNew(newVerticesName, pageSize);
    newVertices.beginPass1(verticesToBeKept.size());
    const uint64_t batchCount = 10000;
    setupLoadBalancing(verticesToBeKept.size(), batchCount);
    runThreads(&MarkerGraph::removeVerticesThreadFunction1, threadCount);
    newVertices.beginPass2();
    setupLoadBalancing(verticesToBeKept.size(), batchCount);
    runThreads(&MarkerGraph::removeVerticesThreadFunction2, threadCount);
    newVertices.endPass2(false, true);



    // Replace the old vertices with the new ones.
    vertices().remove();
    destructVertices();
    newVertices.rename(verticesName);
    verticesPointer = removeVerticesData.newVerticesPointer;
    removeVerticesData.newVerticesPointer = 0;


    // Update the vertexTable, in place.
    fill(vertexTable.begin(), vertexTable.end(), invalidCompressedVertexId);
    setupLoadBalancing(vertexCount(), batchCount);
    runThreads(&MarkerGraph::removeVerticesThreadFunction3, threadCount);



    // Remove everything else.
    if(reverseComplementVertex.isOpen) {
        reverseComplementVertex.remove();
    }
    if(edges.isOpen) {
        edges.remove();
    }
    if(edgeMarkerIntervals.isOpen()) {
        edgeMarkerIntervals.remove();
    }
    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }
    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }
    if(reverseComplementEdge.isOpen) {
        reverseComplementEdge.remove();
    }
    if(vertexRepeatCounts.isOpen) {
        vertexRepeatCounts.remove();
    }
    if(edgeConsensus.isOpen()) {
        edgeConsensus.remove();
    }
    if(edgeConsensusOverlappingBaseCount.isOpen) {
        edgeConsensusOverlappingBaseCount.remove();
    }
    if(vertexCoverageData.isOpen()) {
        vertexCoverageData.remove();
    }
    if(edgeCoverageData.isOpen()) {
        edgeCoverageData.remove();
    }
}



void MarkerGraph::removeVerticesThreadFunction1(size_t threadId)
{
    const MemoryMapped::Vector<VertexId>& verticesToBeKept =
        *removeVerticesData.verticesToBeKept;
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& newVertices
        = *removeVerticesData.newVerticesPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this thread.
        for(VertexId newVertexId=begin; newVertexId!=end; newVertexId++) {
            const VertexId oldVertexId = verticesToBeKept[newVertexId];
            newVertices.incrementCount(newVertexId, vertexCoverage(oldVertexId));
        }
    }
}



void MarkerGraph::removeVerticesThreadFunction2(size_t threadId)
{
    const MemoryMapped::Vector<VertexId>& verticesToBeKept =
        *removeVerticesData.verticesToBeKept;
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& newVertices
        = *removeVerticesData.newVerticesPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this thread.
        for(VertexId newVertexId=begin; newVertexId!=end; newVertexId++) {
            const VertexId oldVertexId = verticesToBeKept[newVertexId];
            copy(vertices().begin(oldVertexId), vertices().end(oldVertexId),
                newVertices.begin(newVertexId));
        }
    }
}



void MarkerGraph::removeVerticesThreadFunction3(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this thread.
        for(VertexId vertexId=begin; vertexId!=end; vertexId++) {
            const CompressedVertexId compressedVertexId = vertexId;
            const span<MarkerId> vertexMarkerIds = vertices()[vertexId];
            for(const MarkerId markerId: vertexMarkerIds) {
                vertexTable[markerId] = compressedVertexId;
            }
        }
    }
}



// This renumbers the vertex table to make sure that
// vertices are numbered contiguously starting at 0.
// This must be called after the vertexTable is changed,
// as in Assembler::cleanupDuplicateMarkers.
// After this is called, all other data structures
// are inconsistent and need to be recreated.
MarkerGraph::VertexId MarkerGraph::renumberVertexTable(size_t threadCount)
{
    // Sanity check.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);
    SHASTA_ASSERT(vertexTable.size() > 0);

    // Find the maximum vertex id.
    const VertexId maxVertexId = findMaxVertexTableEntry(threadCount);

    // Call the lower level version.
    return renumberVertexTable(threadCount, maxVertexId);

}



// This second version can be called if the maximum vertex id
// present in the vertex table is already known, and is faster.
MarkerGraph::VertexId MarkerGraph::renumberVertexTable(size_t threadCount, VertexId maxVertexId)
{
    const bool debug = false;

    // Sanity check.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);
    SHASTA_ASSERT(vertexTable.size() > 0);

    if(debug) {
        ofstream csv("VertexTable-Before.csv");
        csv << "MarkerId,VertexId\n";
        for(MarkerId markerId=0; markerId<vertexTable.size(); markerId++) {
            csv << markerId << "," << vertexTable[markerId] << "\n";
        }
    }

    cout << timestamp << "Renumbering the marker graph vertex table." << endl;

    // Create a vector of bools that tells us which VertexId's are present.
    const string vertexTableName = vertexTable.fileName;
    renumberVertexTableData.isPresent.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-isPresent"),
        vertexTable.getPageSize());
    renumberVertexTableData.isPresent.resize(maxVertexId + 1);
    fill(
        renumberVertexTableData.isPresent.begin(),
        renumberVertexTableData.isPresent.end(),
        false);
    const uint64_t batchSize = 100000;
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::renumberVertexTableThreadFunction1, threadCount);

    // Now we know what VertexId's are present, so we can compute the new VertexId
    // corresponding to each old VertexId.
    renumberVertexTableData.newVertexId.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-newVertexId"),
        vertexTable.getPageSize());
    renumberVertexTableData.newVertexId.resize(maxVertexId + 1);
    VertexId newVertexId = 0;
    for(VertexId oldVertexId=0; oldVertexId<=maxVertexId; oldVertexId++) {
        if(renumberVertexTableData.isPresent[oldVertexId]) {
            renumberVertexTableData.newVertexId[oldVertexId] = newVertexId;
            ++newVertexId;
        } else {
            renumberVertexTableData.newVertexId[oldVertexId] = invalidVertexId;
        }
    }

    // Now we can renumber the vertex table.
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::renumberVertexTableThreadFunction2, threadCount);

    // Clean up.
    renumberVertexTableData.newVertexId.remove();
    renumberVertexTableData.isPresent.remove();

    if(debug) {
        ofstream csv("VertexTable-After.csv");
        csv << "MarkerId,VertexId\n";
        for(MarkerId markerId=0; markerId<vertexTable.size(); markerId++) {
            csv << markerId << "," << vertexTable[markerId] << "\n";
        }
    }

    cout << timestamp << "Done renumbering the marker graph vertex table." << endl;
    return newVertexId - 1;
}



void MarkerGraph::renumberVertexTableThreadFunction1(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];
            if(compressedVertexId != invalidCompressedVertexId) {
                renumberVertexTableData.isPresent[VertexId(compressedVertexId)] = true;
            }
        }
    }
}



void MarkerGraph::renumberVertexTableThreadFunction2(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];
            if(compressedVertexId != invalidCompressedVertexId) {
                const VertexId oldVertexId = VertexId(compressedVertexId);
                const VertexId newVertexId = renumberVertexTableData.newVertexId[oldVertexId];
                vertexTable[markerId] = CompressedVertexId(newVertexId);
            }
        }
    }
}



MarkerGraph::VertexId MarkerGraph::findMaxVertexTableEntry(size_t threadCount)
{
    // Sanity checks.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);

    // Initialize the maximum VertexId found by each thread.
    findMaxVertexTableEntryData.threadMaxVertexId.resize(threadCount);
    fill(
        findMaxVertexTableEntryData.threadMaxVertexId.begin(),
        findMaxVertexTableEntryData.threadMaxVertexId.end(),
        0);

    // Each thread finds the maximum of a subset of the vertex table.
    const uint64_t batchSize = 100000;
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::findMaxVertexTableEntryThreadFunction, threadCount);

    // Return the maximum value found by all threads.
    return *std::max_element(
        findMaxVertexTableEntryData.threadMaxVertexId.begin(),
        findMaxVertexTableEntryData.threadMaxVertexId.end());
}



void MarkerGraph::findMaxVertexTableEntryThreadFunction(size_t threadId)
{
    VertexId maxVertexId = 0;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];
            if(compressedVertexId != invalidCompressedVertexId) {
                maxVertexId = max(maxVertexId, VertexId(compressedVertexId));
            }
        }
    }

    findMaxVertexTableEntryData.threadMaxVertexId[threadId] = maxVertexId;
}




// Recreate the vertices from the vertexTable.
// This assumes that valid VertexId's in the vertex table
// are numbered contiguously starting at 0 (call renumberVertexTable to ensure that).
void MarkerGraph::createVerticesFromVertexTable(size_t threadCount, VertexId maxVertexId)
{
    SHASTA_ASSERT(vertexTable.isOpen);
    const string vertexTableName = vertexTable.fileName;

    // Create our copy of the vertices.
    createVerticesFromVertexTableData.vertices.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-isPresent"),
        vertexTable.getPageSize());

    // Pass 1: count the number of markers in each vertex.
    createVerticesFromVertexTableData.vertices.beginPass1(maxVertexId+1);
    uint64_t batchSize = 100000;
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction1, threadCount);

    // Pass 2: fill in the markers of each vertex.
    createVerticesFromVertexTableData.vertices.beginPass2();
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction2, threadCount);
    createVerticesFromVertexTableData.vertices.endPass2(false);

    // Pass 3: sort the markers of each vertex.
    batchSize = 1000;
    setupLoadBalancing(maxVertexId+1, batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction3, threadCount);

    // Finally, we copy to the main copy of the vertices.
    vertices().clear();
    const uint64_t vertexCount = createVerticesFromVertexTableData.vertices.size();
    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        vertices().appendVector(createVerticesFromVertexTableData.vertices.size(vertexId));
    }
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction4, threadCount);

    // Cleanup.
    createVerticesFromVertexTableData.vertices.remove();
}



void MarkerGraph::createVerticesFromVertexTableThreadFunction1(size_t threadId)
{
    auto& vertices = createVerticesFromVertexTableData.vertices;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];

            if(compressedVertexId != invalidCompressedVertexId) {
                vertices.incrementCountMultithreaded(VertexId(compressedVertexId));
            }
        }
    }

}



void MarkerGraph::createVerticesFromVertexTableThreadFunction2(size_t threadId)
{
    auto& vertices = createVerticesFromVertexTableData.vertices;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];

            if(compressedVertexId != invalidCompressedVertexId) {
                vertices.storeMultithreaded(VertexId(compressedVertexId), markerId);
            }
        }
    }

}



void MarkerGraph::createVerticesFromVertexTableThreadFunction3(size_t threadId)
{
    auto& vertices = createVerticesFromVertexTableData.vertices;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over verteices in this batch.
        for(VertexId vertexId=begin; vertexId!=end; vertexId++) {
            auto vertexMarkers = vertices[vertexId];
            sort(vertexMarkers.begin(), vertexMarkers.end());
        }
    }

}



void MarkerGraph::createVerticesFromVertexTableThreadFunction4(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices in this batch.
        for(VertexId vertexId=begin; vertexId!=end; vertexId++) {

            // Make the copy.
            auto vertexMarkers = createVerticesFromVertexTableData.vertices[vertexId];
            copy(vertexMarkers.begin(), vertexMarkers.end(), verticesPointer->begin(vertexId));
        }
    }

}



// Find the common KmerId for all the markers of a marker graph vertex.
KmerId MarkerGraph::getVertexKmerId(
    MarkerGraphVertexId vertexId,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
    ) const
{
    // Get it from the first marker on this vertex.
    const MarkerId markerId = getVertexMarkerIds(vertexId)[0];

    // Find the OrientedReadId.
    // This is slow as it requires a binary search in the markers toc.
    OrientedReadId orientedReadId;
    uint32_t ordinal;
    tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);

    return getOrientedReadMarkerKmerId(
        orientedReadId,
        ordinal,
        k,
        reads,
        markers
        );
}



// Find the edge that contains a given MarkerInterval.
MarkerGraphEdgeId MarkerGraph::locateMarkerInterval(
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerInterval& markerInterval) const
{
    const OrientedReadId orientedReadId = markerInterval.orientedReadId;
    const uint64_t firstOrientedReadMarkerId =
        markers.begin(orientedReadId.getValue()) - markers.begin();

    // Now locate this marker interval.
    const uint64_t markerId0 = firstOrientedReadMarkerId + markerInterval.ordinals[0];
    const uint64_t markerId1 = firstOrientedReadMarkerId + markerInterval.ordinals[1];
    const MarkerGraphVertexId vertexId0 = vertexTable[markerId0];
    const MarkerGraphVertexId vertexId1 = vertexTable[markerId1];

    for(const auto edgeId: edgesBySource[vertexId0]) {
        if(edges[edgeId].target != vertexId1) {
            continue;
        }
        const auto markerIntervals = edgeMarkerIntervals[edgeId];
        if(find(markerIntervals.begin(), markerIntervals.end(), markerInterval) !=
            markerIntervals.end()) {
            return edgeId;
        }
    }

    return invalid<MarkerGraphEdgeId>;
}


// Apply an ordinal offset in the specified direction to a given MarkerInterval
// and find the edge that contains the offset MarkerInterval.
// This assumes that we have the complete marker graph.
MarkerGraphEdgeId MarkerGraph::locateMarkerIntervalWithOffset(
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MarkerInterval markerInterval,
    uint32_t ordinalOffset,
    uint64_t direction // 0=forward, 1=backward.
    ) const
{
    const OrientedReadId orientedReadId = markerInterval.orientedReadId;
    const uint64_t firstOrientedReadMarkerId =
        markers.begin(orientedReadId.getValue()) - markers.begin();

    // Construct the offset MarkerInterval.
    // If we end up outside the oriented read, return invalid<MarkerGraphEdgeId>.
    if(direction == 0) {
        markerInterval.ordinals[0] += ordinalOffset;
        markerInterval.ordinals[1] += ordinalOffset;
        if(markerInterval.ordinals[1] >= markers.size(orientedReadId.getValue())) {
            return invalid<MarkerGraphEdgeId>;
        }
    } else {
        if(ordinalOffset > markerInterval.ordinals[0]) {
            return invalid<MarkerGraphEdgeId>;
        }
        markerInterval.ordinals[0] -= ordinalOffset;
        markerInterval.ordinals[1] -= ordinalOffset;
    }
    SHASTA_ASSERT(markerInterval.ordinals[1] == markerInterval.ordinals[0] + 1);


    // Now locate this marker interval.
    const uint64_t markerId0 = firstOrientedReadMarkerId + markerInterval.ordinals[0];
    const uint64_t markerId1 = firstOrientedReadMarkerId + markerInterval.ordinals[1];
    const MarkerGraphVertexId vertexId0 = vertexTable[markerId0];
    const MarkerGraphVertexId vertexId1 = vertexTable[markerId1];

    for(const auto edgeId: edgesBySource[vertexId0]) {
        if(edges[edgeId].target != vertexId1) {
            continue;
        }
        const auto markerIntervals = edgeMarkerIntervals[edgeId];
        if(find(markerIntervals.begin(), markerIntervals.end(), markerInterval) !=
            markerIntervals.end()) {
            return edgeId;
        }
    }

    // If this happens, we don't have a complete marker graph.
    SHASTA_ASSERT(0);
}



// Find out if an edge has duplicate oriented reads
// in its MarkerIntervals.
bool MarkerGraph::edgeHasDuplicateOrientedReadIds(EdgeId edgeId) const
{
    const auto markerIntervals = edgeMarkerIntervals[edgeId];
    if(markerIntervals.size() < 2) {
        return false;
    }
    for(uint64_t i=1; i<markerIntervals.size(); i++) {
        if(markerIntervals[i-1].orientedReadId == markerIntervals[i].orientedReadId) {
            return true;
        }
    }

    return false;
}



// Find out if a vertex has more than one marker on the same oriented read.
bool MarkerGraph::vertexHasDuplicateOrientedReadIds(
    VertexId vertexId,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) const
{
    const span<const MarkerId> vertexMarkerIds = vertices()[vertexId];
    if(vertexMarkerIds.size() < 2) {
        return false;
    }

    // The markers are sorted, so we only have to check each marker
    // against the previous one.
    // This could be done faster but is not performance critical.
    for(uint64_t i=1; i<vertexMarkerIds.size(); i++) {
        const MarkerId markerId0 = vertexMarkerIds[i-1];
        const MarkerId markerId1 = vertexMarkerIds[i];
        OrientedReadId orientedReadId0;
        OrientedReadId orientedReadId1;
        tie(orientedReadId0, ignore) = findMarkerId(markerId0, markers);
        tie(orientedReadId1, ignore) = findMarkerId(markerId1, markers);
        if(orientedReadId0 == orientedReadId1) {
            return true;
        }
    }

    return false;
}



// Flag primary edges (only used for Mode 3 assembly).
void MarkerGraph::flagPrimaryEdges(
    uint64_t minPrimaryEdgeCoverage,
    uint64_t maxPrimaryEdgeCoverage,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MemoryMapped::Vector< pair<uint64_t, uint64_t> >& disjointSetsHistogram,
    uint64_t threadCount)
{
    // If minPrimaryEdgeCoverage and maxPrimaryEdgeCoverage are both 0,
    // use the disjoint sets histogram and simple heuristics to choose
    // appropriate values.
    if((minPrimaryEdgeCoverage == 0) and (maxPrimaryEdgeCoverage == 0)) {

        // Set minPrimaryEdgeCoverage to the first value where the
        // disjointSetsHistogram starts increasing.
        bool done = false;
        uint64_t frequencyAtMinPrimaryEdgeCoverage = 0;
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
                minPrimaryEdgeCoverage = coverage;
                frequencyAtMinPrimaryEdgeCoverage = frequency;
                done = true;
                break;
            }
        }
        SHASTA_ASSERT(done);

        // Set maxPrimaryEdgeCoverage to the last coverage with frequency
        // at least equal to frequencyAtMinPrimaryEdgeCoverage.
        done = false;
        for(uint64_t i=disjointSetsHistogram.size()-1; i>0; i--) {
            const uint64_t coverage = disjointSetsHistogram[i].first;
            const uint64_t frequency = disjointSetsHistogram[i].second;
            if(frequency >= frequencyAtMinPrimaryEdgeCoverage) {
                maxPrimaryEdgeCoverage = coverage;
                done= true;
                break;
            }
        }
        SHASTA_ASSERT(done);

        cout << "Automatically set: minPrimaryEdgeCoverage = " << minPrimaryEdgeCoverage <<
            ", maxPrimaryEdgeCoverage = " << maxPrimaryEdgeCoverage << endl;
    }



    // Store the arguments so the threads can see them.
    flagPrimaryEdgesData.minPrimaryEdgeCoverage = minPrimaryEdgeCoverage;
    flagPrimaryEdgesData.maxPrimaryEdgeCoverage = maxPrimaryEdgeCoverage;
    flagPrimaryEdgesData.markersPointer = &markers;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Clear the flags on all edges.
    for(Edge& edge: edges) {
        edge.isPrimary = 0;
    }

    // Multithreaded code to flag primary edges.
    const uint64_t batchCount = 10000;
    setupLoadBalancing(edges.size(), batchCount);
    runThreads(&MarkerGraph::flagPrimaryEdgesThreadFunction, threadCount);

    uint64_t primaryEdgeCount = 0;
    for(Edge& edge: edges) {
        if(edge.isPrimary == 1) {
            ++primaryEdgeCount;
        }
    }
    cout << "Found " << primaryEdgeCount <<
        " primary marker graph edges out of " << edges.size() << " total." << endl;
}



void MarkerGraph::flagPrimaryEdgesThreadFunction(uint64_t threadId)
{
    const uint64_t minPrimaryEdgeCoverage = flagPrimaryEdgesData.minPrimaryEdgeCoverage;
    const uint64_t maxPrimaryEdgeCoverage = flagPrimaryEdgesData.maxPrimaryEdgeCoverage;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers =
        *flagPrimaryEdgesData.markersPointer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(EdgeId edgeId=begin; edgeId!=end; ++edgeId) {
            if(isPrimaryEdge(edgeId, minPrimaryEdgeCoverage, maxPrimaryEdgeCoverage, markers)) {
                edges[edgeId].isPrimary = 1;
            }
        }
    }
}



// Find out if a marker graph edge is a primary edge.
// Only used for Mode 3 assembly.
bool MarkerGraph::isPrimaryEdge(
    MarkerGraphEdgeId edgeId,
    uint64_t minPrimaryEdgeCoverage,
    uint64_t maxPrimaryEdgeCoverage,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) const
{
    // Check coverage.
    const uint64_t coverage = edgeCoverage(edgeId);
    if(coverage < minPrimaryEdgeCoverage) {
        return false;
    }
    if(coverage > maxPrimaryEdgeCoverage) {
        return false;
    }

    // Check for duplicate oriented reads on the edge.
    if(edgeHasDuplicateOrientedReadIds(edgeId)) {
        return false;
    }

    // Check for duplicate oriented reads on its vertices.
    const MarkerGraph::Edge& edge = edges[edgeId];
    if(
        vertexHasDuplicateOrientedReadIds(edge.source, markers) or
        vertexHasDuplicateOrientedReadIds(edge.target, markers)) {
        return false;
    }

    // If all above checks passed, this is a primary edge.
    return true;
}



#if 0
void MarkerGraph::createPrimaryJourneys(
    uint64_t orientedReadCount,
    uint64_t threadCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    primaryJourneys.clear();

    const uint64_t batchCount = 100;

    primaryJourneys.beginPass1(orientedReadCount);
    setupLoadBalancing(edges.size(), batchCount);
    runThreads(&MarkerGraph::createPrimaryJourneysThreadFunction1, threadCount);
    primaryJourneys.beginPass2();
    setupLoadBalancing(edges.size(), batchCount);
    runThreads(&MarkerGraph::createPrimaryJourneysThreadFunction2, threadCount);
    primaryJourneys.endPass2(false, true);
    setupLoadBalancing(orientedReadCount, 1);
    runThreads(&MarkerGraph::createPrimaryJourneysThreadFunction3, threadCount);

    cout << "Found " << primaryJourneys.totalSize() <<
        " marker graph primary journey entries for " << orientedReadCount <<
        " oriented reads." << endl;
    cout << "Average number of marker graph primary journey entries per oriented read is " <<
        double(primaryJourneys.totalSize()) / double(orientedReadCount) << endl;

    writePrimaryJourneys();
}



void MarkerGraph::writePrimaryJourneys()
{
    const uint64_t orientedReadCount = primaryJourneys.size();

    ofstream csv("MarkerGraphPrimaryJourneys.csv");

    for(ReadId readId=0; readId<orientedReadCount/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            csv << orientedReadId << ",";
            for(const auto& primaryJourneyEntry: primaryJourneys[orientedReadId.getValue()]) {
                csv << primaryJourneyEntry.edgeId << ",";
            }
            csv << "\n";
        }
    }
}



void MarkerGraph::createPrimaryJourneysThreadFunction1(uint64_t threadId)
{
    createPrimaryJourneysThreadFunction12(1);
}



void MarkerGraph::createPrimaryJourneysThreadFunction2(uint64_t threadId)
{
    createPrimaryJourneysThreadFunction12(2);
}



void MarkerGraph::createPrimaryJourneysThreadFunction12(uint64_t  pass)
{
    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(EdgeId edgeId=begin; edgeId!=end; ++edgeId) {
            const Edge& edge = edges[edgeId];

            // If this is not a primary edge, skip it.
            if(edge.isPrimary == 0) {
                continue;
            }

            PrimaryJourneyEntry primaryJourneyEntry;
            primaryJourneyEntry.edgeId = edgeId;

            // Loop over the MarkerIntervals of this edge.
            span<MarkerInterval> markerIntervals = edgeMarkerIntervals[edgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                const uint64_t orientedReadIdValue = markerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    primaryJourneys.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    primaryJourneyEntry.ordinals = markerInterval.ordinals;
                    primaryJourneys.storeMultithreaded(orientedReadIdValue, primaryJourneyEntry);
                }

            }

        }
    }
}



void MarkerGraph::createPrimaryJourneysThreadFunction3(uint64_t threadId)
{
    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(uint64_t orientedReadIdValue=begin; orientedReadIdValue!=end; orientedReadIdValue++) {
            auto journey = primaryJourneys[orientedReadIdValue];
            sort(journey.begin(), journey.end());
        }
    }
}



// Starting from a primary marker graph edge, follow the primary journeys
// of all oriented reads on the edge, moving forward.
// Find the set of MarkerGraphEdgeIds that were encountered in this way,
// and for each the number of times it was encountered.
void MarkerGraph::followPrimaryJourneysForward(
    MarkerGraphEdgeId edgeId0,
    vector<MarkerGraphEdgeId>& edgeIds,
    vector<uint64_t>& count) const
{
    edgeIds.clear();
    count.clear();

    // Loop over the oriented reads in edgeId0.
    for(const MarkerInterval& markerInterval: edgeMarkerIntervals[edgeId0]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto primaryJourney = primaryJourneys[orientedReadId.getValue()];

        // Loop over the primary journey backward, stopping when we encounter edgeId0.
        for(uint64_t j=primaryJourney.size(); /* Check later */; --j) {
            const auto& primaryJourneyEntry = primaryJourney[j];
            const MarkerGraphEdgeId edgeId1 = primaryJourneyEntry.edgeId;
            if(edgeId1 == edgeId0) {
                break;
            }
            edgeIds.push_back(edgeId1);
            if(j == 0) {
                break;
            }
        }
    }

    deduplicateAndCount(edgeIds, count);
    SHASTA_ASSERT(edgeIds.size() == count.size());

}



// Same, but moving backward.
void MarkerGraph::followPrimaryJourneysBackward(
    MarkerGraphEdgeId edgeId0,
    vector<MarkerGraphEdgeId>& edgeIds,
    vector<uint64_t>& count) const
{
    edgeIds.clear();
    count.clear();

    // Loop over the oriented reads in edgeId0.
    for(const MarkerInterval& markerInterval: edgeMarkerIntervals[edgeId0]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto primaryJourney = primaryJourneys[orientedReadId.getValue()];

        // Loop over the primary journey, stopping when we encounter edgeId0.
        for(const auto& primaryJourneyEntry: primaryJourney) {
            const MarkerGraphEdgeId edgeId1 = primaryJourneyEntry.edgeId;
            if(edgeId1 == edgeId0) {
                break;
            }
            edgeIds.push_back(edgeId1);
        }
    }

    deduplicateAndCount(edgeIds, count);
    SHASTA_ASSERT(edgeIds.size() == count.size());

}
#endif

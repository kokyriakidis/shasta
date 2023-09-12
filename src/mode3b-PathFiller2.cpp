// Shasta.
#include "mode3b-PathFiller2.hpp"
#include "approximateTopologicalSort.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "globalMsa.hpp"
#include "orderPairs.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>


PathFiller2::PathFiller2(
    const Assembler& assembler,
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    const PathFiller2DisplayOptions& options) :
    assembler(assembler),
    edgeIdA(edgeIdA),
    edgeIdB(edgeIdB),
    options(options),
    html(options.html)

{
    // PARAMETERS THAT SHOULD BE EXPOSED WHEN CODE STABILIZES.

    // The estimated offset gets extended by this ratio to
    // decide how much to extend reads that only appear in edgeIdA or edgeIdB.
    double estimatedOffsetRatio = 1.1;

    // The minimum coverage for a vertex to be created.
    // This gets reduced if necessary to enforce the maximum MSA length.
    uint64_t minVertexCoverage = 20;

    // Control vertex splitting.
    const int64_t maxBaseSkip = 300;

    // The maximum allowed length for a single MSA.
    const uint64_t maxMsaLength = 10000;

    // Reduce minVertexCoverage until we succeed.
    while(true) {
        if(html) {
            html << "<br>Using minVertexCoverage = " << minVertexCoverage;
        }
        const bool success = run(
            estimatedOffsetRatio,
            minVertexCoverage,
            maxBaseSkip,
            maxMsaLength);
        if(success) {
            break;
        } else {
            clearAll();
            --minVertexCoverage;
            if(minVertexCoverage == 0) {
                cout << "Assertion occurred while assembling between marker graph edges " <<
                    edgeIdA << " " << edgeIdB << endl;
                SHASTA_ASSERT(0);
            }
        }
    }
}



void PathFiller2::clearAll()
{
    PathFiller2BaseClass::clear();
    SHASTA_ASSERT(num_vertices(*this) == 0);
    SHASTA_ASSERT(num_edges(*this) == 0);
    orientedReadInfos.clear();
    assemblyPath.clear();
}



bool PathFiller2::run(
    double estimatedOffsetRatio,
    uint64_t minVertexCoverage,
    int64_t maxBaseSkip,
    uint64_t maxMsaLength)
{

    // Store the source target of edgeIdA and the source vertex of edgeIdB.
    const MarkerGraph::Edge& edgeA = assembler.markerGraph.edges[edgeIdA];
    const MarkerGraph::Edge& edgeB = assembler.markerGraph.edges[edgeIdB];
    vertexIdA = edgeA.target;
    vertexIdB = edgeB.source;


    // Oriented reads.
    checkAssumptions();
    gatherOrientedReads();
    writeOrientedReads();
    estimateOffset();

    // Vertices.
    createVertices(estimatedOffsetRatio);
    removeLowCoverageVertices(minVertexCoverage);
    splitVertices(maxBaseSkip);

    // Edges.
    createEdges();
    writeGraph("Initial assembly graph");

    // Remove strongly connected components, then regenerate
    // edges from scratch with the remaining vertices.
    removeStrongComponents();
    removeAllEdges();
    createEdges();
    writeGraph("Assembly graph after removal of strong connected components");

    // Remove vertices that are not accessible from vertexIdA
    // or from which vertexidB is not accessible.
    const uint64_t inaccessibleVertexCount = removeInaccessibleVertices();
    if(inaccessibleVertexCount>0 and html and options.showDebugInformation) {
        writeGraph("Assembly graph after removal of inaccessible vertices");
    }

    if(html and options.showDebugInformation) {
        vertexCoverageHistogram();
        edgeCoverageHistogram();
    }

    // Assemble.
    findAssemblyPath();
    if(not assembleAssemblyPathEdges(maxMsaLength)) {
        return false;
    }
    writeGraph("Assembly graph after assembly");

    // Write assembled sequence.
    if(html) {
        vector<Base> sequence;
        getSecondarySequence(sequence);

        html <<
            "<h2>Assembled sequence</h2>"
            "Assembled sequence not including the first and last edge is " <<
            sequence.size() << " bases long."
            "<pre style='font-family:monospace'>\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));
        html << "</pre>";

        ofstream fasta("PathFiller2.fasta");
        fasta << ">PathFiller2 " << sequence.size() << endl;
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));

        writeVerticesCsv();
    }

    return true;
}



void PathFiller2::checkAssumptions() const
{
    SHASTA_ASSERT(edgeIdA != edgeIdB);
    SHASTA_ASSERT(assembler.assemblerInfo->assemblyMode == 3);
    SHASTA_ASSERT(assembler.getReads().representation == 0);
    SHASTA_ASSERT(not assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA));
    SHASTA_ASSERT(not assembler.markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB));

    const MarkerGraph& markerGraph = assembler.markerGraph;
    const auto& markers = assembler.markers;

    // edgeIdA and edgeIdB cannot have duplicate oriented reads.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA)) {
        throw runtime_error("Duplicated oriented read on edgeIdA.");
    }
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB)) {
        throw runtime_error("Duplicated oriented read on edgeIdB.");
    }

    // Neither can their source and target vertices.
    if(markerGraph.vertexHasDuplicateOrientedReadIds(vertexIdA, markers)) {
        throw runtime_error("Duplicated oriented read on target vertex of edgeIdA.");
    }
    if(markerGraph.vertexHasDuplicateOrientedReadIds(vertexIdB, markers)) {
        throw runtime_error("Duplicated oriented read on source vertex of edgeIdB.");
    }
}



#if 0
// Gather oriented reads and fill in the orientedReadInfos.
// We have already checked that no reads appears twice on edgeIdA or edgeIdB.
void PathFiller2::gatherOrientedReads()
{
    // Gather the OrientedReadIds that appear on edgeIdA or edgeIdB.
    vector<OrientedReadId> orientedReadIds;
    for(const MarkerInterval& markerInterval: assembler.markerGraph.edgeMarkerIntervals[edgeIdA]) {
        orientedReadIds.push_back(markerInterval.orientedReadId);
    }
    for(const MarkerInterval& markerInterval: assembler.markerGraph.edgeMarkerIntervals[edgeIdB]) {
        orientedReadIds.push_back(markerInterval.orientedReadId);
    }
    deduplicate(orientedReadIds);

    // Store the OrientedReadIds.
    orientedReadInfos.clear();
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        orientedReadInfos.push_back(OrientedReadInfo(orientedReadId));
    }

    // Fill in the OrdinalAndPositions of the oriented reads that appear in edgeIdA.
    for(const MarkerInterval& markerInterval: assembler.markerGraph.edgeMarkerIntervals[edgeIdA]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const uint64_t i = getOrientedReadIndex(orientedReadId);
        OrientedReadInfo& info = orientedReadInfos[i];

        const uint32_t ordinal = markerInterval.ordinals[1];    // Because vertexIdA is the target of edgeIdA
        const MarkerId markerId = assembler.getMarkerId(orientedReadId, ordinal);
        info.ordinalAndPositionA = {ordinal, assembler.markers.begin()[markerId].position};
    }

    // Fill in the OrdinalAndPositions of the oriented reads that appear in edgeIdB.
    for(const MarkerInterval& markerInterval: assembler.markerGraph.edgeMarkerIntervals[edgeIdB]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const uint64_t i = getOrientedReadIndex(orientedReadId);
        OrientedReadInfo& info = orientedReadInfos[i];

        const uint32_t ordinal = markerInterval.ordinals[0];    // Because vertexIdB is the source of edgeIdB
        const MarkerId markerId = assembler.getMarkerId(orientedReadId, ordinal);
        info.ordinalAndPositionB = {ordinal, assembler.markers.begin()[markerId].position};
    }
}
#endif



// Gather oriented reads and fill in the orientedReadInfos.
// We have already checked that no reads appears twice on edgeIdA or edgeIdB.
void PathFiller2::gatherOrientedReads()
{
    // Joint loop over marker intervals that appear in edgeIdA and/or edgeIdB.
    const auto markerIntervalsA = assembler.markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = assembler.markerGraph.edgeMarkerIntervals[edgeIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();
    auto itA = beginA;
    auto itB = beginB;
    while(true) {
        if((itA == endA) and (itB == endB)) {
            break;
        }

        // Oriented reads that appear only in edgeIdA.
        if((itB == endB) or (itA != endA and itA->orientedReadId < itB->orientedReadId)) {
            const MarkerInterval& markerIntervalA = *itA;
            const OrientedReadId orientedReadIdA = markerIntervalA.orientedReadId;
            const uint32_t ordinalA = markerIntervalA.ordinals[1];    // Because vertexIdA is the target of edgeIdA
            const MarkerId markerIdA = assembler.getMarkerId(orientedReadIdA, ordinalA);

            OrientedReadInfo info(orientedReadIdA);
            info.ordinalAndPositionA = {ordinalA, assembler.markers.begin()[markerIdA].position};
            orientedReadInfos.push_back(info);

            ++itA;
        }



        // Oriented reads that appear only in edgeIdB.
        else if((itA == endA) or (itB != endB and itB->orientedReadId < itA->orientedReadId)) {
            const MarkerInterval& markerIntervalB = *itB;
            const OrientedReadId orientedReadIdB = markerIntervalB.orientedReadId;
            const uint32_t ordinalB = markerIntervalB.ordinals[0];    // Because vertexIdB is the source of edgeIdB
            const MarkerId markerIdB = assembler.getMarkerId(orientedReadIdB, ordinalB);

            OrientedReadInfo info(orientedReadIdB);
            info.ordinalAndPositionB = {ordinalB, assembler.markers.begin()[markerIdB].position};
            orientedReadInfos.push_back(info);

            ++itB;
        }

        // Oriented reads that appear in both edgeIdA and edgeIdB.
        else {
            SHASTA_ASSERT(itA != endA);
            SHASTA_ASSERT(itB != endB);

            const MarkerInterval& markerIntervalA = *itA;
            const OrientedReadId orientedReadIdA = markerIntervalA.orientedReadId;

            const MarkerInterval& markerIntervalB = *itB;
            const OrientedReadId orientedReadIdB = markerIntervalB.orientedReadId;

            SHASTA_ASSERT(orientedReadIdA == orientedReadIdB);
            const OrientedReadId orientedReadId = orientedReadIdA;

            const uint32_t ordinalA = markerIntervalA.ordinals[1];    // Because vertexIdA is the target of edgeIdA
            const uint32_t ordinalB = markerIntervalB.ordinals[0];    // Because vertexIdB is the source of edgeIdB

            // Only use it if the ordinal offset is not negative.
            if(ordinalB >= ordinalA) {
                const MarkerId markerIdA = assembler.getMarkerId(orientedReadIdA, ordinalA);
                const MarkerId markerIdB = assembler.getMarkerId(orientedReadIdB, ordinalB);

                OrientedReadInfo info(orientedReadId);
                info.ordinalAndPositionA = {ordinalA, assembler.markers.begin()[markerIdA].position};
                info.ordinalAndPositionB = {ordinalB, assembler.markers.begin()[markerIdB].position};
                orientedReadInfos.push_back(info);
            }

            ++itA;
            ++itB;
        }

    }
}



void PathFiller2::writeOrientedReads() const
{
    if(not html) {
        return;
    }

    html <<
        "<h2>Oriented reads</h2>"
        "<table>"
        "<tr>"
        "<th>Index"
        "<th>Oriented<br>read"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Ordinal<br>offset"
        "<th>PositionA"
        "<th>PositionB"
         "<th>Position<br>offset"
        ;

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered>" << info.orientedReadId;

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalAndPositionA.ordinal;
        }

        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalAndPositionB.ordinal;
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            html << info.ordinalOffset();
        }

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalAndPositionA.position;
        }

        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalAndPositionB.position;
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            html << info.positionOffset();
        }
     }

    html << "</table>";
}



// The index of an OrientedReadId is its index in the orientedReadInfos vector.
uint64_t PathFiller2::getOrientedReadIndex(OrientedReadId orientedReadId) const
{
    // Look it up.
    auto it = std::lower_bound(
        orientedReadInfos.begin(),
        orientedReadInfos.end(),
        OrientedReadInfo(orientedReadId));

    // Check that we found it.
    SHASTA_ASSERT(it != orientedReadInfos.end());
    SHASTA_ASSERT(it->orientedReadId == orientedReadId);

    // Return its index.
    return it - orientedReadInfos.begin();
}



void PathFiller2::estimateOffset()
{
    int64_t n = 0;
    int64_t sum = 0;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        if(info.isOnA() and info.isOnB()) {
            n++;
            sum += info.positionOffset();
        }
    }
    estimatedABOffset = int64_t(std::round(double(sum) / double(n)));

    if(html and options.showDebugInformation) {
        html << "<br>Estimated offset is " << estimatedABOffset << " bases.";
    }
}



// Create vertices by following the portion of each oriented read we will use for assembly.
void PathFiller2::createVertices(double estimatedOffsetRatio)
{
    PathFiller2& graph = *this;
    const int64_t offsetThreshold = int64_t(estimatedOffsetRatio * double(estimatedABOffset));

    // During this phase there is at most one vertex corresponding to
    // each marker graph vertex.
    std::map<MarkerGraphVertexId, vertex_descriptor> vertexMap;


    // Loop over our oriented reads.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        OrientedReadInfo& info = orientedReadInfos[i];
        const OrientedReadId orientedReadId = info.orientedReadId;
        info.vertices.clear();

        // Oriented reads that appear on both edgeIdA and edgeIdB.
        if(info.isOnA() and info.isOnB()) {
            info.firstOrdinal = info.ordinalAndPositionA.ordinal;
            info.lastOrdinal =  info.ordinalAndPositionB.ordinal;
            for(int64_t ordinal=info.firstOrdinal; ordinal<=info.lastOrdinal; ordinal++) {
                createVerticesHelper(i, ordinal, vertexMap);
            }
        }

        // Oriented reads that appear on edgeIdA but not on edgeIdB.
        else if(info.isOnA() and not info.isOnB()) {
            info.firstOrdinal = info.ordinalAndPositionA.ordinal;
            const int64_t maxPosition = info.ordinalAndPositionA.position + offsetThreshold;
            const int64_t markerCount = int64_t(assembler.markers.size(orientedReadId.getValue()));
            for(int64_t ordinal=info.firstOrdinal; ordinal<markerCount; ordinal++) {
                const MarkerId markerId = assembler.getMarkerId(orientedReadId, uint32_t(ordinal));
                const int64_t position = int64_t(assembler.markers.begin()[markerId].position);
                if(position > maxPosition) {
                    break;
                }
                createVerticesHelper(i, ordinal, vertexMap);
                info.lastOrdinal = ordinal;
            }
        }

        // Oriented reads that appear on edgeIdB but not on edgeIdA.
        else if(info.isOnB() and not info.isOnA()) {
            info.lastOrdinal = info.ordinalAndPositionB.ordinal;
            const int64_t minPosition = info.ordinalAndPositionB.position - offsetThreshold;
            for(int64_t ordinal=info.lastOrdinal; ordinal>=0; ordinal--) {
                const MarkerId markerId = assembler.getMarkerId(orientedReadId, uint32_t(ordinal));
                const int64_t position = int64_t(assembler.markers.begin()[markerId].position);
                if(position < minPosition) {
                    break;
                }
                createVerticesHelper(i, ordinal, vertexMap);
                info.firstOrdinal = ordinal;
            }
            reverse(info.vertices.begin(), info.vertices.end());
        }

        else {
            SHASTA_ASSERT(0);
        }
    }

    if(html and options.showDebugInformation) {
        html << "<br>After initial creation, there are " << num_vertices(graph) << " vertices." << endl;
    }
}



void PathFiller2::createVerticesHelper(
    uint64_t i,
    int64_t ordinal,
    std::map<MarkerGraphVertexId, vertex_descriptor>& vertexMap)
{
    PathFiller2& graph = *this;
    OrientedReadInfo& info = orientedReadInfos[i];
    const OrientedReadId orientedReadId = info.orientedReadId;

    // Find the marker graph vertex that contains this ordinal.
    const MarkerGraphVertexId vertexId =
        assembler.getGlobalMarkerGraphVertex(orientedReadId, uint32_t(ordinal));
    SHASTA_ASSERT(vertexId != invalid<MarkerGraphVertexId>);

    // Get the corresponding vertex descriptor, creating the vertex
    // if necessary.
    vertex_descriptor v;
    auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        v = add_vertex(PathFiller2Vertex(vertexId, orientedReadInfos.size()), graph);
        vertexMap.insert(make_pair(vertexId, v));
    } else {
        v = it->second;
    }

    // Store this ordinal in the vertex.
    PathFiller2Vertex& vertex = graph[v];
    vertex.ordinals[i].push_back(ordinal);

    // Store this vertex in the OrientedReadInfo.
    info.vertices.push_back(v);
}



PathFiller2Vertex::PathFiller2Vertex(
    MarkerGraphVertexId vertexId,
    uint64_t orientedReadCount) :
    vertexId(vertexId),
    ordinals(orientedReadCount)
{
}



void PathFiller2::removeLowCoverageVertices(uint64_t minVertexCoverage)
{
    PathFiller2& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        auto& vertex = graph[v];
        if((vertex.vertexId == vertexIdA) or (vertex.vertexId == vertexIdB)) {
            continue;
        }
        if(graph[v].coverage() >= minVertexCoverage) {
            continue;
        }
        verticesToBeRemoved.push_back(v);
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        removeVertex(v);
    }

    if(html and options.showDebugInformation) {
        html << "<br>After removing vertices with coverage less than " <<
            minVertexCoverage <<
            ", there are " << num_vertices(graph) << " vertices." << endl;
    }
}



void PathFiller2::removeVertex(vertex_descriptor v)
{
    PathFiller2& graph = *this;

    // Before removing the vertex, we have to record this fact
    // in the oriented reads that visit this vertex.
    const PathFiller2Vertex& vertex = graph[v];
    for(uint64_t i=0; i<vertex.ordinals.size(); i++) {
        OrientedReadInfo& info = orientedReadInfos[i];
        for(const uint64_t ordinal: vertex.ordinals[i]) {
            info.vertexAtOrdinal(ordinal) = null_vertex();
        }
    }

    // Now we can remove the vertex.
    clear_vertex(v, graph);
    remove_vertex(v, graph);

}



// Return the total number of ordinals.
uint64_t PathFiller2Vertex::coverage() const
{
    uint64_t c = 0;
    for(const auto& v: ordinals) {
        c += v.size();
    }
    return c;
}



// Return the total number of marker intervals.
uint64_t PathFiller2Edge::coverage() const
{
    uint64_t c = 0;
    for(const auto& v: markerIntervals) {
        c += v.size();
    }
    return c;
}



void PathFiller2::splitVertices(int64_t maxBaseSkip)
{

    PathFiller2& graph = *this;

    // Gather the vertices we have now so we can safely iterate over them.
    vector<vertex_descriptor> initialVertices;
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        initialVertices.push_back(v);
    }

    class OrdinalInfo {
    public:
        uint64_t i;
        int64_t ordinal;
        int64_t estimatedOffset;
        bool operator<(const OrdinalInfo& that) const
        {
            return estimatedOffset < that.estimatedOffset;
        }
    };



    // Loop over our initial vertices.
    vector<OrdinalInfo> ordinalInfos;
    vector<int64_t> splitPositions;
    uint64_t splitCount = 0;
    for(const vertex_descriptor v: initialVertices) {
        const PathFiller2Vertex& vertex = graph[v];
        SHASTA_ASSERT(vertex.ordinals.size() == orientedReadInfos.size());

        // Don't split vertexIdA or vertexIdB
        const MarkerGraphVertexId vertexId = vertex.vertexId;
        if( vertexId == vertexIdA or
            vertexId == vertexIdB) {
            continue;
        }

        // Loop over all ordinals of all oriented reads in this vertex.
        ordinalInfos.clear();
        for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
            const OrientedReadInfo& info = orientedReadInfos[i];
            const OrientedReadId orientedReadId = info.orientedReadId;

            for(const int64_t ordinal: vertex.ordinals[i]) {
                // Get the position of this marker.
                const MarkerId markerId =assembler.getMarkerId(orientedReadId, uint32_t(ordinal));
                const uint64_t position = assembler.markers.begin()[markerId].position;

                OrdinalInfo ordinalInfo;
                ordinalInfo.i = i;
                ordinalInfo.ordinal = ordinal;

                // Estimate the base offset from A0.
                if(info.isOnA() and info.isOnB()) {

                    // Base offsets from A0 and to B1.
                    const int64_t offsetFromA = position - info.ordinalAndPositionA.position;
                    const int64_t offsetToB = info.ordinalAndPositionB.position - position;

                    // Average both estimates.
                    ordinalInfo.estimatedOffset = (offsetFromA + estimatedABOffset - offsetToB) / 2;

                } else if(info.isOnA() and not info.isOnB()) {

                    const int64_t offsetFromA = position - info.ordinalAndPositionA.position;
                    ordinalInfo.estimatedOffset = offsetFromA;

                } else if(info.isOnB() and not info.isOnA()) {

                    const int64_t offsetToB = info.ordinalAndPositionB.position - position;
                    ordinalInfo.estimatedOffset = estimatedABOffset - offsetToB;

                } else {
                    SHASTA_ASSERT(0);
                }

                // Store this OrdinalInfo.
                ordinalInfos.push_back(ordinalInfo);
            }
        }

        // Sort them by estimated offset.
        sort(ordinalInfos.begin(), ordinalInfos.end());

        /*
        cout << vertex.vertexId << ":";
        for(const auto& info: ordinalInfos) {
            cout << " " << info.estimatedOffset;
        }
        cout << endl;
        */

       // Look for skips larger than maxBaseSkip.
        splitPositions.clear();
        splitPositions.push_back(0);
        for(uint64_t i=1; i<ordinalInfos.size(); i++) {
            const int64_t baseSkip = ordinalInfos[i].estimatedOffset - ordinalInfos[i-1].estimatedOffset;
            if(baseSkip > maxBaseSkip) {
                splitPositions.push_back(i);
            }
        }
        splitPositions.push_back(ordinalInfos.size());

        if(splitPositions.size() == 2) {
            // cout << "No skips found." << endl;
            continue;
        }
        ++splitCount;
        // cout << "Found " << splitPositions.size() - 2 << " skips." << endl;
        // cout << "This vertex will be split in " << splitPositions.size() - 1 << "." << endl;

        // Create the new vertices.
        for(uint64_t i=0; i<splitPositions.size()-1; i++) {
            const vertex_descriptor vNew = add_vertex(graph);
            PathFiller2Vertex& vertexNew = graph[vNew];
            vertexNew.vertexId = vertex.vertexId;
            vertexNew.replicaIndex = i + 1;
            vertexNew.ordinals.resize(orientedReadInfos.size());

            // Loop over the ordinals that go in this new vertex.
            const uint64_t begin = splitPositions[i];
            const uint64_t end = splitPositions[i+1];
            for(uint64_t j=begin; j!=end; j++) {
                const OrdinalInfo& ordinalInfo = ordinalInfos[j];
                const uint64_t i = ordinalInfo.i;
                vertexNew.ordinals[i].push_back(ordinalInfo.ordinal);

                // Store this vertex in the OrientedReadInfo.
                OrientedReadInfo& info = orientedReadInfos[i];
                info.vertexAtOrdinal(ordinalInfo.ordinal) = vNew;
            }
        }

        // Remove the initial vertex.
        remove_vertex(v, graph);
    }

    if(html and options.showDebugInformation) {
        html << "<br>" << splitCount << " vertices were split into one or more new vertices." << endl;
        html << "<br>After splitting vertices, there are " <<
            num_vertices(graph) << " vertices." << endl;
    }
}



void PathFiller2::createEdges()
{
    PathFiller2& graph = *this;

    // Loop over oriented reads.
    // Follow each oriented read over the vertices it visits.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        vertex_descriptor vPrevious = null_vertex();
        int64_t ordinalPrevious = invalid<uint32_t>;

        // Loop over the vertices visited by this oriented read.
        for(int64_t ordinal=info.firstOrdinal; ordinal<=info.lastOrdinal; ordinal++) {
            const vertex_descriptor v = info.vertexAtOrdinal(ordinal);

            // No vertex. Skip.
            if(v == null_vertex()) {
                continue;
            }

            // We found a vertex, and we don't have a previous vertex.
            // Store this vertex as the previous vertex.
            if(vPrevious == null_vertex()) {
                vPrevious = v;
                ordinalPrevious = ordinal;
                continue;
            }

            // We found a vertex, and we also have a previous vertex.
            // Add an edge between these two vertices if necessary,
            // then add a MarkerInterval to it to describe the
            // transition of this oriented read from vPrevious to v.
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(vPrevious, v, graph);
            if(not edgeExists) {
                bool edgeWasAdded = false;
                tie(e, edgeWasAdded) = add_edge(vPrevious, v, graph);
                SHASTA_ASSERT(edgeWasAdded);
                graph[e].markerIntervals.resize(orientedReadInfos.size());
            }
            graph[e].markerIntervals[i].push_back({ordinalPrevious, ordinal});

            // Replace the previous vertex with this one.
            vPrevious = v;
            ordinalPrevious = ordinal;

        }
    }

    if(html and options.showDebugInformation) {
        html << "<br>There are " << num_edges(graph) << " edges." << endl;
    }
}



// Approximate topological sort is only used to improve the
// display of the graph and make it faster.
// It sets the vertex ranks and isDagEdge flags in the edges.
void PathFiller2::approximateTopologicalSort()
{
    PathFiller2& graph = *this;

    vector< pair<edge_descriptor, uint64_t> > edgesWithCoverage;
    BGL_FORALL_EDGES(e, graph, PathFiller2) {
        edgesWithCoverage.push_back({e, graph[e].coverage()});
    }
    sort(edgesWithCoverage.begin(), edgesWithCoverage.end(),
        OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());

    vector<edge_descriptor> sortedEdges;
    for(const auto& p: edgesWithCoverage) {
        sortedEdges.push_back(p.first);
    }

    shasta::approximateTopologicalSort(graph, sortedEdges);
}



void PathFiller2::writeGraph(const string& title)
{
    PathFiller2& graph = *this;

    if(html and options.showGraph) {
        html << "<h2>" << title << "</h2>";
        html << "<p>The assembly graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges.";
        approximateTopologicalSort();
        writeGraph();
    }
}



void PathFiller2::writeGraph() const
{
    // Write out the graph in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    {
        ofstream dotFile(dotFileName);
        writeGraphviz(dotFile);
    }

    // Compute layout in svg format.
    const string command = "dot -O -T svg " + dotFileName;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    const double timeout = 600;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    if(returnCode!=0 or signalOccurred) {
        throw runtime_error("An error occurred while running the following command: " + command);
    }
    if(timeoutTriggered) {
        std::filesystem::remove(dotFileName);
        throw runtime_error("Timeout during graph layout computation.");
    }

    // Remove the .dot file.
    std::filesystem::remove(dotFileName);

    // Copy the svg file to html.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<p>" << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);
}



void PathFiller2::writeGraphviz(ostream& out) const
{
    const PathFiller2& graph = *this;

    // S and V for edges HSV.
    const double S = 0.7;
    const double V = 1.;

    // Gather assembly path edges.
    vector<edge_descriptor> sortedAssemblyPathEdges = assemblyPath;
    sort(sortedAssemblyPathEdges.begin(), sortedAssemblyPathEdges.end());


    out <<
        "digraph PathFiller2 {\n"
        "mclimit=0.01;\n"       // For layout speed
        "edge [penwidth=6];\n"
        "node [fontname=\"Courier New\"];\n"
        "edge [fontname=\"Courier New\"];\n";

    if(options.showVertices) {
        if(options.showVertexLabels) {
            out << "node [shape=rectangle];\n";
        } else {
            out << "node [shape=point width=0.2];\n";
        }
    } else {
        out << "node [shape=point style=invis];\n";
    }


    // To help Graphviz compute the layout, write vertices in rank order.
    vector< pair<vertex_descriptor, uint64_t> > verticesWithRank;
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        verticesWithRank.push_back({v, graph[v].rank});
    }
    sort(verticesWithRank.begin(), verticesWithRank.end(),
        OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());
    for(const auto& p: verticesWithRank) {
        const vertex_descriptor v = p.first;
        const PathFiller2Vertex& vertex = graph[v];
        out << "\"" << vertex.stringId() << "\" [";
        out << "tooltip=\"" << vertex.stringId() << "\\n" << vertex.coverage() << "\"";
        if(options.showVertexLabels) {
            out << " label=\"" << vertex.stringId() << "\\n" << vertex.coverage() << "\"";
        }
        if(vertex.vertexId == vertexIdA) {
            out << "style=filled fillcolor=lightgreen";
        }
        if(vertex.vertexId == vertexIdB) {
            out << "style=filled fillcolor=pink";
        }
        out << "];\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, graph, PathFiller2) {

        const PathFiller2Edge& edge = graph[e];
        const uint64_t coverage = edge.coverage();
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);

        // Compute the hue based on coverage.
        double H;
        if(coverage >= orientedReadInfos.size()) {
            H = 1./3.;
        } else {
            H = (double(coverage - 1) / (3. * double(orientedReadInfos.size() - 1)));
        }
        const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";

        out << "\"" << graph[v0].stringId() << "\"->\"" << graph[v1].stringId() << "\" [";

        // Color is based on coverage.
        out << " color=" << colorString;

        // Tooltip.
        out << " tooltip=\"";
        out << "Coverage " << coverage << "\\n";
        out << "\"";

        // Label.
        if(options.showEdgeLabels) {
            out << " label=\"";
            out << coverage;
            /*
            if(not edge.sequence.empty()) {
                out << "\\n";
                copy(edge.sequence.begin(), edge.sequence.end(), ostream_iterator<Base>(out));
            }
            */
            out << "\"";
        }

        if(not edge.isDagEdge) {
            out << " constraint=false";
        }

        // If we have an assembly path and this edge is not on the assembly path,
        // draw it dashed.
        if(not assemblyPath.empty()) {
            if(not std::binary_search(sortedAssemblyPathEdges.begin(), sortedAssemblyPathEdges.end(), e)) {
                out << " style=dashed";
            }
        }

        out << "];\n";

    }

    out << "}\n";
}



string PathFiller2Vertex::stringId() const
{
    string s = to_string(vertexId);
    if(replicaIndex) {
        s += "." + to_string(replicaIndex);
    }
    return s;
}



void PathFiller2::vertexCoverageHistogram() const
{
    if(not html) {
        return;
    }

    const PathFiller2& graph = *this;
    vector<uint64_t> histogram;
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        const uint64_t coverage = graph[v].coverage();

        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }

    html <<
        "<h2>Vertex coverage histogram</h2>"
        "<table>"
        "<tr><th>Coverage<th>Frequency";
    for(uint64_t coverage=0; coverage< histogram.size(); coverage++) {
        const uint64_t frequency = histogram[coverage];
        if(frequency) {
            html <<
                "<tr>"
                "<td class=centered>" << coverage <<
                "<td class=centered>" << frequency;
        }
    }
    html << "</table>";
}



void PathFiller2::edgeCoverageHistogram() const
{
    if(not html) {
        return;
    }

    const PathFiller2& graph = *this;
    vector<uint64_t> histogram;
    BGL_FORALL_EDGES(e, graph, PathFiller2) {
        const uint64_t coverage = graph[e].coverage();

        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }

    html <<
        "<h2>Edge coverage histogram</h2>"
        "<table>"
        "<tr><th>Coverage<th>Frequency";
    for(uint64_t coverage=0; coverage< histogram.size(); coverage++) {
        const uint64_t frequency = histogram[coverage];
        if(frequency) {
            html <<
                "<tr>"
                "<td class=centered>" << coverage <<
                "<td class=centered>" << frequency;
        }
    }
    html << "</table>";
}



void PathFiller2::removeAllEdges()
{
    PathFiller2& graph = *this;
    vector<edge_descriptor> allEdges;
    BGL_FORALL_EDGES(e, graph, PathFiller2) {
        allEdges.push_back(e);
    }

    for(const edge_descriptor e: allEdges) {
        boost::remove_edge(e, graph);
    }
}



void PathFiller2::removeStrongComponents()
{
    PathFiller2& graph = *this;
    uint64_t removedCount = 0;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        vertexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        graph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }



    // Keep the non-trivial ones.
    // A non-trivial strong component has at least one internal edge.
    // This means that it either has more than one vertex,
    // or it consists of a single vertex with a self-edge.
    for(const auto& p: componentVertices) {

        // Figure out if it is non-trivial.
        bool isNonTrivial;
        if(p.second.size() > 1) {

            // More than one vertex. Certainly non-trivial.
            isNonTrivial = true;
        } else if (p.second.size() == 1) {

            // Only one vertex. Non-trivial if self-edge present.
            const vertex_descriptor v = p.second.front();
            bool selfEdgeExists = false;
            tie(ignore, selfEdgeExists) = edge(v, v, graph);
            isNonTrivial = selfEdgeExists;
        } else {

            // Empty. This should never happen.
            SHASTA_ASSERT(0);
        }

        // If non-trivial, remove all of its vertices.
        // But don't remove vertexIdA or vertexIdB.
        if(isNonTrivial) {
            for(const vertex_descriptor v: p.second) {
                const PathFiller2Vertex& vertex = graph[v];
                if(vertex.vertexId == vertexIdA or vertex.vertexId == vertexIdB) {
                    continue;
                }
                 removeVertex(v);
                ++removedCount;
            }
        }
    }

    if(html and options.showDebugInformation) {
        html <<
            "<br>Removed " << removedCount <<
            " vertices in non-trivial strongly connected components."
            "<br>The graph has now " << num_vertices(graph) <<
            " vertices.";

    }
}



void PathFiller2::findAssemblyPath()
{
    const PathFiller2& graph = *this;
    assemblyPath.clear();


    // Find the first and last vertex of the path we are looking for.
    vertex_descriptor vA = null_vertex();
    vertex_descriptor vB = null_vertex();
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        const PathFiller2Vertex& vertex = graph[v];
        if(vertex.vertexId == vertexIdA) {
            SHASTA_ASSERT(vertex.replicaIndex == 0);
            SHASTA_ASSERT(vA == null_vertex());
            vA = v;
        }
        if(vertex.vertexId == vertexIdB) {
            SHASTA_ASSERT(vertex.replicaIndex == 0);
            SHASTA_ASSERT(vB == null_vertex());
            vB = v;
        }
    }
    SHASTA_ASSERT(vA != null_vertex());
    SHASTA_ASSERT(vB != null_vertex());


    // Main iteration loop.
    vertex_descriptor v = vA;
    while(v != vB) {

        // Find the edge with the most coverage.
        edge_descriptor eNext;
        uint64_t bestCoverage = 0;
        BGL_FORALL_OUTEDGES(v, e, graph, PathFiller2) {
            const uint64_t coverage = graph[e].coverage();
            if(coverage > bestCoverage) {
                eNext = e;
                bestCoverage = coverage;
            }
        }
        if(bestCoverage == 0) {
            cout << "PathFiller2: at " << graph[v].stringId() <<
                ": no out-edge found when filling path from " <<
                edgeIdA << " to " << edgeIdB << endl;
        }
        SHASTA_ASSERT(bestCoverage > 0);

        // Store this edge.
        assemblyPath.push_back(eNext);
        v = target(eNext, graph);
    }

    if(html and options.showDebugInformation) {
        html << "<br>The assembly path has " << assemblyPath.size() << " edges.";
    }
}



bool PathFiller2::assembleAssemblyPathEdges(uint64_t maxMsaLength)
{
    PathFiller2& graph = *this;

    for(const edge_descriptor e: assemblyPath) {

        if(html and options.showAssemblyDetails) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            html << "<h2>Assembly details for edge " <<
                graph[v0].stringId() << " to " << graph[v1].stringId() << "</h2>";
        }

        if(not assembleEdge(e, maxMsaLength)) {
            return false;
        }

    }



    if(html and options.showAssemblyDetails) {
        html <<
            "<h2>Assembly path</h2>"
            "<table>"
            "<tr><th>Vertex0<th>Vertex1<th>Begin<th>End<th>Length<th class=left>Sequence";

        uint64_t position = 0;
        for(const edge_descriptor e: assemblyPath) {
            const auto& edge = graph[e];
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);

            const uint64_t length = edge.sequence.size();
            const uint64_t end = position + length;
            html <<
                "<tr>"
                "<td class=centered>" << graph[v0].stringId() <<
                "<td class=centered>" << graph[v1].stringId() <<
                "<td class=centered>" << position <<
                "<td class=centered>" << end <<
                "<td class=centered>" << length <<
                "<td style='font-family:monospace'>";
            copy(edge.sequence.begin(), edge.sequence.end(), ostream_iterator<Base>(html));

            position += length;
        }

        html << "</table>";
    }

    return true;
}



bool PathFiller2::assembleEdge(edge_descriptor e, uint64_t maxMsaLength)
{
    const bool debug = false;
    PathFiller2& graph = *this;
    PathFiller2Edge& edge = graph[e];

    const uint64_t k = assembler.assemblerInfo->k;
    SHASTA_ASSERT((k % 2) == 0);
    const uint64_t kHalf = k / 2;

    if(debug) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        cout << "PathFiller2::assembleEdge begins for edge " <<
            graph[v0].stringId() << " to " << graph[v1].stringId() << endl;
    }


    // Gather the sequences of the contributing oriented reads.
    // Each sequence is stored with the number of distinct reads that
    // have that sequence.
    vector< pair<vector<Base>, uint64_t> > orientedReadSequences;

    // Loop over oriented reads.
    vector<Base> orientedReadSequence;
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadInfos[i].orientedReadId;
        const vector< pair<int64_t, int64_t> >& orientedReadMarkerIntervals = edge.markerIntervals[i];
        if(orientedReadMarkerIntervals.empty()) {
            continue;
        }

        // Find the first and last ordinals of this oriented read on this edge.
        const int64_t ordinal0 = orientedReadMarkerIntervals.front().first;
        const int64_t ordinal1 = orientedReadMarkerIntervals.back().second;

        // Find the corresponding MarkerIds.
        const MarkerId markerId0 = assembler.getMarkerId(orientedReadId, uint32_t(ordinal0));
        const MarkerId markerId1 = assembler.getMarkerId(orientedReadId, uint32_t(ordinal1));

        // And the corresponding positions in the oriented reads.
        const uint64_t position0 = assembler.markers.begin()[markerId0].position + kHalf;
        const uint64_t position1 = assembler.markers.begin()[markerId1].position + kHalf;

        // Now we can get the sequence contributed by this oriented read.
        orientedReadSequence.clear();
        for(uint64_t position=position0; position!=position1; position++) {
            const Base base = assembler.getReads().getOrientedReadBase(orientedReadId, uint32_t(position));
            orientedReadSequence.push_back(base);
        }

        if(debug) {
            copy(orientedReadSequence.begin(), orientedReadSequence.end(),
                ostream_iterator<Base>(cout));
            cout << " " << orientedReadId << endl;
        }

        // Store it.
        bool found = false;
        for(auto& p: orientedReadSequences) {
            if(p.first == orientedReadSequence) {
                ++p.second;
                found = true;
                break;
            }
        }
        if(not found) {
            orientedReadSequences.push_back(make_pair(orientedReadSequence, 1));
        }
    }

    // Do the MSA.
    vector<Base> consensus;
    if(not globalMsaSpoa(orientedReadSequences, consensus, maxMsaLength)) {
        if(html) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            html << "<br>MSA is too long for edge " <<
                graph[v0].stringId() << " to " << graph[v1].stringId() << endl;
        }
        return false;
    }
    if(debug) {
        copy(consensus.begin(), consensus.end(), ostream_iterator<Base>(cout));
        cout << " Consensus" << endl;
    }

    // Store the consensus in the edge.
    edge.sequence = consensus;

    if(debug) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        cout << "PathFiller1::assembleEdge ends for edge " <<
            graph[v0].stringId() << " " << graph[v1].stringId() << endl;
    }

    return true;
}



// Get the sequence between edgeIdA and edgeIdB.
// This does not include the sequences od edgeIdA and edgeIdB themselves.
void PathFiller2::getSecondarySequence(
    vector<Base>& sequence) const
{
    const PathFiller2& graph = *this;

    sequence.clear();
    for(const edge_descriptor e: assemblyPath) {
        const vector<Base>& edgeSequence = graph[e].sequence;
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
    }
}


// Get the complete sequence, including the sequences of edgeIdA and edgeIdB.
void PathFiller2::getCompleteSequence(vector<Base>& sequence) const
{
    SHASTA_ASSERT(assemblyPath.size() >= 2);
    const PathFiller2& graph = *this;

    sequence.clear();

    const auto edgeASequence = assembler.markerGraph.edgeSequence[edgeIdA];
    copy(edgeASequence.begin(), edgeASequence.end(), back_inserter(sequence));

    for(const edge_descriptor e: assemblyPath) {
        const vector<Base>& edgeSequence = graph[e].sequence;
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
    }

    const auto edgeBSequence = assembler.markerGraph.edgeSequence[edgeIdB];
    copy(edgeBSequence.begin(), edgeBSequence.end(), back_inserter(sequence));

}



// Remove vertices that are not accessible from vertexIdA
// or from which vertexIdB is not accessible.
// Returns the number of vertices that were removed.
uint64_t PathFiller2::removeInaccessibleVertices()
{
    PathFiller2& graph = *this;

    // Find the vertices corresponding to vertexIdA and vertexIdB.
    vertex_descriptor vA = null_vertex();
    vertex_descriptor vB = null_vertex();
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        const PathFiller2Vertex& vertex = graph[v];
        if(vertex.vertexId == vertexIdA) {
            SHASTA_ASSERT(vertex.replicaIndex == 0);
            SHASTA_ASSERT(vA == null_vertex());
            vA = v;
        }
        if(vertex.vertexId == vertexIdB) {
            SHASTA_ASSERT(vertex.replicaIndex == 0);
            SHASTA_ASSERT(vB == null_vertex());
            vB = v;
        }
    }
    SHASTA_ASSERT(vA != null_vertex());
    SHASTA_ASSERT(vB != null_vertex());



    // Use a forward BFS to find the vertices that are accessible from vertexIdA,
    // moving forward. Those vertices get their isAccessibleA flag set.
    {
        std::queue<vertex_descriptor> q;
        q.push(vA);
        graph[vA].isAccessibleA = true;
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();

            BGL_FORALL_OUTEDGES(v0, e, graph, PathFiller2) {
                const vertex_descriptor v1 = target(e, graph);
                auto& vertex1 = graph[v1];
                if(not vertex1.isAccessibleA) {
                    vertex1.isAccessibleA = true;
                    q.push(v1);
                }
            }
        }
        SHASTA_ASSERT(graph[vB].isAccessibleA);
    }



    // Use a backward BFS to find the vertices that are accessible from vertexIdB,
    // moving backward. Those vertices get their isAccessibleB flag set.
    {
        std::queue<vertex_descriptor> q;
        q.push(vB);
        graph[vB].isAccessibleB = true;
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();

            BGL_FORALL_INEDGES(v0, e, graph, PathFiller2) {
                const vertex_descriptor v1 = source(e, graph);
                auto& vertex1 = graph[v1];
                if(not vertex1.isAccessibleB) {
                    vertex1.isAccessibleB = true;
                    q.push(v1);
                }
            }
        }
        SHASTA_ASSERT(graph[vA].isAccessibleB);
    }


    // Gather the vertices to be removed.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        const auto& vertex = graph[v];
        if(not (vertex.isAccessibleA and vertex.isAccessibleB)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    // Remove them.
    for(const vertex_descriptor v: verticesToBeRemoved) {
        removeVertex(v);
    }

    return verticesToBeRemoved.size();
}



void PathFiller2::writeVerticesCsv() const
{
    const PathFiller2& graph = *this;

    ofstream csv("PathFiller2-Vertices.csv");
    csv << "Vertex,OrientedRead,Ordinal\n";

    BGL_FORALL_VERTICES(v, graph, PathFiller2) {
        const PathFiller2Vertex& vertex = graph[v];
        const string vertexStringId = vertex.stringId();

        SHASTA_ASSERT(vertex.ordinals.size() == orientedReadInfos.size());
        for(uint64_t i=0; i<vertex.ordinals.size(); i++) {
            const OrientedReadId orientedReadId = orientedReadInfos[i].orientedReadId;
            const vector<int64_t>& orientedReadOrdinals = vertex.ordinals[i];

            for(const int64_t ordinal: orientedReadOrdinals) {
                csv << vertexStringId << ",";
                csv << orientedReadId << ",";
                csv << ordinal << "\n";
            }
        }

    }
}

#include "mode3b-PathFiller2.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3b;

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



    checkAssumptions();
    gatherOrientedReads();
    writeOrientedReads();
    estimateOffset();
    createVertices(estimatedOffsetRatio);
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
    const MarkerGraph::Edge& edgeA = markerGraph.edges[edgeIdA];
    const MarkerGraph::Edge& edgeB = markerGraph.edges[edgeIdB];
    const MarkerGraphVertexId vertexIdA0 = edgeA.source;
    const MarkerGraphVertexId vertexIdA1 = edgeA.target;
    const MarkerGraphVertexId vertexIdB0 = edgeB.source;
    const MarkerGraphVertexId vertexIdB1 = edgeB.target;
    if(markerGraph.vertexHasDuplicateOrientedReadIds(vertexIdA0, markers)) {
        throw runtime_error("Duplicated oriented read on source vertex of edgeIdA.");
    }
    if(markerGraph.vertexHasDuplicateOrientedReadIds(vertexIdA1, markers)) {
        throw runtime_error("Duplicated oriented read on target vertex of edgeIdA.");
    }
    if(markerGraph.vertexHasDuplicateOrientedReadIds(vertexIdB0, markers)) {
        throw runtime_error("Duplicated oriented read on source vertex of edgeIdB.");
    }
    if(markerGraph.vertexHasDuplicateOrientedReadIds(vertexIdB1, markers)) {
        throw runtime_error("Duplicated oriented read on target vertex of edgeIdB.");
    }
}



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

        const uint32_t ordinal0 = markerInterval.ordinals[0];
        const MarkerId markerId0 = assembler.getMarkerId(orientedReadId, ordinal0);
        info.ordinalAndPositionA0 = {ordinal0, assembler.markers.begin()[markerId0].position};

        const uint32_t ordinal1 = markerInterval.ordinals[1];
        const MarkerId markerId1 = assembler.getMarkerId(orientedReadId, ordinal1);
        info.ordinalAndPositionA1 = {ordinal1, assembler.markers.begin()[markerId1].position};

        SHASTA_ASSERT(ordinal1 == ordinal0 + 1);
    }

    // Fill in the OrdinalAndPositions of the oriented reads that appear in edgeIdB.
    for(const MarkerInterval& markerInterval: assembler.markerGraph.edgeMarkerIntervals[edgeIdB]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const uint64_t i = getOrientedReadIndex(orientedReadId);
        OrientedReadInfo& info = orientedReadInfos[i];

        const uint32_t ordinal0 = markerInterval.ordinals[0];
        const MarkerId markerId0 = assembler.getMarkerId(orientedReadId, ordinal0);
        info.ordinalAndPositionB0 = {ordinal0, assembler.markers.begin()[markerId0].position};

        const uint32_t ordinal1 = markerInterval.ordinals[1];
        const MarkerId markerId1 = assembler.getMarkerId(orientedReadId, ordinal1);
        info.ordinalAndPositionB1 = {ordinal1, assembler.markers.begin()[markerId1].position};

        SHASTA_ASSERT(ordinal1 == ordinal0 + 1);
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
        "<th>OrdinalA0"
        "<th>OrdinalA1"
        "<th>OrdinalB0"
        "<th>OrdinalB1"
        "<th>Ordinal<br>offset"
        "<th>PositionA0"
        "<th>PositionA1"
        "<th>PositionB0"
        "<th>PositionB1"
        "<th>Position<br>offset0"
        "<th>Position<br>offset1"
        ;

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered>" << info.orientedReadId;

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalAndPositionA0.ordinal;
        }

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalAndPositionA1.ordinal;
        }
        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalAndPositionB0.ordinal;
        }

        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalAndPositionB1.ordinal;
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            html << info.ordinalOffset();
        }

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalAndPositionA0.position;
        }

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalAndPositionA1.position;
        }
        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalAndPositionB0.position;
        }

        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalAndPositionB1.position;
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            html << info.positionOffset0();
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            html << info.positionOffset1();
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
            sum += info.positionOffset0();
            sum += info.positionOffset1();
        }
    }
    estimatedOffset = sum / (2*n);

    if(html) {
        html << "<p>Estimated offset is " << estimatedOffset << " bases.";
    }
}



// Create vertices by following the portion of each oriented read we will use for assembly.
void PathFiller2::createVertices(double estimatedOffsetRatio)
{
    PathFiller2& graph = *this;
    const int64_t offsetThreshold = int64_t(estimatedOffsetRatio * double(estimatedOffset));

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
            info.firstOrdinal = info.ordinalAndPositionA0.ordinal;
            info.lastOrdinal =  info.ordinalAndPositionB1.ordinal;
            for(int64_t ordinal=info.firstOrdinal; ordinal<=info.lastOrdinal; ordinal++) {
                createVerticesHelper(i, ordinal, vertexMap);
            }
        }

        // Oriented reads that appear on edgeIdA but not on edgeIdB.
        else if(info.isOnA() and not info.isOnB()) {
            info.firstOrdinal = info.ordinalAndPositionA0.ordinal;
            const int64_t maxPosition = info.ordinalAndPositionA1.position + offsetThreshold;
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
            info.lastOrdinal = info.ordinalAndPositionB1.ordinal;
            const int64_t minPosition = info.ordinalAndPositionB0.position - offsetThreshold;
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

    if(html) {
        html << "<p>Found " << num_vertices(graph) << " vertices." << endl;
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

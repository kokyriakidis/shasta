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
    checkAssumptions();
    gatherOrientedReads();
    writeOrientedReads();

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
        "<p>"
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


// Shasta.
#include "mode3b-PathFiller3.hpp"
#include "Assembler.hpp"
#include "MarkerGraph.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3b;



PathFiller3::PathFiller3(
    const Assembler& assembler,
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    const PathFiller3DisplayOptions& options) :
    assembler(assembler),
    edgeIdA(edgeIdA),
    edgeIdB(edgeIdB),
    options(options),
    html(options.html)
{

    // Store the source target of edgeIdA and the source vertex of edgeIdB.
    const MarkerGraph::Edge& edgeA = assembler.markerGraph.edges[edgeIdA];
    const MarkerGraph::Edge& edgeB = assembler.markerGraph.edges[edgeIdB];
    vertexIdA = edgeA.target;
    vertexIdB = edgeB.source;

    // Check assumptions here as this used vertexIdA and vertexIdB.
    checkAssumptions();

    // Oriented reads.
    gatherOrientedReads();
    writeOrientedReads();

    // Use the oriented reads that appear both on vertexIdA and vertexIdB
    // to estimate the base offset between vertexIdA and vertexIdB.
    estimateOffset();

}



void PathFiller3::checkAssumptions() const
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



void PathFiller3::gatherOrientedReads()
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

            OrientedReadInfo info(orientedReadIdA);
            info.ordinalA = ordinalA;
            orientedReadInfos.push_back(info);

            ++itA;
        }



        // Oriented reads that appear only in edgeIdB.
        else if((itA == endA) or (itB != endB and itB->orientedReadId < itA->orientedReadId)) {

            const MarkerInterval& markerIntervalB = *itB;
            const OrientedReadId orientedReadIdB = markerIntervalB.orientedReadId;
            const uint32_t ordinalB = markerIntervalB.ordinals[0];    // Because vertexIdB is the source of edgeIdB

            OrientedReadInfo info(orientedReadIdB);
            info.ordinalB = ordinalB;
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

                OrientedReadInfo info(orientedReadId);
                info.ordinalA = ordinalA;
                info.ordinalB = ordinalB;
                orientedReadInfos.push_back(info);
            }

            ++itA;
            ++itB;
        }

    }
}



void PathFiller3::writeOrientedReads() const
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
        "<th>Ordinal<br>offset"
        ;

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered>" << info.orientedReadId;

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << info.ordinalA;
        }

        html << "<td class=centered>";
        if(info.isOnB()) {
            html << info.ordinalB;
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            html << info.ordinalOffset();
        }

        html << "<td class=centered>";
        if(info.isOnA()) {
            html << basePosition(info.orientedReadId, info.ordinalA);
        }

        html << "<td class=centered>";
        if(info.isOnB()) {
            html << basePosition(info.orientedReadId, info.ordinalB);
        }

        html << "<td class=centered>";
        if(info.isOnA() and info.isOnB()) {
            const int64_t baseOffset =
                basePosition(info.orientedReadId, info.ordinalB) -
                basePosition(info.orientedReadId, info.ordinalA);
            SHASTA_ASSERT(baseOffset >= 0);
            html << baseOffset;
        }
     }

    html << "</table>";
}



// Get the base position of a marker in an oriented read
// given the ordinal.
int64_t PathFiller3::basePosition(OrientedReadId orientedReadId, int64_t ordinal) const
{
    const MarkerId markerId = assembler.getMarkerId(orientedReadId, uint32_t(ordinal));
    const int64_t position = int64_t(assembler.markers.begin()[markerId].position);
    return position;

}



void PathFiller3::estimateOffset()
{
    int64_t n = 0;
    int64_t sum = 0;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        if(info.isOnA() and info.isOnB()) {
            const OrientedReadId orientedReadId = info.orientedReadId;
            const int64_t positionA = basePosition(orientedReadId, info.ordinalA);
            const int64_t positionB = basePosition(orientedReadId, info.ordinalB);
            const int64_t baseOffset = positionB - positionA;
            SHASTA_ASSERT(baseOffset >= 0);

            sum += baseOffset;
            ++n;
        }
    }
    estimatedABOffset = int64_t(std::round(double(sum) / double(n)));

    if(html) {
        html << "<br>Estimated offset is " << estimatedABOffset << " bases.";
    }
}





// Shasta.
#include "mode3b-PathFiller3.hpp"
#include "Assembler.hpp"
#include "markerAccessFunctions.hpp"
#include "MarkerGraph.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3b;

// Seqan.
#include <seqan/align.h>

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"



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

    // PARAMETERS THAT SHOULD BE EXPOSED WHEN CODE STABILIZES.

    // The estimated offset gets extended by this ratio to
    // decide how much to extend reads that only appear in edgeIdA or edgeIdB.
    double estimatedOffsetRatio = 1.1;

    // Alignment parameters.
    int64_t matchScore = 6;
    int64_t mismatchScore = -1;
    int64_t gapScore = -1;

    const uint64_t minVertexCoverage = 12;


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

    // Markers.
    gatherMarkers(estimatedOffsetRatio);
    writeMarkers();

    // Marker graph.
    alignAndDisjointSets(matchScore, mismatchScore, gapScore);
    createVertices(minVertexCoverage);
    createEdges();
    writeGraph("Initial assembly graph");
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



// Fill in the markerInfos vector of each read.
void PathFiller3::gatherMarkers(double estimatedOffsetRatio)
{
    const int64_t offsetThreshold = int64_t(estimatedOffsetRatio * double(estimatedABOffset));


    // Loop over our oriented reads.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        OrientedReadInfo& info = orientedReadInfos[i];
        const OrientedReadId orientedReadId = info.orientedReadId;
        info.markerInfos.clear();

        // Oriented reads that appear on both edgeIdA and edgeIdB.
        if(info.isOnA() and info.isOnB()) {
            for(int64_t ordinal=info.ordinalA;
                ordinal<=info.ordinalB; ordinal++) {
                addMarkerInfo(i, ordinal);
            }
        }

        // Oriented reads that appear on edgeIdA but not on edgeIdB.
        else if(info.isOnA() and not info.isOnB()) {
            const int64_t maxPosition = basePosition(orientedReadId, info.ordinalA) + offsetThreshold;
            const int64_t markerCount = int64_t(assembler.markers.size(orientedReadId.getValue()));

            for(int64_t ordinal=info.ordinalA;
                ordinal<markerCount; ordinal++) {
                const int64_t position = basePosition(orientedReadId, ordinal);
                if(position > maxPosition) {
                    break;
                }
                addMarkerInfo(i, ordinal);
            }
        }

        // Oriented reads that appear on edgeIdB but not on edgeIdA.
        else if(info.isOnB() and not info.isOnA()) {
            const int64_t minPosition = basePosition(orientedReadId, info.ordinalB) - offsetThreshold;

            for(int64_t ordinal=info.ordinalB; ordinal>=0; ordinal--) {
                const int64_t position = basePosition(orientedReadId, ordinal);
                if(position < minPosition) {
                    break;
                }
                addMarkerInfo(i, ordinal);
            }

            // We added the MarkerInfos in reverse order, so we have to reverse them.
            reverse(info.markerInfos.begin(), info.markerInfos.end());
        }

        else {
            SHASTA_ASSERT(0);
        }
    }

}



// Add the marker at given ordinal to the i-th oriented read.
void PathFiller3::addMarkerInfo(uint64_t i, int64_t ordinal)
{
    OrientedReadInfo& info = orientedReadInfos[i];

    MarkerInfo markerInfo;
    markerInfo.ordinal = ordinal;
    markerInfo.position = basePosition(info.orientedReadId, ordinal);
    markerInfo.kmerId = getOrientedReadMarkerKmerId(
        info.orientedReadId,
        uint32_t(ordinal),
        assembler.assemblerInfo->k,
        assembler.getReads(),
        assembler.markers);

    info.markerInfos.push_back(markerInfo);
}



void PathFiller3::writeMarkers()
{
    if(not (html and options.showDebugInformation)) {
        return;
    }

    const uint64_t k = assembler.assemblerInfo->k;

    html <<
        "<h2>Markers used in this assembly step</h2>"
        "<table>"
        "<tr>"
        "<th>Oriented<br>read<br>index"
        "<th>Oriented<br>read"
        "<th>Ordinal"
        "<th>Position"
        "<th>KmerId"
        "<th>Kmer";

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        for(const MarkerInfo& markerInfo: info.markerInfos) {
            const Kmer kmer(markerInfo.kmerId, k);

            html <<
                "<tr>"
                "<td class=centered>" << i <<
                "<td class=centered>" << info.orientedReadId <<
                "<td class=centered>" << markerInfo.ordinal <<
                "<td class=centered>" << markerInfo.position <<
                "<td class=centered>" << markerInfo.kmerId <<
                "<td class=centered style='font-family:monospace'>";
            kmer.write(html, k);
        }
    }

    html << "</table>";
}



// Compute alignments and use them to create the disjoint set data structure,
// from which the marker graph will be created.
void PathFiller3::alignAndDisjointSets(
    uint64_t matchScore,
    uint64_t mismatchScore,
    uint64_t gapScore
    )
{

    // SeqAn types we need.
    using TSequence = seqan::String<KmerId>;
    using TStringSet = seqan::StringSet<TSequence>;
    using TDepStringSet = seqan::StringSet< TSequence, seqan::Dependent<> >;
    using TAlignGraph = seqan::Graph< seqan::Alignment<TDepStringSet> >;

    // Assign ids to markers.
    uint64_t markerId = 0;
    for(OrientedReadInfo& info: orientedReadInfos) {
        for(MarkerInfo& markerInfo: info.markerInfos) {
            markerInfo.id = markerId++;
        }
    }

    // Initialize the disjoint sets data structure.
    const uint64_t markerCount = markerId;
    vector<uint64_t> rank(markerCount);
    vector<uint64_t> parent(markerCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t markerId=0; markerId<markerCount; markerId++) {
        disjointSets.make_set(markerId);
    }

    // Construct a Seqan sequence containing the KmerIds for each oriented read.
    // Add 100 to each KmerId because Seqan uses 45 to represent a gap.
    vector<TSequence> seqanSequences(orientedReadInfos.size());
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        TSequence& seqanSequence = seqanSequences[i];
        for(const MarkerInfo& markerInfo: info.markerInfos) {
            seqan::appendValue(seqanSequence, markerInfo.kmerId + 100);
        }
    }



    // Loop over pairs of reads.
    for(uint64_t i0=0; i0<orientedReadInfos.size()-1; i0++) {
        const OrientedReadInfo& info0 = orientedReadInfos[i0];
        const uint64_t length0 = info0.markerInfos.size();
        const TSequence& seqanSequence0 = seqanSequences[i0];
        for(uint64_t i1=i0+1; i1<orientedReadInfos.size(); i1++) {
            const OrientedReadInfo& info1 = orientedReadInfos[i1];
            const uint64_t length1 = info1.markerInfos.size();
            const TSequence& seqanSequence1 = seqanSequences[i1];

            // Figure the constraints for this alignment.
            const bool constrainedA = info0.isOnA() and info1.isOnA();
            const bool constrainedB = info0.isOnB() and info1.isOnB();

            // Only do alignments that are constrained on at least one side.
            if(not (constrainedA or constrainedB)) {
                continue;
            }

            // Align the KmerIds of these oriented reads.
            // For now we do a full blown alignment, but later
            // we should use banded alignments instead.
            // Store them in a SeqAn string set.
            TStringSet sequences;
            appendValue(sequences, seqanSequence0);
            appendValue(sequences, seqanSequence1);

            // Compute the alignment.
            using namespace seqan;
            TAlignGraph graph(sequences);
            if(constrainedA and constrainedB) {
                globalAlignment(
                    graph,
                    Score<int, Simple>(int(matchScore), int(mismatchScore), int(gapScore)),
                    AlignConfig<false, false, false, false>(),
                    LinearGaps());
            } else  if(constrainedA and not constrainedB) {
                globalAlignment(
                    graph,
                    Score<int, Simple>(int(matchScore), int(mismatchScore), int(gapScore)),
                    AlignConfig<false, false, true, true>(),
                    LinearGaps());
            } else  if(constrainedB and not constrainedA) {
                globalAlignment(
                    graph,
                    Score<int, Simple>(int(matchScore), int(mismatchScore), int(gapScore)),
                    AlignConfig<true, true, false, false>(),
                    LinearGaps());
            } else {
                SHASTA_ASSERT(0);
            }

            // Extract the alignment from the graph.
            // This creates a single sequence consisting of the two rows
            // of the alignment, concatenated.
            TSequence align;
            convertAlignment(graph, align);
            const uint64_t totalAlignmentLength = seqan::length(align);
            SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
            const uint64_t alignmentLength = totalAlignmentLength / 2;

            // Use the alignment to update the disjoint sets data structure.
            uint64_t j0 = 0;
            uint64_t j1 = 0;
            const uint64_t seqanGapValue = 45;
            for(uint64_t positionInAlignment=0; positionInAlignment<alignmentLength; positionInAlignment++) {
                const KmerId kmerId0 = align[positionInAlignment];
                const KmerId kmerId1 = align[positionInAlignment + alignmentLength];

                if(kmerId0 == seqanGapValue) {
                    if(kmerId1 == seqanGapValue) {
                        // Both 0 and 1 are gaps.
                        SHASTA_ASSERT(0);
                    } else {
                        // 0 is gap, 1 is not gap.
                        ++j1;
                    }
                } else {
                    if(kmerId1 == seqanGapValue) {
                        // 0 is not gap, 1 is gap.
                        ++j0;
                    } else {
                        // Neither 0 nor 1 is a gap.
                        if(kmerId0 == kmerId1) {
                            // If a match, merge the disjoint sets containing these two markers.
                            disjointSets.union_set(info0.markerInfos[j0].id, info1.markerInfos[j1].id);
                        }
                        ++j0;
                        ++j1;
                    }

                }
            }
            SHASTA_ASSERT(j0 == length0);
            SHASTA_ASSERT(j1 == length1);
       }
    }

    // Store in each MarkerInfo the id of the disjoint set it belongs to.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        OrientedReadInfo& info = orientedReadInfos[i];
        for(MarkerInfo& markerInfo: info.markerInfos) {
            markerInfo.disjointSetId = disjointSets.find_set(markerInfo.id);
        }
    }

    // Fill in the disjoint sets map.
    disjointSetsMap.clear();
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        for(uint64_t j=0; j<info.markerInfos.size(); j++) {
            const MarkerInfo& markerInfo = info.markerInfos[j];
            disjointSetsMap[markerInfo.disjointSetId].push_back({i, j});
        }
    }

    // Write a histogram of disjoint sets sizes.
    if(html and options.showDebugInformation) {
        vector<uint64_t> histogram;
        for(const auto& p: disjointSetsMap) {
            const uint64_t disjointSetSize = p.second.size();
            if(disjointSetSize >= histogram.size()) {
                histogram.resize(disjointSetSize + 1, 0);
            }
            ++histogram[disjointSetSize];
        }

        html <<
            "<h2>Disjoint sets size histogram</h2>"
            "<table>"
            "<tr>"
            "<th>Size"
            "<th>Frequency"
            "<th>Markers";

        for(uint64_t disjointSetSize=0; disjointSetSize<histogram.size(); disjointSetSize++) {
            const uint64_t frequency = histogram[disjointSetSize];
            if(frequency) {
                html <<
                    "<tr>"
                    "<td class=centered>" << disjointSetSize <<
                    "<td class=centered>" << frequency <<
                    "<td class=centered>" << frequency * disjointSetSize;
            }
        }

        html << "</table>";
    }

}



// Create vertices. Each disjoint set with at least minVertexCoverage markers
// generates a vertex.
void PathFiller3::createVertices(uint64_t minVertexCoverage)
{
    PathFiller3& graph = *this;

    // Remove all vertices and edges, just in case.
    PathFiller3BaseClass::clear();
    vertexMap.clear();

    // Find the disjoint sets corresponding to vertexIdA and vertexIdB.
    // Those will always generate a vertex regardless of coverage.
    disjointSetIdA = invalid<uint64_t>;
    disjointSetIdB = invalid<uint64_t>;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        if(info.isOnA()) {
            const MarkerInfo& markerInfoA = info.markerInfos.front();
            if(disjointSetIdA == invalid<uint64_t>) {
                disjointSetIdA = markerInfoA.disjointSetId;
            } else {
                SHASTA_ASSERT(disjointSetIdA = markerInfoA.disjointSetId);
            }
        }
        if(info.isOnB()) {
            const MarkerInfo& markerInfoB = info.markerInfos.back();
            if(disjointSetIdB == invalid<uint64_t>) {
                disjointSetIdB = markerInfoB.disjointSetId;
            } else {
                SHASTA_ASSERT(disjointSetIdB = markerInfoB.disjointSetId);
            }
        }
    }

    // Loop over disjoint sets that are large enough.
    for(const auto& p: disjointSetsMap) {
        const auto& disjointSet = p.second;
        if(disjointSet.size() < minVertexCoverage) {
            continue;
        }

        const uint64_t disjointSetId = p.first;
        const vertex_descriptor v = add_vertex({disjointSetId}, graph);
        vertexMap.insert(make_pair(disjointSetId, v));
    }

    if(html and options.showDebugInformation) {
        html << "<br>The assembly graph has " << num_vertices(graph) << " vertices.";
    }
}



// Create edges by following the reads.
void PathFiller3::createEdges()
{
    PathFiller3& graph = *this;

    removeAllEdges();

    // Loop over all reads.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        // Follow this read, finding the vertices it reaches.
        vertex_descriptor v0 = null_vertex();
        PathFiller3MarkerIndexes indexes0;
        for(uint64_t j=0; j<info.markerInfos.size(); j++) {
            const MarkerInfo& markerInfo = info.markerInfos[j];
            const uint64_t disjointSetId = markerInfo.disjointSetId;
            const auto it = vertexMap.find(disjointSetId);

            if(it != vertexMap.end()) {
                const vertex_descriptor v1 = it->second;
                const PathFiller3MarkerIndexes indexes1 = {i, j};
                if(v0 != null_vertex()) {

                    // Get the edge v0->v1, creating it if necessary.
                    edge_descriptor e;
                    bool edgeExists = false;
                    tie(e, edgeExists) = edge(v0, v1, graph);
                    if(not edgeExists) {
                        bool edgeWasAdded = false;
                        tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
                        SHASTA_ASSERT(edgeWasAdded);
                    }
                    PathFiller3Edge& edge = graph[e];

                    edge.markerIntervals.push_back({indexes0, indexes1});
                }

                // v1 becomes the previous vertex.
                v0 = v1;
                indexes0 = indexes1;

            }
        }
    }
    if(html and options.showDebugInformation) {
        html << "<br>The assembly graph has " << num_edges(graph) << " edges.";
    }
}



void PathFiller3::removeAllEdges()
{
    PathFiller3& graph = *this;
    BGL_FORALL_VERTICES(v, graph, PathFiller3) {
        clear_vertex(v, graph);
    }
}



void PathFiller3::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void PathFiller3::writeGraphviz(ostream& s) const
{
    const PathFiller3& graph = *this;

    // S and V for edges HSV.
    const double S = 0.7;
    const double V = 1.;

    s <<
        "digraph PathFiller3 {\n"
        "mclimit=0.01;\n"       // For layout speed
        "edge [penwidth=6];\n"
        "node [fontname=\"Courier New\"];\n"
        "edge [fontname=\"Courier New\"];\n";

    if(options.showVertices) {
        if(options.showVertexLabels) {
            s << "node [shape=rectangle style=filled color=black fillcolor=gray80];\n";
        } else {
            s << "node [shape=point width=0.2];\n";
        }
    } else {
        s << "node [shape=point style=invis];\n";
    }

    // Vertices.
    BGL_FORALL_VERTICES(v, graph, PathFiller3) {
        const uint64_t disjointSetId = graph[v].disjointSetId;
        auto it = disjointSetsMap.find(disjointSetId);
        SHASTA_ASSERT(it != disjointSetsMap.end());
        const uint64_t coverage = it->second.size();

        const bool isA = (graph[v].disjointSetId == disjointSetIdA);
        const bool isB = (graph[v].disjointSetId == disjointSetIdB);

        s << disjointSetId << "[";

        // Label.
        s << "label=\"";
        if(isA) {
            s << "A\\n";
        }
        if(isB) {
            s << "B\\n";
        }
        s << graph[v].disjointSetId << "\\n" << coverage;
        s << "\" ";

        // Special drawing of the begin/end vertices.
        if(isA or isB) {
            s << "shape=rectangle style=filled color=black fillcolor=cyan";
        }

        s << "];\n";
    }

    // Edges.
    BGL_FORALL_EDGES(e, graph, PathFiller3) {
        const PathFiller3Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const uint64_t coverage = edge.markerIntervals.size();

        // Compute the hue based on coverage.
        double H;
        if(coverage >= orientedReadInfos.size()) {
            H = 1./3.;
        } else {
            H = (double(coverage - 1) / (3. * double(orientedReadInfos.size() - 1)));
        }
        const string colorString = "\"" + to_string(H) + " " + to_string(S) + " " + to_string(V) + "\"";

        s <<
            graph[v0].disjointSetId << "->" <<
            graph[v1].disjointSetId << " [";

        if(options.showEdgeLabels) {
            s << "label=\"" << coverage << "\"";
        }
        s << " color=" << colorString;

        // Tooltip.
        s << " tooltip=\"";
        s << "Coverage " << coverage << "\\n";
        s << "\"";

        s << "];\n";
    }

    s << "}\n";
}



void PathFiller3::writeGraph(const string& title)
{
    PathFiller3& graph = *this;

    if(html and options.showGraph) {
        html << "<h2>" << title << "</h2>";
        html << "<p>The assembly graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges.";
        writeGraph();
    }
}



void PathFiller3::writeGraph() const
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


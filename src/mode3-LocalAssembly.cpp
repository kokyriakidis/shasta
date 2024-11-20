// Shasta.
#include "mode3-LocalAssembly.hpp"
#include "globalMsa.hpp"
#include "Marker.hpp"
#include "markerAccessFunctions.hpp"
#include "mode3-Anchor.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Seqan.
#include <seqan/align.h>
namespace seqan = seqan2;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"



// The oriented reads common between edgeIdA and edgeIdB are always
// used for assembly. The oriented reads that appear only
// on edgeIdA or edgeIdB are used for assembly under control
// of useA and useB.
// So, if useA and useB are both true (the default), the assembly uses the
// union of the oriented reads on edgeIdA and edgeIdB.
// If they are both false, the assembly uses the
// intersection of the oriented reads on edgeIdA and edgeIdB.
// If useA is true and useB is false, the assembly uses the
// oriented reads on edgeIdA, regardless of whether they appear on edgeIdB.
// If useA is false and useB is true, the assembly uses the
// oriented reads on edgeIdB, regardless of whether they appear on edgeIdA.
LocalAssembly::LocalAssembly(
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    uint64_t minVertexCoverage, // 0 = automatic
    const LocalAssemblyDisplayOptions& displayOptions,
    const Mode3AssemblyOptions::LocalAssemblyOptions& options,
    bool useA,
    bool useB) :
    k(k),
    reads(reads),
    markers(markers),
    anchors(anchors),
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB),
    options(displayOptions),
    html(displayOptions.html)
{

    SHASTA_ASSERT((k % 2) == 0);
    kHalf = k / 2;


    // If the edges are adjacent, stop here, leaving the AssemblyPath empty.
    // This results in empty secondary sequence.
    if(anchors.areAdjacentAnchors(anchorIdA, anchorIdB)) {
        if(html) {
            html << "<br>The two anchors are adjacent. Intervening sequence is empty.";
        }
        return;
    }

    // Check assumptions here as this used vertexIdA and vertexIdB.
    checkAssumptions();

    // Oriented reads.
    gatherOrientedReads(useA, useB);

    // Use the oriented reads that appear both on vertexIdA and vertexIdB
    // to estimate the base offset between vertexIdA and vertexIdB.
    estimateOffset();

    // If the offset is negative or cannot be estimated, stop here.
    // This is pathological and results in empty assembled sequence.
    if((estimatedABOffset == invalid<int64_t>) or (estimatedABOffset <= 0)) {
        if(html) {
            html << "<br>The estimated offset is not positive." << endl;
        }
        throw std::runtime_error("The estimated offset is not positive.");
    }

    // Markers.
    gatherMarkers(options.estimatedOffsetRatio);
    writeOrientedReads();
    writeOrientedReadsSequences();

    // Assembly graph.
    alignAndDisjointSets(
        options.matchScore, options.mismatchScore, options.gapScore,
        options.maxSkipBases,
        options.maxDrift, options.minHalfBand, options.minScoreRatio);
    writeMarkers();

    // Iteration to reduce minVertexCoverage if a long MSA is encountered.
    while(true) {

        minVertexCoverage = createVertices(minVertexCoverage, options.vertexSamplingRate);
        createEdges();
        writeGraph("Initial assembly graph");

        // Remove strongly connected components, then regenerate
        // edges from scratch with the remaining vertices.
        if(removeStrongComponents() > 0) {
            removeAllEdges();
            createEdges();
            writeGraph("Assembly graph after removal of strong connected components");
        }

        // This must be done after removing strongly connected components.
        // Otherwise we can have inaccessible vertices that cause the
        // assembly path to encounter dead ends.
        if(removeInaccessibleVertices()) {
            writeGraph("Assembly graph after removal of inaccessible vertices.");
        }

        // Assemble.
        findAssemblyPath();
        if(minVertexCoverage > 2) {
            try {
                assembleAssemblyPathEdges(options.maxMsaLength, LongMsaPolicy::throwException);
            } catch(...) {
                --minVertexCoverage;
                clear();
                if(html and displayOptions.showDebugInformation) {
                    html << "<br>minVertexCoverage reduced to " << minVertexCoverage;
                }
                continue;
            }

        } else {
            assembleAssemblyPathEdges(options.maxMsaLength, LongMsaPolicy::assembleAtLowCoverage);
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

            ofstream fasta("LocalAssembly.fasta");
            fasta << ">LocalAssembly " << sequence.size() << endl;
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));

            getCompleteSequence(sequence);

            html <<
                "Assembled sequence including the first and last edge is " <<
                sequence.size() << " bases long."
                "<pre style='font-family:monospace'>\n";
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));
            html << "</pre>";
        }

        break;
    }

}



void LocalAssembly::checkAssumptions() const
{
    SHASTA_ASSERT(anchorIdA != anchorIdB);
    SHASTA_ASSERT(reads.representation == 0);
}



void LocalAssembly::gatherOrientedReads(bool useA, bool useB)
{
    // Joint loop over marker intervals that appear in edgeIdA and/or edgeIdB.
    const auto markerIntervalsA = anchors[anchorIdA];
    const auto markerIntervalsB = anchors[anchorIdB];
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

            if(useA) {
                const AnchorMarkerInterval& markerIntervalA = *itA;
                const OrientedReadId orientedReadIdA = markerIntervalA.orientedReadId;
                const uint32_t ordinalA = markerIntervalA.ordinal0 + anchors.ordinalOffset(anchorIdA);    // Because vertexIdA is the target of anchorIdA

                OrientedReadInfo info(orientedReadIdA);
                info.ordinalA = ordinalA;
                orientedReadInfos.push_back(info);
            }

            ++itA;
        }



        // Oriented reads that appear only in edgeIdB.
        else if((itA == endA) or (itB != endB and itB->orientedReadId < itA->orientedReadId)) {

            if(useB) {
                const AnchorMarkerInterval& markerIntervalB = *itB;
                const OrientedReadId orientedReadIdB = markerIntervalB.orientedReadId;
                const uint32_t ordinalB = markerIntervalB.ordinal0;    // Because vertexIdB is the source of anchorIdB

                OrientedReadInfo info(orientedReadIdB);
                info.ordinalB = ordinalB;
                orientedReadInfos.push_back(info);
            }

            ++itB;
        }



        // Oriented reads that appear in both edgeIdA and edgeIdB.
        // They are always used for assembly regardless of the settings of useA and useB.
        else {
            SHASTA_ASSERT(itA != endA);
            SHASTA_ASSERT(itB != endB);

            const AnchorMarkerInterval& markerIntervalA = *itA;
            const OrientedReadId orientedReadIdA = markerIntervalA.orientedReadId;

            const AnchorMarkerInterval& markerIntervalB = *itB;
            const OrientedReadId orientedReadIdB = markerIntervalB.orientedReadId;

            SHASTA_ASSERT(orientedReadIdA == orientedReadIdB);
            const OrientedReadId orientedReadId = orientedReadIdA;

            const uint32_t ordinalA = markerIntervalA.ordinal0 + anchors.ordinalOffset(anchorIdA);    // Because vertexIdA is the target of anchorIdA
            const uint32_t ordinalB = markerIntervalB.ordinal0;    // Because vertexIdB is the source of anchorIdB

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



void LocalAssembly::writeOrientedReads() const
{
    if(not html) {
        return;
    }
    if(not options.showOrientedReads) {
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
        "<th>First<br>ordinal"
        "<th>Last<br>ordinal"
        "<th>First<br>position"
        "<th>Last<br>position"
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

        SHASTA_ASSERT(not info.markerInfos.empty());
        const MarkerInfo& firstMarkerInfo = info.markerInfos.front();
        const MarkerInfo& lastMarkerInfo = info.markerInfos.back();
        html <<
            "<td class=centered>" << firstMarkerInfo.ordinal <<
            "<td class=centered>" << lastMarkerInfo.ordinal <<
            "<td class=centered>" << firstMarkerInfo.position <<
            "<td class=centered>" << lastMarkerInfo.position;
     }

    html << "</table>";


    // Count reads.
    uint64_t commonCount = 0;
    uint64_t onlyACount = 0;
    uint64_t onlyBCount = 0;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        const bool isOnA = info.isOnA();
        const bool isOnB = info.isOnB();
        if(isOnA) {
            if(isOnB) {
                ++commonCount;
            } else {
                ++onlyACount;
            }
        } else {
            if(isOnB) {
                ++onlyBCount;
            } else {
                SHASTA_ASSERT(0);
            }

        }
    }
    html <<
        "<p><table>"
        "<tr><th class=left>Common<td class=centered>" << commonCount <<
        "<tr><th class=left>On A only<td class=centered>" << onlyACount <<
        "<tr><th class=left>On B only<td class=centered>" << onlyBCount <<
        "<tr><th class=left>Total<td class=centered>" << orientedReadInfos.size() <<
        "</table>";

}



// Get the base position of a marker in an oriented read
// given the ordinal.
int64_t LocalAssembly::basePosition(OrientedReadId orientedReadId, int64_t ordinal) const
{
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];
    const int64_t position = int64_t(orientedReadMarkers[ordinal].position);
    return position;
}



void LocalAssembly::estimateOffset()
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
    if(n == 0) {
        estimatedABOffset = invalid<int64_t>;

        if(html) {
            html << "<br>The offset cannot be estimated because there are no common oriented reads between " <<
                anchorIdA << " and " << anchorIdB;
        }
    } else {
        estimatedABOffset = int64_t(std::round(double(sum) / double(n)));

        if(html) {
            html << "<br>Estimated position offset is " << estimatedABOffset << " bases.";
        }
    }
}



// Fill in the markerInfos vector of each read.
void LocalAssembly::gatherMarkers(double estimatedOffsetRatio)
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
            const int64_t markerCount = int64_t(markers.size(orientedReadId.getValue()));

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
void LocalAssembly::addMarkerInfo(uint64_t i, int64_t ordinal)
{
    OrientedReadInfo& info = orientedReadInfos[i];

    MarkerInfo markerInfo;
    markerInfo.ordinal = ordinal;
    markerInfo.position = basePosition(info.orientedReadId, ordinal);
    markerInfo.kmerId = getOrientedReadMarkerKmerId(
        info.orientedReadId,
        uint32_t(ordinal),
        k,
        reads,
        markers);

    info.markerInfos.push_back(markerInfo);
}



void LocalAssembly::writeMarkers()
{
    if(not (html and options.showMarkers)) {
        return;
    }

    html <<
        "<h2>Markers used in this assembly step</h2>"
        "<table>"
        "<tr>"
        "<th>Oriented<br>read<br>index"
        "<th>Oriented<br>read"
        "<th>Ordinal"
        "<th>Ordinal<br>offset<br>from A"
        "<th>Ordinal<br>offset<br>to B"
        "<th>Position"
        "<th>KmerId"
        "<th>Kmer"
        "<th>Vertex"
        "<th>Coverage";

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        for(const MarkerInfo& markerInfo: info.markerInfos) {
            const Kmer kmer(markerInfo.kmerId, k);

            html <<
                "<tr>"
                "<td class=centered>" << i <<
                "<td class=centered>" << info.orientedReadId <<
                "<td class=centered>" << markerInfo.ordinal;

            // Ordinal offset from A.
            html << "<td class=centered>";
            if(info.isOnA()) {
                html << markerInfo.ordinal - info.markerInfos.front().ordinal;
            }

            // Ordinal offset to B.
            html << "<td class=centered>";
            if(info.isOnB()) {
                html << info.markerInfos.back().ordinal - markerInfo.ordinal;
            }

            html <<
                "<td class=centered>" << markerInfo.position <<
                "<td class=centered>" << markerInfo.kmerId <<
                "<td class=centered style='font-family:monospace'>";
            kmer.write(html, k);
            html <<
                "<td class=centered>" << markerInfo.disjointSetId <<
                "<td class=centered>" << disjointSetsMap[markerInfo.disjointSetId].size();
        }
    }

    html << "</table>";
}



// Compute alignments and use them to create the disjoint set data structure,
// from which the marker graph will be created.
// maxDrift is the maximum tolerated length drift of each read.
// Used to compute the band for banded alignments.
void LocalAssembly::alignAndDisjointSets(
    uint64_t matchScore,
    uint64_t mismatchScore,
    uint64_t gapScore,
    uint64_t /* maxSkipBases */,
    double maxDrift,
    uint64_t minHalfBand,
    double minScoreRatio
    )
{

    // SeqAn types we need.
    using TSequence = seqan::String<KmerId>;
    using TStringSet = seqan::StringSet<TSequence>;
    using TDepStringSet = seqan::StringSet< TSequence, seqan::Dependent<> >;
    using TAlignGraph = seqan::Graph< seqan::Alignment<TDepStringSet> >;

    const bool detailedDebugOutput = false;
    ofstream dot;
    ofstream csv;
    if(detailedDebugOutput) {
        dot.open("LocalAssembly-AlignmentDetails.dot");
        dot << "graph PathFiler3lignments {\n";
        csv.open("LocalAssembly-AlignmentDetails.csv");
        csv << "OrientedReadId0,Ordinal0,OrientedReadId1,Ordinal1\n";
    }

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
            // cout << "*** " << info0.orientedReadId << " " << info1.orientedReadId << endl;
            const uint64_t length1 = info1.markerInfos.size();
            const TSequence& seqanSequence1 = seqanSequences[i1];

            // Figure the constraints for this alignment.
            const bool constrainedA = info0.isOnA() and info1.isOnA();
            const bool constrainedB = info0.isOnB() and info1.isOnB();

            // If constrained on A, merge the first markers of the two reads,
            // as the alignment does not guarantee that.
            // If constrained on B, merge the last markers of the two reads,
            // as the alignment does not guarantee that.
            if(constrainedA) {
                disjointSets.union_set(info0.markerInfos.front().id, info1.markerInfos.front().id);
            }
            if(constrainedB) {
                disjointSets.union_set(info0.markerInfos.back().id, info1.markerInfos.back().id);
            }

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


#if 0
            // Old code that used geneal alignment, not banded alignments.
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
#endif


            // Banded alignment, allowing for the specified maxDrift.
            // This is necessary to prevent large cycles in the graph.
            // It is also good for performance.
            using namespace seqan;
            TAlignGraph graph(sequences);
            int score = 0;
            if(constrainedA and constrainedB) {
                const int64_t diagonalA = 0;
                const int64_t diagonalB = int64_t(length0) - int64_t(length1);
                const int64_t totalDrift = int64_t(maxDrift * 0.5 * double(min(length0, length1)));
                const int64_t halfBand = totalDrift + int64_t(minHalfBand);
                const int64_t minBand = min(diagonalA, diagonalB) - halfBand;
                const int64_t maxBand = max(diagonalA, diagonalB) + halfBand;
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(int(matchScore), int(mismatchScore), int(gapScore)),
                    AlignConfig<false, false, false, false>(),
                    int(minBand), int(maxBand),
                    LinearGaps());
            } else  if(constrainedA and not constrainedB) {
                const int64_t diagonalA = 0;
                const int64_t totalDrift = int64_t(maxDrift * double(min(length0, length1)));
                const int64_t halfBand = totalDrift + int64_t(minHalfBand);
                const int64_t minBand = diagonalA - halfBand;
                const int64_t maxBand = diagonalA + halfBand;
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(int(matchScore), int(mismatchScore), int(gapScore)),
                    AlignConfig<false, false, true, true>(),
                    int(minBand), int(maxBand),
                    LinearGaps());
            } else  if(constrainedB and not constrainedA) {
                const int64_t diagonalB = int64_t(length0) - int64_t(length1);
                const int64_t totalDrift = int64_t(maxDrift * double(min(length0, length1)));
                const int64_t halfBand = totalDrift + int64_t(minHalfBand);
                const int64_t minBand = diagonalB - halfBand;
                const int64_t maxBand = diagonalB + halfBand;
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(int(matchScore), int(mismatchScore), int(gapScore)),
                    AlignConfig<true, true, false, false>(),
                    int(minBand), int(maxBand),
                    LinearGaps());
            } else {
                SHASTA_ASSERT(0);
            }

            // If SeqAn was not able to compute the banded aignment, ignore it.
            if(score == MinValue<int>::VALUE) {
                if(html and options.showDebugInformation) {
                    html << "<br>Alignment between " << info0.orientedReadId <<
                        " and " << info1.orientedReadId <<
                        " ignored.";
                }
                continue;
            }

            // Check that the score is sufficiently good.
            const uint64_t bestPossibleScore = matchScore * min(length0, length1);
            const double scoreRatio = double(score) / double(bestPossibleScore);
            if(scoreRatio < minScoreRatio) {
                if(html and options.showDebugInformation) {
                    html << "<br>Alignment between " << info0.orientedReadId <<
                        " and " << info1.orientedReadId << ": lengths " << length0 << " " << length1 <<
                        ", score " << score << "/" << bestPossibleScore << " " <<
                        double(score) / double(bestPossibleScore) <<
                        " discarded due to low score.";
                }
                continue;
            }



            // Extract the alignment from the graph.
            // This creates a single sequence consisting of the two rows
            // of the alignment, concatenated.
            TSequence align;
            convertAlignment(graph, align);
            const uint64_t totalAlignmentLength = seqan::length(align);
            SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
            const uint64_t alignmentLength = totalAlignmentLength / 2;
            const uint64_t seqanGapValue = 45;

#if 0
            // This is not needed when doing banded alignments.
            // If the alignment has large base skips, don't use it.
            bool hasLargeSkip = false;
            uint64_t j0 = 0;
            uint64_t j1 = 0;
            uint64_t previousPosition0 = invalid<uint64_t>;
            uint64_t previousPosition1 = invalid<uint64_t>;
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
                            // Check for large base skips.
                            const uint64_t position0 = info0.markerInfos[j0].position;
                            const uint64_t position1 = info1.markerInfos[j1].position;
                            // cout << "***A " << position0 << " " << position1 << endl;
                            if(previousPosition0 != invalid<uint64_t>) {
                                const int64_t offset = int64_t(position0) - int64_t(position1);
                                const int64_t previousOffset = int64_t(previousPosition0) - int64_t(previousPosition1);
                                if(abs(offset - previousOffset) > int64_t(maxSkipBases)) {
                                    hasLargeSkip = true;
                                    // cout << "Skip" << endl;
                                    break;
                                }
                            }
                            previousPosition0 = position0;
                            previousPosition1 = position1;
                        }
                        ++j0;
                        ++j1;
                    }

                }
            }
            if(hasLargeSkip) {
                if(html and options.showDebugInformation) {
                    html << "<br>Alignment between " << info0.orientedReadId <<
                        " and " << info1.orientedReadId <<
                        " suppressed.";
                }
                continue;
            }
#endif


            // Use the alignment to update the disjoint sets data structure.
            uint64_t j0 = 0;
            uint64_t j1 = 0;
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
                            if(detailedDebugOutput) {
                                dot << "\"" << info0.orientedReadId << "-";
                                dot << info0.markerInfos[j0].ordinal << "\"--\"";
                                dot << info1.orientedReadId << "-";
                                dot << info1.markerInfos[j1].ordinal << "\";\n";
                                csv <<
                                    info0.orientedReadId << "," <<
                                    info0.markerInfos[j0].ordinal << "," <<
                                    info1.orientedReadId << "," <<
                                    info1.markerInfos[j1].ordinal << "\n";
                            }
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

    // Histogram disjoint sets sizes.
    disjointSetsSizeHistogram.clear();
    for(const auto& p: disjointSetsMap) {
        const uint64_t disjointSetSize = p.second.size();
        if(disjointSetSize >= disjointSetsSizeHistogram.size()) {
            disjointSetsSizeHistogram.resize(disjointSetSize + 1, 0);
        }
        ++disjointSetsSizeHistogram[disjointSetSize];
    }


    // Write the histogram of disjoint sets sizes.
    if(html and options.showDebugInformation) {

        html <<
            "<h2>Disjoint sets size histogram</h2>"
            "<table>"
            "<tr>"
            "<th>Size"
            "<th>Frequency"
            "<th>Markers";

        for(uint64_t disjointSetSize=0; disjointSetSize<disjointSetsSizeHistogram.size(); disjointSetSize++) {
            const uint64_t frequency = disjointSetsSizeHistogram[disjointSetSize];
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

    if(detailedDebugOutput) {
        dot << "}\n";
    }
}



// Create vertices. Each disjoint set with at least minVertexCoverage markers
// generates a vertex.
uint64_t LocalAssembly::createVertices(
    uint64_t minVertexCoverage,
    double vertexSamplingRate)  // Only used if minVertexCoverage is 0
{
    LocalAssembly& graph = *this;

    // Remove all vertices and edges, just in case.
    LocalAssemblyBaseClass::clear();
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
                SHASTA_ASSERT(disjointSetIdA == markerInfoA.disjointSetId);
            }
        }
        if(info.isOnB()) {
            const MarkerInfo& markerInfoB = info.markerInfos.back();
            if(disjointSetIdB == invalid<uint64_t>) {
                disjointSetIdB = markerInfoB.disjointSetId;
            } else {
                SHASTA_ASSERT(disjointSetIdB == markerInfoB.disjointSetId);
            }
        }
    }

    if(html) {
        html << "<br>Start vertex " << disjointSetIdA << ", end vertex " << disjointSetIdB;
    }



    // If minVertexCoverage is 0, select a value automatically.
    // Select a value that gives a number of vertices approximately correct given
    // the estimated offset.
    if(minVertexCoverage == 0) {

        // Estimate the desired number of vertices given the estimated offset.
        const uint64_t totalBaseCount = reads.getTotalBaseCount() * 2; // Both strands.
        const uint64_t totalMarkerCount = markers.totalSize();
        const double markerDensity = double(totalMarkerCount) / double(totalBaseCount);
        const uint64_t desiredVertexCount = uint64_t(
            vertexSamplingRate *  markerDensity * double(estimatedABOffset));

        // Use the disjointSetsSizeHistogram to choose a value of minVertexCoverage
        // that will give us approximately this number of vertices.
        // Never reduce minVertexCoverage below 2.
        uint64_t cumulativeDisjointSetsCount = 0;
        for(minVertexCoverage = disjointSetsSizeHistogram.size()-1; minVertexCoverage>2; --minVertexCoverage) {
            cumulativeDisjointSetsCount += disjointSetsSizeHistogram[minVertexCoverage];
#if 0
            if(html and options.showDebugInformation) {
                html << "<br>minVertexCoverage " << minVertexCoverage <<
                    " would generate " << cumulativeDisjointSetsCount <<
                    " vertices and we want " << desiredVertexCount;
            }
#endif
            if(cumulativeDisjointSetsCount >= desiredVertexCount) {
                break;
            }
        }

        if(html and options.showDebugInformation) {
            html << "<br>Set minVertexCoverage to " << minVertexCoverage <<
                " based on marker density " << markerDensity <<
                ", vertex sampling rate " << vertexSamplingRate <<
                ", desired number of vertices " << desiredVertexCount;
        }
    }



    // Loop over disjoint sets that are large enough.
    // Also always include disjointSetIdA and disjointSetIdB.
    for(const auto& p: disjointSetsMap) {
        const uint64_t disjointSetId = p.first;
        const auto& disjointSet = p.second;
        if(disjointSet.size() >= minVertexCoverage or
            disjointSetId==disjointSetIdA or
            disjointSetId==disjointSetIdB) {

            const vertex_descriptor v = add_vertex({disjointSetId}, graph);
            vertexMap.insert(make_pair(disjointSetId, v));
        }
    }

    if(html and options.showDebugInformation) {
        html << "<br>The assembly graph has " << num_vertices(graph) << " vertices.";
    }

    return minVertexCoverage;
}



// Create edges by following the reads.
void LocalAssembly::createEdges()
{
    LocalAssembly& graph = *this;

    removeAllEdges();

    // Loop over all reads.
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];

        // Follow this read, finding the vertices it reaches.
        vertex_descriptor v0 = null_vertex();
        LocalAssemblyMarkerIndexes indexes0;
        for(uint64_t j=0; j<info.markerInfos.size(); j++) {
            const MarkerInfo& markerInfo = info.markerInfos[j];
            const uint64_t disjointSetId = markerInfo.disjointSetId;
            const auto it = vertexMap.find(disjointSetId);

            if(it != vertexMap.end()) {
                const vertex_descriptor v1 = it->second;
                const LocalAssemblyMarkerIndexes indexes1 = {i, j};
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
                    LocalAssemblyEdge& edge = graph[e];

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



void LocalAssembly::removeAllEdges()
{
    LocalAssembly& graph = *this;
    BGL_FORALL_VERTICES(v, graph, LocalAssembly) {
        clear_vertex(v, graph);
    }
}



void LocalAssembly::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void LocalAssembly::writeGraphviz(ostream& s) const
{
    const LocalAssembly& graph = *this;

    // S and V for edges HSV.
    const double S = 0.7;
    const double V = 1.;

    // Gather assembly path edges.
    vector<edge_descriptor> sortedAssemblyPathEdges = assemblyPath;
    sort(sortedAssemblyPathEdges.begin(), sortedAssemblyPathEdges.end());

    s <<
        "digraph LocalAssembly {\n"
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
    BGL_FORALL_VERTICES(v, graph, LocalAssembly) {
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
    BGL_FORALL_EDGES(e, graph, LocalAssembly) {
        const LocalAssemblyEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const uint64_t coverage = edge.coverage();

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

        // If we have an assembly path and this edge is not on the assembly path,
        // draw it dashed.
        if(not assemblyPath.empty()) {
            if(not std::binary_search(sortedAssemblyPathEdges.begin(), sortedAssemblyPathEdges.end(), e)) {
                s << " style=dashed";
            }
        }

        s << "];\n";
    }

    s << "}\n";
}



void LocalAssembly::writeGraph(const string& title)
{
    LocalAssembly& graph = *this;

    if(html and options.showGraph) {
        html << "<h2>" << title << "</h2>";
        html << "<p>The assembly graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges.";
        writeGraph();
    }
}



void LocalAssembly::writeGraph() const
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



uint64_t LocalAssembly::removeStrongComponents()
{
    LocalAssembly& graph = *this;
    uint64_t removedCount = 0;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, graph, LocalAssembly) {
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
                const LocalAssemblyVertex& vertex = graph[v];
                if(vertex.disjointSetId == disjointSetIdA or vertex.disjointSetId == disjointSetIdB) {
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

    return removedCount;
}



void LocalAssembly::removeVertex(vertex_descriptor v)
{
    LocalAssembly& graph = *this;

    vertexMap.erase(graph[v].disjointSetId);

    clear_vertex(v, graph);
    remove_vertex(v, graph);

}



void LocalAssembly::findAssemblyPath()
{
    const LocalAssembly& graph = *this;
    assemblyPath.clear();


    // Find the first and last vertex of the path we are looking for.
    vertex_descriptor vA = null_vertex();
    vertex_descriptor vB = null_vertex();
    BGL_FORALL_VERTICES(v, graph, LocalAssembly) {
        const LocalAssemblyVertex& vertex = graph[v];
        if(vertex.disjointSetId == disjointSetIdA) {
            SHASTA_ASSERT(vA == null_vertex());
            vA = v;
        }
        if(vertex.disjointSetId == disjointSetIdB) {
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
        BGL_FORALL_OUTEDGES(v, e, graph, LocalAssembly) {
            // Ignore a self-edge A->A.
            // This can exist because we did not allow vertex A (and B)
            // to be removed when removing strong components.
            if(v == vA and target(e, graph) == vA) {
                continue;
            }
            const uint64_t coverage = graph[e].coverage();
            if(coverage > bestCoverage) {
                eNext = e;
                bestCoverage = coverage;
            }
        }
        if(bestCoverage == 0) {
            cout << "LocalAssembly: at " << graph[v].disjointSetId <<
                ": no out-edge found when filling path from " <<
                anchorIdA << " to " << anchorIdB << endl;
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




void LocalAssembly::assembleAssemblyPathEdges(
    uint64_t maxMsaLength,
    LongMsaPolicy longMsaPolicy)
{
    const LocalAssembly& graph = *this;

    for(const edge_descriptor e: assemblyPath) {
        assembleEdge(maxMsaLength, longMsaPolicy, e);
    }



    // Write a table containing a summary of edge sequences with coverage,
    // and their position in assembled sequence.
    if(html and options.showAssemblyDetails) {
        html <<
            "<br><table>"
            "<tr>"
            "<th>Source"
            "<th>Target"
            "<th>Begin"
            "<th>End"
            "<th>length"
            "<th>Sequence"
            ;

        uint64_t position = 0;
        for(const edge_descriptor e: assemblyPath) {
            const LocalAssemblyEdge& edge = graph[e];
            const vector<Base>& sequence = edge.consensusSequence;
            const vector<uint64_t>& coverage = edge.consensusCoverage;
            SHASTA_ASSERT(sequence.size() == coverage.size());

            html <<
                "<tr>"
                "<td class=centered>" << graph[source(e, graph)].disjointSetId <<
                "<td class=centered>" << graph[target(e, graph)].disjointSetId <<
                "<td class=centered>" << position <<
                "<td class=centered>" << position + sequence.size() <<
                "<td class=centered>" << sequence.size() <<
                "<td class=centered style='font-family:monospace'>";
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));
            html << "<br>";
            for(const uint64_t c: coverage) {
                writeCoverageCharacterToHtml(c);
            }

            position += sequence.size();
        }

        html << "</table>";
    }
}



void LocalAssembly::assembleEdge(
    uint64_t maxMsaLength,
    LongMsaPolicy longMsaPolicy,
    edge_descriptor e)
{
    LocalAssembly& graph = *this;
    LocalAssemblyEdge& edge = graph[e];

    if(html and options.showAssemblyDetails) {
        html << "<h2>Assembly details for edge " <<
            graph[source(e, graph)].disjointSetId << "->" <<
            graph[target(e, graph)].disjointSetId << "</h2>"
            "<table>"
            "<tr><th>Oriented<br>read<th>Sequence<br>length<th>Sequence";
    }

    // Gather the sequences of the contributing oriented reads.
    // Each sequence is stored with the number of distinct oriented reads that
    // have that sequence.
    vector< pair<vector<Base>, uint64_t> > orientedReadSequences;

    // Loop over marker intervals of this edge.
    vector<Base> orientedReadSequence;
    for(const auto& p: edge.markerIntervals) {

        // Locate the two markers of this marker interval.
        const LocalAssemblyMarkerIndexes indexes0 = p.first;
        const LocalAssemblyMarkerIndexes indexes1 = p.second;
        const uint64_t i0 = indexes0.i;
        const uint64_t i1 = indexes1.i;
        const uint64_t j0 = indexes0.j;
        const uint64_t j1 = indexes1.j;

        // They must belong to the same oriented read.
        SHASTA_ASSERT(i0 == i1);
        const uint64_t i = i0;
        const OrientedReadInfo& info = orientedReadInfos[i];
        const OrientedReadId orientedReadId = info.orientedReadId;

        const MarkerInfo& markerInfo0 = info.markerInfos[j0];
        const MarkerInfo& markerInfo1 = info.markerInfos[j1];

        // Now we can get the contributing sequence.
        const uint64_t position0 = markerInfo0.position + kHalf;
        const uint64_t position1 = markerInfo1.position + kHalf;

        // Now we can get the sequence contributed by this oriented read.
        orientedReadSequence.clear();
        for(uint64_t position=position0; position!=position1; position++) {
            const Base base = reads.getOrientedReadBase(orientedReadId, uint32_t(position));
            orientedReadSequence.push_back(base);
        }

        if(html and options.showAssemblyDetails) {
            html <<
                "<tr><td class=centered>" << orientedReadId <<
                "<td class=centered>" << orientedReadSequence.size() <<
                "<td class=centered style='font-family:monospace'>";
            copy(orientedReadSequence.begin(), orientedReadSequence.end(),
                ostream_iterator<Base>(html));
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

    // Sort the sequences by decreasing number of supporting reads.
    sort(orientedReadSequences.begin(), orientedReadSequences.end(),
        OrderPairsBySecondOnlyGreater<vector<Base>, uint64_t>());

    if(html and options.showAssemblyDetails) {
        html << "</table>";

        html << "<p><table>"
            "<tr><th>Coverage<th>Sequence<br>length<th>Sequence";
        for(const auto& p: orientedReadSequences) {
            const vector<Base>& sequence = p.first;
            const uint64_t coverage = p.second;
            html <<
                "<tr>"
                "<td class=centered>" << coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=centered style='font-family:monospace'>";
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));

        }
        html << "</table>";
    }

    // If there is only one distinct sequence (all reads agree),
    // store that one sequence as the consensus.
    // This is the most common case.
    if(orientedReadSequences.size() == 1) {
        const auto& p = orientedReadSequences.front();
        const vector<Base>& sequence = p.first;
        const uint64_t coverage = p.second;
        edge.consensusSequence = sequence;
        edge.consensusCoverage.clear();
        edge.consensusCoverage.resize(sequence.size(), coverage);
        return;
    }


    // If getting here, we have more than one sequence, and we must
    // compute a consensus via multiple sequence alignment (MSA).

    // If any of the sequences are too long, react according to longMsaPolicy.
    // This can be problematic.
    if(orientedReadSequences.size() > 1) {

        // Find the length of the longest sequence.
        uint64_t maxLength = 0;
        for(const auto& p: orientedReadSequences) {
            const vector<Base>& sequence = p.first;
            maxLength = max(sequence.size(), maxLength);
        }
        if(html and options.showDebugInformation) {
            html << "<br>Maximum sequence length " << maxLength;
            html << "<br>Maximum allowed sequence length " << maxMsaLength;
        }

        if(maxLength > maxMsaLength) {
            if(html and options.showDebugInformation) {
                html << "<br>MSA length " << maxLength << " at " <<
                    graph[source(e, graph)].disjointSetId << "->" <<
                    graph[target(e, graph)].disjointSetId;
            }
            if(longMsaPolicy == LongMsaPolicy::throwException) {
                throw runtime_error("Long MSA.");
            } else {
                orientedReadSequences.resize(1);
                if(html and options.showDebugInformation) {
                    html << "<br>Assembling this edge at coverage " << orientedReadSequences.front().second;
                }
            }
        }
    }

    // Compute the MSA.
    vector< vector<AlignedBase> > alignment;
    globalMsaSpoa(orientedReadSequences, alignment);
    SHASTA_ASSERT(alignment.size() == orientedReadSequences.size());

    // Compute coverage at each alignment position for each of the 5 AlignedBases.
    const uint64_t alignmentLength = alignment.front().size();
    vector< array<uint64_t, 5> > coverage(alignmentLength, {0, 0, 0, 0, 0});
    for(uint64_t i=0; i<orientedReadSequences.size(); i++) {
        const vector<AlignedBase>& alignmentRow = alignment[i];
        SHASTA_ASSERT(alignmentRow.size() == alignmentLength);
        for(uint64_t position=0; position<alignmentLength; position++) {
            const AlignedBase b = alignmentRow[position];
            coverage[position][b.value] += orientedReadSequences[i].second;
        }
    }

    // Compute coverage-based consensus at each alignment position.
    vector<AlignedBase> alignedConsensus;
    vector<uint64_t> alignmentConsensusCoverage;
    for(const auto& c: coverage) {
        const uint64_t iBase = std::max_element(c.begin(), c.end()) - c.begin();
        alignedConsensus.push_back(AlignedBase::fromInteger(iBase));
        alignmentConsensusCoverage.push_back(c[iBase]);
    }
    SHASTA_ASSERT(alignedConsensus.size() == alignmentLength);

    // Store in the edge the consensus and its coverage, excluding the gaps.
    edge.consensusSequence.clear();
    edge.consensusCoverage.clear();
    for(uint64_t position=0; position<alignedConsensus.size(); position++) {
        const AlignedBase b = alignedConsensus[position];
        if(not b.isGap()) {
            edge.consensusSequence.push_back(Base(b));
            edge.consensusCoverage.push_back(alignmentConsensusCoverage[position]);
        }
    }

    if(html and options.showAssemblyDetails) {

        html << "<p><table>"
            "<tr><th>Coverage<th>Sequence<br>length<th>Aligned<br>sequence";

        // Write one row for each distinct sequence.
        for(uint64_t i=0; i<orientedReadSequences.size(); i++) {
            const auto& p = orientedReadSequences[i];
            const vector<Base>& sequence = p.first;
            const uint64_t coverage = p.second;
            const vector<AlignedBase>& alignedSequence = alignment[i];
            html <<
                "<tr>"
                "<td class=centered>" << coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=centered style='font-family:monospace'>";
            for(uint64_t position=0; position<alignedSequence.size(); position++) {
                const AlignedBase b = alignedSequence[position];
                const bool isDiscordant = (b != alignedConsensus[position]);
                if(isDiscordant) {
                    html << "<span style='background-color:LightCoral'>";
                }
                html << alignedSequence[position];
                if(isDiscordant) {
                    html << "</span>";
                }
            }
        }

        // Write one row with aligned consensus.
        html <<
            "<tr>"
            "<td class=centered colspan=2>Consensus"
            "<td class=centered style='font-family:monospace'>";
        copy(alignedConsensus.begin(), alignedConsensus.end(),
            ostream_iterator<AlignedBase>(html));

        // Write one row with aligned consensus coverage.
        html <<
            "<tr>"
            "<td class=centered colspan=2>Consensus coverage"
            "<td class=centered style='font-family:monospace'>";
        for(uint64_t position=0; position<coverage.size(); position++) {
            writeCoverageCharacterToHtml(alignmentConsensusCoverage[position]);
        }

        // Write one row with aligned discordant coverage.
        html <<
            "<tr>"
            "<td class=centered colspan=2>Discordant coverage"
            "<td class=centered style='font-family:monospace'>";
        for(uint64_t position=0; position<coverage.size(); position++) {
            writeCoverageCharacterToHtml(edge.coverage() - alignmentConsensusCoverage[position]);
        }

        // Write one row with coverage for each of the 5 AlignedBases.
        for(uint64_t b=0; b<5; b++) {
            html <<
                "<tr><td colspan=2 class=centered>" << AlignedBase::fromInteger(b) << " coverage"
                "<td class=centered style='font-family:monospace'>";
            for(uint64_t position=0; position<coverage.size(); position++) {
                writeCoverageCharacterToHtml(coverage[position][b]);
            }
        }
        html << "</table>";

        // Write another table with the final, ungapped consensus and its coverage.
        html <<
            "<p>Consensus length is " << edge.consensusSequence.size() <<
            "<br><table>"
            "<tr><th>Consensus<td class=centered style='font-family:monospace'>";
        copy(edge.consensusSequence.begin(), edge.consensusSequence.end(),
            ostream_iterator<Base>(html));
        html << "<tr><th>Consensus coverage<td class=centered style='font-family:monospace'>";
        for(const uint64_t coverage: edge.consensusCoverage) {
            writeCoverageCharacterToHtml(coverage);
        }
        html << "<tr><th>Discordant coverage<td class=centered style='font-family:monospace'>";
        for(const uint64_t coverage: edge.consensusCoverage) {
            writeCoverageCharacterToHtml(edge.coverage() - coverage);
        }
        html << "</table>";
    }

}



void LocalAssembly::writeCoverageCharacterToHtml(uint64_t coverage) const
{
    if(coverage == 0) {
        html << "&nbsp;";
    } else if(coverage < 10) {
        html << coverage;
    } else if(coverage < 36) {
        html << char((coverage - 10) + 'A');
    } else {
        html << "*";
    }

}


// Get the sequence between edgeIdA and edgeIdB.
// This does not include the sequences of edgeIdA and edgeIdB themselves.
void LocalAssembly::getSecondarySequence(
    vector<Base>& sequence) const
{
    const LocalAssembly& graph = *this;

    sequence.clear();
    for(const edge_descriptor e: assemblyPath) {
        const vector<Base>& edgeSequence = graph[e].consensusSequence;
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
    }

}



// Get the complete sequence, including the sequences of edgeIdA and edgeIdB.
void LocalAssembly::getCompleteSequence(
    vector<Base>& sequence) const
{
    const LocalAssembly& graph = *this;

    sequence.clear();

    const auto edgeASequence = anchors.anchorSequences[anchorIdA];
    copy(edgeASequence.begin(), edgeASequence.end(), back_inserter(sequence));

    for(const edge_descriptor e: assemblyPath) {
        const vector<Base>& edgeSequence = graph[e].consensusSequence;
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
    }

    const auto edgeBSequence = anchors.anchorSequences[anchorIdB];
    copy(edgeBSequence.begin(), edgeBSequence.end(), back_inserter(sequence));


}



// Remove vertices that are not accessible from vertexIdA
// or from which vertexIdB is not accessible.
// Returns the number of vertices that were removed.
uint64_t LocalAssembly::removeInaccessibleVertices()
{
    LocalAssembly& graph = *this;

    // Find the vertices corresponding to vertexIdA and vertexIdB.
    vertex_descriptor vA = null_vertex();
    vertex_descriptor vB = null_vertex();
    BGL_FORALL_VERTICES(v, graph, LocalAssembly) {
        const LocalAssemblyVertex& vertex = graph[v];
        if(vertex.disjointSetId == disjointSetIdA) {
            SHASTA_ASSERT(vA == null_vertex());
            vA = v;
        }
        if(vertex.disjointSetId == disjointSetIdB) {
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

            BGL_FORALL_OUTEDGES(v0, e, graph, LocalAssembly) {
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

            BGL_FORALL_INEDGES(v0, e, graph, LocalAssembly) {
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
    BGL_FORALL_VERTICES(v, graph, LocalAssembly) {
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



// Remove all vertices and edges and clear the vertexMap and assemblyPath.
// All other data are left alone.
void LocalAssembly::clear()
{
    LocalAssemblyBaseClass::clear();
    vertexMap.clear();
    assemblyPath.clear();
}



void LocalAssembly::writeOrientedReadsSequences() const
{
    if(not html) {
        return;
    }
    if(not options.showOrientedReads) {
        return;
    }

    ofstream fasta("LocalAssembly-OrientedReadSequences.fasta");

    for(const OrientedReadInfo& info: orientedReadInfos) {

        SHASTA_ASSERT(not info.markerInfos.empty());
        const uint64_t position0 = uint64_t(info.markerInfos.front().position) + kHalf;
        const uint64_t position1 = uint64_t(info.markerInfos.back().position) + kHalf;

        fasta <<
            ">" << info.orientedReadId << " " <<
            position0 << ":" << position1 <<
            " length " << position1-position0 << "\n";
        for(uint64_t position=position0; position!=position1; position++) {
            const Base base = reads.getOrientedReadBase(info.orientedReadId, uint32_t(position));
            fasta << base;
        }
        fasta << "\n";
    }
}

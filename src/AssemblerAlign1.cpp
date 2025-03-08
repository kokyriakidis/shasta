// Alternative alignment functions with 1 suffix (SeqAn).
#include "Assembler.hpp"
using namespace shasta;

// Standard library.
#include "chrono.hpp"

// Seqan.
#include <seqan/align.h>
namespace seqan = seqan2;

void Assembler::alignOrientedReads1(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    int matchScore,
    int mismatchScore,
    int gapScore)
{
    alignOrientedReads1(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        matchScore, mismatchScore, gapScore);
}



void Assembler::alignOrientedReads1(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore)
{
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    alignOrientedReads1(
        orientedReadId0,
        orientedReadId1,
        matchScore,
        mismatchScore,
        gapScore,
        alignment,
        alignmentInfo);
}



void Assembler::alignOrientedReads1(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo)
{

    // Use SeqAn to compute an alignment free at both ends.
    // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
    using namespace seqan;

    // Hide shasta::Alignment.
    using seqan::Alignment;

    // An oriented read is represented as a sequence of KmerId
    // (the KmerId's of its markers). We want to align a pair of
    // such sequences.
    using TSequence = String<KmerId>;

    // Other SeqAn types we need.
    using TStringSet = StringSet<TSequence>;
    using TDepStringSet = StringSet<TSequence, Dependent<> >;
    using TAlignGraph = Graph<Alignment<TDepStringSet> >;

#if 0
    // Access the markers of our oriented reads.
    const span<CompressedMarker> markers0 =
        markers[orientedReadId0.getValue()];
    const span<CompressedMarker> markers1 =
        markers[orientedReadId1.getValue()];
#endif

    // Get the marker KmerIds for the two oriented reads.
    array<span<KmerId>, 2> allMarkerKmerIds;
    array<vector<KmerId>, 2> allMarkerKmerIdsVectors;
    if(markerKmerIds.isOpen()) {
        allMarkerKmerIds[0] = markerKmerIds[orientedReadId0.getValue()];
        allMarkerKmerIds[1] = markerKmerIds[orientedReadId1.getValue()];
    } else {
        // This is slower and will happen if markerKmerIds is not available.
        // Resize the vectors and make the spans point to the vectors.
        // Then call getOrientedReadMarkerKmerIds to fill them in.
        allMarkerKmerIdsVectors[0].resize(markers.size(orientedReadId0.getValue()));
        allMarkerKmerIdsVectors[1].resize(markers.size(orientedReadId1.getValue()));
        allMarkerKmerIds[0] = span<KmerId>(allMarkerKmerIdsVectors[0]);
        allMarkerKmerIds[1] = span<KmerId>(allMarkerKmerIdsVectors[1]);
        getOrientedReadMarkerKmerIds(orientedReadId0, allMarkerKmerIds[0]);
        getOrientedReadMarkerKmerIds(orientedReadId1, allMarkerKmerIds[1]);
    }

    // Construct the sequences of KmerId's we want to align.
    // SeqAn uses 45 to represent gaps, so we add 100 to the KmerIds passed to SeqAn.
    TSequence seq0;
    for(const KmerId& kmerId: allMarkerKmerIds[0]) {
        appendValue(seq0, kmerId + 100);
    }
    TSequence seq1;
    for(const KmerId& kmerId: allMarkerKmerIds[1]) {
        appendValue(seq1, kmerId + 100);
    }

    // Store them in a SeqAn string set.
    TStringSet sequences;
    appendValue(sequences, seq0);
    appendValue(sequences, seq1);

    // Compute the alignment.
    TAlignGraph graph(sequences);
    // const auto t0 = std::chrono::steady_clock::now();
    /* const int score = */ globalAlignment(
        graph,
        Score<int, Simple>(matchScore, mismatchScore, gapScore),
        AlignConfig<true, true, true, true>(),
        LinearGaps());
    // const auto t1 = std::chrono::steady_clock::now();
    /*
    cout << "Number of markers in these oriented reads: " <<
        markers0.size() << " " << markers1.size() << endl;
    cout << "Alignment score is " << score << endl;
    cout << "Alignment computation took " << seconds(t1-t0) << " s." << endl;
    */

    // Extract the alignment from the graph.
    // This creates a single sequence consisting of the two rows
    // of the alignment, concatenated.
    TSequence align;
    convertAlignment(graph, align);
    const int totalAlignmentLength = int(seqan::length(align));
    SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
    const int alignmentLength = totalAlignmentLength / 2;
    // cout << "Alignment length " << alignmentLength << endl;



    // Fill in the alignment.
    alignment.clear();
    uint32_t ordinal0 = 0;
    uint32_t ordinal1 = 0;
    const uint32_t seqanGapValue = 45;
    for(int i=0;
        i<alignmentLength and
        ordinal0<allMarkerKmerIds[0].size() and
        ordinal1<allMarkerKmerIds[1].size();
        i++) {
        if( align[i] != seqanGapValue and
            align[i + alignmentLength] != seqanGapValue and
            allMarkerKmerIds[0][ordinal0] == allMarkerKmerIds[1][ordinal1]) {
            alignment.ordinals.push_back(array<uint32_t, 2>{ordinal0, ordinal1});
        }
        if(align[i] != seqanGapValue) {
            ++ordinal0;
        }
        if(align[i + alignmentLength] != seqanGapValue) {
            ++ordinal1;
        }
    }

    // Store the alignment info.
    alignmentInfo.create(alignment, uint32_t(allMarkerKmerIds[0].size()), uint32_t(allMarkerKmerIds[1].size()));


    // Debugging.
#if 0
    {
        ofstream debugOut("AlignDebug.txt");

        // Extract the two rows of the alignment.
        array<vector<uint32_t>, 2> alignment;
        alignment[0].resize(alignmentLength);
        alignment[1].resize(alignmentLength);
        for(int i=0; i<alignmentLength; i++) {
            alignment[0][i] = align[i] - 100;
            alignment[1][i] = align[i + alignmentLength] - 100;
        }



        // Write out the alignment.
        for(int i=0; i<alignmentLength; i++) {
            debugOut << i << " ";
            if(alignment[0][i] == seqanGapValue) {
                debugOut << "-";
            } else {
                debugOut << alignment[0][i];
            }
            debugOut << " ";
            if(alignment[1][i] == seqanGapValue) {
                debugOut << "-";
            } else {
                debugOut << alignment[1][i];
            }
            if(
                alignment[0][i]!=seqanGapValue and
                alignment[1][i]!=seqanGapValue and
                alignment[0][i]==alignment[1][i]) {
                debugOut << "***";
            }
            debugOut << "\n";
        }
        debugOut << flush;
    }
#endif

}


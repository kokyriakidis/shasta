#pragma once

/******************************************************************

Class ProjectedSAlignment describes the "projection" of
an Alignment in marker space to base space.

******************************************************************/

#include "array.hpp"
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class ProjectedAlignment;
    class ProjectedAlignmentSegment;

    class Alignment;
    class Assembler;
    class Base;
    class CompressedMarker;
    class LongBaseSequenceView;
    class OrientedReadId;
}



// Each ProjectedAlignmentSegment describes a pair of consecutive
// aligned markers in the Alignment in marker space.
class shasta::ProjectedAlignmentSegment {
public:

    ProjectedAlignmentSegment(
        uint32_t kHalf,
        const array<uint32_t, 2>& ordinalA,
        const array<uint32_t, 2>& ordinalB,
        const array< span<const CompressedMarker>, 2>& markers);

    // In the arrays below, the index can be 0 or 1 and correspond
    // to the first and second oriented reads in the Alignment.
    // The A and B suffix refer to the left and right markers
    // in the pair of consecutive aligned markers.

    // The ordinals of the pair of consecutive aligned markers
    // described by this ProjectedAlignmentSegment.
    array<uint32_t, 2> ordinalsA;
    array<uint32_t, 2> ordinalsB;

    // The markers for the two oriented reads in this alignment.
    const array< span<const CompressedMarker>, 2>& markers;

    // The begin/end positions in base space, taken
    // at the midpoints of the two consecutive aligned markers.
    array<uint32_t, 2> positionsA;
    array<uint32_t, 2> positionsB;

    // The corresponding Base sequences.
    array<vector<Base>, 2> sequences;

    // The alignment between the two sequences.
    // See seqan.hpp for its meaning.
    int64_t editDistance;
    vector< pair<bool, bool> > alignment;
    void computeAlignment(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore);

    // The Base sequences in RLE represenation.
    array<vector<Base>, 2> rleSequences;
    void fillRleSequences();

    // The alignment between the two RLE sequences.
    // See seqan.hpp for its meaning.
    int64_t rleEditDistance;
    vector< pair<bool, bool> > rleAlignment;
    void computeRleAlignment(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore);


    void writeAlignmentHtml(ostream&) const;
    void writeRleAlignmentHtml(ostream&) const;

    void writeHtml(ostream&) const;
};



// Class ProjectedSAlignment describes the "projection" of
// an Alignment in marker space to base space.
// It is a sequence of ProjectedAlignmentSegments as defined above.
class shasta::ProjectedAlignment {
public:
    vector<ProjectedAlignmentSegment> segments;

    ProjectedAlignment(
        const Assembler&,
        const array<OrientedReadId, 2>&,
        const Alignment&);

    ProjectedAlignment(
        uint32_t k,
        const array<OrientedReadId, 2>&,
        const array<LongBaseSequenceView, 2>&,
        const Alignment&,
        const array< span<const CompressedMarker>, 2>& markers);

    // Marker length and its half.
    uint32_t k;
    uint32_t kHalf;

    // The two OrientedReadIds in this alignment.
    const array<OrientedReadId, 2>& orientedReadIds;

    // The base sequences of the reads.
    // These require reverse complementing for OrientedReadIds on strand 1.
    const array<LongBaseSequenceView, 2>& sequences;

    void fillSequences(ProjectedAlignmentSegment&) const;

    // Get the Base at a given position in one of the two oriented reads,
    // doing a reverse complement if necessary.
    Base getBase(uint64_t i, uint32_t position) const;

    // The markers for the two oriented reads in this alignment.
    const array< span<const CompressedMarker>, 2>& markers;



    // Statistics for the entire ProjectedAlignment.

    // Number of aligned bases (raw and RLE) in the aligned portions of the two oriented reads.
    array<uint64_t, 2> totalLength;
    array<uint64_t, 2> totalLengthRle;

    // Total edit distance (raw and RLE).
    int64_t totalEditDistance;
    int64_t totalEditDistanceRle;

    void computeStatistics();
    double errorRate() const;
    double errorRateRle() const;
    double Q() const;
    double QRle() const;

    void writeStatisticsHtml(ostream&) const;
    void writeHtml(ostream&, bool brief) const;
};


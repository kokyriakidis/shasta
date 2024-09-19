#include "ProjectedAlignment.hpp"
#include "Alignment.hpp"
#include "Base.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "seqan.hpp"
using namespace shasta;



ProjectedAlignment::ProjectedAlignment(
    uint32_t k,
    const array<OrientedReadId, 2>& orientedReadIds,
    const array<LongBaseSequenceView, 2>& sequences,
    const Alignment& alignment,
    const array< span<const CompressedMarker>, 2>& markers) :
    k(k),
    kHalf(k / 2),
    orientedReadIds(orientedReadIds),
    sequences(sequences),
    markers(markers)
{
    SHASTA_ASSERT((k % 2) == 0);

    // Scoring scheme for edit distance.
    const int64_t matchScore = 0;
    const int64_t mismatchScore = -1;
    const int64_t gapScore = -1;

    // Loop over pairs of consecutive aligned markers (A, B).
    for(uint64_t iB=1; iB<alignment.ordinals.size(); iB++) {
        const uint64_t iA = iB - 1;

        // Get the ordinals of these pair of consecutive aligned markers.
        const array<uint32_t, 2>& ordinalsA = alignment.ordinals[iA];
        const array<uint32_t, 2>& ordinalsB = alignment.ordinals[iB];

        segments.push_back(ProjectedAlignmentSegment(
            kHalf,
            ordinalsA,
            ordinalsB,
            markers));
        ProjectedAlignmentSegment& segment = segments.back();

        // Fill in the base sequences.
        fillSequences(segment);

        // Align them.
        segment.computeAlignment(matchScore, mismatchScore, gapScore);

        // Same, in RLE.
        segment.fillRleSequences();
        segment.computeRleAlignment(matchScore, mismatchScore, gapScore);
    }

    computeStatistics();
}



ProjectedAlignmentSegment::ProjectedAlignmentSegment(
    uint32_t kHalf,
    const array<uint32_t, 2>& ordinalsA,
    const array<uint32_t, 2>& ordinalsB,
    const array< span<const CompressedMarker>, 2>& markers) :
    ordinalsA(ordinalsA),
    ordinalsB(ordinalsB),
    markers(markers)
{
    for(uint64_t i=0; i<2; i++) {
        positionsA[i] = markers[i][ordinalsA[i]].position + kHalf;
        positionsB[i] = markers[i][ordinalsB[i]].position + kHalf;
    }
}



void ProjectedAlignmentSegment::computeAlignment(
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore)
{
    const vector<uint8_t>& sequence0 = reinterpret_cast< const vector<uint8_t>& >(sequences[0]);
    const vector<uint8_t>& sequence1 = reinterpret_cast< const vector<uint8_t>& >(sequences[1]);

    editDistance =  -seqanAlign(
        sequence0.begin(), sequence0.end(),
        sequence1.begin(), sequence1.end(),
        matchScore,
        mismatchScore,
        gapScore,
        false,
        false,
        alignment);

}



void ProjectedAlignmentSegment::computeRleAlignment(
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore)
{
    const vector<uint8_t>& sequence0 = reinterpret_cast< const vector<uint8_t>& >(rleSequences[0]);
    const vector<uint8_t>& sequence1 = reinterpret_cast< const vector<uint8_t>& >(rleSequences[1]);

    rleEditDistance =  -seqanAlign(
        sequence0.begin(), sequence0.end(),
        sequence1.begin(), sequence1.end(),
        matchScore,
        mismatchScore,
        gapScore,
        false,
        false,
        rleAlignment);

}



void ProjectedAlignment::writeHtml(ostream& html, bool brief) const
{
    html <<
        "<table>"
        "<tr>"
        "<th colspan=6>" << orientedReadIds[0] <<
        "<th colspan=6>" << orientedReadIds[1] <<
        "<th rowspan=2>Edit<br>distance"
        "<th rowspan=2>RLE<br>edit<br>distance"
        "<th rowspan=2 class=left>Alignments"
        "<tr>"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Skip"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Skip"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length";


    for(const ProjectedAlignmentSegment& segment: segments) {
        if(brief and (segment.editDistance == 0)) {
            continue;
        }
        segment.writeHtml(html);
    }

    html << "</table>";
}



void ProjectedAlignmentSegment::writeHtml(ostream& html) const
{
    html << "<tr>";

    for(uint64_t i=0; i<2; i++) {
        html << "<td class=centered>" << ordinalsA[i];
        html << "<td class=centered>" << ordinalsB[i];
        html << "<td class=centered>" << ordinalsB[i] - ordinalsA[i] - 1;
        html << "<td class=centered>" << positionsA[i];
        html << "<td class=centered>" << positionsB[i];
        html << "<td class=centered>" << positionsB[i] - positionsA[i];
    }

    html << "<td class=centered>" << editDistance;
    html << "<td class=centered>" << rleEditDistance;

    html << "<td class=left style='font-family:courier'>";
    writeAlignmentHtml(html);
    html << "<br><br>";
    writeRleAlignmentHtml(html);
}



void ProjectedAlignmentSegment::writeAlignmentHtml(ostream& html) const
{
    const vector<Base>& sequence0 = sequences[0];
    const vector<Base>& sequence1 = sequences[1];

    uint64_t position0 = 0;
    uint64_t position1 = 0;
    std::ostringstream alignment0;
    std::ostringstream alignment1;

    for(const pair<bool, bool>& p: alignment) {
        const bool hasBase0 = p.first;
        const bool hasBase1 = p.second;

        if(hasBase0) {
            alignment0 << sequence0[position0++];
        } else {
            alignment0 << "-";
        }

        if(hasBase1) {
            alignment1 << sequence1[position1++];
        } else {
            alignment1 << "-";
        }

    }

    SHASTA_ASSERT(position0 == sequence0.size());
    SHASTA_ASSERT(position1 == sequence1.size());

    const string alignment0String = alignment0.str();
    const string alignment1String = alignment1.str();

    for(uint64_t i=0; i<alignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style='background-color:pink'>";
        }
        html << alignment0String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

    html << "<br>";

    for(uint64_t i=0; i<alignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style=background-color:pink'>";
        }
        html << alignment1String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

}



void ProjectedAlignmentSegment::writeRleAlignmentHtml(ostream& html) const
{
    const vector<Base>& sequence0 = rleSequences[0];
    const vector<Base>& sequence1 = rleSequences[1];

    uint64_t position0 = 0;
    uint64_t position1 = 0;
    std::ostringstream alignment0;
    std::ostringstream alignment1;

    for(const pair<bool, bool>& p: rleAlignment) {
        const bool hasBase0 = p.first;
        const bool hasBase1 = p.second;

        if(hasBase0) {
            alignment0 << sequence0[position0++];
        } else {
            alignment0 << "-";
        }

        if(hasBase1) {
            alignment1 << sequence1[position1++];
        } else {
            alignment1 << "-";
        }

    }

    SHASTA_ASSERT(position0 == sequence0.size());
    SHASTA_ASSERT(position1 == sequence1.size());

    const string alignment0String = alignment0.str();
    const string alignment1String = alignment1.str();

    for(uint64_t i=0; i<rleAlignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style='background-color:pink'>";
        }
        html << alignment0String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

    html << "<br>";

    for(uint64_t i=0; i<rleAlignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style=background-color:pink'>";
        }
        html << alignment1String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

}



void ProjectedAlignment::fillSequences(ProjectedAlignmentSegment& segment) const
{
    for(uint64_t i=0; i<2; i++) {
        vector<Base>& sequence = segment.sequences[i];
        sequence.clear();
        for(uint32_t position=segment.positionsA[i]; position!=segment.positionsB[i]; position++) {
            sequence.push_back(getBase(i, position));
        }
        SHASTA_ASSERT(not sequence.empty());
    }
}




Base ProjectedAlignment::getBase(uint64_t i, uint32_t position) const
{
    const LongBaseSequenceView& sequence = sequences[i];

    if(orientedReadIds[i].getStrand() == 0) {
        return sequence[position];
    } else {
        return sequence[sequence.baseCount - 1 - position].complement();
    }
}



void ProjectedAlignmentSegment::fillRleSequences()
{
    for(uint64_t i=0; i<2; i++) {
        const vector<Base>& sequence = sequences[i];
        vector<Base>& rleSequence = rleSequences[i];
        rleSequence.clear();

        for(const Base b: sequence) {
            if(rleSequence.empty()) {
                rleSequence.push_back(b);
            } else {
                if(rleSequence.back() != b) {
                    rleSequence.push_back(b);
                }
            }
        }
    }
}



void ProjectedAlignment::computeStatistics()
{
    // Compute total lengths.
    totalLength = {0, 0};
    totalLengthRle = {0, 0};
    for(const ProjectedAlignmentSegment& segment: segments) {
        for(uint64_t i=0; i<2; i++) {
            totalLength[i] += segment.sequences[i].size();
            totalLengthRle[i] += segment.rleSequences[i].size();
        }
    }

    // Sanity check on the total lengths.
    for(uint64_t i=0; i<2; i++) {
        SHASTA_ASSERT(totalLength[i] == segments.back().positionsB[i] - segments.front().positionsA[i]);
    }

    // Compute total edit distances.
    totalEditDistance = 0;
    totalEditDistanceRle = 0;
    for(const ProjectedAlignmentSegment& segment: segments) {
        totalEditDistance += segment.editDistance;
        totalEditDistanceRle += segment.rleEditDistance;
    }
}



double ProjectedAlignment::errorRate() const
{
    return double(totalEditDistance) / double(totalLength[0] + totalLength[1]);
}



double ProjectedAlignment::errorRateRle() const
{
    return double(totalEditDistanceRle) /double(totalLengthRle[0] + totalLengthRle[1]);
}



double ProjectedAlignment::Q() const
{
    return -10. * log10(errorRate());
}



double ProjectedAlignment::QRle() const
{
    return -10. * log10(errorRateRle());
}



void ProjectedAlignment::writeStatisticsHtml(ostream& html) const
{
    using std::fixed;
    using std::setprecision;

    html << "<table>";

    // Header line.
    html <<
        "<tr>"
        "<th>"
        "<th>" << orientedReadIds[0] <<
        "<th>" << orientedReadIds[1] <<
        "<th>Total";

    // Length.
    html <<
        "<tr>"
        "<th class=left>Length of aligned portion" <<
        "<td class=centered>" << totalLength[0] <<
        "<td class=centered>" << totalLength[1] <<
        "<td class=centered>" << totalLength[0] + totalLength[1];

    // RLE length.
    html <<
        "<tr>"
        "<th class=left>RLE Length of aligned portion" <<
        "<td class=centered>" << totalLengthRle[0] <<
        "<td class=centered>" << totalLengthRle[1] <<
        "<td class=centered>" << totalLengthRle[0] + totalLengthRle[1];

    // Edit distance.
    html <<
        "<tr>"
        "<th class=left>Edit distance" <<
        "<td colspan=3 class=centered>" << totalEditDistance;

    // RLE edit distance.
    html <<
        "<tr>"
        "<th class=left>RLE edit distance" <<
        "<td colspan=3 class=centered>" << totalEditDistanceRle;

    // Error rate.
    html <<
        "<tr>"
        "<th class=left>Error rate" <<
        "<td colspan=3 class=centered>" << errorRate();

    // RLE error rate.
    html <<
        "<tr>"
        "<th class=left>RLE error rate" <<
        "<td colspan=3 class=centered>" << errorRateRle();

    // Q.
    html <<
        "<tr>"
        "<th class=left>Q (dB)" <<
        "<td colspan=3 class=centered>" << fixed << setprecision(1) << Q();

    // Q.
    html <<
        "<tr>"
        "<th class=left>RLE Q (dB)" <<
        "<td colspan=3 class=centered>" << fixed << setprecision(1) << QRle();

    html << "</table>";
}

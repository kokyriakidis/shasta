#include "Assembler.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;



void Assembler::createReadGraph4(
    uint32_t maxAlignmentCount)
{
    const bool debug = false;

    // QRle threshold to use an alignment in the read graph.
    const double minQRle = 35.;

    const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

    cout << timestamp << "createReadGraph4 begins, maxAlignmentCount " << maxAlignmentCount << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);

    // The computation of projected alignments is expensive, and so it should be done once for each alignment 
    // in an initial step and store the Error Rates of the projected alignment
    vector<double> alignmentErrorRateRle(alignmentCount);

    // This will hold the decomepressed Alignment.
    // Defined here to reduce memory allocation activity.
    Alignment alignment;

    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if((alignmentId % 1000) == 0) {
            cout << timestamp << alignmentId << "/" << alignmentCount << endl;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];
        thisAlignmentData.info.isInReadGraph = 0;

        // The alignment is stored as an alignment between readId0 on strand 0
        // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
        // The reverse complement alignment also exists, but is not stored explicitly.
        const ReadId readId0 = thisAlignmentData.readIds[0];
        const ReadId readId1 = thisAlignmentData.readIds[1];
        const bool isSameStrand = thisAlignmentData.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
        const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

        // The alignment is stored in compressed form as a string,
        // so we have to decompress it.
        span<const char> compressedAlignment = compressedAlignments[alignmentId];
        shasta::decompress(compressedAlignment, alignment);

        // Project this alignment to base space.
        const ProjectedAlignment projectedAlignment(
            *this,
            {orientedReadId0, orientedReadId1},
            alignment,
            true);

        const double errorRateRle = projectedAlignment.errorRateRle();

        alignmentErrorRateRle[alignmentId] = errorRateRle;

    }


    // Vector to keep the alignments for each read,
    // with their number of markers.
    // Contains pairs(errorRateRle, alignment id).
    vector< pair<double, uint32_t> > readAlignments;

    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Loop over reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        if(debug) {
            cout << "Working on read " << readId << endl;
        }

        // Gather the alignments for this read, each with its number of markers.
        readAlignments.clear();

        for(const uint32_t alignmentId: alignmentTable[OrientedReadId(readId, 0).getValue()]) {

            const double errorRateRle = alignmentErrorRateRle[alignmentId];

            if(errorRateRle<= maxErrorRateRle) {
                readAlignments.push_back(make_pair(errorRateRle, alignmentId));
            }

        }

        if(debug) {
            cout << "Found " << readAlignments.size() << " alignments." << endl;
        }

        // Keep the best maxAlignmentCount.
        if(readAlignments.size() > maxAlignmentCount) {
            std::nth_element(
                readAlignments.begin(),
                readAlignments.begin() + maxAlignmentCount,
                readAlignments.end(),
                std::less< pair<double, uint32_t> >());
            readAlignments.resize(maxAlignmentCount);
        }
        if(debug) {
            cout << "Kept " << readAlignments.size() << " alignments." << endl;
        }

        // Mark the surviving alignments as to be kept.
        for(const auto& p: readAlignments) {
            // Get information for this alignment.
            const uint32_t alignmentId = p.second;
            AlignmentData& thisAlignmentData = alignmentData[alignmentId];
            keepAlignment[alignmentId] = true;
            thisAlignmentData.info.isInReadGraph = 1;
            if(debug) {
                const AlignmentData& alignment = alignmentData[alignmentId];
                cout << "Marked alignment " << alignment.readIds[0] << " " <<
                    alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
            }
        }
    }


    cout << timestamp << "Done processing alignments." << endl;

    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);

    cout << timestamp << "createReadGraph4 ends." << endl;
}

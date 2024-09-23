#include "Assembler.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;



void Assembler::createReadGraph4()
{
    // QRle threshold to use an alignment in the read graph.
    const double minQRle = 25.;

    cout << timestamp << "createReadGraph4 begins." << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);

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
        // We call this with quick=true, so it stores only the following:
        // - RLE sequences and RLE alignments for segments for which the RLE sequences
        //   of the two oriented reads are different.
        // - Total RLE edit distance and total RLE lengths.
        const ProjectedAlignment projectedAlignment(
            *this,
            {orientedReadId0, orientedReadId1},
            alignment,
            true);

        // If the RLE Q is large enough, flag this alignment as to be kept.
        if(projectedAlignment.QRle() >= minQRle) {
            keepAlignment[alignmentId] = true;
            thisAlignmentData.info.isInReadGraph = 1;
        } else {
            thisAlignmentData.info.isInReadGraph = 0;
        }
    }
    cout << timestamp << "Done processing alignments." << endl;

    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);

    cout << timestamp << "createReadGraph4 ends." << endl;
}

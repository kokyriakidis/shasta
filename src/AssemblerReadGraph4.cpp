#include "Assembler.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;

//class AlignmentStats{public: double errorRateRle; uint32_t alignedRange; uint32_t rightUnaligned; uint32_t leftUnaligned; uint32_t alignmentId;};



void Assembler::createReadGraph4(
    uint32_t maxAlignmentCount)
{
    const bool debug = false;

    // QRle threshold to use an alignment in the read graph.
    const double minQRle = 42.;

    const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

    cout << timestamp << "createReadGraph4 begins, maxAlignmentCount " << maxAlignmentCount << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);

    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if((alignmentId % 1000) == 0) {
            cout << timestamp << alignmentId << "/" << alignmentCount << endl;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];
        // keepAlignment[alignmentId] = true;
        // thisAlignmentData.info.isInReadGraph = 1;

        // The alignment is stored as an alignment between readId0 on strand 0
        // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
        // The reverse complement alignment also exists, but is not stored explicitly.
        const ReadId readId0 = thisAlignmentData.readIds[0];
        const ReadId readId1 = thisAlignmentData.readIds[1];
        const bool isSameStrand = thisAlignmentData.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
        const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

        const double errorRateRle = thisAlignmentData.info.errorRateRle;

        // If the RLE Q is large enough, flag thus alignment as to be kept.
        if((errorRateRle <= maxErrorRateRle)) {
            keepAlignment[alignmentId] = true;
            thisAlignmentData.info.isInReadGraph = 1;
        } else {
            keepAlignment[alignmentId] = false;
            thisAlignmentData.info.isInReadGraph = 0;
        }

    }

    cout << timestamp << "Done processing alignments." << endl;

    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);

    cout << timestamp << "createReadGraph4 ends." << endl;

    removeReadGraphBridges(10);

}

















































































// void Assembler::createReadGraph4(
//     uint32_t maxAlignmentCount)
// {
//     const bool debug = false;

//     // QRle threshold to use an alignment in the read graph.
//     const double minQRle = 100000.;

//     const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

//     cout << timestamp << "createReadGraph4 begins, maxAlignmentCount skata " << maxAlignmentCount << endl;

//     // Get the total number of stored alignments.
//     const uint64_t alignmentCount = alignmentData.size();
//     SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

//     // Flag all alignments as not to be kept.
//     vector<bool> keepAlignment(alignmentCount, false);

//     // The computation of projected alignments is expensive, and so it should be done once for each alignment 
//     // in an initial step and store the Error Rates of the projected alignment
//     vector<double> alignmentErrorRateRle(alignmentCount);

//     // This will hold the decomepressed Alignment.
//     // Defined here to reduce memory allocation activity.
//     Alignment alignment;

//     // Loop over all alignments.
//     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         if((alignmentId % 1000) == 0) {
//             cout << timestamp << alignmentId << "/" << alignmentCount << endl;
//         }

//         // Get information for this alignment.
//         AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         thisAlignmentData.info.isInReadGraph = 0;

//         // The alignment is stored as an alignment between readId0 on strand 0
//         // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         // The reverse complement alignment also exists, but is not stored explicitly.
//         const ReadId readId0 = thisAlignmentData.readIds[0];
//         const ReadId readId1 = thisAlignmentData.readIds[1];
//         const bool isSameStrand = thisAlignmentData.isSameStrand;
//         SHASTA_ASSERT(readId0 < readId1);
//         const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
//         const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

//         // The alignment is stored in compressed form as a string,
//         // so we have to decompress it.
//         span<const char> compressedAlignment = compressedAlignments[alignmentId];
//         shasta::decompress(compressedAlignment, alignment);

//         // Project this alignment to base space.
//         const ProjectedAlignment projectedAlignment(
//             *this,
//             {orientedReadId0, orientedReadId1},
//             alignment,
//             true);

//         const double errorRateRle = projectedAlignment.errorRateRle();

//         alignmentErrorRateRle[alignmentId] = errorRateRle;

//     }


//     // Vector to keep the alignments for each read,
//     // with their number of markers.
//     // Contains pairs(errorRateRle, alignment id).
//     vector<AlignmentStats> readAlignments;

    
//     // Find the number of reads and oriented reads.
//     const ReadId orientedReadCount = uint32_t(markers.size());
//     SHASTA_ASSERT((orientedReadCount % 2) == 0);
//     const ReadId readCount = orientedReadCount / 2;

//     // Loop over reads.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         if(debug) {
//             cout << "Working on read " << readId << endl;
//         }

//         OrientedReadId alignmentOrientedReadId0(readId, 0);

//         // Loop over the alignments that this oriented read is involved in, with the proper orientation.
//         const vector< pair<OrientedReadId, AlignmentInfo> > correctOrientedAlignments =
//             findOrientedAlignments(alignmentOrientedReadId0, false);



//         // Gather the alignments for this read, considering alignment range and right unaligned portion.
//         readAlignments.clear();

//         for(uint32_t i=0; i<correctOrientedAlignments.size(); i++) {

//             const uint32_t alignmentId = alignmentTable[alignmentOrientedReadId0.getValue()][i];

//             const double errorRateRle = alignmentErrorRateRle[alignmentId];

//             // Calculate alignment range
//             const uint32_t alignmentRange = correctOrientedAlignments[i].second.markerCount;

//             // Calculate right unaligned portion
//             const uint32_t rightUnaligned0 = correctOrientedAlignments[i].second.rightTrim(0);

//             // Calculate left unaligned portion
//             const uint32_t leftUnaligned0 = correctOrientedAlignments[i].second.leftTrim(0);

//             // Calculate right unaligned portion
//             const uint32_t rightUnaligned1 = correctOrientedAlignments[i].second.rightTrim(1);

//             // Calculate left unaligned portion
//             const uint32_t leftUnaligned1 = correctOrientedAlignments[i].second.leftTrim(1);

//             // if(rightUnaligned1 != 0 and leftUnaligned1 != 0 and errorRateRle<= maxErrorRateRle) {
//             //     readAlignments.push_back(AlignmentStats{errorRateRle, alignmentRange, rightUnaligned1, leftUnaligned1, alignmentId});
//             // }

//             if(errorRateRle<= maxErrorRateRle) {
//                 readAlignments.push_back(AlignmentStats{errorRateRle, alignmentRange, rightUnaligned1, leftUnaligned1, alignmentId});
//             }

//         }
        
//         cout << "Working on read " << readId << endl;


//         // // Find the alignment with the highest alignedRange
//         // auto bestAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return a.alignedRange < b.alignedRange;
//         //     });

//         // if (bestAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestAlignment = *bestAlignmentIt;

//         //     // Clear existing alignments and keep only the best one
//         //     readAlignments.clear();
//         //     readAlignments.push_back(bestAlignment);

//         //     if(debug) {
//         //         cout << "Kept alignment with the highest alignedRange: " << bestAlignment.alignedRange << endl;
//         //     }
//         // } else {
//         //     if(debug) {
//         //         cout << "No alignments found for this read." << endl;
//         //     }
//         // }


//         // // Find the 10 alignments with the lowest errorRateRle
//         // std::partial_sort(readAlignments.begin(), 
//         //                   readAlignments.begin() + std::min(10UL, readAlignments.size()), 
//         //                   readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return a.errorRateRle < b.errorRateRle;
//         //     });

//         // // Find the 10 alignments with the highest sum of rightUnaligned, leftUnaligned, and alignedRange
//         // std::partial_sort(readAlignments.begin(), 
//         //                   readAlignments.begin() + std::min(10UL, readAlignments.size()), 
//         //                   readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return (a.rightUnaligned + a.leftUnaligned + a.alignedRange) >
//         //                (b.rightUnaligned + b.leftUnaligned + b.alignedRange);
//         //     });

//         cout << "read " << readId << " has that many alignments skata:" << readAlignments.size() << endl;

//         // // Keep only the top 10 alignments (or fewer if there are less than 10)
//         // readAlignments.resize(std::min(10UL, readAlignments.size()));

//         // cout << "read " << readId << " has that many alignments:" << readAlignments.size() << endl;

//         // if(debug) {
//         //     cout << "Top 5 alignments (or fewer) based on sum of unaligned portions and aligned range:" << endl;
//         //     for(const auto& alignment : readAlignments) {
//         //         cout << "AlignmentId: " << alignment.alignmentId 
//         //              << ", Sum: " << (alignment.rightUnaligned + alignment.leftUnaligned + alignment.alignedRange)
//         //              << " (Right: " << alignment.rightUnaligned 
//         //              << ", Left: " << alignment.leftUnaligned 
//         //              << ", Aligned: " << alignment.alignedRange << ")" << endl;
//         //     }
//         // }


//         // // Find the alignment with the highest alignedRange which also has the highest rightUnaligned
//         // auto bestRightAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         // if (a.alignedRange == b.alignedRange) {
//         //         //     return a.rightUnaligned > b.rightUnaligned;
//         //         // }
//         //         // return a.alignedRange < b.alignedRange;
//         //         return a.rightUnaligned+a.leftUnaligned+a.alignedRange < b.rightUnaligned+b.leftUnaligned+b.alignedRange;
//         //     });

//         // if (bestRightAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestRightAlignment = *bestRightAlignmentIt;

//         //     readAlignmentsRightTrimmed.clear();
//         //     readAlignmentsRightTrimmed.push_back(bestRightAlignment);
//         // }


//         // // Find the alignment with the highest alignedRange which also has the highest leftUnaligned
//         // auto bestLeftAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         // if (a.alignedRange == b.alignedRange) {
//         //         //     return a.leftUnaligned > b.leftUnaligned;
//         //         // }
//         //         // return a.alignedRange < b.alignedRange;
//         //         return a.leftUnaligned < b.leftUnaligned;
//         //     });

//         // if (bestLeftAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestLeftAlignment = *bestLeftAlignmentIt;

//         //     readAlignmentsLeftTrimmed.clear();
//         //     readAlignmentsLeftTrimmed.push_back(bestLeftAlignment);
//         // }

        

//         // if (!readAlignmentsRightTrimmed.empty()){
//         // cout << "The best rightTrimm alignment is " << readAlignmentsRightTrimmed[0].alignmentId << " alignmentId." << endl;
//         // cout << "The best rightTrimm has " << readAlignmentsRightTrimmed[0].rightUnaligned << " rightUnaligned." << endl;
//         // }

//         // if (!readAlignmentsLeftTrimmed.empty()){
//         //     cout << "The best leftTrimm alignment is " << readAlignmentsLeftTrimmed[0].alignmentId << " alignmentId." << endl;
//         //     cout << "The best leftTrimm has " << readAlignmentsRightTrimmed[0].leftUnaligned << " leftUnaligned." << endl;
//         // }
        

//         // if(debug) {
//         //     cout << "Found " << readAlignments.size() << " alignments." << endl;
//         // }

//         // // Keep the best maxAlignmentCount.
//         // if(readAlignments.size() > maxAlignmentCount) {
//         //     std::nth_element(
//         //         readAlignments.begin(),
//         //         readAlignments.begin() + maxAlignmentCount,
//         //         readAlignments.end(),
//         //         std::less< pair<double, uint32_t> >());
//         //     readAlignments.resize(maxAlignmentCount);
//         // }
//         // if(debug) {
//         //     cout << "Kept " << readAlignments.size() << " alignments." << endl;
//         // }

//         // Mark the surviving alignments as to be kept.
//         for(const auto& p: readAlignments) {
//             // Get information for this alignment.
//             const uint32_t alignmentId = p.alignmentId;
//             AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//             keepAlignment[alignmentId] = true;
//             thisAlignmentData.info.isInReadGraph = 1;
//             if(debug) {
//                 const AlignmentData& alignment = alignmentData[alignmentId];
//                 cout << "Marked alignment " << alignment.readIds[0] << " " <<
//                     alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
//             }
//         }

//         // // Mark the surviving alignments as to be kept.
//         // for(const auto& p: readAlignmentsLeftTrimmed) {
//         //     // Get information for this alignment.
//         //     const uint32_t alignmentId = p.alignmentId;
//         //     AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         //     keepAlignment[alignmentId] = true;
//         //     thisAlignmentData.info.isInReadGraph = 1;
//         //     if(debug) {
//         //         const AlignmentData& alignment = alignmentData[alignmentId];
//         //         cout << "Marked alignment " << alignment.readIds[0] << " " <<
//         //             alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
//         //     }
//         // }
//     }


//     cout << timestamp << "Done processing alignments." << endl;

//     const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

//     // Create the read graph using the alignments we selected.
//     createReadGraphUsingSelectedAlignments(keepAlignment);

//     cout << timestamp << "createReadGraph4 ends." << endl;
// }

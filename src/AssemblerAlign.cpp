// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"

// Shasta.
#include "Assembler.hpp"
#include "Alignment.hpp"
#include "AlignmentGraph.hpp"
#include "Align4.hpp"
#include "Align6.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "performanceLog.hpp"
#include "ProjectedAlignment.hpp"
#include "Reads.hpp"
#include "span.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard libraries.
#include "chrono.hpp"
#include <cmath>
#include <filesystem>
#include "iterator.hpp"
#include "tuple.hpp"



// Compute a marker alignment of two oriented reads.
void Assembler::alignOrientedReads(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    size_t maxSkip,     // Maximum ordinal skip allowed.
    size_t maxDrift,    // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    alignOrientedReads(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        maxSkip, maxDrift, maxMarkerFrequency
        );
}
void Assembler::alignOrientedReads(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxDrift,    // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    reads->checkReadsAreOpen();
    reads->checkReadNamesAreOpen();
    checkMarkersAreOpen();


    // Get the markers sorted by kmerId.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);
    getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

    // Call the lower level function.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

    // Compute the AlignmentInfo.
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
    cout << orientedReadId0 << " has " << reads->getRead(orientedReadId0.getReadId()).baseCount;
    cout << " bases and " << markersSortedByKmerId[0].size() << " markers." << endl;
    cout << orientedReadId1 << " has " << reads->getRead(orientedReadId1.getReadId()).baseCount;
    cout << " bases and " << markersSortedByKmerId[1].size() << " markers." << endl;
    cout << "The alignment has " << alignmentInfo.markerCount;
    cout << " markers. Left trim " << leftTrim;
    cout << " markers, right trim " << rightTrim << " markers." << endl;

#if 0
    // For convenience, also write the two oriented reads.
    ofstream fasta("AlignedOrientedReads.fasta");
    writeOrientedRead(orientedReadId0, fasta);
    writeOrientedRead(orientedReadId1, fasta);
#endif
}



// This lower level version takes as input vectors of
// markers already sorted by kmerId.
void Assembler::alignOrientedReads(
    const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxDrift,            // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    // Compute the alignment.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
}



void Assembler::alignOrientedReads(
    const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxDrift,            // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    bool debug,
    AlignmentGraph& graph,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo
)
{
    align(markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
}



// Compute marker alignments of an oriented read with all reads
// for which we have an alignment.
void Assembler::alignOverlappingOrientedReads(
    ReadId readId, Strand strand,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxDrift,                // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    alignOverlappingOrientedReads(
        OrientedReadId(readId, strand),
        maxSkip, maxDrift, maxMarkerFrequency, minAlignedMarkerCount, maxTrim);
}



void Assembler::alignOverlappingOrientedReads(
    OrientedReadId orientedReadId0,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxDrift,                // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    // Check that we have what we need.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Get the markers for orientedReadId0.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);

    // Loop over all alignments involving this oriented read.
    size_t goodAlignmentCount = 0;
    for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
        const AlignmentData& ad = alignmentData[i];

        // Get the other oriented read involved in this alignment.
        const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);

        // Get the markers for orientedReadId1.
        getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

        // Compute the alignment.
        AlignmentGraph graph;
        Alignment alignment;
        AlignmentInfo alignmentInfo;
        const bool debug = false;
        alignOrientedReads(
            markersSortedByKmerId,
            maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();

        cout << orientedReadId0 << " " << orientedReadId1 << " " << alignmentInfo.markerCount;
        if(alignmentInfo.markerCount) {
            cout << " " << leftTrim << " " << rightTrim;
            if(alignmentInfo.markerCount >= minAlignedMarkerCount && leftTrim<=maxTrim && rightTrim<=maxTrim) {
                cout << " good";
                ++goodAlignmentCount;
            }
        }
        cout << endl;

    }
    cout << "Found " << goodAlignmentCount << " alignments out of ";
    cout << alignmentTable[orientedReadId0.getValue()].size() << "." << endl;

}



// Compute an alignment for each alignment candidate.
// Store the alignments the satisfy our criteria.
void Assembler::computeAlignments(

    const AlignOptions& alignOptions,
    bool computeProjectedAlignmentMetrics,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    size_t threadCount
)
{

    const auto tBegin = steady_clock::now();
    performanceLog << timestamp << "Begin computing alignments for ";
    performanceLog << alignmentCandidates.candidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    SHASTA_ASSERT(kmerChecker);
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.alignOptions = &alignOptions;
    data.computeProjectedAlignmentMetrics = computeProjectedAlignmentMetrics;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // For alignment method 4, compute sorted markers.
    if(alignOptions.alignMethod == 4) {
        cout << timestamp << "Computing sorted markers." << endl;
        computeSortedMarkers(threadCount);
    }

    // For alignment method 5, compute low frequency markers.
    if(alignOptions.alignMethod == 5) {
        cout << timestamp << "Computing unique markers." << endl;
        computeLowFrequencyMarkers(1, threadCount);
    }

    // For alignment method 6, compute Align6Markers.
    if(alignOptions.alignMethod == 6) {
        performanceLog << timestamp << "Computing Align6Markers." << endl;
        computeAlign6Markers(threadCount);
        performanceLog << timestamp << "Done computing Align6Markers." << endl;
    }

    // Pick the batch size for computing alignments.
    size_t batchSize = 10;
    if(batchSize > alignmentCandidates.candidates.size()/threadCount) {
        batchSize = alignmentCandidates.candidates.size()/threadCount;
    }
    if(batchSize == 0) {
        batchSize = 1;
    }


    // Compute the alignments.
    data.threadAlignmentData.resize(threadCount);
    data.threadCompressedAlignments.resize(threadCount);
    
    performanceLog << timestamp << "Alignment computation begins." << endl;
    cout << timestamp << "Alignment computation begins." << endl;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::computeAlignmentsThreadFunction, threadCount);
    performanceLog << timestamp << "Alignment computation completed." << endl;
    cout << timestamp << "Alignment computation completed." << endl;

    // Store the alignments found by each thread.
    performanceLog << timestamp << "Storing the alignment found by each thread." << endl;
    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize);
    compressedAlignments.createNew(largeDataName("CompressedAlignments"), largeDataPageSize);
    
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        for(const AlignmentData& ad: threadAlignmentData) {
            alignmentData.push_back(ad);
        }

        const auto threadCompressedAlignments = data.threadCompressedAlignments[threadId];
        const auto size = threadCompressedAlignments->size();
        for(size_t i=0; i<size; i++) {
            compressedAlignments.appendVector(
                (*threadCompressedAlignments)[i].begin(),
                (*threadCompressedAlignments)[i].end()
            );
        }

        // Clean up thread storage.
        data.threadCompressedAlignments[threadId]->remove();
    }

    // Release unused allocated memory.
    alignmentData.unreserve();
    compressedAlignments.unreserve();

    // Cleanup.
    if(alignOptions.alignMethod == 4) {
        sortedMarkers.remove();
    }
    if(alignOptions.alignMethod == 5) {
        lowFrequencyMarkers.remove();
    }
    if(alignOptions.alignMethod == 6) {
        align6Markers.remove();
    }

    cout << "Found and stored " << alignmentData.size() << " good alignments." << endl;
    performanceLog << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();
    enforceAlignmentTransitivity(alignOptions);
    

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Computation of alignments ";
    performanceLog << "completed in " << tTotal << " s." << endl;

    performanceLog << timestamp;
}


// XXX
// --- START OF: Enforce transitivity of the alignments.
// If read r0 is aligned with r1 and r2, but r1 and r2 are not aligned,
// we attempt to compute an alignment between r1 and r2.
// If a good alignment is found, it is added to the alignment data.
// 1. Iterate through all reads and their neighbors in the alignmentTable.
// 2. For each pair of neighbors (readId1, readId2) of a read (readId0), check if an alignment already exists between them.
// 3. If no alignment exists, compute one using the same method as computeAlignmentsThreadFunction.
// 4. If the new alignment is good (i.e., meets the criteria), add it to a temporary vector of new AlignmentData.
// 5. After checking all reads, add the new alignments to the main alignmentData vector.
// 6. Finally, recompute the alignmentTable to include the newly added alignments.
void Assembler::enforceAlignmentTransitivity(const AlignOptions& alignOptions)
{
    const auto tBegin = steady_clock::now();
    performanceLog << timestamp << "Begin enforcing alignment transitivity." << endl;

    // --- Gather pairs that need alignment in parallel ---
    const ReadId readCount = reads->readCount();
    size_t threadCount = std::thread::hardware_concurrency();
    enforceTransitivityData.threadPairsToAlign.assign(threadCount, {});
    setupLoadBalancing(readCount, 1);
    runThreads(&Assembler::enforceAlignmentTransitivityGatherThreadFunction, threadCount);

    // --- Consolidate the pairs from all threads ---
    vector<pair<OrientedReadId, OrientedReadId>> pairsToAlign;
    for(const auto& threadPairs : enforceTransitivityData.threadPairsToAlign) {
        pairsToAlign.insert(pairsToAlign.end(), threadPairs.begin(), threadPairs.end());
    }
    enforceTransitivityData.threadPairsToAlign.clear();
    
    // Sort and unique the oriented pairs first to remove exact duplicates
    sort(pairsToAlign.begin(), pairsToAlign.end());
    pairsToAlign.erase(unique(pairsToAlign.begin(), pairsToAlign.end()), pairsToAlign.end());

    // Now, unique by the pair of unoriented ReadIds.
    // To do this efficiently without a set, we create a temporary vector
    // of canonical unoriented pairs, sort it, then unique it.
    vector<pair<pair<ReadId, ReadId>, pair<OrientedReadId, OrientedReadId>>> tempPairs;
    tempPairs.reserve(pairsToAlign.size());
    for(const auto& orientedPair : pairsToAlign) {
        ReadId r1 = orientedPair.first.getReadId();
        ReadId r2 = orientedPair.second.getReadId();
        if(r1 > r2) {
            std::swap(r1, r2);
        }
        tempPairs.push_back({{r1, r2}, orientedPair});
    }
    std::sort(tempPairs.begin(), tempPairs.end(),
        [](const auto& a, const auto& b) {
            return a.first < b.first;
    });
    auto last = std::unique(tempPairs.begin(), tempPairs.end(),
        [](const auto& a, const auto& b) {
        return a.first == b.first;
    });
    tempPairs.erase(last, tempPairs.end());

    // Now extract the final list of oriented pairs to align.
    vector<pair<OrientedReadId, OrientedReadId>> uniquePairsToAlign;
    uniquePairsToAlign.reserve(tempPairs.size());
    for(const auto& p : tempPairs) {
        uniquePairsToAlign.push_back(p.second);
    }
    
    performanceLog << timestamp << "Found " << uniquePairsToAlign.size() << " unique pairs to align." << endl;
    
    // // DEBUG: Print the first 20 pairs
    // for(size_t i=0; i<20; i++) {
    //     performanceLog << uniquePairsToAlign[i].first.getReadId() << " " << " Strand " << uniquePairsToAlign[i].first.getStrand() << " " << uniquePairsToAlign[i].second.getReadId() << " " << " Strand " << uniquePairsToAlign[i].second.getStrand() << endl;
    // }

    // --- Align pairs in parallel ---
    enforceTransitivityData.alignOptions = &alignOptions;
    enforceTransitivityData.pairsToAlign = std::move(uniquePairsToAlign);
    enforceTransitivityData.threadAlignmentData.resize(threadCount);
    enforceTransitivityData.threadCompressedAlignments.resize(threadCount);
    for(size_t threadId=0; threadId<threadCount; ++threadId) {
        enforceTransitivityData.threadCompressedAlignments[threadId] =
            make_shared<MemoryMapped::VectorOfVectors<char, uint64_t>>();
        enforceTransitivityData.threadCompressedAlignments[threadId]->createNew(
            largeDataName("tmp-enforceTransitivity-CompressedAlignments-" + to_string(threadId)),
            largeDataPageSize);
    }

    setupLoadBalancing(enforceTransitivityData.pairsToAlign.size(), 1);
    runThreads(&Assembler::enforceAlignmentTransitivityThreadFunction, threadCount);

    // --- Store new alignments ---
    vector<AlignmentData> newAlignmentData;
    for (size_t threadId = 0; threadId < threadCount; threadId++) {
        for (const auto& ad : enforceTransitivityData.threadAlignmentData[threadId]) {
            newAlignmentData.push_back(ad);
        }
        enforceTransitivityData.threadAlignmentData[threadId].clear();
    }

    if (!newAlignmentData.empty()) {
        performanceLog << timestamp << "Adding " << newAlignmentData.size() << " new alignments to enforce transitivity." << endl;

        for(const auto& ad : newAlignmentData) {
            alignmentData.push_back(ad);
        }
        for (size_t threadId = 0; threadId < threadCount; threadId++) {
            const auto& threadCompressedAlignments = *enforceTransitivityData.threadCompressedAlignments[threadId];
            for (size_t i = 0; i < threadCompressedAlignments.size(); i++) {
                compressedAlignments.appendVector(threadCompressedAlignments[i].begin(), threadCompressedAlignments[i].end());
            }
            enforceTransitivityData.threadCompressedAlignments[threadId]->remove();
        }

        alignmentData.unreserve();
        compressedAlignments.unreserve();

        performanceLog << timestamp << "Recomputing alignment table." << endl;
        computeAlignmentTable();
    }

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Enforcing alignment transitivity completed in " << tTotal << " s." << endl;
}

void Assembler::enforceAlignmentTransitivityGatherThreadFunction(size_t threadId)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId0 = begin; readId0 < end; ++readId0) {
            const Strand strand0 = 0;
            OrientedReadId orientedReadId0(readId0, strand0);

            const span<uint32_t> alignmentIds = alignmentTable[orientedReadId0.getValue()];
            vector<pair<OrientedReadId, uint32_t>> neighbors0;
            neighbors0.reserve(alignmentIds.size());
            for (const auto alignmentId : alignmentIds) {
                neighbors0.emplace_back(alignmentData[alignmentId].getOther(orientedReadId0), alignmentId);
            }

            for (size_t i = 0; i < neighbors0.size(); i++) {
                const OrientedReadId orientedReadId1 = neighbors0[i].first;
                if (orientedReadId1 < orientedReadId0) continue;

                for (size_t j = i + 1; j < neighbors0.size(); j++) {
                    const OrientedReadId orientedReadId2 = neighbors0[j].first;
                    bool areNeighbors = false;
                    const span<uint32_t> neighbors1AlignmentIds = alignmentTable[orientedReadId1.getValue()];
                    for (const auto alignmentId1 : neighbors1AlignmentIds) {
                        if (alignmentData[alignmentId1].getOther(orientedReadId1) == orientedReadId2) {
                            areNeighbors = true;
                            break;
                        }
                    }
                    if (!areNeighbors) {
                        // Canonicalize the pair before adding it to ensure uniqueness.
                        auto id1 = orientedReadId1;
                        auto id2 = orientedReadId2;
                        if(id2 < id1) {
                            std::swap(id1, id2);
                        }
                        enforceTransitivityData.threadPairsToAlign[threadId].push_back({id1, id2});
                    }
                }
            }
        }
    }
}

void Assembler::enforceAlignmentTransitivityThreadFunction(size_t threadId)
{
    const auto& alignOptions = *enforceTransitivityData.alignOptions;

    // Work areas used inside the loop
    Alignment alignment;
    AlignmentInfo alignmentInfo;

    // Get the batch of pairs to process
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for (uint64_t i = begin; i < end; ++i) {
            const auto& pair = enforceTransitivityData.pairsToAlign[i];
            OrientedReadId orientedReadId1 = pair.first;
            OrientedReadId orientedReadId2 = pair.second;
            
            try {
                ofstream nullStream;
                alignOrientedReads5(orientedReadId1, orientedReadId2,
                    alignOptions.matchScore, alignOptions.mismatchScore, alignOptions.gapScore,
                    alignOptions.align5DriftRateTolerance, alignOptions.align5MinBandExtend,
                    alignment, alignmentInfo,
                    nullStream);

            } catch (std::exception& e) {
                continue; // Skip if alignment fails
            }

            const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
            const uint32_t positionStart1 = orientedReadMarkers1[alignmentInfo.data[0].firstOrdinal].position;
            const uint32_t positionEnd1 = orientedReadMarkers1[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assemblerInfo->k);
            const uint32_t overlapLength1 = positionEnd1 - positionStart1;

            const auto orientedReadMarkers2 = markers[orientedReadId2.getValue()];
            const uint32_t positionStart2 = orientedReadMarkers2[alignmentInfo.data[1].firstOrdinal].position;
            const uint32_t positionEnd2 = orientedReadMarkers2[alignmentInfo.data[1].lastOrdinal].position + uint32_t(assemblerInfo->k);
            const uint32_t overlapLength2 = positionEnd2 - positionStart2;

            // XXX
            // --- START OF: Require at least 1000 bases of overlap. ---
            //
            const uint32_t minOverlap = 1000;
            if (overlapLength1 < minOverlap || overlapLength2 < minOverlap) {
                continue;
            }
            // XXX
            // --- END OF: Require at least 1000 bases of overlap. ---
            //

            const uint64_t readLength1 = reads->getRead(orientedReadId1.getReadId()).baseCount;
            const uint64_t readLength2 = reads->getRead(orientedReadId2.getReadId()).baseCount;

            const uint32_t leftTrimBases1 = positionStart1;
            const uint32_t rightTrimBases1 = readLength1 - positionEnd1;

            const uint32_t leftTrimBases2 = positionStart2;
            const uint32_t rightTrimBases2 = readLength2 - positionEnd2;

            // XXX
            // --- START OF: Calculate the "effective" 5end and 3end hanging lengths. ---
            // This is the minimum of the left trim bases and minimum of the right trim bases.
            // This basically calculates the overlaping portion of the two reads that is not covered by the alignment.
            uint32_t min5endHangingLength = 0;
            uint32_t min3endHangingLength = 0;

            if (leftTrimBases1 < leftTrimBases2) {
                min5endHangingLength = leftTrimBases1;
            } else {
                min5endHangingLength = leftTrimBases2;
            }

            if (rightTrimBases1 < rightTrimBases2) {
                min3endHangingLength = rightTrimBases1;
            } else {
                min3endHangingLength = rightTrimBases2;
            }

            // If the hanging lengths are too long, skip the alignment.
            const uint32_t maxHangingLength = 1000;

            if (min5endHangingLength > maxHangingLength || min3endHangingLength > maxHangingLength) {
                // cout << orientedReadId1 << " " << orientedReadId2 << " long hanging lengths." << endl;
                continue;
            }
            // XXX
            // --- END OF: Calculate the "effective" 5end and 3end hanging lengths. ---
            //


            // // XXX
            // // --- START OF: At least 80% of the overlap should be covered by the alignment. ---
            // //
            // const float minOverlapFraction = 0.80;

            // if (overlapLength1 < minOverlapFraction * (overlapLength1 + min5endHangingLength + min3endHangingLength) || 
            //     overlapLength2 < minOverlapFraction * (overlapLength2 + min5endHangingLength + min3endHangingLength)) 
            // {
            //     // cout << orientedReadId1 << " " << orientedReadId2 << " low overlap fraction." << endl;
            //     continue;
            // }
            // // XXX
            // // --- END OF: At least 80% of the overlap should be covered by the alignment. ---
            // //

            array<OrientedReadId, 2> orientedReadIds = {orientedReadId1, orientedReadId2};
            const ProjectedAlignment projectedAlignment(
                *this,
                orientedReadIds,
                alignment,
                ProjectedAlignment::Method::QuickRaw);


            // XXX
            // --- START OF: Total Overlap Error Rate filtering    
            //

            float alignmentErrorRate = projectedAlignment.errorRate();
            
            // Skip alignments with error rate greater than 0.07.
            if (alignmentErrorRate > 0.07) {
                continue;
            }

            alignmentInfo.errorRate = alignmentErrorRate;
            alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);

            // XXX
            // --- END OF: Total Overlap Error Rate filtering    
            //



            if(orientedReadId1.getReadId() > orientedReadId2.getReadId()){
                swap(orientedReadId1, orientedReadId2);
                alignmentInfo.swap();
            }
            if(orientedReadId1.getStrand() == 1){
                orientedReadId1.flipStrand();
                orientedReadId2.flipStrand();
                alignmentInfo.reverseComplement();
            }

            OrientedReadPair candidate;
            candidate.readIds[0] = orientedReadId1.getReadId();
            candidate.readIds[1] = orientedReadId2.getReadId();
            candidate.isSameStrand = (orientedReadId1.getStrand() == orientedReadId2.getStrand());
            
            enforceTransitivityData.threadAlignmentData[threadId].push_back(AlignmentData(candidate, alignmentInfo));
            string compressedAlignment;
            shasta::compress(alignment, compressedAlignment);
            enforceTransitivityData.threadCompressedAlignments[threadId]->appendVector(
                compressedAlignment.begin(), compressedAlignment.end());
        }
    }
}

// XXX
// --- END OF: Enforce transitivity of the alignments.
//



void Assembler::computeAlignmentsThreadFunction(size_t threadId)
{

    array<OrientedReadId, 2> orientedReadIds;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    string compressedAlignment;

    const bool debug = false;
    auto& data = computeAlignmentsData;
    const size_t alignmentMethod = data.alignOptions->alignMethod;
    const uint32_t maxMarkerFrequency = data.alignOptions->maxMarkerFrequency;
    const size_t maxSkip = data.alignOptions->maxSkip;
    const size_t maxDrift = data.alignOptions->maxDrift;
    const size_t minAlignedMarkerCount = data.alignOptions->minAlignedMarkerCount;
    const double minAlignedFraction = data.alignOptions->minAlignedFraction;
    const size_t maxTrim = data.alignOptions->maxTrim;
    const int matchScore = data.alignOptions->matchScore;
    const int mismatchScore = data.alignOptions->mismatchScore;
    const int gapScore = data.alignOptions->gapScore;
    const double downsamplingFactor = data.alignOptions->downsamplingFactor;
    const int bandExtend = data.alignOptions->bandExtend;
    const int maxBand = data.alignOptions->maxBand;
    const bool suppressContainments = data.alignOptions->suppressContainments;
    const double align5DriftRateTolerance = data.alignOptions->align5DriftRateTolerance;
    const uint64_t align5MinBandExtend = data.alignOptions->align5MinBandExtend;
    const bool computeProjectedAlignmentMetrics = data.computeProjectedAlignmentMetrics;


    // Align4-specific items.
    Align4::Options align4Options;
    MemoryMapped::ByteAllocator byteAllocator;
    if(alignmentMethod == 4) {
        align4Options.deltaX = data.alignOptions->align4DeltaX;
        align4Options.deltaY = data.alignOptions->align4DeltaY;
        align4Options.minEntryCountPerCell = data.alignOptions->align4MinEntryCountPerCell;
        align4Options.maxDistanceFromBoundary = data.alignOptions->align4MaxDistanceFromBoundary;
        align4Options.minAlignedMarkerCount = minAlignedMarkerCount;
        align4Options.minAlignedFraction = minAlignedFraction;
        align4Options.maxSkip = maxSkip;
        align4Options.maxDrift = maxDrift;
        align4Options.maxTrim = maxTrim;
        align4Options.maxBand = maxBand;
        align4Options.matchScore = matchScore;
        align4Options.mismatchScore = mismatchScore;
        align4Options.gapScore = gapScore;
        byteAllocator.createNew(
            largeDataName("tmp-ByteAllocator-" + to_string(threadId)),
            largeDataPageSize, 2ULL * 1024 * 1024 * 1024);
    }

    ofstream nullStream;
    Align6 align6(
        assemblerInfo->k,
        maxSkip,
        maxDrift,
        data.alignOptions->align6Options,
        assemblerInfo->kmerDistributionInfo,
        nullStream);
    if((threadId == 0) and (alignmentMethod == 6)) {
        std::lock_guard<std::mutex> lock(mutex);
        align6.writeGlobalFrequencyCriteria(cout);
    }

    vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
    
    shared_ptr< MemoryMapped::VectorOfVectors<char, uint64_t> > thisThreadCompressedAlignmentsPointer =
        make_shared< MemoryMapped::VectorOfVectors<char, uint64_t> >();
    data.threadCompressedAlignments[threadId] = thisThreadCompressedAlignmentsPointer;
    auto& thisThreadCompressedAlignments = *thisThreadCompressedAlignmentsPointer;
    thisThreadCompressedAlignments.createNew(
        largeDataName("tmp-ThreadGlobalCompressedAlignments-" + to_string(threadId)),
        largeDataPageSize);

#if 0
    // A vector to store the time taken to compute each alignment.
    vector< pair<uint64_t, double> > elapsedTime;
#endif

    const uint64_t messageFrequency = min(1000000UL, alignmentCandidates.candidates.size()/20);

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for(size_t i=begin; i!=end; i++) {
#if 0
            const auto steadyClock0 = std::chrono::steady_clock::now();
#endif

            const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
            SHASTA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

            // Get the oriented read ids, with the first one on strand 0.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

            if((i % messageFrequency) == 0){
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "Working on alignment " << i;
                performanceLog << " of " << alignmentCandidates.candidates.size();
                // performanceLog << ": " << orientedReadIds[0] << " " << orientedReadIds[1];
                performanceLog << endl;
            }


            // Compute the alignment.
            try {
                if(alignmentMethod == 0) {

                    // Get the markers for the two oriented reads in this candidate.
                    for(size_t j=0; j<2; j++) {
                        getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
                    }

                    // Compute the Alignment.
                    alignOrientedReads(
                        markersSortedByKmerId,
                        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

                } else if(alignmentMethod == 1) {
                    alignOrientedReads1(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        alignment, alignmentInfo);
                } else if(alignmentMethod == 3) {
                    alignOrientedReads3(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        downsamplingFactor, bandExtend, maxBand,
                        alignment, alignmentInfo);
                } else if(alignmentMethod == 4) {
                    alignOrientedReads4(orientedReadIds[0], orientedReadIds[1],
                        align4Options,
                        byteAllocator,
                        alignment, alignmentInfo,
                        false);
                    SHASTA_ASSERT(byteAllocator.isEmpty());
                } else if(alignmentMethod == 5) {
                    ofstream nullStream;
                    alignOrientedReads5(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        align5DriftRateTolerance, align5MinBandExtend,
                        alignment, alignmentInfo,
                        nullStream);
                } else if(alignmentMethod == 6) {
                    ofstream nullStream;
                    alignOrientedReads6(orientedReadIds[0], orientedReadIds[1],
                        alignment, alignmentInfo, align6);
                } else {
                    SHASTA_ASSERT(0);
                }
            } catch (std::exception& e) {
                std::lock_guard<std::mutex> lock(mutex);
                cout <<
                    "An error occurred while computing a marker alignment "
                    " of oriented reads " << orientedReadIds[0] << " and " << orientedReadIds[1] <<
                    ". This alignment candidate will be skipped. Error description is: " <<
                    e.what() << endl;
                continue;
            } catch(...) {
                std::lock_guard<std::mutex> lock(mutex);
                cout <<
                    "An error occurred while computing a marker alignment "
                    " of oriented reads " << orientedReadIds[0] << " and " << orientedReadIds[1] <<
                    ". This alignment candidate will be skipped. " << endl;
                continue;
            }

#if 0
            const auto steadyClock1 = std::chrono::steady_clock::now();
            const auto deltaT = 1.e-9 * double(std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0).count());
            elapsedTime.push_back({i, deltaT});
#endif

            // // If the alignment has too few markers, skip it.
            // if(alignment.ordinals.size() < minAlignedMarkerCount) {
            //     // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " too few markers." << endl;
            //     continue;
            // }

            // // If the aligned fraction is too small, skip it.
            // if(min(alignmentInfo.alignedFraction(0), alignmentInfo.alignedFraction(1)) < minAlignedFraction) {
            //     continue;
            // }

            // // If the alignment has too much trim, skip it.
            // uint32_t leftTrim;
            // uint32_t rightTrim;
            // tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
            // if(leftTrim>maxTrim || rightTrim>maxTrim) {
            //     // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " too much trim." << endl;
            //     continue;
            // }

            // // For alignment methods other than method 0, we also need to check for
            // // maxSip and maxDrift. Method 0 does that automatically.
            // if(alignmentMethod != 0) {
            //     if(alignment.maxSkip() > maxSkip) {
            //         continue;
            //     }
            //     if(alignment.maxDrift() > maxDrift) {
            //         continue;
            //     }
            // }

            // // Skip containing alignments, if so requested.
            // if(suppressContainments and alignmentInfo.isContaining(uint32_t(maxTrim))) {
            //     continue;
            // }


            const auto orientedReadMarkers0 = markers[orientedReadIds[0].getValue()];
            const uint32_t positionStart0 = orientedReadMarkers0[alignmentInfo.data[0].firstOrdinal].position;
            const uint32_t positionEnd0 = orientedReadMarkers0[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assemblerInfo->k);
            const uint32_t overlapLength0 = positionEnd0 - positionStart0;

            const auto orientedReadMarkers1 = markers[orientedReadIds[1].getValue()];
            const uint32_t positionStart1 = orientedReadMarkers1[alignmentInfo.data[1].firstOrdinal].position;
            const uint32_t positionEnd1 = orientedReadMarkers1[alignmentInfo.data[1].lastOrdinal].position + uint32_t(assemblerInfo->k);
            const uint32_t overlapLength1 = positionEnd1 - positionStart1;

            // XXX
            // --- START OF: Require at least 1000 bases of overlap. ---
            //
            const uint32_t minOverlap = 1000;
            if (overlapLength0 < minOverlap || overlapLength1 < minOverlap) {
                continue;
            }
            // XXX
            // --- END OF: Require at least 1000 bases of overlap. ---
            //


            const uint64_t readLength0 = reads->getRead(orientedReadIds[0].getReadId()).baseCount;
            const uint64_t readLength1 = reads->getRead(orientedReadIds[1].getReadId()).baseCount;

            const uint32_t leftTrimBases0 = positionStart0;
            const uint32_t rightTrimBases0 = readLength0 - positionEnd0;

            const uint32_t leftTrimBases1 = positionStart1;
            const uint32_t rightTrimBases1 = readLength1 - positionEnd1;

            
            // XXX
            // --- START OF: Calculate the "effective" 5end and 3end hanging lengths. ---
            // This is the minimum of the left trim bases and minimum of the right trim bases.
            // This basically calculates the overlaping portion of the two reads that is not covered by the alignment.
            uint32_t min5endHangingLength = 0;
            uint32_t min3endHangingLength = 0;

            if (leftTrimBases0 < leftTrimBases1) {
                min5endHangingLength = leftTrimBases0;
            } else {
                min5endHangingLength = leftTrimBases1;
            }

            if (rightTrimBases0 < rightTrimBases1) {
                min3endHangingLength = rightTrimBases0;
            } else {
                min3endHangingLength = rightTrimBases1;
            }

            // If the hanging lengths are too long, skip the alignment.
            const uint32_t maxHangingLength = 1000;

            if (min5endHangingLength > maxHangingLength || min3endHangingLength > maxHangingLength) {
                // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " long hanging lengths." << endl;
                continue;
            }
            // XXX
            // --- END OF: Calculate the "effective" 5end and 3end hanging lengths. ---
            //




            // XXX
            // --- START OF: At least 80% of the overlap should be covered by the alignment. ---
            //
            const float minOverlapFraction = 0.80;

            if (overlapLength0 < minOverlapFraction * (overlapLength0 + min5endHangingLength + min3endHangingLength) || 
                overlapLength1 < minOverlapFraction * (overlapLength1 + min5endHangingLength + min3endHangingLength)) 
            {
                // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " low overlap fraction." << endl;
                continue;
            }
            // XXX
            // --- END OF: At least 80% of the overlap should be covered by the alignment. ---
            //




            // Compute projected alignment metrics, if requested.
            // if(computeProjectedAlignmentMetrics) {
            //     const ProjectedAlignment projectedAlignment(
            //         *this,
            //         orientedReadIds,
            //         alignment,
            //         ProjectedAlignment::Method::QuickRle);
            //     alignmentInfo.errorRate = float(projectedAlignment.errorRateRle());
            //     alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            // }
            if(computeProjectedAlignmentMetrics) {
                const ProjectedAlignment projectedAlignment(
                    *this,
                    orientedReadIds,
                    alignment,
                    ProjectedAlignment::Method::QuickRaw);


                // XXX
                // --- START OF: Total Overlap Error Rate filtering    
                //

                float alignmentErrorRate = projectedAlignment.errorRate();
                
                // Skip alignments with error rate greater than 0.07.
                if (alignmentErrorRate > 0.07) {
                    continue;
                }

                // XXX
                // --- END OF: Total Overlap Error Rate filtering    
                //

                

                // XXX
                // --- START OF: Sliding-window error filtering (read-0)    
                //
                

                // Window length.
                const uint64_t windowLength = 375;

                // Determine the number of windows.
                const uint64_t span0 = positionEnd0 - positionStart0 + 1;   // real number of bases

                // How many full windows do we need? – ceiling division
                const size_t numberOfWindows = (span0 + windowLength - 1) / windowLength;

                // How many bases fall into the last window?
                const uint64_t remainderBases = span0 % windowLength;              // 0 … L-1

                vector<uint32_t> windowBases(numberOfWindows, windowLength);        // start with L everywhere

                if (remainderBases != 0) {
                    windowBases.back() = remainderBases;               // shrink the last one
                }   



                // XXX
                // --- START OF: Preparation for the chemical window filtering.
                //

                uint64_t chemicalWindowLength = 384;
                uint64_t dn = 2000;
                double dr = 0.1;
                bool runChemicalWindowFiltering = false;
                if (span0 > chemicalWindowLength) {
                    
                    if (dn > (span0 * dr)) {
                        dn = span0 * dr;
                    }

                    if ((dn > chemicalWindowLength) && (span0 > dn)) {
                        runChemicalWindowFiltering = true;
                    }

                }

                uint64_t errorsInFirstDnBases = 0;
                uint64_t errorsInLastDnBases = 0;

                // XXX
                // --- END OF: Preparation for the chemical window filtering.
                //


                vector<uint32_t> windowErrors(numberOfWindows, 0);
                
                // Loop over the RAW segments of the projectedAlignment and find windows with high error rate.
                for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
                    
                    // Get the RAW sequences
                    const vector<Base>& sequence0 = segment.sequences[0];
                    const vector<Base>& sequence1 = segment.sequences[1];

                    const uint32_t segmentPositionStart0 = segment.positionsA[0];
                    const uint32_t segmentPositionEnd0 = segment.positionsB[0];
                    

                    // Align them base by base to get the right sequence alignment representation.
                    uint64_t localPosition0 = 0;
                    uint64_t localPosition1 = 0;

                    // Retreive the RAW alignment sequences.
                    for (const pair<bool, bool>& p: segment.alignment) {
                        const bool hasBase0 = p.first;
                        const bool hasBase1 = p.second;

                        // Find the window that contains this position.
                        uint64_t windowIndex = (segmentPositionStart0 + localPosition0 - positionStart0) / windowLength;

                        if (hasBase0 && hasBase1) {
                            if (sequence0[localPosition0] != sequence1[localPosition1]) {
                                // Increment the error count for this window.
                                ++windowErrors[windowIndex];

                                // Check if the position is in the first or last dn bases.
                                if ( ((segmentPositionStart0 + localPosition0 - positionStart0) <= dn) ) {
                                    ++errorsInFirstDnBases;
                                } else if ( ((segmentPositionStart0 + localPosition0) >= (positionEnd0 - dn))) {
                                    ++errorsInLastDnBases;
                                }

                            }
                            // Advance the position.
                            ++localPosition0;
                            ++localPosition1;
                            continue;
                        }

                        if (hasBase0 && !hasBase1) {
                            // Increment the error count for this window.
                            ++windowErrors[windowIndex];

                            // Check if the position is in the first or last dn bases.
                            if ( ((segmentPositionStart0 + localPosition0 - positionStart0) <= dn) ) {
                                ++errorsInFirstDnBases;
                            } else if ( ((segmentPositionStart0 + localPosition0) >= (positionEnd0 - dn))) {
                                ++errorsInLastDnBases;
                            }

                            // Advance the position.
                            ++localPosition0;    
                            continue;
                        }

                        if (!hasBase0 && hasBase1) {       
                            // Increment the error count for this window.
                            ++windowErrors[windowIndex];

                            // Check if the position is in the first or last dn bases.
                            if ( ((segmentPositionStart0 + localPosition0 - positionStart0) <= dn) ) {
                                ++errorsInFirstDnBases;
                            } else if ( ((segmentPositionStart0 + localPosition0) >= (positionEnd0 - dn))) {
                                ++errorsInLastDnBases;
                            }

                            // Advance the position.
                            ++localPosition1;
                            continue;
                        }
                    }
                }

                uint64_t basesThatManagedToAlignAcrossWindows = 0;

                // Find windows with high error rate.
                for (uint64_t i = 0; i < numberOfWindows; i++) {
                    
                    // Calculate the number of errors per window threshold using a 7% error rate.
                    uint64_t numberOfErrorsPerWindowThreshold = windowBases[i] * 0.07;
                    
                    // If the window is small, set the threshold to 1.
                    if (numberOfErrorsPerWindowThreshold == 0 && windowBases[i] > 4) {
                        numberOfErrorsPerWindowThreshold = 1;
                    }
                    
                    // Cap the number of errors per window at 31.
                    if (numberOfErrorsPerWindowThreshold > 31) {
                        numberOfErrorsPerWindowThreshold = 31;
                    }
                    
                    // Edit distance threshold filtering.
                    if (windowErrors[i] > numberOfErrorsPerWindowThreshold) {
                        // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " window " << i << " has " << windowErrors[i] << " errors." << endl;
                        continue;
                    }

                    basesThatManagedToAlignAcrossWindows += windowBases[i];

                }

                // Check if the number of bases that managed to align across windows
                // that passed the sliding-window error filtering cover at least 90% of the total span.
                if (basesThatManagedToAlignAcrossWindows < 0.9 * span0) {
                    // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " too few bases that managed to align across windows." << endl;
                    continue;
                }

                // ---------------------------------------------------
                //   END OF: Sliding-window error filtering (read-0)    
                // ---------------------------------------------------





                


                
                

                // XXX
                // --- START OF: Max Gap filtering
                //   - Filter out alignments with too many consecutive gaps of size >= 6 bases length (alignmentMaxGaps = 64).
                //   - Filter out alignments with too high gap rate (alignmentMaxGapRate = 0.06).
                //

                // Add the hanging lengths to the total edit distance.
                // Treat the hanging lengths as unaligned bases.
                const uint64_t totalLength0 = projectedAlignment.totalLength[0];
                const uint64_t totalLength1 = projectedAlignment.totalLength[1];
                
                const uint64_t alignmentMaxGaps = 64;
                const double alignmentMaxGapRate = 0.06;

                uint64_t totalSVGaps = projectedAlignment.largeGapSum;
                uint64_t longestSVGap = std::max(projectedAlignment.longestGap[0], projectedAlignment.longestGap[1]);

                if(totalSVGaps > alignmentMaxGaps) {
                    // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " have a long SVGap: " << longestSVGap << endl;
                    // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " total SVGap sum (gaps >= 6): " << totalSVGaps << endl;
                    continue;
                }

                if((totalSVGaps > (alignmentMaxGapRate * double(totalLength0))) || (totalSVGaps > (alignmentMaxGapRate * double(totalLength1)))) {
                    // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " too high SVGap rate: " << totalSVGaps << " " << (alignmentMaxGapRate * double(totalLength0)) << " " << (alignmentMaxGapRate * double(totalLength1)) << endl;
                    continue;
                }

                // XXX
                // --- END OF: Max Gap filtering    
                //



                // XXX
                // --- START OF: Chemical Window Filtering
                //   Count for errors in the terminal regions of the read.
                //   - ONT reads often have one bad end
                //   - Chimeric reads may have one good end and one bad end
                //   - Terminal degradation is typically assymetric
                //   Check if in either end of the read:
                //   - We have more than min_err errors (min_err = 128)
                //   - The number of errors is greater than the error rate threshold (er = 0.36)
                //   If both conditions are met, the alignment is discarded.
                //

                double er = 0.36;
                uint64_t min_err = 128;

                if (runChemicalWindowFiltering) {
                    if ( (errorsInFirstDnBases > min_err) && (errorsInFirstDnBases > (dn * er)) ) {
                        // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " too many errors in first " << dn << " bases." << endl;
                        // cout << errorsInFirstDnBases << " " << dn * er << endl;
                        continue;
                    }

                    if ( (errorsInLastDnBases > min_err) && (errorsInLastDnBases > (dn * er)) ) {
                        // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " too many errors in last " << dn << " bases." << endl;
                        // cout << errorsInLastDnBases << " " << dn * er << endl;
                        continue;
                    }
                }

                // XXX
                // --- END OF: Chemical Window Filtering
                //
                
            
            
                alignmentInfo.errorRate = alignmentErrorRate;
                alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            }

            // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " good." << endl;
            threadAlignmentData.push_back(AlignmentData(candidate, alignmentInfo));

            // Store the alignment in compressed form.
            shasta::compress(alignment, compressedAlignment);
            thisThreadCompressedAlignments.appendVector(
                compressedAlignment.c_str(),
                compressedAlignment.c_str() + compressedAlignment.size()
            );
        }
    }

    if(alignmentMethod == 4) {
        std::lock_guard<std::mutex> lock(mutex);
        cout << "Thread " << threadId << " byte allocator: " <<
            byteAllocator.getMaxAllocatedByteCount() << "/" <<
            2ULL * 1024 * 1024 * 1024 << endl;
    }

    thisThreadCompressedAlignments.unreserve();

#if 0
    // Write the elapsed time taken by each alignment.
    ofstream csv("AlignTime-Thread" + to_string(threadId) + ".csv");
    for(const auto& p: elapsedTime) {
        const uint64_t i = p.first;
        const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
        SHASTA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);
        const OrientedReadId orientedReadId0 = OrientedReadId(candidate.readIds[0], 0);
        const OrientedReadId orientedReadId1 = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);
        const double t = p.second;
        csv << i << ",";
        csv << orientedReadId0 << ",";
        csv << orientedReadId1 << ",";
        csv << t << "\n";
    }
#endif
}



void Assembler::accessCompressedAlignments()
{
    compressedAlignments.accessExistingReadOnly(
        largeDataName("CompressedAlignments"));
}



// Compute alignmentTable from alignmentData.
// This could be made multithreaded if it becomes a bottleneck.
void Assembler::computeAlignmentTable()
{
    if(alignmentTable.isOpen()) {
        alignmentTable.remove();
    }
    alignmentTable.createNew(largeDataName("AlignmentTable"), largeDataPageSize);
    alignmentTable.beginPass1(ReadId(2 * reads->readCount()));
    for(const AlignmentData& ad: alignmentData) {
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
    }
    alignmentTable.beginPass2();
    for(uint32_t i=0; i<alignmentData.size(); i++) {
        const AlignmentData& ad = alignmentData[i];
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
    }
    alignmentTable.endPass2();



    // Sort each section of the alignment table by OrientedReadId.
    vector< pair<OrientedReadId, uint32_t> > v;
    for(ReadId readId0=0; readId0<reads->readCount(); readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);

            // Access the section of the alignment table for this oriented read.
            const span<uint32_t> alignmentTableSection =
                alignmentTable[orientedReadId0.getValue()];

            // Store pairs(OrientedReadId, alignmentIndex).
            v.clear();
            for(uint32_t alignmentIndex: alignmentTableSection) {
                const AlignmentData& alignment = alignmentData[alignmentIndex];
                const OrientedReadId orientedReadId1 = alignment.getOther(orientedReadId0);
                v.push_back(make_pair(orientedReadId1, alignmentIndex));
            }

            // Sort.
            sort(v.begin(), v.end());

            // Store the sorted alignmentIndex.
            for(size_t i=0; i<v.size(); i++) {
                alignmentTableSection[i] = v[i].second;
            }
        }
    }

    alignmentTable.unreserve();

}



void Assembler::accessAlignmentData()
{
    alignmentData.accessExistingReadOnly(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadOnly(largeDataName("AlignmentTable"));
}
void Assembler::accessAlignmentDataReadWrite()
{
    alignmentData.accessExistingReadWrite(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadWrite(largeDataName("AlignmentTable"));
}



void Assembler::checkAlignmentDataAreOpen() const
{
    if(!alignmentData.isOpen || !alignmentTable.isOpen()) {
        throw runtime_error("Alignment data are not accessible.");
    }
}



// Find in the alignment table the alignments involving
// a given oriented read, and return them with the correct
// orientation (this may involve a swap and/or reverse complement
// of the AlignmentInfo stored in the alignmentTable).
vector< pair<OrientedReadId, shasta::AlignmentInfo> >
    Assembler::findOrientedAlignments(
        OrientedReadId orientedReadId0Argument,
        bool inReadGraphOnly) const
{
    const ReadId readId0 = orientedReadId0Argument.getReadId();
    const ReadId strand0 = orientedReadId0Argument.getStrand();

    vector< pair<OrientedReadId, AlignmentInfo> > result;

    // Loop over alignment involving this read, as stored in the
    // alignment table.
    const auto alignmentTable0 = alignmentTable[orientedReadId0Argument.getValue()];
    for(const auto i: alignmentTable0) {
        const AlignmentData& ad = alignmentData[i];

        // Get the oriented read ids that the AlignmentData refers to.
        OrientedReadId orientedReadId0(ad.readIds[0], 0);
        OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        AlignmentInfo alignmentInfo = ad.info;

        // Skip it if it is not in the read graph and only alignments
        // in the read graph were requested.
        if(inReadGraphOnly and (not alignmentInfo.isInReadGraph)) {
            continue;
        }

        // Swap oriented reads, if necessary.
        if(orientedReadId0.getReadId() != readId0) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }
        SHASTA_ASSERT(orientedReadId0.getReadId() == readId0);

        // Reverse complement, if necessary.
        if(orientedReadId0.getStrand() != strand0) {
            orientedReadId0.flipStrand();
            orientedReadId1.flipStrand();
            alignmentInfo.reverseComplement();
        }
        SHASTA_ASSERT(orientedReadId0.getStrand() == strand0);
        SHASTA_ASSERT(orientedReadId0 == orientedReadId0Argument);

        result.push_back(make_pair(orientedReadId1, alignmentInfo));
    }
    return result;
}



// Flag palindromic reads.
void Assembler::flagPalindromicReads(
    uint32_t maxSkip,
    uint32_t maxDrift,
    uint32_t maxMarkerFrequency,
    double alignedFractionThreshold,
    double nearDiagonalFractionThreshold,
    uint32_t deltaThreshold,
    size_t threadCount)
{
    performanceLog << timestamp << "Finding palindromic reads." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the parameters so all threads can see them.
    flagPalindromicReadsData.maxSkip = maxSkip;
    flagPalindromicReadsData.maxDrift = maxDrift;
    flagPalindromicReadsData.maxMarkerFrequency = maxMarkerFrequency;
    flagPalindromicReadsData.alignedFractionThreshold = alignedFractionThreshold;
    flagPalindromicReadsData.nearDiagonalFractionThreshold = nearDiagonalFractionThreshold;
    flagPalindromicReadsData.deltaThreshold = deltaThreshold;

    // Reset all palindromic flags.
    reads->assertReadsAndFlagsOfSameSize();
    const ReadId readCount = reads->readCount();
    for(ReadId readId=0; readId<readCount; readId++) {
        reads->setPalindromicFlag(readId, false);
    }

    // Do it in parallel.
    setupLoadBalancing(readCount, 1000);
    runThreads(&Assembler::flagPalindromicReadsThreadFunction, threadCount);

    // Count the reads flagged as palindromic.
    size_t palindromicReadCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(reads->getFlags(readId).isPalindromic) {
            ++palindromicReadCount;
        }
    }
    assemblerInfo->palindromicReadCount = palindromicReadCount;
    cout << "Flagged " << palindromicReadCount <<
        " reads as palindromic out of " << readCount << " total." << endl;
    cout << "Palindromic fraction is " <<
        double(palindromicReadCount)/double(readCount) << endl;

}



void Assembler::flagPalindromicReadsThreadFunction(uint64_t)
{

    // Work areas used inside the loop and defined here
    // to reduce memory allocation activity.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;

    // Make local copies of the parameters.
    const uint32_t maxSkip = flagPalindromicReadsData.maxSkip;
    const uint32_t maxDrift = flagPalindromicReadsData.maxDrift;
    const uint32_t maxMarkerFrequency = flagPalindromicReadsData.maxMarkerFrequency;
    const double alignedFractionThreshold = flagPalindromicReadsData.alignedFractionThreshold;
    const double nearDiagonalFractionThreshold = flagPalindromicReadsData.nearDiagonalFractionThreshold;
    const uint32_t deltaThreshold = flagPalindromicReadsData.deltaThreshold;


    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    reads->assertReadsAndFlagsOfSameSize();

    while(getNextBatch(begin, end)) {

        // Loop over all reads in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Get markers sorted by KmerId for this read and its reverse complement.
            for(Strand strand=0; strand<2; strand++) {
                getMarkersSortedByKmerId(OrientedReadId(readId, strand), markersSortedByKmerId[strand]);
            }

            // Compute a marker alignment of this read versus its reverse complement.
            alignOrientedReads(markersSortedByKmerId, maxSkip, maxDrift, maxMarkerFrequency, false,
                graph, alignment, alignmentInfo);

            // If the alignment has too few markers, skip it.
            const size_t alignedMarkerCount = alignment.ordinals.size();
            const size_t totalMarkerCount = markersSortedByKmerId[0].size();
            const double alignedFraction = double(alignedMarkerCount)/double(totalMarkerCount);
            if(alignedFraction < alignedFractionThreshold) {
                continue;
            }

            // If the alignment has too few markers near the diagonal, skip it.
            size_t nearDiagonalMarkerCount = 0;
            for(size_t i=0; i<alignment.ordinals.size(); i++) {
                const array<uint32_t, 2>& ordinals = alignment.ordinals[i];
                const int32_t ordinal0 = int32_t(ordinals[0]);
                const int32_t ordinal1 = int32_t(ordinals[1]);
                const uint32_t delta = abs(ordinal0 - ordinal1);
                if(delta < deltaThreshold) {
                    nearDiagonalMarkerCount++;
                }
            }
            const double nearDiagonalFraction = double(nearDiagonalMarkerCount)/double(totalMarkerCount);
            if(nearDiagonalFraction < nearDiagonalFractionThreshold) {
                continue;
            }

            // If we got here, mark the read as palindromic.
            reads->setPalindromicFlag(readId, true);

        }
    }
}



void Assembler::analyzeAlignmentMatrix(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    // Get the oriented reads.
    const OrientedReadId orientedReadId0(readId0, strand0);
    const OrientedReadId orientedReadId1(readId1, strand1);

    // Get the markers sorted by kmerId.
    vector<MarkerWithOrdinal> markers0;
    vector<MarkerWithOrdinal> markers1;
    getMarkersSortedByKmerId(orientedReadId0, markers0);
    getMarkersSortedByKmerId(orientedReadId1, markers1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();

    // The number of markers in each oriented read.
    const int64_t n0 = end0 - begin0;
    const int64_t n1 = end1 - begin1;

    // We will use coordinates
    // x = ordinal0 + ordinal1
    // y = ordinal0 - ordinal1 (offset)
    // In these coordinates, diagonals in the alignment matrix
    // are lines of constant y and so they become horizontal.
    const int64_t xMin = 0;
    const int64_t xMax = n0 + n1 - 2;
    const int64_t yMin = -n1;
    const int64_t yMax = n0 - 1;
    const int64_t nx= xMax -xMin + 1;
    const int64_t ny= yMax -yMin + 1;

    // Create a histogram in cells of size (dx, dy).
    const int64_t dx = 100;
    const int64_t dy = 20;
    const int64_t nxCells = (nx-1)/dx + 1;
    const int64_t nyCells = (ny-1)/dy + 1;
    vector< vector<uint64_t> > histogram(nxCells, vector<uint64_t>(nyCells, 0));
    cout << "nxCells " << nxCells << endl;
    cout << "nyCells " << nyCells << endl;



    // Joint loop over the markers, looking for common k-mer ids.
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }


            // Loop over pairs in the streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                const int64_t ordinal0 = jt0->ordinal;
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int64_t ordinal1 = int64_t(jt1->ordinal);

                    const int64_t x = ordinal0 + ordinal1;
                    const int64_t y = ordinal0 - ordinal1;
                    const int64_t ix = (x-xMin) / dx;
                    const int64_t iy = (y-yMin) / dy;
                    SHASTA_ASSERT(ix >= 0);
                    SHASTA_ASSERT(iy >= 0);
                    SHASTA_ASSERT(ix < nxCells);
                    SHASTA_ASSERT(iy < nyCells);

                    ++histogram[ix][iy];
                }
            }


            // Continue joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }
    }

    ofstream csv("Histogram.csv");
    PngImage image = PngImage(int(nxCells), int(nyCells));
    uint64_t activeCellCount = 0;
    for(int64_t iy=0; iy<nyCells; iy++) {
        for(int64_t ix=0; ix<nxCells; ix++) {
            const int64_t frequency = histogram[ix][iy];
            if(frequency > 0) {
                csv << frequency;
            }
            if(frequency >= 10) {
                ++activeCellCount;
            }
            // const int r = (frequency <= 10) ? 0 : min(int(255), int(10*frequency));
            const int r = (frequency>=10) ? 255 : 0;
            // const int r = min(int(255), int(10*frequency));
            const int g = r;
            const int b = r;
            image.setPixel(int(ix), int(iy), r, g, b);
            csv << ",";
        }
        csv << "\n";
    }
    image.write("Histogram.png");
    cout << activeCellCount << " active cells out of " << nxCells*nyCells << endl;


}


// Count the common marker near a given ordinal offset for
// two oriented reads. This can be used to check
// whether an alignmnent near the specified ordinal offset exists.
uint32_t Assembler::countCommonMarkersNearOffset(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int32_t offset,
    int32_t offsetTolerance
)
{
    const int32_t minOffset = offset - offsetTolerance;
    const int32_t maxOffset = offset + offsetTolerance;
    return countCommonMarkersWithOffsetIn(
        orientedReadId0, orientedReadId1,
        minOffset, maxOffset);
}
uint32_t Assembler::countCommonMarkersWithOffsetIn(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int32_t minOffset,
    int32_t maxOffset
)
{
#if 1
    const bool debug = false;
#else
    const bool debug =
        orientedReadId0.getReadId() == 5 and
        orientedReadId0.getStrand() == 0 and
        orientedReadId1.getReadId() == 32 and
        orientedReadId1.getStrand() == 1;
#endif
    if(debug) {
        cout << "countCommonMarkersWithOffsetIn" << endl;
    }

    // Get the markers sorted by kmerId.
    checkMarkersAreOpen();
    vector<MarkerWithOrdinal> markers0;
    vector<MarkerWithOrdinal> markers1;
    getMarkersSortedByKmerId(orientedReadId0, markers0);
    getMarkersSortedByKmerId(orientedReadId1, markers1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();



    // Main loop, looking for common markers.
    uint32_t count = 0;
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common marker.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in markers0 and markers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }

            // Loop over pairs of markers in the streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int32_t ordinalOffset = int32_t(jt0->ordinal) - int32_t(jt1->ordinal);
                    if(debug) {
                        cout << jt0->ordinal << " " <<
                            jt1->ordinal << " " <<
                            ordinalOffset<< endl;
                    }
                    if(ordinalOffset >= minOffset and ordinalOffset <= maxOffset) {
                        ++count;
                    }

                }
            }
            // Continue main loop.
            it0 = it0End;
            it1 = it1End;
        }

    }

    return count;
}



// Check if an alignment between two reads should be suppressed,
// bases on the setting of command line option
// --Align.sameChannelReadAlignment.suppressDeltaThreshold.
bool Assembler::suppressAlignment(
    ReadId readId0,
    ReadId readId1,
    uint64_t delta)
{

    // If the ch meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    // Check the channel first for efficiency,
    // so we can return faster in most cases.
    const auto ch0 = reads->getMetaData(readId0, "ch");
    if(ch0.empty()) {
        return false;
    }
    const auto ch1 = reads->getMetaData(readId1, "ch");
    if(ch1.empty()) {
        return false;
    }
    if(ch0 not_eq ch1) {
        return false;
    }



    // If the sampleid meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    const auto sampleid0 = reads->getMetaData(readId0, "sampleid");
    if(sampleid0.empty()) {
        return false;
    }
    const auto sampleid1 = reads->getMetaData(readId1, "sampleid");
    if(sampleid1.empty()) {
        return false;
    }
    if(sampleid0 not_eq sampleid1) {
        return false;
    }



    // If the runid meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    const auto runid0 = reads->getMetaData(readId0, "runid");
    if(runid0.empty()) {
        return false;
    }
    const auto runid1 = reads->getMetaData(readId1, "runid");
    if(runid1.empty()) {
        return false;
    }
    if(runid0 not_eq runid1) {
        return false;
    }

    // cout << "Checking " << readId0 << " " << readId1 << endl;


    // If the read meta data fields of the two reads are missing,
    // don't suppress the alignment.
    const auto read0 = reads->getMetaData(readId0, "read");
    if(read0.empty()) {
        return false;
    }
    const auto read1 = reads->getMetaData(readId1, "read");
    if(read1.empty()) {
        return false;
    }
    // cout << read0 << " " << read1 << endl;



    // Convert the read meta data fields to integers.
    // Keep in mind the span<char> is not null-terminated.
    const int64_t r0 = int64_t(atoul(read0));
    const int64_t r1 = int64_t(atoul(read1));
    // cout << r0 << " " << r1 << endl;



    // Suppress the alignment if the absolute difference of the
    // read meta data fields is less than delta.
    return std::abs(r0 - r1) < int64_t(delta);

}



// Remove all alignment candidates for which suppressAlignment
// returns false.
void Assembler::suppressAlignmentCandidates(
    uint64_t delta,
    size_t threadCount)
{
    performanceLog << timestamp << "Suppressing alignment candidates." << endl;

    // Allocate memory for flags to keep track of which alignments
    // should be suppressed.
    suppressAlignmentCandidatesData.suppress.createNew(
        largeDataName("tmp-suppressAlignmentCandidates"), largeDataPageSize);
    const uint64_t candidateCount = alignmentCandidates.candidates.size();
    suppressAlignmentCandidatesData.suppress.resize(candidateCount);

    // Figure out which candidates should be suppressed.
    suppressAlignmentCandidatesData.delta = delta;
    const uint64_t batchSize = 10000;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::suppressAlignmentCandidatesThreadFunction, threadCount);

    ofstream csv("SuppressedAlignmentCandidates.csv");
    csv << "ReadId0,ReadId1,SameStrand,Name0,Name1,MetaData0,MetaData1" << endl;

    // Suppress the alignment candidates we flagged.
    cout << "Number of alignment candidates before suppression is " << candidateCount << endl;
    uint64_t j = 0;
    uint64_t suppressCount = 0;
    for(uint64_t i=0; i<candidateCount; i++) {
        if(suppressAlignmentCandidatesData.suppress[i]) {
            ++suppressCount;
            const ReadId readId0 = alignmentCandidates.candidates[i].readIds[0];
            const ReadId readId1 = alignmentCandidates.candidates[i].readIds[1];
            csv << readId0 << "," << readId1 << ","
                << (alignmentCandidates.candidates[i].isSameStrand ? "Yes" : "No") << ","
                << reads->getReadName(readId0) << "," << reads->getReadName(readId1) << ","
                << reads->getReadMetaData(readId0) << "," << reads->getReadMetaData(readId1) << endl;
        } else {
            alignmentCandidates.candidates[j++] =
                alignmentCandidates.candidates[i];
        }
    }
    SHASTA_ASSERT(j + suppressCount == candidateCount);
    alignmentCandidates.candidates.resize(j);
    cout << "Suppressed " << suppressCount << " alignment candidates." << endl;
    cout << "Number of alignment candidates after suppression is " << j << endl;


    // Clean up.
    suppressAlignmentCandidatesData.suppress.remove();

    performanceLog << timestamp << "Done suppressing alignment candidates." << endl;
}



void Assembler::suppressAlignmentCandidatesThreadFunction(uint64_t)
{
    const uint64_t delta = suppressAlignmentCandidatesData.delta;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate alignments in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const OrientedReadPair& p = alignmentCandidates.candidates[i];
            // cout << "Checking " << p.readIds[0] << " " << p.readIds[1] <<  " " << int(p.isSameStrand) << endl;
            suppressAlignmentCandidatesData.suppress[i] =
                suppressAlignment(p.readIds[0], p.readIds[1], delta);
        }

    }
}


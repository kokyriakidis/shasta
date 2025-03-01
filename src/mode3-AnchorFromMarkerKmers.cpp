#include "mode3-Anchor.hpp"
#include "MarkerKmers.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

using namespace shasta;
using namespace mode3;



Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    shared_ptr<MarkerKmers> markerKmers,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    kHalf = k / 2;

    cout << timestamp << "Anchor creation from marker kmers begins." << endl;

    // Store arguments so all threads can see them.
    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;
    data.markerKmers = markerKmers;

    // Use the MarkerKmers to construct anchors.
    // Each thread stores the anchors it finds in a separate vector.
    data.threadAnchors.clear(); // Just in case
    data.threadAnchors.resize(threadCount);
    const uint64_t batchSize = 1000;
    setupLoadBalancing(markerKmers->size(), batchSize);
    runThreads(&Anchors::constructFromMarkerKmersThreadFunction, threadCount);

    // Gather the anchors found by all threads.
    performanceLog << timestamp << "Gathering the anchors found by all threads." << endl;
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadAnchorsPointer = data.threadAnchors[threadId];
        auto& threadAnchors = *threadAnchorsPointer;

        // Loop over the anchors found by this thread.
        for(uint64_t i=0; i<threadAnchors.size(); i++) {
            const auto threadAnchor = threadAnchors[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& markerInfo: threadAnchor) {
                anchorMarkerIntervals.append(AnchorMarkerInterval(markerInfo.orientedReadId, markerInfo.ordinal));
            }
        }
        threadAnchors.remove();
        threadAnchorsPointer = 0;
    }

    // The anchor sequences are all empty.
    // The ordinal offsets are all 0.
    const uint64_t anchorCount = anchorMarkerIntervals.size();
    anchorInfos.resize(anchorCount);
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        anchorSequences.appendVector();
        anchorInfos[anchorId].ordinalOffset = 0;
    }

    cout << "Generated " << anchorCount << " anchors." << endl;
    cout << timestamp << "Anchor creation from marker kmers ends." << endl;
}



void Anchors::constructFromMarkerKmersThreadFunction(uint64_t threadId)
{
    using MarkerInfo = MarkerKmers::MarkerInfo;

    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;
    const MarkerKmers& markerKmers = *data.markerKmers;

    // Initialize the anchors that will be found by this thread.
    auto& threadAnchorsPointer = data.threadAnchors[threadId];
    threadAnchorsPointer = make_shared<MemoryMapped::VectorOfVectors<ConstructFromMarkerKmersData::MarkerInfo, uint64_t> >();
    auto& threadAnchors = *threadAnchorsPointer;
    threadAnchors.createNew(largeDataName("tmp-threadAnchors-" + to_string(threadId)), largeDataPageSize);

    // A vector used below and defined here to reduce memory allocation activity.
    // It will contain the MarkerInfos for a marker Kmer, ecluding
    // the ones for which the same ReadId appears more than once in the same Kmer.
    // There are the ones that will be used to generate anchors.
    vector<MarkerInfo> usableMarkerInfos;

    // Loop over batches of marker Kmers assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker Kmers assigned to this batch.
        for(uint64_t i=begin; i!=end; i++) {

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[i];

            // Check for high coverage using all of the marker infos.
            if(markerInfos.size() > maxAnchorCoverage) {
                continue;
            }



            // Gather the usable MarkerInfos.
            // These are the ones for which the ReadId is different from the ReadId
            // of the previous and next MarkerInfo.
            usableMarkerInfos.clear();
            for(uint64_t i=0; i<markerInfos.size(); i++) {
                const MarkerInfo& markerInfo = markerInfos[i];
                bool isUsable = true;

                // Check if same ReadId of previous MarkerInfo.
                if(i != 0) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i-1].orientedReadId.getReadId());
                }

                // Check if same ReadId of next MarkerInfo.
                if(i != markerInfos.size() - 1) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i+1].orientedReadId.getReadId());
                }

                if(isUsable) {
                    usableMarkerInfos.push_back(markerInfo);
                }
            }



            // Check for low coverage using the usable marker infos.
            if(usableMarkerInfos.size() < minAnchorCoverage) {
                continue;
            }



            // If getting here, we will generate a pair of Anchors corresponding
            // to this Kmer.


            // Generate the first anchor of the pair (no reverse complementing).
            threadAnchors.appendVector(usableMarkerInfos);


            // Reverse complement the usableMarkerInfos, then
            // generate the second anchor in the pair.
            for(MarkerInfo& markerInfo: usableMarkerInfos) {
                markerInfo = markerKmers.reverseComplement(markerInfo);
            }
            threadAnchors.appendVector(usableMarkerInfos);
        }
    }

}


#if 0
// Alignment-free generation of anchors from marker k-mers (old code).
#include "mode3-Anchor.hpp"
#include "Base.hpp"
#include "extractKmer.hpp"
#include "Marker.hpp"
#include "markerAccessFunctions.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ShortBaseSequenceEdit.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

#include "fstream.hpp"


Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    kHalf = k / 2;

    cout << timestamp << "Anchor creation from marker kmers begins." << endl;

    // Store coverage thresholds so all threads can see them.
    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    data.minPrimaryCoverage = minPrimaryCoverage;
    data.maxPrimaryCoverage = maxPrimaryCoverage;

    // We want to construct a hash table that will contain a MarkerInfo object
    // for each marker in all oriented reads.
    // Indexed by bucketId.
    // The bucket is computed by hashing the k-mer of each marker,
    // so all markers with the same k-mer end up in the same bucket.

    // Figure out the number of buckets in the hash table.
    // It will contain one entry for each marker.
    const uint64_t totalMarkerCount = markers.totalSize();
    const uint64_t bucketCount = (std::bit_ceil(totalMarkerCount) >> 4);
    data.mask = bucketCount - 1;
    cout << "Total number of markers " << totalMarkerCount << endl;
    cout << "Number of buckets " << bucketCount << endl;

    // Initialize the hash table.
    auto& buckets = data.buckets;
    buckets.createNew(largeDataName("tmp-constructFromMarkerKmersData-Buckets"), largeDataPageSize);

    // Gather markers and store them in buckets.
    buckets.beginPass1(bucketCount);
    size_t batchSize = 10;
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&Anchors::constructFromMarkerKmersGatherMarkersPass1, threadCount);
    buckets.beginPass2();
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&Anchors::constructFromMarkerKmersGatherMarkersPass2, threadCount);
    buckets.endPass2(true, true);

    // Compute k-mer frequencies.
    data.kmerInfo.createNew(largeDataName("tmp-constructFromMarkerKmersData-KmerInfo"), largeDataPageSize);
    data.kmerInfo.beginPass1(bucketCount);
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&Anchors::constructFromMarkerKmersComputeKmerFrequencyPass1, threadCount);
    data.kmerInfo.beginPass2();
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&Anchors::constructFromMarkerKmersComputeKmerFrequencyPass2, threadCount);
    data.kmerInfo.endPass2(true, true);

#if 1
    // Write out the high frequency k-mers.
    {
        vector< pair<Kmer, uint64_t> > kmerFrequency;
        for(uint64_t bucketId=0; bucketId<bucketCount; bucketId++) {
            const auto bucket = data.kmerInfo[bucketId];
            for(const auto& kmerInfo: bucket) {
                if(kmerInfo.frequency > maxPrimaryCoverage) {
                    kmerFrequency.push_back(
                        make_pair(kmerInfo.kmer, kmerInfo.frequency));
                }
            }
        }

        sort(kmerFrequency.begin(), kmerFrequency.end(),
            OrderPairsBySecondOnlyGreater<Kmer, uint64_t>());

        ofstream csv("HighFrequencyKmers.csv");
        for(const auto& p: kmerFrequency) {
            const Kmer& kmer = p.first;
            const uint64_t frequency = p.second;
            kmer.write(csv, k);
            csv << "," << frequency << "\n";
        }
    }
#endif

    // Flag forbidden k-mers.
    performanceLog << timestamp << "Flagging forbidden k-mers." << endl;
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&Anchors::constructFromMarkerKmersFlagForbiddenKmers, threadCount);


    // Process each bucket to create anchors.
    // Each thread stores for each anchor a vector of (OrientedReadId, ordinal)
    // (data.threadAnchors).
    performanceLog << timestamp << "Creating anchors." << endl;
    data.threadAnchors.resize(threadCount);
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&Anchors::constructFromMarkerKmersCreateAnchors, threadCount);

    // We no longer need the hash tables.
    buckets.remove();
    data.kmerInfo.remove();



    // Gather the anchors found by all threads.
    performanceLog << timestamp << "Gathering the anchors found by all threads." << endl;
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadAnchorsPointer = data.threadAnchors[threadId];
        auto& threadAnchors = *threadAnchorsPointer;

        // Loop over the anchors found by this thread.
        for(uint64_t i=0; i<threadAnchors.size(); i++) {
            const auto threadAnchor = threadAnchors[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& markerInfo: threadAnchor) {
                anchorMarkerIntervals.append(AnchorMarkerInterval(markerInfo.orientedReadId, markerInfo.ordinal));
            }
        }
        threadAnchors.remove();
        threadAnchorsPointer = 0;
    }

    // The anchor sequences are all empty.
    // The ordinal offsets are all 0.
    const uint64_t anchorCount = anchorMarkerIntervals.size();
    anchorInfos.resize(anchorCount);
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        anchorSequences.appendVector();
        anchorInfos[anchorId].ordinalOffset = 0;
    }

    cout << "Generated " << anchorCount << " anchors." << endl;
    cout << timestamp << "Anchor creation from marker kmers ends." << endl;
}



void Anchors::constructFromMarkerKmersGatherMarkersPass1(uint64_t /* threadId */)
{
    constructFromMarkerKmersGatherMarkersPass12(1);
}



void Anchors::constructFromMarkerKmersGatherMarkersPass2(uint64_t /* threadId */)
{
    constructFromMarkerKmersGatherMarkersPass12(2);
}



void Anchors::constructFromMarkerKmersGatherMarkersPass12(uint64_t pass)
{
    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    const uint64_t mask = data.mask;
    auto& buckets = data.buckets;

    ConstructFromMarkerKmersData::MarkerInfo markerInfo;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); ++readId) {
            const LongBaseSequenceView readSequence = reads.getRead(readId);

            // Get the markers on this read, without reverse complementing (strand 0).
            OrientedReadId orientedReadId(readId, 0);
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const uint32_t orientedReadMarkerCount = uint32_t(orientedReadMarkers.size());

            // Loop  over the markers.
            for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
                const CompressedMarker& marker = orientedReadMarkers[ordinal];
                const uint32_t position = marker.position;

                // Extract the Kmer at this position.
                Kmer kmer;
                extractKmer(readSequence, position, k, kmer);

                // Find the bucket this marker goes to.
                const uint64_t hashValue = MurmurHash64A(&kmer, sizeof(kmer), 1241);
                const uint64_t bucketId = hashValue & mask;

                // Increment the bucket count (pass 1) or store this marker (pass 2).
                if(pass == 1) {
                    buckets.incrementCountMultithreaded(bucketId);
                } else {
                    markerInfo.orientedReadId = orientedReadId;
                    markerInfo.ordinal = ordinal;
                    buckets.storeMultithreaded(bucketId, markerInfo);
                }

                // Do the same for the reverse complemented k-mer.
                const Kmer kmerRc = kmer.reverseComplement(k);
                const uint64_t hashValueRc = MurmurHash64A(&kmerRc, sizeof(kmerRc), 1241);
                const uint64_t bucketIdRc = hashValueRc & mask;
                if(pass == 1) {
                    buckets.incrementCountMultithreaded(bucketIdRc);
                } else {
                    markerInfo.orientedReadId.flipStrand();
                    markerInfo.ordinal = orientedReadMarkerCount - 1 - markerInfo.ordinal;
                    buckets.storeMultithreaded(bucketIdRc, markerInfo);
                }
            }
        }
    }
}


void Anchors::constructFromMarkerKmersComputeKmerFrequencyPass1(uint64_t /* threadId */)
{
    constructFromMarkerKmersComputeKmerFrequencyPass12(1);
}



void Anchors::constructFromMarkerKmersComputeKmerFrequencyPass2(uint64_t /* threadId */)
{
    constructFromMarkerKmersComputeKmerFrequencyPass12(2);
}



// Compute the frequency of every k-mer.
void Anchors::constructFromMarkerKmersComputeKmerFrequencyPass12(uint64_t pass)
{
    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    auto& buckets = data.buckets;

    // Vector to contain the MarkerInfos found in a bucket,
    // together with the corresponding marker k-mers.
    // Sorted by k-mer.
    class MarkerInfo : public ConstructFromMarkerKmersData::MarkerInfo {
    public:
        Kmer kmer;
        MarkerInfo(const ConstructFromMarkerKmersData::MarkerInfo& base, const Kmer& kmer) :
            ConstructFromMarkerKmersData::MarkerInfo(base), kmer(kmer) {}
        bool operator<(const MarkerInfo& that) const
        {
            return tie(kmer, orientedReadId, ordinal) < tie(that.kmer, that.orientedReadId, that.ordinal);
        }
    };
    vector<MarkerInfo> bucketInfo;
    vector< pair<Kmer, uint64_t> > kmersWithFrequency;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all buckets assigned to this batch.
        for(uint64_t bucketId=begin; bucketId!=end; ++bucketId) {

            // Gather the markers in this bucket with the corresponding Kmers.
            auto bucket = buckets[bucketId];
            bucketInfo.clear();
            for(const ConstructFromMarkerKmersData::MarkerInfo& markerInfo: bucket) {
                const Kmer kmer = getOrientedReadMarkerKmer(
                    markerInfo.orientedReadId,
                    markerInfo.ordinal,
                    k, reads, markers);
                bucketInfo.push_back(MarkerInfo(markerInfo, kmer));
            }

            // Sort them by Kmer.
            sort(bucketInfo.begin(), bucketInfo.end());



            // Process each streak with the same k-mer.
            kmersWithFrequency.clear();
            for(uint64_t streakBegin=0; streakBegin<bucketInfo.size(); /* Increment later */) {
                const Kmer& kmer = bucketInfo[streakBegin].kmer;

                // Find the end of the streak with this k-mer.
                uint64_t streakEnd = streakBegin + 1;
                while(true) {
                    if((streakEnd == bucketInfo.size()) or ( bucketInfo[streakEnd].kmer != kmer)) {
                        break;
                    }
                    ++streakEnd;
                }

                kmersWithFrequency.push_back(make_pair(kmer, streakEnd - streakBegin));

                // Prepare to process the next streak.
                streakBegin = streakEnd;
            }

            if(pass == 1) {
                data.kmerInfo.incrementCountMultithreaded(bucketId, kmersWithFrequency.size());
            } else {
                for(const auto& p: kmersWithFrequency) {
                    const Kmer& kmer = p.first;
                    const uint64_t frequency = p.second;
                    data.kmerInfo.storeMultithreaded(bucketId,
                        ConstructFromMarkerKmersData::KmerInfo(kmer, frequency));
                }
            }
        }
    }

}



void Anchors::constructFromMarkerKmersFlagForbiddenKmers(uint64_t /* threadId */)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t veryHighFrequencyThreshold = 1000;

    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    const uint64_t mask = data.mask;


    vector<Kmer> editedKmers;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all buckets assigned to this batch.
        for(uint64_t bucketId0=begin; bucketId0!=end; ++bucketId0) {
            const auto bucket0 = constructFromMarkerKmersData.kmerInfo[bucketId0];

            // Loop over all kmers in this bucket that have very high frequency.
            for(const auto& kmerInfo: bucket0) {
                if(kmerInfo.frequency < veryHighFrequencyThreshold) {
                    continue;
                }

                const Kmer& kmer0 = kmerInfo.kmer;

                // Apply 1-base edits.
                applySingleEdit(kmer0, k, editedKmers);

                // The edited k-mers are marked as forbidden, if present.
                for(const Kmer& kmer1: editedKmers) {

                    // Locate the bucket where kmer1 is, if present.
                    const uint64_t hashValue = MurmurHash64A(&kmer1, sizeof(kmer1), 1241);
                    const uint64_t bucketId1 = hashValue & mask;
                    const auto bucket1 = constructFromMarkerKmersData.kmerInfo[bucketId1];

                    // If kmer1 is present, mark it as forbidden.
                    for(auto& kmerInfo1: bucket1) {
                        if(kmerInfo1.kmer == kmer1) {
                            kmerInfo1.isForbidden = true;
                            break;
                        }
                    }

                }
            }

        }
    }

}



// This loops over buckets and creates anchors from the MarkerInfo
// objects stored in each bucket.
void Anchors::constructFromMarkerKmersCreateAnchors(uint64_t threadId )
{

    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    const uint64_t minPrimaryCoverage = data.minPrimaryCoverage;
    const uint64_t maxPrimaryCoverage = data.maxPrimaryCoverage;
    auto& buckets = data.buckets;

    // Initialize the anchors that will be found by this thread.
    auto& threadAnchorsPointer = data.threadAnchors[threadId];
    threadAnchorsPointer = make_shared<MemoryMapped::VectorOfVectors<ConstructFromMarkerKmersData::MarkerInfo, uint64_t> >();
    auto& threadAnchors = *threadAnchorsPointer;
    threadAnchors.createNew(largeDataName("tmp-threadAnchors-" + to_string(threadId)), largeDataPageSize);

    // Vector to contain the MarkerInfos found in a bucket,
    // together with the corresponding marker k-mers.
    // Sorted by k-mer.
    class MarkerInfo : public ConstructFromMarkerKmersData::MarkerInfo {
    public:
        Kmer kmer;
        MarkerInfo(const ConstructFromMarkerKmersData::MarkerInfo& base, const Kmer& kmer) :
            ConstructFromMarkerKmersData::MarkerInfo(base), kmer(kmer) {}
        bool operator<(const MarkerInfo& that) const
        {
            return tie(kmer, orientedReadId, ordinal) < tie(that.kmer, that.orientedReadId, that.ordinal);
        }
    };
    vector<MarkerInfo> bucketInfo;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all buckets assigned to this batch.
        for(uint64_t bucketId=begin; bucketId!=end; ++bucketId) {

            // Gather the markers in this bucket with the corresponding Kmers.
            // If the bucket is too small it cannot generate any anchors.
            auto bucket = buckets[bucketId];
            if(bucket.size() < minPrimaryCoverage) {
                continue;
            }
            bucketInfo.clear();
            for(const ConstructFromMarkerKmersData::MarkerInfo& markerInfo: bucket) {
                const Kmer kmer = getOrientedReadMarkerKmer(
                    markerInfo.orientedReadId,
                    markerInfo.ordinal,
                    k, reads, markers);
                bucketInfo.push_back(MarkerInfo(markerInfo, kmer));
            }

            // Sort them by Kmer.
            sort(bucketInfo.begin(), bucketInfo.end());



            // Each streak with the same k-mer generates a pair of anchors, but only if:
            // - Its coverage is in [minPrimaryCoverage, maxPrimaryCoverage].
            // - The first OrientedReadId is on strand 0, AND
            // - There are no repeated ReadIds.
            for(uint64_t streakBegin=0; streakBegin<bucketInfo.size(); /* Increment later */) {
                const Kmer& kmer = bucketInfo[streakBegin].kmer;

                // Find the end of the streak with this k-mer.
                uint64_t streakEnd = streakBegin + 1;
                while(true) {
                    if((streakEnd == bucketInfo.size()) or ( bucketInfo[streakEnd].kmer != kmer)) {
                        break;
                    }
                    ++streakEnd;
                }

#if 0
                // If this k-mer is marked as forbidden, skip it.
                const auto& kmerInfos = data.kmerInfo[bucketId];
                bool isForbidden = false;
                for(const auto& kmerInfo: kmerInfos) {
                    if(kmerInfo.isForbidden) {
                        isForbidden = true;
                        break;
                    }
                }
                if(isForbidden) {
                    // Prepare to process the next streak.
                    streakBegin = streakEnd;
                    continue;
                }
#endif

                // If the first oriented read is not on strand 0, skip it.
                if(bucketInfo[streakBegin].orientedReadId.getStrand() != 0) {
                    // Prepare to process the next streak.
                    streakBegin = streakEnd;
                    continue;
                }

                // If coverage is not in the desired range, skip it.
                const uint64_t coverage = streakEnd - streakBegin;
                if((coverage < minPrimaryCoverage) or (coverage > maxPrimaryCoverage)) {
                    // Prepare to process the next streak.
                    streakBegin = streakEnd;
                    continue;
                }

                // If there are repeated ReadIds, skip it.
                bool repeatedReadIdsExist = false;
                for(uint64_t i=streakBegin+1; i<streakEnd; i++) {
                    if(bucketInfo[i-1].orientedReadId.getReadId() == bucketInfo[i].orientedReadId.getReadId()) {
                        repeatedReadIdsExist = true;
                        break;
                    }
                }
                if(repeatedReadIdsExist) {
                    // Prepare to process the next streak.
                    streakBegin = streakEnd;
                    continue;
                }

#if 0
                // EXPOSE THE CONSTANTS WHEN CODE STABILIZES
                // If this k-mer has long repeats, skip it.
                // Do this check last because it is expensive.
                if(
                    (kmer.maxHomopolymerLength(k) > 6) or
                    (kmer.countExactRepeatCopies<2>(k) > 6) or
                    (kmer.countExactRepeatCopies<3>(k) > 4) or
                    (kmer.countExactRepeatCopies<4>(k) > 3) or
                    (kmer.countExactRepeatCopies<5>(k) > 3) or
                    (kmer.countExactRepeatCopies<6>(k) > 3)
                    ) {
                    // Prepare to process the next streak.
                    streakBegin = streakEnd;
                    continue;
                }
#endif

                // Generate the first anchor of the pair (no reverse complementing).
                threadAnchors.appendVector();
                for(uint64_t i=streakBegin; i!=streakEnd; i++) {
                    const MarkerInfo& markerInfo = bucketInfo[i];
                    threadAnchors.append(markerInfo);
                }

                // Generate the second anchor of the pair (with reverse complementing).
                threadAnchors.appendVector();
                for(uint64_t i=streakBegin; i!=streakEnd; i++) {
                    MarkerInfo markerInfo = bucketInfo[i];
                    const uint32_t markerCount = uint32_t(markers[markerInfo.orientedReadId.getValue()].size());
                    markerInfo.orientedReadId.flipStrand();
                    markerInfo.ordinal = markerCount - 1 - markerInfo.ordinal;
                    threadAnchors.append(markerInfo);
                }


                // Prepare to process the next streak.
                streakBegin = streakEnd;
            }
        }
    }


}
#endif

// Alignment-free generation of anchors from marker k-mers.
#include "mode3-Anchor.hpp"
#include "Base.hpp"
#include "extractKmer.hpp"
#include "Marker.hpp"
#include "markerAccessFunctions.hpp"
#include "MurmurHash2.hpp"
#include "Reads.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;



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
    buckets.createNew(largeDataName("tmp-constructFromMarkerKmersData"), largeDataPageSize);

    // Pass 1 counts the markers that go to each bucket.
    buckets.beginPass1(bucketCount);
    size_t batchSize = 10;
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&Anchors::constructFromMarkerKmersThreadFunction1, threadCount);

    // In pass 2 we store each marker in its assigned bucket.
    buckets.beginPass2();
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&Anchors::constructFromMarkerKmersThreadFunction2, threadCount);
    buckets.endPass2(true, true);

    // In pass 3 we process each bucket to create anchors.
    // Each thread stores for each anchor a vector of (OrientedReadId, ordinal)
    // (data.threadAnchors).
    data.threadAnchors.resize(threadCount);
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&Anchors::constructFromMarkerKmersThreadFunction3, threadCount);

    // We no longer need the hash table.
    buckets.remove();



    // Gather the anchors found by all threads.
    cout << timestamp << "Gathering the anchors found by all threads." << endl;
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
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {

    }

    SHASTA_ASSERT(0);

   SHASTA_ASSERT(0);

    cout << timestamp << "Anchor creation from marker kmers ends." << endl;

    SHASTA_ASSERT(0);
}



void Anchors::constructFromMarkerKmersThreadFunction1(uint64_t /* threadId */)
{
    constructFromMarkerKmersThreadFunction12(1);
}



void Anchors::constructFromMarkerKmersThreadFunction2(uint64_t /* threadId */)
{
    constructFromMarkerKmersThreadFunction12(2);
}



void Anchors::constructFromMarkerKmersThreadFunction12(uint64_t pass)
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

            // Loop  over the markers.
            for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
                const CompressedMarker& marker = orientedReadMarkers[ordinal];
                const uint32_t position = marker.position;

                // Extract the Kmer at this position.
                Kmer kmer;
                extractKmer(readSequence, position, k, kmer);

                // Find the bucket this bucket goes to.
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
            }
        }
    }
}



// This loops over buckets and creates anchors from the MarkerInfo
// objects stored in each bucket.
void Anchors::constructFromMarkerKmersThreadFunction3(uint64_t threadId )
{
    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
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



            // Each streak with the same k-mer generates a pair of anchors, but only if the
            // first OrientedReadId is on strand 0.
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

                // Only generate a pair of reverse complementary anchors
                // if the first oriented read is on strand 0.
                if(bucketInfo.front().orientedReadId.getStrand() == 0) {

                    // Generate the first anchor of the pair (no reverse complementing).
                    threadAnchors.appendVector();
                    for(const MarkerInfo& markerInfo: bucketInfo) {
                        threadAnchors.append(markerInfo);
                    }

                    // Generate the second anchor of the pair (with reverse complementing).
                    threadAnchors.appendVector();
                    for(ConstructFromMarkerKmersData::MarkerInfo markerInfo: bucketInfo) {
                        const uint32_t markerCount = uint32_t(markers[markerInfo.orientedReadId.getValue()].size());
                        markerInfo.orientedReadId.flipStrand();
                        markerInfo.ordinal = markerCount - 1 - markerInfo.ordinal;
                        threadAnchors.append(markerInfo);
                    }

                }

                // Prepare to process the next streak.
                streakBegin = streakEnd;
            }
        }
    }


}

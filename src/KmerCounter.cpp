#include "KmerCounter.hpp"
#include "deduplicate.hpp"
#include "extractKmer.hpp"
#include "Marker.hpp"
#include "MurmurHash2.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
using namespace shasta;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<KmerCounter>;



KmerCounter::KmerCounter(
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MappedMemoryOwner& mappedMemoryOwner,
    uint64_t threadCount
    ) :
    MultithreadedObject(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    reads(reads),
    markers(markers)
{
    // Initial message.
    performanceLog << timestamp << "Markers counting begins." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Get the number of reads and sanity check.
    const uint64_t readCount = reads.readCount();
    SHASTA_ASSERT(readCount == markers.size() / 2);

    const uint64_t readBatchSize = 4;
    const uint64_t bucketBatchSize = 100;


    // Create a temporary hash table to hold the marker KmerIds of all the reads.
    // (one strand only).
    const uint64_t kmerCount = markers.totalSize() / 2;
    const uint64_t averageBucketSize = 16;
    const uint64_t desiredBucketCount = kmerCount / averageBucketSize;
    const uint64_t bucketCount = std::bit_ceil(desiredBucketCount);
    hashMask = bucketCount - 1;
    cout << "Total number of marker k-mer on one strand " << kmerCount << endl;
    cout << "Number of buckets for k-mer counting " << bucketCount << " " << hex << bucketCount << endl;
    cout << "Hashing mask " << hex << hashMask << endl;
    kmerIds.createNew(largeDataName("tmp-KmerCounterKmerIds"), largeDataPageSize);

    // Gather marker KmerIds of all the reads (one strand only).
    kmerIds.beginPass1(bucketCount);
    setupLoadBalancing(readCount, readBatchSize);
    runThreads(&KmerCounter::threadFunction1, threadCount);
    kmerIds.beginPass2();
    setupLoadBalancing(readCount, readBatchSize);
    runThreads(&KmerCounter::threadFunction2, threadCount);
    kmerIds.endPass2(false, true);

    // Now we can create the final hash table to contain in each bucket
    // pairs(KmerId, frequency).
    kmerIdFrequencies.createNew(largeDataName("KmerFrequencies"), largeDataPageSize);
    kmerIdFrequencies.beginPass1(bucketCount);
    setupLoadBalancing(bucketCount, bucketBatchSize);
    runThreads(&KmerCounter::threadFunction3, threadCount);
    kmerIdFrequencies.beginPass2();
    setupLoadBalancing(bucketCount, bucketBatchSize);
    runThreads(&KmerCounter::threadFunction4, threadCount);
    kmerIdFrequencies.endPass2(false, true);

    // Remove the temporary hash table.
    kmerIds.remove();

    // Write a histogram of k-mer frequencies.
    writeHistogram();

    // Final message.
    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Marker counting completed in " << tTotal << " s." << endl;

}



void KmerCounter::threadFunction1(uint64_t /* threadId */)
{
    threadFunction12(1);
}



void KmerCounter::threadFunction2(uint64_t /* threadId */)
{
    threadFunction12(2);
}



void KmerCounter::threadFunction12(uint64_t pass)
{
    Kmer kmer;
    Kmer rcKmer;
    KmerId kmerId;
    KmerId rcKmerId;
    KmerId canonicalKmerId;

    // Loop over all batches of reads assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin,  end)) {

        // Loop over all reads assigned to this batch.
        for(uint64_t readId=begin; readId!=end; readId++) {
            const OrientedReadId orientedReadId(ReadId(readId), 0);

            // Access the sequence and markers for this read.
            const LongBaseSequenceView readSequence = reads.getRead(ReadId(readId));
            const auto readMarkers = markers[orientedReadId.getValue()];

            // Loop over the markers for this read.
            for(const CompressedMarker& compressedMarker: readMarkers) {
                const uint64_t position = compressedMarker.position;

                // Extract the Kmer at this position.
                 extractKmer(readSequence, position, k, kmer);

                 // Compute its reverse complement.
                 rcKmer = kmer.reverseComplement(k);

                 // Compute the KmerIds.
                 kmerId = kmer.id(k);
                 rcKmerId = rcKmer.id(k);
                 canonicalKmerId = min(kmerId, rcKmerId);

                 // Hash it.
                 const uint64_t hashValue = MurmurHash64A(&canonicalKmerId, sizeof(canonicalKmerId), hashSeed);
                 const uint64_t bucketId = hashValue & hashMask;

                 if(pass == 1) {
                     kmerIds.incrementCountMultithreaded(bucketId);
                 } else {
                     kmerIds.storeMultithreaded(bucketId, canonicalKmerId);
                 }
            }
        }

    }
}



void KmerCounter::threadFunction3(uint64_t /* threadId */)
{
    threadFunction34(3);
}



void KmerCounter::threadFunction4(uint64_t /* threadId */)
{
    threadFunction34(4);
}



void KmerCounter::threadFunction34(uint64_t pass)
{
    vector<KmerId> bucketKmerIds;
    vector<uint64_t> bucketFrequencies;

    // Loop over all batches of buckets assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin,  end)) {

        // Loop over all buckets assigned to this batch.
        for(uint64_t bucketId=begin; bucketId!=end; bucketId++) {
            auto inputBucket = kmerIds[bucketId];

            if(pass == 3) {
                sort(inputBucket.begin(), inputBucket.end());
            }

            // Deduplicate and count.
            bucketKmerIds.resize(inputBucket.size());
            copy(inputBucket.begin(), inputBucket.end(), bucketKmerIds.begin());
            deduplicateAndCount(bucketKmerIds, bucketFrequencies, true);
            const uint64_t n = bucketKmerIds.size();
            SHASTA_ASSERT(bucketFrequencies.size() == n);

            // Update the kmerIdFrequencies has table.
            if(pass == 3) {
                kmerIdFrequencies.incrementCount(bucketId, n);
            } else {
                auto outputBucket = kmerIdFrequencies[bucketId];
                SHASTA_ASSERT(outputBucket.size() == n);
                for(uint64_t i=0; i<n; i++) {
                    outputBucket[i] = {bucketKmerIds[i], bucketFrequencies[i]};
                }
            }
        }

    }
}


// Write a histogram of k-mer frequencies.
void KmerCounter::writeHistogram()
{

    vector<uint64_t> histogram;
    for(uint64_t bucketId=0; bucketId<kmerIdFrequencies.size(); bucketId++) {
        const auto bucket = kmerIdFrequencies[bucketId];

        for(const auto& p: bucket) {
            const uint64_t frequency = p.second;

            if(frequency >= histogram.size()) {
                histogram.resize(frequency + 1, 0);
            }

            ++histogram[frequency];
        }
    }

    ofstream csv("KmerFrequencyHistogram.csv");
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        const uint64_t frequency = histogram[coverage];
        if(frequency) {
            csv << coverage << "," << frequency << "\n";
        }
    }

}

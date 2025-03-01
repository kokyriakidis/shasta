#include "MarkerKmers.hpp"
#include "extractKmer.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
using namespace shasta;

#include "fstream.hpp"

// Expplicit instantiationn.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<MarkerKmers>;



MarkerKmers::MarkerKmers(
    uint64_t k,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t threadCount) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<MarkerKmers>(*this),
    k(k),
    reads(reads),
    markers(markers)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Figure out the number of buckets in the hash table.
    // It will contain one entry for each pairs of marker,
    // because we only store canonical k-mers.
    const uint64_t totalMarkerCount = markers.totalSize() / 2;
    const uint64_t bucketCount = (std::bit_ceil(totalMarkerCount) >> 4);
    mask = bucketCount - 1;

    // Initialize the hash table.
    markerInfos.createNew(largeDataName("MarkerKmers-MarkerInfos"), largeDataPageSize);

    // Gather markers k-mers.
    markerInfos.beginPass1(bucketCount);
    size_t batchSize = 10;
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&MarkerKmers::gatherMarkersPass1, threadCount);
    markerInfos.beginPass2();
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&MarkerKmers::gatherMarkersPass2, threadCount);
    markerInfos.endPass2(true, true);

    // Sort each bucket by Kmer.
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&MarkerKmers::sortMarkers, threadCount);

    // Create the KmerInfos.
    kmerInfos.createNew(largeDataName("MarkerKmers-KmerInfos"), largeDataPageSize);
    kmerInfos.beginPass1(bucketCount);
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&MarkerKmers::fillKmerInfosPass1, threadCount);
    kmerInfos.beginPass2();
    setupLoadBalancing(bucketCount, batchSize);
    runThreads(&MarkerKmers::fillKmerInfosPass2, threadCount);
    kmerInfos.endPass2(false, true);

    writeFrequencyHistogram();
    // writeCsv();
}



MarkerKmers::MarkerKmers(
    uint64_t k,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<MarkerKmers>(*this),
    k(k),
    reads(reads),
    markers(markers)
{
    markerInfos.accessExistingReadOnly(largeDataName("MarkerKmers-MarkerInfos"));
    kmerInfos.accessExistingReadOnly(largeDataName("MarkerKmers-KmerInfos"));
    mask = markerInfos.size() - 1;
}



void MarkerKmers::gatherMarkersPass1(uint64_t /* threadId */)
{
    gatherMarkersPass12(1);
}



void MarkerKmers::gatherMarkersPass2(uint64_t /* threadId */)
{
    gatherMarkersPass12(2);
}



void MarkerKmers::gatherMarkersPass12(uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); ++readId) {

            // Get the sequence for this read (without reverse complementing).
            const LongBaseSequenceView readSequence = reads.getRead(readId);

            // Get the two OrientedReadIds corresponding to this ReadId.
            const OrientedReadId orientedReadId0(readId, 0);
            const OrientedReadId orientedReadId1(readId, 1);

            // Get the markers on this read, without reverse complementing (strand 0).
            const auto orientedReadMarkers = markers[orientedReadId0.getValue()];
            const uint32_t orientedReadMarkerCount = uint32_t(orientedReadMarkers.size());

            // Loop  over the markers.
            for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
                const CompressedMarker& marker = orientedReadMarkers[ordinal];
                const uint32_t position = marker.position;

                // Extract the Kmer at this position and its reverse complement.
                Kmer kmer;
                extractKmer(readSequence, position, k, kmer);
                const Kmer kmerRc = kmer.reverseComplement(k);

                // Only store the lowest of the two, the canonical one.
                if(kmer <= kmerRc) {
                    const uint64_t bucketId = findBucket(kmer);
                    if(pass == 1) {
                        markerInfos.incrementCountMultithreaded(bucketId);
                    } else {
                        // This marker occurs on orientedReadId0.
                        MarkerInfo markerInfo(orientedReadId0, ordinal);
                        markerInfos.storeMultithreaded(bucketId, markerInfo);
                    }
                } else {
                    const uint64_t bucketId = findBucket(kmerRc);
                    if(pass == 1) {
                        markerInfos.incrementCountMultithreaded(bucketId);
                    } else {
                        // This marker occurs on orientedReadId1.
                        MarkerInfo markerInfo(orientedReadId1, orientedReadMarkerCount - 1 - ordinal);
                        markerInfos.storeMultithreaded(bucketId, markerInfo);
                    }
                }
            }
        }

    }
}


// Get the Kmer corresponding to a given MarkerInfo.
Kmer MarkerKmers::getKmer(const MarkerInfo& markerInfo) const
{
    // Get the OrientedReadId and marker ordinal.
    const OrientedReadId orientedReadId = markerInfo.orientedReadId;
    const uint32_t ordinal = markerInfo.ordinal;

    // Get the ReadId and Strand.
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    // Get the sequence for this read (without reverse complementing).
    const LongBaseSequenceView readSequence = reads.getRead(readId);

    if(strand == 0) {
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const CompressedMarker& marker = orientedReadMarkers[ordinal];
        const uint32_t position = marker.position;
        Kmer kmer;
        extractKmer(readSequence, position, k, kmer);
        return kmer;
    } else {

        // Find the corresponding marker on strand 0, then reverse complement it.
        const OrientedReadId orientedReadId0(readId, 0);
        const auto orientedRead0Markers = markers[orientedReadId0.getValue()];
        const uint32_t markerCount = uint32_t(orientedRead0Markers.size());
        const uint32_t ordinal0 = markerCount - 1 - ordinal;
        const CompressedMarker& marker0 = orientedRead0Markers[ordinal0];
        const uint32_t position0 = marker0.position;
        Kmer kmer0;
        extractKmer(readSequence, position0, k, kmer0);
        return kmer0.reverseComplement(k);
    }
}



// This sorts markers in buckets by their Kmer.
void MarkerKmers::sortMarkers(uint64_t /* threadId */)
{
    // Loop over all buckets assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for(uint64_t bucketId=begin; bucketId<end; bucketId++) {
            span<MarkerInfo> bucket = markerInfos[bucketId];
            sort(bucket.begin(), bucket.end(), MarkerInfoSorter(*this));
        }
    }
}



void MarkerKmers::fillKmerInfosPass1(uint64_t /* threadId */)
{
    // Loop over all buckets assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for(uint64_t bucketId=begin; bucketId<end; bucketId++) {
            span<MarkerInfo> bucket = markerInfos[bucketId];

            // The number of distinct k-mers in this bucket.
            uint64_t kmerCount = 0;

            // The bucket is sorted by Kmer.
            // Find streaks with the same Kmer.
            for(uint64_t streakBegin=0; streakBegin<bucket.size(); /* Increment later */) {
                const Kmer kmer = getKmer(bucket[streakBegin]);

                uint64_t streakEnd = streakBegin + 1;
                for(; streakEnd<bucket.size(); streakEnd++) {
                    if(getKmer(bucket[streakEnd]) != kmer) {
                        break;
                    }
                }

                ++kmerCount;

                // Prepare to process the next streak
                streakBegin = streakEnd;
            }

            // This kmerInfo bucket wil contain kmerCount entries.
            kmerInfos.incrementCountMultithreaded(bucketId, kmerCount);
        }
    }

}



void MarkerKmers::fillKmerInfosPass2(uint64_t /* threadId */)
{
    // Loop over all buckets assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for(uint64_t bucketId=begin; bucketId<end; bucketId++) {
            span<MarkerInfo> markerInfoBucket = markerInfos[bucketId];
            span<KmerInfo> kmerInfoBucket = kmerInfos[bucketId];
            uint64_t kmerInfoBucketIndex = 0;

            // The bucket is sorted by Kmer.
            // Find streaks with the same Kmer.
            for(uint64_t streakBegin=0; streakBegin<markerInfoBucket.size(); /* Increment later */) {
                const Kmer kmer = getKmer(markerInfoBucket[streakBegin]);

                uint64_t streakEnd = streakBegin + 1;
                for(; streakEnd<markerInfoBucket.size(); streakEnd++) {
                    if(getKmer(markerInfoBucket[streakEnd]) != kmer) {
                        break;
                    }
                }

                KmerInfo& kmerInfo = kmerInfoBucket[kmerInfoBucketIndex++];
                kmerInfo.markerInfo = markerInfoBucket[streakBegin];
                kmerInfo.begin = &markerInfoBucket[streakBegin] - markerInfos.begin();
                kmerInfo.end = &markerInfoBucket[streakEnd] - markerInfos.begin();

                // Prepare to process the next streak
                streakBegin = streakEnd;
            }

            SHASTA_ASSERT(kmerInfoBucketIndex == kmerInfoBucket.size());
        }
    }


}




void MarkerKmers::writeCsv() const
{
    writeFrequencyHistogram();
    writeMarkerInfosCsv2();
    writeKmerInfosCsv();
}



void MarkerKmers::writeMarkerInfosCsv1() const
{
    ofstream csv("MarkerKmers-MarkerInfos1.csv");
    csv << "Kmer,Global index,Bucket,Index in bucket,OrientedReadId,Ordinal,Position,\n";

    // Loop over all buckets.
    for(uint64_t bucketId=0; bucketId<markerInfos.size(); bucketId++) {
        const span<const MarkerInfo> bucket = markerInfos[bucketId];

        // Loop over MarkerInfos in this bucket.
        for(uint64_t i=0; i<bucket.size(); i++) {
            const MarkerInfo& markerInfo = bucket[i];

            const OrientedReadId orientedReadId = markerInfo.orientedReadId;
            const uint32_t ordinal= markerInfo.ordinal;

            const Kmer kmer = getKmer(markerInfo);

            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const CompressedMarker& marker = orientedReadMarkers[ordinal];
            const uint32_t position = marker.position;

            kmer.write(csv, k);
            csv << ",";

            csv << &markerInfo - markerInfos.begin() << ",";
            csv << bucketId << ",";
            csv << i << ",";
            csv << orientedReadId << ",";
            csv << ordinal << ",";
            csv << position << ",";

            csv << "\n";

        }
    }
}



void MarkerKmers::writeMarkerInfosCsv2() const
{
    ofstream csv("MarkerKmers-MarkerInfos2.csv");
    csv << "Kmer,Frequency,OrientedReadId,Ordinal,Position,\n";

    // Loop over all buckets.
    for(uint64_t bucketId=0; bucketId<markerInfos.size(); bucketId++) {
        const span<const KmerInfo> bucket = kmerInfos[bucketId];

        // Loop over KmerInfos in this bucket.
        for(uint64_t i=0; i<bucket.size(); i++) {
            const KmerInfo& kmerInfo = bucket[i];

            const Kmer kmer = getKmer(kmerInfo.markerInfo);

            for(uint64_t i=kmerInfo.begin; i!=kmerInfo.end; i++) {
                const MarkerInfo& markerInfo = markerInfos.begin()[i];
                const OrientedReadId orientedReadId = markerInfo.orientedReadId;
                const uint32_t ordinal = markerInfo.ordinal;

                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                const CompressedMarker& marker = orientedReadMarkers[ordinal];
                const uint32_t position = marker.position;

                kmer.write(csv, k);
                csv << ",";
                csv << kmerInfo.end - kmerInfo.begin << ",";

                csv << orientedReadId << ",";
                csv << ordinal << ",";
                csv << position << ",";

                csv << "\n";
            }
        }
    }

}


void MarkerKmers::writeKmerInfosCsv() const
{
    ofstream csv("MarkerKmers-KmerInfos.csv");
    csv << "Kmer,Frequency,GlobalIndexBegin,GlobalIndexEnd,\n";

    // Loop over all buckets.
    for(uint64_t bucketId=0; bucketId<markerInfos.size(); bucketId++) {
        const span<const KmerInfo> bucket = kmerInfos[bucketId];

        // Loop over KmerInfos in this bucket.
        for(uint64_t i=0; i<bucket.size(); i++) {
            const KmerInfo& kmerInfo = bucket[i];

            const Kmer kmer = getKmer(kmerInfo.markerInfo);

            kmer.write(csv, k);
            csv << ",";

            csv << kmerInfo.end - kmerInfo.begin << ",";
            csv << kmerInfo.begin << ",";
            csv << kmerInfo.end << ",";

            csv << "\n";
        }
    }

}



void MarkerKmers::writeFrequencyHistogram() const
{
    // Create the histogram.
    vector<uint64_t> histogram;
    for(uint64_t bucketId=0; bucketId<markerInfos.size(); bucketId++) {
        const span<const KmerInfo> bucket = kmerInfos[bucketId];

        for(uint64_t i=0; i<bucket.size(); i++) {
            const KmerInfo& kmerInfo = bucket[i];
            const uint64_t coverage = kmerInfo.end - kmerInfo.begin;

            if(histogram.size() <= coverage) {
                histogram.resize(coverage + 1, 0);
            }
            ++histogram[coverage];
        }
    }

    // Write it out.
    ofstream csv("MarkerKmers-Histogram.csv");
    csv << "Coverage,Frequency,\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        const uint64_t frequency = histogram[coverage];
        if(frequency) {
            csv << coverage << "," << frequency << ",\n";
        }
    }

}



uint64_t MarkerKmers::getFrequency(const Kmer& kmer) const
{
    const Kmer kmerRc = kmer.reverseComplement(k);
    if(kmer <= kmerRc) {
        return getMarkerInfos(kmer).size();
    } else {
        return getMarkerInfos(kmerRc).size();
    }
}



// Return a span of the MarkerInfos for a given canonical k-mer.
// If this is called for a non-canonical k-mer, it returns an empty span.
span<const MarkerKmers::MarkerInfo> MarkerKmers::getMarkerInfos(const Kmer& kmer) const
{
    if(not kmer.isCanonical(k)) {
        return span<const MarkerInfo>();
    }

    const uint64_t bucketId = findBucket(kmer);
    const span<const KmerInfo> bucket = kmerInfos[bucketId];
    for(const KmerInfo& kmerInfo: bucket) {
        if(getKmer(kmerInfo.markerInfo) == kmer) {
            const MarkerInfo* const begin = markerInfos.begin() + kmerInfo.begin;
            const MarkerInfo* const end   = markerInfos.begin() + kmerInfo.end;
            return span<const MarkerInfo>(begin, end);
        }
    }

    // We did not find this Kmer in the bucket where it would have been.
    return span<const MarkerInfo>();
}



// Get MarkerInfo objects for a given Kmer.
void MarkerKmers::get(
    const Kmer& kmer,
    vector<MarkerInfo>& v) const
{
    v.clear();
    const Kmer kmerRc = kmer.reverseComplement(k);

    if(kmer <= kmerRc) {

        // Kmer is canonical.
        span<const MarkerInfo> s = getMarkerInfos(kmer);
        for(const MarkerInfo& markerInfo: s) {
            v.push_back(markerInfo);
        }
    } else {

        // Kmer is not canonical but kmerRc is.
        span<const MarkerInfo> s = getMarkerInfos(kmerRc);
        for(const MarkerInfo& markerInfo: s) {
            v.push_back(reverseComplement(markerInfo));
        }
        sort(v.begin(), v.end(), MarkerInfoSorter(*this));
    }
}



// Construct the reverse complement of a MarkerInfo.
MarkerKmers::MarkerInfo MarkerKmers::reverseComplement(const MarkerInfo& markerInfo) const
{
    MarkerInfo markerInfoRc = markerInfo;
    const uint32_t orientedReadMarkerCount = uint32_t(markers[markerInfo.orientedReadId.getValue()].size());

    markerInfoRc.orientedReadId.flipStrand();
    markerInfoRc.ordinal = orientedReadMarkerCount - 1 - markerInfoRc.ordinal;

    return markerInfoRc;
}

#pragma once

#include "MemoryMappedVectorOfVectors.hpp"
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "MurmurHash2.hpp"
#include "ReadId.hpp"

#include "tuple.hpp"



// The MarkerKmers object keeps track of the locations in the oriented reads
// where each marker k-mer appears. It is only used for alignment-free assembly
// (--Assembly.mode 3 --Assembly.mode3.anchorCreationMethod FromMarkerKmers)
// but can also be created using CreateMarkerKmers.py.

namespace shasta {
    class MarkerKmers;

    class CompressedMarker;
    class Reads;
}



class shasta::MarkerKmers :
    public MappedMemoryOwner,
    public MultithreadedObject<MarkerKmers> {
public:

    MarkerKmers(
        uint64_t k,
        const MappedMemoryOwner&,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        uint64_t threadCount);

    MarkerKmers(
        uint64_t k,
        const MappedMemoryOwner&,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    bool isOpen() const
    {
        return markerInfos.isOpen() and kmerInfos.isOpen();
    }

    uint64_t getFrequency(const Kmer&) const;

    // Class to describe a single marker of an oriented read.
    class MarkerInfo {
    public:
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        MarkerInfo() {}
        MarkerInfo(OrientedReadId orientedReadId, uint32_t ordinal) :
            orientedReadId(orientedReadId), ordinal(ordinal) {}
    };

    // Get MarkerInfo objects for a given Kmer.
    void get(const Kmer&, vector<MarkerInfo>&) const;

    // Return the number of Kmers stored.
    uint64_t size() const
    {
        return kmerInfos.totalSize();
    }

    // Get MarkerInfo objects for the i-th k-mer stored.
    span<const MarkerInfo> operator[](uint64_t i) const
    {
        const KmerInfo& kmerInfo = kmerInfos.begin()[i];
        return span<const MarkerInfo>(
            markerInfos.begin() + kmerInfo.begin,
            markerInfos.begin() + kmerInfo.end);
    }

    // Construct the reverse complement of a MarkerInfo.
    MarkerInfo reverseComplement(const MarkerInfo&) const;

private:

    // Constructor arguments.
    uint64_t k;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

    // Get the Kmer corresponding to a given MarkerInfo.
    Kmer getKmer(const MarkerInfo&) const;



    // A function object class that sorts MarkerInfo objects by Kmer.
    class MarkerInfoSorter {
    public:
        MarkerInfoSorter(const MarkerKmers& markerKmers) : markerKmers(markerKmers) {}
        bool operator()(const MarkerInfo& x, const MarkerInfo& y) const
        {
            const Kmer xKmer = markerKmers.getKmer(x);
            const Kmer yKmer = markerKmers.getKmer(y);
            return
                tie(xKmer, x.orientedReadId, x.ordinal) <
                tie(yKmer, y.orientedReadId, y.ordinal);
        }
    private:
        const MarkerKmers& markerKmers;
    };

    // A hash table that will contain a MarkerInfo object
    // for each marker in all oriented reads.
    // Indexed by bucketId.
    // The bucket is computed by hashing the k-mer of each marker,
    // so all markers with the same k-mer end up in the same bucket.
    // Only canonical k-mers are stored (the ones for which
    // Kmers::isCanonical(k) returns true).
    // After construction, the markerInfos in each bucket
    // are sorted by Kmer, so the ones for each Kmer are
    // all together.
    MemoryMapped::VectorOfVectors<MarkerInfo, uint64_t> markerInfos;

    void gatherMarkersPass1(uint64_t threadId);
    void gatherMarkersPass2(uint64_t threadId);
    void gatherMarkersPass12(uint64_t pass);
    void sortMarkers(uint64_t threadId);

    // The number of buckets is chosen equal to a power of 2,
    // so the bucketId can be obtained with a simple bitwise and
    // with a mask equal to the number of buckets minus 1.
    uint64_t mask;

    // Find the bucket that a canonical k-mer goes to.
    uint64_t findBucket(const Kmer& kmer) const
    {
        const uint64_t hashValue = MurmurHash64A(&kmer, sizeof(kmer), 1241);
        return hashValue & mask;
    }



    // After sortMarkers is called, all the MarkerInfos for each marker Kmer
    // are contiguous and in a single bucket.
    // Index them here.
    // KmerInfos are organized in a hash table similarly to MarkerInfos.
    class KmerInfo {
    public:
        // A MarkerInfo object used to identify the Kmer.
        // We do this to avoid storing the Kmer.
        MarkerInfo markerInfo;

        // Begin/end for this Kmer, relative to markerInfos.begin().
        uint64_t begin;
        uint64_t end;
    };
    MemoryMapped::VectorOfVectors<KmerInfo, uint64_t> kmerInfos;
    void fillKmerInfosPass1(uint64_t threadId);
    void fillKmerInfosPass2(uint64_t threadId);

    // Return a span of the MarkerInfos for a given canonical k-mer.
    // If this is called for a non-canonical k-mer, it returns an empty span.
    span<const MarkerInfo> getMarkerInfos(const Kmer&) const;


    void writeCsv() const;
    void writeFrequencyHistogram() const;
    void writeMarkerInfosCsv1() const;
    void writeMarkerInfosCsv2() const;
    void writeKmerInfosCsv() const;

};

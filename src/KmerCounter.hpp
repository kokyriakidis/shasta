#pragma once

// KmerCounter is a class that holds a hash table capable of returning
// the number of time a marker KmerId appears in
// the markers of this assembly.

// Shasta.
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"



namespace shasta {
    class KmerCounter;

    class CompressedMarker;
    class KmerDistributionInfo;
    class Reads;
}


class shasta::KmerCounter :
    public MultithreadedObject<KmerCounter>,
    public MappedMemoryOwner {
public:

    // This constructor creates the KmerIdFrequencies hash table.
    KmerCounter(
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MappedMemoryOwner& mappedMemoryOwner,
        uint64_t threadCount
        );

    // This constructor accesses an existing KmerIdFrequencies hash table.
    KmerCounter(
        uint64_t k,
        const MappedMemoryOwner& mappedMemoryOwner
        );

    bool isAvailable() const
    {
        return kmerIdFrequencies.isOpen();
    }

    uint64_t getFrequency(KmerId) const;
    uint64_t getFrequency(const Kmer&) const;

    // Override the frequencies stored in this KmerCounter
    // with the ones obtained from another KmerCounter.
    void overrideFrequencies(const KmerCounter&);

private:

    // Data passed in to the constructor.
    uint64_t k;
    Reads const* readsPointer = 0;
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t> const* markersPointer = 0;

    // Hashing.
    const uint32_t hashSeed = 12771;
    uint64_t hashMask;

    // A temporary hash table that stores KmerIds of each marker k-mer we encounter.
    // KmerIds are stored in canonical form. Canonical KmerId is the
    // lowest of the KmerId and the KmerId of its reverse complement.
    MemoryMapped::VectorOfVectors<KmerId, uint64_t> kmerIds;

    // This is the persistent hash table that contains in each bucket
    // pairs(KmerId, frequency).
    MemoryMapped::VectorOfVectors<pair<KmerId, uint64_t>, uint64_t> kmerIdFrequencies;


    // Passes 1 and 2 gather marker KmerIds in the kmerIds hash table.
    void threadFunction1(uint64_t threadId);
    void threadFunction2(uint64_t threadId);
    void threadFunction12(uint64_t pass);

    // Passes 3 and 4 gather fill in the KmerIdFrequencies hash table.
    void threadFunction3(uint64_t threadId);
    void threadFunction4(uint64_t threadId);
    void threadFunction34(uint64_t pass);

public:

    // The frequency histogram.
    // It consists of pairs (coverage, frequency).
    void createHistogram();
    void writeHistogram(ostream&) const;
    MemoryMapped::Vector< pair<uint64_t, uint64_t> > histogram;

    void getHistogramInfo(KmerDistributionInfo&) const;
};

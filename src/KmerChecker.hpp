#ifndef SHASTA_KMER_CHECKER_HPP
#define SHASTA_KMER_CHECKER_HPP

// As we transition to longer markers, we will no longer be able to store a k-mer
// table. Instead, the KmerChecker will be used to find out if a given k-mer
// is a marker. The initial implementation of the KmerChecker is table based,
// but later we will switch to hashing.

#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"
#include "MultithreadedObject.hpp"

#include "utility.hpp"

namespace shasta {
    class KmerChecker;
    class KmersOptions;
    class Reads;
}

class shasta::KmerChecker :
    public MultithreadedObject<KmerChecker>,
    public MappedMemoryOwner {
public:

    bool isMarker(KmerId kmerId) const
    {
        return kmerTable[kmerId].isMarker;
    }

    // Initial creation.
    KmerChecker(
        const KmersOptions&,
        uint64_t threadCount,
        const Reads&,
        const MappedMemoryOwner&);

    // Create from binary data.
    KmerChecker(
        uint64_t k,
        uint64_t generationMethod,
        const Reads&,
        const MappedMemoryOwner&);

private:
    uint64_t k;
    const Reads& reads;
    MemoryMapped::Vector<KmerInfo> kmerTable;
    void initializeKmerTable();

    // Creation functions for each generation method.

    // Randomly select k-mers.
    void create0(
        double probability, // The probability that a k-mer is selected as a marker.
        int seed            // For random number generator.
    );

    // Randomly select the k-mers to be used as markers, but
    // excluding those that are globally overenriched in the input reads,
    // as measured by total frequency in all reads.
    void create1(
        double probability,
        int seed,
        double enrichmentThreshold,
        const Reads&,
        uint64_t threadCount);

    // Randomly select the k-mers to be used as markers, but
    // excluding those that are overenriched even in a single oriented read.
    void create2(
        double probability,
        int seed,
        double enrichmentThreshold,
        uint64_t threadCount);

    // Read the k-mers to be used as markers from a file.
    void create3(const string& fileName);

    // Randomly select the k-mers to be used as markers, but
    // excluding those that appear in two copies close to each other
    // even in a single oriented read.
    void create4(
        double probability,
        int seed,
        uint64_t distanceThreshold,
        uint64_t threadCount);

    void computeKmerFrequency(size_t threadId);

    class SelectKmers2Data {
    public:

        double enrichmentThreshold;

        // The number of times each k-mer appears in an oriented read.
        // Indexed by KmerId.
        MemoryMapped::Vector<uint64_t> globalFrequency;

        // The number of oriented reads that each k-mer is
        // over-enriched in by more than a factor enrichmentThreshold.
        // Indexed by KmerId.
        MemoryMapped::Vector<ReadId> overenrichedReadCount;

    };
    SelectKmers2Data selectKmers2Data;
    void selectKmers2ThreadFunction(size_t threadId);


    void selectKmers4ThreadFunction(size_t threadId);
    class SelectKmers4Data {
    public:

        // The number of times each k-mer appears in an oriented read.
        // Indexed by KmerId.
        MemoryMapped::Vector<uint64_t> globalFrequency;

        // The minimum distance at which two copies of each k-mer
        // appear in any oriented read.
        // Indexed by KmerId.
        MemoryMapped::Vector< pair<std::mutex, uint32_t> > minimumDistance;

    };
    SelectKmers4Data selectKmers4Data;

};


#endif

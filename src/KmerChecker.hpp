#ifndef SHASTA_KMER_CHECKER_HPP
#define SHASTA_KMER_CHECKER_HPP

// Shasta.
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"
#include "MultithreadedObject.hpp"

// Standard library.
#include "memory.hpp"
#include "utility.hpp"

namespace shasta {

    class KmerChecker;
    class KmerCheckerFactory;
    class HashedKmerChecker;

    class KmerTable;
    class KmerTable0;
    class KmerTable1;
    class KmerTable2;
    class KmerTable3;
    class KmerTable4;

    class KmersOptions;
    class Reads;
}



// The KmerChecker is an abstract class that knows how to find
// out if a k-mer is a marker.
class shasta::KmerChecker {
public:
    virtual bool isMarker(KmerId kmerId) const = 0;
};



// The KmerCheckerFactory knows how to create the appropriate
// type of KmerChecker for the options used.
class shasta::KmerCheckerFactory {
public:

    static shared_ptr<KmerChecker> createNew(
        const KmersOptions&,
        uint64_t threadCount,
        const Reads&,
        const MappedMemoryOwner&);

    static shared_ptr<KmerChecker> createFromBinaryData(
        uint64_t k,
        uint64_t generationMethod,
        const Reads&,
        const MappedMemoryOwner&);
};



// Old implementations are table based.
// There are derived classes to support all 5 marker generation methods
// but they are limited to k-mer lengths k<16.
class shasta::KmerTable :
    public KmerChecker,
    public MappedMemoryOwner {
public:

    bool isMarker(KmerId kmerId) const
    {
        return kmerTable[kmerId].isMarker;
    }

    KmerTable(uint64_t k, bool createNew, const MappedMemoryOwner&);

protected:

    uint64_t k;
    MemoryMapped::Vector<KmerInfo> kmerTable;

private:
    void createKmerTable();
    void accessKmerTable();

};



// Marker k-mer generation method 0 (used when --Kmers.generationMethod 0).
class shasta::KmerTable0 : public KmerTable {
public:

    // Construct from scratch.
    KmerTable0(
        uint64_t k,
        double probability,
        int seed,
        const MappedMemoryOwner&);

    // Construct from binary data.
    KmerTable0(
        uint64_t k,
        const MappedMemoryOwner&);

};



// Marker k-mer generation method 1 (used when --Kmers.generationMethod 1).
class shasta::KmerTable1 :
    public KmerTable,
    public MultithreadedObject<KmerTable1> {
public:

    // Construct from scratch.
    KmerTable1(
        uint64_t k,
        double probability,
        int seed,
        double enrichmentThreshold,
        const Reads&,
        uint64_t threadCount,
        const MappedMemoryOwner&);

    // Construct from binary data.
    KmerTable1(
        uint64_t k,
        const Reads&,
        const MappedMemoryOwner&);

private:
    const Reads& reads;
    void computeKmerFrequency(size_t threadId);
};



// Marker k-mer generation method 2 (used when --Kmers.generationMethod 2).
class shasta::KmerTable2 :
    public KmerTable,
    public MultithreadedObject<KmerTable2> {
public:

    // Construct from scratch.
    KmerTable2(
        uint64_t k,
        double probability,
        int seed,
        double enrichmentThreshold,
        const Reads&,
        uint64_t threadCount,
        const MappedMemoryOwner&);

    // Construct from binary data.
    KmerTable2(
        uint64_t k,
        const Reads&,
        const MappedMemoryOwner&);

private:
    const Reads& reads;
    double enrichmentThreshold;

    // The number of times each k-mer appears in an oriented read.
    // Indexed by KmerId.
    MemoryMapped::Vector<uint64_t> globalFrequency;

    // The number of oriented reads that each k-mer is
    // over-enriched in by more than a factor enrichmentThreshold.
    // Indexed by KmerId.
    MemoryMapped::Vector<ReadId> overenrichedReadCount;

    void threadFunction(size_t threadId);
};



// Marker k-mer generation method 3 (used when --Kmers.generationMethod 3).
class shasta::KmerTable3: public KmerTable {
public:

    // Construct from scratch.
    KmerTable3(
        uint64_t k,
        uint64_t readRepresentation,
        const string& fileName,
        const MappedMemoryOwner&);

    // Construct from binary data.
    KmerTable3(
        uint64_t k,
        const MappedMemoryOwner&);

};



// Marker k-mer generation method 4 (used when --Kmers.generationMethod 4).
class shasta::KmerTable4 :
    public KmerTable,
    public MultithreadedObject<KmerTable4> {
public:

    // Construct from scratch.
    KmerTable4(
        uint64_t k,
        double probability,
        int seed,
        uint64_t distanceThreshold,
        const Reads&,
        uint64_t threadCount,
        const MappedMemoryOwner&);

    // Construct from binary data.
    KmerTable4(
        uint64_t k,
        const Reads&,
        const MappedMemoryOwner&);

public:
    const Reads& reads;

    void threadFunction(size_t threadId);

    // The number of times each k-mer appears in an oriented read.
    // Indexed by KmerId.
    MemoryMapped::Vector<uint64_t> globalFrequency;

    // The minimum distance at which two copies of each k-mer
    // appear in any oriented read.
    // Indexed by KmerId.
    MemoryMapped::Vector< pair<std::mutex, uint32_t> > minimumDistance;
};



// The new implementation of the KmerChecker is not table based
// and uses hashing instead.
// It only supports marker generation method 0 (random generation)
// but allow marker lengths k<32.
// This class is not yet implemented.
class shasta::HashedKmerChecker : public KmerChecker {
public:
    bool isMarker(KmerId kmerId) const;
};

#endif

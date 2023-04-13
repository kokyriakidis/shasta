#include "KmerChecker.hpp"
#include "AssemblerOptions.hpp"
using namespace shasta;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<KmerChecker>;



// Initial creation.
KmerChecker::KmerChecker(
    const KmersOptions& kmersOptions,
    uint64_t threadCount,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<KmerChecker>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(kmersOptions.k),
    reads(reads)
{
    // Sanity check on the value of k.
    if(k > Kmer::capacity) {
        throw runtime_error("K-mer capacity exceeded.");
    }

    initializeKmerTable();
    const int seed = 231;

    switch(kmersOptions.generationMethod) {
    case 0:
        // Randomly select kmers.
        create0(
            kmersOptions.probability,
            seed);
        break;

    case 1:
        // Randomly select the k-mers to be used as markers, but
        // excluding those that are globally overenriched in the input reads,
        // as measured by total frequency in all reads.
        create1(
            kmersOptions.probability,
            seed,
            kmersOptions.enrichmentThreshold,
            reads,
            threadCount);
        break;

    case 2:
        // Randomly select the k-mers to be used as markers, but
        // excluding those that are overenriched even in a single oriented read.
        create2(
            kmersOptions.probability,
            seed,
            kmersOptions.enrichmentThreshold,
            threadCount);
        break;

    case 3:
        // Read the k-mers to be used as markers from a file.
        create3(
            kmersOptions.file);
        break;

    case 4:
        // Randomly select the k-mers to be used as markers, but
        // excluding those that appear in two copies close to each other
        // even in a single oriented read.
        create4(
            kmersOptions.probability,
            seed,
            kmersOptions.distanceThreshold,
            threadCount);
        break;

    default:
        throw runtime_error("Invalid --Kmers generationMethod. "
            "Specify a value between 0 and 4, inclusive.");
    }
}



// Create from binary data.
KmerChecker::KmerChecker(
    uint64_t k,
    uint64_t generationMethod,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<KmerChecker>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    reads(reads)
{
    kmerTable.accessExistingReadOnly(largeDataName("Kmers"));
}

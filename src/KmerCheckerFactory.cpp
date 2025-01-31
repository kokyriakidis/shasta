#include "KmerCheckerFactory.hpp"
#include "KmerCheckerFromFile.hpp"
#include "Kmer.hpp"
#include "KmerTable.hpp"
#include "HashedKmerChecker.hpp"
#include "AssemblerOptions.hpp"
#include "Reads.hpp"
using namespace shasta;



std::shared_ptr<KmerChecker> KmerCheckerFactory::createNew(
    const KmersOptions& kmersOptions,
    uint64_t threadCount,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner)
{
    // For generation method 0, always use the HashedKmerChecker.
    if(kmersOptions.generationMethod == 0) {
        return make_shared<HashedKmerChecker>(
            kmersOptions.k,
            kmersOptions.probability,
            mappedMemoryOwner);
    }

    // For generation method 3, always use the KmerCheckerFromFile.
    if(kmersOptions.generationMethod == 3) {
        return make_shared<KmerCheckerFromFile>(
            kmersOptions.k,
            kmersOptions.file,
            mappedMemoryOwner);
    }

    // In all other cases, we are limited to k<=16.
    if(kmersOptions.k > int(Kmer16::capacity)) {
        throw runtime_error("Kmer generation method " +
            to_string(kmersOptions.generationMethod) +
            " is only supported for a maximum marker length of 15.");
    }

    const int seed = 231;
    switch(kmersOptions.generationMethod) {
        case 0:
        return make_shared<KmerTable0>(
            kmersOptions.k,
             kmersOptions.probability,
             seed,
             mappedMemoryOwner);

        case 1:
        return make_shared<KmerTable1>(
             kmersOptions.k,
             kmersOptions.probability,
             seed,
             kmersOptions.enrichmentThreshold,
             reads,
             threadCount,
             mappedMemoryOwner);

        case 2:
        return make_shared<KmerTable2>(
             kmersOptions.k,
             kmersOptions.probability,
             seed,
             kmersOptions.enrichmentThreshold,
             reads,
             threadCount,
             mappedMemoryOwner);

        case 3:
            // This case was handled above.
            SHASTA_ASSERT(0);

        case 4:
        return make_shared<KmerTable4>(
            kmersOptions.k,
            kmersOptions.probability,
            seed,
            kmersOptions.distanceThreshold,
            reads,
            threadCount,
            mappedMemoryOwner);

        default:
            throw runtime_error("Invalid --Kmers generationMethod. "
                "Specify a value between 0 and 4, inclusive.");
     }
}



std::shared_ptr<shasta::KmerChecker> KmerCheckerFactory::createFromBinaryData(
    uint64_t k,
    uint64_t generationMethod,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner)
{
    // For generation method 0, always use the HashedKmerChecker.
    if(generationMethod == 0) {
        return make_shared<HashedKmerChecker>(mappedMemoryOwner);
    }

    // For generation method 0, always use the KmerCheckerFromFile.
    if(generationMethod == 3) {
        return make_shared<KmerCheckerFromFile>(k, mappedMemoryOwner);
    }

    switch(generationMethod) {
    case 0:
         return make_shared<KmerTable0>(k, mappedMemoryOwner);

    case 1:
         return make_shared<KmerTable1>(k, reads, mappedMemoryOwner);

    case 2:
         return make_shared<KmerTable2>(k, reads, mappedMemoryOwner);

    case 3:
        // This case was handled above.
        SHASTA_ASSERT(0);

    case 4:
         return make_shared<KmerTable4>(k, reads, mappedMemoryOwner);


    default:
        throw runtime_error("Invalid --Kmers generationMethod. "
            "Specify a value between 0 and 4, inclusive.");
    }
}

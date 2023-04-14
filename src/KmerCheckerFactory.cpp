#include "KmerCheckerFactory.hpp"
#include "KmerTable.hpp"
#include "AssemblerOptions.hpp"
#include "Reads.hpp"
using namespace shasta;



std::shared_ptr<KmerChecker> KmerCheckerFactory::createNew(
    const KmersOptions& kmersOptions,
    uint64_t threadCount,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner)
{
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
        return make_shared<KmerTable3>(
            kmersOptions.k,
            reads.representation,
            kmersOptions.file,
            mappedMemoryOwner);

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
    switch(generationMethod) {
    case 0:
         return make_shared<KmerTable0>(k, mappedMemoryOwner);

    case 1:
         return make_shared<KmerTable1>(k, reads, mappedMemoryOwner);

    case 2:
         return make_shared<KmerTable2>(k, reads, mappedMemoryOwner);

    case 3:
         return make_shared<KmerTable3>(k, mappedMemoryOwner);

    case 4:
         return make_shared<KmerTable4>(k, reads, mappedMemoryOwner);


    default:
        throw runtime_error("Invalid --Kmers generationMethod. "
            "Specify a value between 0 and 4, inclusive.");
    }
}

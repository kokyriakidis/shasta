#ifndef SHASTA_HASHED_KMER_CHECKER_HPP
#define SHASTA_HASHED_KMER_CHECKER_HPP

#include "KmerChecker.hpp"
#include "MappedMemoryOwner.hpp"

namespace shasta {
    class HashedKmerChecker;
}


// The new implementation of the KmerChecker is not table based
// and uses hashing instead.
// It only supports marker generation method 0 (random generation)
// but allow marker lengths k<32.
class shasta::HashedKmerChecker :
    public KmerChecker,
    public MappedMemoryOwner {
public:
    bool isMarker(KmerId) const;

    // Initial creation.
    HashedKmerChecker(uint64_t k, double markerDensity, const MappedMemoryOwner&);

    // Creation from binary data.
    HashedKmerChecker(const MappedMemoryOwner&);

private:
    uint64_t k;
    uint32_t hashThreshold;

    // This is used to store the hashThreshold in binary data.
    class HashedKmerCheckerData {
    public:
        uint64_t k;
        uint32_t hashThreshold;
    };
};



#endif


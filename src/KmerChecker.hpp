#ifndef SHASTA_KMER_CHECKER_HPP
#define SHASTA_KMER_CHECKER_HPP

// Shasta.
#include "shastaTypes.hpp"

namespace shasta {
    class KmerChecker;
    class HashedKmerChecker;
}



// The KmerChecker is an abstract class that knows how to find
// out if a k-mer is a marker.
class shasta::KmerChecker {
public:
    virtual bool isMarker(KmerId) const = 0;
};



// The new implementation of the KmerChecker is not table based
// and uses hashing instead.
// It only supports marker generation method 0 (random generation)
// but allow marker lengths k<32.
// This class is not yet implemented.
class shasta::HashedKmerChecker : public KmerChecker {
public:
    bool isMarker(KmerId) const;
};

#endif

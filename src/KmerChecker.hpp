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
// All implementations must guarantee that if a KmerId if a marker
// its reverse complement is also a marker.
class shasta::KmerChecker {
public:
    virtual bool isMarker(KmerId) const = 0;
};

#endif

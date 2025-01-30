#pragma once

#include "KmerChecker.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"
#include "shastaTypes.hpp"

namespace shasta {
    class KmerCheckerFromFile;
}


// A KmerChecker that reads the marker k-mers from an input file.
class shasta::KmerCheckerFromFile :
    public KmerChecker,
    public MappedMemoryOwner {
public:
    bool isMarker(KmerId) const;

    // Initial creation.
    KmerCheckerFromFile(uint64_t k, const string& fileName, const MappedMemoryOwner&);

    // Creation from binary data.
    KmerCheckerFromFile(uint64_t k, const MappedMemoryOwner&);

private:
    uint64_t k;

    MemoryMapped::Vector<KmerId> markerKmers;

};

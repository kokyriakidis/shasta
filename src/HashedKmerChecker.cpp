// Shasta.
#include "HashedKmerChecker.hpp"
#include "Kmer.hpp"
#include "MemoryMappedObject.hpp"
using namespace shasta;

// MurmurHash.
#include "MurmurHash2.hpp"

// Standard library.
#include <cmath>



// We must guarantee that if a KmerId if a marker
// its reverse complement is also a marker.
// To do this we check both.
// This will usually require two calls to MurmurHash2,
// but this is probably still faster than two cache misses
// in the old k-mer table.
bool HashedKmerChecker::isMarker(KmerId kmerId) const
{
    // Check the KmerId.
    if(MurmurHash2(&kmerId, sizeof(kmerId), 267457831) < hashThreshold) {
        return true;
    }

    // Check its reverse complement.
    const Kmer kmer(kmerId, k);
    const Kmer kmerRc = kmer.reverseComplement(k);
    const KmerId kmerIdRc = KmerId(kmerRc.id(k));
    return MurmurHash2(&kmerIdRc, sizeof(kmerId), 267457831) < hashThreshold;
}



// Initial creation.
HashedKmerChecker::HashedKmerChecker(
    uint64_t k,
    double markerDensity,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    k(k)
{
    // Sanity check on the marker density.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }



    // Compute the hash threshold that achieves the required marker density.

    // In this computation, we neglect self-complementary k-mers,
    // which are a small minority of total.

    // Call:
    // - hashMax the maximum possible value of a hash
    // - hashValue the hash value for a given KmerId
    // - hashValueRc the hash value for its reverse complement (for length k)
    // - p = hashThreshold / hashMax

    // A KmerId is a marker if
    // Event A: hashValue < hashThreshold
    // OR
    // Event B: hashValueRc < hashThreshold
    // Event A occurs with probability P(A) = hashThreshold / hashMax = p.
    // Event B also occurs with probability P(B) = hashThreshold / hashMax = p.

    // If we use a good hash function, we can consider A and B uncorrelated.
    // Therefore we can use the standard formula:
    // P(A or B) = 1 - P(not(A or B)) =
    // 1 - P((not A) and (not B)) =
    // 1 - (P(not A)) * P(not B)) =
    // 1 - (1 - P(A)) * (1 - P(B))
    // It can also be verified by simple algebra that this is equal to the standard formula
    // P(A or B) =
    // P(A) + P(B) - P(A and B) =
    // P(A) + P(B) - P(A) * P(B)
    // but we don't need this part.

    // Using the above we get:
    // markerDensity =
    // P(A or B) =
    // 1 - (1 - P(A)) * (1 - P(B)) ==
    // 1 - (1 - p)^2
    // (Because P(A) = P(B) = p).
    // From
    // markerDensity = 1 - (1 - p)^2
    // we get
    // p = 1 - sqrt(1 - markerDensity)
    // And finally hashThreshold = hashMax * p.

    const double p = 1. - std::sqrt(1. - markerDensity);
    const double hashMax = std::numeric_limits<uint32_t> :: max();
    hashThreshold = uint32_t(std::round(double(hashMax) * p));

    // Store k and the hash threshold in binary data.
    MemoryMapped::Object<HashedKmerCheckerData> data;
    data.createNew(largeDataName("HashedKmerChecker"), largeDataPageSize);
    data->k = k;
    data->hashThreshold = hashThreshold;

}



// Creation from binary data.
HashedKmerChecker::HashedKmerChecker(
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    MemoryMapped::Object<HashedKmerCheckerData> data;
    data.accessExistingReadOnly(largeDataName("HashedKmerChecker"));
    k = data->k;
    hashThreshold = data->hashThreshold;
}

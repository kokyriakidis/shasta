#ifndef SHASTA_KMER_HPP
#define SHASTA_KMER_HPP

#include "shastaTypes.hpp"
#include "ShortBaseSequence.hpp"
#include <limits>

namespace shasta {

    // Type used to represent a k-mer.
    // This limits the maximum k-mer length that can be used.
    // If this changes, KmerId must also be changed.
    using Kmer = ShortBaseSequence16;
    static_assert(
        std::numeric_limits<KmerId>::digits == 2*Kmer::capacity,
        "Kmer and KmerId types are inconsistent.");
}

#endif

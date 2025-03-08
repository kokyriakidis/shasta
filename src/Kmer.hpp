#ifndef SHASTA_KMER_HPP
#define SHASTA_KMER_HPP

#include "BitCounter.hpp"
#include "shastaTypes.hpp"
#include "ShortBaseSequence.hpp"

namespace shasta {

    // Type used to represent a k-mer.
    // This limits the maximum k-mer length that can be used.
    // If this changes, KmerId must also be changed.
    using Kmer16 = ShortBaseSequence16;
    using Kmer32 = ShortBaseSequence32;
    using Kmer64 = ShortBaseSequence64;

#ifdef SHASTA_LONG_MARKERS
    using Kmer128 = ShortBaseSequence128;
    using Kmer = Kmer128;
    static_assert(
        BitCounter<KmerId>::numberOfBits == 256,
        "Unexpected number of bits in KmerId.");
    static_assert(
        Kmer::capacity == 128,
        "Unexpected Kmer capacity.");
#else
    using Kmer = Kmer64;
    static_assert(
        BitCounter<KmerId>::numberOfBits == 128,
        "Unexpected number of bits in KmerId.");
    static_assert(
        Kmer::capacity == 64,
        "Unexpected Kmer capacity.");
#endif

    static_assert(
        BitCounter<KmerId>::numberOfBits == 2 * Kmer::capacity,
        "Kmer and KmerId types are inconsistent.");
}

#endif

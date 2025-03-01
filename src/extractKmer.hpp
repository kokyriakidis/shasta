#ifndef SHASTA_EXTRACT_KMER_HPP
#define SHASTA_EXTRACT_KMER_HPP

// Extract the kmer at a given position in a LongBaseSequenceView.
// Uses bit operations for speed, so it can be used in
// performance critical code.

#include "SHASTA_ASSERT.hpp"

#include <concepts>
#include "cstdint.hpp"



namespace shasta {
    class LongBaseSequenceView;
    template<class Int> class ShortBaseSequence;

    template<class Int> void extractKmer(
        const LongBaseSequenceView&,
        uint64_t position,
        uint64_t length,
        ShortBaseSequence<Int>&);

    // Extract n bits from x, starting at position xPosition,
    // and store them in y, starting at position yPosition,
    // leaving the remaining bits of y unchanged.
    // Here, xPosition and yPosition are counted with
    // 0 at the most significant bit and moving towards the least
    // significant bit.
    // The above is repeated for x[0]->y[0], x[1]->y[1].
    template<class Int> inline void extractBits(
        const uint64_t* x,
        uint64_t xPosition, // 0 = MSB
        uint64_t n,
        Int* y,
        uint64_t yPosition  // 0 = MSB
        );

    // Explicit instantiations are in extractKmer.cpp.
    extern template void extractKmer(
        const LongBaseSequenceView&,
        uint64_t position,
        uint64_t length,
        ShortBaseSequence<uint8_t>&);

    extern template void extractKmer(
        const LongBaseSequenceView&,
        uint64_t position,
        uint64_t length,
        ShortBaseSequence<uint16_t>&);

    extern template void extractKmer(
        const LongBaseSequenceView&,
        uint64_t position,
        uint64_t length,
        ShortBaseSequence<uint32_t>&);

    extern template void extractKmer(
        const LongBaseSequenceView&,
        uint64_t position,
        uint64_t length,
        ShortBaseSequence<uint64_t>&);

}




#endif


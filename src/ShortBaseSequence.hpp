#pragma once

// Shasta.
#include "Base.hpp"
#include "BitCounter.hpp"
#include "bitReversal.hpp"

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include "iostream.hpp"
#include <limits>
#include "stdexcept.hpp"
#include "string.hpp"



namespace shasta {

    // A short sequence of bases.
    // Uses only two integers, so its capacity is limited
    // by the length of the integers used.
    // This class does not keep track of the number of bases
    // actually stored. All unused positions are left set at "A".
    template<class Int> class ShortBaseSequence;
    using ShortBaseSequence8 = ShortBaseSequence<uint8_t>;
    using ShortBaseSequence16 = ShortBaseSequence<uint16_t>;
    using ShortBaseSequence32 = ShortBaseSequence<uint32_t>;
    using ShortBaseSequence64 = ShortBaseSequence<uint64_t>;
    using ShortBaseSequence128 = ShortBaseSequence<__uint128_t>;
    template<class Int> inline ostream& operator<<(ostream&, const ShortBaseSequence<Int>&);

    void testShortBaseSequence();
}



// A short sequence of bases.
// Uses only two integers, so its capacity is limited
// by the length of the integers used.
// Position 0: the LSB bit of the bases (with base 0 corresponding to the MSB bit).
// Position 1: the MSB bit of the bases (with base 0 corresponding to the MSB bit).
// This class does not keep track of the number of bases
// actually stored. All unused positions are left set at "A".
template<class Int> class shasta::ShortBaseSequence {
public:

    // The number of bases that can be represented equals the number of bits
    // in the Int type.
    static constexpr size_t capacity = BitCounter<Int>::numberOfBits;
    static constexpr size_t capacityMinus1 = capacity - 1;

    // The constructor fills the data with 0, which corresponds to all A's.
    ShortBaseSequence()
    {
        std::fill(data.begin(), data.end(), Int(0));
    }

    // Return the base at a given position.
    Base operator[](uint64_t i) const
    {
        const uint64_t bitIndex = capacityMinus1 - (i & capacityMinus1);
        const uint64_t bit0 = (data[0] >> bitIndex) & 1ULL;
        const uint64_t bit1 = (data[1] >> bitIndex) & 1ULL;
        const uint8_t value = uint8_t((bit1 << 1ULL) + bit0);
        return Base::fromInteger(value);
    }

    // Set the base at a given position.
    void set(uint64_t i, Base base) {
        const uint64_t bitIndex = capacityMinus1 - (i & capacityMinus1);
        const Int mask = Int(1ULL << bitIndex);
        const Int maskComplement = Int(~mask);

        const uint64_t bit0 = (base.value) & 1ULL;
        if(bit0 == 0) {
            data[0] &= maskComplement;
        } else {
            data[0] |= mask;
        }

        const uint64_t bit1 = (base.value >> 1ULL) & 1ULL;
        if(bit1 == 0) {
            data[1] &= maskComplement;
        } else {
            data[1] |= mask;
        }

    }

    // Return an integer consisting of the concatenation
    // of the base bits corresponding to the first n bases.
    using Int2 = typename BitCounter<Int>::doubleSizeType;
    Int2 id(uint64_t n) const
    {
        const uint64_t shift = capacity - n;
        const Int2 lsb = data[0] >> shift;
        const Int2 msb = data[1] >> shift;
        return (msb << n) | lsb;
    }

    // Opposite of the above: construct the sequence given the id.
    ShortBaseSequence(Int2 id, uint64_t n)
    {
        const Int2 mask = (Int2(1) << n) - Int2(1);
        const uint64_t shift = capacity - n;
        data[0] = Int((id & mask) << shift);
        data[1] = Int(((id >> n) & mask) << shift);
    }

    // Return the reverse complement of the first n bases.
    ShortBaseSequence<Int> reverseComplementSlow(uint64_t n) const
    {
        ShortBaseSequence<Int> reverseComplementedSequence;
        for(size_t i=0; i<n; i++) {
            const Base b = (*this)[i].complement();
            reverseComplementedSequence.set(n-i-1, b);
        }
        return reverseComplementedSequence;
    }



    // Return the reverse complement of the first n bases.
    // Use bit reversal for speed. This avoids a loop over the n bases.
    ShortBaseSequence<Int> reverseComplement(uint64_t n) const
    {
        const Int shift = Int(capacity - n);
        const Int mask = Int(1ULL << n) - Int(1);
        ShortBaseSequence<Int> reverseComplementedSequence;
        reverseComplementedSequence.data[0] = Int(((~bitReversal(data[0])) & mask) << shift);
        reverseComplementedSequence.data[1] = Int(((~bitReversal(data[1])) & mask) << shift);

#if 0
        // Testing.
        SHASTA_ASSERT(reverseComplementedSequence == reverseComplementSlow(n));
        SHASTA_ASSERT(reverseComplementedSequence.reverseComplementSlow(n) == *this);
#endif

        return reverseComplementedSequence;
    }



    bool operator==(const ShortBaseSequence<Int>& that) const
    {
        return data == that.data;
    }

    bool operator<(const ShortBaseSequence<Int>& that) const
    {
        return data < that.data;
    }

    bool operator<=(const ShortBaseSequence<Int>& that) const
    {
        return data < that.data;
    }

    // Write the first n bases.
    ostream& write(ostream& s, uint64_t n) const
    {
        for(uint64_t i=0; i<n; i++) {
            s << (*this)[i];
        }
        return s;
    }

    void shiftLeft() {
        data[0] = Int(data[0] << 1);
        data[1] = Int(data[1] << 1);
    }

    void shiftRight() {
        data[0] = Int(data[0] >> 1);
        data[1] = Int(data[1] >> 1);
    }



    // A k-mer is canonical if it is <= than its reverse complement.
    bool isCanonical(uint64_t k) const
    {
        return *this <= reverseComplement(k);
    }

    // A k-mer is palindromic if it is equal to its reverse complement.
    // A palindromic k-mer is also canonical.
    bool isPalindromic(uint64_t k) const
    {
        return *this == reverseComplement(k);
    }

    void classify(uint64_t k, bool& isCanonical, bool& isPalindromic) const
    {
        const ShortBaseSequence<Int> rc = reverseComplement(k);
        isCanonical   = ( (*this) <= rc );
        isPalindromic = ( (*this) == rc );
    }



    // Return the longest homopolymer length in the first k bases.
    uint64_t maxHomopolymerLength(uint64_t k) const
    {
        uint64_t maxLength = 0;

        for(uint64_t i=0; i<k-1; /* Increment later */) {
            const Base base = (*this)[i];

            // Find the first base not equal to this one.
            uint64_t j = i + 1;
            while(j < k and (*this)[j] == base) {
                ++j;
            }
            maxLength = max(maxLength, j - i);
            i = j;
        }

        return maxLength;
    }



    // Looks for exact repeats of period N in the first k bases.
    // Returns the longest number of complete copies found.
    // For example, for N=3 and k=10, given the following:
    // ACTCATCATCG
    // it returns 2, because of the two complete copies TCATCA = (TCA)2
    template<uint64_t N> uint64_t countExactRepeatCopies(uint64_t k) const
    {
        // Handle trivial cases.
        if(k < N) {
            return 0;
        }
        if(k == N) {
            return 1;
        }

        // Gather the bases.
        array<Base, capacity> sequence;
        for(uint64_t i=0; i<k; i++) {
            sequence[i] = (*this)[i];
        }

        // Gather N-mers starting at each position.
        array<array<Base, N>, capacity> nMers;
        for(uint64_t i=0; i <= k-N; i++) {
            for(uint64_t j=0; j<N; j++) {
                nMers[i][j] = sequence[i+j];
            }
        }

        // Starting at each position, count the number of consecutive copies.
        uint64_t maxOffset = 0;
        for(uint64_t i=0; i <= k-N; i++) {
            uint64_t j = i + N;
            for(; j<= k-N; j+=N) {
                if(nMers[j] != nMers[i]) {
                    break;
                }
            }
            const uint64_t offset = j - i;
            maxOffset = max(maxOffset, offset);
        }

        return maxOffset / N ;
    }



    // The data are left public to facilitate low level custom code.
    array<Int, 2> data;
};



template<class Int> inline std::ostream& shasta::operator<<(
    std::ostream& s,
    const shasta::ShortBaseSequence<Int>& sequence)
{
    for(size_t i=0; i<sequence.capacity; i++) {
        s << sequence[i];
    }
    return s;
}


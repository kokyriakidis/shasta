#pragma once


// BitCounter is a traits class that can be used to define the
// number of bits in an integer type.
// It can be used for built-in unsigned integer types and
// also for unsigned integers defined by the Boost Multiprecision library.
// std::numeric_limits::digits only works for built-in types
// and cannot be extended.

#ifdef SHASTA_LONG_MARKERS
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include <cstdint.hpp>

namespace shasta {

    template<class Int> class BitCounter {
    public:
        static constexpr int numberOfBits = 0;
        using doubleSizeType = void;
    };

    template<> class BitCounter<uint8_t> {
    public:
        static constexpr int numberOfBits = 8;
        using doubleSizeType = uint16_t;
    };

    template<> class BitCounter<uint16_t> {
    public:
        static constexpr int numberOfBits = 16;
        using doubleSizeType = uint32_t;
    };

    template<> class BitCounter<uint32_t> {
    public:
        static constexpr int numberOfBits = 32;
        using doubleSizeType = uint64_t;
    };

    template<> class BitCounter<uint64_t> {
    public:
        static constexpr int numberOfBits = 64;
        using doubleSizeType = __uint128_t;
    };

    template<> class BitCounter<__uint128_t> {
    public:
        static constexpr int numberOfBits = 128;
#ifdef SHASTA_LONG_MARKERS
        using doubleSizeType = boost::multiprecision::uint256_t;
#endif
    };

#ifdef SHASTA_LONG_MARKERS
    template<> class BitCounter<boost::multiprecision::uint256_t> {
    public:
        static constexpr int numberOfBits = 256;
        using doubleSizeType = boost::multiprecision::uint512_t;
    };
#endif
}

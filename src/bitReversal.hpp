#pragma once

// See https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel

#include "cstdint.hpp"

namespace shasta {
    inline uint16_t bitReversal(uint16_t);
    inline uint32_t bitReversal(uint32_t);
    inline uint64_t bitReversal(uint64_t);
    inline __uint128_t bitReversal(__uint128_t);

    void testBitReversal();
}



inline uint16_t shasta::bitReversal(uint16_t x)
{
    const uint16_t m1 = uint16_t(0x5555);
    const uint16_t m2 = uint16_t(0x3333);
    const uint16_t m4 = uint16_t(0x0F0F);

    x = ((x >> 1) & m1) | ((x & m1) << 1);
    x = ((x >> 2) & m2) | ((x & m2) << 2);
    x = ((x >> 4) & m4) | ((x & m4) << 4);
    x = (x >> 8) | (x << 8);
    return x;
}



inline uint32_t shasta::bitReversal(uint32_t x)
{
    x = ((x >> 1) & 0x55555555) | ((x & 0x55555555) << 1);
    x = ((x >> 2) & 0x33333333) | ((x & 0x33333333) << 2);
    x = ((x >> 4) & 0x0F0F0F0F) | ((x & 0x0F0F0F0F) << 4);
    x = ((x >> 8) & 0x00FF00FF) | ((x & 0x00FF00FF) << 8);
    x = ( x >> 16) | ( x << 16);
    return x;
}



inline uint64_t shasta::bitReversal(uint64_t x)
{
    x = ((x >> 1)  & 0x5555555555555555UL) | ((x & 0x5555555555555555UL) << 1 );
    x = ((x >> 2)  & 0x3333333333333333UL) | ((x & 0x3333333333333333UL) << 2 );
    x = ((x >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((x & 0x0F0F0F0F0F0F0F0FUL) << 4 );
    x = ((x >> 8)  & 0x00FF00FF00FF00FFUL) | ((x & 0x00FF00FF00FF00FFUL) << 8 );
    x = ((x >> 16) & 0x0000FFFF0000FFFFUL) | ((x & 0x0000FFFF0000FFFFUL) << 16);
    x = (x >> 32) | (x << 32);
    return x;
}



inline __uint128_t shasta::bitReversal(__uint128_t x)
{
    static const __uint128_t x1  = 0x5555555555555555UL;
    static const __uint128_t x2  = 0x3333333333333333UL;
    static const __uint128_t x4  = 0x0F0F0F0F0F0F0F0FUL;
    static const __uint128_t x8  = 0x00FF00FF00FF00FFUL;
    static const __uint128_t x16 = 0x0000FFFF0000FFFFUL;
    static const __uint128_t x32 = 0x00000000FFFFFFFFUL;

    static const __uint128_t y1  = (x1  << 64) | x1;
    static const __uint128_t y2  = (x2  << 64) | x2;
    static const __uint128_t y4  = (x4  << 64) | x4;
    static const __uint128_t y8  = (x8  << 64) | x8;
    static const __uint128_t y16 = (x16 << 64) | x16;
    static const __uint128_t y32 = (x32 << 64) | x32;

    x = ((x >> 1 ) & y1 ) | ((x & y1 ) << 1 );
    x = ((x >> 2 ) & y2 ) | ((x & y2 ) << 2 );
    x = ((x >> 4 ) & y4 ) | ((x & y4 ) << 4 );
    x = ((x >> 8 ) & y8 ) | ((x & y8 ) << 8 );
    x = ((x >> 16) & y16) | ((x & y16) << 16);
    x = ((x >> 32) & y32) | ((x & y32) << 32);
    x = (x >> 64) | (x << 64);

    return x;
}


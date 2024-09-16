#pragma once

// Simple template class to select an integer type give its size in bytes.

namespace shasta {

    template<int N> struct UintBySize{};

    template<> struct UintBySize<1>  {using type =   uint8_t  ;};
    template<> struct UintBySize<2>  {using type =   uint16_t ;};
    template<> struct UintBySize<4>  {using type =   uint32_t ;};
    template<> struct UintBySize<8>  {using type =   uint64_t ;};
    template<> struct UintBySize<16> {using type = __uint128_t;};

}

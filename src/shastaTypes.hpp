#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include <cstdint.hpp>

namespace shasta {

    // A KmerId requires twice the number of bits of the integer type
    // used to define ShortBaseSequence.
    using KmerId16 = uint32_t;
    using KmerId32 = uint64_t;
    using KmerId64 = __uint128_t;
    using KmerId128 = boost::multiprecision::uint256_t;
    using KmerId = KmerId64;

    using ReadId = uint32_t;
    using Strand = ReadId;

    using MarkerId = uint64_t;
    using MarkerGraphVertexId = uint64_t;
    using MarkerGraphEdgeId = uint64_t;

    using AssemblyGraphVertexId = uint64_t;
    using AssemblyGraphEdgeId = uint64_t;
}


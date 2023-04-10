#ifndef SHASTA_MARKER_ACCESS_FUNCTIONS_HPP

#include "Kmer.hpp"
#include "ReadId.hpp"

namespace shasta {

    class CompressedMarker;
    class Reads;
    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }

    // Access functions for markers Kmers and KmerIds.
    // There are similar member functions in class Assembler,
    // but these are accessible anywhere else.

    Kmer getOrientedReadMarkerKmer(
        OrientedReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
        );

    Kmer getOrientedReadMarkerKmerStrand0(
        ReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
        );

    Kmer getOrientedReadMarkerKmerStrand1(
        ReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
        );

    // Get the marker KmerId for an oriented read and ordinal.
    KmerId getOrientedReadMarkerKmerId(
        OrientedReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
        );

}

#endif

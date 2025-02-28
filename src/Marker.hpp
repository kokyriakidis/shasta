#ifndef SHASTA_MARKER_HPP
#define SHASTA_MARKER_HPP


/*******************************************************************************

Among all 4^k k-mers of length k, we choose a subset that we call "markers".
The markers are selected at the beginning of an assembly
and never changed, and selected in such a way that,
if (and only if) a k-mer is a marker, its reverse complement
is also a marker.

*******************************************************************************/

#include "Kmer.hpp"
#include "ReadId.hpp"
#include "Uint.hpp"

#include "utility.hpp"

namespace shasta {

    // Classes that will be used to represent markers
    // when the restructuring of marker storage is complete.
    class CompressedMarker;
    class Marker;
    class MarkerWithOrdinal;

    // A pair the identifies a marker by its OrientedRead and ordinal.
    using MarkerDescriptor = pair<OrientedReadId, uint32_t>;
}



// Markers in shared memory are stored using class CompressedMarker.
class shasta::CompressedMarker {
public:

    // The position of this marker in the oriented read.
    // This limits the length of a read to 2^24=16Mib bases.
    uint32_t position;
};



// This stores the same information as CompressedMarker,
// but using built-in, aligned integers.
class shasta::Marker {
public:

    // The id of the k-mer for this marker.
    KmerId kmerId;

    // The position of this marker in the oriented read.
    uint32_t position;

    Marker(KmerId kmerId, uint32_t position) :
        kmerId(kmerId), position(position) {}

    // Default constructor.
    Marker() {}
};



// This also stores the ordinal, that is the index
// of the marker in the oriented read, when the markers
// are sorted by position in the read.
class shasta::MarkerWithOrdinal : public Marker {
public:
    uint32_t ordinal;

    MarkerWithOrdinal(KmerId kmerId, uint32_t position, uint32_t ordinal) :
        Marker(kmerId, position), ordinal(ordinal) {}

    // Default constructor.
    MarkerWithOrdinal() {}

    // Order by kmerId.
    bool operator<(const MarkerWithOrdinal& that) const
    {
        return kmerId < that.kmerId;
    }
};



#endif

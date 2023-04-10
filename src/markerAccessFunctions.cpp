#include "markerAccessFunctions.hpp"
#include "extractKmer.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
using namespace shasta;



Kmer shasta::getOrientedReadMarkerKmer(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
    )
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        return getOrientedReadMarkerKmerStrand0(readId, ordinal, k, reads, markers);
    } else {
        return getOrientedReadMarkerKmerStrand1(readId, ordinal, k, reads, markers);
    }

}



Kmer shasta::getOrientedReadMarkerKmerStrand0(
    ReadId readId,
    uint32_t ordinal0,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
    )
{
    const auto read = reads.getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];

    Kmer kmer0;
    extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);

    return kmer0;
}



Kmer shasta::getOrientedReadMarkerKmerStrand1(
    ReadId readId,
    uint32_t ordinal1,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
    )
{

    // We only have the read stored without reverse complement, so get it from there...
    const auto read = reads.getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    const uint64_t ordinal0 = readMarkerCount - 1 - ordinal1;
    Kmer kmer0;
    extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);

    // ... then do the reverse complement.
    const Kmer kmer1 = kmer0.reverseComplement(k);
    return kmer1;
}



// Get the marker KmerId for an oriented read and ordinal.
KmerId shasta::getOrientedReadMarkerKmerId(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers
    )
{
    const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal, k, reads, markers);
    return KmerId(kmer.id(k));
}

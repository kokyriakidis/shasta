#pragma once

#include "MarkerInterval.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"

#include "cstdint.hpp"
#include "span.hpp"

namespace shasta {

    class MarkerGraph;
    class MarkerInterval;


    // The main input to mode 3 assembly is a set of anchors.
    // Each anchor consists of a span of MarkerIntervals, with the following requirements:
    // - All MarkerIntervals correspond to exactly the same sequence in the corresponding oriented reads, and:
    //      * Those portions of the oriented reads are believed to be aligned.
    //      * They apear in a low number of copies in the genome being sequenced.
    // - There are no duplicate oriented reads in an anchor.
    // - The anchor coverage (number of oriented reads) is in [minPrimaryCoverage, maxPrimaryCoverage].
    // For now the anchors are simply a reference to assembler.markerGraph.edgeMarkerIntervals,
    // but it might be possible to construct the anchors by other means.

    namespace mode3 {

        using AnchorId = uint64_t;
        class Anchor;
        class Anchors;

    }
}



// An Anchor defines a set of MarkerIntervals accessible via this public interface.
// Internals are kept private to facilitate restructuring.
class shasta::mode3::Anchor : private span<const MarkerInterval> {
private:
    using BaseClass = span<const MarkerInterval>;
public:

    Anchor(const span<const MarkerInterval>& s) : span<const MarkerInterval>(s) {}

    // Operator[] returns a MarkerInterval by copy, not by reference.
    // This way we can change the internal representation of the Anchor.
    // But this also means we can't use BaseClass::operator[] because that
    // returns a reference to a const MarkerInterval.
    MarkerInterval operator[](uint64_t i) const
    {
        return BaseClass::operator[](i);
    }

    using BaseClass::size;
    using BaseClass::begin;
    using BaseClass::end;
    using BaseClass::front;
    using BaseClass::empty;

    uint64_t coverage() const
    {
        return size();
    }

    // Return the number of common oriented reads with another Anchor.
    uint64_t countCommon(const Anchor& that) const;
};



class shasta::mode3::Anchors :
    public MultithreadedObject<Anchors>,
    public MappedMemoryOwner {
public:

    // This constructor creates the Anchors from marker graph edges.
    Anchors(
        const MappedMemoryOwner&,
        const MarkerGraph&);

    // This constructor access existing Anchors.
    Anchors(const MappedMemoryOwner&);

    Anchor operator[](AnchorId anchorId) const;
    uint64_t size() const;

    // Return the number of common oriented reads between two Anchors.
    uint64_t countCommon(AnchorId, AnchorId) const;

private:
    MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> anchorMarkerIntervals;
};

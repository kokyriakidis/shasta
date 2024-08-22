#pragma once

#include "cstdint.hpp"
#include "span.hpp"

namespace shasta {

    class MarkerGraph;
    class MarkerInterval;

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }

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
        using Anchor = span<const MarkerInterval>;
        class Anchors;

    }
}



class shasta::mode3::Anchors {
public:

    Anchors(const MarkerGraph&);
    Anchor operator[](AnchorId anchorId) const;
    uint64_t size() const;

private:
    const MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t>& anchorMarkerIntervals;
};

#pragma once

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "utility.hpp"

namespace shasta {
    class Align6;

    class Alignment;
    class AlignmentInfo;
    class KmerCounter;
}



class shasta::Align6 {
public:

    Align6(
        const array<span< pair<KmerId, uint32_t> >, 2>& orientedReadSortedMarkersSpans,
        uint64_t k,
        const KmerCounter&,
        uint64_t maxSkip,
        uint64_t maxDrift,
        Alignment& alignment,
        AlignmentInfo& alignmentInfo,
        ostream& html);
};

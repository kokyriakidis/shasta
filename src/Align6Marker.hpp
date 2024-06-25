#pragma once

// Class Align6Marker is used to pass marker information
// to the Align6 code.

#include "Kmer.hpp"

#include "cstdint.hpp"

namespace shasta {
    class Align6Marker;
}

class shasta::Align6Marker {
public:
    KmerId kmerId;
    uint32_t ordinal;

    // Global frequency as obtained from the KmerCounter.
    // If the KmewrCounter returns a value not representable by a uint32_t,
    // it is replaced with largest uint32_t.
    uint32_t globalFrequency;
    void setGlobalFrequency(uint64_t globalFrequencyLong);

    // Ordering by KmerId.
    bool operator<(const Align6Marker& that) const
    {
        return kmerId < that.kmerId;
    }
};

#pragma once

#include <invalid.hpp>

namespace shasta {
    class KmerDistributionInfo;
}



// A class to contain some information about the coverage distribution
// of the marker k-mers.
class shasta::KmerDistributionInfo {
public:

    // The coverage at which the histogram of the k-mer distribution starts to go up.
    uint64_t coverageLow = invalid<uint64_t>;

    // The coverage greater than coverageLow at which the histogram
    // of the k-mer distribution reaches its maximum.
    uint64_t coveragePeak = invalid<uint64_t>;

    // The highest coverage greater than coveragePeak at which the histogram
    // had a larger value than at coverageLow.
    uint64_t coverageHigh = invalid<uint64_t>;

    // The following holds:
    // coverageLow <= coveragePeak <= coverageHigh
};

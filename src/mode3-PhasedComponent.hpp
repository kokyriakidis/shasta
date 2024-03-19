#pragma once

#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class PhasedComponent;
    }
}



// A PhasedComponent is a set of phased diploid bubbles
// in a BubbleChain.
// It is a vector of (bubble position in bubble chain, phase),
// sorted by bubble position in bubble chain.
// The phase can be -1 or +1.
// PhasedComponents are created in such a way that their position ranges
// in the bubble chain are not overlapping.
class shasta::mode3::PhasedComponent : public vector< pair<uint64_t, int64_t> > {
public:
    uint64_t minPositionInBubbleChain;
    uint64_t maxPositionInBubbleChain;
    void sort();
    void computePositionRange();
};


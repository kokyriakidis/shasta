#include "mode3-PhasedComponent.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3;

#include "algorithm.hpp"
#include "limits"



void PhasedComponent::sort()
{
    SHASTA_ASSERT(size() > 1);
    std::sort(begin(), end(), OrderPairsByFirstOnly<uint64_t, int64_t>());
    minPositionInBubbleChain = front().first;
    maxPositionInBubbleChain = back().first;
}



void PhasedComponent::computePositionRange()
{
    minPositionInBubbleChain = std::numeric_limits<uint64_t>::max();
    maxPositionInBubbleChain = 0;
    for(const auto& p: *this) {
        const uint64_t positionInBubbleChain = p.first;
        minPositionInBubbleChain = min(minPositionInBubbleChain, positionInBubbleChain);
        maxPositionInBubbleChain = max(maxPositionInBubbleChain, positionInBubbleChain);
    }
}

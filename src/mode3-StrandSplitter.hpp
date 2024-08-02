#pragma once

// Shasta.
#include "mode3-Anchor.hpp"
#include "Mode3Assembler.hpp"

// Standard library.
#include "vector.hpp"

// Class StrandSplitter is used to separate strands in self-complementary components.

namespace shasta {
    class OrientedReadId;
    namespace mode3 {
        class StrandSplitter;
    }
}



class shasta::mode3::StrandSplitter {
public:

    StrandSplitter(
        vector<OrientedReadId>&,
        vector<AnchorId>&,
        const Anchors&,
        const vector<AnchorId>& reverseComplementAnchor
        );

private:

    vector<OrientedReadId>& orientedReadIds;
    vector<AnchorId>& anchorIds;

    uint64_t getReverseComplementOrientedReadIndex(uint64_t orientedReadIndex) const;

    // A vector that given the anchor index for the reverse complement
    // of a given anchor index.
    vector<uint64_t> reverseComplementAnchorIndex;

    uint64_t getOrientedReadIndex(OrientedReadId) const;
    uint64_t getAnchorIndex(AnchorId) const;
};

#pragma once

// This file contains declarations of concrete classes derived from mode3::Detangler.

#include "mode3-Detangler.hpp"

namespace shasta {
    namespace mode3 {
        class ChainPermutationDetangler;
    }
}



// This uses a chi squared test to detangle superbubbles
// where the number of incoming Chains equal the number of outgoing Chains
// and is greater than 1. It applies the chi squared test to all
// possible permutations.
class shasta::mode3::ChainPermutationDetangler : public ChainDetangler {
public:

    ChainPermutationDetangler(
        bool debug,
        AssemblyGraph&,
        uint64_t nMax,
        double epsilon,
        double maxLogP,
        double minLogPDelta,
        uint64_t minDetangledCoverage);

    bool operator()(const vector<vertex_descriptor>& superbubble);

private:
    bool debug;
    uint64_t nMax;
    double epsilon;
    double maxLogP;
    double minLogPDelta;
    uint64_t minDetangledCoverage;
};

#pragma once

// This file contains declarations of concrete classes derived from mode3::Detangler.

#include "mode3-Detangler.hpp"

namespace shasta {
    namespace mode3 {
        class ChainDetanglerNbyNPermutation;
    }
}



class shasta::mode3::ChainDetanglerNbyNPermutation : public ChainDetangler {
public:

    ChainDetanglerNbyNPermutation(
        bool debug,
        AssemblyGraph&,
        uint64_t nMax,
        double epsilon,
        double maxLogP,
        double minLogPDelta);

    bool operator()(const vector<vertex_descriptor>& superbubble);

private:
    bool debug;
    uint64_t nMax;
    double epsilon;
    double maxLogP;
    double minLogPDelta;
};

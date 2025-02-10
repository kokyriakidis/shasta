#pragma once

// This file contains declarations of concrete classes derived from mode3::Detangler.

#include "mode3-Detangler.hpp"

namespace shasta {
    namespace mode3 {
        class ChainDetangler2by2Permutation;
        class ChainDetanglerNbyNPermutation;
    }
}



class shasta::mode3::ChainDetangler2by2Permutation : public ChainDetangler {
public:

    ChainDetangler2by2Permutation(
        bool debug,
        AssemblyGraph&,
        double epsilon,
        double chiSquareThreshold);

    bool operator()(const vector<vertex_descriptor>& superbubble);

private:
    bool debug;
    double epsilon;
    double chiSquareThreshold;
};



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

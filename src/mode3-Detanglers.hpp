#pragma once

// This file contains declarations of concrete classes derived from mode3::Detangler.

#include "mode3-Detangler.hpp"

namespace shasta {
    namespace mode3 {
        class Detangler2by2;
    }
}



class shasta::mode3::Detangler2by2 : public ChainDetangler {
public:

    Detangler2by2(
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

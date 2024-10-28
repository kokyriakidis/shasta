#pragma once

#include "mode3-AssemblyGraph.hpp"


namespace shasta {
    namespace mode3 {
        class AssemblyGraphPostprocessor;
    }
}



// AssemblyGraph functionality needed only during postprocessing.
class shasta::mode3::AssemblyGraphPostprocessor : public AssemblyGraph {
public:
    AssemblyGraphPostprocessor(
        const string& assemblyStage,
        uint64_t componentIdArgument,
        const Anchors& anchors,
        const Mode3AssemblyOptions& options);
};

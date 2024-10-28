#include "mode3-AssemblyGraphPostprocessor.hpp"
using namespace shasta;
using namespace mode3;


AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const string& assemblyStage,
    uint64_t componentId,
    const Anchors& anchors,
    const Mode3AssemblyOptions& options) :
    AssemblyGraph(assemblyStage, componentId, anchors, options)
{
}

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
        span<const OrientedReadId> orientedReadIds,
        span<const AnchorId> anchorIds,
        const Mode3AssemblyOptions& options);

    // Access functions.
    edge_descriptor getEdge(uint64_t edgeId) const;
    const BubbleChain& getBubbleChain(uint64_t edgeId) const;
    const Bubble& getBubble(
        uint64_t edgeId,
        uint64_t positionInBubbleChain) const;
    const Chain& getChain(
        uint64_t edgeId,
        uint64_t positionInBubbleChain,
        uint64_t indexInBubble,
        uint64_t& ploidy) const;
    const Chain& getChain(const string& chainStringId) const;
    const Chain& getChain(const ChainIdentifier&) const;
    ChainIdentifier getChainIdentifier(const string& chainStringId) const;
    string getChainStringId(const ChainIdentifier&) const;

    static void parseChainStringId(
        const string&,
        uint64_t& componentId,
        uint64_t& bubbleChainId,
        uint64_t& positionInBubbleChain,
        uint64_t& indexInBubble,
        uint64_t& ploidy);

private:

    // Map edge ids (that is, bubble chain ids) to edges.
    std::map<uint64_t, edge_descriptor> edgeIdMap;
    void createEdgeIdMap();

};

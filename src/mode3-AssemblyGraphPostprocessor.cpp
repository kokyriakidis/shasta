// Shasta.
#include "mode3-AssemblyGraphPostprocessor.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>



AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const string& assemblyStage,
    uint64_t componentId,
    const Anchors& anchors,
    span<const OrientedReadId> orientedReadIds,
    span<const AnchorId> anchorIds,
    const Mode3AssemblyOptions& options) :
    AssemblyGraph(assemblyStage, componentId, orientedReadIds, anchorIds, anchors, options)
{
    createEdgeIdMap();
    annotateAnchors();
}



void AssemblyGraphPostprocessor::parseChainStringId(
    const string& s,
    uint64_t& componentId,
    uint64_t& bubbleChainId,
    uint64_t& positionInBubbleChain,
    uint64_t& indexInBubble,
    uint64_t& ploidy)
{
    using Separator = boost::char_separator<char>;
    using Tokenizer = boost::tokenizer<Separator>;
    const Separator separator("-");

    Tokenizer tokenizer(s, separator);
    vector<string> tokens;
    tokens.insert(tokens.begin(), tokenizer.begin(), tokenizer.end());

    if(tokens.size() != 5) {
        throw runtime_error("Invalid chain string id " + s +
            ". Must be of the form a-b-c-d-Pn where a, b, c, d, and n are integers.");
    }

    const string& ploidyToken = tokens[4];
    if((ploidyToken.size() < 2) or (ploidyToken[0] != 'P')) {
        throw runtime_error("Invalid chain string id " + s +
            ". Must be of the form a-b-c-d-Pn where a, b, c, d, and n are integers.");
    }

    try {
        componentId = atoul(tokens[0]);
        bubbleChainId = atoul(tokens[1]);
        positionInBubbleChain = atoul(tokens[2]);
        indexInBubble = atoul(tokens[3]);
        ploidy = atoul(tokens[4].substr(1));
    } catch(std::exception&) {
        throw runtime_error("Invalid chain string id " + s +
            ". Must be of the form a-b-c-d-Pn where a, b, c, d, and n are integers.");
    }
}



void AssemblyGraphPostprocessor::createEdgeIdMap()
{
    const AssemblyGraph& assemblyGraph = *this;

    edgeIdMap.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgeIdMap.insert({assemblyGraph[e].id, e});
    }
}



AssemblyGraph::edge_descriptor AssemblyGraphPostprocessor::getEdge(
    uint64_t edgeId) const
{

    // Look it up in our map.
    auto it = edgeIdMap.find(edgeId);
    if(it == edgeIdMap.end()) {
        throw runtime_error("Invalid bubble chain id " + to_string(edgeId));
    }

    return it->second;
}



const BubbleChain& AssemblyGraphPostprocessor::getBubbleChain(
    uint64_t edgeId) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Look it up in our map.
    auto it = edgeIdMap.find(edgeId);
    if(it == edgeIdMap.end()) {
        throw runtime_error("Invalid bubble chain id " + to_string(edgeId));
    }

    const edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];
    SHASTA_ASSERT(edge.id == edgeId);

    const BubbleChain& bubbleChain = edge;
    return bubbleChain;
}



const Bubble& AssemblyGraphPostprocessor::getBubble(
    uint64_t edgeId,
    uint64_t positionInBubbleChain) const
{
    const BubbleChain& bubbleChain = getBubbleChain(edgeId);
    if(positionInBubbleChain >= bubbleChain.size()) {
        throw runtime_error("Invalid position in bubble chain " + to_string(positionInBubbleChain) +
            ": bubble chain " + to_string(edgeId) + " has " + to_string(bubbleChain.size()) + " bubbles." );
    }

    return bubbleChain[positionInBubbleChain];
}



const Chain& AssemblyGraphPostprocessor::getChain(
    uint64_t edgeId,
    uint64_t positionInBubbleChain,
    uint64_t indexInBubble,
    uint64_t& ploidy) const
{
    const Bubble& bubble = getBubble(edgeId, positionInBubbleChain);
    ploidy = bubble.size();
    if(indexInBubble >= ploidy) {
        throw runtime_error("Invalid index in bubble " + to_string(indexInBubble) +
            ": bubble " + to_string(edgeId) + "-" + to_string(positionInBubbleChain) +
            " has ploidy " + to_string(bubble.size()) + "." );
    }

    return bubble[indexInBubble];

}



const Chain& AssemblyGraphPostprocessor::getChain(const string& chainStringId) const
{
    uint64_t chainComponentId;
    uint64_t edgeId;
    uint64_t positionInBubbleChain;
    uint64_t indexInBubble;
    uint64_t bubblePloidy;
    parseChainStringId(
        chainStringId,
        chainComponentId,
        edgeId,
        positionInBubbleChain,
        indexInBubble,
        bubblePloidy);
    SHASTA_ASSERT(chainComponentId == componentId);

    uint64_t ploidy;
    const Chain& chain = getChain(edgeId, positionInBubbleChain, indexInBubble, ploidy);

    // Check the ploidy.
    const BubbleChain& bubbleChain = getBubbleChain(edgeId);
    const Bubble& bubble = getBubble(edgeId, positionInBubbleChain);
    if((bubbleChain.size() == 1) and (bubble.size() == 1)) {
        if(bubblePloidy != 0) {
            throw runtime_error("Invalid chain " + chainStringId + ": ploidy should be 0.");
        }
    } else {
        if(ploidy != bubblePloidy) {
            throw runtime_error("Invalid chain " + chainStringId + ": ploidy should be " + to_string(ploidy));
        }
    }


    return chain;
}



const Chain& AssemblyGraphPostprocessor::getChain(const ChainIdentifier& chainIdentifier) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const BubbleChain& bubbleChain = assemblyGraph[chainIdentifier.e];
    const Bubble& bubble = bubbleChain[chainIdentifier.positionInBubbleChain];
    return bubble[chainIdentifier.indexInBubble];
}



ChainIdentifier AssemblyGraphPostprocessor::getChainIdentifier(const string& chainStringId) const
{
    uint64_t chainComponentId;
    uint64_t edgeId;
    uint64_t positionInBubbleChain;
    uint64_t indexInBubble;
    uint64_t bubblePloidy;
    parseChainStringId(
        chainStringId,
        chainComponentId,
        edgeId,
        positionInBubbleChain,
        indexInBubble,
        bubblePloidy);
    SHASTA_ASSERT(chainComponentId == componentId);

    uint64_t ploidy;
    getChain(edgeId, positionInBubbleChain, indexInBubble, ploidy);

    // Check the ploidy.
    const BubbleChain& bubbleChain = getBubbleChain(edgeId);
    const Bubble& bubble = getBubble(edgeId, positionInBubbleChain);
    if((bubbleChain.size() == 1) and (bubble.size() == 1)) {
        if(bubblePloidy != 0) {
            throw runtime_error("Invalid chain " + chainStringId + ": ploidy should be 0.");
        }
    } else {
        if(ploidy != bubblePloidy) {
            throw runtime_error("Invalid chain " + chainStringId + ": ploidy should be " + to_string(ploidy));
        }
    }



    // Now we can construct the ChainIdentifier.
    ChainIdentifier chainIdentifier;
    auto it = edgeIdMap.find(edgeId);
    if(it == edgeIdMap.end()) {
        throw runtime_error("Invalid chain " + chainStringId + ": bubble chain does nto exist.");
    }
    chainIdentifier.e = it->second;
    chainIdentifier.positionInBubbleChain = positionInBubbleChain;
    chainIdentifier.indexInBubble = indexInBubble;



    return chainIdentifier;
}



string AssemblyGraphPostprocessor::getChainStringId(const ChainIdentifier& chainIdentifier) const
{
    return chainStringId(
        chainIdentifier.e,
        chainIdentifier.positionInBubbleChain,
        chainIdentifier.indexInBubble
        );
}


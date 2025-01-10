// Shasta.
#include "mode3-AssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"

// Standard library.
#include <queue>
#include <set>


#if 0
void AssemblyGraph::detangle2()
{
    AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    uint64_t minCommonCount = 4;
    double minJaccard = 0;
    double minCorrectedJaccard = 0.7;

    // Check that all the BubbleChains consist of a single Chain.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());
    }

    // Gather the first and last internal AnchorId for each segment.
    // Segments without internal AnchorIds are discarded.
    std::map<AnchorId, vector<edge_descriptor> > firstAnchorMap;
    std::map<AnchorId, vector<edge_descriptor> > lastAnchorMap;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        if(chain.size() > 2) {
            const AnchorId anchorId0 = chain[1];
            const AnchorId anchorId1 = chain[chain.size() - 2];
            firstAnchorMap[anchorId0].push_back(e);
            lastAnchorMap[anchorId1].push_back(e);
        }
    }


    ofstream dot("Detangle2.dot");
    dot << "digraph Detangle2 {\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        if(chain.size() > 2) {
            dot << "\"" << bubbleChainStringId(e) << "\""
                " [label=\"" <<
                bubbleChainStringId(e) << "\\n" << chainOffset(chain) <<
                "\"];\n";
        }
    }


    // Loop over pairs (lastAnchorId, firstAnchorId).
    for(const auto& pLast: lastAnchorMap) {

        // Skip if ambiguous (it is the last internal AnchorId of more than one Chain).
        if(pLast.second.size() != 1) {
            continue;
        }
        const edge_descriptor e0 = pLast.second.front();
        const AnchorId anchorId0 = pLast.first;

        for(const auto& pFirst: firstAnchorMap) {

            // Skip if ambiguous (it is the first internal AnchorId of more than one Chain).
            if(pLast.second.size() != 1) {
                continue;
            }
            const edge_descriptor e1 = pFirst.second.front();
            const AnchorId anchorId1 = pFirst.first;

            // Analyze this pair of AnchorIds.
            AnchorPairInfo info;
            anchors.analyzeAnchorPair(anchorId0, anchorId1, info);

            if(info.common == 0) {
                continue;
            }

            cout <<
                bubbleChainStringId(e0) << " " <<
                bubbleChainStringId(e1) << " " <<
                info.common << " " <<
                info.correctedJaccard() << " " <<
                info.offsetInBases << endl;

            if(info.common < minCommonCount) {
                continue;
            }

            if(info.offsetInBases < 0) {
                continue;
            }

            if(info.jaccard() < minJaccard) {
                continue;
            }

            if(info.correctedJaccard() < minCorrectedJaccard) {
                continue;
            }

            dot <<
                "\"" << bubbleChainStringId(e0) << "\"->\"" <<
                bubbleChainStringId(e1) << "\""
                "["
                "label=\"" <<
                info.common << " " << info.correctedJaccard() <<
                "\""
                "]"
                ";\n";
        }
    }

    dot << "}\n";
}
#else



// Detangling with path following.
void AssemblyGraph::detangle2()
{
    AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t chainLengthThreshold = 100000;

    cout << "AssemblyGraph::detangle2 called." << endl;
    cout << "Component " << componentId << endl;
    cout << orientedReadIds.size() << " oriented reads." << endl;
    cout << anchorIds.size() << " anchors." << endl;
    cout << num_vertices(assemblyGraph) << " vertices." << endl;
    cout << num_edges(assemblyGraph) << " edges." << endl;

    // Check that all the BubbleChains consist of a single Chain.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());
    }

    // Make sure we have anchor annotations.
    annotateAnchors();

    // Gather the long chains.
    std::set<edge_descriptor> longChains;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        if(chainOffset(chain) >= chainLengthThreshold) {
            longChains.insert(e);
        }
    }


    // Main loop over the long chains.
    for(const edge_descriptor e0: longChains) {

        // Debugging.
        /*
        const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
        if(edge0.id != 138) {
            continue;
        }
        */

        edge_descriptor e1;
        edge_descriptor eNew;
        detangle2Chain(e0, e1, eNew, longChains);
    }
}



// This attempts to detangle the e0 Chain by finding a forward path
// that connects it to another Chain e1.
// If successful, it stitches e0 and e1 together to create
// a new Chain eNew and returns true.
// Otherwise, it returns false and e1 and eNew are not set.
bool AssemblyGraph::detangle2Chain(
    edge_descriptor e0,
    edge_descriptor& /* e1 */,
    edge_descriptor& /* eNew */,
    const std::set<edge_descriptor>& longChains)
{
    // cout << "Working on " << bubbleChainStringId(e0) << endl;


    // Forward path following.
    std::set<edge_descriptor> longChainsFound0;
    detangle2PathFollowing(e0, 0, longChains, longChainsFound0);
    // cout << "Forward path following from " << bubbleChainStringId(e0) << " found:";
    if(longChainsFound0.empty()) {
        // cout << " nothing";
    } else {
        for(const edge_descriptor e: longChainsFound0) {
            // cout << " " << bubbleChainStringId(e);
            cout << "\"" << bubbleChainStringId(e0) << "\"->\"" << bubbleChainStringId(e) << "\";\n";
        }
    }
    // cout << endl;

    // Backward path following.
    std::set<edge_descriptor> longChainsFound1;
    detangle2PathFollowing(e0, 1, longChains, longChainsFound1);
    // cout << "Backward path following from " << bubbleChainStringId(e0) << " found:";
    if(longChainsFound1.empty()) {
        // cout << " nothing";
    } else {
        for(const edge_descriptor e: longChainsFound1) {
            // cout << " " << bubbleChainStringId(e);
            cout << "\"" << bubbleChainStringId(e) << "\"->\"" << bubbleChainStringId(e0) << "\";\n";
        }
    }
    // cout << endl;


    return false;
}



void AssemblyGraph::detangle2PathFollowing(
    edge_descriptor e0,
    uint64_t direction,
    const std::set<edge_descriptor>& longChains,
    std::set<edge_descriptor>& longChainsFound
    ) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const Chain& chain = assemblyGraph[e0].getOnlyChain();

    cout << bubbleChainStringId(e0) << ((direction==0) ? " forward" : " backward") << endl;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t seedAnchorCount = 10;
    uint64_t minCommonCount = 4;
    double minJaccard = 0;
    double minCorrectedJaccard = 0.8;

    class AnchorInfo {
    public:
        AnchorId anchorId;
        uint64_t offset;
        AnchorInfo(AnchorId anchorId, uint64_t offset) : anchorId(anchorId), offset(offset) {}
        bool operator<(const AnchorInfo& that) const {
            return offset > that.offset;
        }
    };

    // We do a recursive search starting from the last/first seedAnchorCount internal anchors of e0.
    std::vector<AnchorInfo> h;
    std::set<AnchorId> anchorIdsEncountered;
    if(direction == 0) {
        for(uint64_t i=0; i<seedAnchorCount; i++) {
            const uint64_t positionInChain = chain.size() - 2 - i;
            if(positionInChain == 0) {
                break;
            }
            const AnchorId anchorId = chain[positionInChain];
            h.push_back(AnchorInfo(anchorId, 0));
            std::push_heap(h.begin(), h.end());
            anchorIdsEncountered.insert(anchorId);
            // cout << "Seed " << anchorIdToString(anchorId) << endl;
        }
    } else {
        for(uint64_t i=0; i<seedAnchorCount; i++) {
            const uint64_t positionInChain = i + 1;
            if(positionInChain == chain.size() - 1) {
                break;
            }
            const AnchorId anchorId = chain[positionInChain];
            h.push_back(AnchorInfo(anchorId, 0));
            std::push_heap(h.begin(), h.end());
            anchorIdsEncountered.insert(anchorId);
            // cout << "Seed " << anchorIdToString(anchorId) << endl;
        }
    }
    // cout << "Starting with " << h.size() << " seeds." << endl;


    // Main recursive loop.
    vector< pair<AnchorId, AnchorPairInfo> > anchorIds1;
    while(not h.empty()) {
        std::pop_heap(h.begin(), h.end());
        const AnchorInfo& anchorInfo0 = h.back();
        const AnchorId anchorId0 = anchorInfo0.anchorId;
        const uint64_t offset0 = anchorInfo0.offset;
        // cout << "Dequeued " << anchorIdToString(anchorId0) << " at offset " << offset0 << endl;
        h.pop_back();

        // Path following starting at anchorId0;
        anchors. followOrientedReads(anchorId0, direction,
            minCommonCount, minJaccard, minCorrectedJaccard, anchorIds1);
        // cout << "followOrientedReads found " << anchorIds1.size() << " anchors." << endl;

        bool foundLongChain = false;
        for(const auto& p: anchorIds1) {
            const AnchorId anchorId1 = p.first;
            const AnchorPairInfo& info1 = p.second;
            // cout << "Found " << anchorIdToString(anchorId1) << " at relative offset " << info1.offsetInBases << endl;
            if(anchorIdsEncountered.contains(anchorId1)) {
                // cout << "Skipped because it has already been encountered." << endl;
            }
            if(not anchorIdsEncountered.contains(anchorId1)) {
                anchorIdsEncountered.insert(anchorId1);
                // cout << "Added " << anchorIdToString(anchorId1) << " to anchorIdsEncountered." << endl;

                // Check the annotation for this anchor.
                const uint64_t localAnchorId1 = anchors.getLocalAnchorIdInComponent(anchorId1);
                const auto& internalChainInfo = anchorAnnotations[localAnchorId1].internalChainInfo;
                if(internalChainInfo.size() == 1)  {
                    const pair<ChainIdentifier, uint64_t>& p = internalChainInfo.front();
                    const ChainIdentifier& chainIdentifier = p.first;
                    SHASTA_ASSERT(chainIdentifier.positionInBubbleChain == 0);
                    SHASTA_ASSERT(chainIdentifier.indexInBubble == 0);
                    const edge_descriptor e1 = chainIdentifier.e;
                    if((e1 != e0) and longChains.contains(e1)) {
                        longChainsFound.insert(e1);
                        // cout << "Reached " << bubbleChainStringId(e1) << " at offset " << offset0 + info1.offsetInBases << endl;
                        if(longChainsFound.size() >= 1) {
                            foundLongChain = true;
                            break;
                        }
                    }
                }

                h.push_back(AnchorInfo(anchorId1, offset0 + info1.offsetInBases));
                std::push_heap(h.begin(), h.end());
                // cout << "Queued  " << anchorIdToString(anchorId1) << " at offset " << info1.offsetInBases << endl;

            }
            if(foundLongChain) {
                break;
            }
        }
        if(foundLongChain) {
            break;
        }
    }

}

#endif

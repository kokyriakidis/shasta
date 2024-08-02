// Shasta.
#include "mode3-StrandSplitter.hpp"
#include "MarkerInterval.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <iostream.hpp>
#include <map>
#include <random>



StrandSplitter::StrandSplitter(
    vector<OrientedReadId>& orientedReadIds,
    vector<AnchorId>& anchorIds,
    const Anchors& anchors,
    const vector<AnchorId>& reverseComplementAnchor) :
    orientedReadIds(orientedReadIds),
    anchorIds(anchorIds)
{
    cout << "Splitting a self-complementary component with " <<
        orientedReadIds.size() << " oriented reads and " <<
        anchorIds.size() << " anchors." << endl;

    // Make sure the orientedReadIds and anchorIds vectors are sorted.
    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    sort(anchorIds.begin(), anchorIds.end());

    // Fill in the reverseComplementAnchorIndex vector.
    reverseComplementAnchorIndex.resize(anchorIds.size());
    for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
        const AnchorId anchorId = anchorIds[anchorIndex];
        const AnchorId anchorIdRc = reverseComplementAnchor[anchorId];
        reverseComplementAnchorIndex[anchorIndex] = getAnchorIndex(anchorIdRc);
    }

    // Loop over anchors to store pairs
    // (Index in orientedReadIds vector, Index in anchorIds vector).
    vector< pair<uint64_t, uint64_t> > pairs;
    for(const AnchorId anchorId: anchorIds) {
        const int64_t anchorIndex = getAnchorIndex(anchorId);
        const Anchor anchor = anchors[anchorId];
        for(const MarkerInterval& markerInterval: anchor) {
            pairs.push_back({
                getOrientedReadIndex(markerInterval.orientedReadId),
                anchorIndex
                });
        }
    }
    cout << "There are " << pairs.size() << " (oriented read, anchor) pairs." << endl;

    // Random source for shuffles below.
    const uint32_t seed = 231;
    std::mt19937 randomSource(seed);

    // Vector to contain the new OrientedReadIds and AnchorIds as obtained from the best cut.
    vector<OrientedReadId> newOrientedReadIds;
    vector<AnchorId> newAnchorIds;
    uint64_t bestCutSize = std::numeric_limits<uint64_t>::max();

    const uint64_t iterationCount = 100;
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {

        // At each iteration we process the pairs in a different, random order.
        std::shuffle(pairs.begin(), pairs.end(), randomSource);

        // There is one "vertex" for each oriented read and one "vertex" for each anchor.
        // Initialize a disjoint sets data structure with all vertices isolated.
        const uint64_t n = orientedReadIds.size() + anchorIds.size();
        vector<uint64_t> rank(n);
        vector<uint64_t> parent(n);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<n; i++) {
            disjointSets.make_set(i);
        }



        // Loop over pairs (orientedReadIndex, anchorIndex).
        // For each pair, "join" the two vertices corresponding to the pair
        // and the two vertices corresponding to the reverse complement of the pair,
        // but only if doing so keeps the oriented read and its reverse complement
        // in seperate connected components.
        const uint64_t offset = orientedReadIds.size();
        uint64_t cutSize = 0;
        for(const auto& p: pairs) {
            const uint64_t orientedReadIndex = p.first;
            const uint64_t anchorIndex = p.second;

            const uint64_t orientedReadIndexRc = getReverseComplementOrientedReadIndex(orientedReadIndex);
            const uint64_t anchorIndexRc = reverseComplementAnchorIndex[anchorIndex];

            const uint64_t orientedReadComponent = disjointSets.find_set(orientedReadIndex);
            const uint64_t orientedReadComponentRc = disjointSets.find_set(orientedReadIndexRc);
            SHASTA_ASSERT(orientedReadComponent != orientedReadComponentRc);

            const uint64_t anchorComponent = disjointSets.find_set(anchorIndex + offset);
            const uint64_t anchorComponentRc = disjointSets.find_set(anchorIndexRc + offset);

            // We want to "join":
            // - orientedReadIndex with anchorIndex
            // - orientedReadIdIndexRc with anchorIndexRc.
            // But we can only do it if this keeps orientedReadIndex and orientedReadIndexRc
            // in different components.
            if(
                orientedReadComponent != anchorComponentRc and
                orientedReadComponentRc != anchorComponent and
                anchorComponent != anchorComponentRc
                ) {
                disjointSets.union_set(orientedReadIndex, anchorIndex + offset);
                disjointSets.union_set(orientedReadIndexRc, anchorIndexRc + offset);

                SHASTA_ASSERT(disjointSets.find_set(orientedReadIndex) != disjointSets.find_set(orientedReadIndexRc));
            } else {
                cutSize++;
            }

        }
        cout << timestamp << "Cut size at iteration " << iteration << " of " << iterationCount <<
            " is " << cutSize << endl;

        // If this is not better than the best cut we have so far, ignore it.
        if(cutSize > bestCutSize) {
            continue;
        }

        cout << "Storing this as the best cut size." << endl;
        bestCutSize = cutSize;

        // Store the OrientedReadIds on one side of this cut.
        const uint64_t componentId = disjointSets.find_set(0);
        newOrientedReadIds.clear();
        for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadIds.size(); orientedReadIndex++) {
            if(disjointSets.find_set(orientedReadIndex) == componentId) {
                newOrientedReadIds.push_back(orientedReadIds[orientedReadIndex]);
            }
        }
        SHASTA_ASSERT(newOrientedReadIds.size() == orientedReadIds.size()/2);

        // Store the AnchorIds on one side of this cut.
        newAnchorIds.clear();
        for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
            if(disjointSets.find_set(anchorIndex + offset) == componentId) {
                newAnchorIds.push_back(anchorIds[anchorIndex]);
            }
        }
        SHASTA_ASSERT(newAnchorIds.size() == anchorIds.size()/2);

    }
    cout << timestamp << "Completed " << iterationCount << " iterations." << endl;
    cout << "Best cut size is " << bestCutSize << endl;

    // Replace the orientedReadIds and anchorIds with the ones we found from the best cut.
    orientedReadIds.swap(newOrientedReadIds);
    anchorIds.swap(newAnchorIds);
}



uint64_t StrandSplitter::getOrientedReadIndex(OrientedReadId orientedReadId) const
{
    const auto it = lower_bound(orientedReadIds.begin(), orientedReadIds.end(), orientedReadId);
    SHASTA_ASSERT(it != orientedReadIds.end());
    SHASTA_ASSERT(*it == orientedReadId);
    return it - orientedReadIds.begin();
}



uint64_t StrandSplitter::getAnchorIndex(AnchorId anchorId) const
{
    const auto it = lower_bound(anchorIds.begin(), anchorIds.end(), anchorId);
    SHASTA_ASSERT(it != anchorIds.end());
    SHASTA_ASSERT(*it == anchorId);
    return it - anchorIds.begin();

}


uint64_t StrandSplitter::getReverseComplementOrientedReadIndex(uint64_t orientedReadIndex) const
{
    if(orientedReadIndex & 1UL) {
        return orientedReadIndex - 1UL;
    } else {
        return orientedReadIndex + 1UL;
    }
}


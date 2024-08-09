// Shasta.
#include "mode3-StrandSplitter.hpp"
#include "deduplicate.hpp"
#include "MarkerInterval.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include "iostream.hpp"
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



    // Write it out in graphviz format.
    {
        const string graphName = "StrandSplitter_" + to_string(anchorIds.front());
        ofstream dot(graphName + ".dot");
        dot << "graph " << graphName << " {\n";
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            dot << "\"" << orientedReadId << "\" [color=red];\n";
        }
        for(const AnchorId anchorId: anchorIds) {
            dot << anchorId << " [color=green];\n";
        }
        for(const auto& p: pairs) {
            const uint64_t orientedReadIndex = p.first;
            const uint64_t anchorIndex = p.second;
            const OrientedReadId orientedReadId = orientedReadIds[orientedReadIndex];
            dot << "\"" << orientedReadId <<
                "\"--" << anchorIds[anchorIndex] << ";\n";
        }
        dot << "}\n";
    }

    // For connected components of the bipartite graph.
    // Oriented reads start at 0.
    // Anchors start at offset.
    const uint64_t offset = orientedReadIds.size();
    const uint64_t n = orientedReadIds.size() + anchorIds.size();


    // Verify that it has a single connected component.
    {
        vector<uint64_t> rank(n);
        vector<uint64_t> parent(n);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<n; i++) {
            disjointSets.make_set(i);
        }
        for(const auto& p: pairs) {
            const uint64_t orientedReadIndex = p.first;
            const uint64_t anchorIndex = p.second;
            disjointSets.union_set(orientedReadIndex, anchorIndex + offset);
        }
        vector<uint64_t> componentIds;
        for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadIds.size(); orientedReadIndex++) {
            componentIds.push_back(disjointSets.find_set(orientedReadIndex));
        }
        for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
            componentIds.push_back(disjointSets.find_set(anchorIndex + offset));
        }
        deduplicate(componentIds);
        SHASTA_ASSERT(componentIds.size() == 1);
    }


    // Random source for shuffles below.
    const uint32_t seed = 231;
    std::mt19937 randomSource(seed);

    // Vector to contain the new OrientedReadIds and AnchorIds as obtained from the best cut.
    vector<OrientedReadId> newOrientedReadIds;
    vector<AnchorId> newAnchorIds;

    // Vector to contain the cut "edges" at the current iteration
    // and the best cut obtained so far.
    // Stores pairs
    // (Index in orientedReadIds vector, Index in anchorIds vector).
    vector< pair<uint64_t, uint64_t> > cut;
    vector< pair<uint64_t, uint64_t> > bestCut;
    vector<uint64_t> bestCutComponents(n);

    const uint64_t iterationCount = 10;
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {
        cut.clear();

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
                cut.push_back(p);
            }

        }
        cout << timestamp << "Cut size at iteration " << iteration << " of " << iterationCount <<
            " is " << cut.size() << endl;

        // If this is not better than the best cut we have so far, ignore it.
        if((iteration == 0) or (cut.size() < bestCut.size())) {
            cout << "Storing this as the best cut size." << endl;
            swap(cut, bestCut);
            for(uint64_t i=0; i<n; i++) {
                bestCutComponents[i] = disjointSets.find_set(i);
            }
        }
    }


#if 0
    // Store the OrientedReadIds on one side of this cut.
    const uint64_t componentId = bestCutComponents[0];
    newOrientedReadIds.clear();
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadIds.size(); orientedReadIndex++) {
        if(bestCutComponents[orientedReadIndex] == componentId) {
            newOrientedReadIds.push_back(orientedReadIds[orientedReadIndex]);
        }
    }
    SHASTA_ASSERT(newOrientedReadIds.size() == orientedReadIds.size()/2);

    // Store the AnchorIds on one side of this cut.
    newAnchorIds.clear();
    for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
        if(bestCutComponents[anchorIndex + offset] == componentId) {
            newAnchorIds.push_back(anchorIds[anchorIndex]);
        }
    }
    SHASTA_ASSERT(newAnchorIds.size() == anchorIds.size()/2);
#else

    // We only store oriented reads and anchors that are on one side of the cut
    // and that don't participate in any cut edge.
    vector<bool> keepOrientedRead(orientedReadIds.size(), false);
    vector<bool> keepAnchor(anchorIds.size(), false);
    const uint64_t componentId = bestCutComponents[0];
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadIds.size(); orientedReadIndex++) {
        if(bestCutComponents[orientedReadIndex] == componentId) {
            keepOrientedRead[orientedReadIndex] = true;
        }
    }
    for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
        if(bestCutComponents[anchorIndex + offset] == componentId) {
            keepAnchor[anchorIndex] = true;
        }
    }
    for(const auto& p: cut) {
        const uint64_t orientedReadIndex = p.first;
        const uint64_t anchorIndex = p.second;
        keepOrientedRead[orientedReadIndex] = false;
        keepAnchor[anchorIndex] = false;
    }
    newOrientedReadIds.clear();
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadIds.size(); orientedReadIndex++) {
        if(keepOrientedRead[orientedReadIndex]) {
            newOrientedReadIds.push_back(orientedReadIds[orientedReadIndex]);
        }
    }
    newAnchorIds.clear();
    for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
        if(keepAnchor[anchorIndex]) {
            newAnchorIds.push_back(anchorIds[anchorIndex]);
        }
    }
    cout << "Best cut so far keeps " << newOrientedReadIds.size() <<
        " oriented reads and " << newAnchorIds.size() << " anchors." << endl;
    cout << "Oriented read loss fraction " << 1. - double(newOrientedReadIds.size()) / (0.5 * double(orientedReadIds.size())) << endl;
    cout << "Anchor loss fraction " << 1. - double(newAnchorIds.size()) / (0.5 * double(anchorIds.size())) << endl;
#endif

    // Write out the best cut
    ofstream csv("BestCut.csv");
    for(const auto& p: cut) {
        const uint64_t orientedReadIndex = p.first;
        const uint64_t anchorIndex = p.second;
        csv << orientedReadIds[orientedReadIndex] << ",";
        csv << anchorIds[anchorIndex] << "\n";
    }

    cout << timestamp << "Completed " << iterationCount << " iterations." << endl;
    cout << "Best cut size is " << bestCut.size() << endl;

    // Replace the orientedReadIds and anchorIds with the ones we found from the best cut.
    this->orientedReadIds.swap(newOrientedReadIds);
    this->anchorIds.swap(newAnchorIds);

    // Recompute the pairs.
    // Loop over anchors to store pairs
    // (Index in orientedReadIds vector, Index in anchorIds vector).
    pairs.clear();
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

    // Recompute the connected components.
    {
        const uint64_t offset = orientedReadIds.size();
        const uint64_t n = orientedReadIds.size() + anchorIds.size();
        vector<uint64_t> rank(n);
        vector<uint64_t> parent(n);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<n; i++) {
            disjointSets.make_set(i);
        }
        for(const auto& p: pairs) {
            const uint64_t orientedReadIndex = p.first;
            const uint64_t anchorIndex = p.second;
            disjointSets.union_set(orientedReadIndex, anchorIndex + offset);
        }
        vector<uint64_t> componentIds;
        for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadIds.size(); orientedReadIndex++) {
            componentIds.push_back(disjointSets.find_set(orientedReadIndex));
        }
        for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
            componentIds.push_back(disjointSets.find_set(anchorIndex + offset));
        }
        deduplicate(componentIds);
        const uint64_t componentCount = componentIds.size();
        cout << "After the cut, there are " << componentCount << " components." << endl;
    }

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


// Shasta.

#include "bits/stdint-uintn.h"
#include "mode3-AssemblyGraph.hpp"
#include "mode3-PhasingTable.hpp"
#include "MarkerInterval.hpp"
#include "orderPairs.hpp"
#include "PngImage.hpp"
#include "shastaTypes.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/multi_index/detail/bidir_node_iterator.hpp>
#include <boost/multi_index/detail/ord_index_impl.hpp>
#include <boost/operators.hpp>

// Standard library.
#include "algorithm.hpp"
#include "filesystem.hpp"
#include "iostream.hpp"
#include <limits>
#include <set>
#include "stdexcept.hpp"
#include "string.hpp"
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"



void AssemblyGraph::writeBubbleChainsPhasingTables(
    const string& fileNamePrefix,
    double phaseErrorThreshold) const
{
    const AssemblyGraph& cGraph = *this;

    const string directoryName = fileNamePrefix + "-PhasingTables";
    if(not std::filesystem::create_directory(directoryName)) {
        throw runtime_error("Could not create directory " + directoryName);
    }


    // Loop over all BubbleChains.
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = cGraph[ce];
        const BubbleChain& bubbleChain = edge;

        // Create the PhasingTable for this bubble chain.
        PhasingTable phasingTable(bubbleChain, anchors, phaseErrorThreshold);

        if(phasingTable.empty()) {
            continue;
        }
        if(phasingTable.bubbleCount() < 2) {
            continue;
        }

        cout << "Phasing table for " << bubbleChainStringId(ce) <<
            " has " << phasingTable.entryCount() <<
            " entries (of which " << phasingTable.ambiguousEntryCount() <<
            " ambiguous) for " <<
            phasingTable.bubbleCount() << " bubbles and " <<
            phasingTable.orientedReadCount() << " oriented reads." << endl;

        const string fileNamePrefix = directoryName + "/" + bubbleChainStringId(ce);
        phasingTable.writeCsv(fileNamePrefix);
        phasingTable.writePng(fileNamePrefix + "-RelativePhase.png",
            PhasingTable::ColoringMethod::byRelativePhase);
        phasingTable.writePng(fileNamePrefix + "-DiscreteRelativePhase.png",
            PhasingTable::ColoringMethod::byDiscreteRelativePhase);

        phasingTable.greedyPhasing();
        phasingTable.writePng(fileNamePrefix + "-Consistency.png",
            PhasingTable::ColoringMethod::byConsistency);

#if 0
        for(uint64_t i=0; i<6; i++) {
            cout << "Discordant count before sweep " << i << " = " << phasingTable.discordantCount() << endl;
            phasingTable.flipSweep();
        }
        cout << "Final discordant count = " << phasingTable.discordantCount() << endl;
        phasingTable.writePng(directoryName + "/" + bubbleChainStringId(ce) + "-sweep.png", false);
        phasingTable.writePng(directoryName + "/" + bubbleChainStringId(ce) + "-sweep-byType.png", true);
#endif
    }
}


PhasingTable::PhasingTable(
    const BubbleChain& bubbleChain,
    const Anchors& anchors,
    double phaseErrorThreshold)
{
    fill(bubbleChain, anchors, phaseErrorThreshold);
    gatherOrientedReads();
    gatherBubbles();
    fillIndexes();
}



void PhasingTable::fill(
    const BubbleChain& bubbleChain,
    const Anchors& anchors,
    double phaseErrorThreshold)
{
    clear();

    // Loop over the bubbles in this bubble chain.
    for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
        const mode3::Bubble& bubble = bubbleChain[positionInBubbleChain];

        // If this bubble is not diploid, skip it.
        if(not bubble.isDiploid()) {
            continue;
        }

        // Loop over the two chains of this diploid bubble.
        for(uint64_t chainIndexInBubble=0; chainIndexInBubble<bubble.size(); chainIndexInBubble++) {
            SHASTA_ASSERT(chainIndexInBubble < 2);
            const Chain& chain = bubble[chainIndexInBubble];


            // Loop over marker graph edges of this chain, excluding the terminal ones.
            SHASTA_ASSERT(chain.size() >= 2);
            for(uint64_t i=1; i<chain.size()-1; i++) {
                const AnchorId anchorId = chain[i];
                const span<const AnchorMarkerInterval> anchor = anchors[anchorId];

                // Loop over AnchorMarkerInterval of this anchor.
                for(const AnchorMarkerInterval& markerInterval: anchor) {
                    const OrientedReadId orientedReadId = markerInterval.orientedReadId;

                    // Access the PhasingTableEntry for this OrientedReadId and
                    // position in the bubble chain, creating it if necessary.
                    auto it = indexByBoth().find(make_tuple(orientedReadId, positionInBubbleChain));
                    if(it == indexByBoth().end()) {
                        tie(it, ignore) = insert(PhasingTableEntry(orientedReadId, positionInBubbleChain));
                    }
                    // Access it as non-const so we can update the frequency array.
                    // We can do a const_cast because we only update the frequency,
                    // which does not participate in any field used to index the PhasingTable.
                    PhasingTableEntry& entry = const_cast<PhasingTableEntry&>(*it);

                    // Increment the PhasingTableEntry for this OrientedReadId and positionInBubbleChain.
                    ++entry.frequency[chainIndexInBubble];
                }
            }
        }
    }

    // Compute the relative phase of all PhasingTableEntries.
    for(const PhasingTableEntry& phasingTableEntry: indexByBoth()) {
        PhasingTableEntry& nonConstPhasingTableEntry = const_cast<PhasingTableEntry&>(phasingTableEntry);
        nonConstPhasingTableEntry.storeRelativePhase(phaseErrorThreshold);
    }
}



void PhasingTable::gatherOrientedReads()
{

    // Gather the distinct OrientedReadIds that appear in this PhasingTable.
    std::set<OrientedReadId> orientedReadIds;
    for(const PhasingTableEntry& phasingTableEntry: indexByBoth()) {
        orientedReadIds.insert(phasingTableEntry.orientedReadId);
    }

    // Store them in the orientedReads vector.
    orientedReads.clear();
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        OrientedRead orientedRead;
        orientedRead.id = orientedReadId;
        orientedReads.push_back(orientedRead);
    }

    // Fill in the min/max positions in the bubble chain.
    for(OrientedRead& orientedRead: orientedReads) {
        orientedRead.minPositionInBubbleChain = std::numeric_limits<uint64_t>::max();
        orientedRead.maxPositionInBubbleChain = 0;
        for(auto it=indexByOrientedReadId().find(orientedRead.id);
            it!=indexByOrientedReadId().end() and it->orientedReadId == orientedRead.id; ++it) {
            const uint64_t positionInBubbleChain = it->positionInBubbleChain;
            orientedRead.minPositionInBubbleChain = min(orientedRead.minPositionInBubbleChain, positionInBubbleChain);
            orientedRead.maxPositionInBubbleChain = max(orientedRead.maxPositionInBubbleChain, positionInBubbleChain);
        }
    }

    // Sort the orientedReads vector by average position.
    vector< pair<uint64_t, uint64_t> > orientedReadsTable; // (index, minPosition + maxPosition)
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        const OrientedRead& orientedRead = orientedReads[i];
        orientedReadsTable.push_back({i, orientedRead.minPositionInBubbleChain + orientedRead.maxPositionInBubbleChain});
    }
    sort(orientedReadsTable.begin(), orientedReadsTable.end(),
        OrderPairsBySecondOnly<uint64_t, uint64_t>());
    vector<OrientedRead> sortedOrientedReads;
    for(const auto& p: orientedReadsTable) {
        sortedOrientedReads.push_back(orientedReads[p.first]);
    }
    orientedReads.swap(sortedOrientedReads);

    // Fill in the orientedReadIdsMap map.
    orientedReadsMap.clear();
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        orientedReadsMap.insert({orientedReads[i].id, i});
    }
}



void PhasingTable::gatherBubbles()
{

    // Gather the positions in the bubble chains of the diploid bubbles
    // that the oriented reads appear in.
    std::set<uint64_t> positionsInBubbleChain;
    for(const PhasingTableEntry& phasingTableEntry: indexByBoth()) {
        positionsInBubbleChain.insert(phasingTableEntry.positionInBubbleChain);
    }

    // Store them in the bubbles vector.
    bubbles.clear();
    for(const uint64_t positionInBubbleChain: positionsInBubbleChain) {
        bubbles.push_back({positionInBubbleChain});
    }

    // Check that the bubbles are sorted by position.
    for(uint64_t i1=1; i1<bubbles.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const Bubble& bubble0 = bubbles[i0];
        const Bubble& bubble1 = bubbles[i1];
        SHASTA_ASSERT(bubble0.positionInBubbleChain < bubble1.positionInBubbleChain);
    }

    // Fill in the bubble map.
    bubblesMap.clear();
    for(uint64_t i=0; i<bubbles.size(); i++) {
        bubblesMap.insert({bubbles[i].positionInBubbleChain, i});
    }

}



// Fill the orientedReadIndex and bubbleIndex in all PhasingTableEntries.
// This can only be done after gatherOrientedReads and gatherBubbles
// have been called.
void PhasingTable::fillIndexes()
{
    for(const PhasingTableEntry& phasingTableEntry: indexByBoth()) {

        // Access the PhasingTableEntry as a non-const reference.
        // This is ok because we will not modify the fields that participate
        // in the PhasingTable indexes.
        PhasingTableEntry& entry = const_cast<PhasingTableEntry&>(phasingTableEntry);

        entry.orientedReadIndex = orientedReadsMap[entry.orientedReadId];
        entry.bubbleIndex = bubblesMap[entry.positionInBubbleChain];
    }

}

#if 0

void PhasingTable::write(const string& fileNamePrefix) const
{
    writeCsv(fileNamePrefix);
    writeHtml(fileNamePrefix);
    writePng(fileNamePrefix + ".png", true);
}
#endif



void PhasingTable::writeCsv(const string& fileNamePrefix) const
{
    writeOrientedReadsCsv(fileNamePrefix);
    writeBubblesCsv(fileNamePrefix, true);
    writeDetailsCsv(fileNamePrefix);
}



void PhasingTable::writeOrientedReadsCsv(const string& fileNamePrefix) const
{
    ofstream csv(fileNamePrefix + "-OrientedReads.csv");
    csv << "OrientedReadId,Min position in bubble chain,Max position in bubble chain,"
        "Oriented read index,Min bubble index,Max bubble Index,\n";

    for(uint64_t i=0; i<orientedReads.size(); i++) {
        const OrientedRead& orientedRead = orientedReads[i];
        csv << orientedRead.id << ",";
        csv << orientedRead.minPositionInBubbleChain << ",";
        csv << orientedRead.maxPositionInBubbleChain << ",";
        csv << i << ",";
        csv << bubblesMap.find(orientedRead.minPositionInBubbleChain)->second << ",";
        csv << bubblesMap.find(orientedRead.maxPositionInBubbleChain)->second << ",";
        csv << "\n";
    }
}



void PhasingTable::writeBubblesCsv(
    const string& fileNamePrefix,
    bool writePhasingInformation) const
{
    ofstream csv(fileNamePrefix + "-Bubbles.csv");
    csv << "Position in bubble chain,Bubble index,Unambiguous,Ambiguous,";
    if(writePhasingInformation) {
        csv << "Consistent,Inconsistent,Error rate,";
    }
    csv << "\n";

    for(uint64_t i=0; i<bubbles.size(); i++) {
        csv << bubbles[i].positionInBubbleChain << ",";
        csv << i << ",";

        uint64_t unambiguous;
        uint64_t ambiguous;
        tie(unambiguous, ambiguous) = countEntriesForBubble(bubbles[i].positionInBubbleChain);
        csv << unambiguous << ",";
        csv << ambiguous << ",";

        if(writePhasingInformation) {
            uint64_t consistent;
            uint64_t inconsistent;
            tie(consistent, inconsistent) = countConsistentEntriesForBubble(bubbles[i].positionInBubbleChain);
            csv << consistent << ",";
            csv << inconsistent << ",";
            csv << double(inconsistent) / double(consistent + inconsistent) << ",";
        }

        csv << "\n";
    }
}



void PhasingTable::writeDetailsCsv(const string& fileNamePrefix) const
{
    ofstream csv(fileNamePrefix + "-Details.csv");

    csv << "Position in bubble chain,OrientedReadId,Bubble index,Oriented read index,Frequency0,Frequency1,"
        "Relative phase,DiscreteRelative phase\n";

    for(const OrientedRead& orientedRead: orientedReads) {
        const OrientedReadId orientedReadId = orientedRead.id;
        for(auto it=indexByOrientedReadId().find(orientedReadId);
            it!=indexByOrientedReadId().end() and it->orientedReadId == orientedReadId; ++it) {
            const PhasingTableEntry& phasingTableEntry = *it;
            phasingTableEntry.writeCsv(csv);
            csv << "\n";
        }
    }
}



void PhasingTableEntry::writeCsv(ostream& csv) const
{
    csv << positionInBubbleChain << ",";
    csv << orientedReadId << ",";
    csv << bubbleIndex << ",";
    csv << orientedReadIndex << ",";
    csv << frequency[0] << ",";
    csv << frequency[1] << ",";
    csv << relativePhase << ",";
    csv << discreteRelativePhase << ",";
}



void PhasingTable::writePng(const string& fileName, ColoringMethod coloringMethod) const
{
    PngImage image{int(bubbleCount()), int(orientedReadCount())};
    for(uint64_t x=0; x<bubbleCount(); x++) {
        for(uint64_t y=0; y<orientedReadCount(); y++) {
            image.setPixel(int(x), int(y), 255, 255, 255);
        }
    }

    for(const PhasingTableEntry& entry: indexByBoth()) {

        int r, g, b;
        if(coloringMethod == ColoringMethod::byDiscreteRelativePhase) {
            switch(entry.discreteRelativePhase) {
            case 0:
                // Ambiguous: black
                r = 0;
                g = 0;
                b = 0;
                break;
            case +1:
                // In-phase: red.
                r = 255;
                g = 0;
                b = 0;
                break;
            case -1:
                // Out-of-phase: blue.
                r = 0;
                g = 0;
                b = 255;
                break;
            default:
                SHASTA_ASSERT(0);
            }

        } else if(coloringMethod == ColoringMethod::byRelativePhase) {

            // Compute (r, g, b) values that give:
            // - Green if relativePhase is 1 (in-phase).
            // - Red if relativePhase is -1 (out-of-phase).
            if(entry.relativePhase >= 0.) {
                r = 255;
                g = 0;
                b = int(std::round((1. - entry.relativePhase) * 255.));
            } else {
                r = int(std::round((1. +  entry.relativePhase) * 255.));
                g = 0;
                b = 255;
            }
        } else if(coloringMethod == ColoringMethod::byConsistency) {
            const int64_t state = consistencyState(entry);
            switch(state) {
            case +1:
                r = 0;
                g = 255;
                b = 0;
                break;
            case -1:
                r = 255;
                g = 0;
                b = 0;
                break;
            case 0:
                r = 255;
                g = 255;
                b = 0;
                break;
            default:
                SHASTA_ASSERT(0);
            }

        } else {
            SHASTA_ASSERT(0);
        }

        image.setPixel(int(entry.bubbleIndex), int(entry.orientedReadIndex), r, g, b);
    }

    image.write(fileName);
}



uint64_t PhasingTable::unambiguousEntryCount() const
{
    const auto& indexByBoth = get<0>();

    uint64_t n = 0;
    for(const PhasingTableEntry& entry: indexByBoth) {
        if(entry.discreteRelativePhase != 0) {
            ++n;
        }
    }
    return n;
}



uint64_t PhasingTable::ambiguousEntryCount() const
{
    const auto& indexByBoth = get<0>();

    uint64_t n = 0;
    for(const PhasingTableEntry& entry: indexByBoth) {
        if(entry.discreteRelativePhase == 0) {
            ++n;
        }
    }
    return n;
}



// Compute the consistency state of a PhasingTableEntry relative
// to the current phases of its oriented read and bubble.
// It can be +1 (consistent), -1 (inconsistent), or 0 (unassigned or ambiguous).
int64_t PhasingTable::consistencyState(const PhasingTableEntry& entry) const
{
    if(entry.discreteRelativePhase == 0) {
        return 0;
    }

    const int64_t orientedReadPhase = orientedReads[entry.orientedReadIndex].phase;
    if(orientedReadPhase == 0) {
        return 0;
    }

    const int64_t bubblePhase = bubbles[entry.bubbleIndex].phase;
    if(bubblePhase == 0) {
        return 0;
    }

    if(entry.discreteRelativePhase == 1) {
        if(orientedReadPhase == bubblePhase) {
            return +1;
        } else {
            return -1;
        }
    } else {
        if(orientedReadPhase == bubblePhase) {
            return -1;
        } else {
            return +1;
        }
    }
}



// Count the number of (consistent,inconsistent) PhasingTableEntries
// for an oriented read based on the phases currently assigned
// to bubbles and oriented reads.
pair<uint64_t, uint64_t> PhasingTable::countConsistentEntriesForOrientedRead(
    OrientedReadId orientedReadId) const
{
    uint64_t consistentCount = 0;
    uint64_t inconsistentCount = 0;

    for(auto it=indexByOrientedReadId().find(orientedReadId);
        it!=indexByOrientedReadId().end() and it->orientedReadId == orientedReadId; ++it) {
        const PhasingTableEntry& entry = *it;

        const int64_t s = consistencyState(entry);
        switch(s) {
        case +1:
            ++consistentCount;
            break;
        case -1:
            ++inconsistentCount;
            break;
        case 0:
            break;
        default:
            SHASTA_ASSERT(0);
        }
    }

    return {consistentCount, inconsistentCount};
}



// Count the number of (consistent,inconsistent) PhasingTableEntries
// for the bubble at a given bubble chain position based on the phases currently assigned
// to bubbles and oriented reads.
pair<uint64_t, uint64_t> PhasingTable::countConsistentEntriesForBubble(uint64_t positionInBubbleChain) const
{
    uint64_t consistentCount = 0;
    uint64_t inconsistentCount = 0;

    for(auto it=indexByPositionInBubbleChain().find(positionInBubbleChain);
        it!=indexByPositionInBubbleChain().end() and it->positionInBubbleChain == positionInBubbleChain; ++it) {
        const PhasingTableEntry& entry = *it;

        const int64_t s = consistencyState(entry);
        switch(s) {
        case +1:
            ++consistentCount;
            break;
        case -1:
            ++inconsistentCount;
            break;
        case 0:
            break;
        default:
            SHASTA_ASSERT(0);
        }
    }

    return {consistentCount, inconsistentCount};

}



pair<uint64_t, uint64_t> PhasingTable::countEntriesForBubble(uint64_t positionInBubbleChain) const
{
    uint64_t unambiguous = 0;
    uint64_t ambiguous = 0;

    for(auto it=indexByPositionInBubbleChain().find(positionInBubbleChain);
        it!=indexByPositionInBubbleChain().end() and it->positionInBubbleChain == positionInBubbleChain; ++it) {
        const PhasingTableEntry& entry = *it;

        if(entry.discreteRelativePhase == 0) {
            ++ambiguous;
        } else {
            ++unambiguous;
        }
    }

    return {unambiguous, ambiguous};

}


// Count the number of (consistent,inconsistent) PhasingTableEntries
// based on the phases currently assigned
// to bubbles and oriented reads.
pair<uint64_t, uint64_t> PhasingTable::countConsistentEntries() const
{
    uint64_t consistentCount = 0;
    uint64_t inconsistentCount = 0;

    for(const PhasingTableEntry& entry: indexByBoth()) {

        const int64_t s = consistencyState(entry);
        switch(s) {
        case +1:
            ++consistentCount;
            break;
        case -1:
            ++inconsistentCount;
            break;
        case 0:
            break;
        default:
            SHASTA_ASSERT(0);
        }
    }

    return {consistentCount, inconsistentCount};

}



// Iteratively optimize the phases of the oriented reads and of the bubbles.
// Experimental. Do not use.
void PhasingTable::simpleIterativePhasing1()
{
    // Start with the phases of all oriented reads and bubbles set to +1.
    for(OrientedRead& orientedRead: orientedReads) {
        orientedRead.phase = +1;
    }
    for(Bubble& bubble: bubbles) {
        bubble.phase = +1;
    }


    // Iteration loop.
    uint64_t consistentCount;
    uint64_t inconsistentCount;
    tie(consistentCount, inconsistentCount) = countConsistentEntries();
    const uint64_t unassignedCount = size() - (consistentCount + inconsistentCount);
    uint64_t oldInconsistentCount = inconsistentCount;
    cout << "Initial consistency statistics: consistent " << consistentCount <<
        ", inconsistent " << inconsistentCount <<
        ", unassigned " << unassignedCount << endl;
    for(uint64_t iteration=0; ; iteration++) {

        // Set the oriented read phases based on the current bubble phases.
        for(OrientedRead& orientedRead: orientedReads) {

            // Count the number of consistent/inconsistent PhasingTableEntries
            // for this bubble.
            tie(consistentCount, inconsistentCount) =
                countConsistentEntriesForOrientedRead(orientedRead.id);

            // Set the phase of this oriented read accordingly.
            if(consistentCount >= inconsistentCount) {
                // Do nothing.
            } else {
                // Flip it.
                orientedRead.phase = - orientedRead.phase;
            }
        }

        // Set the bubble phases based on the current oriented read phases.
        for(Bubble& bubble: bubbles) {

            // Count the number of consistent/inconsistent PhasingTableEntries
            // for this bubble.
            tie(consistentCount, inconsistentCount) =
                countConsistentEntriesForBubble(bubble.positionInBubbleChain);

            const double consistentFraction = double(consistentCount) / double(consistentCount + inconsistentCount);

            // Set the phase of this bubble accordingly.
            if(consistentFraction > 0.2) {
                // Do nothing.
            } else {
                // Flip it.
                bubble.phase = - bubble.phase;
            }
        }

        tie(consistentCount, inconsistentCount) = countConsistentEntries();
        const uint64_t unassignedCount = size() - (consistentCount + inconsistentCount);
        cout << "Consistency statistics after phasing iteration " << iteration <<
            ": consistent " << consistentCount <<
            ", inconsistent " << inconsistentCount <<
            ", unassigned " << unassignedCount << endl;
        SHASTA_ASSERT(inconsistentCount <= oldInconsistentCount);
        if(inconsistentCount == oldInconsistentCount) {
            break;
        }
        oldInconsistentCount = inconsistentCount;
    }
}



// Iteratively optimize the phases of the oriented reads and of the bubbles.
// Experimental. Do not use.
void PhasingTable::simpleIterativePhasing2()
{
    // Start with the phases of all oriented reads and bubbles set to +1.
    for(OrientedRead& orientedRead: orientedReads) {
        orientedRead.phase = +1;
    }
    for(Bubble& bubble: bubbles) {
        bubble.phase = +1;
    }


    // Iteration loop.
    uint64_t consistentCount;
    uint64_t inconsistentCount;
    tie(consistentCount, inconsistentCount) = countConsistentEntries();
    const uint64_t unassignedCount = size() - (consistentCount + inconsistentCount);
    cout << "Initial consistency statistics: consistent " << consistentCount <<
        ", inconsistent " << inconsistentCount <<
        ", unassigned " << unassignedCount << endl;
    vector<uint64_t> consistentBubbles;
    vector<uint64_t> inconsistentBubbles;
    for(uint64_t iteration=0; iteration<6; iteration++) {

        // Loop over oriented reads.
        for(OrientedRead& orientedRead: orientedReads) {

            // Gather the bubbles that have a consistent/inconsistent
            // PhasingTableEntry with this oriented read.
            // Gather the bubbles where this oriented read appears with phase +1 or -1
            // (with tolerance equal to phaseError).
            consistentBubbles.clear();
            inconsistentBubbles.clear();
            for(auto it=indexByOrientedReadId().find(orientedRead.id);
                it!=indexByOrientedReadId().end() and it->orientedReadId == orientedRead.id; ++it) {
                const PhasingTableEntry& phasingTableEntry = *it;
                const int64_t s = consistencyState(phasingTableEntry);

                if(s == +1) {
                    consistentBubbles.push_back(phasingTableEntry.bubbleIndex);
                } else if(s == -1) {
                    inconsistentBubbles.push_back(phasingTableEntry.bubbleIndex);
                }
            }

            // If there are more consistentBubbles than inconsistentBubbles, flip the minusBubbles.
            // If there are more inconsistentBubbles than consistentBubbles, flip the plusBubbles.
            if(consistentBubbles.size() == inconsistentBubbles.size()) {
                continue;
            }
            const vector<uint64_t>& bubblesToFlip =
                (consistentBubbles.size() > inconsistentBubbles.size()) ? inconsistentBubbles : consistentBubbles;
            for(const uint64_t bubbleIndex: bubblesToFlip) {
                Bubble& bubble = bubbles[bubbleIndex];
                bubble.phase = -bubble.phase;
            }
            if(inconsistentBubbles.size() > consistentBubbles.size()) {
                orientedRead.phase = - orientedRead.phase;
            }
        }

        tie(consistentCount, inconsistentCount) = countConsistentEntries();
        const uint64_t unassignedCount = size() - (consistentCount + inconsistentCount);
        cout << "Consistency statistics after phasing iteration " << iteration <<
            ": consistent " << consistentCount <<
            ", inconsistent " << inconsistentCount <<
            ", unassigned " << unassignedCount << endl;
    }
}



void PhasingTable::greedyPhasing()
{
    const bool debug = false;

    class OrientedReadInfo {
    public:

        // Index of this oriented read in the orientedReads vector.
        uint64_t orientedReadIndex;

        // The total number of unambiguous PhasingTableEntries for this oriented read.
        uint64_t unambiguousBubbleCount = 0;

        // The number of bubbles that have already been phased and that have an
        // unambiguous PhasingTableEntry with this oriented read.
        uint64_t phasedUnambiguousBubbleCount = 0;

        OrientedReadInfo(uint64_t orientedReadIndex) :
            orientedReadIndex(orientedReadIndex) {}
    };

    // The OrientedReadTable is a container of OrientedReadInfo
    // used to keep track of unphased oriented reads by various criteria.
    class OrientedReadTable : public boost::multi_index_container<OrientedReadInfo,
        boost::multi_index::indexed_by <

        // Index by orientedReadIndex (unique).
        boost::multi_index::ordered_unique<boost::multi_index::member<
            OrientedReadInfo,
            uint64_t,
            &OrientedReadInfo::orientedReadIndex> >,

            // Index by unambiguousBubbleCount (non-unique, largest first).
            boost::multi_index::ordered_non_unique<boost::multi_index::member<
                OrientedReadInfo,
                uint64_t,
                &OrientedReadInfo::unambiguousBubbleCount>,
                std::greater<uint64_t> >,

            // Index by phasedUnambiguousBubbleCount (non-unique, largest first).
            boost::multi_index::ordered_non_unique<boost::multi_index::member<
                OrientedReadInfo,
                uint64_t,
                &OrientedReadInfo::phasedUnambiguousBubbleCount>,
                std::greater<uint64_t> >
        > > {
    };
    OrientedReadTable orientedReadTable;



    // Initialize the OrientedReadTable.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadCount(); orientedReadIndex++) {
        const OrientedReadId orientedReadId = orientedReads[orientedReadIndex].id;

        OrientedReadInfo orientedReadInfo(orientedReadIndex);
        for(auto it=indexByOrientedReadId().find(orientedReadId);
            it!=indexByOrientedReadId().end() and it->orientedReadId == orientedReadId; ++it) {
            const PhasingTableEntry& phasingTableEntry = *it;
            if(phasingTableEntry.discreteRelativePhase != 0) {
                ++orientedReadInfo.unambiguousBubbleCount;
            }
        }
        orientedReadTable.insert(orientedReadInfo);
    }


    // Initialize the phases and phasing components of all oriented reads and bubbles.
    for(OrientedRead& orientedRead: orientedReads) {
        orientedRead.phase = 0;
        orientedRead.phasingComponent = invalid<uint64_t>;
    }
    for(Bubble& bubble: bubbles) {
        bubble.phase = 0;
        bubble.phasingComponent = invalid<uint64_t>;
    }



    // Outer loop is over phasing components.
    for(uint64_t phasingComponent=0; ; phasingComponent++) {
        if(orientedReadTable.empty()) {
            break;
        }

        // Find the starting oriented read for this phasing component.
        const auto it = orientedReadTable.get<1>().begin();
        const OrientedReadInfo& orientedReadInfo = *it;
        OrientedRead& orientedRead = orientedReads[orientedReadInfo.orientedReadIndex];

        const uint64_t minPositionInBubbleChain = orientedRead.minPositionInBubbleChain;
        const uint64_t maxPositionInBubbleChain = orientedRead.maxPositionInBubbleChain;
        const uint64_t minBubbleIndex = bubblesMap[minPositionInBubbleChain];
        const uint64_t maxBubbleIndex = bubblesMap[maxPositionInBubbleChain];

        if(debug) {
            cout << "Begin phasing component " << phasingComponent << endl;
            cout << "Phasing group begins at " << orientedRead.id <<
                ", index " << orientedReadInfo.orientedReadIndex <<
                " with " << orientedReadInfo.unambiguousBubbleCount << " unambiguous bubbles." << endl;
            cout << "Bubble index range for this oriented read is [" <<
                minBubbleIndex << "," << maxBubbleIndex << "]." << endl;
        }

        if(orientedReadInfo.unambiguousBubbleCount == 0) {
            break;
        }

        // Assign phase +1 in this phasing group to this starting read for this phasing component.
        SHASTA_ASSERT(orientedRead.phase == 0);
        SHASTA_ASSERT(orientedRead.phasingComponent ==  invalid<uint64_t>);
        orientedRead.phase = +1;
        orientedRead.phasingComponent = phasingComponent;

        // Assign to all unambiguous bubbles of this oriented read a phase consistent with it.
        for(auto it=indexByOrientedReadId().find(orientedRead.id);
            it!=indexByOrientedReadId().end() and it->orientedReadId == orientedRead.id; ++it) {
            const PhasingTableEntry& phasingTableEntry = *it;
            Bubble& bubble = bubbles[phasingTableEntry.bubbleIndex];
            SHASTA_ASSERT(bubble.phase == 0);
            SHASTA_ASSERT(bubble.phasingComponent == invalid<uint64_t>);

            // Skip it if it is ambiguous.
            if(phasingTableEntry.discreteRelativePhase == 0) {
                continue;
            }

            // Set the phase of this bubble to a phase consistent with the +1 phase
            // we assigned to the starting oriented read.
            bubble.phase = phasingTableEntry.discreteRelativePhase;
            bubble.phasingComponent = phasingComponent;

            // Update the OrientedReadTable to reflect the fact that this bubble was just phased.
            for(auto it=indexByPositionInBubbleChain().find(bubble.positionInBubbleChain);
                it!=indexByPositionInBubbleChain().end() and it->positionInBubbleChain == bubble.positionInBubbleChain; ++it) {
                const PhasingTableEntry& phasingTableEntry = *it;
                if(phasingTableEntry.discreteRelativePhase == 0) {
                    continue;
                }

                auto jt = orientedReadTable.get<0>().find(phasingTableEntry.orientedReadIndex);
                SHASTA_ASSERT(jt != orientedReadTable.get<0>().end());
                OrientedReadInfo info = *jt;
                info.phasedUnambiguousBubbleCount++;
                orientedReadTable.get<0>().replace(jt, info);
            }
        }

        // Remove the starting oriented read from the orientedReadTable.
        orientedReadTable.get<1>().erase(it);



        // The inner loop phases one oriented read at a time, adding it to the current
        // phasing component.
        while(not orientedReadTable.empty()) {

            // Find the oriented read with the most phased bubbles.
            const auto it = orientedReadTable.get<2>().begin();
            const OrientedReadInfo& orientedReadInfo = *it;
            OrientedRead& orientedRead = orientedReads[orientedReadInfo.orientedReadIndex];

            const uint64_t minPositionInBubbleChain = orientedRead.minPositionInBubbleChain;
            const uint64_t maxPositionInBubbleChain = orientedRead.maxPositionInBubbleChain;
            const uint64_t minBubbleIndex = bubblesMap[minPositionInBubbleChain];
            const uint64_t maxBubbleIndex = bubblesMap[maxPositionInBubbleChain];

            if(orientedReadInfo.phasedUnambiguousBubbleCount == 0) {
                // Finish this phasing component.
                break;
            }

            if(debug) {
                cout << "Adding to phasing group " << orientedRead.id <<
                    ", index " << orientedReadInfo.orientedReadIndex <<
                    " with " << orientedReadInfo.unambiguousBubbleCount << " unambiguous bubbles," <<
                    " of which " << orientedReadInfo.phasedUnambiguousBubbleCount << " already phased ." << endl;
                cout << "Bubble index range for this oriented read is [" <<
                    minBubbleIndex << "," << maxBubbleIndex << "]." << endl;
            }

            // Use the bubbles that are already phased to assign a phase to this oriented read.
            uint64_t plusCount = 0;
            uint64_t minusCount = 0;
            for(auto it=indexByOrientedReadId().find(orientedRead.id);
                it!=indexByOrientedReadId().end() and it->orientedReadId == orientedRead.id; ++it) {
                const PhasingTableEntry& phasingTableEntry = *it;
                if(phasingTableEntry.discreteRelativePhase == 0) {
                    continue;
                }
                Bubble& bubble = bubbles[phasingTableEntry.bubbleIndex];
                if(bubble.phase == 0) {
                    continue;
                }
                int64_t phase;
                if(phasingTableEntry.discreteRelativePhase == +1) {
                    phase = bubble.phase;
                } else {
                    phase = - bubble.phase;
                }
                if(phase == +1) {
                    ++plusCount;
                } else if(phase == -1) {
                    ++minusCount;
                }
            }

            SHASTA_ASSERT(plusCount + minusCount == orientedReadInfo.phasedUnambiguousBubbleCount);

            // Phase this oriented read in this phasing component.
            SHASTA_ASSERT(orientedRead.phase == 0);
            SHASTA_ASSERT(orientedRead.phasingComponent ==  invalid<uint64_t>);
            orientedRead.phase = (plusCount >= minusCount) ? +1 : -1;
            orientedRead.phasingComponent = phasingComponent;

            // Assign to all unambiguous bubbles of this oriented read
            // that are not already phased a phase consistent with it.
            for(auto it=indexByOrientedReadId().find(orientedRead.id);
                it!=indexByOrientedReadId().end() and it->orientedReadId == orientedRead.id; ++it) {
                const PhasingTableEntry& phasingTableEntry = *it;

                // Skip it if it is ambiguous.
                if(phasingTableEntry.discreteRelativePhase == 0) {
                    continue;
                }
                Bubble& bubble = bubbles[phasingTableEntry.bubbleIndex];

                // If already phased, skip it.
                if(bubble.phase != 0) {
                    continue;
                }

                // Phase this bubble to a phase consistent with this oriented read.
                bubble.phase = (phasingTableEntry.discreteRelativePhase == +1) ? orientedRead.phase : -orientedRead.phase;
                bubble.phasingComponent = phasingComponent;

                // Update the OrientedReadTable to reflect the fact that this bubble was just phased.
                for(auto it=indexByPositionInBubbleChain().find(bubble.positionInBubbleChain);
                    it!=indexByPositionInBubbleChain().end() and it->positionInBubbleChain == bubble.positionInBubbleChain; ++it) {
                    const PhasingTableEntry& phasingTableEntry = *it;
                    if(phasingTableEntry.discreteRelativePhase == 0) {
                        continue;
                    }

                    auto jt = orientedReadTable.get<0>().find(phasingTableEntry.orientedReadIndex);
                    if(jt == orientedReadTable.get<0>().end()) {
                        continue;
                    }
                    OrientedReadInfo info = *jt;
                    info.phasedUnambiguousBubbleCount++;
                    orientedReadTable.get<0>().replace(jt, info);
                }
            }

            // Remove the oriented read from the orientedReadTable.
            orientedReadTable.get<2>().erase(it);
        }
    }
}



double PhasingTable::bubbleErrorRate(uint64_t positionInBubbleChain) const
{
    // Must be called for a diploid bubble
    auto it = bubblesMap.find(positionInBubbleChain);
    SHASTA_ASSERT(it != bubblesMap.end());
    const Bubble& bubble = bubbles[it->second];

    if(bubble.phase == 0) {
        return 1.;
    }

    // This bubble is diploid and phased.
    uint64_t consistent;
    uint64_t inconsistent;
    tie(consistent, inconsistent) = countConsistentEntriesForBubble(positionInBubbleChain);
    return double(inconsistent) / double(consistent+ inconsistent);
}



// Use the phases stored in the Bubbles to consruct the PhasedComponents.
// The PhasedComponents must be non-overlapping and sorted by position.
void PhasingTable::constructPhasedComponents(bool debug)
{
    phasedComponents.clear();

    // Create an initial version of PhasedComponents without
    // worrying about ordering by position and about overlap between PhasedComponents.
    for(const Bubble& bubble: bubbles) {
        if(bubble.phase == 0) {
            continue;
        }
        const uint64_t phasedComponentId = bubble.phasingComponent;
        if(phasedComponentId >= phasedComponents.size()) {
            for(uint64_t i=phasedComponents.size(); i<=phasedComponentId; i++) {
                phasedComponents.push_back(make_shared<PhasedComponent>());
            }
        }
        phasedComponents[phasedComponentId]->push_back({bubble.positionInBubbleChain, bubble.phase});
    }

    if(debug) {
        uint64_t totalPhasedBubbleCount = 0;
        for(const auto& phasedComponent: phasedComponents) {
            totalPhasedBubbleCount += phasedComponent->size();
        }
        cout << "Created " << phasedComponents.size() << " initial phased components "
            "with a total " << totalPhasedBubbleCount << " phased diploid bubbles." << endl;
    }



    // If there is more than one PhasedComponent, we have to eliminate overlaps.
    // We do this by removing bubbles from overlapping PhasedComponents, giving
    // priority to larger PhasedComponents.
    if(phasedComponents.size() > 1) {

        if(debug) {
            cout << "More than one phased components found. Removing overlaps." << endl;
        }

        // Sort the phased components by decreasing size.
        class SortHelper {
        public:
            bool operator()(
                const shared_ptr<PhasedComponent>& p0,
                const shared_ptr<PhasedComponent>& p1
                ) const
            {
                return p0->size() > p1->size();
            }
        };
        sort(phasedComponents.begin(), phasedComponents.end(), SortHelper());

        for(const auto& phasedComponent: phasedComponents) {
            phasedComponent->computePositionRange();
        }

        // Process the PhasedComponents in order of decreasing size.
        vector< pair<uint64_t, uint64_t> > forbiddenRanges; // (min, max)
        for(auto& phasedComponent: phasedComponents) {

            // See if it overlaps any of the forbidden ranges.
            bool overlaps = false;
            for(const auto& forbiddenRange: forbiddenRanges) {
                const bool disjointLeft  = phasedComponent->maxPositionInBubbleChain < forbiddenRange.first;
                const bool disjointRight = phasedComponent->minPositionInBubbleChain > forbiddenRange.second;
                if(not(disjointLeft or disjointRight)) {
                    overlaps = true;
                    break;
                }
            }

            if(debug) {
                cout << "Phased component at " << phasedComponent->minPositionInBubbleChain << " " <<
                    phasedComponent->maxPositionInBubbleChain;
                if(overlaps) {
                    cout << " overlaps a previous phased component." << endl;
                } else {
                    cout << " has no overlaps with previous phased components." << endl;
                }
            }

            if(not overlaps) {
                forbiddenRanges.push_back(
                    {phasedComponent->minPositionInBubbleChain, phasedComponent->maxPositionInBubbleChain});
                continue;
            }



            // This PhasedComponent overlaps a forbiddenRange.
            // We need to remove the offending bubbles.
            shared_ptr<PhasedComponent> newPhasedComponent = make_shared<PhasedComponent>();
            for(const auto& p: *phasedComponent) {
                const uint64_t positionInBubbleChain = p.first;

                // See if this bubble overlaps any forbidden ranges.
                bool overlaps = false;
                for(const auto& forbiddenRange: forbiddenRanges) {
                    if( positionInBubbleChain >= forbiddenRange.first and
                        positionInBubbleChain <= forbiddenRange.second) {
                        overlaps = true;
                        break;
                    }
                }

                // Only keep it if there is no overlap.
                if(not overlaps) {
                    newPhasedComponent->push_back(p);
                }

            }

            // Replace this phased component with the new one.
            phasedComponent = newPhasedComponent;
            phasedComponent->computePositionRange();
            forbiddenRanges.push_back({phasedComponent->minPositionInBubbleChain, phasedComponent->maxPositionInBubbleChain});

            if(debug) {
                cout << "After removing overlap, this phased component has " << phasedComponent->size() <<
                    " diploid bubbles and position range " << phasedComponent->minPositionInBubbleChain << " " <<
                    phasedComponent->maxPositionInBubbleChain << endl;
            }
        }
    }



    // This could have created empty PhasedComponents.
    // Remove them if they are present.
    {
        vector< shared_ptr<PhasedComponent> > nonEmptyPhasedComponents;
        for(const shared_ptr<PhasedComponent>& phasedComponent: phasedComponents) {
            if(not phasedComponent->empty()) {
                nonEmptyPhasedComponents.push_back(phasedComponent);
            } else {
                if(debug) {
                    cout << "Removing empty phased component." << endl;
                }
            }
        }
        if(nonEmptyPhasedComponents.size() != phasedComponents.size()) {
            phasedComponents.swap(nonEmptyPhasedComponents);
        }
    }



    // Compute the position ranges.
    for(const auto& phasedComponent: phasedComponents) {
        phasedComponent->computePositionRange();
    }

    // Sort the phased components in order of increasing position.
    class SortHelper {
    public:
        bool operator()(
            const shared_ptr<PhasedComponent>& p0,
            const shared_ptr<PhasedComponent>& p1
            ) const
        {
            return p0->minPositionInBubbleChain < p1->minPositionInBubbleChain;
        }
    };
    sort(phasedComponents.begin(), phasedComponents.end(), SortHelper());

    if(debug) {
        cout << phasedComponents.size() << " phased components:" << endl;
        for(const auto& phasedComponent: phasedComponents) {
            cout  << phasedComponent->size() << " diploid bubbles at positions " <<
                phasedComponent->minPositionInBubbleChain << "..." <<
                phasedComponent->maxPositionInBubbleChain << " in bubble chain." << endl;

        }
        // phasingGraph.writeGraphviz("PhasingGraph.dot");
    }
}

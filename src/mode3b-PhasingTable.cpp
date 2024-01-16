// Shasta.
#include "mode3b-PhasingTable.hpp"
#include "mode3b-CompressedPathGraph1B.hpp"
#include "Assembler.hpp"
#include "html.hpp"
#include "orderPairs.hpp"
#include "PngImage.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <filesystem>
#include "tuple.hpp"


void CompressedPathGraph1B::writeBubbleChainsPhasingTables(
    const string& fileNamePrefix,
    double phaseErrorThreshold) const
{
    const CompressedPathGraph1B& cGraph = *this;

    const string directoryName = fileNamePrefix + "-PhasingTables";
    if(not std::filesystem::create_directory(directoryName)) {
        throw runtime_error("Could not create directory " + directoryName);
    }


    // Loop over all BubbleChains.
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1B) {
        const CompressedPathGraph1BEdge& edge = cGraph[ce];
        const BubbleChain& bubbleChain = edge;

        // Create the phasing table for this bubble chain.
        PhasingTable phasingTable(bubbleChain, assembler.markerGraph, phaseErrorThreshold);

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

        phasingTable.simpleIterativePhasing2();
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
    const MarkerGraph& markerGraph,
    double phaseErrorThreshold)
{
    fill(bubbleChain, markerGraph, phaseErrorThreshold);
    gatherOrientedReads();
    gatherBubbles();
    fillIndexes();
}



void PhasingTable::fill(
    const BubbleChain& bubbleChain,
    const MarkerGraph& markerGraph,
    double phaseErrorThreshold)
{
    clear();

    // Loop over the bubbles in this bubble chain.
    for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
        const mode3b::Bubble& bubble = bubbleChain[positionInBubbleChain];

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
                const MarkerGraphEdgeId markerGraphEdgeId = chain[i];

                // Loop over MarkerIntervals of this marker graph edge.
                const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
                for(const MarkerInterval& markerInterval: markerIntervals) {
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
    writeBubblesCsv(fileNamePrefix);
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



void PhasingTable::writeBubblesCsv(const string& fileNamePrefix) const
{
    ofstream csv(fileNamePrefix + "-Bubbles.csv");
    csv << "Position in bubble chain,Bubble index,\n";

    for(uint64_t i=0; i<bubbles.size(); i++) {
        csv << bubbles[i].positionInBubbleChain << ",";
        csv << i << ",";
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



#if 0
void PhasingTable::writeHtml(const string& fileNamePrefix) const
{
    const auto& indexByBoth = get<0>();

    ofstream html(fileNamePrefix + ".html");
    writeHtmlBegin(html, "Phasing Table");

    html << "<body><canvas id='canvas1' width='" << bubbleCount() <<
        "' height='" << orientedReadCount() << "' style='border:1px solid #000000;'></canvas>\n";

    html <<
        "<script>\n"
        "var canvas = document.getElementById('canvas1');\n"
        "var c = canvas.getContext('2d');\n"
        "c.fillStyle = 'white';\n"
        "c.fillRect(0, 0, " << bubbleCount() << ", " << orientedReadCount() << ");\n"
        "c.fillStyle = 'hsl(0, 100%, 50%)';\n";
    for(const PhasingTableEntry& entry: indexByBoth) {
        const uint64_t h = hue(phase(entry));
        html << "c.fillStyle = 'hsl(" << h << ", 100%, 50%)';\n";
        html << "c.fillRect(" << entry.bubbleIndex << ", " << entry.orientedReadIndex << ", 1, 1);\n";
    }
    html << "</script>";

    html << "</body>";
    writeHtmlEnd(html);
}
#endif




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
                // Ambiguous: yellow
                r = 255;
                g = 255;
                b = 0;
                break;
            case +1:
                // In-phase: green.
                r = 0;
                g = 255;
                b = 0;
                break;
            case -1:
                // Out-of-phase: red.
                r = 255;
                g = 0;
                b = 0;
                break;
            default:
                SHASTA_ASSERT(0);
            }

        } else if(coloringMethod == ColoringMethod::byRelativePhase) {

            // Compute (r, g, b) values that give:
            // - Green if relativePhase is 1 (in-phase).
            // - Red if relativePhase is -1 (out-of-phase).
            if(entry.relativePhase >= 0.) {
                r = int(std::round((1. - entry.relativePhase) * 255.));
                g = 255;
                b = 0;
            } else {
                r = 255;
                g = int(std::round((1. +  entry.relativePhase) * 255.));
                b = 0;
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



#if 0
void PhasingTable::flipSweep()
{
    const auto& indexByOrientedReadId = get<1>();

    // Loop over oriented reads.
    vector<uint64_t> plusBubbles;
    vector<uint64_t> minusBubbles;
    for(uint64_t i=0; i<orientedReadCount(); i++) {
        const OrientedReadId orientedReadId = orientedReadInfos[i].orientedReadId;

        // Gather the bubbles where this oriented read appears with phase +1 or -1
        // (with tolerance equal to phaseError).
        plusBubbles.clear();
        minusBubbles.clear();
        for(auto it=indexByOrientedReadId.find(orientedReadId);
            it!=indexByOrientedReadId.end() and it->orientedReadId == orientedReadId; ++it) {
            const PhasingTableEntry& phasingTableEntry = *it;

            // Compute the phase for this entry, taking into account possible flip of the Bubble.
            const double p = phase(phasingTableEntry);

            if(p >= phaseThresholdPlus) {
                plusBubbles.push_back(phasingTableEntry.bubbleIndex);
            } else if(p <= phaseThresholdMinus) {
                minusBubbles.push_back(phasingTableEntry.bubbleIndex);
            }
        }


        // If there are more plusBubbles than minusBubbles, flip the minusBubbles.
        // If there are more minusBubbles than plusBubbles, flip the plusBubbles.
        if(plusBubbles.size() == minusBubbles.size()) {
            continue;
        }
        const vector<uint64_t>& bubblesToFlip =
            (plusBubbles.size() > minusBubbles.size()) ? minusBubbles : plusBubbles;
        for(const uint64_t bubbleIndex: bubblesToFlip) {
            Bubble& bubble = bubbles[bubbleIndex];
            bubble.flip = not bubble.flip;
        }
    }
}



// Count the number of PhaseTableEntries for a given oriented read
// with phase +1 or -1, allowing a phaseError discrepancy up to phase Error.
// That is, countPlus is the number of oriented read entries with phase >= 1 - phaseError,
// and countMinus is the number of oriented read entries with phase <= -1 + phaseError.
void PhasingTable::count(
    OrientedReadId orientedReadId,
    uint64_t& countPlus,
    uint64_t& countMinus) const
{
    countPlus = 0;
    countMinus = 0;

    const auto& indexByOrientedReadId = get<1>();

    for(auto it=indexByOrientedReadId.find(orientedReadId);
        it!=indexByOrientedReadId.end() and it->orientedReadId == orientedReadId; ++it) {
        const PhasingTableEntry& phasingTableEntry = *it;

        // Compute the phase for this entry, taking into account possible flip of the Bubble.
        const double p = phase(phasingTableEntry);

        if(p >= phaseThresholdPlus) {
            ++countPlus;
        } else if(p <= phaseThresholdMinus) {
            ++countMinus;
        }
    }
}



uint64_t PhasingTable::discordantCount() const
{
    uint64_t n = 0;
    for(uint64_t i=0; i<orientedReadCount(); i++) {
        const OrientedReadId orientedReadId = orientedReadInfos[i].orientedReadId;
        uint64_t countPlus;
        uint64_t countMinus;
        count(orientedReadId, countPlus, countMinus);
        n += min(countPlus, countMinus);
    }
    return n;
}
#endif



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



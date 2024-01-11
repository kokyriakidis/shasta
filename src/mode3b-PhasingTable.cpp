// Shasta.
#include "mode3b-PhasingTable.hpp"
#include "mode3b-CompressedPathGraph1B.hpp"
#include "Assembler.hpp"
#include "html.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <filesystem>
#include "tuple.hpp"


void CompressedPathGraph1B::writeBubbleChainsPhasingTables(const string& fileNamePrefix) const
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
        PhasingTable phasingTable(bubbleChain, assembler.markerGraph);

        cout << "Phasing table for " << bubbleChainStringId(ce) <<
            " has " << phasingTable.bubbleCount() << " bubbles and " <<
            phasingTable.orientedReadCount() << " oriented reads." << endl;

        // Write to csv.
        phasingTable.write(directoryName + "/" + bubbleChainStringId(ce));

    }
}



PhasingTable::PhasingTable(
    const BubbleChain& bubbleChain,
    const MarkerGraph& markerGraph)
{
    fill(bubbleChain, markerGraph);
    gatherBubbles();
    gatherOrientedReads();
    fillIndexes();
}



void PhasingTable::fill(
    const BubbleChain& bubbleChain,
    const MarkerGraph& markerGraph)
{
    clear();
    const auto& indexByBoth = get<0>();

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

                    // Increment the PhasingTableEntry for this OrientedReadId and positionInBubbleChain.
                    auto it = indexByBoth.find(make_tuple(orientedReadId, positionInBubbleChain));
                    if(it == indexByBoth.end()) {
                        tie(it, ignore) = insert({orientedReadId, positionInBubbleChain});
                    }

                    // We can do a const_cast because we only update the frequency,
                    // which does not participate in any index.
                    PhasingTableEntry& entry = const_cast<PhasingTableEntry&>(*it);
                    ++entry.frequency[chainIndexInBubble];
                }
            }

        }

    }
}



void PhasingTable::write(const string& fileNamePrefix) const
{
    writeCsv(fileNamePrefix);
    writeHtml(fileNamePrefix);
}



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

    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        csv << info.orientedReadId << ",";
        csv << info.minPositionInBubbleChain << ",";
        csv << info.maxPositionInBubbleChain << ",";
        csv << i << ",";
        csv << bubblesMap.find(info.minPositionInBubbleChain)->second << ",";
        csv << bubblesMap.find(info.maxPositionInBubbleChain)->second << ",";
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
    const auto& indexByOrientedReadId = get<1>();

    ofstream csv(fileNamePrefix + "-Details.csv");

    csv << "Position in bubble chain,OrientedReadId,Bubble index,Oriented read index,Frequency0,Frequency1,\n";

    for(const OrientedReadInfo& info: orientedReadInfos) {
        const OrientedReadId orientedReadId = info.orientedReadId;
        for(auto it=indexByOrientedReadId.find(orientedReadId);
            it!=indexByOrientedReadId.end() and it->orientedReadId == orientedReadId; ++it) {
            const PhasingTableEntry& phasingTableEntry = *it;
            phasingTableEntry.writeCsv(csv);
            csv << "\n";
        }
    }
}



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
        const double fraction = entry.fraction0();
        const double angle = 360. - 120. * fraction;
        uint64_t iAngle = uint64_t(std::round(angle));
        iAngle = min(iAngle, 360UL);
        html << "c.fillStyle = 'hsl(" << iAngle << ", 100%, 50%)';\n";
        html << "c.fillRect(" << entry.bubbleIndex << ", " << entry.orientedReadIndex << ", 1, 1);\n";
    }
    html << "</script>";

    html << "</body>";
    writeHtmlEnd(html);
}



void PhasingTable::gatherBubbles()
{
    const auto& indexByBoth = get<0>();

    std::set<uint64_t> positionsInBubbleChain;
    for(const PhasingTableEntry& phasingTableEntry: indexByBoth) {
        positionsInBubbleChain.insert(phasingTableEntry.positionInBubbleChain);
    }

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



void PhasingTable::gatherOrientedReads()
{
    const auto& indexByBoth = get<0>();
    const auto& indexByOrientedReadId = get<1>();

    std::set<OrientedReadId> orientedReadIds;
    for(const PhasingTableEntry& phasingTableEntry: indexByBoth) {
        orientedReadIds.insert(phasingTableEntry.orientedReadId);
    }

    orientedReadInfos.clear();
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        OrientedReadInfo info;
        info.orientedReadId = orientedReadId;
        orientedReadInfos.push_back(info);
    }

    for(OrientedReadInfo& info: orientedReadInfos) {
        info.minPositionInBubbleChain = std::numeric_limits<uint64_t>::max();
        info.maxPositionInBubbleChain = 0;
        for(auto it=indexByOrientedReadId.find(info.orientedReadId);
            it!=indexByOrientedReadId.end() and it->orientedReadId == info.orientedReadId; ++it) {
            const uint64_t positionInBubbleChain = it->positionInBubbleChain;
            info.minPositionInBubbleChain = min(info.minPositionInBubbleChain, positionInBubbleChain);
            info.maxPositionInBubbleChain = max(info.maxPositionInBubbleChain, positionInBubbleChain);
        }
    }

    // We want the reads ordered by position.
    vector< pair<uint64_t, uint64_t> > orientedReadIdsTable; // (index, minPosition + maxPosition)
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& info = orientedReadInfos[i];
        orientedReadIdsTable.push_back({i, info.minPositionInBubbleChain + info.maxPositionInBubbleChain});
    }
    sort(orientedReadIdsTable.begin(), orientedReadIdsTable.end(),
        OrderPairsBySecondOnly<uint64_t, uint64_t>());
    vector<OrientedReadInfo> sortedInfos;
    for(const auto& p: orientedReadIdsTable) {
        sortedInfos.push_back(orientedReadInfos[p.first]);
    }
    orientedReadInfos.swap(sortedInfos);

    // Fill in the orientedReadIdsMap map.
    orientedReadIdsMap.clear();
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        orientedReadIdsMap.insert({orientedReadInfos[i].orientedReadId, i});
    }
}



void PhasingTable::fillIndexes()
{
    const auto& indexByBoth = get<0>();

    for(const PhasingTableEntry& phasingTableEntry: indexByBoth) {
        PhasingTableEntry& entry = const_cast<PhasingTableEntry&>(phasingTableEntry); // OK because we are not modifying the keys
        entry.bubbleIndex = bubblesMap[entry.positionInBubbleChain];
        entry.orientedReadIndex = orientedReadIdsMap[entry.orientedReadId];
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
}

// Shasta.
#include "mode3-TangleGraph.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <algorithm.hpp>
#include <fstream.hpp>
#include <queue>



TangleGraph::TangleGraph(
    bool debug,
    uint64_t tangleId,
    const Anchors& anchors,
    const vector<AnchorId>& entranceAnchors,
    const vector<AnchorId>& exitAnchors,
    bool bidirectional,
    double maxLoss,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold) :
    debug(debug),
    tangleId(tangleId),
    anchors(anchors),
    bidirectional(bidirectional)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minVertexCoverage = 5;
    const uint64_t anchorCoverageThreshold = 3;


    if(debug) {
        cout << "Creating a tangle graph for tangle " << tangleId <<
            " with " << entranceAnchors.size() <<
            " entrances and " << exitAnchors.size() << " exits." << endl;
        cout << "Entrances:";
        for(const AnchorId anchorId: entranceAnchors) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;
        cout << "Exits:";
        for(const AnchorId anchorId: exitAnchors) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;
    }

    constructEntrances(entranceAnchors);
    constructExits(exitAnchors);
    gatherOrientedReads();

    if(not createVertices(anchorCoverageThreshold, minVertexCoverage)) {
        failure = true;
        return;
    }
    createEdges();
    if(debug) {
        cout << "The initial tangle graph has " << num_vertices(*this) <<
            " vertices and " << num_edges(*this) << " edges." << endl;
    }
    // writeGraphviz("Initial");

    removeWeakEdges(maxLoss);
    removeCrossEdges(lowCoverageThreshold, highCoverageThreshold);
    // writeGraphviz("Cleanedup");

    if(not removeUnreachable()) {
        failure = true;
        if(debug) {
            cout << "Reachability failure." << endl;
        }
        return;
    }

    if(debug) {
        cout << "The final tangle graph has " << num_vertices(*this) <<
            " vertices and " << num_edges(*this) << " edges." << endl;
        writeGraphviz("Final");
    }

}



void TangleGraph::constructEntrances(const vector<AnchorId>& entranceAnchors)
{
    for(const AnchorId anchorId: entranceAnchors) {
        entrances.push_back(Entrance(anchorId));
    }

#if 0
    // Initialize the anchorMarkerIntervals of the entrances.
    for(const AnchorId anchorId: entranceAnchors) {
        entrances.push_back(Entrance(anchorId, anchors[anchorId]));
    }

    // If an OrientedReadId appears in more than one entrance,
    // we want to remove the corresponding AnchorMarkerIntervals
    // from all entrances.

    // Find the duplicate OrientedReadIds.
    vector<OrientedReadId> orientedReadIds;
    for(const Entrance& entrance: entrances) {
        for(const AnchorMarkerInterval& anchorMarkerInterval: entrance.anchorMarkerIntervals) {
            orientedReadIds.push_back(anchorMarkerInterval.orientedReadId);
        }
    }
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(orientedReadIds, count, 2UL);

    if(debug) {
        if(orientedReadIds.empty()) {
            cout << "The entrances contain no duplicate oriented reads." << endl;
        } else {
            cout << "The following oriented reads appear in more than once entrance "
                " will be neglected during read following from the entrances:";
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                cout << " " << orientedReadId;
            }
            cout << endl;
        }
    }

    // If we found any duplicate OrientedReadIds, remove them from all entrances.
    if(not orientedReadIds.empty()) {
        for(Entrance& entrance: entrances) {
            vector<AnchorMarkerInterval> newAnchorMarkerIntervals;
            for(const AnchorMarkerInterval& anchorMarkerInterval: entrance.anchorMarkerIntervals) {
                if(not binary_search(orientedReadIds.begin(), orientedReadIds.end(), anchorMarkerInterval.orientedReadId)) {
                    newAnchorMarkerIntervals.push_back(anchorMarkerInterval);
                }
            }
            entrance.anchorMarkerIntervals.swap(newAnchorMarkerIntervals);
        }
    }



    // Read following for each entrance.
    // This fills in the journeyAnchorIds of each entrance.
    for(Entrance& entrance: entrances) {
        entrance.readFollowing(debug, anchors, bidirectional);
    }



    // Find the AnchorIds that were found in more than one Entrance.
    vector<AnchorId> duplicateAnchorIds;
    for(Entrance& entrance: entrances) {
        copy(entrance.journeyAnchorIds.begin(), entrance.journeyAnchorIds.end(),
            back_inserter(duplicateAnchorIds));
    }
    deduplicateAndCountWithThreshold(duplicateAnchorIds, count, 2UL);

    if(debug) {
        cout << duplicateAnchorIds.size() <<
            " anchors were found by read following from more than one entrance." << endl;
    }


    // Now we can fill in the uniqueJourneyAnchorIds for each Entrance.
    for(Entrance& entrance: entrances) {
        std::set_difference(
            entrance.journeyAnchorIds.begin(), entrance.journeyAnchorIds.end(),
            duplicateAnchorIds.begin(), duplicateAnchorIds.end(),
            back_inserter(entrance.uniqueJourneyAnchorIds));

        if(debug) {
            cout << "Read following for entrance " << anchorIdToString(entrance.anchorId) <<
                " found " << entrance.uniqueJourneyAnchorIds.size() << " anchors unique to this entrance." << endl;
        }
    }
#endif
}



void TangleGraph::constructExits(const vector<AnchorId>& exitAnchors)
{
    for(const AnchorId anchorId: exitAnchors) {
        exits.push_back(Exit(anchorId));
    }

    #if 0
    // Initialize the anchorMarkerIntervals of the entrances.
    for(const AnchorId anchorId: exitAnchors) {
        exits.push_back(Exit(anchorId, anchors[anchorId]));
    }

    // If an OrientedReadId appears in more than one exit,
    // we want to remove the corresponding AnchorMarkerIntervals
    // from all exits.

    // Find the duplicate OrientedReadIds.
    vector<OrientedReadId> orientedReadIds;
    for(const Exit& exit: exits) {
        for(const AnchorMarkerInterval& anchorMarkerInterval: exit.anchorMarkerIntervals) {
            orientedReadIds.push_back(anchorMarkerInterval.orientedReadId);
        }
    }
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(orientedReadIds, count, 2UL);

    if(debug) {
        if(orientedReadIds.empty()) {
            cout << "The exits contain no duplicate oriented reads." << endl;
        } else {
            cout << "The following oriented reads appear in more than once exits "
                " and will be neglected during read following from the exits:";
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                cout << " " << orientedReadId;
            }
            cout << endl;
        }
    }

    // If we found any duplicate OrientedReadIds, remove them from all exits.
    if(not orientedReadIds.empty()) {
        for(Exit& exit: exits) {
            vector<AnchorMarkerInterval> newAnchorMarkerIntervals;
            for(const AnchorMarkerInterval& anchorMarkerInterval: exit.anchorMarkerIntervals) {
                if(not binary_search(orientedReadIds.begin(), orientedReadIds.end(), anchorMarkerInterval.orientedReadId)) {
                    newAnchorMarkerIntervals.push_back(anchorMarkerInterval);
                }
            }
            exit.anchorMarkerIntervals.swap(newAnchorMarkerIntervals);
        }
    }



    // Read following for each exit.
    // This fills in the journeyAnchorIds of each exit.
    for(Exit& exit: exits) {
        exit.readFollowing(debug, anchors, bidirectional);
    }



    // Find the AnchorIds that were found in more than one Exit.
    vector<AnchorId> duplicateAnchorIds;
    for(Exit& exit: exits) {
        copy(exit.journeyAnchorIds.begin(), exit.journeyAnchorIds.end(),
            back_inserter(duplicateAnchorIds));
    }
    deduplicateAndCountWithThreshold(duplicateAnchorIds, count, 2UL);

    if(debug) {
        cout << duplicateAnchorIds.size() <<
            " anchors were found by read following from more than one exit." << endl;
    }


    // Now we can fill in the uniqueJourneyAnchorIds for each Entrance.
    for(Exit& exit: exits) {
        std::set_difference(
            exit.journeyAnchorIds.begin(), exit.journeyAnchorIds.end(),
            duplicateAnchorIds.begin(), duplicateAnchorIds.end(),
            back_inserter(exit.uniqueJourneyAnchorIds));

        if(debug) {
            cout << "Read following for exit " << anchorIdToString(exit.anchorId) <<
                " found " << exit.uniqueJourneyAnchorIds.size() << " anchors unique to this exit." << endl;
        }
    }
#endif
}



TangleGraph::EntranceOrExit::EntranceOrExit(
    AnchorId anchorId) :
    anchorId(anchorId)
{
}


#if 0
// This fills in the journeyAnchorIds.
void TangleGraph::Entrance::readFollowing(
    bool debug, const Anchors& anchors, bool bidirectional)
{
    for(const AnchorMarkerInterval& anchorMarkerInterval: anchorMarkerIntervals) {
        const OrientedReadId orientedReadId = anchorMarkerInterval.orientedReadId;
        const auto journey = anchors.journeys[orientedReadId.getValue()];

        const uint64_t begin = (bidirectional ? 0 : anchorMarkerInterval.positionInJourney);
        const uint64_t end = journey.size();
        for(uint64_t position = begin; position != end; position++) {
            journeyAnchorIds.push_back(journey[position]);
        }
    }

    deduplicate(journeyAnchorIds);

    if(debug) {
        cout << "Read following for entrance " << anchorIdToString(anchorId) <<
            " found " << journeyAnchorIds.size() << " anchors after deduplication." << endl;
    }

}



// This fills in the journeyAnchorIds.
void TangleGraph::Exit::readFollowing(
    bool debug, const Anchors& anchors, bool bidirectional)
{
    for(const AnchorMarkerInterval& anchorMarkerInterval: anchorMarkerIntervals) {
        const OrientedReadId orientedReadId = anchorMarkerInterval.orientedReadId;
        const auto journey = anchors.journeys[orientedReadId.getValue()];

        const uint64_t begin = 0;
        const uint64_t end = (bidirectional ? journey.size() : anchorMarkerInterval.positionInJourney + 1);
        for(uint64_t position = begin; position != end; position++) {
            journeyAnchorIds.push_back(journey[position]);
        }
    }

    deduplicate(journeyAnchorIds);

    if(debug) {
        cout << "Read following for exit " << anchorIdToString(anchorId) <<
            " found " << journeyAnchorIds.size() << " anchors after deduplication." << endl;
    }

}
#endif



void TangleGraph::gatherOrientedReads()
{

    // Gather the OrientedReadIds that appear in one entrance and no more than one.
    vector<OrientedReadId> entranceOrientedReadIds;
    for(const Entrance& entrance: entrances) {
        const Anchor anchor = anchors[entrance.anchorId];
        for(const auto& anchorMarkerInterval: anchor) {
            entranceOrientedReadIds.push_back(anchorMarkerInterval.orientedReadId);
        }
    }
    deduplicateAndCountAndKeepUnique(entranceOrientedReadIds);
    // The entranceOrientedReadIds are now sorted.

    // Gather the OrientedReadIds that appear in one exit and no more than one.
    vector<OrientedReadId> exitOrientedReadIds;
    for(const Exit& exit: exits) {
        const Anchor anchor = anchors[exit.anchorId];
        for(const auto& anchorMarkerInterval: anchor) {
            exitOrientedReadIds.push_back(anchorMarkerInterval.orientedReadId);
        }
    }
    deduplicateAndCountAndKeepUnique(exitOrientedReadIds);
    // The exitOrientedReadIds are now sorted.

    // We will work with the union set of entranceOrientedReadIds and exitOrientedReadIds.
    // This is also stored sorted.
    vector<OrientedReadId> orientedReadIds;
    std::set_union(
        entranceOrientedReadIds.begin(), entranceOrientedReadIds.end(),
        exitOrientedReadIds.begin(), exitOrientedReadIds.end(),
        back_inserter(orientedReadIds));

    // Initialize the OrientedReadInfos.
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        orientedReadInfos.push_back(OrientedReadInfo(orientedReadId));
    }
    // We can now use getOrientedReadInfo.

    // Fill in the OrientedReadInfos.
    for(uint64_t entranceIndex=0; entranceIndex<entrances.size(); entranceIndex++) {
        const Entrance& entrance = entrances[entranceIndex];
        const Anchor anchor = anchors[entrance.anchorId];
        for(const auto& anchorMarkerInterval: anchor) {
            OrientedReadInfo* orientedReadInfo = getOrientedReadInfo(anchorMarkerInterval.orientedReadId);
            if(orientedReadInfo) {
                orientedReadInfo->entranceIndex = entranceIndex;
                orientedReadInfo->entrancePositionInJourney = anchorMarkerInterval.positionInJourney;
            }
        }
    }
    for(uint64_t exitIndex=0; exitIndex<exits.size(); exitIndex++) {
        const Exit& exit = exits[exitIndex];
        const Anchor anchor = anchors[exit.anchorId];
        for(const auto& anchorMarkerInterval: anchor) {
            OrientedReadInfo* orientedReadInfo = getOrientedReadInfo(anchorMarkerInterval.orientedReadId);
            if(orientedReadInfo) {
                orientedReadInfo->exitIndex = exitIndex;
                orientedReadInfo->exitPositionInJourney = anchorMarkerInterval.positionInJourney;
            }
        }
    }



    // Fill in the journeyBegin and journeyEnd fields for each OrientedReadInfo.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        const auto journey = anchors.journeys[orientedReadId.getValue()];

        // Start with the entire journey.
        orientedReadInfo.journeyBegin = 0;
        orientedReadInfo.journeyEnd = journey.size();

        // If bidirectional is false, clip it at the entrances/exits as appropriate.
        if(not bidirectional) {
            if(orientedReadInfo.entranceIndex != invalid<uint64_t>) {
                orientedReadInfo.journeyBegin = orientedReadInfo.entrancePositionInJourney;
            }
            if(orientedReadInfo.exitIndex != invalid<uint64_t>) {
                orientedReadInfo.journeyEnd = orientedReadInfo.exitPositionInJourney + 1;
            }
        }
    }


    // Remove any OrientedReadInfos for which journeyBegin is not less than journeyEnd.
    {
        vector<OrientedReadInfo> newOrientedReadInfos;
        for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
            if(orientedReadInfo.journeyBegin < orientedReadInfo.journeyEnd) {
                newOrientedReadInfos.push_back(orientedReadInfo);
            }
        }
        orientedReadInfos.swap(newOrientedReadInfos);
    }



    if(debug) {
        ofstream csv("TangleOrientedReads-" + to_string(tangleId) + ".csv");
        csv << "OrientedReadId,Journey length,"
            "Entrance index,Entrance,Entrance position in journey,"
            "Exit index,Exit,Exit position in journey,Type,Journey begin,Journey end," << endl;
        for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
            const auto journey = anchors.journeys[orientedReadInfo.orientedReadId.getValue()];
            csv << orientedReadInfo.orientedReadId << ",";
            csv << journey.size() << ",";

            // Entrance information.
            if(orientedReadInfo.entranceIndex == invalid<uint64_t>) {
                csv << ",,,";
            } else {
                const Entrance& entrance = entrances[orientedReadInfo.entranceIndex];
                csv <<
                    orientedReadInfo.entranceIndex << "," <<
                    anchorIdToString(entrance.anchorId) << "," <<
                    orientedReadInfo.entrancePositionInJourney << ",";

            }

            // Exit information.
            if(orientedReadInfo.exitIndex == invalid<uint64_t>) {
                csv << ",,,";
            } else {
                const Exit& exit = exits[orientedReadInfo.exitIndex];
                csv <<
                    orientedReadInfo.exitIndex << "," <<
                    anchorIdToString(exit.anchorId) << "," <<
                    orientedReadInfo.exitPositionInJourney << ",";
            }

            // Type.
            if(orientedReadInfo.entranceIndex == invalid<uint64_t>) {
                SHASTA_ASSERT(orientedReadInfo.exitIndex != invalid<uint64_t>);
                csv << "Exit only,";
            } else {
                if(orientedReadInfo.exitIndex == invalid<uint64_t>) {
                    csv << "Entrance only,";
                } else {
                    csv << "Both,";
                }
            }

            // Journey information.
            csv << orientedReadInfo.journeyBegin << ",";
            csv << orientedReadInfo.journeyEnd << ",";

            csv << endl;
        }
    }


    // We can compute a tangle matrix at tangle boundary by just counting the oriented reads
    // present in each entrance/exit pair.
    vector< vector<uint64_t> > tangleMatrixAtBoundary(entrances.size(), vector<uint64_t>(exits.size(), 0));
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const uint64_t entranceIndex = orientedReadInfo.entranceIndex;
        if(orientedReadInfo.entranceIndex == invalid<uint64_t>) {
            continue;
        }
        const uint64_t exitIndex = orientedReadInfo.exitIndex;
        if(exitIndex == invalid<uint64_t>) {
            continue;
        }
        ++tangleMatrixAtBoundary[entranceIndex][exitIndex];
    }

    if(debug) {
        cout << "Tangle matrix at tangle boundary:" << endl;
        for(uint64_t entranceIndex=0; entranceIndex<entrances.size(); entranceIndex++) {
            const Entrance& entrance = entrances[entranceIndex];
            for(uint64_t exitIndex=0; exitIndex<exits.size(); exitIndex++) {
                const Exit& exit = exits[exitIndex];
                cout <<
                    "(" << entranceIndex << "," << exitIndex << ") " <<
                    anchorIdToString(entrance.anchorId) << " " <<
                    anchorIdToString(exit.anchorId) << " " << tangleMatrixAtBoundary[entranceIndex][exitIndex] << endl;
            }
        }
    }
}



TangleGraph::OrientedReadInfo* TangleGraph::getOrientedReadInfo(OrientedReadId orientedReadId)
{
    const auto it = std::lower_bound(
        orientedReadInfos.begin(), orientedReadInfos.end(),
        OrientedReadInfo(orientedReadId));;

    if(it == orientedReadInfos.end()) {
        return 0;
    }
    if(it->orientedReadId != orientedReadId)
    {
        return 0;
    }

    return &(*it);
}



#if 0
void TangleGraph::computeTangleMatrix()
{
    if(debug) {
        cout << "Tangle matrix:" << endl;
    }

    tangleMatrix.resize(entrances.size(), vector<uint64_t>(exits.size()));

    for(uint64_t entranceIndex=0; entranceIndex<entrances.size(); entranceIndex++) {
        const Entrance& entrance = entrances[entranceIndex];
        for(uint64_t exitIndex=0; exitIndex<exits.size(); exitIndex++) {
            const Exit& exit = exits[exitIndex];

            vector<AnchorId> commonUniqueAnchors;
            std::set_intersection(
                entrance.uniqueJourneyAnchorIds.begin(), entrance.uniqueJourneyAnchorIds.end(),
                exit.uniqueJourneyAnchorIds.begin(), exit.uniqueJourneyAnchorIds.end(),
                back_inserter(commonUniqueAnchors));
            tangleMatrix[entranceIndex][exitIndex] = commonUniqueAnchors.size();

            if(debug) {
                cout <<
                    "(" << entranceIndex << "," << exitIndex << ") " <<
                    anchorIdToString(entrance.anchorId) << " " <<
                    anchorIdToString(exit.anchorId) << " " << tangleMatrix[entranceIndex][exitIndex] << endl;
            }

        }
    }

}
#endif



bool TangleGraph::createVertices(uint64_t anchorCoverageThreshold, uint64_t minVertexCoverage)
{
    TangleGraph& tangleGraph = *this;

    // Gather the AnchorIds in the journeys of oriented reads that appear in each entrance.
    vector< vector<AnchorId> > entranceAnchorIds(entrances.size());
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const uint64_t entranceIndex = orientedReadInfo.entranceIndex;
        if(entranceIndex == invalid<uint64_t>) {
            // This oriented read does not appear in any entrance.
            continue;
        }

        // Gather AnchorIds reached by this oriented read.
        const auto journey = anchors.journeys[orientedReadInfo.orientedReadId.getValue()];
        SHASTA_ASSERT(orientedReadInfo.journeyBegin < orientedReadInfo.journeyEnd);
        copy(
            journey.begin() + orientedReadInfo.journeyBegin,
            journey.begin() + orientedReadInfo.journeyEnd,
            back_inserter(entranceAnchorIds[entranceIndex]));
    }

    // Deduplicate the AnchorIds for each entrance.
    // Only keep the ones that appear in at least anchorCoverageThreshold oriented reads.
    vector<uint64_t> count;
    for(auto& v: entranceAnchorIds) {
        deduplicateAndCountWithThreshold(v, count, anchorCoverageThreshold);
    }


    // Now find the AnchorIds that are in more than one entrance.
    vector<AnchorId> duplicateEntranceAnchorIds;
    for(auto& v: entranceAnchorIds) {
        copy(v.begin(), v.end(), back_inserter(duplicateEntranceAnchorIds));
    }
    deduplicateAndCountWithThreshold(duplicateEntranceAnchorIds, count, 2UL);



    // Do the same for the exits.
    // Gather the AnchorIds in the journeys of oriented reads that appear in each exit.
    vector< vector<AnchorId> > exitAnchorIds(exits.size());
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const uint64_t exitIndex = orientedReadInfo.exitIndex;
        if(exitIndex == invalid<uint64_t>) {
            // This oriented read does not appear in any exit.
            continue;
        }

        // Gather AnchorIds reached by this oriented read.
        const auto journey = anchors.journeys[orientedReadInfo.orientedReadId.getValue()];
        SHASTA_ASSERT(orientedReadInfo.journeyBegin < orientedReadInfo.journeyEnd);
        copy(
            journey.begin() + orientedReadInfo.journeyBegin,
            journey.begin() + orientedReadInfo.journeyEnd,
            back_inserter(exitAnchorIds[exitIndex]));
    }

    // Deduplicate the AnchorIds for each exit.
    for(auto& v: exitAnchorIds) {
        deduplicateAndCountWithThreshold(v, count, anchorCoverageThreshold);
    }

    // Now find the AnchorIds that are in more than one exit.
    vector<AnchorId> duplicateExitAnchorIds;
    for(auto& v: exitAnchorIds) {
        copy(v.begin(), v.end(), back_inserter(duplicateExitAnchorIds));
    }
    deduplicateAndCountWithThreshold(duplicateExitAnchorIds, count, 2UL);



    // The forbiddenAnchorIds are the union set of
    // duplicateEntranceAnchorIds and duplicateExitAnchorIds.
    vector<AnchorId> forbiddenAnchorIds;
    std::set_union(
        duplicateEntranceAnchorIds.begin(), duplicateEntranceAnchorIds.end(),
        duplicateExitAnchorIds.begin(), duplicateExitAnchorIds.end(),
        back_inserter(forbiddenAnchorIds));


    // Generate the set of all allowed anchorIds.
    vector<AnchorId> allAnchorIds;
    for(auto& v: entranceAnchorIds) {
        copy(v.begin(), v.end(), back_inserter(allAnchorIds));
    }
    for(auto& v: exitAnchorIds) {
        copy(v.begin(), v.end(), back_inserter(allAnchorIds));
    }
    deduplicate(allAnchorIds);
    vector<AnchorId> allowedAnchorIds;
    std::set_difference(
        allAnchorIds.begin(), allAnchorIds.end(),
        forbiddenAnchorIds.begin(), forbiddenAnchorIds.end(),
        back_inserter(allowedAnchorIds));

    // If an entrance or exit is not in this allowed set, give up.
    for(const Entrance& entrance: entrances) {
        if(not binary_search(allowedAnchorIds.begin(), allowedAnchorIds.end(), entrance.anchorId)) {
            if(debug) {
                cout << "Giving up because entrance anchor " << anchorIdToString(entrance.anchorId) <<
                    " is not in the allowed set." << endl;
            }
            return false;
        }
    }
    for(const Exit& exit: exits) {
        if(not binary_search(allowedAnchorIds.begin(), allowedAnchorIds.end(), exit.anchorId)) {
            if(debug) {
                cout << "Giving up because exit anchor " << anchorIdToString(exit.anchorId) <<
                    " is not in the allowed set." << endl;
            }
            return false;
        }
    }



    // Now we generate one vertex for each of these AnchorIds.
    for(const AnchorId anchorId: allowedAnchorIds) {
        const vertex_descriptor v = add_vertex(TangleGraphVertex(anchorId), tangleGraph);
        vertexTable.push_back(make_pair(anchorId, v));
        // cout << "Added to vertexTable " << anchorIdToString(anchorId) << " " << v << endl;
    }

    // At this point the vertexTable is valid and we can use getVertex.

    // Sanity checks.
    for(const Entrance& entrance: entrances) {
        SHASTA_ASSERT(getVertex(entrance.anchorId) != null_vertex());
    }
    for(const Exit& exit: exits) {
        SHASTA_ASSERT(getVertex(exit.anchorId) != null_vertex());
    }
    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        SHASTA_ASSERT(getVertex(tangleGraph[v].anchorId) == v);
    }



    // Find the oriented reads that visit each vertex.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        const auto globalJourney = anchors.journeys[orientedReadId.getValue()];

        // Loop over the portion of the global journey we selected for this oriented read.
        const uint64_t begin = orientedReadInfo.journeyBegin;
        const uint64_t end = orientedReadInfo.journeyEnd;
        SHASTA_ASSERT(begin < end);
        for(uint64_t positionInJourney=begin; positionInJourney!=end; positionInJourney++) {
            const AnchorId anchorId = globalJourney[positionInJourney];
            const vertex_descriptor v = getVertex(anchorId);
            if(v != null_vertex()) {
                SHASTA_ASSERT(tangleGraph[v].anchorId == anchorId);
                tangleGraph[v].orientedReadIds.push_back(orientedReadId);
            }
        }
    }



    // Remove low coverage vertices.
    // Don't allow removing entrances and exits.
    {
        vector<vertex_descriptor> verticesToBeRemoved;
        BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
            const TangleGraphVertex& vertex = tangleGraph[v];
            if(vertex.coverage() < minVertexCoverage) {
                const AnchorId anchorId = vertex.anchorId;
                if(not (isEntrance(anchorId) or isExit(anchorId))) {
                    verticesToBeRemoved.push_back(v);
                }
            }
        }
        for(const vertex_descriptor v: verticesToBeRemoved) {
            boost::remove_vertex(v, tangleGraph);
        }

        // We also need to recreate the vertexTable.
        vertexTable.clear();
        BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
            vertexTable.push_back(make_pair(tangleGraph[v].anchorId, v));
        }
        sort(vertexTable.begin(), vertexTable.end());
    }




    // Now we can fill in the tangleJourney of each OrientedReadInfo.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        const auto globalJourney = anchors.journeys[orientedReadId.getValue()];

        // Loop over the portion of the global journey we selected for this oriented read.
        const uint64_t begin = orientedReadInfo.journeyBegin;
        const uint64_t end = orientedReadInfo.journeyEnd;
        SHASTA_ASSERT(begin < end);
        for(uint64_t positionInJourney=begin; positionInJourney!=end; positionInJourney++) {
            const AnchorId anchorId = globalJourney[positionInJourney];
            const vertex_descriptor v = getVertex(anchorId);
            if(v != null_vertex()) {
                orientedReadInfo.tangleJourney.push_back(v);
            }
        }

        if(false) {
            cout << "The tangle journey for " << orientedReadId <<
                " has " << orientedReadInfo.tangleJourney.size() << " anchors." << endl;
        }
    }

    return true;
}



TangleGraph::vertex_descriptor TangleGraph::getVertex(AnchorId anchorId) const
{
    auto it = lower_bound(
        vertexTable.begin(), vertexTable.end(),
        make_pair(anchorId, null_vertex()),
        OrderPairsByFirstOnly<AnchorId, vertex_descriptor>());
    if(it == vertexTable.end()) {
        return null_vertex();
    }
    if(it->first != anchorId) {
        return null_vertex();
    }

    return it->second;
}




void TangleGraph::createEdges()
{
    TangleGraph& tangleGraph = *this;

    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        const vector<vertex_descriptor>& tangleJourney = orientedReadInfo.tangleJourney;

        for(uint64_t i1=1; i1<tangleJourney.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const vertex_descriptor v0 = tangleJourney[i0];
            const vertex_descriptor v1 = tangleJourney[i1];

            // Find this edge, and create it if necessary.
            edge_descriptor e;
            bool edgeWasFound = false;
            tie(e, edgeWasFound) = boost::edge(v0, v1, tangleGraph);
            if(not edgeWasFound) {
                bool edgeWasAdded = false;
                tie(e, edgeWasAdded) = add_edge(v0, v1, tangleGraph);
                SHASTA_ASSERT(edgeWasAdded);
            }

            // Store this OrientedReadId in the edge.
            tangleGraph[e].orientedReadIds.push_back(orientedReadId);
        }
    }
}




void TangleGraph::writeGraphviz(const string& name) const
{
    const TangleGraph& tangleGraph = *this;
    ofstream dot("TangleGraph-" + to_string(tangleId) + "-" + name + ".dot");
    dot << "digraph TangleGraph" << tangleId << "{\n";

    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        const TangleGraphVertex& vertex = tangleGraph[v];
        const AnchorId anchorId = vertex.anchorId;
        dot << "\"" << anchorIdToString(anchorId) << "\"";

        // Begin attributes.
        dot << " [";



        // Label.
        dot << "label=\"";
        dot << anchorIdToString(anchorId) << "\\n" << vertex.coverage() << "\\n";
        dot << "\"";



        // Color.
        const bool vertexIsEntrance = isEntrance(anchorId);
        const bool vertexIsExit = isExit(anchorId);
        if(vertexIsEntrance and not vertexIsExit) {
            dot << " style=filled color=green fillcolor=green";
        }
        if(vertexIsExit and not vertexIsEntrance) {
            dot << " style=filled color=red fillcolor=red";
        }
        if(vertexIsEntrance and vertexIsExit) {
            dot << " style=filled color=magenta fillcolor=magenta";
        }



        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, tangleGraph, TangleGraph) {
        const TangleGraphEdge& edge = tangleGraph[e];
        const vertex_descriptor v0 = source(e, tangleGraph);
        const vertex_descriptor v1 = target(e, tangleGraph);

        const AnchorId anchorId0 = tangleGraph[v0].anchorId;
        const AnchorId anchorId1 = tangleGraph[v1].anchorId;

        dot << "\"" << anchorIdToString(anchorId0) << "\"->";
        dot << "\"" << anchorIdToString(anchorId1) << "\"";

        // Begin attributes.
        dot << " [";

        // Tooltip.
        dot <<
            "tooltip=\"" <<
            anchorIdToString(anchorId0) << "->" <<
            anchorIdToString(anchorId1) << " " <<
            edge.coverage() << " " <<
            edgeLoss(e) << "\"";

        // Thickness.
        dot << "penwidth=" << 0.5 * double(edge.coverage());

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";
    }

    dot << "}\n";
}



// Remove edges for which loss = (commonCount - coverage) / commonCount > maxLoss
// This is similar to AnchorGraph::removeWeakEdges.
void TangleGraph::removeWeakEdges(double maxLoss)
{
    TangleGraph& tangleGraph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, tangleGraph, TangleGraph) {
        if(edgeLoss(e) > maxLoss) {
            edgesToBeRemoved.push_back(e);
        }
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, tangleGraph);
    }

    if(debug) {
        cout << "Removed " << edgesToBeRemoved.size() << " weak edges." << endl;
    }

}




double TangleGraph::edgeLoss(edge_descriptor e) const
{
    const TangleGraph& tangleGraph = *this;

    const TangleGraphEdge& edge = tangleGraph[e];
    const vertex_descriptor v0 = source(e, tangleGraph);
    const vertex_descriptor v1 = target(e, tangleGraph);

    // Find the number of common oriented reads between the two vertices.
    const vector<OrientedReadId>& orientedReadIds0 = tangleGraph[v0].orientedReadIds;
    const vector<OrientedReadId>& orientedReadIds1 = tangleGraph[v1].orientedReadIds;
    vector<OrientedReadId> commonOrientedReadIds;
    std::set_intersection(
        orientedReadIds0.begin(), orientedReadIds0.end(),
        orientedReadIds1.begin(), orientedReadIds1.end(),
        back_inserter(commonOrientedReadIds));
    const uint64_t commonCount = commonOrientedReadIds.size();

    return double(commonCount - edge.coverage()) / double(commonCount);
}



// Remove cross-edges.
// This removes an edge v0->v1 if the following are all true:
// - Its coverage is at most lowCoverageThreshold.
// - v0 has at least one out-edge with coverage at least highCoverageThreshold
// - v1 has at least one in-edge with coverage at least highCoverageThreshold.
void TangleGraph::removeCrossEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    TangleGraph& tangleGraph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, tangleGraph, TangleGraph) {
        const TangleGraphEdge& edge = tangleGraph[e];

        // Check coverage.
        if(edge.coverage() > lowCoverageThreshold) {
            continue;
        }

        // Check out-edges of v0.
        const vertex_descriptor v0 = source(e, tangleGraph);
        bool v0HasStrongOutEdge = false;
        BGL_FORALL_OUTEDGES(v0, e0, tangleGraph, TangleGraph) {
            if(tangleGraph[e0].coverage() >= highCoverageThreshold) {
                v0HasStrongOutEdge = true;
                break;
            }
        }
        if(not v0HasStrongOutEdge) {
            continue;
        }

        // Check in-edges of v1.
        const vertex_descriptor v1 = target(e, tangleGraph);
        bool v1HasStrongOutEdge = false;
        BGL_FORALL_INEDGES(v1, e1, tangleGraph, TangleGraph) {
            if(tangleGraph[e1].coverage() >= highCoverageThreshold) {
                v1HasStrongOutEdge = true;
                break;
            }
        }
        if(not v1HasStrongOutEdge) {
            continue;
        }

        // If all above checks passed, this edge will be removed.
        edgesToBeRemoved.push_back(e);
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(debug) {
            const vertex_descriptor v0 = source(e, tangleGraph);
            const vertex_descriptor v1 = target(e, tangleGraph);
            cout << "Removed cross edge " << anchorIdToString(tangleGraph[v0].anchorId) <<
                " -> " << anchorIdToString(tangleGraph[v1].anchorId) << endl;
        }
        boost::remove_edge(e, tangleGraph);
    }
    if(debug) {
        cout << "Removed " << edgesToBeRemoved.size() << " cross edges." << endl;
    }
}



void TangleGraph::getChains(vector< vector<AnchorId> >& anchorChains) const
{
    const TangleGraph& tangleGraph = *this;

    // Find chains of edges (linear paths).
    vector< vector<edge_descriptor> > edgeChains;
    findLinearChains(tangleGraph, 1, edgeChains);
    const uint64_t chainCount = edgeChains.size();



    // Create the corresponding chains of AnchorIds.
    anchorChains.clear();
    vector<AnchorId> anchorChain;
    vector<uint64_t> splitPoints;
    for(uint64_t i=0; i<chainCount; i++) {
        const vector<edge_descriptor>& edgeChain = edgeChains[i];

        // Fill in the AnchorIds corresponding to this chain.
        anchorChain.clear();
        // Add the source AnchorId for the first edge.
        const vertex_descriptor v = source(edgeChain.front(), tangleGraph);
        anchorChain.push_back(tangleGraph[v].anchorId);
        // Add the target AnchorId for all the edges.
        for(const edge_descriptor e: edgeChain) {
            const vertex_descriptor v = target(e, tangleGraph);
            anchorChain.push_back(tangleGraph[v].anchorId);
        }


        // If any of the internal AnchorIds of this chain correspond to an
        // Entrance or Exit, we have to split it.
        splitPoints.clear();
        for(uint64_t i=1; i<anchorChain.size()-1; i++) {
            const AnchorId anchorId = anchorChain[i];
            if(isEntrance(anchorId) or isExit(anchorId)) {
                splitPoints.push_back(i);
            }
        }

        if(splitPoints.empty()) {
            anchorChains.push_back(anchorChain);
        } else {

            if(debug) {
                cout << "Splitting chain " << anchorChain.front() << "..." << anchorChain.back() << endl;
            }

            // Add a split chain from the beginning to the first split point.
            anchorChains.resize(anchorChains.size() + 1);
            vector<AnchorId>& firstSplitChain = anchorChains.back();
            for(uint64_t j=0; j<=splitPoints.front(); j++) {
                firstSplitChain.push_back(anchorChain[j]);
            }
            // Add split chains in-between split points.
            for(uint64_t i1=1; i1<splitPoints.size(); i1++) {
                const uint64_t i0 = i1 - 1;
                const uint64_t j0 = splitPoints[i0];
                const uint64_t j1 = splitPoints[i1];
                anchorChains.resize(anchorChains.size() + 1);
                vector<AnchorId>& splitChain = anchorChains.back();
                for(uint64_t j=j0; j<=j1; j++) {
                    splitChain.push_back(anchorChain[j]);
                }
            }
            // Add a split chain from the last split point to the end.
            anchorChains.resize(anchorChains.size() + 1);
            vector<AnchorId>& lastSplitChain = anchorChains.back();
            for(uint64_t j=splitPoints.back(); j<anchorChain.size(); j++) {
                lastSplitChain.push_back(anchorChain[j]);
            }
        }
    }

}



// Find out if a given AnchorId is an entrance.
bool TangleGraph::isEntrance(AnchorId anchorId) const
{
    for(const Entrance& entrance: entrances) {
        if(entrance.anchorId == anchorId) {
            return true;
        }
    }
    return false;
}



// Find out if a given AnchorId is an exit.
bool TangleGraph::isExit(AnchorId anchorId) const
{
    for(const Exit& exit: exits) {
        if(exit.anchorId == anchorId) {
            return true;
        }
    }
    return false;
}



bool TangleGraph::removeUnreachable()
{
    TangleGraph& tangleGraph = *this;

    // Clear the BFS flags for all vertices.
    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        tangleGraph[v].wasSeenByForwardBfs = false;
        tangleGraph[v].wasSeenByBackwardBfs = false;
    }

    // Do a forward BFS starting at the entrances.
    std::queue<vertex_descriptor> q;
    for(const Entrance& entrance: entrances) {
        const vertex_descriptor v = getVertex(entrance.anchorId);
        SHASTA_ASSERT(v != null_vertex());
        tangleGraph[v].wasSeenByForwardBfs = true;
        q.push(v);
    }
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        BGL_FORALL_OUTEDGES(v0, e, tangleGraph, TangleGraph) {
            const vertex_descriptor v1 = target(e, tangleGraph);
            TangleGraphVertex& vertex1 = tangleGraph[v1];
            if(not vertex1.wasSeenByForwardBfs) {
                vertex1.wasSeenByForwardBfs = true;
                q.push(v1);
            }
        }
    }



    // Do a backward BFS starting at the exits.
    SHASTA_ASSERT(q.empty());
    for(const Exit& exit: exits) {
        const vertex_descriptor v = getVertex(exit.anchorId);
        SHASTA_ASSERT(v != null_vertex());
        tangleGraph[v].wasSeenByBackwardBfs = true;
        q.push(v);
    }
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        BGL_FORALL_INEDGES(v0, e, tangleGraph, TangleGraph) {
            const vertex_descriptor v1 = source(e, tangleGraph);
            TangleGraphVertex& vertex1 = tangleGraph[v1];
            if(not vertex1.wasSeenByBackwardBfs) {
                vertex1.wasSeenByBackwardBfs = true;
                q.push(v1);
            }
        }
    }



    // Remove vertices that were not seen by both BFSs.
    // If this includes any entrance or exit, give up.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        const TangleGraphVertex& vertex = tangleGraph[v];
        const AnchorId anchorId = vertex.anchorId;
        if(vertex.wasSeenByForwardBfs and vertex.wasSeenByBackwardBfs) {
            continue;
        }
        if(isEntrance(anchorId)) {
            return false;
        }
        if(isExit(anchorId)) {
            return false;
        }
        verticesToBeRemoved.push_back(v);
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, tangleGraph);
        boost::remove_vertex(v, tangleGraph);
    }

    // Recreate the vertex table.
    vertexTable.clear();
    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        vertexTable.push_back(make_pair(tangleGraph[v].anchorId, v));
    }
    deduplicate(vertexTable);


    return true;
}



// Return true if successful, that is, all Entrances are
// connecte to at least one Exit, and all Exits are
// connected to at least one Entrance.
bool TangleGraph::isSuccessful() const
{
    const TangleGraph& tangleGraph = *this;
    if(failure) {
        return false;
    }

    for(const Entrance& entrance: entrances) {
        const AnchorId anchorId = entrance.anchorId;
        const vertex_descriptor v = getVertex(anchorId);
        const TangleGraphVertex& vertex = tangleGraph[v];
        if(vertex.wasSeenByForwardBfs and vertex.wasSeenByBackwardBfs) {
            continue;
        } else {
            return false;
        }
    }

    for(const Exit& exit: exits) {
        const AnchorId anchorId = exit.anchorId;
        const vertex_descriptor v = getVertex(anchorId);
        const TangleGraphVertex& vertex = tangleGraph[v];
        if(vertex.wasSeenByForwardBfs and vertex.wasSeenByBackwardBfs) {
            continue;
        } else {
            return false;
        }
    }

    return true;
}

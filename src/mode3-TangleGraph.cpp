// Shasta.
#include "mode3-TangleGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <algorithm.hpp>
#include <fstream.hpp>



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

    createVertices(minVertexCoverage);
    createEdges();
    if(debug) {
        cout << "The initial tangle graph has " << num_vertices(*this) <<
            " vertices and " << num_edges(*this) << " edges." << endl;
    }

    removeWeakEdges(maxLoss);
    writeGraphviz("A");
    removeCrossEdges(lowCoverageThreshold, highCoverageThreshold);

    if(debug) {
        cout << "The final tangle graph has " << num_vertices(*this) <<
            " vertices and " << num_edges(*this) << " edges." << endl;
        writeGraphviz("B");
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



// Create TangleGraph vertices.
// There is a vertex for each AnchorId that:
// - Appears in at least one Entrance and/or one Exit.
// - Does not appear in more than one Entrance.
// - Does not appear in more than one Exit.
void TangleGraph::createVertices(uint64_t minVertexCoverage)
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
        copy(
            journey.begin() + orientedReadInfo.journeyBegin,
            journey.begin() + orientedReadInfo.journeyEnd,
            back_inserter(entranceAnchorIds[entranceIndex]));
    }

    // Deduplicate the AnchorIds for each entrance.
    for(auto& v: entranceAnchorIds) {
        deduplicate(v);
    }


    // Now find the AnchorIds that are in more than one entrance.
    vector<AnchorId> duplicateEntranceAnchorIds;
    for(auto& v: entranceAnchorIds) {
        copy(v.begin(), v.end(), back_inserter(duplicateEntranceAnchorIds));
    }
    vector<uint64_t> count;
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
        copy(
            journey.begin() + orientedReadInfo.journeyBegin,
            journey.begin() + orientedReadInfo.journeyEnd,
            back_inserter(exitAnchorIds[exitIndex]));
    }

    // Deduplicate the AnchorIds for each exit.
    for(auto& v: exitAnchorIds) {
        deduplicate(v);
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


    // Now we generate one vertex for each of these AnchorIds.
    for(const AnchorId anchorId: allowedAnchorIds) {
        const vertex_descriptor v = add_vertex(TangleGraphVertex(anchorId), tangleGraph);
        vertexTable.push_back(make_pair(anchorId, v));
    }


    // At this point the vertexTable is valid and we can use getVertex.

    // Fill in the entranceIndex and exitIndex of each vertex.
    for(uint64_t entranceIndex=0; entranceIndex<entrances.size(); entranceIndex++) {
        for(const AnchorId anchorId: entranceAnchorIds[entranceIndex]) {
            const vertex_descriptor v = getVertex(anchorId);
            if(v == null_vertex()) {
                continue;
            }
            TangleGraphVertex& vertex = tangleGraph[v];
            SHASTA_ASSERT(vertex.entranceIndex == invalid<uint64_t>);
            vertex.entranceIndex = entranceIndex;
        }
    }
    for(uint64_t exitIndex=0; exitIndex<exits.size(); exitIndex++) {
        for(const AnchorId anchorId: exitAnchorIds[exitIndex]) {
            const vertex_descriptor v = getVertex(anchorId);
            if(v == null_vertex()) {
                continue;
            }
            TangleGraphVertex& vertex = tangleGraph[v];
            SHASTA_ASSERT(vertex.exitIndex == invalid<uint64_t>);
            vertex.exitIndex = exitIndex;
        }
    }



    // Find the oriented reads that visit each vertex.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        const auto globalJourney = anchors.journeys[orientedReadId.getValue()];

        // Loop over the portion of the global journey we selected for this oriented read.
        const uint64_t begin = orientedReadInfo.journeyBegin;
        const uint64_t end = orientedReadInfo.journeyEnd;
        for(uint64_t positionInJourney=begin; positionInJourney!=end; positionInJourney++) {
            const AnchorId anchorId = globalJourney[positionInJourney];
            const vertex_descriptor v = getVertex(anchorId);
            if(v != null_vertex()) {
                tangleGraph[v].orientedReadIds.push_back(orientedReadId);
            }
        }
    }


    // Remove low coverage vertices.
    {
        vector<vertex_descriptor> verticesToBeRemoved;
        BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
            if(tangleGraph[v].coverage() < minVertexCoverage) {
                verticesToBeRemoved.push_back(v);
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
        for(uint64_t positionInJourney=begin; positionInJourney!=end; positionInJourney++) {
            const AnchorId anchorId = globalJourney[positionInJourney];
            const vertex_descriptor v = getVertex(anchorId);
            if(v != null_vertex()) {
                orientedReadInfo.tangleJourney.push_back(v);
            }
        }

        if(debug) {
            cout << "The tangle journey for " << orientedReadId <<
                " has " << orientedReadInfo.tangleJourney.size() << " anchors." << endl;
        }
    }



    // Histogram the vertices by their entranceIndex and exitIndex.
    std::map< pair<uint64_t, uint64_t>, uint64_t> histogram;
    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        const TangleGraphVertex& vertex = tangleGraph[v];
        const auto p = make_pair(vertex.entranceIndex, vertex.exitIndex);
        auto it = histogram.find(p);
        if(it == histogram.end()) {
            histogram.insert({p, 1});
        } else {
            ++(it->second);
        }
    }



    if(debug) {
        cout << "Histogram of vertex count by (entrance, exit):" << endl;
        for(const auto& p: histogram) {
            const auto& q = p.first;
            const uint64_t frequency = p.second;
            const uint64_t entranceIndex = q.first;
            const uint64_t exitIndex = q.second;

            cout << "Entrance: ";
            if(entranceIndex == invalid<uint64_t>) {
                cout << "None";
            } else {
                cout << entranceIndex << " " << anchorIdToString(entrances[entranceIndex].anchorId);
            }

            cout << ", exit: ";
            if(exitIndex == invalid<uint64_t>) {
                cout << "None";
            } else {
                cout << exitIndex << " " << anchorIdToString(exits[exitIndex].anchorId);
            }

            cout << ", frequency " << frequency << "\n";
        }
    }




    // Use the histogram to fill in the tangleMatrix.
    tangleMatrix.resize(entrances.size(), vector<uint64_t>(exits.size(), 0));
    for(const auto& p: histogram) {
        const auto& q = p.first;
        const uint64_t frequency = p.second;
        const uint64_t entranceIndex = q.first;
        if(entranceIndex == invalid<uint64_t>) {
            continue;
        }

        const uint64_t exitIndex = q.second;
        if(exitIndex == invalid<uint64_t>) {
            continue;
        }

        tangleMatrix[entranceIndex][exitIndex] = frequency;
    }



    if(debug) {
        cout << "Tangle matrix:" << endl;
        for(uint64_t entranceIndex=0; entranceIndex<entrances.size(); entranceIndex++) {
            const Entrance& entrance = entrances[entranceIndex];
            for(uint64_t exitIndex=0; exitIndex<exits.size(); exitIndex++) {
                const Exit& exit = exits[exitIndex];
                cout <<
                    "(" << entranceIndex << "," << exitIndex << ") " <<
                    anchorIdToString(entrance.anchorId) << " " <<
                    anchorIdToString(exit.anchorId) << " " << tangleMatrix[entranceIndex][exitIndex] << endl;
            }
        }
    }
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

        if(vertex.entranceIndex == invalid<uint64_t>) {
            dot << "N";
        } else {
            dot << vertex.entranceIndex;
        }
        dot << ",";
        if(vertex.exitIndex == invalid<uint64_t>) {
            dot << "N";
        } else {
            dot << vertex.exitIndex;
        }

        // End of label.
        dot << "\"";



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

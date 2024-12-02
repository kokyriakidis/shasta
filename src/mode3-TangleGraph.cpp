#include "mode3-TangleGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode3;


TangleGraph::TangleGraph(
    bool debug,
    const Anchors& anchors,
    const vector<AnchorId>& entranceAnchors,
    const vector<AnchorId>& exitAnchors,
    bool bidirectional) :
    debug(debug),
    anchors(anchors),
    bidirectional(bidirectional)
{
    if(debug) {
        cout << "Creating a tangle graph with " << entranceAnchors.size() <<
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
    computeTangleMatrix();
}



void TangleGraph::constructEntrances(const vector<AnchorId>& entranceAnchors)
{
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
}



void TangleGraph::constructExits(const vector<AnchorId>& exitAnchors)
{
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
}



TangleGraph::EntranceOrExit::EntranceOrExit(
    AnchorId anchorId,
    const Anchor& anchor) :
    anchorId(anchorId)
{
    copy(anchor.begin(), anchor.end(), back_inserter(anchorMarkerIntervals));
}



// This fills in the journeyAnchorIds.
void TangleGraph::Entrance::readFollowing(
    bool debug, const Anchors& anchors, bool bidirectional)
{
    for(const AnchorMarkerInterval& anchorMarkerInterval: anchorMarkerIntervals) {
        const OrientedReadId orientedReadId = anchorMarkerInterval.orientedReadId;
        const auto journey = anchors.journeys[orientedReadId.getValue()];

        const uint64_t begin = (bidirectional ? 0 : anchorMarkerInterval.positionInJourney + 1);
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
        const uint64_t end = (bidirectional ? journey.size() : anchorMarkerInterval.positionInJourney);
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
                cout << anchorIdToString(entrance.anchorId) << " " <<
                    anchorIdToString(exit.anchorId) << " " << tangleMatrix[entranceIndex][exitIndex] << endl;
            }

        }
    }

}

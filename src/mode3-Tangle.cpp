// Shasta.
#include "mode3-Tangle.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>


Tangle::Tangle(
    bool debug,
    AssemblyGraph& assemblyGraph,
    uint64_t maxOffset,
    const vector<AssemblyGraph::vertex_descriptor>& tangleVerticesArgument) :
    debug(debug),
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVerticesArgument)
{
    if(debug) {
        cout << "Working on a tangle with " << tangleVertices.size() << " vertices." << endl;
    }

    // Sort the tangleVertices so we can do binary searches in them
    // in isTangleVertex.
    sort(tangleVertices.begin(), tangleVertices.end());
    if(debug) {
        writeTangleVertices();
    }

    findTangleEdges(maxOffset);
    if(debug) {
        writeTangleEdges();
    }

    findEntrances();
    findExits();
    if(debug) {
        writeEntrances();
        writeExits();
    }

    readFollowingFromEntrances();
    readFollowingFromExits();

    // Compute a "tangle matrix" using the number of common unique AnchorIds between
    // each entrance/exit pair.
    if(debug) {
        cout << "Tangle matrix:" << endl;
    }
    for(const Entrance& entrance: entrances) {
        for(const Exit& exit: exits) {
            vector<AnchorId> commonUniqueAnchors;
            std::set_intersection(
                entrance.uniqueJourneyAnchorIds.begin(), entrance.uniqueJourneyAnchorIds.end(),
                exit.uniqueJourneyAnchorIds.begin(), exit.uniqueJourneyAnchorIds.end(),
                back_inserter(commonUniqueAnchors));
            if(debug) {
                cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
                    assemblyGraph.bubbleChainStringId(exit.e) <<
                    " " << commonUniqueAnchors.size() << endl;
            }
        }

    }

}



bool Tangle::isTangleVertex(AssemblyGraph::vertex_descriptor v) const
{
    return binary_search(tangleVertices.begin(), tangleVertices.end(), v);
}



void Tangle::writeTangleVertices() const
{
    cout << "Tangle vertices:";
    for(const AssemblyGraph::vertex_descriptor v: tangleVertices) {
        const AnchorId anchorId = assemblyGraph[v].getAnchorId();
        cout << " " << anchorIdToString(anchorId);
    }
    cout << endl;

}



void Tangle::findTangleEdges(uint64_t maxOffset)
{
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;

    // Loop over vertices in the Tangle.
    for(const vertex_descriptor v0: tangleVertices) {

        // Loop over its out-edges.
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {

            // Extract the Chain corresponding to this edge.
            const BubbleChain& bubbleChain = assemblyGraph[e];
            SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());

            // If this is a Tangle edge, store it.
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(isTangleVertex(v1)) {
                const Chain& chain = bubbleChain.getOnlyChain();
                if(assemblyGraph.chainOffset(chain) <= maxOffset) {
                    tangleEdges.push_back(e);
                }
            }
        }
    }

    sort(tangleEdges.begin(), tangleEdges.end());

}



bool Tangle::isTangleEdge(AssemblyGraph::edge_descriptor e) const
{
    return binary_search(tangleEdges.begin(), tangleEdges.end(), e);
}



void Tangle::writeTangleEdges() const
{
    cout << "Tangle edges:";
    for(const AssemblyGraph::edge_descriptor e: tangleEdges) {
        cout << " " << assemblyGraph.bubbleChainStringId(e);
    }
    cout << endl;
}



Tangle::EntranceOrExit::EntranceOrExit(
    AssemblyGraph::edge_descriptor e,
    AnchorId anchorId,
    const Anchor& anchor) :
    e(e),
    anchorId(anchorId)
{
    copy(anchor.begin(), anchor.end(), back_inserter(anchorMarkerIntervals));
}



// The entrances are AssemblyGraph edges that are not in the Tangle
// but whose target vertex is in the Tangle.
void Tangle::findEntrances()
{
    for(const AssemblyGraph::vertex_descriptor v: tangleVertices) {
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not isTangleEdge(e)) {
                const Chain& chain = assemblyGraph[e].getOnlyChain();
                const AnchorId anchorId = chain.secondToLast();
                entrances.push_back(Entrance(e, anchorId, assemblyGraph.anchors[anchorId]));
            }
        }
    }

    // If an OrientedReadId appears in more than one entrance,
    // we want to remove the corresponding AnchorMarkerIntervals
    // from all entrances.
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
                " and will be neglected when following oriented reads forward:";
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
}



// The exits are AssemblyGraph edges that are not in the Tangle
// but whose source vertex is in the Tangle.
void Tangle::findExits()
{
    for(const AssemblyGraph::vertex_descriptor v: tangleVertices) {
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not isTangleEdge(e)) {
                const Chain& chain = assemblyGraph[e].getOnlyChain();
                const AnchorId anchorId = chain.second();
                exits.push_back(Exit(e, anchorId, assemblyGraph.anchors[anchorId]));
            }
        }
    }

    // If an OrientedReadId appears in more than one exit,
    // we want to remove the corresponding AnchorMarkerIntervals
    // from all exits.
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
            cout << "The following oriented reads appear in more than once exit "
                " and will be neglected when following oriented reads forward:";
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
}



bool Tangle::isEntrance(AssemblyGraph::edge_descriptor e) const
{
    for(const Entrance& entrance: entrances) {
        if(entrance.e == e) {
            return true;
        }
    }
    return false;
}



bool Tangle::isExit(AssemblyGraph::edge_descriptor e) const
{
    for(const Exit& exit: exits) {
        if(exit.e == e) {
            return true;
        }
    }
    return false;
}



void Tangle::writeEntrances() const
{
    cout << "Entrances:" << endl;
    for(const Entrance& entrance: entrances) {
        cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
            anchorIdToString(entrance.anchorId) << endl;
    }
}



void Tangle::writeExits() const
{
    cout << "Exits:" << endl;
    for(const Exit& exit: exits) {
        cout << assemblyGraph.bubbleChainStringId(exit.e) << " " <<
            anchorIdToString(exit.anchorId) << endl;
    }
}



void Tangle::readFollowingFromEntrances()
{
    for(Entrance& entrance: entrances) {
        readFollowingFromEntrance(entrance);
    }

    // Find the AnchorIds that were found in more than one Entrance.
    vector<AnchorId> duplicateAnchorIds;
    for(Entrance& entrance: entrances) {
        copy(entrance.journeyAnchorIds.begin(), entrance.journeyAnchorIds.end(),
            back_inserter(duplicateAnchorIds));
    }
    vector<uint64_t> count;
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
            cout << "Read following for entrance " << assemblyGraph.bubbleChainStringId(entrance.e) <<
                " found " << entrance.uniqueJourneyAnchorIds.size() << " anchors unique to this entrance." << endl;
        }
    }
}



void Tangle::readFollowingFromEntrance(Entrance& entrance)
{
    for(const AnchorMarkerInterval& anchorMarkerInterval: entrance.anchorMarkerIntervals) {
        const OrientedReadId orientedReadId = anchorMarkerInterval.orientedReadId;
        const auto journey = assemblyGraph.anchors.journeys[orientedReadId.getValue()];
        for(uint64_t position=anchorMarkerInterval.positionInJourney+1; position<journey.size(); position++) {
            entrance.journeyAnchorIds.push_back(journey[position]);
        }
    }

    if(debug) {
        cout << "Read following for entrance " << assemblyGraph.bubbleChainStringId(entrance.e) <<
            " found " << entrance.journeyAnchorIds.size() << " anchors before deduplication." << endl;
    }

    deduplicate(entrance.journeyAnchorIds);
    if(debug) {
        cout << "After deduplication, read following found " << entrance.journeyAnchorIds.size() << " anchors." << endl;
    }
}



void Tangle::readFollowingFromExits()
{
    for(Exit& exit: exits) {
        readFollowingFromExit(exit);
    }

    // Find the AnchorIds that were found in more than one Exit.
    vector<AnchorId> duplicateAnchorIds;
    for(Exit& exit: exits) {
        copy(exit.journeyAnchorIds.begin(), exit.journeyAnchorIds.end(),
            back_inserter(duplicateAnchorIds));
    }
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(duplicateAnchorIds, count, 2UL);

    if(debug) {
        cout << duplicateAnchorIds.size() <<
            " anchors were found by read following from more than one exit." << endl;
    }


    // Now we can fill in the uniqueJourneyAnchorIds for each Exit.
    for(Exit& exit: exits) {
        std::set_difference(
            exit.journeyAnchorIds.begin(), exit.journeyAnchorIds.end(),
            duplicateAnchorIds.begin(), duplicateAnchorIds.end(),
            back_inserter(exit.uniqueJourneyAnchorIds));

        if(debug) {
            cout << "Read following for exit " << assemblyGraph.bubbleChainStringId(exit.e) <<
                " found " << exit.uniqueJourneyAnchorIds.size() << " anchors unique to this exit." << endl;
        }
    }
}



void Tangle::readFollowingFromExit(Exit& exit)
{
    for(const AnchorMarkerInterval& anchorMarkerInterval: exit.anchorMarkerIntervals) {
        const OrientedReadId orientedReadId = anchorMarkerInterval.orientedReadId;
        const auto journey = assemblyGraph.anchors.journeys[orientedReadId.getValue()];
        for(uint64_t position=0; position<anchorMarkerInterval.positionInJourney; position++) {
            exit.journeyAnchorIds.push_back(journey[position]);
        }
    }

    if(debug) {
        cout << "Read following for exit " << assemblyGraph.bubbleChainStringId(exit.e) <<
            " found " << exit.journeyAnchorIds.size() << " anchors before deduplication." << endl;
    }

    deduplicate(exit.journeyAnchorIds);
    if(debug) {
        cout << "After deduplication, read following found " << exit.journeyAnchorIds.size() << " anchors." << endl;
    }
}
// Shasta.
#include "mode3-Detangler.hpp"
#include "mode3-AssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



Detangler::Detangler(
    bool debug,
    AssemblyGraph& assemblyGraph) :
    debug(debug),
    assemblyGraph(assemblyGraph)
{
}



void Detangler::writeInitialMessage(const vector<vertex_descriptor>& superbubble) const
{
    if(debug) {
        cout << "Found a superbubble with " << superbubble.size() <<
            " vertices:";
        for(const vertex_descriptor v: superbubble) {
            cout << " " << anchorIdToString(assemblyGraph[v].getAnchorId());
        }
        cout << endl;
    }
}



ChainDetangler::ChainDetangler(bool debug, AssemblyGraph& assemblyGraph) :
    Detangler(debug, assemblyGraph)
{
}



// This computes the tangle matrix using the second to last Anchor
// of each incoming chain and the second Anchor of each outgoing Chain.
void ChainDetangler::prepare(const vector<vertex_descriptor>& superbubble)
{
    // Find the Entrances.
    entrances.clear();
    for(const vertex_descriptor v: superbubble) {
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not assemblyGraph.isInternalToSuperbubble(e)) {
                entrances.push_back(Entrance(e, assemblyGraph));
            }
        }
    }

    // Find the ExitChains.
    exits.clear();
    for(const vertex_descriptor v: superbubble) {
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not assemblyGraph.isInternalToSuperbubble(e)) {
                exits.push_back(Exit(e, assemblyGraph));
            }
        }
    }

}



ChainDetangler::Entrance::Entrance(
    edge_descriptor e,
    const AssemblyGraph& assemblyGraph) :
    e(e),
    anchorId(assemblyGraph[e].getOnlyChain().secondToLast())
{
}



ChainDetangler::Exit::Exit(
    edge_descriptor e,
    const AssemblyGraph& assemblyGraph) :
    e(e),
    anchorId(assemblyGraph[e].getOnlyChain().second())
{
}



void ChainDetangler::writeEntrancesAndExits() const
{
    if(debug) {
        cout << entrances.size() << " entrances:";
        for(const Entrance& entrance: entrances) {
            cout << assemblyGraph.bubbleChainStringId(entrance.e) <<
                " " << anchorIdToString(entrance.anchorId) << endl;
        }

        cout << exits.size() << " exit chains:";
        for(const Exit& exit: exits) {
            cout << assemblyGraph.bubbleChainStringId(exit.e) <<
                " " << anchorIdToString(exit.anchorId) << endl;
        }
    }

}



// Return true if there is one or more Entrance/Exit pair
// with the same edge_descriptor.
bool ChainDetangler::commonChainsBetweenEntrancesAndExitsExists() const
{
    for(const Entrance& entrance: entrances) {
        for(const Exit& exit: exits) {
            if(entrance.e == exit.e) {
                return true;
            }
        }
    }

    return false;
}



// Return true if there is one or more Entrance/Exit pair
// with the same AnchorId.
bool ChainDetangler::commonAnchorsBetweenEntrancesAndExitsExists() const
{
    for(const Entrance& entrance: entrances) {
        for(const Exit& exit: exits) {
            if(entrance.anchorId == exit.anchorId) {
                return true;
            }
        }
    }

    return false;
}



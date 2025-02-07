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

    computeTangleMatrix();
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



void ChainDetangler::computeTangleMatrix()
{
    const uint64_t inDegree = entrances.size();
    const uint64_t outDegree = exits.size();

    totalCommonCoverage = 0;

    for(Entrance& entrance: entrances) {
        entrance.commonCoverage = 0;
    }
    for(Exit& exit: exits) {
        exit.commonCoverage = 0;
    }

    tangleMatrix.clear();
    tangleMatrix.resize(inDegree, vector<uint64_t>(outDegree));

    for(uint64_t iEntrance=0; iEntrance<inDegree; iEntrance++) {
        Entrance& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<outDegree; iExit++) {
            Exit& exit = exits[iExit];
            const uint64_t n = assemblyGraph.anchors.countCommon(entrance.anchorId, exit.anchorId, true);
            tangleMatrix[iEntrance][iExit] = n;
            totalCommonCoverage += n;
            entrance.commonCoverage += n;
            exit.commonCoverage += n;
        }
    }

    // Sanity checks.
    {
        uint64_t sum = 0;
        for(Entrance& entrance: entrances) {
            sum += entrance.commonCoverage;
        }
        SHASTA_ASSERT(sum == totalCommonCoverage);
    }
    {
        uint64_t sum = 0;
        for(Exit& exit: exits) {
            sum += exit.commonCoverage;
        }
        SHASTA_ASSERT(sum == totalCommonCoverage);
    }
}



void ChainDetangler::writeTangleMatrix() const
{
    if(debug) {
        cout << "Tangle matrix with total coverage " << totalCommonCoverage << ":" << endl;

        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            const Entrance& entrance = entrances[iEntrance];

            for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                const Exit& exit = exits[iExit];

                cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
                    assemblyGraph.bubbleChainStringId(exit.e) <<
                    " " << tangleMatrix[iEntrance][iExit];

                cout << endl;
            }
        }

        for(const Entrance& entrance: entrances) {
            cout << "Total common coverage for entrance " <<
                assemblyGraph.bubbleChainStringId(entrance.e) <<
                " is " << entrance.commonCoverage << endl;
        }
        for(const Exit& exit: exits) {
            cout << "Total common coverage for exit " <<
                assemblyGraph.bubbleChainStringId(exit.e) <<
                " is " << exit.commonCoverage << endl;
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



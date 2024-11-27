// Shasta.
#include "mode3-Tangle.hpp"
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



// The entrances are AssemblyGraph edges that are not in the Tangle
// but whose target vertex is in the Tangle.
void Tangle::findEntrances()
{
    for(const AssemblyGraph::vertex_descriptor v: tangleVertices) {
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not isTangleEdge(e)) {
                const Chain& chain = assemblyGraph[e].getOnlyChain();
                const AnchorId anchorId = chain.secondToLast();
                entrances.push_back(Entrance(e, anchorId));
            }
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
                exits.push_back(Exit(e, anchorId));
            }
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

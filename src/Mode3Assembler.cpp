// Shasta.
#include "Mode3Assembler.hpp"
#include "Assembler.hpp"
#include "dset64-gccAtomic.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard library.
#include "iostream.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Mode3Assembler>;



Mode3Assembler::Mode3Assembler(
    const Assembler& assembler,
    bool debug) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(assembler),
    assembler(assembler),
    debug(debug)
{
    performanceLog << timestamp << "Mode 3 assembly begins." << endl;

    gatherPrimaryMarkerGraphEdgeIds();
    computeConnectedComponents();

    performanceLog << timestamp << "Mode 3 assembly ends." << endl;
}



void Mode3Assembler::gatherPrimaryMarkerGraphEdgeIds()
{
    const auto& markerGraphEdges = assembler.markerGraph.edges;

    primaryMarkerGraphEdgeIds.clear();
    for(MarkerGraphEdgeId edgeId=0; edgeId<markerGraphEdges.size(); edgeId++) {
        if(markerGraphEdges[edgeId].isPrimary) {
            primaryMarkerGraphEdgeIds.push_back(edgeId);
        }
    }
    cout << "Of " << markerGraphEdges.size() << " marker graph edges, " <<
        primaryMarkerGraphEdgeIds.size() << " are primary." << endl;
}



// The oriented reads present in each primary marker graph edge
// define a bipartite graph. We want to compute connected components
// of this bipartite graph and process them one at a time.
void Mode3Assembler::computeConnectedComponents()
{
    performanceLog << timestamp << "Mode3Assembler::computeConnectedComponents begins." << endl;

    // Compute connected components of oriented reads.
    // Here oriented reads are indexed by OrientedReadId::getValue().
    const uint64_t orientedReadCount = assembler.markers.size();
    vector<DisjointSets::Aint> disjointSetsData(orientedReadCount);
    DisjointSets disjointSets(&disjointSetsData[0], orientedReadCount);

    // Loop over all primary marker graph edges.
    // This will be multithreaded.
    for(const MarkerGraphEdgeId edgeId: primaryMarkerGraphEdgeIds) {
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        SHASTA_ASSERT(not markerIntervals.empty());
        const OrientedReadId orientedReadId0 = markerIntervals.front().orientedReadId;
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId1 = markerInterval.orientedReadId;
            disjointSets.unite(orientedReadId0.getValue(), orientedReadId1.getValue());
        }
    }

    performanceLog << timestamp << "Mode3Assembler::computeConnectedComponents ends." << endl;
}

// Shasta.
#include "Mode3Assembler.hpp"
#include "Assembler.hpp"
#include "dset64-gccAtomic.hpp"
#include "orderPairs.hpp"
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

    // Compute connected components of the oriented reads portion
    // of the bipartite graph.
    // Here oriented reads are indexed by OrientedReadId::getValue().
    const uint64_t orientedReadCount = assembler.markers.size();
    vector<DisjointSets::Aint> disjointSetsData(orientedReadCount);
    DisjointSets disjointSets(&disjointSetsData[0], orientedReadCount);

    // Loop over all primary marker graph edges.
    // This could be multithreaded but runs at decent speed as is.
    for(const MarkerGraphEdgeId edgeId: primaryMarkerGraphEdgeIds) {
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        SHASTA_ASSERT(not markerIntervals.empty());
        const OrientedReadId orientedReadId0 = markerIntervals.front().orientedReadId;
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId1 = markerInterval.orientedReadId;
            disjointSets.unite(orientedReadId0.getValue(), orientedReadId1.getValue());
        }
    }

    // Gather the oriented reads in each connected component.
    vector< vector<OrientedReadId> > componentsOrientedReads(orientedReadCount);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const uint64_t componentId = disjointSets.find(i);
        componentsOrientedReads[componentId].push_back(OrientedReadId::fromValue(ReadId(i)));
    }

    // Gather the primary marker graph edges in each connected component.
    // This stores PrimaryIds, not MarkerGraphEdgeIds.
    vector< vector<uint64_t> > componentsPrimaryIds(orientedReadCount);
    for(uint64_t primaryId=0; primaryId<primaryMarkerGraphEdgeIds.size(); primaryId++) {
        const MarkerGraphEdgeId edgeId = primaryMarkerGraphEdgeIds[primaryId];
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        SHASTA_ASSERT(not markerIntervals.empty());
        const OrientedReadId orientedReadId0 = markerIntervals.front().orientedReadId;
        const uint64_t componentId = disjointSets.find(orientedReadId0.getValue());

        // Check that all MarkerIntervals are in the same component.
        // THIS CHECK CAN BE REMOVED FOR PERFORMANCE.
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId1 = markerInterval.orientedReadId;
            SHASTA_ASSERT(disjointSets.find(orientedReadId1.getValue()) == componentId);
        }
        componentsPrimaryIds[componentId].push_back(primaryId);
    }



    disjointSetsData.clear();



    // Gather the components with more than one read and their sizes.
    // The connected components cannot be self-complementary because
    // we are using read strand separation method 2.
    // This means that the ReadIds must be all distinct (and increasing).
    // For each complementary pair, only keep the one
    // that has the first oriented read on strand 0.
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(uint64_t componentId=0; componentId<orientedReadCount; componentId++) {
        const vector<OrientedReadId>& component = componentsOrientedReads[componentId];
        const uint64_t componentSize = component.size();
        if(componentSize < 2) {
            continue;
        }
        if(component.front().getStrand() != 0) {
            continue;
        }

        // Verify that the ReadIds are all distinct.
        // THIS CHECK CAN BE REMOVED FOR PERFORMANCE.
        for(uint64_t i1=1; i1<component.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            SHASTA_ASSERT(component[i0].getReadId() < component[i1].getReadId());
        }

        // Store this component in the componentTable.
        componentTable.push_back({componentId, componentSize});
    }

    // Sort the component table by decreasing size.
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Store the connected components we kept.
    connectedComponents.resize(componentTable.size());
    for(uint64_t i=0; i<connectedComponents.size(); i++) {
        const uint64_t componentId = componentTable[i].first;
        connectedComponents[i].orientedReadIds.swap(componentsOrientedReads[componentId]);
        connectedComponents[i].primaryIds.swap(componentsPrimaryIds[componentId]);
    }

    // Summarize in a csv file the connected components we kept.
    {
        ofstream csv("Mode3Assembler-OrientedReadIdComponents.csv");
        csv << "Rank,OrientedReads,First OrientedReadId,Primary edges\n";
        for(uint64_t i=0; i<connectedComponents.size(); i++) {
            const ConnectedComponent& connectedComponent = connectedComponents[i];
            csv << i << "," << connectedComponent.orientedReadIds.size() << "," <<
                connectedComponent.orientedReadIds.front() << "," <<
                connectedComponent.primaryIds.size() << "\n";
        }
    }


    performanceLog << timestamp << "Mode3Assembler::computeConnectedComponents ends." << endl;
}

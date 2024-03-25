// Shasta.
#include "Mode3Assembler.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "dset64-gccAtomic.hpp"
#include "mode3-PrimaryGraph.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

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
    assembleConnectedComponents();

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

    // Fill in the orientedReadIdTable.
    orientedReadIdTable.clear();
    orientedReadIdTable.resize(orientedReadCount, {invalid<uint64_t>, invalid<uint64_t>});
    for(uint64_t componentId=0; componentId<connectedComponents.size(); componentId++) {
        const vector<OrientedReadId>& orientedReadIds = connectedComponents[componentId].orientedReadIds;
        for(uint64_t position=0; position<orientedReadIds.size(); position++) {
            const OrientedReadId orientedReadId = orientedReadIds[position];
            orientedReadIdTable[orientedReadId.getValue()] = {componentId, position};
        }
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



void Mode3Assembler::assembleConnectedComponents()
{
    performanceLog << timestamp << "Mode3Assembler::assembleConnectedComponents begins." << endl;
    for(uint64_t componentId=0; componentId<connectedComponents.size(); componentId++) {
        assembleConnectedComponent(componentId);
    }
    performanceLog << timestamp << "Mode3Assembler::assembleConnectedComponents ends." << endl;
}



void Mode3Assembler::assembleConnectedComponent(uint64_t componentId)
{
    performanceLog << timestamp << "Assembling connected component " <<
        componentId << " of " << connectedComponents.size() << endl;
    cout << timestamp << "Assembling connected component " <<
        componentId << " of " << connectedComponents.size() << endl;

    ConnectedComponent connectedComponent = connectedComponents[componentId];
    const vector<OrientedReadId>& orientedReadIds = connectedComponent.orientedReadIds;
    const vector<uint64_t>& primaryIds = connectedComponent.primaryIds;

    cout << "This connected component has " << orientedReadIds.size() <<
        " reads and " << primaryIds.size() << " primary marker graph edges." << endl;



    // We need to compute the primary journey of each oriented read,
    // that is, the sequence of primary edges encountered by each read.
    // We store each journey as a vector of pairs of
    // (ordinal0, localPrimaryId), where localPrimaryId is an index into primaryIds
    // for this connected component.
    vector< vector< pair<uint32_t, uint64_t> > > journeys(orientedReadIds.size());

    performanceLog << timestamp << "Journey computation begins." << endl;
    for(uint64_t localPrimaryId=0; localPrimaryId<primaryIds.size(); localPrimaryId++) {
        const uint64_t primaryId = primaryIds[localPrimaryId];
        const MarkerGraphEdgeId edgeId = primaryMarkerGraphEdgeIds[primaryId];
        const auto markerIntervals = assembler.markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            const auto& p = orientedReadIdTable[orientedReadId.getValue()];
            SHASTA_ASSERT(p.first == componentId);
            journeys[p.second].push_back({ordinal0, localPrimaryId});
        }
    }
    for(vector< pair<uint32_t, uint64_t> >& journey: journeys) {
        sort(journey.begin(), journey.end(), OrderPairsByFirstOnly<uint32_t, uint64_t>());
    }
    performanceLog << timestamp << "Journey computation ends." << endl;

    // Check that the journeys computed in this way are identical to the ones stored in the MarkerGraph.
    // The ones stored in the MarkerGraph will eventually go away.
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto journey = journeys[i];
        const auto storedJourney = assembler.markerGraph.primaryJourneys[orientedReadId.getValue()];
        SHASTA_ASSERT(journey.size() == storedJourney.size());

        for(uint64_t j=0; j<journey.size(); j++) {
            const auto& p = journey[j];
            const uint64_t localPrimaryId = p.second;
            const uint64_t primaryId = primaryIds[localPrimaryId];
            const MarkerGraphEdgeId edgeId = primaryMarkerGraphEdgeIds[primaryId];
            // cout << orientedReadId << " " << storedJourney[j].edgeId << " " << edgeId << endl;
            SHASTA_ASSERT(edgeId == storedJourney[j].edgeId);
        }
    }



    // Now we can create the PrimaryGraph for this connected component.
    PrimaryGraph primaryGraph;

    // Create the vertices first.
    vector<PrimaryGraph::vertex_descriptor> vertexDescriptors;
    for(uint64_t localPrimaryId=0; localPrimaryId<primaryIds.size(); localPrimaryId++) {
        const uint64_t primaryId = primaryIds[localPrimaryId];
        const MarkerGraphEdgeId edgeId = primaryMarkerGraphEdgeIds[primaryId];
        vertexDescriptors.push_back(primaryGraph.addVertex(edgeId));
    }



    // To generate edges of the PrimaryGraph, we need to gather pairs of consecutive
    // journey entries. Each pair (localPrimaryId0, localPrimaryId1) is stored
    // as a localPrimaryId1 in journeyPairs[localPrimaryId0].
    // For now use a simple vector of vector and sequential code, but later
    // switch to MemoryMapped::VectorOfVectors<uint64_t, uint64_t> and multithreaded code.
    vector< vector<uint64_t> > journeyPairs(primaryIds.size());
    performanceLog << timestamp << "PrimaryGraph edge creation begins begins." << endl;
    for(const auto& journey: journeys) {
        for(uint64_t i1=1; i1<journey.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t localPrimaryId0 = journey[i0].second;
            const uint64_t localPrimaryId1 = journey[i1].second;
            journeyPairs[localPrimaryId0].push_back(localPrimaryId1);
        }
     }
     vector<uint64_t> count;
     for(uint64_t localPrimaryId0=0; localPrimaryId0<primaryIds.size(); localPrimaryId0++) {
         const PrimaryGraph::vertex_descriptor v0 = vertexDescriptors[localPrimaryId0];
         const MarkerGraphEdgeId edgeId0 = primaryGraph[v0].edgeId;
         auto journeyPairs0 = journeyPairs[localPrimaryId0];
         deduplicateAndCount(journeyPairs0, count);
         SHASTA_ASSERT(journeyPairs0.size() == count.size());
         for(uint64_t j=0; j<journeyPairs0.size(); j++) {
             const uint64_t localPrimaryId1 = journeyPairs0[j];
             const uint64_t coverage = count[j];
             const PrimaryGraph::vertex_descriptor v1 = vertexDescriptors[localPrimaryId1];
             const MarkerGraphEdgeId edgeId1 = primaryGraph[v1].edgeId;
             MarkerGraphEdgePairInfo info;
             SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
             primaryGraph.addEdgeFromVertexDescriptors(v0, v1, info, coverage);
         }
     }
     performanceLog << timestamp << "PrimaryGraph edge creation ends." << endl;

     cout << "The PrimaryGraph for this connected component has " <<
         num_vertices(primaryGraph) << " and " << num_edges(primaryGraph) << " edges." << endl;
}

// Shasta.
#include "Mode3Assembler.hpp"
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "deduplicate.hpp"
#include "dset64-gccAtomic.hpp"
#include "mode3-AssemblyGraph.hpp"
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
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool debug) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(assembler),
    assembler(assembler),
    debug(debug),
    anchors(assembler.markerGraph.edgeMarkerIntervals)
{
    performanceLog << timestamp << "Mode 3 assembly begins." << endl;

    computeConnectedComponents();
    if(debug) {
        writeConnectedComponents();
    }
    assembleConnectedComponents(threadCount, options, debug);

    performanceLog << timestamp << "Mode 3 assembly ends." << endl;
}



// The oriented reads present in each anchor (primary marker graph edge)
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

    // Loop over all marker graph edges (that is, over all anchors).
    // This could be multithreaded but runs at decent speed as is.
    for(MarkerGraphEdgeId edgeId=0; edgeId<anchors.size(); edgeId++) {
        const auto markerIntervals = anchors[edgeId];
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

    // Gather the anchors (marker graph edges) in each connected component.
    vector< vector<uint64_t> > componentsMarkerGraphEdgeIds(orientedReadCount);
    for(MarkerGraphEdgeId edgeId=0; edgeId<anchors.size(); edgeId++) {
        const auto markerIntervals = anchors[edgeId];
        SHASTA_ASSERT(not markerIntervals.empty());
        const OrientedReadId orientedReadId0 = markerIntervals.front().orientedReadId;
        const uint64_t componentId = disjointSets.find(orientedReadId0.getValue());

        // Check that all MarkerIntervals are in the same component.
        // THIS CHECK CAN BE REMOVED FOR PERFORMANCE.
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId1 = markerInterval.orientedReadId;
            SHASTA_ASSERT(disjointSets.find(orientedReadId1.getValue()) == componentId);
        }
        componentsMarkerGraphEdgeIds[componentId].push_back(edgeId);
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
        connectedComponents[i].markerGraphEdgeIds.swap(componentsMarkerGraphEdgeIds[componentId]);
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

    performanceLog << timestamp << "Mode3Assembler::computeConnectedComponents ends." << endl;
}



void Mode3Assembler::assembleConnectedComponents(
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool debug)
{
    performanceLog << timestamp << "Mode3Assembler::assembleConnectedComponents begins." << endl;

    vector< vector<uint64_t> > assemblyChainLengthsByPValue;
    vector<uint64_t> assemblyBubbleChainLengths;

    ofstream summaryCsv("Components.csv");
    summaryCsv << "Component,Reads,Segments,Sequence,N50,Total Bubble chain length,Bubble chain N50\n";

    vector< shared_ptr<mode3::AssemblyGraph> > assemblyGraphs;
    for(uint64_t componentId=0; componentId<connectedComponents.size(); componentId++) {
        const shared_ptr<AssemblyGraph> assemblyGraph =
            assembleConnectedComponent(componentId, threadCount, options, true, debug);
        assemblyGraphs.push_back(assemblyGraph);

        // Chain length statistics.
        vector< vector<uint64_t> > chainLengths;
        assemblyGraph->getChainLengthsByPValue(chainLengths);

        // Assembly statistics by P-value.
        cout << "Assembly statistics by P-Value for component " << componentId << ":" << endl;
        for(uint64_t pValue=0; pValue<chainLengths.size(); pValue++) {
            uint64_t totalLength, n50;
            tie(totalLength, n50) = AssemblyGraph::n50(chainLengths[pValue]);
            cout << "P-value " << pValue << ": total assembled length " << totalLength <<
                ", N50 " << n50 << endl;
        }

        // Combined chain length statistics for this component.
        vector<uint64_t> allChainLengths;
        for(const auto& v: chainLengths) {
            copy(v.begin(), v.end(), back_inserter(allChainLengths));
        }
        sort(allChainLengths.begin(), allChainLengths.end(), std::greater<uint64_t>());
        uint64_t totalLength, n50;
        tie(totalLength, n50) = AssemblyGraph::n50(allChainLengths);
        cout << "Combined for this component: total assembled length " << totalLength <<
            ", N50 " << n50 << endl;

        // Bubble chain length statistics (non-trivial bubble chains only).
        vector<uint64_t> bubbleChainLengths;
        assemblyGraph->getBubbleChainLengths(bubbleChainLengths);
        uint64_t totalBubbleChainLength, bubbleChainN50;
        tie(totalBubbleChainLength, bubbleChainN50) = AssemblyGraph::n50(bubbleChainLengths);
        copy(bubbleChainLengths.begin(), bubbleChainLengths.end(),
            back_inserter(assemblyBubbleChainLengths));
        cout << "Total non-trivial bubble chain length for this component " << totalBubbleChainLength <<
            ", N50 " << bubbleChainN50 << endl;

        // Write a line to the summaryCsv.
        summaryCsv << componentId << ",";
        summaryCsv << connectedComponents[componentId].orientedReadIds.size() << ",";
        summaryCsv << allChainLengths.size() << ",";
        summaryCsv << totalLength << ",";
        summaryCsv << n50 << ",";
        summaryCsv << totalBubbleChainLength << ",";
        summaryCsv << bubbleChainN50 << "\n";

        // Store the chain lengths.
        if(assemblyChainLengthsByPValue.size() < chainLengths.size()) {
            assemblyChainLengthsByPValue.resize(chainLengths.size());
        }
        for(uint64_t pValue=0; pValue<chainLengths.size(); pValue++) {
            copy(chainLengths[pValue].begin(), chainLengths[pValue].end(),
                back_inserter(assemblyChainLengthsByPValue[pValue]));
        }
    }

    cout << "Global assembly statistics by P-Value:" << endl;
    for(uint64_t pValue=0; pValue<assemblyChainLengthsByPValue.size(); pValue++) {
        sort(assemblyChainLengthsByPValue[pValue].begin(), assemblyChainLengthsByPValue[pValue].end(),
            std::greater<uint64_t>());
        uint64_t totalLength, n50;
        tie(totalLength, n50) = AssemblyGraph::n50(assemblyChainLengthsByPValue[pValue]);
        cout << "P-value " << pValue << ": total assembled length " << totalLength <<
            ", N50 " << n50 << endl;
    }
    vector<uint64_t> allChainLengths;
    for(const auto& v: assemblyChainLengthsByPValue) {
        copy(v.begin(), v.end(), back_inserter(allChainLengths));
    }
    sort(allChainLengths.begin(), allChainLengths.end(), std::greater<uint64_t>());
    uint64_t totalLength, n50;
    tie(totalLength, n50) = AssemblyGraph::n50(allChainLengths);
    cout << "Global assembly statistics, combined for all P-values: total assembled length " << totalLength <<
        ", N50 " << n50 << endl;

    sort(assemblyBubbleChainLengths.begin(), assemblyBubbleChainLengths.end(), std::greater<uint64_t>());
    uint64_t totalBubbleChainLength, bubbleChainN50;
    tie(totalBubbleChainLength, bubbleChainN50) = AssemblyGraph::n50(assemblyBubbleChainLengths);
    cout << "Total non-trivial bubble chain length " << totalBubbleChainLength <<
        ", N50 " << bubbleChainN50 << endl;


    // Create a csv file with one line for each assembled segment.
    // This can also be loaded in Bandage.
    {
        ofstream csv("Assembly.csv");
        csv << "Chain,Connectivity,Component,Bubble chain,Position in bubble chain,Index in bubble,"
            "Sequence length,Primary coverage,P value,Color,"
            "Preceded by,Followed by,"
            "\n";
        for(const shared_ptr<mode3::AssemblyGraph>& assemblyGraph: assemblyGraphs) {;
            assemblyGraph->writeCsvSummary(csv);
        }
    }

    // Create a global FASTA file with output from all the connected components.
    {
        ofstream fasta("Assembly.fasta");
        for(const shared_ptr<mode3::AssemblyGraph>& assemblyGraph: assemblyGraphs) {
            assemblyGraph->writeFastaExpanded(fasta);
        }
    }

    // Create a global GFA file with output from all the connected components.
    {
        ofstream gfa("Assembly.gfa");
        AssemblyGraph::writeGfaHeader(gfa);
        for(const shared_ptr<mode3::AssemblyGraph>& assemblyGraph: assemblyGraphs) {
            assemblyGraph->writeGfaSegmentsExpanded(gfa, true, true);
        }
        for(const shared_ptr<mode3::AssemblyGraph>& assemblyGraph: assemblyGraphs) {
            assemblyGraph->writeGfaLinksExpanded(gfa);
        }
    }

    // Also create a global GFA file without sequence.
    ofstream gfa("Assembly-NoSequence.gfa");
    {
        AssemblyGraph::writeGfaHeader(gfa);
        for(const shared_ptr<mode3::AssemblyGraph>& assemblyGraph: assemblyGraphs) {
            assemblyGraph->writeGfaSegmentsExpanded(gfa, false, true);
        }
        for(const shared_ptr<mode3::AssemblyGraph>& assemblyGraph: assemblyGraphs) {
            assemblyGraph->writeGfaLinksExpanded(gfa);
        }
    }

    performanceLog << timestamp << "Mode3Assembler::assembleConnectedComponents ends." << endl;
}



shared_ptr<AssemblyGraph> Mode3Assembler::assembleConnectedComponent(
    uint64_t componentId,
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool assembleSequence,
    bool debug)
{
    performanceLog << timestamp << "Assembling connected component " <<
        componentId << " of " << connectedComponents.size() << endl;
    cout << timestamp << "Assembling connected component " <<
        componentId << " of " << connectedComponents.size() << endl;

    const ConnectedComponent& connectedComponent = connectedComponents[componentId];
    const vector<OrientedReadId>& orientedReadIds = connectedComponent.orientedReadIds;
    const vector<uint64_t>& markerGraphEdgeIds = connectedComponent.markerGraphEdgeIds;

    cout << "This connected component has " << orientedReadIds.size() <<
        " oriented reads and " << markerGraphEdgeIds.size() << " anchors." << endl;



    // We need to compute the anchor journey of each oriented read,
    // that is, the sequence of anchors encountered by each read.
    // We store each journey as a vector of pairs of
    // (ordinal0, localAnchorId), where localAnchorId is an index into anchorIds (markerGraphEdgeIds)
    // for this connected component.
    vector< vector< pair<uint32_t, uint64_t> > > journeys(orientedReadIds.size());

    performanceLog << timestamp << "Journey computation begins." << endl;
    for(uint64_t localAnchorId=0; localAnchorId<markerGraphEdgeIds.size(); localAnchorId++) {
        const MarkerGraphEdgeId edgeId = markerGraphEdgeIds[localAnchorId];
        const auto anchor = anchors[edgeId];
        for(const MarkerInterval& markerInterval: anchor) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            const auto& p = orientedReadIdTable[orientedReadId.getValue()];
            SHASTA_ASSERT(p.first == componentId);
            journeys[p.second].push_back({ordinal0, localAnchorId});
        }
    }
    for(vector< pair<uint32_t, uint64_t> >& journey: journeys) {
        sort(journey.begin(), journey.end(), OrderPairsByFirstOnly<uint32_t, uint64_t>());
    }
    performanceLog << timestamp << "Journey computation ends." << endl;



    // Write out the journeys.
    if(debug) {
        ofstream csv("Journeys-" + to_string(componentId) + ".csv");
        for(uint64_t i=0; i<orientedReadIds.size(); i++) {
            csv << orientedReadIds[i] << ",";
            const auto& journey = journeys[i];
            for(const auto& p: journey) {
                const uint64_t localAnchorId = p.second;
                const MarkerGraphEdgeId edgeId = markerGraphEdgeIds[localAnchorId];
                csv << edgeId << ",";
            }
            csv << "\n";
        }
    }



    // Now we can create the PrimaryGraph for this connected component.
    PrimaryGraph primaryGraph;

    // Create the vertices first.
    vector<PrimaryGraph::vertex_descriptor> vertexDescriptors;
    for(uint64_t localPrimaryId=0; localPrimaryId<markerGraphEdgeIds.size(); localPrimaryId++) {
        const MarkerGraphEdgeId edgeId = markerGraphEdgeIds[localPrimaryId];
        vertexDescriptors.push_back(primaryGraph.addVertex(edgeId));
    }



    // To generate edges of the PrimaryGraph, we need to gather pairs of consecutive
    // journey entries. Each pair (localPrimaryId0, localPrimaryId1) is stored
    // as a localPrimaryId1 in journeyPairs[localPrimaryId0].
    // For now use a simple vector of vector and sequential code, but later
    // switch to MemoryMapped::VectorOfVectors<uint64_t, uint64_t> and multithreaded code.
    vector< vector<uint64_t> > journeyPairs(markerGraphEdgeIds.size());
    performanceLog << timestamp << "PrimaryGraph edge creation begins." << endl;
    for(const auto& journey: journeys) {
        for(uint64_t i1=1; i1<journey.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t localPrimaryId0 = journey[i0].second;
            const uint64_t localPrimaryId1 = journey[i1].second;
            journeyPairs[localPrimaryId0].push_back(localPrimaryId1);
        }
     }
     vector<uint64_t> count;
     for(uint64_t localPrimaryId0=0; localPrimaryId0<markerGraphEdgeIds.size(); localPrimaryId0++) {
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
         num_vertices(primaryGraph) << " vertices and " << num_edges(primaryGraph) << " edges." << endl;


     // Graphviz output.
     if(debug) {
         PrimaryGraphDisplayOptions options;
         options.showNonTransitiveReductionEdges = true;
         primaryGraph.writeGraphviz(
             "PrimaryGraphInitial" + to_string(componentId), options, assembler.markerGraph);
         options.makeCompact();
         primaryGraph.writeGraphviz(
             "PrimaryGraphCompactInitial" + to_string(componentId), options, assembler.markerGraph);
         primaryGraph.writeEdgeCoverageHistogram("PrimaryGraphInitial" + to_string(componentId) + "-EdgeCoverageHistogram.csv");
     }

     // Remove weak edges..
     primaryGraph.removeWeakEdges(options.primaryGraphOptions.maxLoss, debug);

     // Remove cross-edges.
     primaryGraph.removeCrossEdges(
         options.primaryGraphOptions.crossEdgesLowCoverageThreshold,
         options.primaryGraphOptions.crossEdgesHighCoverageThreshold,
         0,
         debug);

     // Graphviz output.
     if(debug) {
         PrimaryGraphDisplayOptions options;
         options.showNonTransitiveReductionEdges = false;
         primaryGraph.writeGraphviz(
             "PrimaryGraph" + to_string(componentId), options, assembler.markerGraph);
         options.makeCompact();
         primaryGraph.writeGraphviz(
             "PrimaryGraphCompact" + to_string(componentId), options, assembler.markerGraph);
     }

     // Create the assembly graph for this connected component.
     return make_shared<AssemblyGraph>(
         primaryGraph, componentId, assembler, orientedReadIds, markerGraphEdgeIds, threadCount,
         options, assembleSequence, debug);
}



// Debug output of all connected components.
void Mode3Assembler::writeConnectedComponents() const
{
    for(uint64_t componentId=0; componentId<connectedComponents.size(); componentId++) {
        writeConnectedComponent(componentId);
    }
}



// Debug output of a connected components.
void Mode3Assembler::writeConnectedComponent(uint64_t componentId) const
{
    const ConnectedComponent& component = connectedComponents[componentId];
    const vector<OrientedReadId>& orientedReadIds = component.orientedReadIds;
    const vector<uint64_t>& markerGraphEdgeIds = component.markerGraphEdgeIds;

    // Write the oriented reads.
    {
        ofstream csv("OrientedReadIds-" + to_string(componentId) + ".csv");
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            csv << orientedReadId << "\n";
        }
    }



    // Write the marker intervals.
    {
        ofstream csv("MarkerIntervals-" + to_string(componentId) + ".csv");
        csv << "AnchorId,OrientedReadId,Ordinal0,Ordinal1\n";

        for(uint64_t localPrimaryId=0; localPrimaryId<component.markerGraphEdgeIds.size(); localPrimaryId++) {
            const MarkerGraphEdgeId edgeId = markerGraphEdgeIds[localPrimaryId];
            const auto anchor = anchors[edgeId];
            for(const MarkerInterval& markerInterval: anchor) {
                csv << edgeId << ",";
                csv << markerInterval.orientedReadId << ",";
                csv << markerInterval.ordinals[0] << ",";
                csv << markerInterval.ordinals[1] << "\n";
            }
        }

    }
}


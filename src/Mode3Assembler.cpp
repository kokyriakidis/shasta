// Shasta.
#include "Mode3Assembler.hpp"
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "deduplicate.hpp"
#include "dset64-gccAtomic.hpp"
#include "mode3-AssemblyGraph.hpp"
#include "mode3-AnchorGraph.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
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
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    shared_ptr<mode3::Anchors> anchorsPointer,
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool debug) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(assembler),
    assembler(assembler),
    k(k),
    reads(reads),
    markers(markers),
    debug(debug),
    anchorsPointer(anchorsPointer)
{
    SHASTA_ASSERT(anchorsPointer);

    performanceLog << timestamp << "Mode 3 assembly begins." << endl;
    findReverseComplementAnchors();

    computeConnectedComponents();
    if(debug) {
        writeConnectedComponents();
    }
    assembleConnectedComponents(threadCount, options, debug);

    performanceLog << timestamp << "Mode 3 assembly ends." << endl;
}



Mode3Assembler::Mode3Assembler(
    const Assembler& assembler,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    shared_ptr<mode3::Anchors> anchorsPointer) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(assembler),
    assembler(assembler),
    k(k),
    reads(reads),
    markers(markers),
    anchorsPointer(anchorsPointer)
{
    SHASTA_ASSERT(anchorsPointer);
}



void Mode3Assembler::findReverseComplementAnchors()
{
    const uint64_t readCount = assembler.markers.size() / 2;

    // For each Anchor, we look at the lowest numbered ReadId,
    // and store the ReadId and ordinal in anchorInfos[strand][readId].
    class AnchorInfo {
    public:
        AnchorId anchorId;
        uint32_t ordinal;
        bool operator<(const AnchorInfo& that) const
        {
            return ordinal < that.ordinal;
        }
    };
    array<vector< vector<AnchorInfo> >, 2> anchorInfos;
    for(uint64_t i=0; i<2; i++) {
        anchorInfos[i].resize(readCount);
    }
    for(AnchorId anchorId=0; anchorId<anchors().size(); anchorId++) {
        const MarkerInterval& markerInterval = anchors()[anchorId].front();
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const uint32_t ordinal0 = markerInterval.ordinals[0];
        anchorInfos[strand][readId].push_back({anchorId, ordinal0});
    }

    // Sort by ordinal the AnchorInfos for each Strand and Read.
    for(Strand strand=0; strand<2; strand++) {
        for(ReadId readId=0; readId<readCount; readId++) {
            vector<AnchorInfo>& v = anchorInfos[strand][readId];
            sort(v.begin(), v.end());
        }
    }



    // Now we can find the reverse complement of each Anchor
    // by scanning the anchorInfos for the two OrientedReadId of each ReadId.
    reverseComplementAnchor.resize(anchors().size(), invalid<AnchorId>);
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        const uint64_t markerCount = assembler.markers.size(orientedReadId0.getValue());

        const vector<AnchorInfo>& anchorInfos0 = anchorInfos[0][readId];
        const vector<AnchorInfo>& anchorInfos1 = anchorInfos[1][readId];
        const uint64_t n = anchorInfos0.size();
        SHASTA_ASSERT(anchorInfos1.size() == n);

        for(uint64_t i0=0; i0<n; i0++) {
            const uint64_t i1 = n - i0 - 1;

            const AnchorInfo& anchorInfo0 = anchorInfos0[i0];
            const AnchorInfo& anchorInfo1 = anchorInfos1[i1];

            const AnchorId anchorId0 = anchorInfo0.anchorId;
            const AnchorId anchorId1 = anchorInfo1.anchorId;

            const Anchor& anchor0 = anchors()[anchorId0];
            const Anchor& anchor1 = anchors()[anchorId1];

            const MarkerInterval& markerInterval0 = anchor0.front();
            const MarkerInterval& markerInterval1 = anchor1.front();

            SHASTA_ASSERT(markerInterval0.orientedReadId == orientedReadId0);
            SHASTA_ASSERT(markerInterval1.orientedReadId == orientedReadId1);

            SHASTA_ASSERT(markerInterval0.ordinals[0] + markerInterval1.ordinals[1] == markerCount - 1);

            reverseComplementAnchor[anchorId0] = anchorId1;
            reverseComplementAnchor[anchorId1] = anchorId0;
        }
    }

    // Sanity check.
    for(AnchorId anchorId=0; anchorId<anchors().size(); anchorId++) {
        const AnchorId anchorIdRc = reverseComplementAnchor[anchorId];
        SHASTA_ASSERT(anchorIdRc != invalid<AnchorId>);
        SHASTA_ASSERT(anchorIdRc != anchorId);
        SHASTA_ASSERT(reverseComplementAnchor[anchorIdRc] == anchorId);
    }
}



// The oriented reads present in each anchor (primary marker graph edge)
// define a bipartite graph. We want to compute connected components
// of this bipartite graph and process them one at a time.
void Mode3Assembler::computeConnectedComponents()
{
    performanceLog << timestamp << "Mode3Assembler::computeConnectedComponents begins." << endl;

    // Compute connected components of the oriented reads portion
    // of the bipartite graph.
    // OrientedReadIds are indexed by OrientedReadId::getValue().
    const uint64_t orientedReadCount = assembler.markers.size();
    vector<DisjointSets::Aint> disjointSetsData(orientedReadCount);
    DisjointSets disjointSets(&disjointSetsData[0], orientedReadCount);

    // Loop over all Anchors.
    // All oriented reads in an Anchor belong to the same connected component.
    // This could be multithreaded but runs at decent speed as is.
    for(AnchorId anchorId=0; anchorId<anchors().size(); anchorId++) {
        const Anchor anchor = anchors()[anchorId];
        SHASTA_ASSERT(not anchor.empty());
        const OrientedReadId orientedReadId0 = anchor.front().orientedReadId;
        for(const MarkerInterval& markerInterval: anchor) {
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
    for(MarkerGraphEdgeId edgeId=0; edgeId<anchors().size(); edgeId++) {
        const auto markerIntervals = anchors()[edgeId];
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
    // Each component can be self-complementary (contains both strands)
    // or part of a pairs of reverse complementary component.
    // We keep all self-complementary components.
    // For each complementary pair, only keep the one
    // that has the first oriented read on strand 0.
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(uint64_t componentId=0; componentId<orientedReadCount; componentId++) {
        const vector<OrientedReadId>& component = componentsOrientedReads[componentId];
        if(component.size() < 2) {
            continue;
        }

        // Figure out if this is a self-complementary component (contains both strands).
        const bool isSelfComplementary =  component[0].getReadId() == component[1].getReadId();

        if(isSelfComplementary or (component.front().getStrand() == 0)) {
            componentTable.push_back({componentId, component.size()});
        }
    }

    // Sort the component table by decreasing size.
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Store the connected components we kept.
    connectedComponents.resize(componentTable.size());
    for(uint64_t i=0; i<connectedComponents.size(); i++) {
        const uint64_t componentId = componentTable[i].first;
        connectedComponents[i].orientedReadIds.swap(componentsOrientedReads[componentId]);
        connectedComponents[i].anchorIds.swap(componentsMarkerGraphEdgeIds[componentId]);
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

    ofstream orientedReadsCsv("OrientedReadsByComponent.csv");
    orientedReadsCsv << "Component,OrientedReadId,ReadName\n";

    vector< shared_ptr<mode3::AssemblyGraph> > assemblyGraphs;
    for(uint64_t componentId=0; componentId<connectedComponents.size(); componentId++) {
        const shared_ptr<AssemblyGraph> assemblyGraph =
            assembleConnectedComponent(componentId, threadCount, options, true, orientedReadsCsv, debug);
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




bool Mode3Assembler::ConnectedComponent::isSelfComplementary() const
{
    // A self-complementary component must have at least two oriented reads.
    if(orientedReadIds.size() < 2) {
        return false;
    }

    // In a self-complementary component, the number of oriented reads must be even.
    if((orientedReadIds.size() % 2) == 1) {
        return false;
    }

    // Otherwise, check the first two OrientedReadId.
    return orientedReadIds[0].getReadId() == orientedReadIds[1].getReadId();
}



void Mode3Assembler::ConnectedComponent::checkIsValid() const
{
    if(isSelfComplementary()) {

        // For a self-complementary component, all oriented reads come in reverse complemented pairs.
        SHASTA_ASSERT((orientedReadIds.size() %2) == 0);
        for(uint64_t i=0; i<orientedReadIds.size(); i+=2) {
            const OrientedReadId orientedReadId0 = orientedReadIds[i];
            const OrientedReadId orientedReadId1 = orientedReadIds[i+1];
            SHASTA_ASSERT(orientedReadId0.getReadId() == orientedReadId1.getReadId());
            SHASTA_ASSERT(orientedReadId0.getStrand() == 0);
            SHASTA_ASSERT(orientedReadId1.getStrand() == 1);
        }

    } else {

        // For a non-self-complementary component, all ReadIds must be distinct.
        for(uint64_t i1=1; i1<orientedReadIds.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            SHASTA_ASSERT(orientedReadIds[i0].getReadId() < orientedReadIds[i1].getReadId());
        }

    }

}



shared_ptr<AssemblyGraph> Mode3Assembler::assembleConnectedComponent(
    uint64_t componentId,
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool assembleSequence,
    ostream& orientedReadsCsv,
    bool debug)
{
    performanceLog << timestamp << "Assembling connected component " <<
        componentId << " of " << connectedComponents.size() << endl;
    cout << timestamp << "Assembling connected component " <<
        componentId << " of " << connectedComponents.size() << endl;

    const ConnectedComponent& connectedComponent = connectedComponents[componentId];
    const vector<OrientedReadId>& orientedReadIds = connectedComponent.orientedReadIds;
    const vector<uint64_t>& anchorIds = connectedComponent.anchorIds;

    const bool isSelfComplementary = connectedComponent.isSelfComplementary();
    connectedComponent.checkIsValid();

    cout << "This connected component has " << orientedReadIds.size() <<
        " oriented reads and " << anchorIds.size() << " anchors." << endl;
    if(isSelfComplementary) {
        cout << "This connected component is self-complementary and will be assembled double-stranded." << endl;
    }

    // Write to orientedReadsCsv the oriented reads for this component.
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        const auto readName = assembler.getReads().getReadName(orientedReadId.getReadId());
        orientedReadsCsv <<
            componentId << "," <<
            orientedReadId << ",";
        copy(readName.begin(), readName.end(), ostream_iterator<char>(orientedReadsCsv));
        orientedReadsCsv << "\n";
    }


    // We need to compute the anchor journey of each oriented read,
    // that is, the sequence of anchors encountered by each read.
    // We store each journey as a vector of pairs of
    // (ordinal0, localAnchorId), where localAnchorId is an index into anchorIds (markerGraphEdgeIds)
    // for this connected component.
    vector< vector< pair<uint32_t, uint64_t> > > journeys(orientedReadIds.size());

    performanceLog << timestamp << "Journey computation begins." << endl;
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        const AnchorId anchorId = anchorIds[localAnchorId];
        const Anchor anchor = anchors()[anchorId];
        for(const MarkerInterval& markerInterval: anchor) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            const auto& p = orientedReadIdTable[orientedReadId.getValue()];
            if(p.first == componentId) { // Due to StrandSplitter this is not always the case.
                journeys[p.second].push_back({ordinal0, localAnchorId});
            }
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
                const AnchorId anchorId = anchorIds[localAnchorId];
                csv << anchorId << ",";
            }
            csv << "\n";
        }
    }

    // Now we can create the AnchorGraph for this connected component.
    // The constructor generates the vertices.
    AnchorGraph anchorGraph(anchorIds);

    // To generate edges of the AnchorGraph, we need to gather pairs of consecutive
    // journey entries. Each pair (localAnchorId0, localAnchorId1) is stored
    // as a localAnchorId1 in journeyPairs[localAnchorId0].
    // For now use a simple vector of vector and sequential code, but later
    // switch to MemoryMapped::VectorOfVectors<uint64_t, uint64_t> and multithreaded code.
    vector< vector<uint64_t> > journeyPairs(anchorIds.size());
    performanceLog << timestamp << "AnchorGraph edge creation begins." << endl;
    for(const auto& journey: journeys) {
        for(uint64_t i1=1; i1<journey.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t localAnchorId0 = journey[i0].second;
            const uint64_t localAnchorId1 = journey[i1].second;
            journeyPairs[localAnchorId0].push_back(localAnchorId1);
        }
     }
     vector<uint64_t> count;
     for(uint64_t localAnchorId0=0; localAnchorId0<anchorIds.size(); localAnchorId0++) {
         const AnchorId anchorId0 = anchorIds[localAnchorId0];
         auto journeyPairs0 = journeyPairs[localAnchorId0];
         deduplicateAndCount(journeyPairs0, count);
         SHASTA_ASSERT(journeyPairs0.size() == count.size());
         for(uint64_t j=0; j<journeyPairs0.size(); j++) {
             const uint64_t localAnchorId1 = journeyPairs0[j];
             const uint64_t coverage = count[j];
             const AnchorId anchorId1 = anchorIds[localAnchorId1];
             MarkerGraphEdgePairInfo info;
             SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(anchorId0, anchorId1, info));
             anchorGraph.addEdgeFromLocalAnchorIds(localAnchorId0, localAnchorId1, info, coverage);
         }
     }
     performanceLog << timestamp << "AnchorGraph edge creation ends." << endl;

     cout << "The AnchorGraph for this connected component has " <<
         num_vertices(anchorGraph) << " vertices and " << num_edges(anchorGraph) << " edges." << endl;

     // Graphviz output.
     if(debug) {
         AnchorGraphDisplayOptions options;
         options.showNonTransitiveReductionEdges = true;
         anchorGraph.writeGraphviz(
             "AnchorGraphInitial" + to_string(componentId), options, assembler.markerGraph);
         options.makeCompact();
         anchorGraph.writeGraphviz(
             "AnchorGraphCompactInitial" + to_string(componentId), options, assembler.markerGraph);
         anchorGraph.writeEdgeCoverageHistogram("AnchorGraphInitial" + to_string(componentId) + "-EdgeCoverageHistogram.csv");
     }

     // Remove weak edges.
     anchorGraph.removeWeakEdges(options.primaryGraphOptions.maxLoss, debug);
     cout << "After removing weak edges, the AnchorGraph for this connected component has " <<
         num_vertices(anchorGraph) << " vertices and " << num_edges(anchorGraph) << " edges." << endl;

     // Remove cross-edges.
     anchorGraph.removeCrossEdges(
         options.primaryGraphOptions.crossEdgesLowCoverageThreshold,
         options.primaryGraphOptions.crossEdgesHighCoverageThreshold,
         0,
         debug);
     cout << "After removing cross-edges, the AnchorGraph for this connected component has " <<
         num_vertices(anchorGraph) << " vertices and " << num_edges(anchorGraph) << " edges." << endl;

     // Strand separation does not work well ans is skipped.
     // If the component is self-complementary, it will be assembled double-stranded.
     if(false /*isSelfComplementary*/) {
         anchorGraph.findReverseComplementAnchors(anchors(), assembler.markers);
         anchorGraph.findReverseComplementVertices();
         anchorGraph.findReverseComplementEdges();
         anchorGraph.separateStrands();
         cout << "After strand separation, the AnchorGraph for this connected component has " <<
             num_vertices(anchorGraph) << " vertices and " << num_edges(anchorGraph) << " edges." << endl;
     }

     // Graphviz output.
     if(debug) {
         AnchorGraphDisplayOptions options;
         options.showNonTransitiveReductionEdges = false;
         anchorGraph.writeGraphviz(
             "AnchorGraph" + to_string(componentId), options, assembler.markerGraph);
         options.makeCompact();
         anchorGraph.writeGraphviz(
             "AnchorGraphCompact" + to_string(componentId), options, assembler.markerGraph);
     }

     // Create the assembly graph for this connected component.
     return make_shared<AssemblyGraph>(
         anchorGraph, anchors(), componentId, assembler.assemblerInfo->k, orientedReadIds, anchorIds, threadCount,
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
    const vector<uint64_t>& anchorIds = component.anchorIds;

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

        for(uint64_t localAnchorId=0; localAnchorId<component.anchorIds.size(); localAnchorId++) {
            const AnchorId anchorId = anchorIds[localAnchorId];
            const Anchor anchor = anchors()[anchorId];
            for(const MarkerInterval& markerInterval: anchor) {
                csv << anchorId << ",";
                csv << markerInterval.orientedReadId << ",";
                csv << markerInterval.ordinals[0] << ",";
                csv << markerInterval.ordinals[1] << "\n";
            }
        }

    }
}


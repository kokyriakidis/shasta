// Shasta.
#include "Mode3Assembler.hpp"
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



// This constructor runs the assembly.
Mode3Assembler::Mode3Assembler(
    const MappedMemoryOwner& mappedMemoryOwner,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    shared_ptr<mode3::Anchors> anchorsPointer,
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool debug) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    reads(reads),
    markers(markers),
    debug(debug),
    anchorsPointer(anchorsPointer),
    options(options)
{
    SHASTA_ASSERT(anchorsPointer);
    anchors().writeCoverageHistogram();
    anchors().writeAnchorGapsByRead();

    performanceLog << timestamp << "Mode 3 assembly begins." << endl;

    computeConnectedComponents();
    if(debug) {
        writeConnectedComponents();
    }
    assembleConnectedComponents(threadCount, debug);

    performanceLog << timestamp << "Mode 3 assembly ends." << endl;
}



// This constructor just accesses binary data.
Mode3Assembler::Mode3Assembler(
    const MappedMemoryOwner& mappedMemoryOwner,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    shared_ptr<mode3::Anchors> anchorsPointer,
    const Mode3AssemblyOptions& options) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    reads(reads),
    markers(markers),
    anchorsPointer(anchorsPointer),
    options(options)
{
    SHASTA_ASSERT(anchorsPointer);
    componentOrientedReadIds.accessExistingReadOnly(largeDataName("Mode3Assembler-componentOrientedReadIds"));
    componentAnchorIds.accessExistingReadOnly(largeDataName("Mode3Assembler-componentAnchorIds"));
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
    const uint64_t orientedReadCount = markers.size();
    vector<DisjointSets::Aint> disjointSetsData(orientedReadCount);
    DisjointSets disjointSets(&disjointSetsData[0], orientedReadCount);

    // Loop over all Anchors.
    // All oriented reads in an Anchor belong to the same connected component.
    // This could be multithreaded but runs at decent speed as is.
    for(AnchorId anchorId=0; anchorId<anchors().size(); anchorId++) {
        const Anchor anchor = anchors()[anchorId];
        SHASTA_ASSERT(not anchor.empty());
        const OrientedReadId orientedReadId0 = anchor.front().orientedReadId;
        for(const AnchorMarkerInterval& markerInterval: anchor) {
            const OrientedReadId orientedReadId1 = markerInterval.orientedReadId;
            disjointSets.unite(orientedReadId0.getValue(), orientedReadId1.getValue());
        }
    }

    // Gather the oriented reads in each connected component.
    vector< vector<OrientedReadId> > orientedReads(orientedReadCount);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const uint64_t componentId = disjointSets.find(i);
        orientedReads[componentId].push_back(OrientedReadId::fromValue(ReadId(i)));
    }

    // Gather the anchors in each connected component.
    vector< vector<uint64_t> > anchorIds(orientedReadCount);
    for(AnchorId anchorId=0; anchorId<anchors().size(); anchorId++) {
        const auto markerIntervals = anchors()[anchorId];
        SHASTA_ASSERT(not markerIntervals.empty());
        const OrientedReadId orientedReadId0 = markerIntervals.front().orientedReadId;
        const uint64_t componentId = disjointSets.find(orientedReadId0.getValue());
        anchorIds[componentId].push_back(anchorId);
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
        const vector<OrientedReadId>& component = orientedReads[componentId];
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
    componentOrientedReadIds.createNew(largeDataName("Mode3Assembler-componentOrientedReadIds"), largeDataPageSize);
    componentAnchorIds.createNew(largeDataName("Mode3Assembler-componentAnchorIds"), largeDataPageSize);
    for(uint64_t i=0; i<componentTable.size(); i++) {
        const uint64_t componentId = componentTable[i].first;
        componentOrientedReadIds.appendVector(orientedReads[componentId]);
        componentAnchorIds.appendVector(anchorIds[componentId]);
    }

    // Fill in the orientedReadIdTable.
    // This is currently not used.
    orientedReadIdTable.clear();
    orientedReadIdTable.resize(orientedReadCount, {invalid<uint64_t>, invalid<uint64_t>});
    for(uint64_t componentId=0; componentId<componentOrientedReadIds.size(); componentId++) {
        const span<OrientedReadId> orientedReadIds = componentOrientedReadIds[componentId];
        for(uint64_t position=0; position<orientedReadIds.size(); position++) {
            const OrientedReadId orientedReadId = orientedReadIds[position];
            orientedReadIdTable[orientedReadId.getValue()] = {componentId, position};
        }
    }

    performanceLog << timestamp << "Mode3Assembler::computeConnectedComponents ends." << endl;
}



void Mode3Assembler::assembleConnectedComponents(
    uint64_t threadCount,
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
    for(uint64_t componentId=0; componentId<componentOrientedReadIds.size(); componentId++) {
        const shared_ptr<AssemblyGraph> assemblyGraph =
            assembleConnectedComponent(componentId, threadCount, true, orientedReadsCsv, debug);
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
        summaryCsv << componentOrientedReadIds[componentId].size() << ",";
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




bool Mode3Assembler::isSelfComplementaryComponent(uint64_t componentId) const
{
    const span<const OrientedReadId> orientedReadIds = componentOrientedReadIds[componentId];

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



void Mode3Assembler::checkComponentIsValid(uint64_t componentId) const
{
    const span<const OrientedReadId> orientedReadIds = componentOrientedReadIds[componentId];

    if(isSelfComplementaryComponent(componentId)) {

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
    bool assembleSequence,
    ostream& orientedReadsCsv,
    bool debug)
{
    // EXPOSE WHEN CODE STABILIZES.
    // Minimum anchor graph edge coverage, only for alignment-free assembly.
    const uint64_t minEdgeCoverageAlignmentFree = 3;


    performanceLog << timestamp << "Assembling connected component " <<
        componentId << " of " << componentOrientedReadIds.size() << endl;
    cout << timestamp << "Assembling connected component " <<
        componentId << " of " << componentOrientedReadIds.size() << endl;

    const span<const OrientedReadId> orientedReadIds = componentOrientedReadIds[componentId];
    const span<const uint64_t> anchorIds = componentAnchorIds[componentId];

    const bool isSelfComplementary = isSelfComplementaryComponent(componentId);
    checkComponentIsValid(componentId);

    cout << "This connected component has " << orientedReadIds.size() <<
        " oriented reads and " << anchorIds.size() << " anchors." << endl;
    if(isSelfComplementary) {
        cout << "This connected component is self-complementary and will be assembled double-stranded." << endl;
    }



    // Store AnchorInfos for this component.
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        const AnchorId anchorId = anchorIds[localAnchorId];
        anchors().storeAnchorInfo(anchorId, uint32_t(componentId), localAnchorId);
    }



    // Write to orientedReadsCsv the oriented reads for this component.
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        const auto readName = reads.getReadName(orientedReadId.getReadId());
        orientedReadsCsv <<
            componentId << "," <<
            orientedReadId << ",";
        copy(readName.begin(), readName.end(), ostream_iterator<char>(orientedReadsCsv));
        orientedReadsCsv << "\n";
    }



    // Now we can create the AnchorGraph for this connected component.
    // The constructor generates the vertices and edges.
    const uint64_t minEdgeCoverage = (
        options.anchorCreationMethod == "FromMarkerKmers" ? minEdgeCoverageAlignmentFree :
        0);
    AnchorGraph anchorGraph(anchors(), anchorIds, minEdgeCoverage);


     cout << "The AnchorGraph for this connected component has " <<
         num_vertices(anchorGraph) << " vertices and " << num_edges(anchorGraph) << " edges." << endl;

     // Graphviz output.
     if(debug) {
         AnchorGraphDisplayOptions options;
         options.showNonTransitiveReductionEdges = true;
         anchorGraph.writeGraphviz(
             "AnchorGraphInitial" + to_string(componentId), options, anchors());
         options.makeCompact();
         anchorGraph.writeGraphviz(
             "AnchorGraphCompactInitial" + to_string(componentId), options, anchors());
         anchorGraph.writeEdgeCoverageHistogram("AnchorGraphInitial" + to_string(componentId) + "-EdgeCoverageHistogram.csv");
         anchorGraph.writeEdgeDetails("AnchorGraphInitial" + to_string(componentId) + "-EdgeDetails.csv", anchors());
     }

     if(options.anchorCreationMethod != "FromMarkerKmers") {

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

         // Graphviz output.
         if(debug) {
             AnchorGraphDisplayOptions options;
             options.showNonTransitiveReductionEdges = false;
             anchorGraph.writeGraphviz(
                 "AnchorGraph" + to_string(componentId), options, anchors());
             options.makeCompact();
             anchorGraph.writeGraphviz(
                 "AnchorGraphCompact" + to_string(componentId), options, anchors());
         }
     }

     // Create the assembly graph for this connected component.
     return make_shared<AssemblyGraph>(
         anchorGraph, anchors(), componentId, k, orientedReadIds, anchorIds, threadCount,
         options, assembleSequence, debug);
}



// Debug output of all connected components.
void Mode3Assembler::writeConnectedComponents() const
{
    for(uint64_t componentId=0; componentId<componentOrientedReadIds.size(); componentId++) {
        writeConnectedComponent(componentId);
    }
}



// Debug output of a connected components.
void Mode3Assembler::writeConnectedComponent(uint64_t componentId) const
{
    const span<const OrientedReadId> orientedReadIds = componentOrientedReadIds[componentId];
    const span<const uint64_t> anchorIds = componentAnchorIds[componentId];

    // Write the oriented reads.
    {
        ofstream csv("OrientedReadIds-" + to_string(componentId) + ".csv");
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            csv << orientedReadId << "\n";
        }
    }



    // Write the anchors.
    {
        ofstream csv("Anchors-" + to_string(componentId) + ".csv");
        csv << "AnchorId,OrientedReadId,Ordinal0,Ordinal1\n";

        for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
            const AnchorId anchorId = anchorIds[localAnchorId];
            const Anchor anchor = anchors()[anchorId];
            for(const AnchorMarkerInterval& markerInterval: anchor) {
                csv << anchorId << ",";
                csv << markerInterval.orientedReadId << ",";
                csv << markerInterval.ordinal0 << ",";
                csv << markerInterval.ordinal0 + anchors().ordinalOffset(anchorId) << "\n";
            }
        }

    }
}


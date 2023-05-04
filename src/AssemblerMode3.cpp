#include "Assembler.hpp"
#include "mode3.hpp"
#include "mode3-Detangler.hpp"
#include "mode3-PathGraph.hpp"
#include "mode3a.hpp"
#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "mode3b-PathFinder.hpp"
#include "Reads.hpp"
using namespace shasta;



void Assembler::mode3Assembly(
    size_t threadCount)
{
    // EXPOSE WHEN CODE STABILIZES.
    // const uint64_t minClusterSize = 3;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    assemblyGraph3Pointer = std::make_shared<mode3::AssemblyGraph>(
        largeDataFileNamePrefix,
        largeDataPageSize,
        threadCount,
        assemblerInfo->readRepresentation,
        assemblerInfo->k,
        *reads,
        markers,
        markerGraph,
        *consensusCaller);
    auto& assemblyGraph3 = *assemblyGraph3Pointer;
    assemblyGraph3.writeGfa("AssemblyGraph");
    // assemblyGraph3.clusterSegments(threadCount, minClusterSize);
    assemblyGraph3.createJaccardGraph(threadCount);
    // assemblyGraph3.assembleJaccardGraphPaths();
    assemblyGraph3.createDeBruijnGraph();

}



void Assembler::accessMode3AssemblyGraph()
{
    assemblyGraph3Pointer = std::make_shared<mode3::AssemblyGraph>(
        largeDataFileNamePrefix,
        assemblerInfo->readRepresentation,
        assemblerInfo->k,
        *reads, markers, markerGraph, *consensusCaller);
}



void Assembler::analyzeMode3Subgraph(const vector<uint64_t>& segmentIds)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    vector<mode3::AssemblyGraph::AnalyzeSubgraphClasses::Cluster> clusters;
    assemblyGraph3Pointer->analyzeSubgraph(segmentIds, clusters, true);
}



void Assembler::createMode3PathGraph()
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    const mode3::AssemblyGraph& assemblyGraph = *assemblyGraph3Pointer;

    mode3::PathGraph pathGraph(assemblyGraph);

}



void Assembler::createMode3Detangler()
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    const mode3::AssemblyGraph& assemblyGraph = *assemblyGraph3Pointer;

    mode3::Detangler detangler(assemblyGraph);

}



void Assembler::mode3aAssembly(
    size_t threadCount)
{

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    mode3a::Assembler mode3aAssembler(
        threadCount,
        assemblerInfo->k,
        MappedMemoryOwner(*this),
        *reads,
        markers,
        markerGraph);
}



void Assembler::accessMode3aAssemblyData()
{
    mode3aAssemblyData.packedMarkerGraph = make_shared<mode3a::PackedMarkerGraph>(
        MappedMemoryOwner(*this), "Mode3a-PackedMarkerGraph", assemblerInfo->k, *reads, markers, markerGraph, true);

    mode3aAssemblyData.assemblyGraphSnapshots.clear();
    for(uint64_t i=0; ; i++) {
        try {
            const string name = "Mode3a-AssemblyGraphSnapshot-" + to_string(i);
            const auto snapshot =
                make_shared<mode3a::AssemblyGraphSnapshot>(
                    name,
                    MappedMemoryOwner(*this),
                    *mode3aAssemblyData.packedMarkerGraph);
            mode3aAssemblyData.assemblyGraphSnapshots.push_back(snapshot);
        } catch (exception&) {
            break;
        }
    }

    cout << "Found " << mode3aAssemblyData.assemblyGraphSnapshots.size() <<
        " assembly graph snapshots." << endl;
}

#pragma once

// Shasta.
#include "mode3-Anchor.hpp"
#include "AssemblerOptions.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Mode3Assembler;

    class MarkerInterval;
    class Mode3AssemblyOptions;

    namespace mode3 {
        class AssemblyGraph;
        class AssemblyGraphPostprocessor;
    }

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}


class shasta::Mode3Assembler :
    public MultithreadedObject<Mode3Assembler>,
    public MappedMemoryOwner {
public:

    // This constructor runs the assembly.
    Mode3Assembler(
        const MappedMemoryOwner&,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        shared_ptr<mode3::Anchors>,
        uint64_t threadCount,
        const Mode3AssemblyOptions&,
        bool debug);

    // This constructor just accesses binary data.
    Mode3Assembler(
        const MappedMemoryOwner&,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        shared_ptr<mode3::Anchors>,
        const Mode3AssemblyOptions&);

private:
    uint64_t k;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    bool debug;

    // The main input to Mode3Assembler is a set of anchors.
    // Each anchor consists of a span of MarkerIntervals, with the following requirements:
    // - All MarkerIntervals correspond to exactly the same sequence in the corresponding oriented reads, and:
    //      * Those portions of the oriented reads are believed to be aligned.
    //      * They apear in a low number of copies in the genome being sequenced.
    // - There are no duplicate oriented reads in an anchor.
    // - The anchor coverage (number of oriented reads) is isn [minPrimaryCoverage, maxPrimaryCoverage].
    // For now the anchors are simply a reference to assembler.markerGraph.edgeMarkerIntervals,
    // butit might be possible to construct the anchors by other means.
    shared_ptr<mode3::Anchors> anchorsPointer;
public:
    mode3::Anchors& anchors()
    {
        return *anchorsPointer;
    }
    const mode3::Anchors& anchors() const
    {
        return *anchorsPointer;
    }

    const Mode3AssemblyOptions& options;



    // The oriented reads present in each anchor
    // define a bipartite graph. We want to compute connected components
    // of this bipartite graph and process them one at a time.
    // These are also connected components of the global anchor graph
    // (with one vertex for each anchor, and edges created by following the reads).

    // The OrientedReadIds in each component. Indexed by componentId.
    MemoryMapped::VectorOfVectors<OrientedReadId, uint64_t> componentOrientedReadIds;

    // The anchorIds in each component. Indexed by componentId.
    MemoryMapped::VectorOfVectors<mode3::AnchorId, uint64_t> componentAnchorIds;

public:

    bool isSelfComplementaryComponent(uint64_t componentId) const;
    void checkComponentIsValid(uint64_t componentId) const;

    void computeConnectedComponents();

    // Debug output of connected components.
    void writeConnectedComponents() const;
    void writeConnectedComponent(uint64_t componentId) const;



    // For each oriented read, store which ConnectedComponent it belongs to,
    // and at what position.
    // Indexed by OrientedReadId::getValue().
    // For each OrientedReadId we store a pair (componentId, position),
    // where componentId is the index in the connectedComponents vector
    // and position is the index in the orientedReadIds vector
    // for that connected component.
    vector< pair<uint64_t, uint64_t> > orientedReadIdTable;

    void assembleConnectedComponents(
        uint64_t threadCount,
        bool debug);
    shared_ptr<mode3::AssemblyGraph> assembleConnectedComponent(
        uint64_t componentId,
        uint64_t threadCount,
        bool assembleSequence,
        ostream& orientedReadsCsv,
        bool debug);

public:
    // Http server functions for Mode 3 assembly.
    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreReadFollowing(const vector<string>&, ostream&);
    void exploreLocalAssembly(
        const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(
        const vector<string>&, ostream&);
    void exploreAssemblyGraph(
        const vector<string>&, ostream&);
    void exploreSegment(
        const vector<string>&, ostream&);
    void exploreReadFollowingAssemblyGraph(const vector<string>&, ostream&);

    // During postprocessing, we maintain a table of currently loaded versions
    // of the AssemblyGraph, keyed by assembly stage and componentId.
    std::map< pair<string, uint64_t>, shared_ptr<const mode3::AssemblyGraphPostprocessor> >
        assemblyGraphsMap;
    const mode3::AssemblyGraphPostprocessor& getAssemblyGraph(
        const string& assemblyStage,
        uint64_t componentId
        );
};

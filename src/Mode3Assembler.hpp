#pragma once

// Shasta.
#include "mode3-Anchor.hpp"
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

    class Assembler;
    class MarkerInterval;
    class Mode3AssemblyOptions;

    namespace mode3 {
        class AssemblyGraph;
    }

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}


class shasta::Mode3Assembler :
    public MultithreadedObject<Mode3Assembler>,
    public MappedMemoryOwner {
public:
    Mode3Assembler(
        const Assembler&,
        uint64_t threadCount,
        const Mode3AssemblyOptions&,
        bool debug);
private:
    const Assembler& assembler;
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
    const mode3::Anchors& anchors;

    // Keep track of the reverse complement of each anchor:
    // the reverse complement of anchorId is reverseComplementAnchor[anchorId],
    // and therefore for any anchorId:
    // reverseComplementAnchor[reverseComplementAnchor[anchorId]] == anchorId.
    // Anchors are not allowed to be self-complementary, and therefore for any anchorId:
    // reverseComplementAnchor[anchorId] != anchorId.
    vector<mode3::AnchorId> reverseComplementAnchor;
    void findReverseComplementAnchors();

    // The oriented reads present in each anchor
    // define a bipartite graph. We want to compute connected components
    // of this bipartite graph and process them one at a time.
    // These are also connected components of the global primary graph
    // (with one vertex for each anchor, and edges created by following the reads).
    class ConnectedComponent {
    public:
        // The oriented reads in this connected component.
        vector<OrientedReadId> orientedReadIds;

        // The anchors (marker graph edge ids) in this connected component.
        vector<mode3::AnchorId> anchorIds;

        bool isSelfComplementary() const;
        void checkIsValid() const;  // SHASTA_ASSERT if not valid.
    };
    vector<ConnectedComponent> connectedComponents;
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
        const Mode3AssemblyOptions&,
        bool debug);
    shared_ptr<mode3::AssemblyGraph> assembleConnectedComponent(
        uint64_t componentId,
        uint64_t threadCount,
        const Mode3AssemblyOptions&,
        bool assembleSequence,
        ostream& orientedReadsCsv,
        bool debug);
};

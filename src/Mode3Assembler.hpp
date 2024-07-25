#pragma once

// Shasta.
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Mode3Assembler;
    class Assembler;
    class Mode3AssemblyOptions;
    namespace mode3 {
        class AssemblyGraph;
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

    // For Mode 3 assembly we only generate primary marker graph edges,
    // so all marker graph edges participate in Mode 3 assembly.
    // Each marker graph edge becomes an "anchor" for Mode 3 assembly.

    // The oriented reads present in each primary marker graph edge
    // define a bipartite graph. We want to compute connected components
    // of this bipartite graph and process them one at a time.
    // These are also connected components of the global primary graph
    // (with one vertex for each primary marker graph edge,
    // and edges created by following the reads).
    class ConnectedComponent {
    public:
        // The oriented reads in this connected component.
        vector<OrientedReadId> orientedReadIds;

        // The PrimaryIds of the marker graph edges in this connected component.
        // These are indices into primaryMarkerGraphEdgeIds.
        vector<uint64_t> primaryIds;
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
        bool debug);
};

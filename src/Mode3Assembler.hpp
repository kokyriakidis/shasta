#pragma once

// Shasta.
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Mode3Assembler;
    class Assembler;
}


class shasta::Mode3Assembler :
    public MultithreadedObject<Mode3Assembler>,
    public MappedMemoryOwner {
public:
    Mode3Assembler(const Assembler&, uint64_t threadCount, bool debug);
private:
    const Assembler& assembler;
    bool debug;

    // The MarkerGraphEdgeIds of the primary marker graph edges.
    // These are sorted.
    // An index in this vector is called PrimaryId.
    vector<MarkerGraphEdgeId> primaryMarkerGraphEdgeIds;
    void gatherPrimaryMarkerGraphEdgeIds();

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

    // For each oriented read, store which ConnectedComponent it belongs to,
    // and at what position.
    // Indexed by OrientedReadId::getValue().
    // For each OrientedReadId we store a pair (componentId, position),
    // where componentId is the index in the connectedComponents vector
    // and position is the index in the orientedReadIds vector
    // for that connected component.
    vector< pair<uint64_t, uint64_t> > orientedReadIdTable;

    void assembleConnectedComponents(uint64_t threadCount);
    void assembleConnectedComponent(uint64_t componentId, uint64_t threadCount);
};
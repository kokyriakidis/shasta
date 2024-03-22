#pragma once

// Shasta.
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "vector.hpp"

namespace shasta {
    class Mode3Assembler;
    class Assembler;
}


class shasta::Mode3Assembler :
    public MultithreadedObject<Mode3Assembler>,
    public MappedMemoryOwner {
public:
    Mode3Assembler(const Assembler&, bool debug);
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
    void computeConnectedComponents();
};

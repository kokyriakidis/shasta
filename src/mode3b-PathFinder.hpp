#ifndef SHASTA_MODE3B_PATH_FINDER_HPP
#define SHASTA_MODE3B_PATH_FINDER_HPP

#include "shastaTypes.hpp"

namespace shasta {
    namespace mode3b {
        class PathFinder;
    }

    class CompressedMarker;
    class MarkerGraph;
    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
};


// Find and assemble a path in the complete marker graph.
class shasta::mode3b::PathFinder {
public:

    PathFinder(
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        MarkerGraphEdgeId startEdgeId,  // The path starts here.
        uint64_t direction              // 0=forward, 1=backward
        );
private:

    // Thinsg we get from the constructor.
    uint64_t k;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;


    // Move in the specified direction until we find an edge with similar
    // read composition which is suitable to becoe the next primary vertex.
    MarkerGraphEdgeId findNextPrimaryEdge(MarkerGraphEdgeId, uint64_t direction);

    // Given two marker graph edges, figure out if they have
    // similar read compositions.
    bool haveSimilarCompositions(MarkerGraphEdgeId, MarkerGraphEdgeId) const;
};

#endif

#ifndef SHASTA_MODE3B_PATH_FINDER_HPP
#define SHASTA_MODE3B_PATH_FINDER_HPP

#include "shastaTypes.hpp"

namespace shasta {
    namespace mode3b {
        class PathFinder;
    }

    class Assembler;
};


// Find and assemble a path in the complete marker graph.
class shasta::mode3b::PathFinder {
public:

    PathFinder(
        const Assembler&,
        MarkerGraphEdgeId startEdgeId,  // The path starts here.
        uint64_t direction              // 0=forward, 1=backward
        );
private:

    // Things we get from the constructor.
    const Assembler& assembler;


    // Move in the specified direction until we find an edge with similar
    // read composition which is suitable to becoe the next primary vertex.
    MarkerGraphEdgeId findNextPrimaryEdge(
        MarkerGraphEdgeId,
        uint64_t direction,
        uint64_t maxMarkerOffset,
        uint64_t minCommonCount,
        double minCorrectedJaccard) const;

};

#endif

#ifndef SHASTA_MODE3B_ASSEMBLY_PATH_HPP
#define SHASTA_MODE3B_ASSEMBLY_PATH_HPP

#include "MarkerGraphEdgePairInfo.hpp"
#include "shastaTypes.hpp"

#include "iosfwd.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3b {
        class AssemblyPath;
    }

    class Assembler;
    class Base;
}



// An AssemblyPath is path in the marker graph that is used to assemble sequence.
// It is created starting with of N+1 primary vertices, numbered from 0 to N,
// found by mode3b::PathFinder. The primary vertices are not necessarily
// adjacent in the marker graph.
// Two consecutive primary vertices have similar read compositions
// and therefore a number of common oriented reads.
// Each pair of consecutive primary vertices generates an AssemblyStep.
// There are N AssemblySteps, numbered from 0 to N-1.
// Each AssemblyStep consists of a number of secondary vertices
// created by mode3b::PathFiller.
// The primary and secondary vertices are a path in the marker graph and
// can be used to assemble sequence.
class shasta::mode3b::AssemblyPath {
public:
    AssemblyPath(
        const Assembler&,
        MarkerGraphEdgeId,
        uint64_t direction  // 0 = forward, 1 = backward.
        );

    void getSequence(vector<Base>&) const;
    void writeFasta(ostream&) const;
private:
    const Assembler& assembler;

    // There are N+1 primary edges numbered from 0 to N
    // and N steps numbered from 0 to N-1.
    // The i-th step goes from primary edge i
    // to primary edge i+1.
    vector<MarkerGraphEdgeId> primaryEdges;
    class Step {
    public:
        MarkerGraphEdgePairInfo info;
        vector<MarkerGraphEdgeId> secondaryEdges;
        Step(const MarkerGraphEdgePairInfo& info);
    };
    vector<Step> steps;
};

#endif
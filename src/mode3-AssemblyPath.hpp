#pragma once

#include "MarkerGraphEdgePairInfo.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"

#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class AssemblyPath;
    }

    class Assembler;
    class Base;
}



// An AssemblyPath is path in the marker graph that is used to assemble sequence.
// It is created starting with of N+1 primary vertices, numbered from 0 to N,
// found by mode3::PathFinder. The primary vertices are not necessarily
// adjacent in the marker graph.
// Two consecutive primary vertices have similar read compositions
// and therefore a number of common oriented reads.
// Each pair of consecutive primary vertices generates an AssemblyStep.
// There are N AssemblySteps, numbered from 0 to N-1.
// Each AssemblyStep consists of a number of secondary vertices
// created by mode3::PathFiller.
// The primary and secondary vertices are a path in the marker graph and
// can be used to assemble sequence.
class shasta::mode3::AssemblyPath : public MultithreadedObject<AssemblyPath> {
public:

    // Create the assembly path starting from a given primary edge.
    AssemblyPath(
        const Assembler&,
        MarkerGraphEdgeId,
        uint64_t direction,  // 0 = forward, 1 = backward, 2=bidirectional
        bool allowOrientedReadsOnFirst = true,
        bool allowOrientedReadsOnLast = true
        );

    // Create the assembly path given n primary edges and
    // the n-1 MarkerGraphEdgePairInfo between consecutive primary edges.
    AssemblyPath(
        const Assembler&,
        const vector<MarkerGraphEdgeId>&,
        const vector<MarkerGraphEdgePairInfo>&,
        bool allowOrientedReadsOnFirst = true,  // Allow using for assembly the oriented reads on the first edge of the path.
        bool allowOrientedReadsOnLast  = true,  // Allow using for assembly the oriented reads on the last edge of the path.
        uint64_t threadCount = 0);

    void getSequence(vector<Base>&) const;
    void getInternalSequence(vector<Base>&) const;  // Excluding the begin/end edges.
    void writeFasta(ostream&, const string& name) const;
    void writeCsv(ostream&, const string& name) const;
private:
    const Assembler& assembler;

    // If set, allow using for assembly the oriented reads on the first edge of the path.
    bool allowOrientedReadsOnFirst;

    // If set, allow using for assembly the oriented reads on the first edge of the path.
    bool allowOrientedReadsOnLast;

    // There are N+1 primary edges numbered from 0 to N
    // and N steps numbered from 0 to N-1.
    // The i-th step goes from primary edge i
    // to primary edge i+1.
    vector<MarkerGraphEdgeId> primaryEdges;

    // Each of the N steps stores the assembled
    // sequence intervening between the trimary edges for the step.
    class Step {
    public:
        MarkerGraphEdgePairInfo info;
        vector<Base> sequence;
        Step(const MarkerGraphEdgePairInfo& info);
    };
    vector<Step> steps;

    // Create the primaryEdges and the steps.
    void create(
        MarkerGraphEdgeId,
        uint64_t direction  // 0 = forward, 1 = backward, 2=bidirectional
        );

    // Assemble the sequence of each Step.
    void assembleSequential();
    void assembleStep(uint64_t i);
    void assembleParallel(uint64_t threadCount);
    void assembleThreadFunction(uint64_t threadId);
};


#ifndef SHASTA_MODE3B_PATH_GRAPH_HPP
#define SHASTA_MODE3B_PATH_GRAPH_HPP

/*******************************************************************************

In the mode3b::PathGraph, each vertex corresponds to a primary edge of
of the marker graph, which is believed to correspond to a single copy
of sequence. It is characterized as follows:
- minCoverage <= coverage <= maxCoverage
- No duplicate oriented reads on the marker graph edge or its vertices.

A directed edge v0->v1 is generated if a sufficient number of oriented reads
visit the marker graph edge corresponding to v1 after
the marker graph edge corresponding to v0, without visiting other
primary marker graph edges in between.

*******************************************************************************/


namespace shasta {
    class Assembler;
    namespace mode3b {
        class PathGraph;
    }
}


class shasta::mode3b::PathGraph {
public:
    PathGraph(const Assembler&);
private:
    const Assembler& assembler;
};

#endif


#ifndef SHASTA_MODE3B_PATH_FILLER_HPP
#define SHASTA_MODE3B_PATH_FILLER_HPP

// Shasta.
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include <map>
#include "vector.hpp"


namespace shasta {
    namespace mode3b {
        class PathFiller;
    }

    class Assembler;
    class Base;
    class MarkerInterval;
};


// Assemble a step between two primary edges of an assembly path.
// Given the two marker graph that are primary for an assembly path,
// find the secondary edges in between and assemble the sequence.
// This works with mode3b assembly and the complete marker graph.
// This assumes that the two edges have no duplicate oriented read ids.
class shasta::mode3b::PathFiller {
public:

    PathFiller(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        ostream& html);

private:
    const Assembler& assembler;
    MarkerGraphEdgeId edgeIdA;
    MarkerGraphEdgeId edgeIdB;

    void checkAssumptions() const;

    // The OrientedReadIds we will be using for this assembly.
    // These are the ones that are common between edgeIdA and edgeIdB,
    // and that visit edgeIdB after visiting edgeIdA.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The ordinals of this oriented read on the two edges.
        uint32_t ordinalA0;
        uint32_t ordinalA1;
        uint32_t ordinalB0;
        uint32_t ordinalB1;
    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();


    // A local marker graph that sees only the oriented reads used in this assembly.
    class Vertex {
    public:
        MarkerGraphVertexId vertexId;
    };
    class Edge {
    public:
        MarkerGraphEdgeId edgeId;
        vector<MarkerInterval> markerIntervals;
        vector<Base> sequence;
    };
    using GraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex,
        Edge
        >;
    class Graph: public GraphBaseClass {
    public:
        Graph(const Assembler&);

        std::map<MarkerGraphVertexId, vertex_descriptor> vertexMap;
        vertex_descriptor addVertexIfNecessary(MarkerGraphVertexId);

        std::map<MarkerGraphEdgeId, edge_descriptor> edgeMap;
        edge_descriptor addEdgeIfNecessary(MarkerGraphEdgeId);

        void writeGraphviz(ostream&, uint64_t peakCoverage) const;

    private:
        const Assembler& assembler;
    };
    Graph graph;
    void createGraph();
    void writeGraph(ostream& html) const;

    // EXPOSE WHEN CODE STABILIZES.
    static const uint64_t minEdgeCoverage = 6;

};

#endif

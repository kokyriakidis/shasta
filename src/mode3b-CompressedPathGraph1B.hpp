#ifndef SHASTA_MODE3B_COMPRESSED_PATH_GRAPH1B_HPP
#define SHASTA_MODE3B_COMPRESSED_PATH_GRAPH1B_HPP

// Shasta
#include "invalid.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library
#include <map>
#include "string.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3b {

        // Each edge of the CompressedPathGraph1B describes a BubbleChain.

        // A Chain is a sequence of MarkerGraphEdgeIds.
        // It can be used to generate an AssemblyPath.
        using Chain = vector<MarkerGraphEdgeId>;

        // A Bubble is a set of Chains that begin and end at the same MarkerGraphEdgeId.
        // It can consist of one or more Chains.
        class Bubble;

        // A BubbleChain is a sequence of Bubbles.
        using BubbleChain = vector<Bubble>;

        class CompressedPathGraph1B;
        class CompressedPathGraph1BVertex;
        class CompressedPathGraph1BEdge;
        using CompressedPathGraph1BBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            CompressedPathGraph1BVertex,
            CompressedPathGraph1BEdge>;

        class PathGraph1;
    }
    class Assembler;
}



class shasta::mode3b::Bubble : public vector<Chain> {
public:
    bool isHaploid() const
    {
        return size() == 1;
    }
    bool isDiploid() const
    {
        return size() == 2;
    }
    bool isGeneral() const
    {
        return size() > 2;
    }
};



class shasta::mode3b::CompressedPathGraph1BVertex {
public:
    MarkerGraphEdgeId edgeId;
};



class shasta::mode3b::CompressedPathGraph1BEdge : public BubbleChain {
public:
    uint64_t id = invalid<uint64_t>;
};



class shasta::mode3b::CompressedPathGraph1B: public CompressedPathGraph1BBaseClass {
public:
    CompressedPathGraph1B(
        const PathGraph1&,
        uint64_t componentId,
        const Assembler&);
private:
    // Information stored by the constructor.
    const PathGraph1& graph;
    uint64_t componentId;
    const Assembler& assembler;

    // Initial creation from the PathGraph1.
    // Each linear chain of edges in the PathGraph1 after transitive reduction generates
    // a CompressedPathGraph1BEdge (BubbleChain) consisting of a single haploid bubble.
    void create();
    std::map<MarkerGraphEdgeId, vertex_descriptor> vertexMap;
    uint64_t nextEdgeId = 0;

    // Return the vertex corresponding to a given MarkerGraphEdgeId,
    // creating it if necessary.
    vertex_descriptor getVertex(MarkerGraphEdgeId);

    // Compress parallel edges into bubbles, where possible.
    void compressParallelEdges();

    // Output.
    void write(const string& name) const;
    void writeCsv(const string& fileNamePrefix) const;
    void writeBubbleChainsCsv(const string& fileNamePrefix) const;
    void writeBubblesCsv(const string& fileNamePrefix) const;
    void writeChainsCsv(const string& fileNamePrefix) const;
    void writeChainsDetailsCsv(const string& fileNamePrefix) const;
    void writeGraphviz(const string& fileNamePrefix) const;
    void writeGfa(const string& fileNamePrefix) const;

    string bubbleChainStringId(edge_descriptor) const;
    string bubbleStringId(edge_descriptor, uint64_t positionInBubbleChain) const;
    string chainStringId(edge_descriptor, uint64_t positionInBubbleChain, uint64_t indexInBubble) const;

    uint64_t chainOffset(const Chain&) const;
    void bubbleOffset(
        const Bubble&,
        uint64_t& averageOffset,
        uint64_t& minOffset,
        uint64_t& maxOffset
        ) const;
    void bubbleChainOffset(
        const BubbleChain&,
        uint64_t& averageOffset,
        uint64_t& minOffset,
        uint64_t& maxOffset
        ) const;
};

#endif

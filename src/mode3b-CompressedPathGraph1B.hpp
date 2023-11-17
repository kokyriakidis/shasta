#ifndef SHASTA_MODE3B_COMPRESSED_PATH_GRAPH1B_HPP
#define SHASTA_MODE3B_COMPRESSED_PATH_GRAPH1B_HPP

// Shasta
#include "Base.hpp"
#include "invalid.hpp"
#include "shastaTypes.hpp"
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library
#include "array.hpp"
#include <map>
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3b {

        // Each edge of the CompressedPathGraph1B describes a BubbleChain.

        // A Chain is a sequence of MarkerGraphEdgeIds.
        // It can be used to generate an AssemblyPath.
        class Chain;

        // A Bubble is a set of Chains that begin and end at the same MarkerGraphEdgeId.
        // It can consist of one or more Chains.
        class Bubble;

        // A BubbleChain is a sequence of Bubbles.
        class BubbleChain;

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



// A Chain is a sequence of MarkerGraphEdgeIds.
// It can be used to generate an AssemblyPath.
class shasta::mode3b::Chain : public vector<MarkerGraphEdgeId> {
public:
    vector<Base> sequence;
};



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



class shasta::mode3b::BubbleChain : public vector<Bubble> {
public:
    const Bubble& firstBubble() const
    {
        SHASTA_ASSERT(not empty());
        return front();
    }
    const Bubble& lastBubble() const
    {
        SHASTA_ASSERT(not empty());
        return back();
    }

    // Collapse consecutive haploid bubbles.
    void compress();
};



class shasta::mode3b::CompressedPathGraph1BVertex {
public:
    MarkerGraphEdgeId edgeId;

    // Numbering of vertices consecutively starting at zero.
    // This is computed by renumberVertices, and becomes
    // invalid as soon as a vertex is added or removed.
    uint64_t index = invalid<uint64_t>;

    // The id of the Superbubble this vertex belongs to, if any.
    // Stored by class Superbubbles.
    uint64_t superbubbleId = invalid<uint64_t>;

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
        const Assembler&,
        uint64_t threadCount0,
        uint64_t threadCount1);
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
    void renumberEdges();

    // Return the vertex corresponding to a given MarkerGraphEdgeId,
    // creating it if necessary.
    vertex_descriptor getVertex(MarkerGraphEdgeId);

    void removeVertex(vertex_descriptor);

    // Compute vertexIndex for every vertex.
    // This numbers vertices consecutively starting at zero.
    // This numbering becomes invalid as soon as a vertex is added or removed.
    void numberVertices();
    void clearVertexNumbering();

    // Compress parallel edges into bubbles, where possible.
    bool compressParallelEdges();

    // Compress linear sequences of edges (BubbleChains) into longer BubbleChains.
    bool compressSequentialEdges();

    // Call compressParallelEdges and compressSequentialEdges iteratively until nothing changes.
    void compress();

    // Compute the tangle matrix given in-edges and out-edges.
    // The last bubble of each in-edge and the first bubble
    // of each out-edge must be haploid.
    void computeTangleMatrix(
        const vector<edge_descriptor>& inEdges,
        const vector<edge_descriptor>& outEdges,
        vector< vector<uint64_t> >& tangleMatrix
        ) const;

    // Vertex detangling.
    bool detangleVerticesStrict(bool debug);
    bool detangleVertexStrict(vertex_descriptor, bool debug);
    bool detangleVertices(bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);
    bool detangleVertex(
        vertex_descriptor,
        bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);

    // Vertex detangling that can deal with non-haploid bubbles adjacent to the
    // vertex to be detangled.
    bool detangleVerticesGeneral(bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);
    bool detangleVertexGeneral(
        vertex_descriptor,
        bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);

    // Split the first/last bubble of a bubble chain.
    // Used by detangleVertexGeneral to eliminate
    // non-haploid bubble adjacent to a vertex to be detangled.
    void splitBubbleChainAtBeginning(edge_descriptor);
    void splitBubbleChainAtEnd(edge_descriptor);


    // Edge detangling.
    bool detangleEdges(
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);
    bool detangleEdge(
        std::map<uint64_t, edge_descriptor>& edgeMap,
        std::map<uint64_t, edge_descriptor>::iterator&,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);

    // Special treatment to detangle back edges that were too long
    // to be handled by detangleEdges.
    bool detangleBackEdges(
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);
    bool detangleBackEdge(
        std::map<uint64_t, edge_descriptor>& edgeMap,
        std::map<uint64_t, edge_descriptor>::iterator&,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh);



    // Find short superbubbles in the CompressedPathGraph1B.
    class Superbubbles {
    public:
        Superbubbles(
            CompressedPathGraph1B&,
            uint64_t maxOffset1,    // Used to define superbubbles
            uint64_t maxOffset2     // Compared against the offset between entries and exits
            );
        ~Superbubbles();

        // Return the number of superbubbbles.
        uint64_t size() const
        {
            return superbubbles.size();
        }

        // Return the vertices in the specified superbubble.
        const vector<vertex_descriptor>& getSuperbubble(uint64_t superBubbleId) const
        {
            return components[superbubbles[superBubbleId]];
        }

        // Figure out if a vertex is in the specified superbubble.
        bool isInSuperbubble(uint64_t superbubbleId, vertex_descriptor cv)
        {
            return cGraph[cv].superbubbleId == superbubbleId;
        }

    private:

        CompressedPathGraph1B& cGraph;

        // The connected components of the CompressedPathGraph1B, computed
        // using only edges with offset up to maxOffset1.
        vector< vector<vertex_descriptor> > components;

        // The superbubbles are the components with size at least 2.
        // Store the componentIds (index in components vector).
        vector<uint64_t> superbubbles;
    };



    // Remove short superbubbles with one entry and one exit.
    void removeShortSuperbubbles(
        uint64_t maxOffset1,    // Used to define superbubbles
        uint64_t maxOffset2     // Compared against the offset between entry and exit
    );



    // Phasing of bubble chains.
    void phaseBubbleChains(
        bool debug,
        uint64_t lowThreshold,
        uint64_t highThreshold,
        uint64_t longBubbleThreshold);
    void phaseBubbleChain(
        edge_descriptor e,
        uint64_t lowThreshold,
        uint64_t highThreshold,
        uint64_t longBubbleThreshold,
        bool debug);

    // In the phasing graph, each vertex corresponds to a diploid bubble
    // in the BubbleChain being phased.
    class TangleMatrix : public array< array<uint64_t, 2>, 2> {
    public:
        void analyze(
            uint64_t lowThreshold,
            uint64_t highThreshold,
            int64_t& phase,
            uint64_t& minConcordant,
            uint64_t& maxDiscordant,
            uint64_t& total) const;
    };


    // A PhasedComponent is a set of phased diploid bubbles
    // in a BubbleChain.
    // It is a vector of (bubble position in bubble chain, phase),
    // sorted by bubble position in bubble chain.
    // PhasedComponents are created in such a way that their position ranges
    // in the bubble chain are not overlapping.
    class PhasedComponent : public vector< pair<uint64_t, int64_t> > {
    public:
        uint64_t minPositionInBubbleChain;
        uint64_t maxPositionInBubbleChain;
        void sort();
    };

    class PhasingGraphVertex {
    public:
        uint64_t positionInBubbleChain;
        int64_t phase = 0;  // +1 or -1 for phased vertices, 0 otherwise
    };

    class PhasingGraphEdge {
    public:
        int64_t phase;          // +1 (in phase) or -1 (out of phase)

        // Tangle matrix metrics.
        // If phase = +1, minConcordant = min(m00, m11), maxDiscordant = max(m01, m10).
        // If phase = -1, minConcordant = min(m01, m10), maxDiscordant = max(m00, m11).
        uint64_t minConcordant;
        uint64_t maxDiscordant;

        bool operator<(const PhasingGraphEdge& that) const
        {
            if(maxDiscordant < that.maxDiscordant) {
                return true;
            }
            if(maxDiscordant > that.maxDiscordant) {
                return false;
            }
            return minConcordant > that.minConcordant;
        }
        bool isSpanningTreeEdge = false;
    };
    using PhasingGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::undirectedS,
        PhasingGraphVertex,
        PhasingGraphEdge>;
    class PhasingGraph : public PhasingGraphBaseClass {
    public:
        void phase();
        bool isConsistent(edge_descriptor) const;
        vector< shared_ptr<PhasedComponent> > phasedComponents;
    };

    // Sequence assembly.
    void assembleChains(uint64_t threadCount0, uint64_t threadCount1);
    void assembleChain(Chain&, uint64_t threadCount1) const;


    // Output.
    void write(const string& name) const;
    void writeCsv(const string& fileNamePrefix) const;
    void writeBubbleChainsCsv(const string& fileNamePrefix) const;
    void writeBubblesCsv(const string& fileNamePrefix) const;
    void writeChainsCsv(const string& fileNamePrefix) const;
    void writeChainsDetailsCsv(const string& fileNamePrefix) const;
    void writeGraphviz(const string& fileNamePrefix, bool labels) const;
    void writeGfa(const string& fileNamePrefix) const;
    void writeGfaExpanded(const string& fileNamePrefix, bool includeSequence) const;
    void writeFastaExpanded(const string& fileNamePrefix) const;

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

#pragma once

// Shasta
#include "Base.hpp"
#include "invalid.hpp"
#include "mode3-Anchor.hpp"
#include "mode3-Detangler.hpp"
#include "mode3-PhasedComponent.hpp"
#include "mode3-Superbubbles.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/serialization/vector.hpp>

// Standard library
#include "algorithm.hpp"
#include "array.hpp"
#include <map>
#include "memory.hpp"
#include "fstream.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {

    class Assembler;
    class Mode3AssemblyOptions;
    class OrientedReadId;

    namespace mode3 {

        // Each edge of the AssemblyGraph describes a BubbleChain.

        // A Chain is a sequence of AnchorIds.
        class Chain;

        // A Bubble is a set of Chains that begin and end at the same AnchorId.
        // It can consist of one or more Chains.
        class Bubble;

        // A BubbleChain is a sequence of Bubbles.
        class BubbleChain;

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            AssemblyGraphVertex,
            AssemblyGraphEdge>;
        class AssemblyGraphEdgePredicate;
        class AssemblyGraphCrossEdgePredicate;
        class AssemblyGraphNoInternalAnchorsEdgePredicate;

        class ChainIdentifier;

        class AnchorGraph;
        class Anchors;
        class ChainPermutationDetangler;
    }
}



// A Chain is a sequence of AnchorIds.
class shasta::mode3::Chain : public vector<AnchorId> {
public:

    // Flag used to indicate that this Chain needs to be assembled.
    // Used by assembleChainsMultithreaded.
    bool shouldBeAssembled = false;
    bool wasAssembled = false;

    // Assembled sequence, including the sequence of the first and
    // last primary marker graph edges.
    vector<Base> sequence;

    // The internal sequence assembled between consecutive pairs
    // of AnchorIds in the chain.
    // If a local assembly fails, the success flag remains false and the sequence remains empty.
    class StepSequence {
    public:
        vector<Base> sequence;
        bool success = false;
    };
    vector<StepSequence> stepSequences;


    AnchorId second() const
    {
        SHASTA_ASSERT(size() > 1);
        return (*this)[1];
    }
    AnchorId secondToLast() const
    {
        SHASTA_ASSERT(size() > 1);
        return (*this)[size() - 2];
    }

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AnchorId> >(*this);
        ar & shouldBeAssembled;
        ar & wasAssembled;
        ar & sequence;
    }
};



class shasta::mode3::Bubble : public vector<Chain> {
public:
    bool isHaploid() const
    {
        return size() == 1;
    }
    bool isDiploid() const
    {
        return size() == 2;
    }

    // Remove duplicate chains.
    void deduplicate();

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<Chain> >(*this);
    }
};



class shasta::mode3::BubbleChain : public vector<Bubble> {
public:
    const Bubble& firstBubble() const
    {
        SHASTA_ASSERT(not empty());
        return front();
    }
    Bubble& firstBubble()
    {
        SHASTA_ASSERT(not empty());
        return front();
    }
    const Bubble& lastBubble() const
    {
        SHASTA_ASSERT(not empty());
        return back();
    }
    Bubble& lastBubble()
    {
        SHASTA_ASSERT(not empty());
        return back();
    }

    uint64_t diploidBubbleCount() const
    {
        uint64_t n = 0;
        for(const Bubble& bubble: *this) {
            if(bubble.isDiploid()) {
                ++n;
            }
        }
        return n;
    }

    // This returns true if this BubbleChain consists of a single haploid bubble.
    bool isSimpleChain() const
    {
        return size() == 1 and firstBubble().isHaploid();
    }
    Chain& getOnlyChain()
    {
        SHASTA_ASSERT(isSimpleChain());
        return firstBubble().front();
    }
    const Chain& getOnlyChain() const
    {
        SHASTA_ASSERT(isSimpleChain());
        return firstBubble().front();
    }

    // Collapse consecutive haploid bubbles.
    bool compress();


    // Return the AnchorId of the initial Anchor of this BubbleChain.
    AnchorId firstAnchorId() const
    {
        // Sanity check.
        SHASTA_ASSERT(not empty());

        // Access the first Bubble of this BubbleChain.
        const Bubble& firstBubble = front();

        // Get the first AnchorId of the first Chain of this bubble.
        const AnchorId anchorId = firstBubble.front().front();

        // Check that the AnchorId is the same for all the Chains of this bubble.
        for(const Chain& chain: firstBubble) {
            SHASTA_ASSERT(chain.front() == anchorId);
        }

        return anchorId;
    }



    // Return the AnchorId of the final Anchor of this BubbleChain.
    AnchorId lastAnchorId() const
    {
        // Sanity check.
        SHASTA_ASSERT(not empty());

        // Access the last Bubble of this BubbleChain.
        const Bubble& lastBubble = back();

        // Get the last AnchorId of the first Chain of this bubble.
        const AnchorId anchorId = lastBubble.front().back();

        // Check that the AnchorId is the same for all the Chains of this bubble.
        for(const Chain& chain: lastBubble) {
            SHASTA_ASSERT(chain.back() == anchorId);
        }

        return anchorId;
    }



    // Return the total length of this BubbleChain.
    uint64_t totalLength() const;

    // Return the total number of anchors in this BubbleChain.
    uint64_t anchorCount() const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<Bubble> >(*this);
    }

};



class shasta::mode3::ChainIdentifier {
public:
    AssemblyGraphBaseClass::edge_descriptor e;
    uint64_t positionInBubbleChain;
    uint64_t indexInBubble;

    ChainIdentifier() {}
    ChainIdentifier(
        AssemblyGraphBaseClass::edge_descriptor e,
        uint64_t positionInBubbleChain,
        uint64_t indexInBubble) :
        e(e),
        positionInBubbleChain(positionInBubbleChain),
        indexInBubble(indexInBubble)
        {}

};



class shasta::mode3::AssemblyGraphVertex {
public:
    AssemblyGraphVertex(AnchorId anchorId = invalid<AnchorId>) : anchorId(anchorId) {}

    AnchorId getAnchorId() const
    {
        return anchorId;
    }

    // Numbering of vertices consecutively starting at zero.
    // This is computed by renumberVertices, and becomes
    // invalid as soon as a vertex is added or removed.
    uint64_t index = invalid<uint64_t>;

    // The id of the Superbubble this vertex belongs to, if any.
    // Stored by storeSuperbubblesInformation.
    uint64_t superbubbleId = invalid<uint64_t>;

    template<class Archive> void serialize(Archive & ar, unsigned int /* version */)
    {
        ar & anchorId;
    }

    AnchorId anchorId;
};



class shasta::mode3::AssemblyGraphEdge : public BubbleChain {
public:
    uint64_t id = invalid<uint64_t>;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<BubbleChain>(*this);
        ar & id;
    }
};



class shasta::mode3::AssemblyGraph:
    public AssemblyGraphBaseClass,
    public MultithreadedObject<shasta::mode3::AssemblyGraph>,
    public MappedMemoryOwner {
public:

    // Create from a connected component of the AnchorGraph, then call run.
    AssemblyGraph(
        const AnchorGraph&,
        const Anchors&,
        uint64_t componentId,
        uint64_t k,
        span<const OrientedReadId>,
        span<const AnchorId>,
        uint64_t threadCount,
        const Mode3AssemblyOptions& options,
        bool assembleSequence,
        bool debug);

    // Constructor from binary data, for postprocessing.
    // Used by AssemblyGraphPostprocessor in the http server.
    AssemblyGraph(
        const string& assemblyStage,
        uint64_t componentId,
        span<const OrientedReadId>,
        span<const AnchorId>,
        const Anchors&,
        const Mode3AssemblyOptions&);

    // Another constructor from binary data, used in the Python API.
    // This is similar to the previous constructor, but gets the
    // orientedReadIds, anchorIds, and Anchors from the Mode3AssemblyGraph in the Assembler.
    AssemblyGraph(
        const string& assemblyStage,
        uint64_t componentId,
        const Assembler&,
        const Mode3AssemblyOptions&);

    // Constructor from a vector of vectors of AnchorIds representing Chains.
    // Used for detangling with read following.
    AssemblyGraph(
        const Anchors&,
        uint64_t componentId,
        uint64_t k,
        span<const OrientedReadId>,
        span<const AnchorId>,
        const vector< vector<AnchorId> >& anchorChains,
        uint64_t threadCount,
        const Mode3AssemblyOptions& options,
        bool debug);


    // Hide Base defined by the base class.
    using Base = shasta::Base;

    // Information stored by the constructor.
    uint64_t componentId;
    const Anchors& anchors;
    uint64_t k;
    const Mode3AssemblyOptions& options;

    // The OrientedReadIds of the connected component that generated this AssemblyGraph.
    // These are sorted.
    span<const OrientedReadId> orientedReadIds;

    // Get the index of an OrientedReadId in the orientedReadIds sorted vector.
    uint64_t getOrientedReadIndex(OrientedReadId) const;

    // The AnchorIds of the anchors
    // of the connected component that generated this AssemblyGraph.
    // These are sorted.
    // An index in this vector is called a local anchor id.
    span<const AnchorId> anchorIds;

    bool sequenceWasAssembled = false;

private:
    // void computeJourneys(bool debug);

    void run(
        uint64_t threadCount,
        bool assembleSequence,
        bool debug);

    // Alternate version for testing.
    void run3(
        uint64_t threadCount,
        bool assembleSequence,
        bool debug);



    // Initial creation from the AnchorGraph.
    // Each linear chain of edges in the AnchorGraph after transitive reduction generates
    // an AssemblyGraphEdge (BubbleChain) consisting of a single haploid bubble.
    void create(const AnchorGraph&, bool debug);
public:
    uint64_t nextEdgeId = 0;
private:
    void renumberEdges();

    // Return the vertex corresponding to a given AnchorId,
    // creating it if it is not in the given vertexMap.
    // This is only used in create().
    vertex_descriptor getVertex(
        AnchorId,
        std::map<AnchorId, vertex_descriptor>& vertexMap
        );

    // Create a new vertex with a given AnchorId.
    vertex_descriptor createVertex(AnchorId);

    void removeVertex(vertex_descriptor);

    // Compute vertexIndex for every vertex.
    // This numbers vertices consecutively starting at zero.
    // This numbering becomes invalid as soon as a vertex is added or removed.
public:
    void numberVertices();
    void clearVertexNumbering();
private:

    // Create a new edge connecting cv0 and cv1.
    // The new edge will consist of a simple BubbleChain with a single
    // haploid Bubble with a Chain of length 2.
    edge_descriptor connect(vertex_descriptor cv0, vertex_descriptor cv1);

    // Create a new edge consisting of the concatenation of two edges,
    // which must be simple chains. The original edges are not removed.
    // The concatenation is done connecting one of the last n
    // AnchorIds of e0 (excluding the very last) and
    // one of the first n AnchorIds of e1 (excluding the very first),
    // choosing the pair with the largest number of common oriented reads.
    edge_descriptor connect(bool debug, edge_descriptor e0, edge_descriptor e1, uint64_t n);

    // Compress parallel edges into bubbles, where possible.
    bool compressParallelEdges();

    // Compress linear sequences of edges (BubbleChains) into longer BubbleChains.
    bool compressSequentialEdges();

    // Call compress on all BubbleChains to merge adjacent haploid bubbles.
    bool compressBubbleChains();

    // Call compressParallelEdges, compressSequentialEdges, and compressBubbleChains
    // iteratively until nothing changes.
    bool compress();

    // This does the opposite of compress. All bubble chains that
    // consist of more than one simple haploid bubble are expanded into one
    // edge for each edge of each bubble.
    // For optimal results it is best to call compressBubbleChains before expand.
    void expand();

    uint64_t totalChainCount() const;

    // Compute the tangle matrix given in-edges and out-edges.
    // The last bubble of each in-edge and the first bubble
    // of each out-edge must be haploid.
    void computeTangleMatrix(
        const vector<edge_descriptor>& inEdges,
        const vector<edge_descriptor>& outEdges,
        vector< vector<uint64_t> >& tangleMatrix
        ) const;

    // Low level primitives used in detangling.
    // See the implementation for details.
    vertex_descriptor cloneAndTruncateAtEnd(bool debug, edge_descriptor);
    vertex_descriptor cloneAndTruncateAtBeginning(bool debug, edge_descriptor);

    // Vertex detangling.
    bool detangleVertices(bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh,
        bool useBayesianModel,
        double epsilon,
        double minLogP);
    bool detangleVertex(
        vertex_descriptor,
        bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh,
        bool useBayesianModel,
        double epsilon,
        double minLogP);

    // This handles the case of a vertex with one in-edge, one out-edge.
    // and one in-out-edge (cycle).
    bool detangleVertexWithCycle(
        vertex_descriptor,
        bool debug,
        double epsilon,
        double minLogP);

    // Edge detangling using only
    // the second-to-last AnchorId of incoming chains and
    // the second AnchorId of outgoing chains.
    bool detangleEdges(
        bool debug,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh,
        bool useBayesianModel,
        double epsilon,
        double minLogP);
    bool detangleEdge(
        bool debug,
        std::map<uint64_t, edge_descriptor>& edgeMap,
        std::map<uint64_t, edge_descriptor>::iterator&,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh,
        bool useBayesianModel,
        double epsilon,
        double minLogP);

    // Edge detangling using up to n AnchorIds
    // of incoming and outgoing chains.
    // This version only handles the 2 by 2 case and always uses the Bayesian model.
    bool detangleEdges(
        bool debug,
        double epsilon,
        double minLogP,
        uint64_t n);
    bool detangleEdge(
        bool debug,
        std::map<uint64_t, edge_descriptor>& edgeMap,
        std::map<uint64_t, edge_descriptor>::iterator&,
        double epsilon,
        double minLogP,
        uint64_t n);



    // Bubble cleanup, with the purpose of eliminating most bubbles caused by errors.
    // See the code for details of what this does.
    uint64_t cleanupBubbles(
        bool debug,
        uint64_t maxOffset,
        uint64_t chainTerminalCommonThreshold,
        uint64_t threadCount);
    uint64_t cleanupBubbles(
        bool debug,
        edge_descriptor ce,
        uint64_t maxOffset,
        uint64_t chainTerminalCommonThreshold);

public:

    // Store Superbubbles information in the superbubbleId field of the vertices.
    // Vertices that don't belong to any superbubble are assigned invalid<uint64_t>.
    void storeSuperbubblesInformation(const Superbubbles&);

    // Find out if two vertices are in the same superbubble.
    bool areInSameSuperbubble(vertex_descriptor v0, vertex_descriptor v1) const
    {
        const AssemblyGraph& assemblyGraph = *this;
        return assemblyGraph[v0].superbubbleId == assemblyGraph[v1].superbubbleId;
    }

    // Find out if the two vertices of an edge are in the same superbubble.
    bool isInternalToSuperbubble(edge_descriptor e) const
    {
        const AssemblyGraph& assemblyGraph = *this;
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        return areInSameSuperbubble(v0, v1);
    }

private:

    // Remove short superbubbles with one entry and one exit.
    bool removeShortSuperbubbles(
        bool debug,
        uint64_t maxOffset1,    // Used to define superbubbles
        uint64_t maxOffset2     // Compared against the offset between entry and exit
    );

    // Detangle short superbubbles with any number of entrances and exits.
    bool detangleShortSuperbubbles(
        bool debug,
        uint64_t maxOffset1,    // Used to define superbubbles
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh,
        bool useBayesianModel,
        double epsilon,
        double minLogP);
    bool detangleShortSuperbubble(
        bool debug,
        const Superbubbles&,
        uint64_t superbubbleId,
        uint64_t maxOffset1,
        uint64_t detangleToleranceLow,
        uint64_t detangleToleranceHigh,
        bool useBayesianModel,
        double epsilon,
        double minLogP);

    // Detangle with read following.
    enum class SuperbubbleCreationMethod {
        SingleEdges,
        ByLength
    };
    uint64_t detangleSuperbubblesWithReadFollowing(
        bool debug,
        SuperbubbleCreationMethod,
        uint64_t maxOffset,
        double maxLoss,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold);
    bool detangleSuperbubbleWithReadFollowing(
        bool debug,
        const Superbubbles&,
        uint64_t superbubbleId,
        uint64_t maxOffset,
        double maxLoss,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold);

    // Detangling with path following.
public:
    void detangle2(); // Python callable.
private:
    void detangle3();

    // Cleanup/simplify superbubbles that are likely to be caused by errors,
    // completely or in part.
    uint64_t cleanupSuperbubbles(
        bool debug,
        uint64_t maxOffset1,    // Used to define superbubbles
        uint64_t maxOffset2,    // Compared against the offset between entry and exit
        uint64_t chainTerminalCommonThreshold);
    bool cleanupSuperbubble(
        bool debug,
        const Superbubbles&,
        uint64_t superbubbleId,
        uint64_t maxOffset2,    // Compared against the offset between entry and exit
        uint64_t chainTerminalCommonThreshold,
        std::set<vertex_descriptor>& previousSuperbubblesVertices);

    // This version of superbubble cleanup uses dominator trees to define superbubbles,
    // instead of computing connected components using edges of length uo tp maxOffset1.
    uint64_t cleanupSuperbubbles(
        bool debug,
        uint64_t maxOffset2,    // Compared against the offset between entry and exit
        uint64_t chainTerminalCommonThreshold);

    // Split terminal haploid bubbles out of bubble chains, to facilitate detangling.
    void splitTerminalHaploidBubbles();
    void splitTerminalHaploidBubbles(edge_descriptor);

    // Phasing of bubble chains using the PhasingGraph.
    void phaseBubbleChainsUsingPhasingGraph(
        bool debug,
        uint64_t n, // Maximum number of Chain AnchorIds to use when computing tangle matrices.
        uint64_t lowThreshold,
        uint64_t highThreshold,
        bool useBayesianModel,
        double epsilon,
        double minLogP,
        uint64_t longBubbleThreshold);
    void phaseBubbleChainUsingPhasingGraph(
        edge_descriptor e,
        uint64_t n, // Maximum number of Chain AnchorIds to use when computing tangle matrices.
        uint64_t lowThreshold,
        uint64_t highThreshold,
        bool useBayesianModel,
        double epsilon,
        double minLogP,
        uint64_t longBubbleThreshold,
        bool debug);
    void phaseBubbleChainUsingPhasedComponents(
        bool debug,
        edge_descriptor e,
        const vector< shared_ptr<PhasedComponent> >&,
        uint64_t longBubbleThreshold);

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
            uint64_t& total,
            double epsilon,
            double& logPin, // log[P(in-phase)/P(random)] in decibels
            double& logPout // log[P(out-of-phase)/P(random)] in decibels
            ) const;
    };


    // Compute the tangle matrix between two incoming chains
    // and two outgoing chains, taking into account up to
    // n AnchorIds for each Chain.
    void computeTangleMatrix(
        const array<const Chain*, 2> inChains,
        const array<const Chain*, 2> outChains,
        uint64_t n,
        TangleMatrix&) const;

    // Gather OrientedReadIds from up to n AnchorIds
    // near the beginning or end of a chain.
    void gatherOrientedReadIdsAtBeginning(
        const Chain&,
        uint64_t n,
        vector<OrientedReadId>&) const;
    void gatherOrientedReadIdsAtEnd(
        const Chain&,
        uint64_t n,
        vector<OrientedReadId>&) const;



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
        double logPInPhase;
        double logPOutOfPhase;
        double logP() const
        {
            return max(max(logPInPhase, logPOutOfPhase), fabs(logPInPhase - logPOutOfPhase));
        }

#if 0
        bool sortByCounts(const PhasingGraphEdge& that) const
        {
            if(maxDiscordant < that.maxDiscordant) {
                return true;
            }
            if(maxDiscordant > that.maxDiscordant) {
                return false;
            }
            return minConcordant > that.minConcordant;
        }
        bool sortByProbabilities(const PhasingGraphEdge& that) const
        {
            return logP() > that.logP();
        }
#endif
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
        void phase(bool debug);
        void phase1(bool debug, bool useBayesianModel);
        bool isConsistent(edge_descriptor) const;
        void writeGraphviz(const string& fileName) const;
        vector< shared_ptr<PhasedComponent> > phasedComponents;

        // Sort edges in order of decreasing significance:
        // - If using the Bayesian model, logP.
        // - Otherwise, minConcordant/maxDiscordant.
        void sortEdges(vector<edge_descriptor>& sortedEdges, bool useBayesianModel) const;
    };



    // Phasing of bubble chains using the PhasingTable.
    void phaseBubbleChainsUsingPhasingTable(
        const string& debugOutputFileNamePrefix,
        double phaseErrorThreshold,
        double bubbleErrorThreshold,
        uint64_t longBubbleThreshold);
    void phaseBubbleChainUsingPhasingTable(
        const string& debugOutputFileNamePrefix,
        edge_descriptor e,
        double phaseErrorThreshold,
        double bubbleErrorThreshold,
        uint64_t longBubbleThreshold);
    void cleanupBubbleChainUsingPhasingTable(
        const string& debugOutputFileNamePrefix,
        edge_descriptor e,
        double phaseErrorThreshold,
        double bubbleErrorThreshold,
        uint64_t longBubbleThreshold);

public:


    // Anchor annotations.
    // We keep track of:
    // - Which AssemblyGraph vertices contain a given AnchorId.
    // - Which Chains contain a given AnchorId at its beginning or end.
    // - Which Chains contain a given AnchorId internally, and at what position.
    //   "Internally" means that we don't consider the first and last AnchorIds
    //   of each chain (those correspond to an AssemblyGraph vertex).
    class AnchorAnnotation {
    public:

        // The vertices that contain this AnchorId.
        vector<vertex_descriptor> vertices;

        // The Chains that contain this AnchorId at the beginning.
        vector<ChainIdentifier> chainsFirstAnchor;

        // The Chains that contain this AnchorId at the end.
        vector<ChainIdentifier> chainsLastAnchor;

        // The Chains that contain this AnchorId internally, and the corresponding
        // position of the AnchorId in each Chain.
        // This stores pairs(ChaiIdentifier, positionInChain).
        vector< pair<ChainIdentifier, uint64_t> > internalChainInfo;
    };

    // Vector of AnchorAnnotations for all anchors in this component, indexed by local anchor id.
    // This is not maintained. It is only created when needed.
    vector<AnchorAnnotation> anchorAnnotations;
    void annotateAnchors();

private:


    // Optimize chains before assembly, to remove assembly steps with
    // less that minCommon reads.
    void optimizeChains(
        bool debug,
        uint64_t minCommon,
        uint64_t k
        );
    void optimizeChain(
        bool debug,
        Chain&,
        uint64_t minCommon,
        uint64_t k
        );

    // Assemble sequence for a single Chain.
    void assembleChain(
        Chain&,
        uint64_t chainTerminalCommonThreshold);

    // Multithreaded version of sequence assembly.
    // This only assembles the chains that have the shouldBeAssembled flag set.
    void assembleChainsMultithreaded(
        uint64_t chainTerminalCommonThreshold,
        uint64_t threadCount);
    // This sets the shouldBeAssembled flag for all chains, then
    // calls assembleChainsMultithreaded.
    void assembleAllChainsMultithreaded(
        uint64_t chainTerminalCommonThreshold,
        uint64_t threadCount);
    // This clears the shouldBeAssembled flag from all Chains.
    void clearAllShouldBeAssembledFlags();
    void cleanupSequence();

    void assembleChainsMultithreadedTheadFunction(uint64_t threadId);
    void combineStepSequences(Chain&);
    class AssemblyStep {
    public:
        edge_descriptor e;              // This identified the BubbleChain.
        uint64_t positionInBubbleChain; // This identifies the Bubble.
        uint64_t indexInBubble;         // This identifies the Chain.
        uint64_t positionInChain;
        uint64_t offsetInBases;

        // For better load balancing, order them by decreasing offsetInBases.
        bool operator<(const AssemblyStep& that) const
        {
            return offsetInBases > that.offsetInBases;
        }
    };
    void runAssemblyStep(
        uint64_t chainTerminalCommonThreshold,
        const AssemblyStep&);
    void runAssemblyStep(
        Chain& chain,
        uint64_t positionInChain,
        uint64_t chainTerminalCommonThreshold);
    class AssembleChainsMultithreadedData {
    public:
        uint64_t chainTerminalCommonThreshold;
        vector<AssemblyStep> assemblySteps;
    };
    AssembleChainsMultithreadedData assembleChainsMultithreadedData;



    // Get the lengths of Chains assembled sequence for each Chain P-value.
    // On return, chainLengths[pValue] contains the lengths of all
    // Chains with that pValue, sorted in decreasing order.
    // This can be used for N50 statistics.
public:
    void getChainLengthsByPValue(vector< vector<uint64_t> >& chainLengths) const;

    // Get the lengths of all non-trivial bubble chains.
    void getBubbleChainLengths(vector<uint64_t>&) const;

    // Given a vector of lengths in decreasing order, compute the total length and N50.
    static pair<uint64_t, uint64_t> n50(const vector<uint64_t>&);
private:

    // Output.
    void write(const string& name, bool writeSequence = false) const;
    void writeCsv(const string& fileNamePrefix) const;
public:
    void writeCsvSummary(ostream&) const;
private:
    void writeBubbleChainsCsv(const string& fileNamePrefix) const;
    void writeBubbleChainsPhasingTables(const string& fileNamePrefix, double phaseErrorThreshold) const;
    void writeBubblesCsv(const string& fileNamePrefix) const;
    void writeChainsCsv(const string& fileNamePrefix) const;
    void writeChainsDetailsCsv(const string& fileNamePrefix) const;
    void writeChainDetailsCsv(ostream&, edge_descriptor, bool writeHeader) const;
    void writeGraphviz(const string& fileNamePrefix, bool labels) const;
    void writeGfa(const string& fileNamePrefix) const;
    void writeGfaExpanded(
        const string& fileNamePrefix,
        bool includeSequence,
        bool useSequenceLength) const;
    void writeGfaExpanded(
        ostream&,
        bool includeSequence,
        bool useSequenceLength) const;
    void writeAssemblyDetails() const;
public:
    void writeGfaSegmentsExpanded(
        ostream&,
        bool includeSequence,
        bool useSequenceLength) const;
    void writeGfaLinksExpanded(ostream&) const;
    static void writeGfaHeader(ostream&);
    void writeFastaExpanded(ostream&) const;
    void writeFastaExpanded(const string& fileNamePrefix) const;
    void writeSnapshot(uint64_t& snapshotNumber) const;

    string bubbleChainStringId(edge_descriptor) const;
    string bubbleStringId(edge_descriptor, uint64_t positionInBubbleChain) const;
    string chainStringId(edge_descriptor, uint64_t positionInBubbleChain, uint64_t indexInBubble) const;


    // Return average coverage for the internal AnchorIds of a Chain.
    // For chain of length 2, this returns 0.
    double primaryCoverage(const Chain&) const;

    // This returns a "P-value" for a Chain defined as follows:
    // If the Chain is the only chain of a BubbleChain, the P-value is 0.
    // Otherwise, the P-value is the ploidy of the Bubble that the Chain belongs to.
    uint64_t chainPValue(edge_descriptor, uint64_t positionInBubbleChain, uint64_t indexInBubble) const;

    uint64_t chainOffset(const Chain&) const;
    void bubbleOffset(
        const Bubble&,
        uint64_t& averageOffset,
        uint64_t& minOffset,
        uint64_t& maxOffset
        ) const;
    bool bubbleOffsetNoException(
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


private:
    // Code in mode3-AssemblyGraph-test.cpp.
    void haplotizeWronglyPolyploidBubbles(
        bool debug
        );

    void removeChainsInBubblesWithNoInternalAnchors(
        bool debug
        );

    void removeSimpleBubbleChainsWithNoInternalAnchors(
        bool debug
        );

    void prune(
        bool debug,
        uint64_t pruneLength
        );

    void removeCrossEdgesInAssemblyGraph(
        bool debug
        );



    // Code in mode3-AssemblyGraphRun4.cpp.
    void run4(
        uint64_t threadCount,
        bool assembleSequence,
        bool debug);
    uint64_t detangle(const Superbubbles&, Detangler&);
    // uint64_t detangleShortSuperbubbles4(bool debug, const Superbubbles&);
    // bool detangleShortSuperbubble4(bool debug, const vector<vertex_descriptor>& superbubble);
    uint64_t cleanupBubbles();
    uint64_t detangleCrossEdgesIndividually(
        bool debug,
        ChainPermutationDetangler&);

    // This detangles induced subgraphs of the AssemblyGraph
    // that are isomorphic to a given Subgraph.
    using Subgraph = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS>;
    uint64_t detangleInducedSubgraphs(
        bool debug,
        const Subgraph&,
        ChainPermutationDetangler&);



    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AssemblyGraphBaseClass>(*this);
        ar & componentId;
        ar & k;
        ar & sequenceWasAssembled;
        ar & nextEdgeId;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph-Stage-ComponentId.
    void save(const string& stage) const;
    void load(const string& stage, uint64_t componentId);

};



// An AssemblyGraphEdgePredicate is an abstract base class for an object
// whose operator() can be applied to an edge_descriptor and can return
// true or false.
class shasta::mode3::AssemblyGraphEdgePredicate {
public:

    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    virtual bool operator()(edge_descriptor) const = 0;

    AssemblyGraphEdgePredicate(const AssemblyGraph&);

protected:

    // All derived classes will need access to the AssemblyGraph.
    const AssemblyGraph& assemblyGraph;
};



// An AssemblyGraphEdgePredicate that returns true for a cross-edge.
// An edge v0->v1 is defined to be a cross-edge if removing it
// does not cause a forward dead end at v0- and a backward dead end
// at v1. That is:
// out-degree(v0) > 1
// in_degree(v1) > 1.
class shasta::mode3::AssemblyGraphCrossEdgePredicate : public AssemblyGraphEdgePredicate {
public:
    bool operator()(edge_descriptor) const;
    AssemblyGraphCrossEdgePredicate(const AssemblyGraph&);
};



// An AssemblyGraphEdgePredicate that returns true for an edge without internal anchors.
// Must be called for an edge corresponding to a single Chain.
class shasta::mode3::AssemblyGraphNoInternalAnchorsEdgePredicate : public AssemblyGraphEdgePredicate {
public:
    bool operator()(edge_descriptor) const;
    AssemblyGraphNoInternalAnchorsEdgePredicate(const AssemblyGraph&);
};

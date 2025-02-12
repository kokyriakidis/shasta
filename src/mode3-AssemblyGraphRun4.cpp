// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detanglers.hpp"
#include "AssemblerOptions.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;



void AssemblyGraph::run4(
    uint64_t threadCount,
    bool /* assembleSequence */,
    bool debug)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double epsilon = 0.05;
    // const double chiSquareThreshold = 30.;
    const uint64_t superbubbleLengthThreshold = 10000;
    const double maxLogP = 10.;
    const double minLogPDelta = 20.;

    AssemblyGraph& assemblyGraph = *this;

    cout << "AssemblyGraph::run4 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    write("A");
    compress();
    write("B");

    // An iteration loop in which we try at each iteration
    // all operations that can simplify the AssemblyGraph.
    for(uint64_t iteration=0; ; iteration++) {
        cout << "Simplify iteration " << iteration << " begins with " <<
            num_vertices(assemblyGraph) << " vertices and " << num_edges(assemblyGraph) <<
            " edges." << endl;

        uint64_t simplificationCount = 0;

        // Bubbles that are removed or cleaned up don't participate in detangling
        // but can still be assembled correctly if we are able to phase/detangle around them.

        // Clean up Bubbles in which one or more Chains have no internal Anchors.
        compress();
        const uint64_t cleanedUpBubbleCount1 = cleanupBubbles();
        cout << "Cleaned up "<< cleanedUpBubbleCount1 << " bubbles containing chains "
            "without internal anchors." << endl;
        simplificationCount += cleanedUpBubbleCount1;

        compressBubbleChains();
        compress();

        // Bubble cleanup. This cleans up Bubbles that are probably caused by errors.
        const uint64_t cleanedUpBubbleCount2 = cleanupBubbles(
            debug,
            options.assemblyGraphOptions.bubbleCleanupMaxOffset,
            options.assemblyGraphOptions.chainTerminalCommonThreshold,
            threadCount);
        simplificationCount += cleanedUpBubbleCount2;
        cout << "Cleaned up " << cleanedUpBubbleCount2 << " bubbles "
            "probably caused by errors." << endl;

        compressBubbleChains();
        compress();

        // Clean up short superbubbles.
        const uint64_t cleanedUpSuperbubbleCount =
            cleanupSuperbubbles(false,
            options.assemblyGraphOptions.superbubbleLengthThreshold1,
            options.assemblyGraphOptions.chainTerminalCommonThreshold);
        simplificationCount += cleanedUpSuperbubbleCount;
        cout << "Cleaned up " << cleanedUpSuperbubbleCount << " superbubbles." << endl;

        compressBubbleChains();
        compress();

        // For detangling the AssemblyGraph needs to be in expanded form.
        expand();

        // Vertex detangling.
        {
            ChainPermutationDetangler detangler(false, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
            Superbubbles superbubbles(assemblyGraph, Superbubbles::FromTangledVertices{});
            const uint64_t detangledCount = detangle(superbubbles, detangler);
            simplificationCount += detangledCount;
            compressSequentialEdges();
            compressBubbleChains();
            cout << "Detangled " << detangledCount << " vertices. The assembly graph now has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " edges." << endl;
        }

        // Edge detangling.
        {
            ChainPermutationDetangler detangler(false, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
            Superbubbles superbubbles(assemblyGraph, Superbubbles::FromTangledEdges{});
            const uint64_t detangledCount = detangle(superbubbles, detangler);
            simplificationCount += detangledCount;
            compressSequentialEdges();
            compressBubbleChains();
            cout << "Detangled " << detangledCount << " edges. The assembly graph now has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " edges." << endl;
        }


        // Detangle superbubbles defined by cross-edges.
        {
            AssemblyGraphCrossEdgePredicate edgePredicate(assemblyGraph);
            Superbubbles superbubbles(assemblyGraph, edgePredicate);

            ChainPermutationDetangler detangler(false, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
            const uint64_t detangledCount = detangle(superbubbles, detangler);

            simplificationCount += detangledCount;
            compressSequentialEdges();
            compressBubbleChains();
            cout << "Detangled " << detangledCount << " superbubbles defined by cross-edges. The assembly graph now has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " edges." << endl;
        }



        // Superbubble detangling.
        {
            ChainPermutationDetangler detangler(false, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
            Superbubbles superbubbles(assemblyGraph, superbubbleLengthThreshold);
            const uint64_t detangledCount = detangle(superbubbles, detangler);
            simplificationCount += detangledCount;
            compressSequentialEdges();
            compressBubbleChains();
            cout << "Detangled " << detangledCount << " superbubbles. The assembly graph now has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " edges." << endl;
        }


        cout << "Detangle iteration " << iteration << " had " <<
            simplificationCount << " successful detangling operations." << endl;

        // After detangling put the AssemblyGraph back in compressed form.
        compress();

        if(simplificationCount == 0) {
            break;
        }
    }
    cout << "After simplifying iterations, the assembly graph has " <<
        num_vertices(assemblyGraph) << " vertices and " << num_edges(assemblyGraph) <<
        " edges." << endl;

    write("Z");

#if 0
    // Assemble sequence.
    assembleAllChainsMultithreaded(
        options.assemblyGraphOptions.chainTerminalCommonThreshold,
        threadCount);
    writeAssemblyDetails();
    write("Final", true);
#endif
}


#if 0
uint64_t AssemblyGraph::detangleShortSuperbubbles4(
    bool debug,
    const Superbubbles& superbubbles)
{
    uint64_t detangledCount = 0;

    for(const Superbubble& superbubble: superbubbles.superbubbles) {
        if(detangleShortSuperbubble4(debug, superbubble)) {
            ++detangledCount;
        }
    }

    return detangledCount;
}
#endif



// Loop over all Superbubbles and let the Detangler try detangling each one.
uint64_t AssemblyGraph::detangle(const Superbubbles& superbubbles, Detangler& detangler)
{
    uint64_t detangledCount = 0;
    for(const Superbubble& superbubble: superbubbles.superbubbles) {
        if(detangler(superbubble)) {
            ++detangledCount;
        }
    }
    return detangledCount;
}



#if 0
// This only detangles 2 by 2 superbubbles.
// It uses a chi-squared test for phasing.
bool AssemblyGraph::detangleShortSuperbubble4(
    bool debug,
    const vector<vertex_descriptor>& superbubble)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double epsilon = 0.05;
    const double chiSquareThreshold = 30.;

    AssemblyGraph& assemblyGraph = *this;

    if(debug) {
        cout << "Found a superbubble with " << superbubble.size() <<
            " vertices:";
        for(const vertex_descriptor cv: superbubble) {
            cout << " " << anchorIdToString(assemblyGraph[cv].getAnchorId());
        }
        cout << endl;
    }

    // Fill in the in-edges and out-edges.
    // These cannot be computed while constructing the superbubbles
    // as they can change when other superbubbles are detangled.
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> outEdges;
    for(const vertex_descriptor v: superbubble) {
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not assemblyGraph.isInternalToSuperbubble(e)) {
                 inEdges.push_back(e);
            }
        }
    }
    for(const vertex_descriptor v: superbubble) {
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            if(not assemblyGraph.isInternalToSuperbubble(e)) {
                 outEdges.push_back(e);
            }
        }
    }
    const uint64_t inDegree = inEdges.size();
    const uint64_t outDegree = outEdges.size();

    if(debug) {
        cout << inDegree << " in-edges:";
        for(const edge_descriptor e: inEdges) {
            cout << " " << bubbleChainStringId(e);
        }
        cout << endl;
        cout << outDegree << " out-edges:";
        for(const edge_descriptor e: outEdges) {
            cout << " " << bubbleChainStringId(e);
        }
        cout << endl;
    }

    // If an inEdge is also an outEdge, don't do anything.
    for(const edge_descriptor e: inEdges) {
        if(find(outEdges.begin(), outEdges.end(), e) != outEdges.end()) {
            if(debug) {
                cout << "Not detangling because " << bubbleChainStringId(e) <<
                    " is both an in-edge and out-edge." << endl;
            }
            return false;
        }
    }

    // This only detangles superbubbles with 2 entrances and 2 exits.
    if(not ((inDegree == 2) and (outDegree == 2))) {
        if(debug) {
            cout << "Not detangling because it is not a 2 by 2 superbubble." << endl;
        }
        return false;
    }

    // Gather the second to last AnchorId of each inEdge and the second AnchorId
    // of each outEdge.
    vector<AnchorId> inAnchors;
    for(const edge_descriptor e: inEdges) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        inAnchors.push_back(chain.secondToLast());
    }
    vector<AnchorId> outAnchors;
    for(const edge_descriptor e: outEdges) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        outAnchors.push_back(chain.second());
    }

    if(debug) {
        cout << inDegree << " in-anchors:";
        for(const AnchorId anchorId: inAnchors) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;
        cout << outDegree << " out-anchors:";
        for(const AnchorId anchorId: outAnchors) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;
    }

    // If an AnchorId appears both in the inAnchors and in the outAnchors,
    // detangling could generate a chain with two consecutive copies of the same
    // AnchorId. Don't detangle.
    for(const AnchorId anchorId: inAnchors) {
        if(find(outAnchors.begin(), outAnchors.end(), anchorId) != outAnchors.end()) {
            if(debug) {
                cout << "Not detangling because " << anchorIdToString(anchorId) <<
                    " is both an in-anchor and out-anchor." << endl;
            }
            return false;
        }
    }


    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix(2, vector<uint64_t>(2));
    uint64_t N = 0;
    vector<uint64_t> inCoverage(inDegree, 0);
    vector<uint64_t> outCoverage(inDegree, 0);
    for(uint64_t i0=0; i0<inDegree; i0++) {
        const AnchorId anchorId0 = inAnchors[i0];
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const AnchorId anchorId1 = outAnchors[i1];
            const uint64_t n = anchors.countCommon(anchorId0, anchorId1, true);
            tangleMatrix[i0][i1] = n;
            N += n;
            inCoverage[i0] += n;
            outCoverage[i1] += n;
        }
    }

    if(debug) {
        cout << "Tangle matrix with total coverage " << N << ":" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrix[i0][i1];

                cout << endl;
            }
        }
        cout << "In-coverage: " << inCoverage[0] << " " << inCoverage[1] << endl;
        cout << "Out-coverage: " << outCoverage[0] << " " << outCoverage[1] << endl;
        SHASTA_ASSERT(inCoverage[0] + inCoverage[1] == N);
        SHASTA_ASSERT(outCoverage[0] + outCoverage[1] == N);
    }

    // If the inCoverage or outCoverage have zero entries, do nothing.
    for(uint64_t i=0; i<inDegree; i++) {
        if(inCoverage[i] == 0) {
            if(debug) {
                cout << "Not detangling because of zero common coverage on in-edge " << bubbleChainStringId(inEdges[i]) << endl;
            }
            return false;
        }
    }
    for(uint64_t i=0; i<outDegree; i++) {
        if(outCoverage[i] == 0) {
            if(debug) {
                cout << "Not detangling because of zero common coverage on out-edge " << bubbleChainStringId(outEdges[i]) << endl;
            }
            return false;
        }
    }


    // Create expected values for the tangle matrix under various assumptions.

    // Random.
    vector< vector<double> > tangleMatrixRandom(2, vector<double>(2, 0.));
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixRandom[i0][i1] = double(inCoverage[i0]) * double(outCoverage[i1]) / double(N);
        }
    }

    // In phase, ideal.
    vector< vector<double> > tangleMatrixInPhaseIdeal(2, vector<double>(2, 0.));
    tangleMatrixInPhaseIdeal[0][0]    = sqrt(double(inCoverage[0]) * double(outCoverage[0]));
    tangleMatrixInPhaseIdeal[1][1]    = sqrt(double(inCoverage[1]) * double(outCoverage[1]));
    // Normalize.
    double inPhaseSum = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            inPhaseSum += tangleMatrixInPhaseIdeal[i0][i1];
        }
    }
    double inPhaseFactor = double(N) / inPhaseSum;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixInPhaseIdeal[i0][i1] *= inPhaseFactor;
        }
    }

    // Out of phase, ideal.
    vector< vector<double> > tangleMatrixOutOfPhaseIdeal(2, vector<double>(2, 0.));
    tangleMatrixOutOfPhaseIdeal[0][1] = sqrt(double(inCoverage[0]) * double(outCoverage[1]));
    tangleMatrixOutOfPhaseIdeal[1][0] = sqrt(double(inCoverage[1]) * double(outCoverage[0]));
    // Normalize.
    double outOfPhaseSum = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            outOfPhaseSum += tangleMatrixOutOfPhaseIdeal[i0][i1];
        }
    }
    double outOfPhaseFactor = double(N) / outOfPhaseSum;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixOutOfPhaseIdeal[i0][i1] *= outOfPhaseFactor;
        }
    }

    // In phase and out of phase, non-ideal.
    vector< vector<double> > tangleMatrixInPhase(2, vector<double>(2, 0.));
    vector< vector<double> > tangleMatrixOutOfPhase(2, vector<double>(2, 0.));
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixInPhase[i0][i1] = epsilon * tangleMatrixRandom[i0][i1] + (1. - epsilon) * tangleMatrixInPhaseIdeal[i0][i1];
            tangleMatrixOutOfPhase[i0][i1] = epsilon * tangleMatrixRandom[i0][i1] + (1. - epsilon) * tangleMatrixOutOfPhaseIdeal[i0][i1];
        }
    }


    if(false) {
        cout << "Random tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixRandom[i0][i1];

                cout << endl;
            }
        }
        cout << "Ideal in phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixInPhaseIdeal[i0][i1];

                cout << endl;
            }
        }
        cout << "Ideal out of phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixOutOfPhaseIdeal[i0][i1];

                cout << endl;
            }
        }
        cout << "Non-ideal in phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixInPhase[i0][i1];

                cout << endl;
            }
        }
        cout << "Non-ideal out of phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixOutOfPhase[i0][i1];

                cout << endl;
            }
        }
    }


    // Do a chi-square test.
    double chi2InPhase = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const double expected = tangleMatrixInPhase[i0][i1];
            const double delta = double(tangleMatrix[i0][i1]) - expected;
            chi2InPhase += delta * delta / expected;
        }
    }
    double chi2OutOfPhase = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const double expected = tangleMatrixOutOfPhase[i0][i1];
            const double delta = double(tangleMatrix[i0][i1]) - expected;
            chi2OutOfPhase += delta * delta / expected;
        }
    }

    if(debug) {
        cout << "Chi square test: in phase " << chi2InPhase <<
            ", out of phase " << chi2OutOfPhase << endl;
    }

    const bool isInPhase = (chi2InPhase < chiSquareThreshold) and (chi2OutOfPhase > chiSquareThreshold);
    const bool isOutOfPhase = (chi2InPhase > chiSquareThreshold) and (chi2OutOfPhase < chiSquareThreshold);

    if(not (isInPhase or isOutOfPhase)) {
        if(debug) {
            cout << "Not detangling bacause phasing is not reliable." << endl;
        }
        return false;
    }

    if(debug) {
        if(isInPhase) {
            cout << "In-phase." << endl;
        }
        if(isOutOfPhase) {
            cout << "Out-of-phase." << endl;
        }
    }



    // We are in-phase or out-of-phase. Generate two new edges.
    for(uint64_t i=0; i<2; i++) {

        // Get the two edges to be connected.
        const edge_descriptor e0 = inEdges[i];
        const edge_descriptor e1 = (isInPhase ? outEdges[i] : outEdges[1 - i]);

        // Get the corresponding Chains.
        const Chain& chain0 = assemblyGraph[e0].getOnlyChain();
        const Chain& chain1 = assemblyGraph[e1].getOnlyChain();

        // Get the two vertices for the new edge.
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        const vertex_descriptor v1 = target(e1, assemblyGraph);

        // Create the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = boost::add_edge(v0, v1, assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.id = nextEdgeId++;
        BubbleChain& newBubbleChain = edgeNew;
        newBubbleChain.resize(1);   // The new BubbleChain has a single Bubble
        Bubble& newBubble = newBubbleChain.front();
        newBubble.resize(1);        // The new Bubble is haploid.
        Chain& newChain = newBubble.front();

        // Build the new chain.
        copy(chain0.begin(), chain0.end() - 1, back_inserter(newChain));
        copy(chain1.begin() + 1, chain1.end(), back_inserter(newChain));
    }

    // Now we can remove  all the vertices inside the superbubble
    // and their edges. This includes the inEdges and outEdges.
    for(const vertex_descriptor v: superbubble) {
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }

    return true;
}
#endif



// This cleans up non-haploid Bubbles in which one or more Chains have no internal Anchors.
uint64_t AssemblyGraph::cleanupBubbles()
{
    AssemblyGraph& assemblyGraph = *this;

    /// Loop over all BubbleChains.
    uint64_t n = 0;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Loop over non-haploid Bubbles of this BubbleChain.
        for(Bubble& bubble: bubbleChain) {
            if(bubble.isHaploid()) {
                continue;
            }

            // Look for a Chain in this Bubble that has no internal Anchors.
            uint64_t j = invalid<uint64_t>;
            for(uint64_t i=0; i<bubble.size(); i++) {
                const Chain& chain = bubble[i];
                if(chain.size() == 2) {
                    j = i;
                    break;
                }
            }

            // If did not find any such Chains, do nothing.
            if(j == invalid<uint64_t>) {
                continue;
            }

            // Only keep the Chain without internal Anchor.
            const Chain chain = bubble[j];
            bubble.clear();
            bubble.push_back(chain);
            ++n;
        }
    }

    return n;
}

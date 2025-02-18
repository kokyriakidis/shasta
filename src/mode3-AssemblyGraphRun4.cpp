// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detanglers.hpp"
#include "AssemblerOptions.hpp"
#include "deduplicate.hpp"
#include "inducedSubgraphIsomorphisms.hpp"
#include "orderPairs.hpp"
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

    // write("A");
    compress();
    // write("B");

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

        // Detangle of cross-edges individually.
        {
            ChainPermutationDetangler detangler(iteration==4, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
            const uint64_t detangledCount = detangleCrossEdgesIndividually(iteration==4, detangler);
            simplificationCount += detangledCount;
            compressSequentialEdges();
            compressBubbleChains();
            cout << "Detangled " << detangledCount << " cross-edges individually. The assembly graph now has " <<
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

        // Detangle superbubbles defined by edges without internal anchors.
        {
            AssemblyGraphNoInternalAnchorsEdgePredicate edgePredicate(assemblyGraph);
            Superbubbles superbubbles(assemblyGraph, edgePredicate);

            ChainPermutationDetangler detangler(false, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
            const uint64_t detangledCount = detangle(superbubbles, detangler);

            simplificationCount += detangledCount;
            compressSequentialEdges();
            compressBubbleChains();
            cout << "Detangled " << detangledCount << " superbubbles defined by edges without internal anchors. The assembly graph now has " <<
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

    // write("Z");



    // Test detangleInducedSubgraphs.
    {
        expand();

        // A diploid bubble preceded and followed by haploid segments.
        Subgraph subgraph(4);
        add_edge(0, 1, subgraph);
        add_edge(1, 2, subgraph);
        add_edge(1, 2, subgraph);
        add_edge(2, 3, subgraph);
        ChainPermutationDetangler detangler(false, assemblyGraph, 6, epsilon, maxLogP, minLogPDelta);
        write("U");
        while(true) {
            const uint64_t detangledCount =  detangleInducedSubgraphs(false, subgraph, detangler);
            cout << "Detangled " << detangledCount << " induced subgraphs." << endl;
            if(detangledCount == 0) {
                break;
            }
        }
        write("V");
        compress();
    }



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



// This detangles cross-edges individually one by one.
// For each cross-edge, we construt a Superbubbles object
// containing just a single Superbubble with the two vertices
// of the cross-edge. We can't use a single Superbubbles object
// with one Superbubble for each cross-edge because
// one vertex can only belong to one superbubble.
uint64_t AssemblyGraph::detangleCrossEdgesIndividually(
    bool debug,
    ChainPermutationDetangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    // While detangling, some edges are removed.
    // To avoid attempting to detangle edges that were removed
    // during detangling, we maintain a map that contains all
    // the unprocessed edges. The map gives the edge_descriptor
    // corresponding to each edgeId. We remove map entries as edges
    // are removed during detangling.
    // We cannot use edge_descriptors are reliable keys because
    // an edge_descriptor can be reused after an edge is removed.
    std::map<uint64_t, edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
         edgeMap.insert(make_pair(assemblyGraph[e].id, e));
    }

    // Clear the superbubbleId of all vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        assemblyGraph[v].superbubbleId = invalid<uint64_t>;
    }

    // Create a Superbubbles object containing a single Superbubble
    // with two vertices. We will these vertices in for each cross edge separately.
    Superbubbles superbubbles(assemblyGraph, Superbubbles::Empty());
    superbubbles.superbubbles.resize(1);
    Superbubble& superbubble = superbubbles.superbubbles.front();
    superbubble.resize(2);



    // Main loop over unprocessed edges.
    AssemblyGraphCrossEdgePredicate crossEdgePredicate(assemblyGraph);
    uint64_t detangledCount = 0;
    vector<uint64_t> detangledEdges;
    while(not edgeMap.empty()) {

        // Get the first edge in the edgeMap and remove it from the edgeMap.
        auto it = edgeMap.begin();
        // const uint64_t edgeId = it->first;
        const edge_descriptor e = it->second;
        edgeMap.erase(it);

        // If not a cross-edge, do nothing.
        if(not crossEdgePredicate(e)) {
            continue;
        }

        if(debug) {
            cout << "Detangling cross-edge " << bubbleChainStringId(e) << endl;
        }

        // Get the vertices of our cross-edge.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        // If the detangling is successful, all the in-edges and out-edges
        // of v0 and  v1 will be removed, and we would have to remove them from
        // the edgeMap. Find out what those edges are.
        detangledEdges.clear();
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            detangledEdges.push_back(assemblyGraph[e].id);
        }
        BGL_FORALL_INEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            detangledEdges.push_back(assemblyGraph[e].id);
        }
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            detangledEdges.push_back(assemblyGraph[e].id);
        }
        BGL_FORALL_OUTEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            detangledEdges.push_back(assemblyGraph[e].id);
        }
        deduplicate(detangledEdges);

        // Fill in our Superbubble.
        superbubble[0] = v0;
        superbubble[1] = v1;
        assemblyGraph[v0].superbubbleId = 0;
        assemblyGraph[v1].superbubbleId = 0;

        // Attempt to detangle it.
        const uint64_t success = detangle(superbubbles, detangler);


        if(success) {
            ++detangledCount;
            for(const uint64_t edgeId: detangledEdges) {
                edgeMap.erase(edgeId);
            }
        } else {
            assemblyGraph[v0].superbubbleId = invalid<uint64_t>;
            assemblyGraph[v1].superbubbleId = invalid<uint64_t>;
        }
    }


    return detangledCount;
}



// This detangles induced subgraphs of the AssemblyGraph
// that are isomorphic to a given Subgraph.
uint64_t AssemblyGraph::detangleInducedSubgraphs(
    bool debug,
    const Subgraph& subgraph,
    ChainPermutationDetangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    vector< vector<vertex_descriptor> > isomorphisms;
    inducedSubgraphIsomorphisms(assemblyGraph, subgraph, isomorphisms);



    if(debug) {
        cout << "Found " << isomorphisms.size() << " isomorphisms." << endl;

        // Find sets of parallel edges in the Subgraph.
        // They are needed below.
        vector< pair<Subgraph::vertex_descriptor, Subgraph::vertex_descriptor> > parallelEdgesSet;
        BGL_FORALL_EDGES(e, subgraph, Subgraph) {
            const Subgraph::vertex_descriptor v0 = source(e, subgraph);
            const Subgraph::vertex_descriptor v1 = target(e, subgraph);
            parallelEdgesSet.push_back({v0, v1});
        }
        deduplicate(parallelEdgesSet);

        // Loop over all isomorphisms.
        for(uint64_t i=0; i<isomorphisms.size(); i++) {
            const vector<vertex_descriptor>& isomorphism = isomorphisms[i];

            // Write the vertex isomorphism.
            cout << "Isomorphism " << i << endl;
            cout << "Vertex isomorphism:" << endl;
            for(Subgraph::vertex_descriptor v0s=0; v0s<isomorphism.size(); v0s++) {
                const vertex_descriptor v0 = isomorphism[v0s];
                cout << v0s << ": " << anchorIdToString(assemblyGraph[v0].anchorId) << endl;
            }

            // The edge isomorphism is more complicated because we have to
            // account correctly for parallel edges.


            // Find the corresponding sets of parallel edges in the AssemblyGraph.
            cout << "Edge isomorphism:" << endl;
            for(const pair<Subgraph::vertex_descriptor, Subgraph::vertex_descriptor>& p: parallelEdgesSet) {
                const Subgraph::vertex_descriptor v0s = p.first;
                const Subgraph::vertex_descriptor v1s = p.second;
                const vertex_descriptor v0 = isomorphism[v0s];
                const vertex_descriptor v1 = isomorphism[v1s];
                cout << v0s << "->" << v1s << ":";
                BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                    if(target(e, assemblyGraph) == v1) {
                        cout << " " << bubbleChainStringId(e);
                    }
                }
                cout << endl;
            }
        }
    }



    // The induced subgraphs can overlap each other.
    // If we successfully detangle one induced subgraph, we have to discard
    // all of the overlapping induced subgraphs.
    // We create an undirected graph in which each vertex represents
    // one of the isomorphisms we found (that is, one of the induced subgragphs).
    // We add edges between vertices corresponding to isomorphisms that share
    // AssemblyGraph vertices.
    using OverlapGraph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;
    OverlapGraph overlapGraph(isomorphisms.size());
    {
        // Find the isomorphisms that each AssemblyGraph vertex is involved in.
        std::map<vertex_descriptor, vector<uint64_t> > vertexIsomorphisms;
        for(uint64_t i=0; i<isomorphisms.size(); i++) {
            const vector<vertex_descriptor>& isomorphism = isomorphisms[i];
            for(Subgraph::vertex_descriptor v0s=0; v0s<isomorphism.size(); v0s++) {
                const vertex_descriptor v0 = isomorphism[v0s];
                vertexIsomorphisms[v0].push_back(i);
            }
        }

        // Add edges to the OverlapGraph.
        for(const auto& p: vertexIsomorphisms) {
            const vector<uint64_t>& v = p.second;
            for(uint64_t i=0; i<v.size()-1; i++) {
                for(uint64_t j=i+1; j<v.size(); j++) {
                    add_edge(v[i], v[j], overlapGraph);
                }
            }
        }

        if(debug) {
            ofstream dot("OverlapGraph.dot");
            dot << "graph OverlapGraph{\n";
            BGL_FORALL_VERTICES(v, overlapGraph, OverlapGraph) {
                dot << v << ";\n";
            }
            BGL_FORALL_EDGES(e, overlapGraph, OverlapGraph) {
                const OverlapGraph::vertex_descriptor v0 = source(e, overlapGraph);
                const OverlapGraph::vertex_descriptor v1 = target(e, overlapGraph);
                dot << v0 << "--" << v1 << ";\n";
            }
            dot << "}\n";
        }

    }


    // Each induced subgraph isomorphism corresponds to a superbubble.
    // Each of these superbubbles is detangled individually.
    // However if an induced subgraph is detangled successfully we
    // have to skip all other induced subgraphs that overlap with it.
    vector<bool> wasDetangledSuccessfully(isomorphisms.size(), false);

    // Clear the superbubbleId of all vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        assemblyGraph[v].superbubbleId = invalid<uint64_t>;
    }

    // Create a Superbubbles object containing a single Superbubble.
    Superbubbles superbubbles(assemblyGraph, Superbubbles::Empty());
    superbubbles.superbubbles.resize(1);
    Superbubble& superbubble = superbubbles.superbubbles.front();



    // Main loop over induced subgraph isomorphisms.
    uint64_t detangledCount = 0;
    for(uint64_t i=0; i<isomorphisms.size(); i++) {
        // cout << "Detangling " << i << " of " << isomorphisms.size() << endl;

        // If any overlapping isomorphisms were successfully detangled,
        // we have to skip this one.
        bool skip = false;
        BGL_FORALL_OUTEDGES(i, e, overlapGraph, OverlapGraph) {
            SHASTA_ASSERT(source(e, overlapGraph) == i);
            const uint64_t j = target(e, overlapGraph);
            if(wasDetangledSuccessfully[j]) {
                skip = true;
                break;
            }
        }
        if(skip) {
            // cout << "Skipped " << i << endl;
            continue;
        }


        // Fill in our Superbubble.
        superbubble.clear();
        const vector<vertex_descriptor>& isomorphism = isomorphisms[i];
        for(const vertex_descriptor v: isomorphism) {
            superbubble.push_back(v);
            assemblyGraph[v].superbubbleId = 0;
        }

        // Attempt to detangle it.
        wasDetangledSuccessfully[i] = detangle(superbubbles, detangler);


        if(wasDetangledSuccessfully[i]) {
            ++detangledCount;
            // cout << "Success " << i << endl;
        } else {
            // cout << "Failure " << i << endl;
            for(const vertex_descriptor v: isomorphism) {
                superbubble.push_back(v);
                assemblyGraph[v].superbubbleId = invalid<uint64_t>;
            }
        }
    }

    return detangledCount;
}

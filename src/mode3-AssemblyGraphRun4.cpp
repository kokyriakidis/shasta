// Shasta.
#include "mode3-AssemblyGraph.hpp"
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
    AssemblyGraph& assemblyGraph = *this;

    cout << "AssemblyGraph::run4 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    // write("A");

    // Bubble cleanup.
    compress();
    for(uint64_t iteration=0; ; iteration ++) {
        performanceLog << timestamp << "Iteration " << iteration <<
            " of bubble cleanup begins." << endl;
        const uint64_t cleanedUpBubbleCount = cleanupBubbles(
            debug,
            options.assemblyGraphOptions.bubbleCleanupMaxOffset,
            options.assemblyGraphOptions.chainTerminalCommonThreshold,
            threadCount);
        if(cleanedUpBubbleCount == 0) {
            break;
        }
        cout << "Cleaned up " << cleanedUpBubbleCount << " bubbles." << endl;
        compressBubbleChains();
        compress();
    }

    // Clean up short superbubbles.
    cleanupSuperbubbles(false,
        options.assemblyGraphOptions.superbubbleLengthThreshold1,
        options.assemblyGraphOptions.chainTerminalCommonThreshold);
    compressBubbleChains();
    compress();

    // For detangling the Assemblyraph needs to be in expanded form.
    expand();

    // Edge detangling.
    Superbubbles superbubbles(assemblyGraph, Superbubbles::FromEdges{});
    detangleShortSuperbubbles4(false, superbubbles);

    // Before assembly sequence we put the AssemblyGraph back to compressed form.
    compress();

    // Assemble sequence.
    assembleAllChainsMultithreaded(
        options.assemblyGraphOptions.chainTerminalCommonThreshold,
        threadCount);
    writeAssemblyDetails();
    write("Final", true);
}



uint64_t AssemblyGraph::detangleShortSuperbubbles4(
    bool debug,
    const Superbubbles& superbubbles)
{
    uint64_t n = 0;

    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(detangleShortSuperbubble4(debug, superbubbles, superbubbleId)) {
            ++n;
        }
    }

    return n;
}



// This only detangles 2 by 2 superbubbles.
// It uses a chi-squared test for phasing.
bool AssemblyGraph::detangleShortSuperbubble4(
    bool debug,
    const Superbubbles& superbubbles,
    uint64_t superbubbleId)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double epsilon = 0.05;
    const double chiSquareThreshold = 10.;

    AssemblyGraph& assemblyGraph = *this;
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

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
    for(const vertex_descriptor v0: superbubble) {
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = source(e, assemblyGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, v1)) {
                 inEdges.push_back(e);
            }
        }
    }
    for(const vertex_descriptor v0: superbubble) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, v1)) {
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
            const uint64_t n = anchors.countCommon(anchorId0, anchorId1);
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

    return false;
}


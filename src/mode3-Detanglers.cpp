// This file contains implementations of concrete classes derived from mode3::Detangler.

// Shasta.
#include "mode3-Detanglers.hpp"
#include "mode3-Anchor.hpp"
#include "mode3-AssemblyGraph.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <iostream.hpp>



Detangler2by2::Detangler2by2(
    bool debug,
    AssemblyGraph& assemblyGraph,
    double epsilon,
    double chiSquareThreshold) :
    ChainDetangler(debug, assemblyGraph),
    debug(debug),
    epsilon(epsilon),
    chiSquareThreshold(chiSquareThreshold)
{
}



bool Detangler2by2::operator()(const vector<vertex_descriptor>& superbubble)
{
    writeInitialMessage(superbubble);
    prepare(superbubble);
    writeEntrancesAndExits();

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
    const uint64_t inDegree = entrances.size();
    const uint64_t outDegree = exits.size();

    writeEntrancesAndExits();

    // This only detangles superbubbles with 2 entrances and 2 exits.
    if(not ((inDegree == 2) and (outDegree == 2))) {
        if(debug) {
            cout << "Not detangling because it is not a 2 by 2 superbubble." << endl;
        }
        return false;
    }

    // If there are common Chains between the entrance and exits, don't do anything.
    if(commonChainsBetweenEntrancesAndExitsExists()) {
        if(debug) {
            cout << "Not detangling because of a cycle between the entrances and exits." << endl;
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

    // If there are common Anchors between the entrance and exits, don't do anything.
    // Detangling could generate Chain with consecutive identical AnchorIds.
    if(commonAnchorsBetweenEntrancesAndExitsExists()) {
        if(debug) {
            cout << "Not detangling because of a common anchors between the entrances and exits." << endl;
        }
        return false;
    }


    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix(inDegree, vector<uint64_t>(outDegree));
    uint64_t N = 0;
    for(uint64_t iEntrance=0; iEntrance<inDegree; iEntrance++) {
        Entrance& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<outDegree; iExit++) {
            Exit& exit = exits[iExit];
            const uint64_t n = assemblyGraph.anchors.countCommon(entrance.anchorId, exit.anchorId, true);
            tangleMatrix[iEntrance][iExit] = n;
            N += n;
            entrance.commonCoverage += n;
            exit.commonCoverage += n;
        }
    }

    if(debug) {
        cout << "Tangle matrix with total coverage " << N << ":" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << assemblyGraph.bubbleChainStringId(inEdge) << " " <<
                    assemblyGraph.bubbleChainStringId(outEdge) << " " << tangleMatrix[i0][i1];

                cout << endl;
            }
        }
        cout << "In-coverage: " << entrances[0].commonCoverage << " " << entrances[1].commonCoverage << endl;
        cout << "Out-coverage: " << exits[0].commonCoverage << " " << exits[1].commonCoverage << endl;
        SHASTA_ASSERT(entrances[0].commonCoverage + entrances[1].commonCoverage == N);
        SHASTA_ASSERT(exits[0].commonCoverage + exits[1].commonCoverage == N);
    }

    // If any Entrance or Exit has zero common coverage, do nothing.
    for(const Entrance& entrance: entrances) {
        if(entrance.commonCoverage == 0) {
            if(debug) {
                cout << "Not detangling because of zero common coverage on entrance " <<
                    assemblyGraph.bubbleChainStringId(entrance.e) << endl;
            }
            return false;
        }
    }
    for(const Exit& exit: exits) {
        if(exit.commonCoverage == 0) {
            if(debug) {
                cout << "Not detangling because of zero common coverage on exit " <<
                    assemblyGraph.bubbleChainStringId(exit.e) << endl;
            }
            return false;
        }
    }


    // Create expected values for the tangle matrix under various assumptions.

    // Random.
    vector< vector<double> > tangleMatrixRandom(2, vector<double>(2, 0.));
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixRandom[i0][i1] =
                double(entrances[i0].commonCoverage) * double(exits[i1].commonCoverage) / double(N);
        }
    }

    // In phase, ideal.
    vector< vector<double> > tangleMatrixInPhaseIdeal(2, vector<double>(2, 0.));
    tangleMatrixInPhaseIdeal[0][0] = sqrt(double(entrances[0].commonCoverage) * double(exits[0].commonCoverage));
    tangleMatrixInPhaseIdeal[1][1] = sqrt(double(entrances[1].commonCoverage) * double(exits[1].commonCoverage));
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
    tangleMatrixOutOfPhaseIdeal[0][1] = sqrt(double(entrances[0].commonCoverage) * double(exits[1].commonCoverage));
    tangleMatrixOutOfPhaseIdeal[1][0] = sqrt(double(entrances[1].commonCoverage) * double(exits[0].commonCoverage));
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

                cout << assemblyGraph.bubbleChainStringId(inEdge) << " " <<
                    assemblyGraph.bubbleChainStringId(outEdge) << " " << tangleMatrixRandom[i0][i1];

                cout << endl;
            }
        }
        cout << "Ideal in phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << assemblyGraph.bubbleChainStringId(inEdge) << " " <<
                    assemblyGraph.bubbleChainStringId(outEdge) << " " << tangleMatrixInPhaseIdeal[i0][i1];

                cout << endl;
            }
        }
        cout << "Ideal out of phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << assemblyGraph.bubbleChainStringId(inEdge) << " " <<
                    assemblyGraph.bubbleChainStringId(outEdge) << " " << tangleMatrixOutOfPhaseIdeal[i0][i1];

                cout << endl;
            }
        }
        cout << "Non-ideal in phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << assemblyGraph.bubbleChainStringId(inEdge) << " " <<
                    assemblyGraph.bubbleChainStringId(outEdge) << " " << tangleMatrixInPhase[i0][i1];

                cout << endl;
            }
        }
        cout << "Non-ideal out of phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << assemblyGraph.bubbleChainStringId(inEdge) << " " <<
                    assemblyGraph.bubbleChainStringId(outEdge) << " " << tangleMatrixOutOfPhase[i0][i1];

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
        edgeNew.id = assemblyGraph.nextEdgeId++;
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

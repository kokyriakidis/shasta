// This file contains implementations of concrete classes derived from mode3::Detangler.

// Shasta.
#include "mode3-Detanglers.hpp"
#include "mode3-Anchor.hpp"
#include "mode3-AssemblyGraph.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// Standard library.
#include <iostream.hpp>



ChainPermutationDetangler::ChainPermutationDetangler(
    bool debug,
    AssemblyGraph& assemblyGraph,
    uint64_t nMax,
    double epsilon,
    double maxLogP,
    double minLogPDelta,
    uint64_t minDetangledCoverage) :
    ChainDetangler(debug, assemblyGraph),
    debug(debug),
    nMax(nMax),
    epsilon(epsilon),
    maxLogP(maxLogP),
    minLogPDelta(minLogPDelta),
    minDetangledCoverage(minDetangledCoverage)
{
}



bool ChainPermutationDetangler::operator()(const vector<vertex_descriptor>& superbubble)
{

    writeInitialMessage(superbubble);
    prepare(superbubble);
    writeEntrancesAndExits();
    writeTangleMatrix();

    const uint64_t inDegree = entrances.size();
    const uint64_t outDegree = exits.size();

    // Check the degrees.
    if(inDegree != outDegree) {
        if(debug) {
            cout << "Not detangling because the numbers of entrances and exits are not equal." << endl;
        }
        return false;
    }
    const uint64_t n = inDegree;
    if(n < 2) {
        if(debug) {
            cout << "Not detangling because the number of entrances and exits are too low." << endl;
        }
        return false;
    }
    if(n > nMax) {
        if(debug) {
            cout << "Not detangling because the number of entrances and exits are too high." << endl;
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

    // If there are common Anchors between the entrance and exits, don't do anything.
    // Detangling could generate Chain with consecutive identical AnchorIds.
    if(commonAnchorsBetweenEntrancesAndExitsExists()) {
        if(debug) {
            cout << "Not detangling because of a common anchors between the entrances and exits." << endl;
        }
        return false;
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



    // Chi squared distribution used below to test each permutation.
    boost::math::chi_squared_distribution chi2Distribution(double(n * n - 1));

    // Permutations and their logP.
    // Only store the ones for which logP is not inf.
    vector< pair<vector<uint64_t>, double> > permutationTable;



    // Try all possible permutations of the Exits, leaving the Entrances in the same order.
    // For each permutation, compute the permuted tangle matrix and perform a chi square test on it.
    vector<uint64_t> permutation(n);
    std::iota(permutation.begin(), permutation.end(), 0);
    vector<Exit> permutedExits(n);
    vector< vector<uint64_t> > permutedTangleMatrix(n, vector<uint64_t>(n));
    do {
        if(debug) {
            cout << "Trying permutation ";
            copy(permutation.begin(), permutation.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
        }

        // Store the permuted Exits.
        for(uint64_t iExit=0; iExit<n; iExit++) {
            permutedExits[iExit] = exits[permutation[iExit]];
        }
        if(debug) {
            cout << "Permuted exits:";
            for(const Exit& exit: permutedExits) {
                cout << " " << assemblyGraph.bubbleChainStringId(exit.e);
            }
            cout << endl;
        }

        // Compute the permuted tangle matrix.
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                permutedTangleMatrix[iEntrance][iExit] = tangleMatrix[iEntrance][permutation[iExit]];
            }
        }

        if(debug) {
            cout << "Permuted tangle matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                const Entrance& entrance = entrances[iEntrance];

                for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                    const Exit& exit = permutedExits[iExit];

                    cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
                        assemblyGraph.bubbleChainStringId(exit.e) <<
                        " " << permutedTangleMatrix[iEntrance][iExit];

                    cout << endl;
                }
            }
        }


        // If any diagonal element of the permuted tangle matrix is
        // less than minDetangledCoverage, ignore this permutation.
        bool hasSmallDiagonalElements = false;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                if(permutedTangleMatrix[iEntrance][iExit] < minDetangledCoverage) {
                    hasSmallDiagonalElements = true;
                }
            }
        }
        if(hasSmallDiagonalElements) {
            if(debug) {
                cout << "This permutation was discarded because of small diagonal elements." << endl;
            }
            continue;
        }


        // Assuming this permutation, create expected values for the tangle matrix
        // under various assumptions.

        // Random.
        vector< vector<double> > tangleMatrixRandom(n, vector<double>(n, 0.));
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                tangleMatrixRandom[iEntrance][iExit] =
                    double(entrances[iEntrance].commonCoverage) *
                    double(permutedExits[iExit].commonCoverage) /
                    double(totalCommonCoverage);
            }
        }

        // Assuming this permutation, the ideal tangle matrix is diagonal.
        vector< vector<double> > tangleMatrixIdeal(n, vector<double>(n, 0.));
        for(uint64_t i=0; i<n; i++) {
            tangleMatrixIdeal[i][i] = 0.5 * (double(entrances[i].commonCoverage) + double(permutedExits[i].commonCoverage));
        }

        // Compute the expected tangle matrix under non-ideal assumptions.
        vector< vector<double> > tangleMatrixNonIdeal(n, vector<double>(n, 0.));
        for(uint64_t iEntrance=0; iEntrance<inDegree; iEntrance++) {
            for(uint64_t iExit=0; iExit<outDegree; iExit++) {
                tangleMatrixNonIdeal[iEntrance][iExit] =
                    epsilon * tangleMatrixRandom[iEntrance][iExit] +
                    (1. - epsilon) * tangleMatrixIdeal[iEntrance][iExit];
            }
        }

        if(debug) {
            cout << "Random tangle matrix for this permutation:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                const Entrance& entrance = entrances[iEntrance];

                for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                    const Exit& exit = permutedExits[iExit];

                    cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
                        assemblyGraph.bubbleChainStringId(exit.e) <<
                        " " << tangleMatrixRandom[iEntrance][iExit];

                    cout << endl;
                }
            }

            cout << "Expected ideal tangle matrix for this permutation:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                const Entrance& entrance = entrances[iEntrance];

                for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                    const Exit& exit = permutedExits[iExit];

                    cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
                        assemblyGraph.bubbleChainStringId(exit.e) <<
                        " " << tangleMatrixIdeal[iEntrance][iExit];

                    cout << endl;
                }
            }

            cout << "Expected non-ideal tangle matrix for this permutation:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                const Entrance& entrance = entrances[iEntrance];

                for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                    const Exit& exit = permutedExits[iExit];

                    cout << assemblyGraph.bubbleChainStringId(entrance.e) << " " <<
                        assemblyGraph.bubbleChainStringId(exit.e) <<
                        " " << tangleMatrixNonIdeal[iEntrance][iExit];

                    cout << endl;
                }
            }
        }

        // Do a chi-square test.
        double chi2 = 0.;
        for(uint64_t iEntrance=0; iEntrance<inDegree; iEntrance++) {
            for(uint64_t iExit=0; iExit<outDegree; iExit++) {
                const double expected = tangleMatrixNonIdeal[iEntrance][iExit];
                const double delta = double(permutedTangleMatrix[iEntrance][iExit]) - expected;
                chi2 += delta * delta / expected;
            }
        }
        const double cumulativeChi2 = cdf(complement(chi2Distribution, chi2));
        if(cumulativeChi2 == 0.) {
            if(debug) {
                cout << "Chi square cumulative is zero." << endl;
            }
        }

        if(cumulativeChi2 > 0.) {
            // Compute the corresponding logP in decibels.
            const double logP = - 10. * log10(cdf(complement(chi2Distribution, chi2)));

            // Store this permutation and its logP.
            permutationTable.push_back(make_pair(permutation, logP));

            if(debug) {
                cout << "Chi square for this permutation: " << chi2 <<
                    ", logP " << logP << " dB." << endl;
            }
        }

    } while(std::next_permutation(permutation.begin(), permutation.end()));


    // If we did not find any usable permutation, do nothing.
    if(permutationTable.empty()) {
        if(debug) {
            cout << "The permutation table is empty." << endl;
        }
        return false;
    }


    // Sort the permutation table by increasing logP.
    sort(permutationTable.begin(), permutationTable.end(),
        OrderPairsBySecondOnly<vector<uint64_t>, double>());

    if(debug) {
        cout << "Permutation table:" << endl;
        for(const auto& p: permutationTable) {
            const vector<uint64_t>& permutation = p.first;
            const double logP = p.second;

            for(uint64_t i=0; i<n; i++) {
                cout << "(" <<
                    assemblyGraph.bubbleChainStringId(entrances[i].e) << "," <<
                    assemblyGraph.bubbleChainStringId(exits[permutation[i]].e) << ") ";
            }
            cout << logP << endl;
        }
    }

    // If the best logP is not good enough, do nothing.
    const double bestLogP = permutationTable.front().second;
    if(bestLogP > maxLogP) {
        if(debug) {
            cout << "Not detangling because the best logP is too large." << endl;
        }
        return false;
    }

    // Also check the second best logP.
    if(permutationTable.size() >= 2) {
        const double secondBestLogP = permutationTable[1].second;
        const double logPDelta = secondBestLogP - bestLogP;
        if(logPDelta < minLogPDelta) {
            if(debug) {
                cout << "Not detangling because the second best logP is too small." << endl;
            }
            return false;
        }

    }

    if(debug) {
        cout << "This superbubble will be detangled." << endl;
    }

    // Create new Chains by connecting each Entrance with the
    // corresponging Exit according to the bets permutation.
    const vector<uint64_t>& bestPermutation = permutationTable.front().first;
    for(uint64_t i=0; i<n; i++) {
        const Entrance& entrance = entrances[i];
        const Exit& exit = exits[bestPermutation[i]];
        connect(entrance, exit);
    }

    // Now we can remove all the superbubble vertices.
    removeAllSuperbubbleVertices(superbubble);

    return true;    // Detangling was successful.
}

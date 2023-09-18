#ifndef SHASTA_GLOBAL_MSA_HPP
#define SHASTA_GLOBAL_MSA_HPP

/*******************************************************************************

Global multiple sequence alignment.
Global means constrained on both sides, aka Needleman–Wunsch.

This supports sequences of arbitrary length.
If all the sequences are at most maxSpoaLength long,
this invokes spoa.

Otherwise it finds a common subsequence of length kmerLength
and splits the MSA at that location, invoking itself recursively
to solve the two MSAs.

Each of the input sequences is passed in as a pair.
The second member of the pair is the "weight" of the sequence
(that is, typically the number of reads with that sequence).

*******************************************************************************/

#include "cstdint.hpp"
#include "utility.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta {

    class Base;
    class AlignedBase;

    void globalMsa(
        const vector< pair<vector<Base>, uint64_t> >& sequences,
        uint64_t maxSpoaLength,
        uint64_t kmerLength,
        vector<Base>& consensus
        );

    void globalMsaSpoa(
        const vector< pair<vector<Base>, uint64_t> >& sequences,
        vector<Base>& consensus
        );
    bool globalMsaSpoa(
        const vector< pair<vector<Base>, uint64_t> >& sequences,
        vector<Base>& consensus,
        uint64_t maximumMsaLength
        );
    void globalMsaSpoa(
        const vector< pair<vector<Base>, uint64_t> >& sequences,
        vector< vector<AlignedBase> >& alignment
        );

    // Python-callable version.
    string globalMsaPython(
        const vector< pair<string, uint64_t> >& sequenceStrings,
        uint64_t maxSpoaLength,
        uint64_t kmerLength
    );
}

#endif

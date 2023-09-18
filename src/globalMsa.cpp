// Shasta.
#include "globalMsa.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
#include "ShortBaseSequence.hpp"

// Spoa.
#include "spoa/spoa.hpp"

// Standard library.
#include "algorithm.hpp"
#include <map>
#include "tuple.hpp"

// See the comments in globalMsa.hpp.



void shasta::globalMsa(
    const vector< pair<vector<Base>, uint64_t> >& sequences,
    uint64_t maxSpoaLength,
    uint64_t kmerLength,
    vector<Base>& consensus
    )
{
    const bool debug = true;
    if(debug) {
        cout << "globalMsa called with " << sequences.size() << " sequences with (length,weight):" << endl;
        uint64_t totalWeight = 0;
        for(const auto& p: sequences) {
            cout << "(" << p.first.size() << "," << p.second << ") ";
            totalWeight += p.second;
        }
        cout << endl;
        cout << "Total weight is " << totalWeight << endl;
    }

    // Sanity check.
    SHASTA_ASSERT(not sequences.empty());

    using Kmer = ShortBaseSequence64;
    SHASTA_ASSERT(kmerLength <= Kmer::capacity);

    // Trivial case.
    if(sequences.size() == 1) {
        consensus = sequences.front().first;
        return;
    }

    // Compute the maximum length of the input sequences.
    uint64_t maxLength = 0;
    for(const auto& p: sequences) {
        maxLength = max(maxLength, p.first.size());
    }

    // If short enough, use spoa.
    if(maxLength <= maxSpoaLength) {
        if(debug) {
            cout << "Using spoa." << endl;
        }
        globalMsaSpoa(sequences, consensus);
        return;
    }



    // Create a table of unique k-mers for each of the sequences.
    class KmerInfo {
    public:
        Kmer kmer;
        uint64_t position;
        bool operator<(const KmerInfo& that) const
        {
            return kmer.data < that.kmer.data;
        }
        bool operator==(const KmerInfo& that) const
        {
            return kmer.data == that.kmer.data;
        }
    };
    vector< vector<KmerInfo> > kmerTable1(sequences.size());

    for(uint64_t i=0; i<sequences.size(); i++) {
        const vector<Base>& sequence = sequences[i].first;
        if(false) {
            cout << "Finding unique k-mers for sequence of length " << sequence.size() << endl;
        }
        vector<KmerInfo>& kmerInfos = kmerTable1[i];

        Kmer kmer;
        for(uint64_t position=0; position<kmerLength; position++) {
            kmer.set(position, sequence[position]);
        }

        for(uint64_t position=0; /* Check later */; position++) {
            kmerInfos.push_back({kmer, position});

            if(position + kmerLength == sequence.size()) {
                break;
            }

            // Update the k-mer.
            kmer.shiftLeft();
            kmer.set(kmerLength - 1, sequence[position + kmerLength]);
        }
        SHASTA_ASSERT(kmerInfos.size() == sequence.size() - kmerLength + 1);

        // Only keep the k-mers that appear once.
        if(false) {
            cout << kmerInfos.size() << " total kmers." << endl;
        }
        deduplicateAndCountAndKeepUnique(kmerInfos);
        if(false) {
            cout << kmerInfos.size() << " unique kmers." << endl;
        }
    }



    // Create a global table of unique k-mers in all the sequences.
    class KmerData {
    public:
        Kmer kmer;
        uint64_t sequenceIndex;
        uint64_t position;
        bool operator<(const KmerData& that) const
        {
            return tie(kmer.data, sequenceIndex) < tie(that.kmer.data, that.sequenceIndex);
        }
    };
    vector<KmerData> kmerTable2;
    for(uint64_t sequenceIndex=0; sequenceIndex<sequences.size(); sequenceIndex++) {
        const vector<KmerInfo>& kmerInfos = kmerTable1[sequenceIndex];
        for(const KmerInfo& kmerInfo: kmerInfos) {
            kmerTable2.push_back({kmerInfo.kmer, sequenceIndex, kmerInfo.position});
        }
    }
    sort(kmerTable2.begin(), kmerTable2.end());



    // Now construct a third table that for each unique k-mer
    // gives the sequence indexes and positions the k-mer appears in.
    class UniqueKmerInfo {
    public:
        Kmer kmer;
        uint64_t totalWeight = 0;
        uint64_t minDistanceFromEnds = invalid<uint64_t>;
        class Occurrence {
        public:
            uint64_t sequenceIndex;
            uint64_t position;
        };
        vector<Occurrence> occurrences;
        bool operator<(const UniqueKmerInfo& that) const
        {
            return tie(totalWeight, minDistanceFromEnds) > tie(that.totalWeight, that.minDistanceFromEnds);
        }
        void write(ostream& s, uint64_t kmerLength) const
        {
            kmer.write(s, kmerLength);
            s << " " << totalWeight;
            s << " " << minDistanceFromEnds;
            for(const auto& occurrence: occurrences) {
                s << " (" << occurrence.sequenceIndex << "," <<
                    occurrence.position << ")";
            }
            s << endl;
        }
    };
    vector<UniqueKmerInfo> kmerTable3;
    for(auto it=kmerTable2.begin(); it!= kmerTable2.end(); /* Increment later */) {
        const Kmer kmer = it->kmer;

        // Find the end of the streak for the same kmer.
        auto jt = it;
        while(true) {
            if(jt == kmerTable2.end()) {
                break;
            }
            if(jt->kmer != kmer) {
                break;
            }
            ++jt;
        }

        // Store this streak in kmerTable3.
        UniqueKmerInfo uniqueKmerInfo;
        uniqueKmerInfo.kmer = kmer;
        for(; it!=jt; it++) {
            const uint64_t sequenceIndex = it->sequenceIndex;
            const uint64_t sequenceLength = sequences[sequenceIndex].first.size();
            const uint64_t position = it->position;
            const uint64_t distanceFromLeft = position;
            const uint64_t distanceFromRight = sequenceLength - kmerLength - position;
            const uint64_t distanceFromEnds = min(distanceFromLeft, distanceFromRight);
            uniqueKmerInfo.occurrences.push_back({sequenceIndex, it->position});
            uniqueKmerInfo.totalWeight += sequences[it->sequenceIndex].second;
            uniqueKmerInfo.minDistanceFromEnds = min(uniqueKmerInfo.minDistanceFromEnds, distanceFromEnds);
        }
        kmerTable3.push_back(uniqueKmerInfo);
    }
    sort(kmerTable3.begin(), kmerTable3.end());


    if(false) {
        for(const auto& uniqueKmerInfo: kmerTable3) {
            uniqueKmerInfo.write(cout, kmerLength);
        }
    }

    // The first entry in kmerTable3 gives the optimal splitting kmer,
    // the sequences involves (all of them, in most cases),
    // and the position of the splitting k-mer in each of the sequences.
    SHASTA_ASSERT(not kmerTable3.empty());
    const UniqueKmerInfo& optimalSplitting = kmerTable3.front();
    if(debug) {
        cout << "Splitting at ";
        optimalSplitting.write(cout, kmerLength);
    }


    // Prepare the sequences for the left and right MSA.
    vector< pair<vector<Base>, uint64_t> > leftSequences;
    vector< pair<vector<Base>, uint64_t> > rightSequences;
    vector<Base> leftConsensus;
    vector<Base> rightConsensus;
    for(const auto& occurrence: optimalSplitting.occurrences) {
        const uint64_t sequenceIndex = occurrence.sequenceIndex;
        const auto& p = sequences[sequenceIndex];
        const vector<Base>& sequence = p.first;
        const uint64_t weight = p.second;
        const uint64_t position = occurrence.position;
        leftSequences.push_back(make_pair(vector<Base>(), weight));
        rightSequences.push_back(make_pair(vector<Base>(), weight));
        vector<Base>& leftSequence = leftSequences.back().first;
        vector<Base>& rightSequence = rightSequences.back().first;
        copy(sequence.begin(), sequence.begin() + position,
            back_inserter(leftSequence));
        copy(sequence.begin() + position + kmerLength, sequence.end(),
            back_inserter(rightSequence));
    }

    // Recursive call to do the left and right MSA.
    globalMsa(leftSequences , maxSpoaLength, kmerLength, leftConsensus);
    globalMsa(rightSequences, maxSpoaLength, kmerLength, rightConsensus);

    // Now stitch the pieces together.
    consensus = leftConsensus;
    for(uint64_t position=0; position<kmerLength; position++) {
        consensus.push_back(optimalSplitting.kmer[position]);
    }
    copy(rightConsensus.begin(), rightConsensus.end(),
        back_inserter(consensus));
}



// This just uses spoa.
// It cannot be used for very long sequences due to quadratic
// memory and time. Practical limit is a few thousand bases.
void shasta::globalMsaSpoa(
    const vector< pair<vector<Base>, uint64_t> >& sequences,
    vector<Base>& consensus
    )
{
    // Sanity check.
    SHASTA_ASSERT(not sequences.empty());

    // Trivial case.
    if(sequences.size() == 1) {
        consensus = sequences.front().first;
        return;
    }

    // We want to enter the sequences in order of decreasing weight.
    // Create a table of pairs (sequenceIndex, weight)
    // where sequenceIndex is the index in the sequences vector.
    // Then sort by decreasing weight.
    vector< pair<uint64_t, uint64_t> > sequencesTable;
    for(uint64_t sequenceIndex=0; sequenceIndex<sequences.size(); sequenceIndex++) {
        const auto& p = sequences[sequenceIndex];
        const uint64_t weight = p.second;
        sequencesTable.push_back(make_pair(sequenceIndex, weight));
    }
    sort(sequencesTable.begin(), sequencesTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Create the spoa alignment engine and alignment graph.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto spoaAlignmentEngine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);
    spoa::Graph spoaAlignmentGraph;

    // Add the sequences to the MSA in order of decreasing weight.
    string sequenceString;
    for(uint64_t indexByWeight=0; indexByWeight<sequencesTable.size(); indexByWeight++) {
        const auto& p = sequencesTable[indexByWeight];
        const uint64_t sequenceIndex = p.first;
        const uint64_t weight = p.second;
        const auto& q = sequences[sequenceIndex];
        SHASTA_ASSERT(q.second == weight);
        const vector<Base>& sequence = q.first;

        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString += base.character();
        }
        auto alignment = spoaAlignmentEngine->Align(sequenceString, spoaAlignmentGraph);
        spoaAlignmentGraph.AddAlignment(alignment, sequenceString, uint32_t(weight));
    }

    // Get the MSA alignment.
    // The true argument causes a final alignment entry equal to the consensus.
    vector<string> alignment = spoaAlignmentGraph.GenerateMultipleSequenceAlignment(false);
    SHASTA_ASSERT(alignment.size() == sequencesTable.size());

    // Compute coverage at each alignment position for each of the 5 AlignedBases.
    const uint64_t alignmentLength = alignment.front().size();
    vector< array<uint64_t, 5> > coverage(alignmentLength, {0, 0, 0, 0, 0});
    for(uint64_t indexByWeight=0; indexByWeight<sequencesTable.size(); indexByWeight++) {
        const string& alignmentRow = alignment[indexByWeight];
        SHASTA_ASSERT(alignmentRow.size() == alignmentLength);
        for(uint64_t position=0; position<alignmentLength; position++) {
            const AlignedBase b = AlignedBase::fromCharacter(alignmentRow[position]);
            coverage[position][b.value] += sequencesTable[indexByWeight].second;
        }
    }

    // Compute coverage-based consensus at each alignment position.
    vector<AlignedBase> alignedConsensus;
    for(const auto& c: coverage) {
        const uint64_t iBase = std::max_element(c.begin(), c.end()) - c.begin();
        alignedConsensus.push_back(AlignedBase::fromInteger(iBase));
    }
    SHASTA_ASSERT(alignedConsensus.size() == alignmentLength);

    // Take out the gaps.
    consensus.clear();
    for(const AlignedBase b: alignedConsensus) {
        if(not b.isGap()) {
            consensus.push_back(Base(b));
        }
    }
}



// This just uses spoa.
// It cannot be used for very long sequences due to quadratic
// memory and time. Practical limit is a few thousand bases.
// Version that returns the alignment.
// THE SEQUENCES MUST BE PASSED IN ORDER OF DECREASING WEIGHT.
void shasta::globalMsaSpoa(
    const vector< pair<vector<Base>, uint64_t> >& sequences,
    vector< vector<AlignedBase> >& alignmentArgument
    )
{
    // Sanity check.
    SHASTA_ASSERT(not sequences.empty());

    // Check that the sequences are ordered by decreasing weight.
    for(uint64_t i=1; i<sequences.size(); i++) {
        SHASTA_ASSERT(sequences[i-1].second >= sequences[i].second);
    }

    // Create the spoa alignment engine and alignment graph.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto spoaAlignmentEngine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);
    spoa::Graph spoaAlignmentGraph;

    // Add the sequences to the MSA in order of decreasing weight.
    string sequenceString;
    for(uint64_t i=0; i<sequences.size(); i++) {
        const auto& p = sequences[i];
        const vector<Base>& sequence = p.first;
        const uint64_t weight = p.second;

        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString += base.character();
        }
        auto alignment = spoaAlignmentEngine->Align(sequenceString, spoaAlignmentGraph);
        spoaAlignmentGraph.AddAlignment(alignment, sequenceString, uint32_t(weight));
    }

    // Get the MSA alignment.
    // The true argument causes a final alignment entry equal to the consensus.
    vector<string> alignment = spoaAlignmentGraph.GenerateMultipleSequenceAlignment(false);
    SHASTA_ASSERT(alignment.size() == sequences.size());

    // Copy it to alignmentArgument.
    alignmentArgument.clear();
    alignmentArgument.resize(alignment.size());
    for(uint64_t i=0 ; i<alignment.size(); i++) {
        const string& alignmentRow = alignment[i];
        vector<AlignedBase>& alignmentArgumentRow = alignmentArgument[i];
        alignmentArgumentRow.resize(alignmentRow.size());
        for(uint64_t j=0; j<alignmentRow.size(); j++) {
            alignmentArgumentRow[j] = AlignedBase::fromCharacter(alignmentRow[j]);
        }
    }

}



// Version that enforces a maximum MSA length and returns false if it is exceeded.
bool shasta::globalMsaSpoa(
    const vector< pair<vector<Base>, uint64_t> >& sequences,
    vector<Base>& consensus,
    uint64_t maximumMsaLength
    )
{
    if(sequences.size() > 1) {
        uint64_t maxLength = 0;
        for(const auto& sequence: sequences) {
            maxLength = max(maxLength, sequence.first.size());
        }
        if(maxLength > maximumMsaLength) {
            return false;
        }
    }

    // If getting here, the MSA is no longer than the specified maximum length
    // (or it is trivial, consisting of just one sequence).
    globalMsaSpoa(sequences, consensus);
    return true;
}



// Python-callable version.
std::string shasta::globalMsaPython(
    const vector< pair<string, uint64_t> >& sequenceStrings,
    uint64_t maxSpoaLength,
    uint64_t kmerLength)
{
    // Extract the sequences.
    vector< pair<vector<Base>, uint64_t> > sequences;
    sequences.reserve(sequenceStrings.size());
    for(const auto& p: sequenceStrings) {
        sequences.resize(sequences.size() + 1);
        sequences.back().second = p.second;
        const string& sequenceString = p.first;
        vector<Base>& sequence = sequences.back().first;
        for(const char c: sequenceString) {
            sequence.push_back(Base::fromCharacter(c));
        }
    }

    // Do the MSA.
    vector<Base> consensus;
    globalMsa(sequences, maxSpoaLength, kmerLength, consensus);

    // Construct the consensus string and return it.
    string consensusString;
    for(const Base b: consensus) {
        consensusString.push_back(b.character());
    }
    return consensusString;
}


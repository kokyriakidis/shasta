#pragma once

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>

// Standard libraries.
#include "array.hpp"
#include <map>
#include <cmath>
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class PhasingComponent;
        class PhasingTable;
        class PhasingTableEntry;

        class Anchors;
        class BubbleChain;
    }
}



// A PhasingTableEntry describes the appearances of one oriented read
// on one or both sides of a diploid Bubble of a BubbleChain.
// The frequency array contains the number of times the oriented read
// appears on non-terminal marker graph edges of the two Chains of the diploid Bubble.
class shasta::mode3::PhasingTableEntry {
public:

    PhasingTableEntry(
        OrientedReadId orientedReadId,
        uint64_t positionInBubbleChain) :
        orientedReadId(orientedReadId),
        positionInBubbleChain(positionInBubbleChain)
        {}

    // The OrientedReadId this PhasingTableEntry refers to,
    // and its index in the PhasingTable::orientedReads vector.
    OrientedReadId orientedReadId;
    uint64_t orientedReadIndex = invalid<uint64_t>;

    // The position in the bubble chain of the diploid bubble
    // this PhasingTableEntry refers to,
    // and its index in the PhasingTable::orientedReads vector.
    uint64_t positionInBubbleChain;
    uint64_t bubbleIndex = invalid<uint64_t>;

    // The number of times this oriented read
    // appears on non-terminal marker graph edges of the two Chains of the diploid Bubble.
    // The two entries in the array corresponds to the two chains of the diploid Bubble.
    array<uint64_t, 2> frequency = {0, 0};

    // The phase of this oriented read relative to this bubble
    // is computed from the frequency array.

    // The relative phase varies continuously between -1 and 1 and is:
    // * +1 if this oriented read always appears in Chain 0 (that is, frequency[1] is 0).
    // * -1 if this oriented read always appears in Chain 1 (that is, frequency[0] is 0).
    // * 0 if this oriented appears with equal frequency on Chain 0 and Chain 1
    //   (that is, frequency[0] = frequency[1]).
    double relativePhase = invalid<double>;

    // The discrete relative phase can be:
    // +1 if relativePhase > +1. - phaseErrorThreshold.
    // -1 if relativePhase < -1. + phaseErrorThreshold.
    // 0 otherwise.
    int64_t discreteRelativePhase = invalid<int64_t>;

    // Compute and store the relativePhase and discreteRelativePhase.
    void storeRelativePhase(double phaseErrorThreshold)
    {
        relativePhase = 2. * double(frequency[0]) / double(frequency[0] + frequency[1]) - 1.;
        if(relativePhase > 1. - phaseErrorThreshold) {
            discreteRelativePhase = +1;
        } else if(relativePhase < -1. + phaseErrorThreshold) {
            discreteRelativePhase = -1;
        } else {
            discreteRelativePhase = 0;
        }
    }

    void writeCsv(ostream&) const;
};



// A PhasingTable is a set of PhasingTableEntry objects,
// randomly accessible by orientedReadId and by positionInBubbleChain.
class shasta::mode3::PhasingTable: public boost::multi_index_container<PhasingTableEntry,
    boost::multi_index::indexed_by <

        // Index by (orientedReadId, positionInBubbleChain) (unique).
        boost::multi_index::ordered_unique<
            boost::multi_index::composite_key<
                PhasingTableEntry,
                boost::multi_index::member<PhasingTableEntry, OrientedReadId ,&PhasingTableEntry::orientedReadId>,
                boost::multi_index::member<PhasingTableEntry, uint64_t, &PhasingTableEntry::positionInBubbleChain>
                > >,

        // Index by orientedReadId (non-unique).
        boost::multi_index::ordered_non_unique<boost::multi_index::member<
                PhasingTableEntry,
                OrientedReadId,
                &PhasingTableEntry::orientedReadId> >,

        // Index by positionInBubbleChain (non-unique).
        boost::multi_index::ordered_non_unique<boost::multi_index::member<
            PhasingTableEntry,
            uint64_t,
            &PhasingTableEntry::positionInBubbleChain> >
    > > {
public:

    PhasingTable(
        const BubbleChain&,
        const Anchors&,
        double phaseErrorThreshold);

    uint64_t entryCount() const
    {
        return size();
    }
    uint64_t unambiguousEntryCount() const;
    uint64_t ambiguousEntryCount() const;

    uint64_t bubbleCount() const
    {
        return bubbles.size();
    }

    uint64_t orientedReadCount() const
    {
        return orientedReads.size();
    }

    // Experimental. Do not use.
    void simpleIterativePhasing1();
    void simpleIterativePhasing2();

    // Optimize the phases of the oriented reads and of the bubbles.
    void greedyPhasing();

    void writeCsv(const string& fileNamePrefix) const;
    enum class ColoringMethod {
        byRelativePhase,
        byDiscreteRelativePhase,
        byConsistency
    };
    void writePng(const string& fileName, ColoringMethod) const;

    double bubbleErrorRate(uint64_t positionInBubbleChain) const;

    vector< shared_ptr<PhasedComponent> > phasedComponents;
    void constructPhasedComponents(bool debug);

private:
    const auto& indexByBoth() const {return get<0>();}
    const auto& indexByOrientedReadId() const {return get<1>();}
    const auto& indexByPositionInBubbleChain() const {return get<2>();}

    void fill(
        const BubbleChain&,
        const Anchors&,
        double phaseErrorThreshold);



    // Information about the orientedReads that appears in the PhasingTable.
    class OrientedRead {
    public:
        OrientedReadId id;
        uint64_t minPositionInBubbleChain;
        uint64_t maxPositionInBubbleChain;
        int64_t phase = 0;  // -1, 0 or +1
        uint64_t phasingComponent = invalid<uint64_t>;
    };
    void gatherOrientedReads();
    vector<OrientedRead> orientedReads;

    // Map OrientedReadId to an index in the orientedReadInfos vector.
    std::map<OrientedReadId, uint64_t> orientedReadsMap;



    // Information about the diploid bubbles in this PhasingTable.
    class Bubble {
    public:
        uint64_t positionInBubbleChain;
        int64_t phase = 0;  // -1, 0 or +1
        uint64_t phasingComponent = invalid<uint64_t>;
    };
    vector<Bubble> bubbles;
    void gatherBubbles();

    // Map a positionInBubbleChain to an index in the bubbles vector.
public:
    std::map<uint64_t, uint64_t> bubblesMap;
private:



    // Fill the orientedReadIndex and bubbleIndex in all PhasingTableEntries.
    // This can only be done after gatherOrientedReads and gatherBubbles
    // have been called.
    void fillIndexes();

    // Compute the consistency state of a PhasingTableEntry relative
    // to the current phases of its oriented read and bubble.
    // It can be +1 (consistent), -1 (inconsistent), or 0 (unassigned or ambiguous).
    // See the implementation for details.
    int64_t consistencyState(const PhasingTableEntry&) const;

    // Count the number of (consistent,inconsistent) PhasingTableEntries
    // for an oriented read based on the phases currently assigned
    // to bubbles and oriented reads.
    pair<uint64_t, uint64_t> countConsistentEntriesForOrientedRead(OrientedReadId) const;

    // Count the number of (consistent,inconsistent) PhasingTableEntries
    // for the bubble at a given bubble chain position based on the phases currently assigned
    // to bubbles and oriented reads.
    pair<uint64_t, uint64_t> countConsistentEntriesForBubble(uint64_t positionInBubbleChain) const;

    // Count the number of (unambiguous, ambiguous) PhasingTableEntries
    // for the bubble at a given bubble chain position based on the phases currently assigned
    // to bubbles and oriented reads.
    pair<uint64_t, uint64_t> countEntriesForBubble(uint64_t positionInBubbleChain) const;

public:
    // Count the number of (consistent,inconsistent) PhasingTableEntries
    // based on the phases currently assigned
    // to bubbles and oriented reads.
    pair<uint64_t, uint64_t> countConsistentEntries() const;

private:
    void writeOrientedReadsCsv(const string& fileNamePrefix) const;
    void writeBubblesCsv(const string& fileNamePrefix, bool writePhasingInformation) const;
    void writeDetailsCsv(const string& fileNamePrefix) const;
};



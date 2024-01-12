#ifndef SHASTA_MODE3B_PHASING_TABLE_HPP
#define SHASTA_MODE3B_PHASING_TABLE_HPP

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
    namespace mode3b {
        class PhasingTable;
        class PhasingTableEntry;

        class BubbleChain;
    }
    class MarkerGraph;
}



// A PhasingTableEntry describes the appearances of one oriented read
// on one or both sides of a diploid Bubble of a BubbleChain.
// The frequency array contains the number of times the oriented read
// appears on non-terminal marker graph edges of the two Chains of the diploid Bubble.
class shasta::mode3b::PhasingTableEntry {
public:
    OrientedReadId orientedReadId;

    // The position in the bubble chain of the diploid bubble
    // this PhasingTableEntry refers to.
    uint64_t positionInBubbleChain;

    uint64_t orientedReadIndex = invalid<uint64_t>;
    uint64_t bubbleIndex = invalid<uint64_t>;

    // The number of times this oriented read
    // appears on non-terminal marker graph edges of the two Chains of the diploid Bubble.
    // The two entries in the array corresponds to the two chains of the diploid Bubble.
    array<uint64_t, 2> frequency = {0, 0};

    // The phase is:
    // * +1 if this oriented read always appears in Chain 0 (that is, frequency[1] is 0).
    // * -1 if this oriented read always appears in Chain 1 (that is, frequency[0] is 0).
    // * 0 if this oriented appears with equal frequency on Chain 0 and Chain 1
    //   (that is, frequency[0] = frequency[1]).
    // This does not take into accoubnt possible flipping of the Bubble
    // as stored in the PhasingTable.
    double phase() const
    {
        return 2. * double(frequency[0]) / double(frequency[0] + frequency[1]) - 1.;
    }

    pair<OrientedReadId, uint64_t> key() const
    {
        return *reinterpret_cast< const pair<OrientedReadId, uint64_t>* >(&orientedReadId);
    }

    void writeCsv(ostream&) const;
};



// A PhasingTable is a set of PhasingTableEntry objects,
// randomly accessible by orientedReadId and by positionInBubbleChain.
class shasta::mode3b::PhasingTable: public boost::multi_index_container<PhasingTableEntry,
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
        const MarkerGraph&,
        double phaseError);

    void write(const string& fileNamePrefix) const;
    void writePng(const string& fileName, bool colorByType) const;

    uint64_t entryCount() const
    {
        return size();
    }
    uint64_t unambiguousEntryCount() const;

    uint64_t bubbleCount() const
    {
        return bubbles.size();
    }

    uint64_t orientedReadCount() const
    {
        return orientedReadInfos.size();
    }

    void flipSweep();
    uint64_t discordantCount() const;

private:
    double phaseError;
    double phaseThresholdPlus;
    double phaseThresholdMinus;

    void fill(
        const BubbleChain&,
        const MarkerGraph&);

    class Bubble {
    public:
        uint64_t positionInBubbleChain;
        bool flip = false;  // 0 and 1 sides should be swapped.
    };
    vector<Bubble> bubbles;
    void gatherBubbles();

    // The phase of a PhasingTableEntry, taking into account possible flipping
    // of the bubble.
    // The phase is in [-1., 1.].
    double phase(const PhasingTableEntry& phasingTableEntry) const
    {
        double p = phasingTableEntry.phase();
        if(bubbles[phasingTableEntry.bubbleIndex].flip) {
            p = -p;
        }
        return p;
    }

    // The hue corresponding to a phase. This is:
    // 360 (red) if phase is +1 (that is, all frequency is on Chain 0).
    // 300 (magenta) if phase is 0.
    // 240 (blue) if phase is -1 (that is, all frequency is on Chain 1).
    static uint64_t hue(double phase)
    {
        return uint64_t(std::round(300. + 60. * phase));
    }

    // Map a positionInBubbleChain to an index in the bubbles vector.
    std::map<uint64_t, uint64_t> bubblesMap;

    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        uint64_t minPositionInBubbleChain;
        uint64_t maxPositionInBubbleChain;
    };
    void gatherOrientedReads();
    vector<OrientedReadInfo> orientedReadInfos;

    // Count the number of PhaseTableEntries for a given oriented read
    // with phase +1 or -1, allowing a phaseError discrepancy up to phase Error.
    // That is, count0 is the number of oriented read entries with phase >= 1 - phaseError,
    // and count1 is the number of oriented read entries with phase <= -1 + phaseError.
    void count(
        OrientedReadId,
        uint64_t& countPlus,
        uint64_t& countMinus) const;

    void fillIndexes();

    // Map OrientedReadId to an index in the orientedReadInfos vector.
    std::map<OrientedReadId, uint64_t> orientedReadIdsMap;

    void writeCsv(const string& fileNamePrefix) const;
    void writeOrientedReadsCsv(const string& fileNamePrefix) const;
    void writeBubblesCsv(const string& fileNamePrefix) const;
    void writeDetailsCsv(const string& fileNamePrefix) const;
    void writeHtml(const string& fileNamePrefix) const;
};

#endif


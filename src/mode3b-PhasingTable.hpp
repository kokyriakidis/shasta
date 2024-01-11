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
        const MarkerGraph&);

    void writeCsv(const string& fileNamePrefix) const;

    uint64_t bubbleCount() const
    {
        return bubbles.size();
    }

    uint64_t orientedReadCount() const
    {
        return orientedReadInfos.size();
    }

private:
    void fill(
        const BubbleChain&,
        const MarkerGraph&);

    class Bubble {
    public:
        uint64_t positionInBubbleChain;
    };
    vector<Bubble> bubbles;
    void gatherBubbles();

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

    void fillIndexes();

    // Map OrientedReadId to an index in the orientedReadInfos vector.
    std::map<OrientedReadId, uint64_t> orientedReadIdsMap;


    void writeOrientedReadsCsv(const string& fileNamePrefix) const;
    void writeBubblesCsv(const string& fileNamePrefix) const;
    void writeDetailsCsv(const string& fileNamePrefix) const;
};

#endif


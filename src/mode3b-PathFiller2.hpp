#ifndef SHASTA_MODE3B_PATH_FILLER2_HPP
#define SHASTA_MODE3B_PATH_FILLER2_HPP

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.

// Standard library.
#include "iosfwd.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3b {
        class PathFiller2;
        class PathFiller2DisplayOptions;
    }
    class Assembler;
};



class shasta::mode3b::PathFiller2DisplayOptions {
public:

    // If this is not open, no output takes place.
    ostream& html;

    PathFiller2DisplayOptions(ostream& html) : html(html) {}
};



class shasta::mode3b::PathFiller2 {
public:

    // Hide class Base defined in boost::adjacency_list.
    using Base = shasta::Base;

    PathFiller2(
        const Assembler&,
        MarkerGraphEdgeId edgeIdA,
        MarkerGraphEdgeId edgeIdB,
        const PathFiller2DisplayOptions&);

private:

    // Store constructor arguments.
    const Assembler& assembler;
    MarkerGraphEdgeId edgeIdA;
    MarkerGraphEdgeId edgeIdB;
    const PathFiller2DisplayOptions& options;
    ostream& html;

    void checkAssumptions() const;



    // A class used to store an ordinal and the corresponding position
    // of a marker in an oriented read.
    class OrdinalAndPosition {
    public:
        uint32_t ordinal = invalid<uint32_t>;
        uint32_t position  = invalid<uint32_t>;
        bool isValid() const
        {
            return ordinal != invalid<uint32_t>;
        }
    };



    // Information about the portion of an oriented read used in this assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        OrientedReadInfo(OrientedReadId orientedReadId) :
            orientedReadId(orientedReadId)
            {}

        // The position of the source and target vertex of edgeIdA in this oriented read.
        // If this oriented read does not appear in edgeIdA, these are left uninitialized.
        OrdinalAndPosition ordinalAndPositionA0;
        OrdinalAndPosition ordinalAndPositionA1;
        bool isOnA() const
        {
            return ordinalAndPositionA0.isValid();
        }

        // The position of the source and target vertex of edgeIdB in this oriented read.
        // If this oriented read does not appear in edgeIdB, this is left uninitialized.
        OrdinalAndPosition ordinalAndPositionB0;
        OrdinalAndPosition ordinalAndPositionB1;
        bool isOnB() const
        {
            return ordinalAndPositionB0.isValid();
        }

        // Order by OrientedReadId.
        bool operator<(const OrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }

        int32_t ordinalOffset() const
        {
            SHASTA_ASSERT(isOnA() and isOnB());
            const int32_t offset0 =
                int32_t(ordinalAndPositionB0.ordinal) -
                int32_t(ordinalAndPositionA0.ordinal);
            const int32_t offset1 =
                int32_t(ordinalAndPositionB1.ordinal) -
                int32_t(ordinalAndPositionA1.ordinal);
            SHASTA_ASSERT(offset0 == offset1);
            return offset0;
        }

        int32_t positionOffset0() const
        {
            SHASTA_ASSERT(isOnA() and isOnB());
            return
                int32_t(ordinalAndPositionB0.position) -
                int32_t(ordinalAndPositionA0.position);
        }

        int32_t positionOffset1() const
        {
            SHASTA_ASSERT(isOnA() and isOnB());
            return
                int32_t(ordinalAndPositionB1.position) -
                int32_t(ordinalAndPositionA1.position);
        }
    };



    // For assembly, we use the union of the oriented reads
    // that appear in edgeIdA and edgeIdB.
    // In contrast, class PathFiller1 used the intersection of those oriented reads,
    // which resulted in low coverage and low accuracy when edgeIdA and edgeIdB are
    // very distant from each other, and there are not many oriented reads that
    // cover both of them.
    // OrientedReadInfos are stored sorted by OrientedReadId.
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    void writeOrientedReads() const;

    // The index of an OrientedReadId is its index in the orientedReadInfos vector.
    uint64_t getOrientedReadIndex(OrientedReadId) const;
};

#endif

/********************************************************************************

Creation of Anchors from json input.

This expects on input a json file containing a list of candidate anchors,
where each candidate anchor is a list of base intervals, and each base interval
is of the form:

["readName", strand, begin, end]
where:
- readName identifies the read using its name as in the input fasta or fastq.
- strand can be 1 or 0 to indicate whether the read is reverse complemented or not.
- begin is the first base position in the read that belongs to the anchor.
- end if the first base position in the read after the anchor.
If strand is "-", begin and end refer to base positions after reverse complement\ing,
so end is always > begin.

All the anchor read sequences in each candidate anchor must be exactly identical.

Candidate anchors are clipped to the first and last markers entirely contained in the candidate anchor.
If no markers are entirely contained in the candidate anchor, the candidate anchor is discarded.
Otherwise, the clipped candidate anchor is used to generate a pair of referce complemented anchors.

So in summary the input looks something like this (spacing only used for readability):
[
    [
        ["read1", "+". 400, 450],
        ["read2", "-", 2000, 2050],
        ["read3", "+", 6100 ,6150]
    [
        ["read4", "-". 12300, 12400],
        ["read5", "-", 4300, 4400],
        ["read6", "+", 8470 ,8570]
    ]
]

********************************************************************************/

// Shasta.
#include "mode3-Anchor.hpp"
#include "Base.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// Standatd library.
#include "fstream.hpp"



Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const string& jsonFileName,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t /* threadCount */) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    // Open the input json file.
    ifstream inputFile(jsonFileName);
    if(not inputFile) {
        throw runtime_error("Could not open " + jsonFileName);
    }

    // Read the json file into a boost property tree.
    Ptree json;
    boost::property_tree::read_json(inputFile, json);



    // Loop over the candidate anchors.
    uint64_t successCount = 0;
    for(const auto& p: json) {
        const Ptree& candidateAnchor = p.second;

        // If not in the desired coverage range, skip it.
        const uint64_t coverage = candidateAnchor.size();
        if(
            (coverage < minPrimaryCoverage) or
            (coverage > maxPrimaryCoverage)) {
            continue;
        }

        if(processCandidateAnchor(candidateAnchor)) {
            ++successCount;
        } else {
            cout << "Discarded." << endl;
        }
    }
    cout << "Of " << json.size() << " candidate anchors, " << successCount <<
        " were used and " << json.size() - successCount << " were discarded." << endl;


    throw runtime_error("Anchor creation from json is not yet implemented.");
}


// Process a candidate anchor from json input.
bool Anchors::processCandidateAnchor(
    const Ptree& candidateAnchor)
{
    using Ptree = boost::property_tree::basic_ptree<string, string>;

    class Interval {
    public:
        OrientedReadId orientedReadId;
        uint64_t begin;
        uint64_t end;
        uint64_t clippedBegin;
        uint64_t clippedEnd;
        uint32_t ordinal0;
        uint32_t ordinal1;
        uint64_t length() const
        {
            return end - begin;
        }
        uint64_t leftClip() const
        {
            return clippedBegin - begin;
        }
        uint64_t rightClip() const
        {
            return end - clippedEnd;
        }
    };

    // Gather the intervals of this candidate anchor.
    vector<Interval> intervals;
    for(const auto& p: candidateAnchor) {
        const Ptree& candidateAnchorInterval = p.second;

        // The interval must have 4 entries (read name strand, begin, end).
        if(candidateAnchorInterval.size() != 4) {
            std::ostringstream s;
            boost::property_tree::write_json(s, candidateAnchorInterval);
            throw runtime_error("Anchor interval has size " + to_string(candidateAnchorInterval.size()) +
                ", must be 4: " + s.str());
        }

        // Parse the interval.
        try {
            // Get the read name.
            auto it = candidateAnchorInterval.begin();
            const string readName = ((it++)->second).get<string>("");

            // Get the ReadId.
            const ReadId readId = reads.getReadId(readName);
            if(readId == invalidReadId) {
                cout << "Read " << readName << " does not exist." << endl;
                throw runtime_error("Read does not exist.");
            }

            // Get the Strand
            const Strand strand = ((it++)->second).get<Strand>("");
            if((strand != 0) and (strand != 1)) {
                cout << "Invalid strand." << endl;
                throw runtime_error("Invalid strand.");
            }
            const OrientedReadId orientedReadId(readId, strand);

            // Get the begin, end of the interval.
            const uint64_t readLength = reads.getReadRawSequenceLength(readId);
            const uint64_t begin = ((it++)->second).get<uint64_t>("");
            const uint64_t end = ((it++)->second).get<uint64_t>("");
            if((begin >= readLength) or (end > readLength)) {
                cout << "Invalid begin/end. Read length is " << readLength << endl;
                throw runtime_error("Invalid begin/end.");
            }
            intervals.push_back({orientedReadId, begin, end});

        } catch (...) {
            std::ostringstream s;
            boost::property_tree::write_json(s, candidateAnchorInterval);
            throw runtime_error("Invalid anchor interval : " + s.str());
        }
    }


    // Check that the sequences on these intervals are identical.
    const Interval& interval0 = intervals[0];
    const uint64_t length0 = interval0.length();
    if(length0 < k) {
        // cout << "Anchor is too short." << endl;
        return false;
    }
    for(uint64_t i=1; i<intervals.size(); i++) {
        const Interval& interval = intervals[i];

        // Check the lengths.
        const uint64_t length = interval0.length();
        if(length != length0) {
            std::ostringstream s;
            boost::property_tree::write_json(s, candidateAnchor);
            throw runtime_error("Interval lengths must all be identical: " + s.str());
        }

        // Check the bases.
        for(uint64_t j=0; j<length; j++) {
            const Base b0 = reads.getOrientedReadBase(interval0.orientedReadId, uint32_t(interval0.begin + j));
            const Base b = reads.getOrientedReadBase(interval.orientedReadId, uint32_t(interval.begin + j));
            SHASTA_ASSERT(b == b0);
        }
    }

    // Object to compare CompressedMarkers by position.
    class Comparator {
    public:
        bool operator()(
            const CompressedMarker& x,
            const CompressedMarker& y) const
        {
            return x.position < y.position;
        }
    };
    const Comparator comparator;

    // Clip to the first/last marker entirely contained in each interval.
    for(Interval& interval: intervals) {
        // cout << "Clipping " << interval.orientedReadId << " " << interval.begin << " " << interval.end << endl;
        const span<const CompressedMarker> orientedReadMarkers = markers[interval.orientedReadId.getValue()];

        // Find the ordinal of the first marker entirely contained in this interval.
        auto it = std::lower_bound(
            orientedReadMarkers.begin(), orientedReadMarkers.end(),
            CompressedMarker({Uint24(uint32_t(interval.begin))}), comparator);
        if(it == orientedReadMarkers.end()) {
            // cout << "The interval begins after the last marker." << endl;
            return false;
        }

        const uint32_t ordinal0 = uint32_t(it - orientedReadMarkers.begin());
        const uint32_t position0 = it->position;
        if(position0 + k > interval.end) {
            // cout << "No marker is entirely inside this interval." << endl;
            return false;
        }

        // Ok, the marker at ordinal0 is entirely inside the interval.
        // Find the ordinal of the last marker entirely contained in this interval.
        uint32_t ordinal1 = ordinal0;
        uint32_t position1 = position0;
        for(; ordinal1<orientedReadMarkers.size(); ++ordinal1) {
            const uint32_t newPosition1 = orientedReadMarkers[ordinal1].position;
            if(newPosition1 + k > interval.end) {
                break;
            }
            position1 = newPosition1;
        }
        ordinal1 -= 1;
        position1 += uint32_t(k);
        // cout << "Clipped to ordinals " << ordinal0 << " " << ordinal1 <<
        //     ", positions " << position0 << " " << position1 << endl;

        interval.ordinal0 = ordinal0;
        interval.ordinal1 = ordinal1;
        interval.clippedBegin = position0;
        interval.clippedEnd = position1;
    }

    // cout << "Clipped bases: " << intervals.front().leftClip() << " " << intervals.front().rightClip() << endl;

    // Check that the number of bases clipped is the same for all the intervals.
    const uint64_t leftClip0 = intervals[0].leftClip();
    const uint64_t rightClip0 = intervals[0].rightClip();
    for(uint64_t i=1; i<intervals.size(); i++) {
        const Interval& interval = intervals[i];
        SHASTA_ASSERT(interval.leftClip() == leftClip0);
        SHASTA_ASSERT(interval.rightClip() == rightClip0);
    }

    return true;
}


// Shasta.
#include "Assembler.hpp"
#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
// Shasta.
#include "seqan.hpp"
using namespace shasta;
using namespace mode0;

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include <map>
#include <set>



void Assembler::computePseudoPath(
    OrientedReadId orientedReadId,

    // The marker graph path computed using computeOrientedReadMarkerGraphPath.
    // This is computed by this function - it does not need to be filled in
    // in advance.
    vector<MarkerGraph::EdgeId>& path,
    vector< pair<uint32_t, uint32_t> >& pathOrdinals,

    // The pseudo-path computed by this function.
    PseudoPath& pseudoPath) const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    using SegmentId = AssemblyGraphEdgeId;

    // Compute the marker graph path.
    const uint64_t markerCount = markers.size(orientedReadId.getValue());
    if(markerCount < 2) {
        pathOrdinals.clear();
        path.clear();
    } else {
        computeOrientedReadMarkerGraphPath(
            orientedReadId,
            0, uint32_t(markerCount - 1),
            path, pathOrdinals);
        SHASTA_ASSERT(path.size() == pathOrdinals.size());
    }



    // Now compute the pseudo-path.
    pseudoPath.clear();
    pseudoPath.clear();
    PseudoPathEntry pseudoPathEntry;
    pseudoPathEntry.segmentId = std::numeric_limits<SegmentId>::max();
    for(uint64_t i=0; i<path.size(); i++) {
        const MarkerGraph::EdgeId markerGraphEdgeId = path[i];
        const pair<uint32_t, uint32_t>& ordinals = pathOrdinals[i];

        // Get the corresponding assembly graph segments.
        const span<const pair<SegmentId, uint32_t> > v =
            assemblyGraph.markerToAssemblyTable[markerGraphEdgeId];

        // If no segments, skip.
        if(v.size() == 0) {
            continue;
        }

        // If detangling was used, there can be more than one,
        // and we don't want this here.
        SHASTA_ASSERT(v.size() == 1);

        // There is only one segment.
        const SegmentId segmentId = v.front().first;
        const uint32_t positionInSegment = v.front().second;

        // Same as the previous.
        if(segmentId == pseudoPathEntry.segmentId) {
            pseudoPathEntry.lastOrdinal = ordinals.second;
            pseudoPathEntry.lastPosition = positionInSegment;
            ++pseudoPathEntry.markerGraphEdgeCount;
            continue;
        }

        // This is the next segment edge encountered
        // by this oriented read along its marker graph path.
        if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
            pseudoPath.push_back(pseudoPathEntry);
        }
        pseudoPathEntry.segmentId = segmentId;
        pseudoPathEntry.firstOrdinal = ordinals.first;
        pseudoPathEntry.lastOrdinal = ordinals.second;
        pseudoPathEntry.firstPosition = positionInSegment;
        pseudoPathEntry.lastPosition = positionInSegment;
        pseudoPathEntry.markerGraphEdgeCount = 1;
    }

    // Add the last entry.
    if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
        pseudoPath.push_back(pseudoPathEntry);
    }
}



// Write the pseudo-path of an oriented read to a csv file.
// The pseudo-path is the sequence of assembly graph edges
// (not necsssarily all adjacent, so not necessatily a path)
// encountered by the oriented read.
void Assembler::writePseudoPath(ReadId readId, Strand strand) const
{
    // Compute the pseudo path.
    const OrientedReadId orientedReadId(readId, strand);
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    computePseudoPath(orientedReadId, path, pathOrdinals, pseudoPath);

    // Write it out.
    ofstream csv("PseudoPath.csv");
    csv << "Segment id,First ordinal,Last ordinal,"
        "First position in segment,Last position in segment, Marker graph edge count\n";
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        csv << pseudoPathEntry.segmentId << ",";
        csv << pseudoPathEntry.firstOrdinal << ",";
        csv << pseudoPathEntry.lastOrdinal << ",";
        csv << pseudoPathEntry.firstPosition << ",";
        csv << pseudoPathEntry.lastPosition << ",";
        csv << pseudoPathEntry.markerGraphEdgeCount << "\n";
    }
}



// Get the vector of segments corresponding to a PseudoPath.
void Assembler::getPseudoPathSegments(
    const PseudoPath& pseudoPath,
    vector<AssemblyGraph::EdgeId>& segmentIds)
{
    segmentIds.clear();
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        segmentIds.push_back(pseudoPathEntry.segmentId);
    }
}

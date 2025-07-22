

#include "Assembler.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "diploidBayesianPhase.hpp"
#include "extractKmer.hpp"
#include "performanceLog.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "shastaTypes.hpp"
#include "timestamp.hpp"
#include "orderPairs.hpp"
#include "Mode3Assembler.hpp"
#include <iostream>
#include <vector>
#include "deduplicate.hpp"
#include <unordered_set>
#include <boost/icl/interval_map.hpp> // Include Boost Interval Container Library header
#include <set> // Include set for the interval map value
#include <algorithm>

using namespace shasta;

// Standard library.
#include "fstream.hpp"
#include <queue>
#include <queue>
#include <set>
#include <unordered_map>

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graph_utility.hpp>


namespace shasta {
    class ReadGraph5;
    class ReadGraph5Vertex;
    class ReadGraph5Edge;

    using ReadGraph5BaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ReadGraph5Vertex,
        ReadGraph5Edge>;

    class ReadGraph5AllAlignments;
    class ReadGraph5AllAlignmentsVertex;
    class ReadGraph5AllAlignmentsEdge;

    using ReadGraph5AllAlignmentsBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ReadGraph5AllAlignmentsVertex,
        ReadGraph5AllAlignmentsEdge>;
}



class shasta::ReadGraph5Vertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

};

class shasta::ReadGraph5AllAlignmentsVertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

};


class shasta::ReadGraph5Edge {
public:
    uint64_t alignmentId;
    ReadGraph5Edge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};

class shasta::ReadGraph5AllAlignmentsEdge {
public:
    uint64_t alignmentId;
    ReadGraph5AllAlignmentsEdge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};

// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph5: public ReadGraph5BaseClass {
public:


    ReadGraph5(uint64_t n) : ReadGraph5BaseClass(n) {}
    
    bool findPathWithPositiveOffset(
        OrientedReadId start,
        vector<vector<OrientedReadId>>& paths,
        vector<vector<double>>& pathsOffsets,
        vector<OrientedReadId>& currentPath,
        vector<double>& currentPathOffset,
        std::set<vertex_descriptor>& visited,
        uint64_t maxDistance,
        uint64_t currentDistance,
        MemoryMapped::Vector<AlignmentData>& alignmentData,
        ReadGraph5& readGraph);
    void findNeighborsUndirectedGraph(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideRight(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideLeft(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphBothSides(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsSkipSameComponentNodes(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachSameComponentNode(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};

// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph5AllAlignments: public ReadGraph5AllAlignmentsBaseClass {
public:

    ReadGraph5AllAlignments(uint64_t n) : ReadGraph5AllAlignmentsBaseClass(n) {}
    void computeShortPath(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        uint64_t maxDistance,
        vector<uint64_t>& path,
        vector<uint64_t>& distance,
        vector<OrientedReadId>& reachedVertices,
        vector<uint64_t>& parentEdges,
        MemoryMapped::Vector<AlignmentData>& alignmentData,
        ReadGraph5AllAlignments& readGraph);
    void findNeighborsUndirectedGraph(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideRight(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphOneSideLeft(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsDirectedGraphBothSides(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsSkipSameComponentNodes(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachSameComponentNode(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachEndNode(OrientedReadId orientedReadId, vector<bool>& finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, uint64_t maxDistance, vector<OrientedReadId>& neighbors, vector<uint64_t>& neighborsAlignmentIds, ReadGraph5AllAlignments& readGraphAllAlignments);
    void findNeighborsEarlyStopWhenReachEndNode(OrientedReadId orientedReadId, vector<bool>& finalDeadEndReadsWithNoIncomingNodesPlusDistanceNeighbors, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    uint64_t findAllPathsToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<vector<OrientedReadId>>& paths,
        vector<OrientedReadId>& currentPath,
        std::set<vertex_descriptor>& visited,
        uint64_t maxDistance,
        uint64_t currentDistance);
    uint64_t findSinglePathToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<vector<OrientedReadId>>& paths,
        vector<OrientedReadId>& currentPath,
        std::set<vertex_descriptor>& visited,
        uint64_t maxDistance,
        uint64_t currentDistance);
    uint64_t findShortestPathToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<OrientedReadId>& shortestPath);
    uint64_t findShortestPathToNode(
        OrientedReadId start,
        OrientedReadId endNode,
        vector<OrientedReadId>& shortestPath,
        uint64_t maxDistance);

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};



bool isSiteInHomopolymerRegionReadGraph5(
    uint64_t sitePositionInRead,
    const shasta::LongBaseSequenceView& readSequence,
    uint64_t homopolymerThreshold = 3)
{
    const uint64_t readLength = readSequence.baseCount;
    SHASTA_ASSERT(sitePositionInRead < readLength);

    const shasta::Base siteBase = readSequence[sitePositionInRead];
    const uint8_t siteBaseValue = siteBase.value; // Assuming Base::value exists

    // Determine scan boundaries for characters to the left/right of the site.
    // These are the furthest characters (inclusive) that the scan loops will consider.
    uint64_t leftScanLoopEndPos = 0;
    if (sitePositionInRead - homopolymerThreshold < 0) {
        leftScanLoopEndPos = 0;
    } else {
        leftScanLoopEndPos = sitePositionInRead - homopolymerThreshold;
    }

    uint64_t rightScanLoopEndPos = 0;
    if (sitePositionInRead + homopolymerThreshold >= readLength) {
        rightScanLoopEndPos = readLength - 1;
    } else {
        rightScanLoopEndPos = sitePositionInRead + homopolymerThreshold;
    }


    // --- Forward scan (strictly to the right of sitePositionInRead) ---
    uint8_t forwardRunBaseVal = 0; // Will hold the base value of the forward run
    uint64_t forwardRunLength = 0;
    bool forwardBaseFound = false;

    // Start scanning from the base immediately to the right of the site
    if (sitePositionInRead + 1 < readLength) { // Check if there is any base to the right
        for (uint64_t i = sitePositionInRead + 1; i <= rightScanLoopEndPos; ++i) {    
            const shasta::Base currentForwardBase = readSequence[i];
            if (!forwardBaseFound) {
                forwardRunBaseVal = currentForwardBase.value;
                forwardRunLength = 1;
                forwardBaseFound = true;
            } else {
                if (currentForwardBase.value == forwardRunBaseVal) {
                    forwardRunLength++;
                } else {
                    break; // Homopolymer run ended
                }
            }
        }
    }

    // --- Backward scan (strictly to the left of sitePositionInRead) ---
    uint8_t backwardRunBaseVal = 0; // Will hold the base value of the backward run
    uint64_t backwardRunLength = 0;
    bool backwardBaseFound = false;

    // Start scanning from the base immediately to the left of the site
    if (sitePositionInRead > 0) { // Check if there is any base to the left
        for (int64_t i = static_cast<int64_t>(sitePositionInRead) - 1; i >= leftScanLoopEndPos; --i) {
            const shasta::Base currentBackwardBase = readSequence[static_cast<uint64_t>(i)];
            if (!backwardBaseFound) {
                backwardRunBaseVal = currentBackwardBase.value;
                backwardRunLength = 1;
                backwardBaseFound = true;
            } else {
                if (currentBackwardBase.value == backwardRunBaseVal) {
                    backwardRunLength++;
                } else {
                    break; // Homopolymer run ended
                }
            }
        }
    }

    // --- Combine results and check conditions ---
    // 'effectiveForwardLength' and 'effectiveBackwardLength' will store the length of runs
    // including the site base, if the site base extends them.
    uint64_t effectiveForwardLength = forwardRunLength;
    uint64_t effectiveBackwardLength = backwardRunLength;

    if (forwardBaseFound && forwardRunBaseVal == siteBaseValue) {
        effectiveForwardLength++; // Site extends the forward run
    } else if (backwardBaseFound && backwardRunBaseVal == siteBaseValue) {
        // Site extends the backward run (and not the forward one, due to 'else if')
        effectiveBackwardLength++;
    }

    if (effectiveForwardLength >= homopolymerThreshold || effectiveBackwardLength >= homopolymerThreshold) {
        return true;
    }

    // If the site base, the forward run base, and the backward run base are all identical,
    // and the total length of (run_strictly_left + site + run_strictly_right) meets threshold.
    if (forwardBaseFound && backwardBaseFound &&
        siteBaseValue == forwardRunBaseVal && 
        backwardRunBaseVal == forwardRunBaseVal) {
        // All three (site, char of forward run, char of backward run) are identical.
        // The length to check is the sum of the original strict run lengths plus 1 for the site itself.
        if ((forwardRunLength + backwardRunLength + 1) >= homopolymerThreshold) {
            return true;
        }
    }

    return false;
}




bool isSiteInHomopolymerRegion(
    uint64_t sitePositionInRead,
    const shasta::LongBaseSequenceView& readSequence,
    uint64_t homopolymerContextLength = 12,  // HPC_PL (homopolymer flank length)
    uint64_t maxRepeatLength = 4,            // HPC_RR (homopolymer repeat range)
    uint64_t homopolymerCutoff = 2           // HPC_CC (homopolymer cutoff)
) {
    const uint64_t readLength = readSequence.baseCount;
    SHASTA_ASSERT(sitePositionInRead < readLength);

    // Convert to int64_t for compatibility with hifiasm logic
    int64_t p = static_cast<int64_t>(sitePositionInRead);
    int64_t sn = static_cast<int64_t>(readLength);
    int64_t hpc_flk = static_cast<int64_t>(homopolymerContextLength);
    int64_t hpc_rr = static_cast<int64_t>(maxRepeatLength);
    int64_t hpc_cutoff = static_cast<int64_t>(homopolymerCutoff);

    // Calculate boundaries (same as hifiasm)
    int64_t s = ((p >= hpc_flk) ? (p - hpc_flk) : 0);
    int64_t e = (((p + hpc_flk) <= sn) ? (p + hpc_flk) : sn);

    // Safety checks
    if (s < 0) s = 0;
    if (e > sn) e = sn;
    SHASTA_ASSERT(s < e);

    for (int64_t r = 1; r <= hpc_rr; r++) {
        int64_t rc = r * hpc_cutoff;
        
        // Pattern 1: Including site (forward and backward)
        int64_t k, zs, ze;
        
        // Forward scan from position p + r
        for (k = p + r; (k < e) && ((k - r) >= s) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k - r)].value); k++);
        ze = k;
        if (ze > e) ze = e;
        
        // Backward scan from position p - 1
        for (k = p - 1; (k >= s) && ((k + r) < e) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k + r)].value); k--);
        zs = k + 1;
        if (zs < s) zs = s;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }

        // Pattern 2: Excluding site (forward only)
        for (k = p + r + 1; (k < e) && ((k - r) >= s) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k - r)].value); k++);
        zs = p + 1;
        if (zs < s) zs = s;
        ze = k;
        if (ze > e) ze = e;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }

        // Pattern 3: Including site (backward and forward)
        for (k = p - r; (k >= s) && ((k + r) < e) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k + r)].value); k--);
        zs = k + 1;
        if (zs < s) zs = s;
        
        for (k = p + 1; (k < e) && ((k - r) >= s) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k - r)].value); k++);
        ze = k;
        if (ze > e) ze = e;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }

        // Pattern 4: Excluding site (backward only)
        for (k = p - r - 1; (k >= s) && ((k + r) < e) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k + r)].value); k--);
        zs = k + 1;
        if (zs < s) zs = s;
        ze = p;
        if (ze > e) ze = e;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }
    }
    
    return false;
}



// -----------------------------------------------------------------------------
// Overload: isSiteInHomopolymerRegion accepting a vector<Base> sequence
// -----------------------------------------------------------------------------
inline bool isSiteInHomopolymerRegion(
    uint64_t sitePositionInRead,
    const std::vector<shasta::Base>& readSequence,
    uint64_t homopolymerContextLength = 12,  // HPC_PL (homopolymer flank length)
    uint64_t maxRepeatLength        = 4,     // HPC_RR (homopolymer repeat range)
    uint64_t homopolymerCutoff      = 2)     // HPC_CC (homopolymer cutoff)
{
    const uint64_t readLength = readSequence.size();
    SHASTA_ASSERT(sitePositionInRead < readLength);

    // Convert to int64_t for computation identical to the other overload.
    int64_t p          = static_cast<int64_t>(sitePositionInRead);
    int64_t sn         = static_cast<int64_t>(readLength);
    int64_t hpc_flk    = static_cast<int64_t>(homopolymerContextLength);
    int64_t hpc_rr     = static_cast<int64_t>(maxRepeatLength);
    int64_t hpc_cutoff = static_cast<int64_t>(homopolymerCutoff);

    int64_t s = (p >= hpc_flk) ? (p - hpc_flk) : 0;
    int64_t e = ((p + hpc_flk) <= sn) ? (p + hpc_flk) : sn;
    if(s < 0) s = 0;
    if(e > sn) e = sn;
    SHASTA_ASSERT(s < e);

    for(int64_t r = 1; r <= hpc_rr; ++r) {
        int64_t rc = r * hpc_cutoff;
        int64_t k, zs, ze;

        // Pattern 1: Including site (forward & backward)
        for(k = p + r; (k < e) && ((k - r) >= s) && (readSequence[k].value == readSequence[k - r].value); ++k) {}
        ze = (k > e) ? e : k;
        for(k = p - 1; (k >= s) && ((k + r) < e) && (readSequence[k].value == readSequence[k + r].value); --k) {}
        zs = k + 1;
        if(zs < s) zs = s;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;

        // Pattern 2: Excluding site (forward only)
        for(k = p + r + 1; (k < e) && ((k - r) >= s) && (readSequence[k].value == readSequence[k - r].value); ++k) {}
        zs = p + 1;
        if(zs < s) zs = s;
        ze = (k > e) ? e : k;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;

        // Pattern 3: Including site (backward & forward)
        for(k = p - r; (k >= s) && ((k + r) < e) && (readSequence[k].value == readSequence[k + r].value); --k) {}
        zs = k + 1;
        if(zs < s) zs = s;
        for(k = p + 1; (k < e) && ((k - r) >= s) && (readSequence[k].value == readSequence[k - r].value); ++k) {}
        ze = (k > e) ? e : k;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;

        // Pattern 4: Excluding site (backward only)
        for(k = p - r - 1; (k >= s) && ((k + r) < e) && (readSequence[k].value == readSequence[k + r].value); --k) {}
        zs = k + 1;
        if(zs < s) zs = s;
        ze = p;
        if(ze > e) ze = e;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;
    }

    return false;
}


// Function to check for a specific strand bias pattern
// Returns true if the strand bias is detected, and populates outDominantStrandHet1 and outDominantStrandHet2.
// dominantStrand values: 0 (all on strand 0), 1 (all on strand 1), 2 (no dominance or not enough reads).
bool hasSiteStrandBiasReadGraph5(
    const std::set<OrientedReadId>& hetBase1OrientedReads,
    const std::set<OrientedReadId>& hetBase2OrientedReads,
    uint64_t& outDominantStrandHet1,
    uint64_t& outDominantStrandHet2)
{
    // For this specific bias pattern, we need at least one read supporting each allele's strand pattern.
    const uint64_t minReadsPerAlleleForThisBias = 1;

    outDominantStrandHet1 = 2; // Initialize to no dominance
    if (hetBase1OrientedReads.size() >= minReadsPerAlleleForThisBias) {
        uint64_t hetBase1CountStrand0 = 0;
        uint64_t hetBase1CountStrand1 = 0;
        for (const auto& hetBase1OrientedRead : hetBase1OrientedReads) {
            if (hetBase1OrientedRead.getStrand() == 0) {
                hetBase1CountStrand0++;
            } else {
                hetBase1CountStrand1++;
            }
        }
        if (hetBase1CountStrand0 > 0 && hetBase1CountStrand1 == 0) {
            outDominantStrandHet1 = 0; // All reads for hetBase1 are on strand 0
        } else if (hetBase1CountStrand1 > 0 && hetBase1CountStrand0 == 0) {
            outDominantStrandHet1 = 1; // All reads for hetBase1 are on strand 1
        }
    }

    outDominantStrandHet2 = 2; // Initialize to no dominance
    if (hetBase2OrientedReads.size() >= minReadsPerAlleleForThisBias) {
        uint64_t hetBase2CountStrand0 = 0;
        uint64_t hetBase2CountStrand1 = 0;
        for (const auto& hetBase2OrientedRead : hetBase2OrientedReads) {
            if (hetBase2OrientedRead.getStrand() == 0) {
                hetBase2CountStrand0++;
            } else {
                hetBase2CountStrand1++;
            }
        }
        if (hetBase2CountStrand0 > 0 && hetBase2CountStrand1 == 0) {
            outDominantStrandHet2 = 0; // All reads for hetBase2 are on strand 0
        } else if (hetBase2CountStrand1 > 0 && hetBase2CountStrand0 == 0) {
            outDominantStrandHet2 = 1; // All reads for hetBase2 are on strand 1
        }
    }

    // Check for the specific bias pattern:
    // Both alleles must show complete strand dominance, and on opposite strands.
    if (outDominantStrandHet1 != 2 && outDominantStrandHet2 != 2 && outDominantStrandHet1 != outDominantStrandHet2) {
        return true; // Bias detected
    }

    return false; // No bias detected
}








class HaplotypeEvidence{    
    public:
        OrientedReadId orientedReadId; // The oriented read ID
        uint64_t positionInOrientedRead = invalid<uint64_t>; // The position in the oriented read
        AlignedBase baseOfOrientedRead; // Base value (0=A, 1=C, 2=G, 3=T, 4=Gap)
        uint64_t alignmentId;
        bool isInHomopolymerRegion = false;
};


class AlignmentPositionBaseStatsSingleEntry{
    public:

        uint64_t positionInTargetRead = invalid<uint64_t>;
        AlignedBase baseOfTargetRead; // Base value (0=A, 1=C, 2=G, 3=T, 4=Gap)

        uint64_t coverage = 0;

        HaplotypeEvidence haplotypeEvidence;

        bool isInHomopolymerRegion = false;

        // Add a comparison operator to allow sorting by position.
        bool operator<(const AlignmentPositionBaseStatsSingleEntry& other) const
        {
            return positionInTargetRead < other.positionInTargetRead;
        }
};

class AlignmentPositionBaseStatsMergedEntries{
    public:

        uint64_t positionInTargetRead = invalid<uint64_t>;
        AlignedBase baseOfTargetRead; // Base value (0=A, 1=C, 2=G, 3=T, 4=Gap)

        uint64_t haplotypeEvidenceCoverage = 0;

        std::array<std::vector<HaplotypeEvidence>, 5> haplotypeEvidenceByBase;

        bool isInHomopolymerRegion = false;

        // Add a comparison operator to allow sorting by position.
        bool operator<(const AlignmentPositionBaseStatsMergedEntries& other) const
        {
            return positionInTargetRead < other.positionInTargetRead;
        }
};

class SnpStats{
    public:
    
        uint64_t positionInTargetRead = invalid<uint64_t>;

        AlignedBase hetBase0;
        vector<OrientedReadId> hetBase0OrientedReadIds;
        bool isInHomopolymerRegion = false;

        AlignedBase hetBase1;
        vector<HaplotypeEvidence> hetBase1OrientedReadIdsHE;

        uint64_t occ_0 = 0;
        uint64_t occ_1 = 0;

        uint64_t overlap_num = 0;

        bool isHpcVector = false;
        bool hasSufficentCoverage = false;
        int64_t score = -1;

        // Add a comparison operator to allow sorting by position.
        bool operator<(const SnpStats& other) const
        {
            return positionInTargetRead < other.positionInTargetRead;
        }
        
};


class Site {
    public:
        std::set<OrientedReadId> orientedReads;
        OrientedReadId targetOrientedReadId;
        std::set<OrientedReadId> excludedOrientedReads;
    };


// Helper: check quickly if a position is heterozygous for a given oriented read.
// 'vec' must be the per-read sorted vector filled in hetPos. Uses binary_search.
static inline bool isHetSite(const std::vector<uint32_t>& v, uint32_t p)
{
    return std::binary_search(v.begin(), v.end(), p);
}


// orientedReadHetPositions[v]  = sorted list of positions that look heterozygous
//              for oriented read with value v (0 … 2*readCount-1).
// NOTE: positions are always stored in forward-strand coordinates.
//       For strand-1 we just mirror them on first access
//       (posRC = readLen - 1 - posFwd).
std::vector< std::vector<uint32_t> > orientedReadHetPositions;



// Structure to hold data for the phasing threads
struct PhasingThreadDataReadGraph5 {
    // Pointers to assembler data (const access needed in thread)
    const Assembler *assembler; // Use const Assembler*
    uint64_t alignmentCount;
    uint64_t readCount; // Add readCount for checks
    uint64_t orientedReadCount; // Add orientedReadCount for checks

    // Parameters needed by the thread function
    // (Add any other parameters from createReadGraph5withStrandSeparation if needed inside the loop)

    // Thread-local storage for results
    // Each outer vector is indexed by threadId
    // Each inner vector<bool> is indexed by alignmentId
    vector<vector<bool> > threadForbiddenAlignments;
    vector<vector<bool> > threadFirstPassHetAlignments;
    vector<vector<bool> > threadAlignmentsAlreadyConsidered;
    vector<vector<Site> > threadSites;
    vector<bool> isReadIdContained;

    

    // Mutex for thread-safe cout if debugging is needed inside threads
    // std::mutex coutMutex; // Uncomment if needed

    PhasingThreadDataReadGraph5(const Assembler *asmPtr, size_t threadCount) : assembler(asmPtr) {
        // --- Check asmPtr validity early ---
        SHASTA_ASSERT(assembler != nullptr);

        alignmentCount = assembler->alignmentData.size();
        readCount = assembler->getReads().readCount();
        orientedReadCount = readCount * 2;

        // --- Check threadCount validity ---
        SHASTA_ASSERT(threadCount > 0);

        threadForbiddenAlignments.resize(threadCount);
        threadFirstPassHetAlignments.resize(threadCount);
        threadAlignmentsAlreadyConsidered.resize(threadCount);
        threadSites.resize(threadCount);

        for (size_t i = 0; i < threadCount; ++i) {
            // Resize inner vectors
            threadForbiddenAlignments[i].resize(alignmentCount, false);
            threadFirstPassHetAlignments[i].resize(alignmentCount, false);
            threadAlignmentsAlreadyConsidered[i].resize(alignmentCount, false);
        }

        

    }
};

// Global pointer for thread access (similar pattern to computeAlignmentsData)
PhasingThreadDataReadGraph5 *phasingThreadDataReadGraph5 = nullptr;





// Thread function for the phasing analysis part
void Assembler::createReadGraph5ThreadFunction(uint64_t threadId) {

    ofstream debugOut("Debug-" + to_string(threadId) + ".txt");
    debugOut << "Thread ID: " << threadId << " is starting read phasing" << endl;

    // --- Ensure phasingThreadDataReadGraph5 is valid ---
    SHASTA_ASSERT(phasingThreadDataReadGraph5 != nullptr);

    // Get reference to thread-local storage
    SHASTA_ASSERT(threadId < phasingThreadDataReadGraph5->threadForbiddenAlignments.size()); // Check threadId validity
    vector<bool>& forbiddenAlignments = phasingThreadDataReadGraph5->threadForbiddenAlignments[threadId];
    vector<bool>& firstPassHetAlignments = phasingThreadDataReadGraph5->threadFirstPassHetAlignments[threadId];
    vector<bool>& alignmentsAlreadyConsidered = phasingThreadDataReadGraph5->threadAlignmentsAlreadyConsidered[threadId];
    vector<Site>& sites = phasingThreadDataReadGraph5->threadSites[threadId];
    const uint64_t alignmentCount = phasingThreadDataReadGraph5->alignmentCount;
    const uint64_t orientedReadCount = phasingThreadDataReadGraph5->orientedReadCount; // Get from shared data

    // Define types for the interval map
    using Interval = boost::icl::discrete_interval<uint32_t>;
    // Store pairs of (alignmentId, the other OrientedReadId relative to orientedReadId0)
    using AlignmentInfoPair = std::pair<uint64_t, OrientedReadId>;
    // Define a custom comparison for the set to handle pairs
    struct CompareAlignmentInfoPair {
        bool operator()(const AlignmentInfoPair& a, const AlignmentInfoPair& b) const {
            if (a.first != b.first) {
                return a.first < b.first; // Compare alignmentId first
            }
            return a.second < b.second; // Then compare OrientedReadId
        }
    };
    using AlignmentInfoSet = std::set<AlignmentInfoPair, CompareAlignmentInfoPair>;
    using AlignmentIntervalMap = boost::icl::interval_map<uint32_t, AlignmentInfoSet>;


    // Access assembler data via the pointer
    const Assembler& assembler = *(phasingThreadDataReadGraph5->assembler);




    // Get batches of read IDs to process
    uint64_t readIdBegin;
    uint64_t readIdEnd; // variable to receive the end of the batch
    while (getNextBatch(readIdBegin, readIdEnd)) {

        debugOut << "Thread ID: " << threadId << " is processing read IDs from " << readIdBegin << " to " << readIdEnd << endl;

        // --- Check batch validity ---
        SHASTA_ASSERT(readIdBegin < phasingThreadDataReadGraph5->readCount);

        // Ensure readIdEnd doesn't wrap around or go below readIdBegin after clamping
        SHASTA_ASSERT(readIdEnd >= readIdBegin);

        // Loop over read IDs in this batch
        for (ReadId readId = readIdBegin; readId < readIdEnd; readId++) {

            // --- Start of per-read analysis code ---
            // if ((readId != 2228) && (readId != 2229) && (readId != 2230) && (readId != 2231) && (readId != 2232) && (readId != 2233) && (readId != 2234) && (readId != 2235) && (readId != 2236) && (readId != 2237)) { // Keep this for debugging specific reads if needed
            //     continue;
            // }

            // if (((readId < 19300) || (readId > 19420))) { // Keep this for debugging specific reads if needed
            //     if ((readId < 1980) || (readId > 2050)) {
            //         continue;
            //     }
            // }

            // // --- Start of per-read analysis code ---
            // if ((readId != 0)) {
            //     continue;
            // }

            debugOut << "Thread ID: " << threadId << " is processing read ID: " << readId << endl;
            debugOut << timestamp << endl;

            const ReadId targetReadId = readId;
            const Strand targetReadStrand = 0; // Analyze strand 0 arbitrarily, results should be consistent
            const OrientedReadId targetOrientedReadId(targetReadId, targetReadStrand);


            // XXX
            // --- START OF: Build Interval Tree for alignments relative to targetOrientedReadId ---
            //

            AlignmentIntervalMap alignmentIntervals; // Create the interval map for this targetReadId
            
            debugOut << timestamp << "Thread ID: " << threadId << " Building interval tree for read ID: " << targetReadId << endl;

            const auto alignmentTableForIntervalTreeConstruction = alignmentTable[targetOrientedReadId.getValue()];
            for (const auto alignmentId : alignmentTableForIntervalTreeConstruction) {

                const AlignmentData& thisAlignmentData = alignmentData[alignmentId];

                OrientedReadId currentOrientedReadId0(thisAlignmentData.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thisAlignmentData.info;

                Alignment alignment;
                span<const char> compressedAlignment = assembler.compressedAlignments[alignmentId];
                shasta::decompress(compressedAlignment, alignment);
                
                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != targetReadId) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                    alignment.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetReadId);

                uint32_t markerCount0 = uint32_t(markers[currentOrientedReadId0.getValue()].size());
                uint32_t markerCount1 = uint32_t(markers[currentOrientedReadId1.getValue()].size());

                // Reverse complement if necessary
                if (currentOrientedReadId0.getStrand() != targetReadStrand) {
                    // Ensure marker counts are valid before calling reverseComplement
                    alignment.reverseComplement(markerCount0, markerCount1);
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == targetReadStrand);
                SHASTA_ASSERT(currentOrientedReadId0 == targetOrientedReadId); // Check consistency

                const auto orientedReadMarkers = markers[currentOrientedReadId0.getValue()];
                const uint32_t position0 = orientedReadMarkers[alignmentInfo.data[0].firstOrdinal].position;
                const uint32_t position1 = orientedReadMarkers[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assembler.assemblerInfo->k);

                SHASTA_ASSERT(position0 < position1);

                // --- Add interval to the map ---
                if (position1 > position0) { // Ensure a valid interval
                    Interval interval = Interval::closed(position0, position1); // Create a [position0, position1) interval
                    // Create a pair with alignmentId and the correctly oriented other read
                    AlignmentInfoPair infoPair = {alignmentId, currentOrientedReadId1};
                    AlignmentInfoSet currentAlignmentInfoSet = {infoPair}; // Create a set with the pair
                    alignmentIntervals.add({interval, currentAlignmentInfoSet}); // Add to the map (handles overlaps)
                    debugOut << "Thread ID: " << threadId << " Added interval [" << position0 << ", " << position1 << ") for alignment " << alignmentId << " (Other Read: " << currentOrientedReadId1 << ")" << endl;
                } else {
                    // debugOut << "Thread ID: " << threadId << " Skipping invalid interval [" << position0 << ", " << position1 << "] for alignment " << alignmentId << endl;
                }
            }

            debugOut << timestamp << "Thread ID: " << threadId << " Finished building interval tree for read ID: " << targetReadId << endl;

            // XXX
            // --- END OF: Build Interval Tree for alignments relative to targetOrientedReadId ---
            //
            

            debugOut << timestamp << "Thread ID: " << threadId << " Starting looping over alignments for read ID: " << targetReadId << endl;







            // std::map<uint64_t, AlignmentPositionBaseStats> positionStatsOnTargetOrientedReadId;
            vector<AlignmentPositionBaseStatsSingleEntry> positionStatsOnTargetOrientedReadIdSingleEntries;
            vector<AlignmentPositionBaseStatsMergedEntries> positionStatsOnTargetOrientedReadIdMergedEntries;
            vector<SnpStats> potentialHetSitesOnTargetOrientedReadId;
            vector<SnpStats> potentialHetSitesOnTargetOrientedReadIdFiltered;


            // Loop over alignments involving the target read (targetOrientedReadId)
            const auto alignmentTableForTargetOrientedReadId = alignmentTable[targetOrientedReadId.getValue()];

            
            // // XXX
            // // --- Find and forbid alignments between the same read on different strands.
            // //
            // vector <bool> readIdsUsed(phasingThreadDataReadGraph5->readCount, false);
            // vector <bool> forbiddenReadIds(phasingThreadDataReadGraph5->readCount, false);

            // for (const auto alignmentId : alignmentTableForTargetOrientedReadId) {
                
            //     const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
            //     OrientedReadId currentOrientedReadId0(thisAlignmentData.readIds[0], 0);
            //     OrientedReadId currentOrientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
            //     AlignmentInfo alignmentInfo = thisAlignmentData.info;

            
            //     const ReadId currentReadId0 = thisAlignmentData.readIds[0];
            //     const ReadId currentReadId1 = thisAlignmentData.readIds[1];
            //     if (currentReadId0 == targetReadId) {
            //         if(readIdsUsed[currentReadId1]) {
            //             forbiddenAlignments[alignmentId] = true;
            //             forbiddenReadIds[currentReadId1] = true;
            //             continue;
            //         }
            //         readIdsUsed[currentReadId1] = true;
            //     } else if (currentReadId1 == targetReadId) {
            //         if(readIdsUsed[currentReadId0]) {
            //             forbiddenAlignments[alignmentId] = true;
            //             forbiddenReadIds[currentReadId0] = true;
            //             continue;
            //         }
            //         readIdsUsed[currentReadId0] = true;                
            //     }
                
            // }

            // // We need a second loop over the alignments to forbid the second alignment instance in some cases
            // for (const auto alignmentId : alignmentTableForTargetOrientedReadId) {
                
            //     const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
            //     const ReadId currentReadId0 = thisAlignmentData.readIds[0];
            //     const ReadId currentReadId1 = thisAlignmentData.readIds[1];
            //     if (currentReadId0 == targetReadId) {
            //         if(forbiddenReadIds[currentReadId1]) {
            //             forbiddenAlignments[alignmentId] = true;
            //         }
            //     } else if (currentReadId1 == targetReadId) {
            //         if(forbiddenReadIds[currentReadId0]) {
            //             forbiddenAlignments[alignmentId] = true;
            //         }               
            //     }

            // }

            // // XXX
            // // --- END OF: Find and forbid alignments between the same read on different strands.
            // //

            

            // XXX
            // --- We found and forbid alignments between the same read on different strands.
            //     Now we can proceed with the analysis of the rest of alignments.

            for (const auto alignmentId : alignmentTableForTargetOrientedReadId) {

                if (forbiddenAlignments[alignmentId]) {
                    // Skip alignments that are forbidden due to strand issues
                    continue;
                }

                const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
                OrientedReadId currentOrientedReadId0(thisAlignmentData.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thisAlignmentData.info;

                Alignment alignment;
                span<const char> compressedAlignment = assembler.compressedAlignments[alignmentId];
                shasta::decompress(compressedAlignment, alignment);

                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != targetReadId) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                    alignment.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetReadId);

                uint32_t markerCount0 = uint32_t(markers[currentOrientedReadId0.getValue()].size());
                uint32_t markerCount1 = uint32_t(markers[currentOrientedReadId1.getValue()].size());
                
                // Reverse complement if necessary
                if (currentOrientedReadId0.getStrand() != targetReadStrand) {
                    // Ensure marker counts are valid before calling reverseComplement
                    alignment.reverseComplement(markerCount0, markerCount1);
                    // alignment.checkStrictlyIncreasing(); // Can be expensive, maybe remove in production
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == targetReadStrand);
                
                debugOut << timestamp << "Thread ID: " << threadId << " Processing alignment ID: " << alignmentId 
                         << " between ReadId: " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() 
                         << " and ReadId: " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << endl;

                
                // XXX
                // --- Project this alignment to base space.
                //

                ProjectedAlignment projectedAlignment(
                    assembler,
                    {currentOrientedReadId0, currentOrientedReadId1},
                    alignment,
                    ProjectedAlignment::Method::QuickRaw
                );

                // XXX
                // --- END OF: Project this alignment to base space.
                //

                debugOut << timestamp << "Thread ID: " << threadId << " Projected alignment for alignment ID: " << alignmentId 
                         << " between ReadId: " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() 
                         << " and ReadId: " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() 
                         << " with segments count: " << projectedAlignment.segments.size() << endl;

                
                // Loop over the RAW segments of the projectedAlignment.
                for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
                    
                    // Get the RAW sequences
                    const vector<Base>& sequence0 = segment.sequences[0];
                    const vector<Base>& sequence1 = segment.sequences[1];

                    // Align them base by base to get the right sequence alignment representation.
                    uint64_t position0 = 0;
                    uint64_t position1 = 0;
                    vector<AlignedBase> rawAlignmentSequence0;
                    vector<AlignedBase> rawAlignmentSequence1;

                    // Retreive the RAW alignment sequences.
                    for (const pair<bool, bool>& p: segment.alignment) {
                        const bool hasBase0 = p.first;
                        const bool hasBase1 = p.second;
                        
                        AlignedBase baseAtPosition0;
                        AlignedBase baseAtPosition1;
                        uint64_t positionWithDetectedChange0 = 0;
                        uint64_t positionWithDetectedChange1 = 0;
                        if(hasBase0) {
                            SHASTA_ASSERT(position0 < sequence0.size());
                            baseAtPosition0 = AlignedBase(sequence0[position0]);
                            rawAlignmentSequence0.push_back(AlignedBase(sequence0[position0]));
                            positionWithDetectedChange0 = position0;
                            position0++;
                        } else {
                            rawAlignmentSequence0.push_back(AlignedBase::gap());
                            baseAtPosition0 = AlignedBase::gap();
                        }

                        if(hasBase1) {
                            SHASTA_ASSERT(position1 < sequence1.size());
                            baseAtPosition1 = AlignedBase(sequence1[position1]);
                            rawAlignmentSequence1.push_back(AlignedBase(sequence1[position1]));
                            positionWithDetectedChange1 = position1;
                            position1++;
                        } else {
                            rawAlignmentSequence1.push_back(AlignedBase::gap());
                            baseAtPosition1 = AlignedBase::gap();
                        }

                        // We need to keep track of the gaps even though we do not rely on them for phasing.
                        // The number of alignments with gaps in this position is important to infer statistics
                        // about this position later (like the total alignments in that position).
                        if (hasBase0 && (baseAtPosition0 != baseAtPosition1)) {

                            // //
                            // // --- DEBUG: Print the RAW alignment sequences that have differences ---
                            // //            AND the differences are not gaps
                            // //
                            // if(hasBase1) {
                            //     debugOut << endl;
                            //     debugOut << "rawAlignmentSequence0 marker positionsA start: " << segment.positionsA[0] << endl;
                            //     debugOut << "rawAlignmentSequence1 marker positionsA start: " << segment.positionsA[1] << endl;
                            //     debugOut << "sequence0 base change positionsA start: " << segment.positionsA[0] + positionWithDetectedChange0 << endl;
                            //     debugOut << "sequence1 base change positionsA start: " << segment.positionsA[1] + positionWithDetectedChange1 << endl;
                            //     debugOut << "sequence0 base at base change positionsA: " << baseAtPosition0 << endl;
                            //     debugOut << "sequence1 base at base change positionsA: " << baseAtPosition1 << endl;
                            //     debugOut << "We found " << segment.editDistance << " edit distance differences between the RAW sequences" << endl;
                            //     debugOut << "rawAlignmentSequence0: ";
                            //     for(uint64_t i=0; i<rawAlignmentSequence0.size(); i++) {
                            //         const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
                            //         if(isDifferent) {
                            //             debugOut << "[";
                            //         }
                            //         debugOut << rawAlignmentSequence0[i].character();
                            //         if(isDifferent) {
                            //             debugOut << "]";
                            //         }
                            //     }
                            //     debugOut << endl;

                            //     debugOut << "rawAlignmentSequence1: ";
                            //     for(uint64_t i=0; i<rawAlignmentSequence1.size(); i++) {
                            //         const bool isDifferent = (rawAlignmentSequence0[i] != rawAlignmentSequence1[i]);
                            //         if(isDifferent) {
                            //             debugOut << "[";
                            //         }
                            //         debugOut << rawAlignmentSequence1[i].character();
                            //         if(isDifferent) {
                            //             debugOut << "]";
                            //         }
                            //     }
                            //     debugOut << endl;
                            // }
                            // //
                            // // --- END DEBUG ---
                            // //




                            // XXX
                            // --- Update statistics for positionInTheTargetRead ---
                            //

                            uint64_t positionInTheTargetRead = segment.positionsA[0] + positionWithDetectedChange0;
                            uint64_t positionInTheOtherRead = segment.positionsA[1] + positionWithDetectedChange1;

                            HaplotypeEvidence haplotypeEvidence;
                            haplotypeEvidence.orientedReadId = currentOrientedReadId1;
                            haplotypeEvidence.positionInOrientedRead = positionInTheOtherRead;
                            haplotypeEvidence.baseOfOrientedRead = baseAtPosition1;
                            haplotypeEvidence.alignmentId = alignmentId;
                            haplotypeEvidence.isInHomopolymerRegion = false;

                            // XXX
                            // --- START OF: Check if this potential het site in the other read
                            //     is in a homopolymer region.
                            //

                            const LongBaseSequenceView currentOrientedReadId1Sequence = reads->getRead(currentOrientedReadId1.getReadId());

                            uint64_t homopolymerContextLength = 12;
                            uint64_t maxRepeatLength = 4;
                            uint64_t homopolymerCutoff = 2;

                            if (currentOrientedReadId1.getStrand() == 0) {
                                bool isInHomopolymerRegion = isSiteInHomopolymerRegion(positionInTheOtherRead, currentOrientedReadId1Sequence, homopolymerContextLength, maxRepeatLength, homopolymerCutoff);
                                if (isInHomopolymerRegion) {
                                    haplotypeEvidence.isInHomopolymerRegion = true;
                                }
                            }

                            // If the other read is on the reverse strand, we need to convert the position to the forward strand.
                            if (currentOrientedReadId1.getStrand() == 1) {
                                uint64_t positionInTheOtherReadConvertedToForwardStrand = (currentOrientedReadId1Sequence.baseCount - 1) - positionInTheOtherRead;
                                bool isInHomopolymerRegion = isSiteInHomopolymerRegion(positionInTheOtherReadConvertedToForwardStrand, currentOrientedReadId1Sequence, homopolymerContextLength, maxRepeatLength, homopolymerCutoff);
                                if (isInHomopolymerRegion) {
                                    haplotypeEvidence.isInHomopolymerRegion = true;
                                }
                            }
    
                            // XXX
                            // --- END OF: Check if this potential het site in the other read ---
                            //     is in a homopolymer region.
                            //

                            AlignmentPositionBaseStatsSingleEntry thisPositionStatsSingleEntry;
                            thisPositionStatsSingleEntry.positionInTargetRead = positionInTheTargetRead;
                            thisPositionStatsSingleEntry.baseOfTargetRead = baseAtPosition0;
                            thisPositionStatsSingleEntry.haplotypeEvidence = haplotypeEvidence;
                            thisPositionStatsSingleEntry.coverage = 1;
                            positionStatsOnTargetOrientedReadIdSingleEntries.push_back(thisPositionStatsSingleEntry);


                        } // --- End of printing the sequences that have differences and updating the stats of this position ---

                    } // --- End of loop over the segment ---
                    SHASTA_ASSERT(rawAlignmentSequence0.size() == rawAlignmentSequence1.size());

                } // --- End loop over segments ---

                // debugOut << "Loop over the RAW segments of the Projected alignment between ReadId: " << currentOrientedReadId0.getReadId() << "-" << currentOrientedReadId0.getStrand() << " and ReadId: " << currentOrientedReadId1.getReadId() << "-" << currentOrientedReadId1.getStrand() << " completed." << endl;  
            
            } // --- End loop over alignments for targetReadId ---

            debugOut << timestamp << "Thread ID: " << threadId << " Finished looping over alignments for read ID: " << targetReadId << endl;

            // Sort the positionStatsOnTargetOrientedReadId by positionInTargetRead.
            std::sort(positionStatsOnTargetOrientedReadIdSingleEntries.begin(), positionStatsOnTargetOrientedReadIdSingleEntries.end());


            for(auto streakBegin=positionStatsOnTargetOrientedReadIdSingleEntries.begin(); streakBegin!=positionStatsOnTargetOrientedReadIdSingleEntries.end(); ) {
                auto streakEnd = streakBegin + 1;

                // Find the streak of positions with the same positionInTargetRead.
                for(;
                    streakEnd!=positionStatsOnTargetOrientedReadIdSingleEntries.end() && streakEnd->positionInTargetRead==streakBegin->positionInTargetRead;
                    streakEnd++) {
                }

                // Create a merged representation of the streak.
                AlignmentPositionBaseStatsMergedEntries thisPositionStatsMergedEntries;
                thisPositionStatsMergedEntries.positionInTargetRead = streakBegin->positionInTargetRead;
                thisPositionStatsMergedEntries.baseOfTargetRead = streakBegin->baseOfTargetRead;

                for(auto it=streakBegin; it!=streakEnd; it++) {
                    thisPositionStatsMergedEntries.haplotypeEvidenceByBase[it->haplotypeEvidence.baseOfOrientedRead.value].push_back(it->haplotypeEvidence);
                    thisPositionStatsMergedEntries.haplotypeEvidenceCoverage += 1;
                }

                positionStatsOnTargetOrientedReadIdMergedEntries.push_back(thisPositionStatsMergedEntries);

                // Process the next streak.
                streakBegin = streakEnd;

            }


            // // DEBUG: Print the positionStatsOnTargetOrientedReadIdMergedEntries.
            // for(const auto& positionStatsMergedEntries : positionStatsOnTargetOrientedReadIdMergedEntries) {
            //     debugOut << "Position: " << positionStatsMergedEntries.positionInTargetRead << " Base: " << positionStatsMergedEntries.baseOfTargetRead.character() << " Total Coverage: " << positionStatsMergedEntries.totalCoverage << endl;
            // }

            
            
            // XXX
            // --- Analyze potential heterozygous sites ---
            //

            // We finished analyzing all alignments for the target read (targetReadId).
            // Now we need to check each potential site in targetReadId (in positionStatsOnOrientedReadId0)
            // to see if it involves a heterozygous site.
            uint64_t sitesSkippedDueToInsufficientCoverage = 0;
            for (auto positionStatsOnTargetOrientedRead: positionStatsOnTargetOrientedReadIdMergedEntries) {


                // XXX
                // --- Check if this potential het site in the target read ---
                //     is in a homopolymer region. The target reas is always 
                //     on the forward strand.
                //

                // ShortBaseSequence8 kmer;
                ShortBaseSequence32 kmer;

                const LongBaseSequenceView targetReadSequence = reads->getRead(targetReadId);

                // Extract the kmer from the targetReadId sequence that includes
                // the 12 preceding and 12 following bases on that position.
                // Ensure we don't go out of bounds
                uint64_t kmerStart = (positionStatsOnTargetOrientedRead.positionInTargetRead >= 12) ? positionStatsOnTargetOrientedRead.positionInTargetRead - 12 : 0;
                uint64_t kmerLength = min(25UL, targetReadSequence.baseCount - kmerStart);
                extractKmer(targetReadSequence, kmerStart, kmerLength, kmer);
                
                uint64_t homopolymerContextLength = 12;
                uint64_t maxRepeatLength = 4;
                uint64_t homopolymerCutoff = 2;
                bool siteIsInHomopolymerRegion = isSiteInHomopolymerRegion(positionStatsOnTargetOrientedRead.positionInTargetRead, targetReadSequence, homopolymerContextLength, maxRepeatLength, homopolymerCutoff);
                if (siteIsInHomopolymerRegion) {
                    // debugOut << "Position " << positionStatsOnTargetOrientedRead.positionInTargetRead << " is in a homopolymer region." << std::endl;
                    // kmer.write(debugOut, kmerLength) << std::endl;
                    positionStatsOnTargetOrientedRead.isInHomopolymerRegion = true;
                }


                // XXX
                // --- END OF: Check if this potential het site in the target read ---
                //     is in a homopolymer region.
                //
            

                // Get the number of alignment supporting each base.
                const uint64_t totalNumberOfA = positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[0].size();
                const uint64_t totalNumberOfC = positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[1].size();
                const uint64_t totalNumberOfG = positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[2].size();
                const uint64_t totalNumberOfT = positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[3].size();
                const uint64_t totalNumberOfGap = positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[4].size();

                uint64_t totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead =
                    totalNumberOfA +
                    totalNumberOfC +
                    totalNumberOfG +
                    totalNumberOfT +
                    totalNumberOfGap;


                // This is the minimum needed number of alignments that
                // support a change from the base of the target read in that position
                uint64_t minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead = 3;

                bool baseAHasEnoughCoverage = totalNumberOfA >= minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead;
                bool baseCHasEnoughCoverage = totalNumberOfC >= minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead;
                bool baseGHasEnoughCoverage = totalNumberOfG >= minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead;
                bool baseTHasEnoughCoverage = totalNumberOfT >= minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead;
                bool baseGapHasEnoughCoverage = totalNumberOfGap >= minimumNumberOfAlignmentsSupportingAChangeFromTheTargetRead;


                // Skip this position if there is no allele base with enough coverage
                if (!(baseAHasEnoughCoverage || baseCHasEnoughCoverage || baseGHasEnoughCoverage || baseTHasEnoughCoverage)) {
                    continue;
                }


                //
                // This is a potential heterozygous site
                //


                // Find the orientedReadIds supporting a change from the target read
                std::set<OrientedReadId> orientedReadIdsSupportingAChangeFromTheTargetRead;
                for (const auto& haplotypeEvidenceOfBase : positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase) {
                    for (const HaplotypeEvidence& haplotypeEvidence : haplotypeEvidenceOfBase) {
                        orientedReadIdsSupportingAChangeFromTheTargetRead.insert(haplotypeEvidence.orientedReadId);
                    }
                }


                // We don't have the number of alignments that agree with the base of the target read
                // in that specific position. We only have the number of alignments that support a change.
                // So, we need to find the total coverage of that position from the intervalTree.
                uint64_t totalCoverageInThisPosition = 0;
                std::set<OrientedReadId> coveringOrientedReadIds;
                std::set<uint64_t> coveringAlignmentIds;
                
                auto it = alignmentIntervals.find(positionStatsOnTargetOrientedRead.positionInTargetRead);
                if (it != alignmentIntervals.end()) {
                    
                    // Retrieve the covering alignment info pairs of this position,
                    // which include the pairs (alignmentId, orientedReadId of the other read).
                    const AlignmentInfoSet& coveringAlignmentInfo = it->second;
                    
                    // Calculate the total coverage in this position. 
                    totalCoverageInThisPosition = coveringAlignmentInfo.size();

                    // Iterate through the pairs in coveringAlignmentInfo to
                    // extract the alignmentId and the orientedReadId of the other read.
                    for (const AlignmentInfoPair& infoPair : coveringAlignmentInfo) {
                        uint64_t alignmentId = infoPair.first;
                        OrientedReadId storedOrientedReadId1 = infoPair.second;
                        coveringAlignmentIds.insert(alignmentId);
                        coveringOrientedReadIds.insert(storedOrientedReadId1);
                    }

                }

                // Add the target read to the coverage
                totalCoverageInThisPosition += 1; 



                // Calculate the total coverage of the reads supporting the target read
                uint64_t totalNumberOfAlignmentsSupportingTheTargetRead = totalCoverageInThisPosition - totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead;


                // Skip this position if the coverage supporting the target read is less than 3
                uint64_t minimumNumberOfAlignmentsSupportingTheTargetRead = 3;
                if (totalNumberOfAlignmentsSupportingTheTargetRead < minimumNumberOfAlignmentsSupportingTheTargetRead) {
                    sitesSkippedDueToInsufficientCoverage++;
                    // debugOut << "Skipping position " << positionInRead0 << " due to insufficient coverage supporting the target read: " << totalCoverageInThisPosition - totalNumberOfAlignmentsSupportingAChangeFromTheTargetRead << endl;
                    // debugOut << "Position: " << positionInRead0 << " A: " << totalNumberOfA << " C: " << totalNumberOfC << " G: " << totalNumberOfG << " T: " << totalNumberOfT << " Gap: " << totalNumberOfGap << " Total Alignments that support a change: " << totalNumberOfAlignmentsWithoutGapsSupportingAChangeFromTheTargetRead << ". Total alignments: " << totalCoverageInThisPosition << ". Total alignments without gaps: " << totalCoverageInThisPositionWithoutGaps << "." << endl;
                    continue;
                }


                // Now find the orientedReadIds that support the target base by:
                // 1. Taking all covering reads
                // 2. Subtracting reads that support any other alternative bases (C, G, T, Gap)

                std::set<OrientedReadId> orientedReadIdsSupportingTheTargetReadBase;

                // Start with all covering reads
                orientedReadIdsSupportingTheTargetReadBase = coveringOrientedReadIds;

                // Remove reads supporting the specific alternative base
                std::set<OrientedReadId> tempResult;
                std::set_difference(
                    orientedReadIdsSupportingTheTargetReadBase.begin(), orientedReadIdsSupportingTheTargetReadBase.end(),
                    orientedReadIdsSupportingAChangeFromTheTargetRead.begin(), orientedReadIdsSupportingAChangeFromTheTargetRead.end(),
                    std::inserter(tempResult, tempResult.begin())
                );
                orientedReadIdsSupportingTheTargetReadBase.swap(tempResult);

                // Add the target read to the set of reads supporting the target read base
                orientedReadIdsSupportingTheTargetReadBase.insert(OrientedReadId(targetReadId, targetReadStrand));



                // XXX
                // --- Check for strand bias ---
                // 

                double st_rate = 0.05;
                uint64_t st_max = 2;

                // Find how many of the hetBase1OrientedReadIds are on the same strand as the target read (0)
                uint64_t orientedReadIdsSupportingTheTargetReadBaseOnStrand0 = 0;
                for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheTargetReadBase) {
                    if (orientedReadId.getStrand() == 0) {
                        orientedReadIdsSupportingTheTargetReadBaseOnStrand0++;
                    }
                }

                if (((orientedReadIdsSupportingTheTargetReadBaseOnStrand0 + st_max) >= totalNumberOfAlignmentsSupportingTheTargetRead) && ((totalNumberOfAlignmentsSupportingTheTargetRead*(st_rate) + orientedReadIdsSupportingTheTargetReadBaseOnStrand0) >= totalNumberOfAlignmentsSupportingTheTargetRead) ) {
                    // debugOut << "Strand bias detected at position " << positionStatsOnTargetOrientedRead.positionInTargetRead << endl;
                    // debugOut << "Total number of alignments supporting the target read on strand 0: " << orientedReadIdsSupportingTheTargetReadBaseOnStrand0 << endl;
                    // debugOut << "Total number of alignments supporting the target read: " << totalNumberOfAlignmentsSupportingTheTargetRead << endl;
                    continue;
                }

                // XXX
                // --- END OF: Check for strand bias ---
                //



                if (baseAHasEnoughCoverage) {

                    // Populate the SnpStats object
                    SnpStats hetPositionStatsInTargetRead;
                    hetPositionStatsInTargetRead.positionInTargetRead = positionStatsOnTargetOrientedRead.positionInTargetRead;
                    hetPositionStatsInTargetRead.hetBase0 = positionStatsOnTargetOrientedRead.baseOfTargetRead;
                    hetPositionStatsInTargetRead.hetBase1 = AlignedBase::fromCharacter('A');
                    hetPositionStatsInTargetRead.isInHomopolymerRegion = positionStatsOnTargetOrientedRead.isInHomopolymerRegion;
                    hetPositionStatsInTargetRead.occ_0 = totalNumberOfAlignmentsSupportingTheTargetRead; // Occurrence of the target base
                    hetPositionStatsInTargetRead.occ_1 = totalNumberOfA; // Occurrence of the alternative base
                    hetPositionStatsInTargetRead.score = -1;
                    hetPositionStatsInTargetRead.overlap_num = totalCoverageInThisPosition;


                    const double supportingTheTargetReadBaseA = double(totalNumberOfAlignmentsSupportingTheTargetRead) / double(totalCoverageInThisPosition);
                    const double supportingTheHetReadBaseA = double(totalNumberOfA) / double(totalCoverageInThisPosition);


                    // Get the haplotypeEvidence for the orientedReadIds supporting the alternative base
                    for (const HaplotypeEvidence& haplotypeEvidence : positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[AlignedBase::fromCharacter('A').value]) {
                        hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE.push_back(haplotypeEvidence);
                    }

                    for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheTargetReadBase) {
                        hetPositionStatsInTargetRead.hetBase0OrientedReadIds.push_back(orientedReadId);
                    }


                    // XXX
                    // --- START OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //

                    uint64_t occ_1_temp = hetPositionStatsInTargetRead.occ_1;

                    // Loop over the haplotype evidence and check if any of the reads are in a homopolymer region
                    for (const HaplotypeEvidence& hapEvidence : hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE) {
                        
                        if (hapEvidence.isInHomopolymerRegion) {
                            occ_1_temp--;
                            continue;
                        }

                    }

                    if (occ_1_temp < 3) {
                        hetPositionStatsInTargetRead.isHpcVector = true;
                    } 

                    // XXX
                    // --- END OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //




                    // XXX
                    // --- START OF: Check if the position has sufficent coverage ---
                    //

                    // TODO: WE NEED ADJUSTABLE HOM_COV
                    double cut_rate = 0.7;
                    uint64_t cut_bd = 6;
                    int64_t hom_cov = 44;
                    int64_t n_hap = 2;

                    uint64_t cc = (hom_cov/n_hap);
                    cc *= cut_rate;
                    if (cc < cut_bd) {
                        cc = cut_bd;
                    }

                    if (hetPositionStatsInTargetRead.occ_0 >= cc) {
                        hetPositionStatsInTargetRead.hasSufficentCoverage = true;
                    }

                    // XXX
                    // --- END OF: Check if the position has sufficent coverage ---
                    //



                    // This Snp position passed all the filters
                    // We can now populate the potentialHetSitesOnOrientedReadId0
                    potentialHetSitesOnTargetOrientedReadId.push_back(hetPositionStatsInTargetRead);

                    
                }


                if (baseCHasEnoughCoverage) {

                    // Populate the SnpStats object
                    SnpStats hetPositionStatsInTargetRead;
                    hetPositionStatsInTargetRead.positionInTargetRead = positionStatsOnTargetOrientedRead.positionInTargetRead;
                    hetPositionStatsInTargetRead.hetBase0 = positionStatsOnTargetOrientedRead.baseOfTargetRead;
                    hetPositionStatsInTargetRead.hetBase1 = AlignedBase::fromCharacter('C');
                    hetPositionStatsInTargetRead.isInHomopolymerRegion = positionStatsOnTargetOrientedRead.isInHomopolymerRegion;
                    hetPositionStatsInTargetRead.occ_0 = totalNumberOfAlignmentsSupportingTheTargetRead; // Occurrence of the target base
                    hetPositionStatsInTargetRead.occ_1 = totalNumberOfC; // Occurrence of the alternative base
                    hetPositionStatsInTargetRead.score = -1;
                    hetPositionStatsInTargetRead.overlap_num = totalCoverageInThisPosition;


                    const double supportingTheTargetReadBaseC = double(totalNumberOfAlignmentsSupportingTheTargetRead) / double(totalCoverageInThisPosition);
                    const double supportingTheHetReadBaseC = double(totalNumberOfC) / double(totalCoverageInThisPosition);


                    // Get the haplotypeEvidence for the orientedReadIds supporting the alternative base
                    for (const HaplotypeEvidence& haplotypeEvidence : positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[AlignedBase::fromCharacter('C').value]) {
                        hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE.push_back(haplotypeEvidence);
                    }

                    for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheTargetReadBase) {
                        hetPositionStatsInTargetRead.hetBase0OrientedReadIds.push_back(orientedReadId);
                    }

                    // XXX
                    // --- START OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //

                    uint64_t occ_1_temp = hetPositionStatsInTargetRead.occ_1;

                    // Loop over the haplotype evidence and check if any of the reads are in a homopolymer region
                    for (const HaplotypeEvidence& hapEvidence : hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE) {
                        
                        if (hapEvidence.isInHomopolymerRegion) {
                            occ_1_temp--;
                            continue;
                        }

                    }

                    if (occ_1_temp < 3) {
                        hetPositionStatsInTargetRead.isHpcVector = true;
                    } 

                    // XXX
                    // --- END OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //




                    // XXX
                    // --- START OF: Check if the position has sufficent coverage ---
                    //

                    double cut_rate = 0.7;
                    uint64_t cut_bd = 6;
                    int64_t hom_cov = 44;
                    int64_t n_hap = 2;

                    uint64_t cc = (hom_cov/n_hap);
                    cc *= cut_rate;
                    if (cc < cut_bd) {
                        cc = cut_bd;
                    }

                    if (hetPositionStatsInTargetRead.occ_0 >= cc) {
                        hetPositionStatsInTargetRead.hasSufficentCoverage = true;
                    }

                    // XXX
                    // --- END OF: Check if the position has sufficent coverage ---
                    //



                    // This Snp position passed all the filters
                    // We can now populate the potentialHetSitesOnOrientedReadId0
                    potentialHetSitesOnTargetOrientedReadId.push_back(hetPositionStatsInTargetRead);

                    
                }


                if (baseGHasEnoughCoverage) {

                    // Populate the SnpStats object
                    SnpStats hetPositionStatsInTargetRead;
                    hetPositionStatsInTargetRead.positionInTargetRead = positionStatsOnTargetOrientedRead.positionInTargetRead;
                    hetPositionStatsInTargetRead.hetBase0 = positionStatsOnTargetOrientedRead.baseOfTargetRead;
                    hetPositionStatsInTargetRead.hetBase1 = AlignedBase::fromCharacter('G');
                    hetPositionStatsInTargetRead.isInHomopolymerRegion = positionStatsOnTargetOrientedRead.isInHomopolymerRegion;
                    hetPositionStatsInTargetRead.occ_0 = totalNumberOfAlignmentsSupportingTheTargetRead; // Occurrence of the target base
                    hetPositionStatsInTargetRead.occ_1 = totalNumberOfG; // Occurrence of the alternative base
                    hetPositionStatsInTargetRead.score = -1;
                    hetPositionStatsInTargetRead.overlap_num = totalCoverageInThisPosition;


                    const double supportingTheTargetReadBaseG = double(totalNumberOfAlignmentsSupportingTheTargetRead) / double(totalCoverageInThisPosition);
                    const double supportingTheHetReadBaseG = double(totalNumberOfG) / double(totalCoverageInThisPosition);


                    // Get the haplotypeEvidence for the orientedReadIds supporting the alternative base
                    for (const HaplotypeEvidence& haplotypeEvidence : positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[AlignedBase::fromCharacter('G').value]) {
                        hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE.push_back(haplotypeEvidence);
                    }

                    for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheTargetReadBase) {
                        hetPositionStatsInTargetRead.hetBase0OrientedReadIds.push_back(orientedReadId);
                    }

                    // XXX
                    // --- START OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //

                    uint64_t occ_1_temp = hetPositionStatsInTargetRead.occ_1;

                    // Loop over the haplotype evidence and check if any of the reads are in a homopolymer region
                    for (const HaplotypeEvidence& hapEvidence : hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE) {
                        
                        if (hapEvidence.isInHomopolymerRegion) {
                            occ_1_temp--;
                            continue;
                        }

                    }

                    if (occ_1_temp < 3) {
                        hetPositionStatsInTargetRead.isHpcVector = true;
                    } 

                    // XXX
                    // --- END OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //




                    // XXX
                    // --- START OF: Check if the position has sufficent coverage ---
                    //

                    double cut_rate = 0.7;
                    uint64_t cut_bd = 6;
                    int64_t hom_cov = 44;
                    int64_t n_hap = 2;

                    uint64_t cc = (hom_cov/n_hap);
                    cc *= cut_rate;
                    if (cc < cut_bd) {
                        cc = cut_bd;
                    }

                    if (hetPositionStatsInTargetRead.occ_0 >= cc) {
                        hetPositionStatsInTargetRead.hasSufficentCoverage = true;
                    }

                    // XXX
                    // --- END OF: Check if the position has sufficent coverage ---
                    //



                    // This Snp position passed all the filters
                    // We can now populate the potentialHetSitesOnOrientedReadId0
                    potentialHetSitesOnTargetOrientedReadId.push_back(hetPositionStatsInTargetRead);

                    
                }



                if (baseTHasEnoughCoverage) {

                    // Populate the SnpStats object
                    SnpStats hetPositionStatsInTargetRead;
                    hetPositionStatsInTargetRead.positionInTargetRead = positionStatsOnTargetOrientedRead.positionInTargetRead;
                    hetPositionStatsInTargetRead.hetBase0 = positionStatsOnTargetOrientedRead.baseOfTargetRead;
                    hetPositionStatsInTargetRead.hetBase1 = AlignedBase::fromCharacter('T');
                    hetPositionStatsInTargetRead.isInHomopolymerRegion = positionStatsOnTargetOrientedRead.isInHomopolymerRegion;
                    hetPositionStatsInTargetRead.occ_0 = totalNumberOfAlignmentsSupportingTheTargetRead; // Occurrence of the target base
                    hetPositionStatsInTargetRead.occ_1 = totalNumberOfT; // Occurrence of the alternative base
                    hetPositionStatsInTargetRead.score = -1;
                    hetPositionStatsInTargetRead.overlap_num = totalCoverageInThisPosition;


                    const double supportingTheTargetReadBaseT = double(totalNumberOfAlignmentsSupportingTheTargetRead) / double(totalCoverageInThisPosition);
                    const double supportingTheHetReadBaseT = double(totalNumberOfT) / double(totalCoverageInThisPosition);


                    // Get the haplotypeEvidence for the orientedReadIds supporting the alternative base
                    for (const HaplotypeEvidence& haplotypeEvidence : positionStatsOnTargetOrientedRead.haplotypeEvidenceByBase[AlignedBase::fromCharacter('T').value]) {
                        hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE.push_back(haplotypeEvidence);
                    }

                    for (const OrientedReadId& orientedReadId : orientedReadIdsSupportingTheTargetReadBase) {
                        hetPositionStatsInTargetRead.hetBase0OrientedReadIds.push_back(orientedReadId);
                    }

                    // XXX
                    // --- START OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //

                    uint64_t occ_1_temp = hetPositionStatsInTargetRead.occ_1;

                    // Loop over the haplotype evidence and check if any of the reads are in a homopolymer region
                    for (const HaplotypeEvidence& hapEvidence : hetPositionStatsInTargetRead.hetBase1OrientedReadIdsHE) {
                        
                        if (hapEvidence.isInHomopolymerRegion) {
                            occ_1_temp--;
                            continue;
                        }

                    }

                    if (occ_1_temp < 3 || hetPositionStatsInTargetRead.isInHomopolymerRegion) {
                        hetPositionStatsInTargetRead.isHpcVector = true;
                    } 

                    // XXX
                    // --- END OF: Check if the haplotypeEvidence in this position suggests a HPC region ---
                    //




                    // XXX
                    // --- START OF: Check if the position has sufficent coverage ---
                    //

                    double cut_rate = 0.7;
                    uint64_t cut_bd = 6;
                    int64_t hom_cov = 44;
                    int64_t n_hap = 2;

                    uint64_t cc = (hom_cov/n_hap);
                    cc *= cut_rate;
                    if (cc < cut_bd) {
                        cc = cut_bd;
                    }

                    cc = 6;

                    if (hetPositionStatsInTargetRead.occ_0 >= cc) {
                        hetPositionStatsInTargetRead.hasSufficentCoverage = true;
                    }

                    // XXX
                    // --- END OF: Check if the position has sufficent coverage ---
                    //



                    // This Snp position passed all the filters
                    // We can now populate the potentialHetSitesOnOrientedReadId0
                    potentialHetSitesOnTargetOrientedReadId.push_back(hetPositionStatsInTargetRead);

                    
                }




            } // --- End of loop over all potential heterozygous sites ---




            debugOut << "Thread ID: " << threadId << " Finished analyzing potential heterozygous sites for read ID: " << targetReadId << endl;

            // Clear the positionStatsOnOrientedReadId0 map
            // positionStatsOnTargetOrientedReadId.clear();




            // XXX
            // --- Now we need to analyze the potential heterozygous sites and try to find sets of reads that belong to the same haplotype ---
            //

            if (potentialHetSitesOnTargetOrientedReadId.size() > 0) {
                debugOut << endl;
                debugOut << "Found " << potentialHetSitesOnTargetOrientedReadId.size() << " potential heterozygous sites in readId " << targetReadId << endl;
                debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
            } else {
                debugOut << endl;
                debugOut << "Found no potential heterozygous sites in readId " << targetReadId << endl;
                debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
                continue;
            }

            // Loop over the potential heterozygous sites and print them
            char bases[] = {'A', 'C', 'G', 'T', '-'};
            for (const auto& positionStatsInTargetRead : potentialHetSitesOnTargetOrientedReadId) {

                debugOut << "Potential heterozygous site at position " << positionStatsInTargetRead.positionInTargetRead
                    << " in readId " << targetReadId << " in strand " << targetReadStrand << ":" << endl
                    << "  Base of targetReadId: " << bases[positionStatsInTargetRead.hetBase0.value] << endl
                    << "  First variant: " << bases[positionStatsInTargetRead.hetBase0.value] << " (" << positionStatsInTargetRead.occ_0 << " "<< positionStatsInTargetRead.occ_0 / positionStatsInTargetRead.overlap_num * 100 << "%)" << endl
                    << "  Second variant: " << bases[positionStatsInTargetRead.hetBase1.value] << " (" << positionStatsInTargetRead.occ_1 << " " << positionStatsInTargetRead.occ_1 / positionStatsInTargetRead.overlap_num * 100 << "%)" << endl;
                
                debugOut << "ReadIds that support the first variant: " << endl;
                for (const auto& orientedReadId : positionStatsInTargetRead.hetBase0OrientedReadIds) {
                    debugOut << "ReadId: " << orientedReadId.getReadId() << " Strand: " << orientedReadId.getStrand() << endl;
                }
                debugOut << "ReadIds that support the second variant: " << endl;
                for (const auto& haplotypeEvidence : positionStatsInTargetRead.hetBase1OrientedReadIdsHE) {
                    debugOut << "ReadId: " << haplotypeEvidence.orientedReadId.getReadId() << " Strand: " << haplotypeEvidence.orientedReadId.getStrand() << endl;
                }
                
            }

            debugOut << "Thread ID: " << threadId << " Finished analyzing potential heterozygous sites for read ID: " << targetReadId << endl;

            

            // XXX
            // --- Dynamic Programming for Grouping Compatible Sites and Phasing ---
            //
            
            // debugOut << timestamp << "Starting DP for grouping compatible sites for read " << targetOrientedReadId << endl;

            uint64_t N = potentialHetSitesOnTargetOrientedReadId.size();
            if (N == 0) {
                debugOut << "No potential het sites found after filtering for read " << targetOrientedReadId << ", skipping DP and subsequent phasing." << endl;
                continue; // Skip to the next readId if no sites
            } else {
                debugOut << "Found " << N << " potential heterozygous sites for DP in readId " << targetReadId << endl;
                // debugOut << "Skipped " << sitesSkippedDueToInsufficientCoverage << " sites due to insufficient coverage" << endl;
            }

            // Define compatibility check function locally (or move to class scope)
            auto areCompatibleSites =
                [&](const SnpStats& stats_i,
                    const SnpStats& stats_j,
                    const OrientedReadId& targetOrientedReadId) -> bool
            {   

                // Never link a SNP to itself
                if (stats_i.positionInTargetRead == stats_j.positionInTargetRead) {
                    return false;
                }


                // 1a. Get all reads overlapping site i
                std::vector<OrientedReadId> reads_at_i;
                std::vector<OrientedReadId> reads_at_i_0;
                std::vector<OrientedReadId> reads_at_i_1;
                reads_at_i.reserve(stats_i.hetBase0OrientedReadIds.size() + stats_i.hetBase1OrientedReadIdsHE.size());
                reads_at_i.insert(reads_at_i.end(), stats_i.hetBase0OrientedReadIds.begin(), stats_i.hetBase0OrientedReadIds.end());
                reads_at_i_0.insert(reads_at_i_0.end(), stats_i.hetBase0OrientedReadIds.begin(), stats_i.hetBase0OrientedReadIds.end());
                for (const HaplotypeEvidence& haplotypeEvidence : stats_i.hetBase1OrientedReadIdsHE) {
                    reads_at_i.push_back(haplotypeEvidence.orientedReadId);
                    reads_at_i_1.push_back(haplotypeEvidence.orientedReadId);
                }
                std::sort(reads_at_i.begin(), reads_at_i.end());

                // 1b. Get all reads overlapping site j
                std::vector<OrientedReadId> reads_at_j;
                std::vector<OrientedReadId> reads_at_j_0;
                std::vector<OrientedReadId> reads_at_j_1;
                reads_at_j.reserve(stats_j.hetBase0OrientedReadIds.size() + stats_j.hetBase1OrientedReadIdsHE.size());
                reads_at_j.insert(reads_at_j.end(), stats_j.hetBase0OrientedReadIds.begin(), stats_j.hetBase0OrientedReadIds.end());
                reads_at_j_0.insert(reads_at_j_0.end(), stats_j.hetBase0OrientedReadIds.begin(), stats_j.hetBase0OrientedReadIds.end());
                for (const HaplotypeEvidence& haplotypeEvidence : stats_j.hetBase1OrientedReadIdsHE) {
                    reads_at_j.push_back(haplotypeEvidence.orientedReadId);
                    reads_at_j_1.push_back(haplotypeEvidence.orientedReadId);
                }
                std::sort(reads_at_j.begin(), reads_at_j.end());


                // 2. Find common overlapping reads
                std::vector<OrientedReadId> commonOverlappingReads;
                std::set_intersection(reads_at_i.begin(), reads_at_i.end(),
                                    reads_at_j.begin(), reads_at_j.end(),
                                    std::back_inserter(commonOverlappingReads));


                bool found_phase_00 = false; // Flag for reads matching target at both sites (Phase 0-0)
                bool found_phase_11 = false; // Flag for reads differing from target at both sites (Phase 1-1)
                uint64_t found_phase_00_count = 0; // Count of reads matching target at both sites
                uint64_t found_phase_11_count = 0; // Count of reads differing from target at both sites


                // 3. Check phase consistency for each common read relative to target 
                //    and also evidence for both phase patterns
                for (const auto& orientedReadId : commonOverlappingReads) {
                    // Determine phase at site i (0: matches target, 1: differs from target, -1: unknown/not informative)
                    int phase_i = -1;
                    AlignedBase allele_at_i; // Invalid allele
                    if (std::find(reads_at_i_0.begin(), reads_at_i_0.end(), orientedReadId) != reads_at_i_0.end()) {
                        allele_at_i = stats_i.hetBase0;
                    } else if (std::find(reads_at_i_1.begin(), reads_at_i_1.end(), orientedReadId) != reads_at_i_1.end()) {
                        allele_at_i = stats_i.hetBase1;
                    }

                    if (allele_at_i == stats_i.hetBase0) {
                        phase_i = 0; // Matches the base of the target read
                    } else {
                        phase_i = 1; // Differs from the base of the target read
                    }

                    // debugOut << "Read " << read << " at site " << site_pair_i.first << " has phase " << phase_i << endl;

                    // Determine phase at site j (0: matches target, 1: differs from target, -1: unknown/not informative)
                    int phase_j = -1;
                    AlignedBase allele_at_j;
                    if (std::find(reads_at_j_0.begin(), reads_at_j_0.end(), orientedReadId) != reads_at_j_0.end()) {
                        allele_at_j = stats_j.hetBase0;
                    } else if (std::find(reads_at_j_1.begin(), reads_at_j_1.end(), orientedReadId) != reads_at_j_1.end()) {
                        allele_at_j = stats_j.hetBase1;
                    }

                    if (allele_at_j == stats_j.hetBase0) {
                        phase_j = 0; // Matches target
                    } else {
                        phase_j = 1; // Differs from target
                    } 

                    // debugOut << "Read " << read << " at site " << site_pair_j.first << " has phase " << phase_j << endl;

                    // Check for inconsistency or update flags if consistent and informative
                    if (phase_i != -1 && phase_j != -1) { // Both phases defined relative to target for this read
                        if (phase_i != phase_j) {
                            // Found a read with inconsistent phasing relative to the target across the two sites.
                            // debugOut << "Incompatibility: Read " << read << " phase " << phase_i << " at site " << site_pair_i.first << ", phase " << phase_j << " at site " << site_pair_j.first << " relative to target." << endl;
                            return false; // Inconsistent phasing detected
                        } else if (phase_i == 0) { // phase_i == 0 && phase_j == 0
                            found_phase_00 = true; // Found evidence for 0-0 linkage relative to target
                            found_phase_00_count++;
                        } else { // phase_i == 1 && phase_j == 1
                            found_phase_11 = true; // Found evidence for 1-1 linkage relative to target
                            found_phase_11_count++;
                        }
                    } else if (phase_i == -1 || phase_j == -1) {
                        // debugOut << "Read " << read << " is not informative for phasing relative to target at site " << site_pair_i.first << " or site " << site_pair_j.first << "." << endl;
                        return false; // Inconsistent phasing detected
                    }

                }

                // Sites are compatible *relative to the target* only if:
                // 1. No common reads showed inconsistent phasing between positions i and j (phase_i != phase_j).
                // 2. There was at least one read supporting the 0-0 linkage relative to target.
                // 3. There was at least one read supporting the 1-1 linkage relative to target.
                // cout << "Compatibility check for sites " << site_pair_i.first << " and " << site_pair_j.first << " relative to target " << targetReadId << ": found_phase_00=" << found_phase_00 << ", found_phase_11=" << found_phase_11 << endl;
                return found_phase_00 && found_phase_11;

            };


            // 1. Initialize DP table and parent pointers
            std::vector<uint64_t> LCG(N, 1); // LCG[i] = size of largest compatible group ending at i
            std::vector<int64_t> parent(N, -1); // parent[i] = index j that gave the max LCG[i]

            // 2. Fill DP table using recurrence relation: LCG(i) = max_{j<i, S[j]<->S[i]} {LCG(j)} + 1
            for (uint64_t i = 1; i < N; ++i) {
                for (uint64_t j = 0; j < i; ++j) {
                    // Check compatibility S[j] <-> S[i]
                    if (areCompatibleSites(potentialHetSitesOnTargetOrientedReadId[i], potentialHetSitesOnTargetOrientedReadId[j], targetOrientedReadId)) {
                        if (LCG[j] + 1 > LCG[i]) {
                            LCG[i] = LCG[j] + 1;
                            parent[i] = j;
                        }
                    }
                }
            }

            // --- Traceback and Phasing ---
            std::vector<bool> isAssigned(N, false);

            // Create pairs of (LCG value, index) to sort by LCG value descending
            std::vector<std::pair<uint64_t, size_t>> sortedLCGIndices;
            std::vector<std::pair<uint64_t, size_t>> isolatedSitesIndices;
            for(size_t i = 0; i < N; ++i) {
                if (LCG[i] > 1) { // Groupped sites (they have LCG > 1)
                    sortedLCGIndices.push_back({LCG[i], i});
                } else if (LCG[i] == 1) { // Isolated sites
                    isolatedSitesIndices.push_back({LCG[i], i});
                }
            }
            // Sort descending by LCG value
            std::sort(sortedLCGIndices.rbegin(), sortedLCGIndices.rend()); 

    

            // Data structure for storing the chains of compatible sites.
            // Each inner vector contains the indices of SnpStats that form a compatible chain.
            std::vector<std::vector<uint64_t>> siteChains;
            std::vector<std::vector<uint64_t>> siteChainsIsolatedSites;

            // --- Traceback to build the chains ---
            // First, process the multi-site chains from the DP results.
            for (const auto& lcgPair : sortedLCGIndices) {
                size_t endOfChainIndex = lcgPair.second;

                if (!isAssigned[endOfChainIndex]) {
                    std::vector<uint64_t> newChain;
                    int traceIndex = endOfChainIndex;

                    // Trace back from the end of the chain to the beginning using the parent pointers.
                    while (traceIndex != -1 && !isAssigned[traceIndex]) {
                        newChain.push_back(traceIndex);
                        isAssigned[traceIndex] = true;
                        traceIndex = parent[traceIndex];
                    }

                    // The chain was built backwards, so reverse it to store it in positional order.
                    std::reverse(newChain.begin(), newChain.end());
                    siteChains.push_back(newChain);
                }
            }

            // Now, handle the isolated sites (those not part of any larger chain) as chains of size 1.
            for (const auto& isolatedPair : isolatedSitesIndices) {
                uint64_t siteIndex = isolatedPair.second;
                if (!isAssigned[siteIndex]) {
                    siteChainsIsolatedSites.push_back({siteIndex}); // A chain containing a single site.
                }
            }



            // XXX
            // --- START OF: Chain Scoring ---
            //

            // Long chain are inherently more reliable because they have more linked SNPs.
            // So initialy all SNP positions in the chain get a score of 1.
            for (const auto& chain : siteChains) {
                for (const uint64_t siteIndex : chain) {
                    if(potentialHetSitesOnTargetOrientedReadId[siteIndex].hasSufficentCoverage){
                        potentialHetSitesOnTargetOrientedReadId[siteIndex].score = 1;
                    } else {
                        potentialHetSitesOnTargetOrientedReadId[siteIndex].score = -1;
                    }
                }
            }

            for (const auto& chain : siteChainsIsolatedSites) {
                for (const uint64_t siteIndex : chain) {
                    if(potentialHetSitesOnTargetOrientedReadId[siteIndex].hasSufficentCoverage && !potentialHetSitesOnTargetOrientedReadId[siteIndex].isHpcVector){
                        potentialHetSitesOnTargetOrientedReadId[siteIndex].score = 1;
                    } else {
                        potentialHetSitesOnTargetOrientedReadId[siteIndex].score = -1;
                    }
                }
            }

            // XXX
            // --- END OF: Chain Scoring ---
            //




            // Print the chains keeping only the sites that have a score of 1
            debugOut << "DP results for read " << targetOrientedReadId << ": LCG values > 1:" << endl;
            for(const auto& chain : siteChains) {
                bool foundGoodChain = false;
                for(const auto& siteIndex : chain) {
                    if(potentialHetSitesOnTargetOrientedReadId[siteIndex].score == 1){
                        foundGoodChain = true;
                        debugOut << "SCORE 1:  Site Index: " << siteIndex << " (Pos: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].positionInTargetRead << ") LCG: " << LCG[siteIndex] << endl;
                    }
                }
                if (foundGoodChain){
                    debugOut << endl;
                }
            }

            for(const auto& chain : siteChains) {
                bool foundGoodChain = false;
                for(const auto& siteIndex : chain) {
                    debugOut << "SCORE ALL:  Site Index: " << siteIndex << " (Pos: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].positionInTargetRead << ") LCG: " << LCG[siteIndex] << endl;
                }
                
            }

            debugOut << "DP results for read " << targetOrientedReadId << ": LCG values == 1:" << endl;
            for(const auto& chain : siteChainsIsolatedSites) {
                bool foundGoodChain = false;
                for(const auto& siteIndex : chain) {
                    if(potentialHetSitesOnTargetOrientedReadId[siteIndex].score == 1){
                        // if (potentialHetSitesOnTargetOrientedReadId[siteIndex].positionInTargetRead == 22129) {
                        //     debugOut << "  Site Index: " << siteIndex << " hasSufficentCoverage: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].hasSufficentCoverage << " isHpcVector: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].isHpcVector << " occ_0: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].occ_0 << " occ_1: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].occ_1 << endl;
                        // }
                        foundGoodChain = true;
                        debugOut << "SCORE 1:  Site Index: " << siteIndex << " (Pos: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].positionInTargetRead << ") LCG: " << LCG[siteIndex] << endl;
                    }
                }
                if (foundGoodChain){
                    debugOut << endl;
                }
            }

            for(const auto& chain : siteChainsIsolatedSites) {
                for(const auto& siteIndex : chain) {
                    debugOut << "SCORE ALL ISOLATED:  Site Index: " << siteIndex << " (Pos: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].positionInTargetRead << ") LCG: " << LCG[siteIndex] << endl;
                    debugOut << "Has good coverage: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].hasSufficentCoverage << endl;
                    debugOut << "Is HPC vector: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].isHpcVector << endl;
                    debugOut << "Occ_0: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].occ_0 << endl;
                    debugOut << "Occ_1: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].occ_1 << endl;
                    debugOut << "Score: " << potentialHetSitesOnTargetOrientedReadId[siteIndex].score << endl;
                }
                
            }


            // Loop over all the sites of the chains and if they have a score of 1
            // add them to the filtered vector
            for (const auto& chain : siteChains) {
                for (const uint64_t siteIndex : chain) {
                    if (potentialHetSitesOnTargetOrientedReadId[siteIndex].score == 1) {
                        potentialHetSitesOnTargetOrientedReadIdFiltered.push_back(potentialHetSitesOnTargetOrientedReadId[siteIndex]);
                    }
                }
            }

            for (const auto& chain : siteChainsIsolatedSites) {
                for (const uint64_t siteIndex : chain) {
                    if (potentialHetSitesOnTargetOrientedReadId[siteIndex].score == 1) {
                        potentialHetSitesOnTargetOrientedReadIdFiltered.push_back(potentialHetSitesOnTargetOrientedReadId[siteIndex]);
                    }
                }
            }

            // Sort the filtered vector by position in target read
            std::sort(potentialHetSitesOnTargetOrientedReadIdFiltered.begin(), potentialHetSitesOnTargetOrientedReadIdFiltered.end());


            debugOut << "Filtered vector size: " << potentialHetSitesOnTargetOrientedReadIdFiltered.size() << endl;

            // Print the filtered vector
            for (const auto& site : potentialHetSitesOnTargetOrientedReadIdFiltered) {
                debugOut << "  Site Index: " << site.positionInTargetRead << " (Pos: " << site.positionInTargetRead << ")" << endl;
            }

            N = potentialHetSitesOnTargetOrientedReadIdFiltered.size();
            if (N == 0) {
                debugOut << "No potential het sites found after filtering for read " << targetOrientedReadId << ", skipping DP and subsequent phasing." << endl;
                continue; // Skip to the next readId if no sites
            }









            // // XXX
            // // --- START OF: NAIVE PHASING ---
            // // generate_haplotypes_naive_HiFi


            // // Filter out streaks of adjacent site positions.
            // std::vector<SnpStats> potentialHetSitesOnTargetOrientedReadIdFilteredNoAdjacentSnps;
            // if (!potentialHetSitesOnTargetOrientedReadIdFiltered.empty()) {
            //     for (uint64_t i = 0; i < potentialHetSitesOnTargetOrientedReadIdFiltered.size(); ) {
            //         // Start a new streak with the index of the current site.
            //         std::vector<uint64_t> currentStreakIndices;
            //         currentStreakIndices.push_back(i);

            //         // Look ahead to find all adjacent sites for the current streak.
            //         uint64_t j = i + 1;
            //         for (; j < potentialHetSitesOnTargetOrientedReadIdFiltered.size(); ++j) {
            //             // A site is adjacent if its position is the same as the previous
            //             // or exactly one greater.
            //             if ((potentialHetSitesOnTargetOrientedReadIdFiltered[j].positionInTargetRead == potentialHetSitesOnTargetOrientedReadIdFiltered[j-1].positionInTargetRead) || 
            //                 (potentialHetSitesOnTargetOrientedReadIdFiltered[j].positionInTargetRead == potentialHetSitesOnTargetOrientedReadIdFiltered[j-1].positionInTargetRead + 1)) {
            //                 currentStreakIndices.push_back(j);
            //             } else {
            //                 // The streak is broken.
            //                 break;
            //             }
            //         }
                    
            //         // Add the completed streak to the list of all streaks, but only if it contains multiple sites.
            //         if (currentStreakIndices.size() == 1) {
            //             for (const uint64_t siteStreakIndex : currentStreakIndices) {
            //                 potentialHetSitesOnTargetOrientedReadIdFilteredNoAdjacentSnps.push_back(potentialHetSitesOnTargetOrientedReadIdFiltered[siteStreakIndex]);
            //             }
            //         }
                    
            //         // Move the main index to the start of the next potential streak.
            //         i = j;
            //     }
            // }

            // // Print the size of the filtered vector
            // debugOut << "Filtered vector size before filtering out streaks: " << potentialHetSitesOnTargetOrientedReadIdFiltered.size() << endl;
            // debugOut << "Filtered vector size after filtering out streaks: " << potentialHetSitesOnTargetOrientedReadIdFilteredNoAdjacentSnps.size() << endl;

            // // Print the filtered vector
            // for (const auto& site : potentialHetSitesOnTargetOrientedReadIdFilteredNoAdjacentSnps) {
            //     debugOut << "  Site Index: " << site.positionInTargetRead << " (Pos: " << site.positionInTargetRead << ")" << endl;
            // }

            // if (potentialHetSitesOnTargetOrientedReadIdFilteredNoAdjacentSnps.size() == 0) {
            //     debugOut << "No potential het sites found after filtering for read " << targetOrientedReadId << ", skipping DP and subsequent phasing." << endl;
            //     continue; // Skip to the next readId if no sites
            // }


            // // XXX
            // // --- END OF: NAIVE PHASING ---
            // //

            


            // XXX
            // --- START OF: Filter out clusters of SNPs ---
            //
            // Sequencing error and artifacts often appear as clusters of nearby SNPs.
            // To avoid clusters of errors, the informative SNPs need to be well-separated.
            // Only SNPs at least 32bp apart are considered.
            // E.g. If we have SNPs on positions: [9992, 33364, 33368, 40000]
            // We want to keep the SNPs on positions: [9992, 40000]
            std::vector<SnpStats> potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated;
            const auto& inputSnps = potentialHetSitesOnTargetOrientedReadIdFiltered;

            // If there are 0 or 1 SNPs, they are all "well-separated" by definition and should be kept.
            if (inputSnps.size() < 2) {
                potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated = inputSnps;
            } else {
                for (size_t i = 0; i < inputSnps.size(); ++i) {
                    const bool hasPrev = (i > 0);
                    const bool hasNext = (i < inputSnps.size() - 1);

                    // A SNP is "too close" if the distance is > 0 and < 32.
                    // This correctly preserves multi-allelic sites (distance == 0), which should not be filtered.
                    const bool prevIsTooClose = hasPrev &&
                        (inputSnps[i].positionInTargetRead > inputSnps[i-1].positionInTargetRead) &&
                        (inputSnps[i].positionInTargetRead - inputSnps[i-1].positionInTargetRead < 32);

                    // Check distance to next SNP.
                    const bool nextIsTooClose = hasNext &&
                        (inputSnps[i+1].positionInTargetRead > inputSnps[i].positionInTargetRead) &&
                        (inputSnps[i+1].positionInTargetRead - inputSnps[i].positionInTargetRead < 32);

                    // Only keep the SNP if it is not too close to an adjacent, different SNP.
                    if (!prevIsTooClose && !nextIsTooClose) {
                        potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated.push_back(inputSnps[i]);
                    }
                }
            }


            // print the size of the filtered vector and the sites
            debugOut << "Filtered vector size after filtering out well separated SNPs: " << potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated.size() << endl;
            for (const auto& site : potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated) {
                debugOut << "  Site Index: " << site.positionInTargetRead << " (Pos: " << site.positionInTargetRead << ")" << endl;
            }


            // Only if at least 2 filtered and well separated het SNPs are found, continue with the phasing.
            if (potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated.size() < 2) {
                debugOut << "Less than 2 filtered and well separated het SNPs found for read " << targetOrientedReadId << ", skipping DP and subsequent phasing." << endl;
                continue;
            }

            // XXX
            // --- END OF: Filtering out clusters of SNPs ---
            //



            

            



            // XXX
            // --- START OF: Phasing based on the filtered, well-separated SNPs ---
            //

            // For clarity, create a reference to the final list of SNPs to be used for phasing.
            const auto& informativeSnps = potentialHetSitesOnTargetOrientedReadIdFilteredWellSeparated;
            // const auto& informativeSnps = potentialHetSitesOnTargetOrientedReadIdFiltered;

            // 1. Collect all in-phase and out-of-phase reads from all informative sites.
            std::vector<OrientedReadId> inPhaseReads;
            std::vector<OrientedReadId> outOfPhaseReads;

            for (const auto& snp : informativeSnps) {
                // In-phase reads are those that match the reference allele (hetBase0).
                inPhaseReads.insert(inPhaseReads.end(), snp.hetBase0OrientedReadIds.begin(), snp.hetBase0OrientedReadIds.end());

                // Out-of-phase reads are those that match the alternate allele (hetBase1).
                for (const auto& evidence : snp.hetBase1OrientedReadIdsHE) {
                    outOfPhaseReads.push_back(evidence.orientedReadId);
                }
            }

            // 2. Sort and remove duplicates to create unique sets for comparison.
            std::sort(inPhaseReads.begin(), inPhaseReads.end());
            inPhaseReads.erase(std::unique(inPhaseReads.begin(), inPhaseReads.end()), inPhaseReads.end());

            std::sort(outOfPhaseReads.begin(), outOfPhaseReads.end());
            outOfPhaseReads.erase(std::unique(outOfPhaseReads.begin(), outOfPhaseReads.end()), outOfPhaseReads.end());

            // 3. Find reads that are in the in-phase set but NOT in the out-of-phase set.
            // A read is considered consistently in-phase only if it never supports an out-of-phase allele.
            std::vector<OrientedReadId> finalInPhaseReads;
            std::set_difference(
                inPhaseReads.begin(), inPhaseReads.end(),
                outOfPhaseReads.begin(), outOfPhaseReads.end(),
                std::back_inserter(finalInPhaseReads)
            );

            // 4. The target read is always considered in-phase with itself.
            finalInPhaseReads.push_back(targetOrientedReadId);
            // Sort and unique again to handle the case where the target read was already present.
            std::sort(finalInPhaseReads.begin(), finalInPhaseReads.end());
            finalInPhaseReads.erase(std::unique(finalInPhaseReads.begin(), finalInPhaseReads.end()), finalInPhaseReads.end());

            // Debug print the final set of consistently in-phase reads.
            debugOut << "Final in-phase reads for read " << targetOrientedReadId << ": " << endl;
            for (const auto& read : finalInPhaseReads) {
                debugOut << "  ReadId: " << read.getReadId() << " Strand: " << read.getStrand() << endl;
            }

            // XXX
            // --- END OF: Phasing based on the filtered, well-separated SNPs ---
            //





            debugOut << "Thread ID: " << threadId << " Finished grouping compatible sites for read ID: " << targetReadId << endl;

            // XXX
            // --- Keep the correct alignments and forbid the bad alignments based on the compatible site phasing ---
            //

            uint64_t numberOfFirstPassHetAlignments = 0;
            // --- Update thread-local boolean vectors ---
            // Loop over all alignments and mark those involving only reads within finalInPhaseOrientedReads
            // Also forbid alignments involving reads in excludedOutOfPhaseOrientedReads

            // Loop over alignments involving the target read (targetReadId)
            for (const auto alignmentId : alignmentTableForTargetOrientedReadId) {

                const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
                ReadId currentReadId0 = thiAlignmentData.readIds[0];
                ReadId currentReadId1 = thiAlignmentData.readIds[1];
                OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thiAlignmentData.info;

                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != targetReadId) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetReadId);

                // Flip strands, if necessary
                if (currentOrientedReadId0.getStrand() != targetReadStrand) {
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == targetReadStrand);

                // Check if the alignment involves reads that are in-phase or out-of-phase.
                // The vectors are sorted, so we use std::binary_search for efficient lookups.
                const bool read0IsInPhase = std::binary_search(finalInPhaseReads.begin(), finalInPhaseReads.end(), currentOrientedReadId0);
                const bool read1IsInPhase = std::binary_search(finalInPhaseReads.begin(), finalInPhaseReads.end(), currentOrientedReadId1);
                const bool involvesfinalInPhaseOrientedReads = read0IsInPhase && read1IsInPhase;
                
                if (involvesfinalInPhaseOrientedReads) {
                    if (phasingThreadDataReadGraph5->isReadIdContained[currentOrientedReadId0.getReadId()]) {
                        continue;
                    }
                    // Mark firstPassHet in thread-local vector
                    // debugOut << "Marking alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " as first pass (intra-phase)" << endl;
                    firstPassHetAlignments[alignmentId] = true;
                    numberOfFirstPassHetAlignments++;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                } else {
                    if (phasingThreadDataReadGraph5->isReadIdContained[currentOrientedReadId0.getReadId()]) {
                        continue;
                    }
                    // Mark forbidden in thread-local vector
                    // debugOut << "Forbidding alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
                    forbiddenAlignments[alignmentId] = true;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

            }

           
            






            // // // We do not perform the analysis on contained reads because the longer reads will
            // // // make the decision for them (where they best belong) given that they will have more informative
            // // // sites to consider. However, these contained reads should be assigned to one cluster of reads only.
            // // // Otherwise, we observe between-haplotypes alignments that connect to
            // // // these contained reads from both sides of a het region.
            // // // So, we also need to forbid alignments between contained reads and excludedOutOfPhaseOrientedReads.
            // // for (const auto& orientedRead : finalInPhaseOrientedReads) {
            // //     if (phasingThreadDataReadGraph5->isReadIdContained[orientedRead.getReadId()]) {
            // //         const auto alignmentTableForContainedRead = alignmentTable[orientedRead.getValue()];
            // //         // Loop over alignments involving contained read.
            // //         for (const auto alignmentId : alignmentTableForContainedRead) {
            // //
            // //             const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            // //             ReadId currentReadId0 = thiAlignmentData.readIds[0];
            // //             ReadId currentReadId1 = thiAlignmentData.readIds[1];
            // //             OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            // //             OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            // //
            // //             // Swap oriented reads, if necessary.
            // //             if (currentOrientedReadId0.getReadId() != orientedRead.getReadId()) {
            // //                 swap(currentOrientedReadId0, currentOrientedReadId1);
            // //             }
            // //             SHASTA_ASSERT(currentOrientedReadId0.getReadId() == orientedRead.getReadId());
            // //
            // //             // Flip strands, if necessary
            // //             if (currentOrientedReadId0.getStrand() != orientedRead.getStrand()) {
            // //                 currentOrientedReadId0.flipStrand();
            // //                 currentOrientedReadId1.flipStrand();
            // //             }
            // //             SHASTA_ASSERT(currentOrientedReadId0.getStrand() == orientedRead.getStrand());
            // //
            // //             // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
            // //             bool involvesExcludedOrientedReads = (finalInPhaseOrientedReads.count(currentOrientedReadId0) && excludedOutOfPhaseOrientedReads.count(currentOrientedReadId1));
            // //
            // //             if (involvesExcludedOrientedReads) {
            // //                 debugOut << "Forbidding alignment " << alignmentId << " involving contained read: " << currentReadId0 << " and read: " << currentReadId1 << " because the other read is in excludedOutOfPhaseOrientedReads." << endl;
            // //                 forbiddenAlignments[alignmentId] = true;
            // //                 alignmentsAlreadyConsidered[alignmentId] = true;
            // //             }
            // //
            // //         }
            // //
            // //     }
            // // }

            // // cout << "Number of first pass het alignments for readId " << targetReadId << ": " << numberOfFirstPassHetAlignments << endl;

            // // debugOut << timestamp << "Finished DP for grouping compatible sites for read " << targetOrientedReadId << endl;

            // // XXX
            // // --- END OF: Keep the correct alignments and forbid the bad alignments based on the compatible site phasing ---
            // //


            // debugOut << "Thread ID: " << threadId << " Finished processing read ID: " << targetReadId << endl;

            // // XXX
            // // --- Populate the sites vector for this read's haplotype block ---
            // //

            // // Create a new Site object for the identified orientedReadId cluster of the target read
            // Site newSite;
            // newSite.orientedReads.insert(finalInPhaseOrientedReads.begin(), finalInPhaseOrientedReads.end());
            // newSite.targetOrientedReadId = targetOrientedReadId;
            // newSite.excludedOrientedReads.insert(excludedOutOfPhaseOrientedReads.begin(), excludedOutOfPhaseOrientedReads.end());

            // if (!newSite.orientedReads.empty()) {
            //     sites.push_back(newSite);
            // }

            // // Create a site of the reverse strand
            // Site newSiteReverseStrand;
            // for (const auto& orientedRead : finalInPhaseOrientedReads) {
            //     OrientedReadId reverseOrientedRead(orientedRead.getReadId(), orientedRead.getStrand() == 0 ? 1 : 0);
            //     newSiteReverseStrand.orientedReads.insert(reverseOrientedRead);
            // }
            // for (const auto& orientedRead : excludedOutOfPhaseOrientedReads) {
            //     OrientedReadId reverseOrientedRead(orientedRead.getReadId(), orientedRead.getStrand() == 0 ? 1 : 0);
            //     newSiteReverseStrand.excludedOrientedReads.insert(reverseOrientedRead);
            // }
            // newSiteReverseStrand.targetOrientedReadId = OrientedReadId(targetOrientedReadId.getReadId(), targetOrientedReadId.getStrand() == 0 ? 1 : 0);

            // if (!newSiteReverseStrand.orientedReads.empty()) {
            //     sites.push_back(newSiteReverseStrand);
            // }












            // // XXX
            // // --- END OF: Populate the sites vector for this read's haplotype block ---
            // //


            
            // debugOut << "Thread ID: " << threadId << " Finished creating site for read ID: " << targetReadId << endl;

            // // Second loop over alignments to make decisions based on the now-complete site info for this read.
            // for (const auto alignmentId : alignmentTable0) {

            //     if (forbiddenAlignments[alignmentId]) {
            //         continue;
            //     }

            //     const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            //     ReadId currentReadId0 = thiAlignmentData.readIds[0];
            //     ReadId currentReadId1 = thiAlignmentData.readIds[1];
            //     OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            //     OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            //     AlignmentInfo alignmentInfo = thiAlignmentData.info;

            //     // Swap oriented reads, if necessary.
            //     if (currentOrientedReadId0.getReadId() != targetReadId) {
            //         swap(currentOrientedReadId0, currentOrientedReadId1);
            //         alignmentInfo.swap();
            //     }
            //     SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetReadId);

            //     // Flip strands, if necessary
            //     if (currentOrientedReadId0.getStrand() != strand0) {
            //         currentOrientedReadId0.flipStrand();
            //         currentOrientedReadId1.flipStrand();
            //         alignmentInfo.reverseComplement();
            //     }
            //     SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

                

            //     // Get the alignment range on read0.
            //     const auto& markers0 = markers[currentOrientedReadId0.getValue()];
            //     const uint32_t start0 = markers0[alignmentInfo.data[0].firstOrdinal].position;
            //     const uint32_t end0   = markers0[alignmentInfo.data[0].lastOrdinal].position;

            //     uint64_t sharedHet = 0;
            //     uint64_t uniqueHet0 = 0;
            //     uint64_t sharedMatch = 0;

            //     for (const auto& [pos0, stats0] : potentialHetSitesOnOrientedReadId0) {
            //         if (pos0 >= start0 && pos0 <= end0) {
            //             // This site is in the alignment range of read 0.
            //             bool isInHet1 = stats0.hetBase1OrientedReadIds.count(currentOrientedReadId1);
            //             bool isInHet2 = stats0.hetBase2OrientedReadIds.count(currentOrientedReadId1);

            //             if (isInHet1 || isInHet2) {
            //                 sharedHet++;
            //                 // This is a shared site. Check if the alleles match.
            //                 uint64_t alleleOnRead0 = stats0.baseOfReadId0;
            //                 uint64_t alleleOnRead1 = isInHet1 ? stats0.hetBase1 : stats0.hetBase2;
            //                 if (alleleOnRead0 == alleleOnRead1) {
            //                     sharedMatch++;
            //                 }
            //             } else {
            //                 uniqueHet0++;
            //             }
            //         }
            //     }

            //     // We loop over the sites of currentOrientedReadId0 and check in each site
            //     // how many sites are shared and matched with currentOrientedReadId1 and
            //     // how many are shared and not match with currentOrientedReadId1.
            //     // We also count how many are totally shared.
            //     const uint64_t sharedMismatch = sharedHet - sharedMatch;
            //     const uint64_t uniqueHet1 = 0;


            //     if (sharedHet > 0 || uniqueHet0 > 0 || uniqueHet1 > 0) {
            //         debugOut << "Alignment " << alignmentId << " between reads " << currentReadId0 << " and " << currentReadId1
            //                  << " from perspective of " << targetOrientedReadId
            //                  << ": shared=" << sharedHet
            //                  << ", shared_match=" << sharedMatch
            //                  << ", shared_mismatch=" << sharedMismatch
            //                  << ", unique_to_0=" << uniqueHet0
            //                  << ", unique_to_1=" << uniqueHet1 << endl;
            //     }
            // }

            // debugOut << "Thread ID: " << threadId << " Finished second loop over alignments for read ID: " << targetReadId << endl;

            // // --- END OF: Second loop over alignments to make decisions based on the now-complete site info for this read. ---



            // // XXX
            // // --- Target read overlaps filtering ---
            // //

            // uint64_t numberOfAlignmentsToKeep = 0;

            // for (const auto alignmentId : alignmentTable0) {

            //     if (forbiddenAlignments[alignmentId]) {
            //         continue;
            //     }

            //     const AlignmentData& thisAlignmentData = alignmentData[alignmentId];
            //     ReadId currentReadId0 = thisAlignmentData.readIds[0];
            //     ReadId currentReadId1 = thisAlignmentData.readIds[1];
            //     OrientedReadId currentOrientedReadId0(thisAlignmentData.readIds[0], 0);
            //     OrientedReadId currentOrientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
            //     AlignmentInfo alignmentInfo = thisAlignmentData.info;

            //     // Swap oriented reads, if necessary.
            //     if (currentOrientedReadId0.getReadId() != targetReadId) {
            //         swap(currentOrientedReadId0, currentOrientedReadId1);
            //         alignmentInfo.swap();
            //     }
            //     SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetReadId);

            //     // Flip strands, if necessary
            //     if (currentOrientedReadId0.getStrand() != strand0) {
            //         currentOrientedReadId0.flipStrand();
            //         currentOrientedReadId1.flipStrand();
            //         alignmentInfo.reverseComplement();
            //     }
            //     SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

            //     // // TODO: Skip reads that have been marked as forbidden
            //     // if (forbiddenReads[currentReadId0] || forbiddenReads[currentReadId1]) {
            //     //     continue;
            //     // }

            //     // Get the base ranges of the oriented reads.
            //     const uint64_t baseRange0 = alignmentInfo.baseRange(assemblerInfo->k, currentOrientedReadId0, 0, markers);
            //     const uint64_t baseRange1 = alignmentInfo.baseRange(assemblerInfo->k, currentOrientedReadId1, 1, markers);

            //     // Find the alignment range of the oriented reads.
            //     const uint64_t range0 = alignmentInfo.data[0].range();
            //     const uint64_t range1 = alignmentInfo.data[1].range();

            //     // Find the leftTrim and rightTrim of the target read in this alignment.
            //     const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            //     const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();

            //     // Find the leftTrim and rightTrim of the query read in this alignment.
            //     const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            //     const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            //     const uint64_t minimumOfLeftTrims = std::min(leftTrim0, leftTrim1);
            //     const uint64_t minimumOfRightTrims = std::min(rightTrim0, rightTrim1);

            //     // The legth of the aligned segment on the query read must be at least 80% 
            //     // of the total span of the alignment (aligned segment + rightTrim + leftTrim).
            //     // If the aligned part is too small compared to the overhangs, it's likely
            //     // a spurious or low-quality alignment.
            //     const double minimumAlignedFraction = 0.8;

            //     const uint64_t max_hang = 1000;

                
            //     if (minimumOfLeftTrims > max_hang || minimumOfRightTrims > max_hang) {
            //         // If either left or right trim is too large, skip this alignment.
            //         // This is to avoid very long overhangs that might indicate a poor alignment.
            //         // debugOut << "Skipping alignment " << alignmentId << " due to excessive left or right trim." << endl;
            //         forbiddenAlignments[alignmentId] = true;
            //         continue;

            //     }

            //     // The legth of the aligned segment on the target read must be at least [minimumAlignedFraction]% 
            //     // of the total span of the alignment (aligned range + rightTrim + leftTrim).
            //     if (range0 < minimumAlignedFraction * (range0 + minimumOfLeftTrims + minimumOfRightTrims)) {
            //         // If the aligned part is too small compared to the overhangs, it's likely
            //         // a spurious or low-quality alignment.
            //         forbiddenAlignments[alignmentId] = true;
            //         continue;

            //     }

            //     // The legth of the aligned segment on the query read must be at least [minimumAlignedFraction]% 
            //     // of the total span of the alignment (aligned range + rightTrim + leftTrim).
            //     if (range1 < minimumAlignedFraction * (range1 + minimumOfLeftTrims + minimumOfRightTrims)) {
            //         // If the aligned part is too small compared to the overhangs, it's likely
            //         // a spurious or low-quality alignment.
            //         forbiddenAlignments[alignmentId] = true;
            //         continue;

            //     }


            //     // // Classify alignment
            //     // if (alignmentInfo.classify(20) == AlignmentType::read0IsContained) {
            //     //     // 0 is unambiguously contained in 1.
            //     //     forbiddenAlignments[alignmentId] = true;
            //     //     continue;
            //     // }

            //     // if (alignmentInfo.classify(20) == AlignmentType::read1IsContained) {
            //     //     // 1 is unambiguously contained in 0.
            //     //     forbiddenAlignments[alignmentId] = true;
            //     //     continue;
            //     // }




            // }




        } // End loop over read IDs in batch

    } // End while getNextBatch

}





void Assembler::createReadGraph5(
    uint64_t maxAlignmentCount,
    double epsilon,
    double delta,
    double WThreshold,
    double WThresholdForBreaks
    )
{
    cout << timestamp << "createReadGraph5 with strand separation begins" << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Get stats about the reads
    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2*readCount;

    // Initialize result vectors
    vector<bool> forbiddenAlignments(alignmentCount, false);
    vector<bool> firstPassHetAlignments(alignmentCount, false);
    vector<bool> alignmentsAlreadyConsidered(alignmentCount, false);
    vector<Site> sites;

    // --- Parallel Phasing Section ---
    cout << timestamp << "Starting parallel phasing analysis." << endl;
    uint64_t threadCount = assemblerInfo->threadCount;
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Initialize shared data structure for threads
    phasingThreadDataReadGraph5 = new PhasingThreadDataReadGraph5(this, threadCount);

    phasingThreadDataReadGraph5->isReadIdContained.resize(readCount, false);

    // Loop over read IDs to find contained reads
    for (ReadId readId = 0; readId < readCount; readId++) {

        const ReadId readId0 = readId;
        const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
        const OrientedReadId targetOrientedReadId(readId0, strand0);

        // Loop over alignments involving this read
        const auto alignmentTable0 = alignmentTable[targetOrientedReadId.getValue()];
        for (const auto alignmentId : alignmentTable0) {

            const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            ReadId currentReadId0 = thiAlignmentData.readIds[0];
            ReadId currentReadId1 = thiAlignmentData.readIds[1];
            OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = thiAlignmentData.info;

            // Swap oriented reads, if necessary.
            if (currentOrientedReadId0.getReadId() != readId0) {
                swap(currentOrientedReadId0, currentOrientedReadId1);
                alignmentInfo.swap();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

            // Flip strands, if necessary
            if (currentOrientedReadId0.getStrand() != strand0) {
                currentOrientedReadId0.flipStrand();
                currentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

            const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();

            const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            // TODO: this is wrong
            // const uint64_t maxTrim = assemblerInfo->actualMaxTrim;
            // assemblerOptions.alignOptions.maxTrim
            const uint64_t maxTrim = 50;

            // Check for containment.
            const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
            const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);

            if (isContained0 && !isContained1) {
                // 0 is unambiguously contained in 1.
                phasingThreadDataReadGraph5->isReadIdContained[currentOrientedReadId0.getReadId()] = true;
            }
            // if(isContained1 && !isContained0) {
            //     // 1 is unambiguously contained in 0.
            //     phasingThreadDataReadGraph5->isReadIdContained[currentReadId1] = true;
            // }
            // if(isContained0 && isContained1) {
            //     // Near complete overlap.
            //     continue;
            // }

        }
    }

    // Print the contained reads
    uint64_t containedReads = 0;
    cout << "Contained reads: " << endl;
    for (ReadId readId = 0; readId < readCount; readId++) {
        if (phasingThreadDataReadGraph5->isReadIdContained[readId]) {
            cout << "  ReadId: " << readId << endl;
            containedReads++;
        }
    }
    cout << "Finished printing contained reads." << endl;
    cout << "Total number of contained reads is: " << containedReads << endl;

    // Copy the isReadIdContained from threads
    vector<bool> isReadIdContained;
    isReadIdContained.resize(readCount);
    for (ReadId readId = 0; readId < readCount; readId++) {
        isReadIdContained[readId] = phasingThreadDataReadGraph5->isReadIdContained[readId];
    }


    // Initialize batch processing for reads
    const uint64_t requestedBatchSize = 1; // Or some other suitable value
    setupLoadBalancing(readCount, requestedBatchSize);


    runThreads(&Assembler::createReadGraph5ThreadFunction, threadCount);


    cout << timestamp << "Parallel phasing analysis completed." << endl;

    // Aggregate results from threads
    cout << timestamp << "Aggregating results from threads." << endl;
    for(uint64_t threadId=0; threadId<threadCount; ++threadId) {
        for(uint64_t alignmentId=0; alignmentId<alignmentCount; ++alignmentId) {
            if(phasingThreadDataReadGraph5->threadForbiddenAlignments[threadId][alignmentId]) {
                forbiddenAlignments[alignmentId] = true;
            }
            if(phasingThreadDataReadGraph5->threadAlignmentsAlreadyConsidered[threadId][alignmentId]) {
                alignmentsAlreadyConsidered[alignmentId] = true;
            }
            if(phasingThreadDataReadGraph5->threadFirstPassHetAlignments[threadId][alignmentId]) {
                firstPassHetAlignments[alignmentId] = true;
            }
        }
        
        // Aggregate the sites from each thread
        for(const Site& site : phasingThreadDataReadGraph5->threadSites[threadId]) {
            sites.push_back(site);
        }

    }



    // Clean up shared data
    delete phasingThreadDataReadGraph5;
    phasingThreadDataReadGraph5 = nullptr;
    cout << timestamp << "Finished aggregating results." << endl;
    // --- End Parallel Phasing Section ---

    // Count results
    const long forbiddenCount = count(forbiddenAlignments.begin(), forbiddenAlignments.end(), true);
    const long firstPassCount = count(firstPassHetAlignments.begin(), firstPassHetAlignments.end(), true);
    const uint64_t sitesCount = sites.size();
    cout << timestamp << "Total forbidden alignments after phasing: " << forbiddenCount << endl;
    cout << timestamp << "Total first pass het alignments after phasing: " << firstPassCount << endl;
    cout << timestamp << "Total sites after phasing: " << sitesCount << endl;




    vector <bool> isReadIdContained2;
    isReadIdContained2.resize(readCount, false);
    // Loop over read IDs to find contained reads
    for (ReadId readId = 0; readId < readCount; readId++) {

        const ReadId readId0 = readId;
        const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
        const OrientedReadId targetOrientedReadId(readId0, strand0);

        // Loop over alignments involving this read
        const auto alignmentTable0 = alignmentTable[targetOrientedReadId.getValue()];
        for (const auto alignmentId : alignmentTable0) {

            if (readId0 != 3743) {
                continue;
            }

            if (!firstPassHetAlignments[alignmentId]) {
                continue;
            }
            
            if (forbiddenAlignments[alignmentId]) {
                continue;
            }

            const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            ReadId currentReadId0 = thiAlignmentData.readIds[0];
            ReadId currentReadId1 = thiAlignmentData.readIds[1];
            OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = thiAlignmentData.info;

            // Swap oriented reads, if necessary.
            if (currentOrientedReadId0.getReadId() != readId0) {
                swap(currentOrientedReadId0, currentOrientedReadId1);
                alignmentInfo.swap();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

            // Flip strands, if necessary
            if (currentOrientedReadId0.getStrand() != strand0) {
                currentOrientedReadId0.flipStrand();
                currentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

            const uint32_t startingAlignmentBasePosition0 = alignmentInfo.startingAlignmentBasePosition[0];
            const uint32_t startingAlignmentBasePosition1 = alignmentInfo.startingAlignmentBasePosition[1];
            const uint32_t endingAlignmentBasePosition0 = alignmentInfo.endingAlignmentBasePosition[0];
            const uint32_t endingAlignmentBasePosition1 = alignmentInfo.endingAlignmentBasePosition[1];
            const uint32_t leftTrimBases0 = alignmentInfo.leftTrimBases[0];
            const uint32_t leftTrimBases1 = alignmentInfo.leftTrimBases[1];
            const uint32_t rightTrimBases0 = alignmentInfo.rightTrimBases[0];
            const uint32_t rightTrimBases1 = alignmentInfo.rightTrimBases[1];

            cout << "AlignmentId: " << alignmentId << endl;
            cout << "startingAlignmentBasePosition0: " << startingAlignmentBasePosition0 << endl;
            cout << "startingAlignmentBasePosition1: " << startingAlignmentBasePosition1 << endl;
            cout << "endingAlignmentBasePosition0: " << endingAlignmentBasePosition0 << endl;
            cout << "endingAlignmentBasePosition1: " << endingAlignmentBasePosition1 << endl;
            cout << "leftTrimBases0: " << leftTrimBases0 << endl;
            cout << "leftTrimBases1: " << leftTrimBases1 << endl;
            cout << "rightTrimBases0: " << rightTrimBases0 << endl;
            cout << "rightTrimBases1: " << rightTrimBases1 << endl;

            if (leftTrimBases0 > leftTrimBases1 && rightTrimBases0 > rightTrimBases1) {
                cout << "ReadId: " << currentOrientedReadId1.getReadId() << " is contained in readId: " << currentOrientedReadId0.getReadId() << endl;
            }

            if (leftTrimBases0 < leftTrimBases1 && rightTrimBases0 < rightTrimBases1) {
                cout << "ReadId: " << currentOrientedReadId0.getReadId() << " is contained in readId: " << currentOrientedReadId1.getReadId() << endl;
            }



            const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();

            const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            // TODO: this is wrong
            // const uint64_t maxTrim = assemblerInfo->actualMaxTrim;
            // assemblerOptions.alignOptions.maxTrim
            const uint64_t maxTrim = 50;

            // Check for containment.
            const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
            const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);

            if (readId0 == 3743) {
                cout << "isContained0: " << isContained0 << endl;
                cout << "isContained1: " << isContained1 << endl;
                cout << "leftTrim0: " << leftTrim0 << endl;
                cout << "rightTrim0: " << rightTrim0 << endl;
                cout << "leftTrim1: " << leftTrim1 << endl;
                cout << "rightTrim1: " << rightTrim1 << endl;
            }

            if (isContained0 && !isContained1) {
                // 0 is unambiguously contained in 1.
                isReadIdContained2[currentOrientedReadId0.getReadId()] = true;
            }
            // if(isContained1 && !isContained0) {
            //     // 1 is unambiguously contained in 0.
            //     isReadIdContained2[currentReadId1] = true;
            // }
            // if(isContained0 && isContained1) {
            //     // Near complete overlap.
            //     continue;
            // }

        }
    }

    // Print the contained reads
    uint64_t newContainedReads = 0;
    cout << "Contained reads: " << endl;
    for (ReadId readId = 0; readId < readCount; readId++) {
        if (isReadIdContained2[readId]) {
            cout << "  ReadId: " << readId << endl;
            newContainedReads++;
        }
    }
    cout << "Finished printing contained reads." << endl;
    cout << "New total number of contained reads is: " << newContainedReads << endl;
    

    
    // Loop over all sites and check if a read is included in both strands in the site
    vector<bool> sitesThatHaveStrandIssuesForbidden(sites.size(), false);
    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        const Site& site = sites[siteId];
        std::set<ReadId> readIdsInSite;
        for(const OrientedReadId orientedReadId: site.orientedReads) {
            readIdsInSite.insert(orientedReadId.getReadId());
        }
        // Check if both strands of the same read are present in the site
        for(const ReadId readId: readIdsInSite) {
            OrientedReadId orientedReadId0(readId, 0);
            OrientedReadId orientedReadId1(readId, 1);
            if (site.orientedReads.count(orientedReadId0) && site.orientedReads.count(orientedReadId1)) {
                cout << timestamp << "Found both strands of read " << readId << " in site " << siteId << "." << endl;
                // print all reads in the site
                cout << "Reads in site " << siteId << ": ";
                for(const OrientedReadId orientedReadId: site.orientedReads) {
                    cout << orientedReadId.getReadId() << " ";
                }
                cout << endl;
                // Mark the site as having strand issues
                sitesThatHaveStrandIssuesForbidden[siteId] = true;
            }
        }
    }


    // Loop over all sites and check if a read is included in both strands in the site
    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        const Site& site = sites[siteId];

        if (site.targetOrientedReadId.getReadId() != 1055) {
            continue;
        }


        for(const OrientedReadId orientedReadId: site.orientedReads) {
            cout << "Site of " << site.targetOrientedReadId.getReadId() << " has oriented read: " << orientedReadId.getReadId() << endl;
        }

        for(const OrientedReadId orientedReadId: site.excludedOrientedReads) {
            cout << "Site of " << site.targetOrientedReadId.getReadId() << " has excluded this oriented read: " << orientedReadId.getReadId() << endl;
        }

    }
    
    
    // We need to find, for each orientedReadId, the sites that are associated with it
    // and fill in the orientedReadSites
    vector< vector<uint64_t> > orientedReadSites(readCount * 2); // Indexed via OrientedReadId::getValue()
    
    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        // if (sitesThatHaveStrandIssuesForbidden[siteId]) {
        //     continue; // Skip sites that have strand issues
        // }
        const Site& site = sites[siteId];
        OrientedReadId targetOrientedReadId = site.targetOrientedReadId;
        // for(const OrientedReadId orientedReadId: site.orientedReads) {
        //     orientedReadSites[orientedReadId.getValue()].push_back(siteId);
        // }
        orientedReadSites[targetOrientedReadId.getValue()].push_back(siteId);
    }
    
    // for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
    //     // if (sitesThatHaveStrandIssuesForbidden[siteId]) {
    //     //     continue; // Skip sites that have strand issues
    //     // }
    //     const Site& site = sites[siteId];
    //     for(const OrientedReadId orientedReadId: site.orientedReads) {
    //         orientedReadSites[orientedReadId.getValue()].push_back(siteId);
    //     }
    // }
    
    cout << timestamp << "Finished finding for each orientedReadId, the sites that are associated with it." << endl;



    for(uint64_t siteId=0; siteId<sites.size(); siteId++) {
        const Site& site = sites[siteId];
        OrientedReadId targetOrientedReadId = site.targetOrientedReadId;
        if (isReadIdContained[targetOrientedReadId.getReadId()] && !isReadIdContained2[targetOrientedReadId.getReadId()]) {
            // cout << "Trying to analyze readId: " << targetOrientedReadId.getReadId() << endl;
            const auto alignmentTable0 = alignmentTable[targetOrientedReadId.getValue()];
            for (const auto alignmentId : alignmentTable0) {

                const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
                ReadId currentReadId0 = thiAlignmentData.readIds[0];
                ReadId currentReadId1 = thiAlignmentData.readIds[1];
                OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
                OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
                AlignmentInfo alignmentInfo = thiAlignmentData.info;

                // Swap oriented reads, if necessary.
                if (currentOrientedReadId0.getReadId() != targetOrientedReadId.getReadId()) {
                    swap(currentOrientedReadId0, currentOrientedReadId1);
                    alignmentInfo.swap();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getReadId() == targetOrientedReadId.getReadId());

                // Flip strands, if necessary
                if (currentOrientedReadId0.getStrand() != targetOrientedReadId.getStrand()) {
                    currentOrientedReadId0.flipStrand();
                    currentOrientedReadId1.flipStrand();
                    alignmentInfo.reverseComplement();
                }
                SHASTA_ASSERT(currentOrientedReadId0.getStrand() == targetOrientedReadId.getStrand());

                // Check if either read in the alignment is in excludedOutOfPhaseOrientedReads
                bool involvesfinalInPhaseOrientedReads = (site.orientedReads.count(currentOrientedReadId0) && site.orientedReads.count(currentOrientedReadId1));
                bool involvesExcludedOrientedReads = (site.excludedOrientedReads.count(currentOrientedReadId0) && site.excludedOrientedReads.count(currentOrientedReadId1));
                
                if (involvesfinalInPhaseOrientedReads) {
                    //cout << "Added firstPassHetAlignments alignmentId: " << alignmentId << endl;
                    // if (phasingThreadDataReadGraph5->isReadIdContained[currentOrientedReadId0.getReadId()]) {
                    //     continue;
                    // }
                    // Mark firstPassHet in thread-local vector
                    // debugOut << "Marking alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " as first pass (intra-phase)" << endl;
                    firstPassHetAlignments[alignmentId] = true;
                    // numberOfFirstPassHetAlignments++;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

                if (involvesExcludedOrientedReads) {
                    //cout << "Added forbiddenAlignments alignmentId: " << alignmentId << endl;
                    // Mark forbidden in thread-local vector
                    // debugOut << "Forbidding alignment " << alignmentId << " involving reads: " << currentReadId0 << " and " << currentReadId1 << " because one read is in excludedOutOfPhaseOrientedReads." << endl;
                    forbiddenAlignments[alignmentId] = true;
                    alignmentsAlreadyConsidered[alignmentId] = true;
                }

            }
        }
        
    }













    vector <bool> isReadIdContained3;
    isReadIdContained3.resize(readCount, false);
    // Loop over read IDs to find contained reads
    for (ReadId readId = 0; readId < readCount; readId++) {

        const ReadId readId0 = readId;
        const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
        const OrientedReadId targetOrientedReadId(readId0, strand0);

        // Loop over alignments involving this read
        const auto alignmentTable0 = alignmentTable[targetOrientedReadId.getValue()];
        for (const auto alignmentId : alignmentTable0) {

            const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            ReadId currentReadId0 = thiAlignmentData.readIds[0];
            ReadId currentReadId1 = thiAlignmentData.readIds[1];
            OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = thiAlignmentData.info;

            // Swap oriented reads, if necessary.
            if (currentOrientedReadId0.getReadId() != readId0) {
                swap(currentOrientedReadId0, currentOrientedReadId1);
                alignmentInfo.swap();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

            // Flip strands, if necessary
            if (currentOrientedReadId0.getStrand() != strand0) {
                currentOrientedReadId0.flipStrand();
                currentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);


            if (!firstPassHetAlignments[alignmentId]) {
                continue;
            }
            
            if (forbiddenAlignments[alignmentId]) {
                continue;
            }

            if (currentOrientedReadId0.getReadId() == 3737) {
                cout << "AlignmentId: " << alignmentId << endl;
                cout << "FirstPassHetAlignments: " << firstPassHetAlignments[alignmentId] << endl;
                cout << "ForbiddenAlignments: " << forbiddenAlignments[alignmentId] << endl;
            }

            const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();

            const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            // TODO: this is wrong
            // const uint64_t maxTrim = assemblerInfo->actualMaxTrim;
            // assemblerOptions.alignOptions.maxTrim
            const uint64_t maxTrim = 50;

            // Check for containment.
            const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
            const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);

            if (currentOrientedReadId0.getReadId() == 3737) {
                cout << "isContained0: " << isContained0 << endl;
                cout << "isContained1: " << isContained1 << endl;
                cout << "leftTrim0: " << leftTrim0 << endl;
                cout << "rightTrim0: " << rightTrim0 << endl;
                cout << "leftTrim1: " << leftTrim1 << endl;
                cout << "rightTrim1: " << rightTrim1 << endl;
            }

            if (isContained0 && !isContained1) {
                // 0 is unambiguously contained in 1.
                isReadIdContained3[currentOrientedReadId0.getReadId()] = true;
            }
            // if(isContained1 && !isContained0) {
            //     // 1 is unambiguously contained in 0.
            //     isReadIdContained3[currentReadId1] = true;
            // }
            // if(isContained0 && isContained1) {
            //     // Near complete overlap.
            //     continue;
            // }

        }
    }

    // Print the contained reads
    uint64_t finalNewContainedReads = 0;
    cout << "Contained reads: " << endl;
    for (ReadId readId = 0; readId < readCount; readId++) {
        if (isReadIdContained3[readId]) {
            cout << "  ReadId: " << readId << endl;
            finalNewContainedReads++;
        }
    }
    cout << "Finished printing contained reads." << endl;
    cout << "Final new total number of contained reads is: " << finalNewContainedReads << endl;





    




    // Loop over read IDs to find contained reads
    for (ReadId readId = 0; readId < readCount; readId++) {

        const ReadId readId0 = readId;
        const Strand strand0 = 0; // Analyze strand 0 arbitrarily, results should be consistent
        const OrientedReadId targetOrientedReadId(readId0, strand0);

        if (!isReadIdContained3[readId0]) {
            continue;
        }

        vector <uint32_t> readsContainingReadId0;

        // Loop over alignments involving this read
        const auto alignmentTable0 = alignmentTable[targetOrientedReadId.getValue()];
        for (const auto alignmentId : alignmentTable0) {

            if (!firstPassHetAlignments[alignmentId]) {
                continue;
            }
            
            if (forbiddenAlignments[alignmentId]) {
                continue;
            }

            const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
            ReadId currentReadId0 = thiAlignmentData.readIds[0];
            ReadId currentReadId1 = thiAlignmentData.readIds[1];
            OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = thiAlignmentData.info;

            // Swap oriented reads, if necessary.
            if (currentOrientedReadId0.getReadId() != readId0) {
                swap(currentOrientedReadId0, currentOrientedReadId1);
                alignmentInfo.swap();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getReadId() == readId0);

            // Flip strands, if necessary
            if (currentOrientedReadId0.getStrand() != strand0) {
                currentOrientedReadId0.flipStrand();
                currentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(currentOrientedReadId0.getStrand() == strand0);

            // If the other read is also contained, forbid the alignment.
            if (isReadIdContained3[currentOrientedReadId1.getReadId()]) {
                forbiddenAlignments[alignmentId] = true;
                continue;
            }

            // if (currentOrientedReadId0.getReadId() == 330 || currentOrientedReadId0.getReadId() == 333) {
                
            // }

            // // Check if the other reads is the one that contains the current read.
            // const uint64_t leftTrim0 = alignmentInfo.data[0].leftTrim();
            // const uint64_t rightTrim0 = alignmentInfo.data[0].rightTrim();
            // const uint64_t leftTrim1 = alignmentInfo.data[1].leftTrim();
            // const uint64_t rightTrim1 = alignmentInfo.data[1].rightTrim();

            // // TODO: this is wrong
            // // const uint64_t maxTrim = assemblerInfo->actualMaxTrim;
            // // assemblerOptions.alignOptions.maxTrim
            // const uint64_t maxTrim = 50;

            // // Check for containment.
            // const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
            // const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);

            // // If 0 is not contained in 1, so forbid the alignment.
            // if (!(isContained0 && !isContained1)) {
            //     forbiddenAlignments[alignmentId] = true;
            //     continue;
            // }

            // // Reaching here meand that 0 is contained in 1
            // readsContainingReadId0.push_back(currentOrientedReadId1.getReadId());
        }

        // // If the read is contained in more than one read, find the alignment 
        // // with the best alignment score and forbid the other alignments.
        // if (readsContainingReadId0.size() > 1) {
        //     // Find the alignment with the lowest error rate
        //     uint64_t bestAlignmentId = 0;
        //     double bestAlignmentErrorRate = 1.;
        //     for (const uint64_t alignmentId: readsContainingReadId0) {
        //         const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
        //         const double alignmentErrorRate = thiAlignmentData.info.errorRate;
        //         if (alignmentErrorRate < bestAlignmentErrorRate) {
        //             bestAlignmentId = alignmentId;
        //             bestAlignmentErrorRate = alignmentErrorRate;
        //         }
        //     }

        //     // Forbid the other alignments.
        //     for (const uint64_t alignmentId: readsContainingReadId0) {
        //         if (alignmentId != bestAlignmentId) {
        //             forbiddenAlignments[alignmentId] = true;
        //         }
        //     }
        // }
        
    }


    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //     const AlignmentData& thiAlignmentData = alignmentData[alignmentId];
    //     ReadId currentReadId0 = thiAlignmentData.readIds[0];
    //     ReadId currentReadId1 = thiAlignmentData.readIds[1];
    //     OrientedReadId currentOrientedReadId0(thiAlignmentData.readIds[0], 0);
    //     OrientedReadId currentOrientedReadId1(thiAlignmentData.readIds[1], thiAlignmentData.isSameStrand ? 0 : 1);
    //     AlignmentInfo alignmentInfo = thiAlignmentData.info;

    //     // Procede only if both reads are contained.
    //     if (!isReadIdContained[currentReadId0] || !isReadIdContained[currentReadId1]) {
    //         continue;
    //     }

    //     bool foundNonContainedReadThatIsPartOfTheSiteReadSet = false;
    //     for (const uint64_t siteId: orientedReadSites[currentOrientedReadId0.getValue()]) {
    //         const Site &site = sites[siteId];
    //         for (const OrientedReadId &orientedReadId: site.orientedReads) {
    //             if (!isReadIdContained[orientedReadId.getReadId()]) {
    //                 foundNonContainedReadThatIsPartOfTheSiteReadSet = true;
    //                 break;
    //                 // forbiddenAlignments[alignmentId] = true;
    //                 // break;
    //             }
    //         }
    //     }

    //     if (foundNonContainedReadThatIsPartOfTheSiteReadSet) {
    //         for (const uint64_t siteId: orientedReadSites[currentOrientedReadId0.getValue()]) {
    //             const Site &site = sites[siteId];
    //             for (const OrientedReadId &orientedReadId: site.orientedReads) {
    //                 if (isReadIdContained[orientedReadId.getReadId()] && (orientedReadId.getReadId() == currentReadId1)) {
    //                     forbiddenAlignments[alignmentId] = true;
    //                     break;
    //                 }
    //             }
    //         }
    //     }

    // }





























    
    // // print all the reads in each site of the orientedReadId 19364-1
    // OrientedReadId testOrientedReadId(3094, 0);
    // cout << "OrientedReadId: " << testOrientedReadId << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId2(3117, 0);
    // cout << "OrientedReadId: " << testOrientedReadId2 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId2 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId2.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId3(3119, 0);
    // cout << "OrientedReadId: " << testOrientedReadId3 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId3 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId3.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId4(2949, 1);
    // cout << "OrientedReadId: " << testOrientedReadId4 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId4 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId4.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId5(2957, 0);
    // cout << "OrientedReadId: " << testOrientedReadId5 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId5 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId5.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId6(2935, 0);
    // cout << "OrientedReadId: " << testOrientedReadId6 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId6 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId6.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId7(1180, 0);
    // cout << "OrientedReadId: " << testOrientedReadId7 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId7 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId7.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId8(1185, 0);
    // cout << "OrientedReadId: " << testOrientedReadId8 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId8 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId8.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId9(1191, 1);
    // cout << "OrientedReadId: " << testOrientedReadId9 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId9 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId9.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId10(1194, 1);
    // cout << "OrientedReadId: " << testOrientedReadId10 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId10 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId10.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId11(1200, 0);
    // cout << "OrientedReadId: " << testOrientedReadId11 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId11 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId11.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId12(3118, 0);
    // cout << "OrientedReadId: " << testOrientedReadId12 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId12 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId12.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;



    // // print all the reads in each site of the orientedReadId 19364-1
    // OrientedReadId testOrientedReadId(3094, 0);
    // cout << "OrientedReadId: " << testOrientedReadId << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId2(3117, 0);
    // cout << "OrientedReadId: " << testOrientedReadId2 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId2 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId2.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId3(3119, 0);
    // cout << "OrientedReadId: " << testOrientedReadId3 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId3 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId3.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId4(2949, 1);
    // cout << "OrientedReadId: " << testOrientedReadId4 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId4 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId4.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId5(1055, 1);
    // cout << "OrientedReadId: " << testOrientedReadId5 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId5 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId5.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId6(1045, 0);
    // cout << "OrientedReadId: " << testOrientedReadId6 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId6 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId6.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId7(1057, 0);
    // cout << "OrientedReadId: " << testOrientedReadId7 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId7 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId7.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId8(1043, 1);
    // cout << "OrientedReadId: " << testOrientedReadId8 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId8 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId8.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId9(1498, 1);
    // cout << "OrientedReadId: " << testOrientedReadId9 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId9 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId9.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId10(1494, 1);
    // cout << "OrientedReadId: " << testOrientedReadId10 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId10 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId10.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId11(1489, 1);
    // cout << "OrientedReadId: " << testOrientedReadId11 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId11 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId11.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;

    // OrientedReadId testOrientedReadId12(1488, 0);
    // cout << "OrientedReadId: " << testOrientedReadId12 << endl;
    // cout << "Sites and their reads associated with orientedReadId " << testOrientedReadId12 << ": " << endl;
    // for (const uint64_t siteId: orientedReadSites[testOrientedReadId12.getValue()]) {
    //     cout << "Site " << siteId << " reads:" << endl;
    //     const Site &site = sites[siteId];
    //     for (const OrientedReadId &readId: site.orientedReads) {
    //         cout << "  ReadId: " << readId.getReadId() << "-" << readId.getStrand() << endl;
    //     }
    // }
    // cout << endl;


    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // if (!firstPassHetAlignments[alignmentId]) {
        //     continue;
        // }

        if (forbiddenAlignments[alignmentId]) {
            continue;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(thisAlignmentData.readIds[0], 0);
        OrientedReadId orientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);
        AlignmentInfo alignmentInfo = thisAlignmentData.info;

        // Swap them if necessary, depending on the average alignment offset at center.
        if(alignmentInfo.offsetAtCenter() < 0.) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }

        
    }



    //*
    //
    // Create the dynamically adjustable boost readGraph using 
    // all the alignments that are not forbidden.
    //
    //*
    
    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph5AllAlignments readGraphAllAlignments(orientedReadCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // if (!firstPassHetAlignments[alignmentId]) {
        //     continue;
        // }

        if (forbiddenAlignments[alignmentId]) {
            continue;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(thisAlignmentData.readIds[0], 0);
        OrientedReadId orientedReadId1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);
        AlignmentInfo alignmentInfo = thisAlignmentData.info;

        // Swap them if necessary, depending on the average alignment offset at center.
        if(alignmentInfo.offsetAtCenter() < 0.) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph5AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph5AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
    }
    
    cout << "The read graph with all alignments that are not forbidden has " << num_vertices(readGraphAllAlignments) << " vertices and " << num_edges(readGraphAllAlignments) << " edges." << endl;


    //
    //
    // We finished looping over the reads and found the haplotype sets for each one
    // Now we need to create the first pass read graph that involves only those haplotype specific 
    // alignments we found in the het sites
    //
    //


    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);
    // vector<bool> keepAlignment(alignmentCount, true);
    // createReadGraphUsingSelectedAlignments(keepAlignment);
    // return;

    // for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {

    //     AlignmentData& alignment = alignmentData[alignmentId];

    //     if (forbiddenAlignments[alignmentId]) {
    //         continue;
    //     }
        
    //     if (firstPassHetAlignments[alignmentId]) {
    //         // Add the alignment to the read graph.
    //         keepAlignment[alignmentId] = true;
    //         alignment.info.isInReadGraph = 1;
    //     }
    // }

    // createReadGraphUsingSelectedAlignments(keepAlignment);
    // return;

    
    //*
    //
    // Order alignments in order of increasing Q. 
    //
    // Gather in alignmentTable[alignmentID, Q]
    // alignments in order of increasing Q.
    // Q(n) = (1 + δ/2ε)^n * e-δL
    // ε = 1e-4, δ = 5e-4
    // logQ(n) = αn - δL
    //
    //*

    // const double epsilon = 1e-4;
    // const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // const double WThreshold = 1e-8;
    const double logWThreshold = log(WThreshold);

    // const double WThresholdForBreaks = 1e+15;
    const double logWThresholdForBreaks = log(WThresholdForBreaks);

    vector< pair<uint64_t, double> > alignmentTableHetSites;
    vector< pair<uint64_t, double> > alignmentTableHetSitesPlusNotForbidden;
    
    // Keep track of which readIds were used in alignments
    vector<bool> readUsed(readCount, false);



    //
    //
    //
    //
    // Do a first pass in which we allow only het loci in order of increasing Q, 
    // then do another pass to fill the breaks where you also allow all the other alignments
    //
    //
    //
    //



    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {

        // cout << "Processing alignment " << alignmentId << endl;
        
        if (!firstPassHetAlignments[alignmentId]) {
            continue;
        }

        if (forbiddenAlignments[alignmentId]) {
            continue;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // The alignment is stored as an alignment between readId0 on strand 0
        // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
        // The reverse complement alignment also exists, but is not stored explicitly.
        
        ReadId readId0 = thisAlignmentData.readIds[0];
        ReadId readId1 = thisAlignmentData.readIds[1];
        bool isSameStrand = thisAlignmentData.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
        OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.
        AlignmentInfo alignmentInfo = thisAlignmentData.info;

        // // Swap them if necessary, depending on the average alignment offset at center.
        // if(alignmentInfo.offsetAtCenter() < 0.) {
        //     swap(orientedReadId0, orientedReadId1);
        //     alignmentInfo.swap();
        // }

        // if (readId0 == 3118 || readId1 == 3118 || readId0 == 3251 || readId1 == 3251) {
        //     continue;
        // }

        // if (bridgingReads.count(readId0) || bridgingReads.count(readId1)) {
        //     continue;
        // }

        // // Swap them if necessary, depending on the average alignment offset at center.
        // if(thisAlignmentData.info.offsetAtCenter() < 0.) {
        //     swap(orientedReadId0, orientedReadId1);
        // }

        // if (isReadIdContained2[readId0] || isReadIdContained2[readId1]) {
        //     continue;
        // }

        // if (isReadIdContained2[readId1]) {
        //     continue;
        // }

        // if (orientedReadId0.getReadId() == 3579 || orientedReadId1.getReadId() == 3579) {
        //     if (orientedReadId0.getReadId() == 3386 || orientedReadId1.getReadId() == 3386) {
        //         continue;
        //     }
        // }

        double errorRate = alignmentInfo.errorRate;
        if (orientedReadId0.getReadId() == 3579 || orientedReadId1.getReadId() == 3579) {
            cout << "Error rate: " << errorRate << " between reads " << readId0 << " and " << readId1 << endl;
            // if (errorRate > 0.03) {
            //     cout << "Skipping alignment " << alignmentId << " because of error rate " << errorRate << " between reads " << readId0 << " and " << readId1 << endl;
            //     continue;
            // }
        }

        if (errorRate > 0.07) {
            cout << "Skipping alignment " << alignmentId << " because of error rate " << errorRate << " between reads " << readId0 << " and " << readId1 << endl;
            continue;
        }
            

        // Store this pair of edges in our edgeTable.
        const uint64_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
        const uint64_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
        const double L = double(range0 + range1)/2.;
        // const uint64_t n = thisAlignmentData.info.mismatchCountRle;
        
        const double nRLE = errorRate * 2 * L;
        // const double markerCount = thisAlignmentData.info.markerCount;

        // logQ(n) = αn - δL
        const double logQ = alpha * double(nRLE) - delta * L;

        // Add the alignment to the table along with its logQ.
        alignmentTableHetSites.push_back(make_pair(alignmentId, logQ));
        readUsed[readId0] = true;
        readUsed[readId1] = true;

        keepAlignment[alignmentId] = true;
        thisAlignmentData.info.isInReadGraph = 1;
        // cout << "Added alignment " << alignmentId << " to the read graph." << endl;

        // // This time use the regular Bayesian filtering
        // if (logQ <= logWThreshold) {
        //     alignmentTableHetSites.push_back(make_pair(alignmentId, logQ));
        //     alignmentsAlreadyConsidered[alignmentId] = true;
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        //     keepAlignment[alignmentId] = true;
        //     thisAlignmentData.info.isInReadGraph = 1;
        // }
        // } else if(logQ <= logWThresholdForBreaks){
        //     alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
        //     alignmentsAlreadyConsidered[alignmentId] = true;
        //     keepAlignmentsForBreaks[alignmentId] = true;
        // }

    }

    sort(alignmentTableHetSites.begin(), alignmentTableHetSites.end(), OrderPairsBySecondOnly<uint64_t, double>());
    cout << "The alignmentTableHetSites has " << alignmentTableHetSites.size() << " entries." << endl;



    createReadGraphUsingSelectedAlignments(keepAlignment);
    return;


}


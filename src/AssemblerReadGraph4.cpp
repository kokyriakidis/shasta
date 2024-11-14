#include "Assembler.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
#include <queue>
#include <map>
#include "orderPairs.hpp"
#include "Mode3Assembler.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/poisson.hpp>

//class AlignmentStats{public: double errorRateRle; uint32_t alignedRange; uint32_t rightUnaligned; uint32_t leftUnaligned; uint32_t alignmentId;};


void Assembler::createReadGraph4withStrandSeparation(
    uint32_t maxAlignmentCount)
{
    cout << timestamp << "createReadGraph4 with strand separation begins" << endl;

    const double QThreshold = 1e-5;
    const double logQThreshold = log(QThreshold);

    //
    // 1. Order alignments in order of increasing Q. 
    //

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Gather in alignmentTable[alignmentID, Q]
    // alignments in order of increasing Q.
    // Q(n) = (1 + δ/2ε)^n * e-δL
    // ε = 1e-4, δ = 5e-4
    // logQ(n) = αn - δL
    vector< pair<uint64_t, double> > alignmentTable;
    vector< pair<uint64_t, double> > alignmentTableFiltered;
    const double epsilon = 1e-4;
    const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // Loop over all alignments.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if((alignmentId % 10000) == 0) {
            cout << timestamp << alignmentId << "/" << alignmentCount << endl;
        }

        // Get information for this alignment.
        AlignmentData& thisAlignmentData = alignmentData[alignmentId];

        // The alignment is stored as an alignment between readId0 on strand 0
        // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
        // The reverse complement alignment also exists, but is not stored explicitly.
        const ReadId readId0 = thisAlignmentData.readIds[0];
        const ReadId readId1 = thisAlignmentData.readIds[1];
        const bool isSameStrand = thisAlignmentData.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
        const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

        // Store this pair of edges in our edgeTable.
        const uint32_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
        const uint32_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = thisAlignmentData.info.mismatchCountRle;     

        // logQ(n) = αn - δL
        const double logQ = alpha * double(n) - delta * L;


        if (logQ <= logQThreshold) {
            alignmentTableFiltered.push_back(make_pair(alignmentId, logQ));
        }
        alignmentTable.push_back(make_pair(alignmentId, logQ));
    }


    // Sort by increasing Q
    sort(alignmentTable.begin(), alignmentTable.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTableFiltered.begin(), alignmentTableFiltered.end(), OrderPairsBySecondOnly<uint64_t, double>());

    cout << "The alignment table has " << alignmentTable.size() << " entries." << endl; // 123863

    // Print the first 100 entries of the alignment table
    cout << "First 10 entries of the alignment table:" << endl;
    cout << "AlignmentId\tlogQ\tReadId0\tReadId1" << endl;
    for(size_t i=0; i<min(size_t(10), alignmentTable.size()); i++) {
        const uint64_t alignmentId = alignmentTable[i].first;
        const AlignmentData& alignment = alignmentData[alignmentId];
        cout << alignmentId << "\t" 
             << alignmentTable[i].second << "\t"
             << alignment.readIds[0] << "\t" 
             << alignment.readIds[1] << endl;
    }

    // Create and initialize the disjoint sets data structure needed below.
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;

    // Keep track of which readIds were used in alignments
    vector<bool> readUsed(readCount, false);
    for(size_t i=0; i<alignmentTableFiltered.size(); i++) {
        const uint64_t alignmentId = alignmentTableFiltered[i].first;
        const AlignmentData& alignment = alignmentData[alignmentId];
        readUsed[alignment.readIds[0]] = true;
        readUsed[alignment.readIds[1]] = true;
    }

    // Create a map of readId -> vector of alignment stats for alignments involving that read
    std::map<ReadId, vector<pair<uint64_t, double>>> readToAlignments;

    // Loop over alignments in alignmentTable
    for(size_t i=0; i<alignmentTable.size(); i++) {
        const uint64_t alignmentId = alignmentTable[i].first;
        const AlignmentData& alignment = alignmentData[alignmentId];

        // check if the reads were used in alignments
        if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
            continue;
        }

        // check if the first read was used and the second was not in alignments 
        if(readUsed[alignment.readIds[0]] && !readUsed[alignment.readIds[1]]) {
            readToAlignments[alignment.readIds[1]].push_back(make_pair(alignmentId, alignmentTable[i].second));
        }

        // check if the second read was used and the first was not in alignments
        if(!readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
            readToAlignments[alignment.readIds[0]].push_back(make_pair(alignmentId, alignmentTable[i].second));
        }
    }

    // // Loop over all alignments.
    // for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
    //      {
    //         // Get information for this alignment.
    //         if((alignmentId % 10000) == 0) {
    //             cout << timestamp << alignmentId << "/" << alignmentCount << endl;
    //         }
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // check if the reads were used in alignments
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
    //         continue;
    //     }

    //     // check if the first read was used and the second was not in alignments
    //     if(readUsed[alignment.readIds[0]] && !readUsed[alignment.readIds[1]]) {
    //         readToAlignments[alignment.readIds[1]].push_back(make_pair(alignmentId, alignmentTable[alignmentId].second));
    //     }

    //     // check if the second read was used and the first was not in alignments
    //     if(!readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
    //         readToAlignments[alignment.readIds[0]].push_back(make_pair(alignmentId, alignmentTable[alignmentId].second));
    //     }

    // }

    // For each read that was used in alignments, keep only the best alignment
    for(const auto& readEntry : readToAlignments) {
        if(!readEntry.second.empty()) {
            // Find alignment with lowest Q value (best alignment)
            auto bestAlignment = std::min_element(
                readEntry.second.begin(),
                readEntry.second.end(),
                [](const auto& a, const auto& b) {
                    return a.second < b.second;
                });
                
            // Mark all other alignments for this read as not to be kept
            for(const auto& alignment : readEntry.second) {
                if(alignment != *bestAlignment) {
                    alignmentData[alignment.first].info.isInReadGraph = false;
                }
            }
        }
    }

    

    ///
    // 2. Process alignments in order of increasing Q. 
    //
    // i.   Start with no edges in the read graph. 
    // ii.  Process alignments in order of increasing Q. 
    // iii. If the alignment breaks strand separation, it is skipped. 
    // iv.  If both vertices of the potential edge have at least the required minimum number of neighbors, the alignment is also skipped.  
    // v.   Otherwise, the pair of reverse complement edges corresponding to the alignment are added to the read graph.


    

    

    // Maintain a vector containing the degree of each vertex
    // verticesDegree[vertexID] -> degree
    vector<uint64_t> verticesDegree(orientedReadCount, 0);

    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }
    
    cout << "Number of reads: " << readCount << endl; // 6600
    cout << "Number of oriented reads: " << orientedReadCount << endl; // 13200
    cout << "The vertices degree vector has " << verticesDegree.size() << " entries." << endl; // 13200
    // Print the first 10 entries of the verticesDegree table
    cout << "First 10 entries of the vertices degree table:" << endl;
    cout << "VertexId\tDegree" << endl;
    for(size_t i=0; i<min(size_t(10), verticesDegree.size()); i++) {
        cout << i << "\t" 
             << verticesDegree[i] << endl;
    }


    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);

    // Process alignments in order of increasing Q
    uint64_t crossStrandEdgeCount = 0;
    for(auto it=alignmentTable.begin(); it!=alignmentTable.end(); ++it) {
        const pair<uint64_t, double>& p = *it;
        const uint64_t alignmentId = p.first;

        // Get the alignment data
        AlignmentData& alignment = alignmentData[alignmentId];
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        const bool isSameStrand = alignment.isSameStrand;
        SHASTA_ASSERT(readId0 < readId1);
        const OrientedReadId A0 = OrientedReadId(readId0, 0);
        const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
        const OrientedReadId A1 = OrientedReadId(readId0, 1);
        const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

        SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
        SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
        SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
        SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

        // Get the connected components that these oriented reads are in.
        const uint32_t a0 = disjointSets.find_set(A0.getValue());
        const uint32_t b0 = disjointSets.find_set(B0.getValue());
        const uint32_t a1 = disjointSets.find_set(A1.getValue());
        const uint32_t b1 = disjointSets.find_set(B1.getValue());


        // If the alignment breaks strand separation, it is skipped.
        // If A0 and B1 are in the same connected component,
        // A1 and B0 also must be in the same connected component.
        // Adding this pair of edges would create a self-complementary
        // connected component containing A0, B0, A1, and B1,
        // and to ensure strand separation we don't want to do that.
        // So we mark these edges as cross-strand edges
        // and don't use them to update the disjoint set data structure.
        if(a0 == b1) {
            SHASTA_ASSERT(a1 == b0);
            crossStrandEdgeCount += 2;
            continue;
        }

        // If both vertices of the potential edge have at least the required minimum number 
        // of neighbors, the alignment is also skipped. 
        const uint64_t degreeA0 = verticesDegree[A0.getValue()];
        const uint64_t degreeB0 = verticesDegree[B0.getValue()];
        if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
            // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
            continue;
        }

        // Add the alignment to the read graph.
        keepAlignment[alignmentId] = true;
        alignment.info.isInReadGraph = 1;

        // Update vertex degrees
        verticesDegree[A0.getValue()]++;
        verticesDegree[B0.getValue()]++;
        verticesDegree[A1.getValue()]++;
        verticesDegree[B1.getValue()]++;
        

        // Update disjoint sets
        disjointSets.union_set(a0, b0);
        disjointSets.union_set(a1, b1);
    }

    // Print the first 100 vertices and their degrees
    cout << "First 100 vertices and their degrees:" << endl;
    for(uint64_t i=0; i<min(uint64_t(100), uint64_t(verticesDegree.size())); i++) {
        cout << "Vertex " << i << ": degree " << verticesDegree[i] << endl;
    }
    

    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }

    cout << timestamp << "Done processing alignments." << endl;

    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;


    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);

    cout << timestamp << "createReadGraph4 with strand separation ends." << endl;

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << readGraph.edges.size() << " total." << endl;


    // Remove bridges from the read graph.
    // removeReadGraphBridges(5);

    // cout << timestamp << "Bridge removal ends." << endl;


    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    
    cout << "The read graph has " << componentMap.size() << " connected components." << endl;



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // Stop writing when we reach connected components
        // consisting of a single isolated read.
        if(component.size() == 1) {
            break;
        }

        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        // The above process of strand separation should have removed
        // all self-complementary components.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        SHASTA_ASSERT(not isSelfComplementary);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // For Mode 2 and Mode 3 assembly, we will only assemble one connected component
    // of each pair. In each pair, we choose the component in the pair
    // that has the lowest numbered read on strand 0.
    // Then, for each read we store in its ReadFlags the strand
    // that the read appears in in this component.
    // That flag will be used in Mode 2 assembly to
    // select portions of the marker graph that should be assembled.
    uint64_t n = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If the lowest numbered read is on strand 1, this is not one of
        // the connected components we want to use.
        if(component.front().getStrand() == 1) {
            continue;
        }

        // Store the strand for each read in this component.
        for(const OrientedReadId orientedReadId: component) {
            reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
        }
        n += component.size();
    }
    SHASTA_ASSERT(n == readCount);


}
























void Assembler::createReadGraph4(
    uint32_t maxAlignmentCount)
{
    const bool debug = false;

    bool assemblyIsAvailable = false;
    try {
        accessMode3Assembler();
        SHASTA_ASSERT(mode3Assembler);
        const mode3::Anchors& anchors = mode3Assembler->anchors();
        SHASTA_ASSERT(anchors.journeys.isOpen());
        SHASTA_ASSERT(anchors.journeys.size() == markers.size());
        assemblyIsAvailable= true;
    } catch (...) {
    }


    if(assemblyIsAvailable) { 
        // Skip if assembly is available
    } else {
        // QRle threshold to use an alignment in the read graph.
        const double minQRle = 26;

        const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

        cout << timestamp << "createReadGraph4 begins, minQRle " << minQRle << endl;

        // Get the total number of stored alignments.
        const uint64_t alignmentCount = alignmentData.size();
        SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

        // Flag all alignments as not to be kept.
        vector<bool> keepAlignment(alignmentCount, false);

        // Loop over all alignments.
        for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
            if((alignmentId % 1000) == 0) {
                cout << timestamp << alignmentId << "/" << alignmentCount << endl;
            }

            // Get information for this alignment.
            AlignmentData& thisAlignmentData = alignmentData[alignmentId];
            // keepAlignment[alignmentId] = true;
            // thisAlignmentData.info.isInReadGraph = 1;

            // The alignment is stored as an alignment between readId0 on strand 0
            // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
            // The reverse complement alignment also exists, but is not stored explicitly.
            const ReadId readId0 = thisAlignmentData.readIds[0];
            const ReadId readId1 = thisAlignmentData.readIds[1];
            const bool isSameStrand = thisAlignmentData.isSameStrand;
            SHASTA_ASSERT(readId0 < readId1);
            const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
            const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

            const uint32_t range0 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId0, 0, markers);
            const uint32_t range1 = thisAlignmentData.info.baseRange(assemblerInfo->k, orientedReadId1, 1, markers);
            const double L = (range0 + range1)/2;
            const uint64_t n = thisAlignmentData.info.mismatchCountRle;
            const double errorRateRle = double(n)/(2.0*L);;

            // If the RLE Q is large enough, flag thus alignment as to be kept.
            if((errorRateRle <= maxErrorRateRle)) {
                keepAlignment[alignmentId] = true;
                thisAlignmentData.info.isInReadGraph = 1;
            } else {
                keepAlignment[alignmentId] = false;
                thisAlignmentData.info.isInReadGraph = 0;
            }

        }

        cout << timestamp << "Done processing alignments." << endl;

        const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
        cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

        // Create the read graph using the alignments we selected.
        createReadGraphUsingSelectedAlignments(keepAlignment);

        cout << timestamp << "Initial createReadGraph4 ends." << endl;

        // Remove bridges from the read graph.
        removeReadGraphBridges(5);

        cout << timestamp << "createReadGraph4 ends." << endl;
    }

    

}




// Strict strand separation in the read graph based on using Poisson CDF.
// This guarantees that the read graph contains no self-complementary
// connected components.
// In other words, for any ReadId x, the two oriented reads
// x-0 and x-1 are guaranteed to be in distinct components of the
// read graph.
void Assembler::flagCrossStrandReadGraphEdges4()
{
    // Each alignment used in the read graph generates a pair of
    // consecutively numbered edges in the read graph
    // which are the reverse complement of each other.
    SHASTA_ASSERT((readGraph.edges.size() % 2) == 0);

    // Below, each pair is identified by the (even) index of
    // the first edge in the pair.

    // For each number of aligned markers alignedMarkerCount,
    // Gather in edgeTable[alignedMarkerCount] pairs
    // with that number of aligned markers.
    vector< pair<uint64_t, double> > edgeTable;

    // To loop over pairs of edges, increment by 2.
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId+=2) {
        const ReadGraphEdge& edge = readGraph.edges[edgeId];
        SHASTA_ASSERT(not edge.crossesStrands);

        // Skip edges flagged as inconsistent.
        if(edge.hasInconsistentAlignment) {
            continue;
        }

        const uint64_t alignmentId = edge.alignmentId;
        const AlignmentData& alignment = alignmentData[alignmentId];

        // Skip edges involving reads classified as chimeric.
        if(getReads().getFlags(alignment.readIds[0]).isChimeric) {
            continue;
        }
        if(getReads().getFlags(alignment.readIds[1]).isChimeric) {
            continue;
        }

        // Sanity check.
        SHASTA_ASSERT(alignmentData[alignmentId].info.isInReadGraph);

        // Check that the next edge is the reverse complement of
        // this edge.
        {
            const ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
            SHASTA_ASSERT(not nextEdge.crossesStrands);
            array<OrientedReadId, 2> nextEdgeOrientedReadIds = nextEdge.orientedReadIds;
            nextEdgeOrientedReadIds[0].flipStrand();
            nextEdgeOrientedReadIds[1].flipStrand();
            SHASTA_ASSERT(nextEdgeOrientedReadIds == edge.orientedReadIds);
            SHASTA_ASSERT(nextEdge.alignmentId == alignmentId);
        }


        // // Store this pair of edges in our edgeTable.
        // const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[0], 0, markers);
        // const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[1], 1, markers);
        // const double L = (range0 + range1) / 2;
        // const uint64_t n = alignment.info.mismatchCountRle;
        // const double m = L * 0.0001;

        // // Calculate Q using Poisson CDF
        // // Helper function to calculate gamma_q(k+1, λ)
        // auto calculate_gamma_q = [](uint64_t k, double lambda) -> double
        // {
        //     // gamma_q(a,z) is the normalized upper incomplete gamma function
        //     // For Poisson CDF, we need gamma_q(k+1, λ)
        //     double a = static_cast<double>(k + 1);
        //     return boost::math::gamma_q(a, lambda);
        // };

        // const double poissonCDF = calculate_gamma_q(n, m);


        // edgeTable.push_back(make_pair(edgeId, poissonCDF));

        // Store this pair of edges in our edgeTable.
        const double errorRateRle = alignment.info.errorRateRle;
        const double L = (alignment.info.range(0) + alignment.info.range(1)) / 0.04;
        const double m = L * 0.0001;
        const uint64_t n = static_cast<uint64_t>(std::round(errorRateRle * L));
        // Calculate Q using Poisson CDF
        // Helper function to calculate gamma_q(k+1, λ)
        auto calculate_gamma_q = [](uint64_t k, double lambda) -> double
        {
            // gamma_q(a,z) is the normalized upper incomplete gamma function
            // For Poisson CDF, we need gamma_q(k+1, λ)
            double a = static_cast<double>(k + 1);
            return boost::math::gamma_q(a, lambda);
        };
        const double poissonCDF = calculate_gamma_q(n, m);

        edgeTable.push_back(make_pair(edgeId, poissonCDF));
        
    }

    // Sort by increasing Q
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<uint64_t, double>());

    // Print out top 100 alignments by Q value
    cout << "Top 100 alignments by Q value:" << endl;
    cout << "EdgeId\tQ\tMarkerCount\tErrorRateRleNew\tErrorRateRleOld\tL\tn" << endl;
    for(size_t i = 0; i < min(size_t(100), edgeTable.size()); i++) {
        const auto& p = edgeTable[i];
        const uint64_t edgeId = p.first;
        const double q = p.second;
        const AlignmentData& alignment = alignmentData[readGraph.edges[edgeId].alignmentId];
        const uint64_t markerCount = alignment.info.markerCount;
        const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[0], 0, markers);
        const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[1], 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = alignment.info.mismatchCountRle;
        const double errorRateRleNew = double(n)/(2.0*L);;
        const double errorRateRleOld = alignment.info.errorRateRle;
        const double m = L * 0.0001;
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        cout << edgeId << "(" << readId0 << "," << readId1 << ")" << "\t" << q << "\t" << markerCount << "\t" << errorRateRleNew << "\t" << errorRateRleOld << "\t" << L << "\t" << n << endl;
    }


    // Create and initialize the disjoint sets data structure needed below.
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    SHASTA_ASSERT(readGraph.connectivity.size() == orientedReadCount);
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }


    // Now process our edge table in order of decreasing alignedMarkerCount.
    // Each pair consists of two edges:
    // A0--B0
    // A1--B1
    // where A1 is the reverse complement of A0 and B1 is the reverse complement of B0.
    // Maintain a disjoint set data structure for read graph vertices.
    // For each pair, compute the current connected components
    // for A0 B0 A1 B1 from the disjoint set data structure, call them a0 b0 a1 b1.
    // If a1==b0 (in which case also b1==a0), adding these edges
    // would create a self-complementary connected component,
    // so we mark them as cross-strand edges and don't add them to the
    // disjoint set data structure.
    // Otherwise, the two edges are added to the disjoint set data structure.
    uint64_t crossStrandEdgeCount = 0;
    for(auto it=edgeTable.begin(); it!=edgeTable.end(); ++it) {
        const pair<uint64_t, double>& p = *it;
        const uint64_t edgeId = p.first;
        ReadGraphEdge& edge = readGraph.edges[edgeId];
        ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
        const OrientedReadId A0 = edge.orientedReadIds[0];
        const OrientedReadId B0 = edge.orientedReadIds[1];
        const OrientedReadId A1 = nextEdge.orientedReadIds[0];
        const OrientedReadId B1 = nextEdge.orientedReadIds[1];

        SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
        SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
        SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
        SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

        // Get the connected components that these oriented reads are in.
        const uint32_t a0 = disjointSets.find_set(A0.getValue());
        const uint32_t b0 = disjointSets.find_set(B0.getValue());
        const uint32_t a1 = disjointSets.find_set(A1.getValue());
        const uint32_t b1 = disjointSets.find_set(B1.getValue());

        // If A0 and B0 are in the same connected component,
        // A1 and B1 also must be in the same connected component.
        // There is nothing to do as this pair of edges
        // does not affect the connected components.
        if (a0 == b0) {
            SHASTA_ASSERT(a1 == b1);
            continue;
        }

        // If A0 and B1 are in the same connected component,
        // A1 and B0 also must be in the same connected component.
        // Adding this pair of edges would create a self-complementary
        // connected component containing A0, B0, A1, and B1,
        // and to ensure strand separation we don't want to do that.
        // So we mark these edges as cross-strand edges
        // and don't use them to update the disjoint set data structure.
        if(a0 == b1) {
            SHASTA_ASSERT(a1 == b0);
            edge.crossesStrands = 1;
            nextEdge.crossesStrands = 1;
            alignmentData[edge.alignmentId].info.isInReadGraph = false;
            crossStrandEdgeCount += 2;
            continue;
        }

        // Otherwise, just update the disjoint sets data structure
        // with these two edges.
        disjointSets.union_set(a0, b0);
        disjointSets.union_set(a1, b1);
    }

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << readGraph.edges.size() << " total." << endl;



    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }


    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    // cout << "The read graph has " << componentMap.size() << " connected components." << endl;



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // Stop writing when we reach connected components
        // consisting of a single isolated read.
        if(component.size() == 1) {
            break;
        }

        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        // The above process of strand separation should have removed
        // all self-complementary components.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        SHASTA_ASSERT(not isSelfComplementary);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // For Mode 2 assembly, we will only assemble one connected component
    // of each pair. In each pair, we choose the component in the pair
    // that has the lowest numbered read on strand 0.
    // Then, for each read we store in its ReadFlags the strand
    // that the read appears in in this component.
    // That flag will be used in Mode 2 assembly to
    // select portions of the marker graph that should be assembled.
    uint64_t n = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If the lowest numbered read is on strand 1, this is not one of
        // the connected components we want to use.
        if(component.front().getStrand() == 1) {
            continue;
        }

        // Store the strand for each read in this component.
        for(const OrientedReadId orientedReadId: component) {
            reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
        }
        n += component.size();
    }
    SHASTA_ASSERT(n == readCount);

}









// Strict strand separation in the read graph based on using Poisson CDF.
// This guarantees that the read graph contains no self-complementary
// connected components.
// In other words, for any ReadId x, the two oriented reads
// x-0 and x-1 are guaranteed to be in distinct components of the
// read graph.
void Assembler::flagCrossStrandReadGraphEdges5()
{
    // Each alignment used in the read graph generates a pair of
    // consecutively numbered edges in the read graph
    // which are the reverse complement of each other.
    SHASTA_ASSERT((readGraph.edges.size() % 2) == 0);

    // Below, each pair is identified by the (even) index of
    // the first edge in the pair.

    // For each number of aligned markers alignedMarkerCount,
    // Gather in edgeTable[alignedMarkerCount] pairs
    // with that number of aligned markers.
    vector< pair<uint64_t, double> > edgeTable;
    const double epsilon = 1e-4;
    const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // To loop over pairs of edges, increment by 2.
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId+=2) {
        const ReadGraphEdge& edge = readGraph.edges[edgeId];
        SHASTA_ASSERT(not edge.crossesStrands);

        // Skip edges flagged as inconsistent.
        if(edge.hasInconsistentAlignment) {
            continue;
        }

        const uint64_t alignmentId = edge.alignmentId;
        const AlignmentData& alignment = alignmentData[alignmentId];

        // Skip edges involving reads classified as chimeric.
        if(getReads().getFlags(alignment.readIds[0]).isChimeric) {
            continue;
        }
        if(getReads().getFlags(alignment.readIds[1]).isChimeric) {
            continue;
        }

        // Sanity check.
        SHASTA_ASSERT(alignmentData[alignmentId].info.isInReadGraph);

        // Check that the next edge is the reverse complement of
        // this edge.
        {
            const ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
            SHASTA_ASSERT(not nextEdge.crossesStrands);
            array<OrientedReadId, 2> nextEdgeOrientedReadIds = nextEdge.orientedReadIds;
            nextEdgeOrientedReadIds[0].flipStrand();
            nextEdgeOrientedReadIds[1].flipStrand();
            SHASTA_ASSERT(nextEdgeOrientedReadIds == edge.orientedReadIds);
            SHASTA_ASSERT(nextEdge.alignmentId == alignmentId);
        }
        

        // Q(n) = (1 + δ/2ε)^n * e-δL
        // ε = 1e-4, δ = 5e-4
        // Store this pair of edges in our edgeTable.
        const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[0], 0, markers);
        const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, edge.orientedReadIds[1], 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = alignment.info.mismatchCountRle;     

        // logQ(n) = αn - δL
        const double logQ_n = alpha * double(n) - delta * L;

        edgeTable.push_back(make_pair(edgeId, logQ_n));
        
    }

    // Sort by increasing Q
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnly<uint64_t, double>());

    // Print out top 100 alignments by Q value
    cout << "Top 100 alignments by logQ value:" << endl;
    cout << "EdgeId\tlogQ\tMarkerCount\tErrorRateRle\tL\tn" << endl;
    for(size_t i = 0; i < min(size_t(100), edgeTable.size()); i++) {
        const auto& p = edgeTable[i];
        const uint64_t edgeId = p.first;
        const double logQ = p.second;
        const AlignmentData& alignment = alignmentData[readGraph.edges[edgeId].alignmentId];
        const uint64_t markerCount = alignment.info.markerCount;
        const uint32_t range0 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[0], 0, markers);
        const uint32_t range1 = alignment.info.baseRange(assemblerInfo->k, readGraph.edges[edgeId].orientedReadIds[1], 1, markers);
        const double L = (range0 + range1)/2;
        const uint64_t n = alignment.info.mismatchCountRle;
        const double errorRateRle = double(n)/(2.0*L);;
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        cout << edgeId << " ( " << readId0 << "," << readId1 << " )" << "\t" << logQ << "\t" << markerCount << "\t" << errorRateRle << "\t" << L << "\t" << n << endl;
    }


    // Create and initialize the disjoint sets data structure needed below.
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    SHASTA_ASSERT(readGraph.connectivity.size() == orientedReadCount);
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }


    // Now process our edge table in order of decreasing alignedMarkerCount.
    // Each pair consists of two edges:
    // A0--B0
    // A1--B1
    // where A1 is the reverse complement of A0 and B1 is the reverse complement of B0.
    // Maintain a disjoint set data structure for read graph vertices.
    // For each pair, compute the current connected components
    // for A0 B0 A1 B1 from the disjoint set data structure, call them a0 b0 a1 b1.
    // If a1==b0 (in which case also b1==a0), adding these edges
    // would create a self-complementary connected component,
    // so we mark them as cross-strand edges and don't add them to the
    // disjoint set data structure.
    // Otherwise, the two edges are added to the disjoint set data structure.
    uint64_t crossStrandEdgeCount = 0;
    for(auto it=edgeTable.begin(); it!=edgeTable.end(); ++it) {
        const pair<uint64_t, double>& p = *it;
        const uint64_t edgeId = p.first;
        ReadGraphEdge& edge = readGraph.edges[edgeId];
        ReadGraphEdge& nextEdge = readGraph.edges[edgeId + 1];
        const OrientedReadId A0 = edge.orientedReadIds[0];
        const OrientedReadId B0 = edge.orientedReadIds[1];
        const OrientedReadId A1 = nextEdge.orientedReadIds[0];
        const OrientedReadId B1 = nextEdge.orientedReadIds[1];

        SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
        SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
        SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
        SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

        // Get the connected components that these oriented reads are in.
        const uint32_t a0 = disjointSets.find_set(A0.getValue());
        const uint32_t b0 = disjointSets.find_set(B0.getValue());
        const uint32_t a1 = disjointSets.find_set(A1.getValue());
        const uint32_t b1 = disjointSets.find_set(B1.getValue());

        // If A0 and B0 are in the same connected component,
        // A1 and B1 also must be in the same connected component.
        // There is nothing to do as this pair of edges
        // does not affect the connected components.
        if (a0 == b0) {
            SHASTA_ASSERT(a1 == b1);
            continue;
        }

        // If A0 and B1 are in the same connected component,
        // A1 and B0 also must be in the same connected component.
        // Adding this pair of edges would create a self-complementary
        // connected component containing A0, B0, A1, and B1,
        // and to ensure strand separation we don't want to do that.
        // So we mark these edges as cross-strand edges
        // and don't use them to update the disjoint set data structure.
        if(a0 == b1) {
            SHASTA_ASSERT(a1 == b0);
            edge.crossesStrands = 1;
            nextEdge.crossesStrands = 1;
            alignmentData[edge.alignmentId].info.isInReadGraph = false;
            crossStrandEdgeCount += 2;
            continue;
        }

        // Otherwise, just update the disjoint sets data structure
        // with these two edges.
        disjointSets.union_set(a0, b0);
        disjointSets.union_set(a1, b1);
    }

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << readGraph.edges.size() << " total." << endl;



    // Verify that for any read the two oriented reads are in distinct
    // connected components.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(readId, 0);
        const OrientedReadId orientedReadId1(readId, 1);
        SHASTA_ASSERT(
            disjointSets.find_set(orientedReadId0.getValue()) !=
            disjointSets.find_set(orientedReadId1.getValue())
        );
    }


    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    // cout << "The read graph has " << componentMap.size() << " connected components." << endl;



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    performanceLog << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // Stop writing when we reach connected components
        // consisting of a single isolated read.
        if(component.size() == 1) {
            break;
        }

        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        // The above process of strand separation should have removed
        // all self-complementary components.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        SHASTA_ASSERT(not isSelfComplementary);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // For Mode 2 assembly, we will only assemble one connected component
    // of each pair. In each pair, we choose the component in the pair
    // that has the lowest numbered read on strand 0.
    // Then, for each read we store in its ReadFlags the strand
    // that the read appears in in this component.
    // That flag will be used in Mode 2 assembly to
    // select portions of the marker graph that should be assembled.
    uint64_t n = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If the lowest numbered read is on strand 1, this is not one of
        // the connected components we want to use.
        if(component.front().getStrand() == 1) {
            continue;
        }

        // Store the strand for each read in this component.
        for(const OrientedReadId orientedReadId: component) {
            reads->setStrandFlag(orientedReadId.getReadId(), orientedReadId.getStrand());
        }
        n += component.size();
    }
    SHASTA_ASSERT(n == readCount);

}











// Add this helper method to properly remove the read graph
void Assembler::removeReadGraph()
{
    // Close existing read graph if it's open
    if (readGraph.edges.isOpen) {
        readGraph.edges.close();
    }
    if (readGraph.connectivity.isOpen()) {
        readGraph.connectivity.close();
    }

    cout << timestamp << "Removed previous read graph files" << endl;
}






















































































// void Assembler::createReadGraph4(
//     uint32_t maxAlignmentCount)
// {
//     const bool debug = false;

//     // QRle threshold to use an alignment in the read graph.
//     const double minQRle = 100000.;

//     const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

//     cout << timestamp << "createReadGraph4 begins, maxAlignmentCount skata " << maxAlignmentCount << endl;

//     // Get the total number of stored alignments.
//     const uint64_t alignmentCount = alignmentData.size();
//     SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

//     // Flag all alignments as not to be kept.
//     vector<bool> keepAlignment(alignmentCount, false);

//     // The computation of projected alignments is expensive, and so it should be done once for each alignment 
//     // in an initial step and store the Error Rates of the projected alignment
//     vector<double> alignmentErrorRateRle(alignmentCount);

//     // This will hold the decomepressed Alignment.
//     // Defined here to reduce memory allocation activity.
//     Alignment alignment;

//     // Loop over all alignments.
//     for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
//         if((alignmentId % 1000) == 0) {
//             cout << timestamp << alignmentId << "/" << alignmentCount << endl;
//         }

//         // Get information for this alignment.
//         AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         thisAlignmentData.info.isInReadGraph = 0;

//         // The alignment is stored as an alignment between readId0 on strand 0
//         // and readId1 on strand 0 or 1 depending on the value of isSameStrand.
//         // The reverse complement alignment also exists, but is not stored explicitly.
//         const ReadId readId0 = thisAlignmentData.readIds[0];
//         const ReadId readId1 = thisAlignmentData.readIds[1];
//         const bool isSameStrand = thisAlignmentData.isSameStrand;
//         SHASTA_ASSERT(readId0 < readId1);
//         const OrientedReadId orientedReadId0(readId0, 0);   // On strand 0.
//         const OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);   // On strand 0 or 1.

//         // The alignment is stored in compressed form as a string,
//         // so we have to decompress it.
//         span<const char> compressedAlignment = compressedAlignments[alignmentId];
//         shasta::decompress(compressedAlignment, alignment);

//         // Project this alignment to base space.
//         const ProjectedAlignment projectedAlignment(
//             *this,
//             {orientedReadId0, orientedReadId1},
//             alignment,
//             true);

//         const double errorRateRle = projectedAlignment.errorRateRle();

//         alignmentErrorRateRle[alignmentId] = errorRateRle;

//     }


//     // Vector to keep the alignments for each read,
//     // with their number of markers.
//     // Contains pairs(errorRateRle, alignment id).
//     vector<AlignmentStats> readAlignments;

    
//     // Find the number of reads and oriented reads.
//     const ReadId orientedReadCount = uint32_t(markers.size());
//     SHASTA_ASSERT((orientedReadCount % 2) == 0);
//     const ReadId readCount = orientedReadCount / 2;

//     // Loop over reads.
//     for(ReadId readId=0; readId<readCount; readId++) {
//         if(debug) {
//             cout << "Working on read " << readId << endl;
//         }

//         OrientedReadId alignmentOrientedReadId0(readId, 0);

//         // Loop over the alignments that this oriented read is involved in, with the proper orientation.
//         const vector< pair<OrientedReadId, AlignmentInfo> > correctOrientedAlignments =
//             findOrientedAlignments(alignmentOrientedReadId0, false);



//         // Gather the alignments for this read, considering alignment range and right unaligned portion.
//         readAlignments.clear();

//         for(uint32_t i=0; i<correctOrientedAlignments.size(); i++) {

//             const uint32_t alignmentId = alignmentTable[alignmentOrientedReadId0.getValue()][i];

//             const double errorRateRle = alignmentErrorRateRle[alignmentId];

//             // Calculate alignment range
//             const uint32_t alignmentRange = correctOrientedAlignments[i].second.markerCount;

//             // Calculate right unaligned portion
//             const uint32_t rightUnaligned0 = correctOrientedAlignments[i].second.rightTrim(0);

//             // Calculate left unaligned portion
//             const uint32_t leftUnaligned0 = correctOrientedAlignments[i].second.leftTrim(0);

//             // Calculate right unaligned portion
//             const uint32_t rightUnaligned1 = correctOrientedAlignments[i].second.rightTrim(1);

//             // Calculate left unaligned portion
//             const uint32_t leftUnaligned1 = correctOrientedAlignments[i].second.leftTrim(1);

//             // if(rightUnaligned1 != 0 and leftUnaligned1 != 0 and errorRateRle<= maxErrorRateRle) {
//             //     readAlignments.push_back(AlignmentStats{errorRateRle, alignmentRange, rightUnaligned1, leftUnaligned1, alignmentId});
//             // }

//             if(errorRateRle<= maxErrorRateRle) {
//                 readAlignments.push_back(AlignmentStats{errorRateRle, alignmentRange, rightUnaligned1, leftUnaligned1, alignmentId});
//             }

//         }
        
//         cout << "Working on read " << readId << endl;


//         // // Find the alignment with the highest alignedRange
//         // auto bestAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return a.alignedRange < b.alignedRange;
//         //     });

//         // if (bestAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestAlignment = *bestAlignmentIt;

//         //     // Clear existing alignments and keep only the best one
//         //     readAlignments.clear();
//         //     readAlignments.push_back(bestAlignment);

//         //     if(debug) {
//         //         cout << "Kept alignment with the highest alignedRange: " << bestAlignment.alignedRange << endl;
//         //     }
//         // } else {
//         //     if(debug) {
//         //         cout << "No alignments found for this read." << endl;
//         //     }
//         // }


//         // // Find the 10 alignments with the lowest errorRateRle
//         // std::partial_sort(readAlignments.begin(), 
//         //                   readAlignments.begin() + std::min(10UL, readAlignments.size()), 
//         //                   readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return a.errorRateRle < b.errorRateRle;
//         //     });

//         // // Find the 10 alignments with the highest sum of rightUnaligned, leftUnaligned, and alignedRange
//         // std::partial_sort(readAlignments.begin(), 
//         //                   readAlignments.begin() + std::min(10UL, readAlignments.size()), 
//         //                   readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         return (a.rightUnaligned + a.leftUnaligned + a.alignedRange) >
//         //                (b.rightUnaligned + b.leftUnaligned + b.alignedRange);
//         //     });

//         cout << "read " << readId << " has that many alignments skata:" << readAlignments.size() << endl;

//         // // Keep only the top 10 alignments (or fewer if there are less than 10)
//         // readAlignments.resize(std::min(10UL, readAlignments.size()));

//         // cout << "read " << readId << " has that many alignments:" << readAlignments.size() << endl;

//         // if(debug) {
//         //     cout << "Top 5 alignments (or fewer) based on sum of unaligned portions and aligned range:" << endl;
//         //     for(const auto& alignment : readAlignments) {
//         //         cout << "AlignmentId: " << alignment.alignmentId 
//         //              << ", Sum: " << (alignment.rightUnaligned + alignment.leftUnaligned + alignment.alignedRange)
//         //              << " (Right: " << alignment.rightUnaligned 
//         //              << ", Left: " << alignment.leftUnaligned 
//         //              << ", Aligned: " << alignment.alignedRange << ")" << endl;
//         //     }
//         // }


//         // // Find the alignment with the highest alignedRange which also has the highest rightUnaligned
//         // auto bestRightAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         // if (a.alignedRange == b.alignedRange) {
//         //         //     return a.rightUnaligned > b.rightUnaligned;
//         //         // }
//         //         // return a.alignedRange < b.alignedRange;
//         //         return a.rightUnaligned+a.leftUnaligned+a.alignedRange < b.rightUnaligned+b.leftUnaligned+b.alignedRange;
//         //     });

//         // if (bestRightAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestRightAlignment = *bestRightAlignmentIt;

//         //     readAlignmentsRightTrimmed.clear();
//         //     readAlignmentsRightTrimmed.push_back(bestRightAlignment);
//         // }


//         // // Find the alignment with the highest alignedRange which also has the highest leftUnaligned
//         // auto bestLeftAlignmentIt = std::max_element(readAlignments.begin(), readAlignments.end(),
//         //     [](const AlignmentStats& a, const AlignmentStats& b) {
//         //         // if (a.alignedRange == b.alignedRange) {
//         //         //     return a.leftUnaligned > b.leftUnaligned;
//         //         // }
//         //         // return a.alignedRange < b.alignedRange;
//         //         return a.leftUnaligned < b.leftUnaligned;
//         //     });

//         // if (bestLeftAlignmentIt != readAlignments.end()) {
//         //     // Keep only this best alignment
//         //     AlignmentStats bestLeftAlignment = *bestLeftAlignmentIt;

//         //     readAlignmentsLeftTrimmed.clear();
//         //     readAlignmentsLeftTrimmed.push_back(bestLeftAlignment);
//         // }

        

//         // if (!readAlignmentsRightTrimmed.empty()){
//         // cout << "The best rightTrimm alignment is " << readAlignmentsRightTrimmed[0].alignmentId << " alignmentId." << endl;
//         // cout << "The best rightTrimm has " << readAlignmentsRightTrimmed[0].rightUnaligned << " rightUnaligned." << endl;
//         // }

//         // if (!readAlignmentsLeftTrimmed.empty()){
//         //     cout << "The best leftTrimm alignment is " << readAlignmentsLeftTrimmed[0].alignmentId << " alignmentId." << endl;
//         //     cout << "The best leftTrimm has " << readAlignmentsRightTrimmed[0].leftUnaligned << " leftUnaligned." << endl;
//         // }
        

//         // if(debug) {
//         //     cout << "Found " << readAlignments.size() << " alignments." << endl;
//         // }

//         // // Keep the best maxAlignmentCount.
//         // if(readAlignments.size() > maxAlignmentCount) {
//         //     std::nth_element(
//         //         readAlignments.begin(),
//         //         readAlignments.begin() + maxAlignmentCount,
//         //         readAlignments.end(),
//         //         std::less< pair<double, uint32_t> >());
//         //     readAlignments.resize(maxAlignmentCount);
//         // }
//         // if(debug) {
//         //     cout << "Kept " << readAlignments.size() << " alignments." << endl;
//         // }

//         // Mark the surviving alignments as to be kept.
//         for(const auto& p: readAlignments) {
//             // Get information for this alignment.
//             const uint32_t alignmentId = p.alignmentId;
//             AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//             keepAlignment[alignmentId] = true;
//             thisAlignmentData.info.isInReadGraph = 1;
//             if(debug) {
//                 const AlignmentData& alignment = alignmentData[alignmentId];
//                 cout << "Marked alignment " << alignment.readIds[0] << " " <<
//                     alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
//             }
//         }

//         // // Mark the surviving alignments as to be kept.
//         // for(const auto& p: readAlignmentsLeftTrimmed) {
//         //     // Get information for this alignment.
//         //     const uint32_t alignmentId = p.alignmentId;
//         //     AlignmentData& thisAlignmentData = alignmentData[alignmentId];
//         //     keepAlignment[alignmentId] = true;
//         //     thisAlignmentData.info.isInReadGraph = 1;
//         //     if(debug) {
//         //         const AlignmentData& alignment = alignmentData[alignmentId];
//         //         cout << "Marked alignment " << alignment.readIds[0] << " " <<
//         //             alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
//         //     }
//         // }
//     }


//     cout << timestamp << "Done processing alignments." << endl;

//     const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
//     cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

//     // Create the read graph using the alignments we selected.
//     createReadGraphUsingSelectedAlignments(keepAlignment);

//     cout << timestamp << "createReadGraph4 ends." << endl;
// }

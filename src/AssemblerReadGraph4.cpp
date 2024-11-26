#include "Assembler.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
#include "orderPairs.hpp"
#include "Mode3Assembler.hpp"

using namespace shasta;

// Standard library.
#include "fstream.hpp"
#include "chrono.hpp"
#include "iterator.hpp"
#include <numeric>
#include <queue>
#include <random>
#include <stack>
#include <queue>

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
    class ReadGraph4;
    class ReadGraph4Vertex;
    class ReadGraph4Edge;

    using ReadGraph4BaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::undirectedS,
        ReadGraph4Vertex,
        ReadGraph4Edge>;

    class ReadGraph4AllAlignments;
    class ReadGraph4AllAlignmentsVertex;
    class ReadGraph4AllAlignmentsEdge;

    using ReadGraph4AllAlignmentsBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::undirectedS,
        ReadGraph4AllAlignmentsVertex,
        ReadGraph4AllAlignmentsEdge>;
}



class shasta::ReadGraph4Vertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

    // The distances from a starting vertex and its reverse complement.
    uint64_t distance0 = invalid<uint64_t>;
    uint64_t distance1 = invalid<uint64_t>;
};

class shasta::ReadGraph4AllAlignmentsVertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

    // The distances from a starting vertex and its reverse complement.
    uint64_t distance0 = invalid<uint64_t>;
    uint64_t distance1 = invalid<uint64_t>;
};


class shasta::ReadGraph4Edge {
public:
    uint64_t alignmentId;
    ReadGraph4Edge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};

class shasta::ReadGraph4AllAlignmentsEdge {
public:
    uint64_t alignmentId;
    ReadGraph4AllAlignmentsEdge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};

// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph4: public ReadGraph4BaseClass {
public:

    ReadGraph4(uint64_t n) : ReadGraph4BaseClass(n) {}
    void findNeighbors(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsSkipSameComponentNodes(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachSameComponentNode(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};

// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph4AllAlignments: public ReadGraph4AllAlignmentsBaseClass {
public:

    ReadGraph4AllAlignments(uint64_t n) : ReadGraph4AllAlignmentsBaseClass(n) {}
    void findNeighbors(OrientedReadId orientedReadId, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsSkipSameComponentNodes(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);
    void findNeighborsEarlyStopWhenReachSameComponentNode(OrientedReadId orientedReadId, boost::disjoint_sets<ReadId*, ReadId*>& disjointSets, uint64_t maxDistance, vector<OrientedReadId>& neighbors);

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};



void ReadGraph4::findNeighbors(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if(!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighbors(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors)
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        if (distance > 0) { // Don't add the starting vertex
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if(!visitedVertices.contains(targetVertex)) {
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}



void ReadGraph4::findNeighborsSkipSameComponentNodes(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId != sourceComponentId) {
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
            }
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    // Only explore vertices in different components
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId != sourceComponentId) {
                        visitedVertices.insert(targetVertex);
                        q.push(make_pair(targetVertex, distance + 1));
                    }
                }
            }
        }
    }
}



void ReadGraph4AllAlignments::findNeighborsSkipSameComponentNodes(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId != sourceComponentId) {
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
            }
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    // Only explore vertices in different components
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId != sourceComponentId) {
                        visitedVertices.insert(targetVertex);
                        q.push(make_pair(targetVertex, distance + 1));
                    }
                }
            }
        }
    }
}









void ReadGraph4::findNeighborsEarlyStopWhenReachSameComponentNode(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId == sourceComponentId) {
                // If we find a vertex in the same component, add it to the neighbors and stop
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId == sourceComponentId) {
                        // If we find a vertex in the same component, add it to the neighbors and stop
                        neighbors.push_back(OrientedReadId::fromValue(targetVertex));
                        return;
                    }
                    // Only explore vertices in different components
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}


void ReadGraph4AllAlignments::findNeighborsEarlyStopWhenReachSameComponentNode(
    OrientedReadId orientedReadId,
    boost::disjoint_sets<ReadId*, ReadId*>& disjointSets,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) 
{
    neighbors.clear();

    // Keep track of visited vertices
    std::set<vertex_descriptor> visitedVertices;

    // Queue for BFS
    std::queue<pair<vertex_descriptor, uint64_t>> q; // vertex and distance
    
    // Get the strong component ID of the starting vertex
    uint64_t sourceComponentId = disjointSets.find_set(orientedReadId.getValue());
    
    // Start BFS from orientedReadId
    q.push(make_pair(orientedReadId.getValue(), 0));
    visitedVertices.insert(orientedReadId.getValue());

    while (!q.empty()) {
        vertex_descriptor currentVertex = q.front().first;
        uint64_t distance = q.front().second;
        q.pop();

        // If this is not the starting vertex and it's in a different component,
        // add it to neighbors
        if (distance > 0) {
            uint64_t currentComponentId = disjointSets.find_set(currentVertex);
            if (currentComponentId == sourceComponentId) {
                // If we find a vertex in the same component, add it to the neighbors and stop
                neighbors.push_back(OrientedReadId::fromValue(currentVertex));
                return;
            }
            neighbors.push_back(OrientedReadId::fromValue(currentVertex));
        }

        // Only continue exploring if we haven't hit the max distance
        if (distance < maxDistance) {
            // Iterate over adjacent vertices
            BGL_FORALL_OUTEDGES(currentVertex, edge, *this, ReadGraph4AllAlignments) {
                vertex_descriptor targetVertex = target(edge, *this);
                if (!visitedVertices.contains(targetVertex)) {
                    uint64_t targetComponentId = disjointSets.find_set(targetVertex);
                    if (targetComponentId == sourceComponentId) {
                        // If we find a vertex in the same component, add it to the neighbors and stop
                        neighbors.push_back(OrientedReadId::fromValue(targetVertex));
                        return;
                    }
                    // Only explore vertices in different components
                    visitedVertices.insert(targetVertex);
                    q.push(make_pair(targetVertex, distance + 1));
                }
            }
        }
    }
}























//class AlignmentStats{public: double errorRateRle; uint32_t alignedRange; uint32_t rightUnaligned; uint32_t leftUnaligned; uint32_t alignmentId;};


void Assembler::createReadGraph4withStrandSeparation(
    uint32_t maxAlignmentCount)
{
    cout << timestamp << "createReadGraph4 with strand separation begins" << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // const double minQRle = 30.;
    // const double maxErrorRateRle = std::pow(10.0, - minQRle / 10.0);

    const double minQRle60 = 60.;
    const double maxErrorRateRle60 = std::pow(10.0, - minQRle60 / 10.0);

    const double minQRle50 = 50.;
    const double maxErrorRateRle50 = std::pow(10.0, - minQRle50 / 10.0);

    const double minQRle40 = 40.;
    const double maxErrorRateRle40 = std::pow(10.0, - minQRle40 / 10.0);

    const double minQRle30 = 30.;
    const double maxErrorRateRle30 = std::pow(10.0, - minQRle30 / 10.0);

    

    // const double QThreshold = 1e-3;
    // const double logQThreshold = log(QThreshold);

    const double QThreshold6 = 1e-6;
    const double logQThreshold6 = log(QThreshold6);

    const double QThreshold3 = 1e-3;
    const double logQThreshold3 = log(QThreshold3);


    // const double QThreshold1 = 1e-14;
    // const double logQThreshold1 = log(QThreshold1);

    // const double QThreshold2 = 1e-12;
    // const double logQThreshold2 = log(QThreshold2);

    // const double QThreshold3 = 1e-10;
    // const double logQThreshold3 = log(QThreshold);

    // const double QThreshold4 = 1e-8;
    // const double logQThreshold4 = log(QThreshold4);

    // const double QThreshold5 = 1e-6;
    // const double logQThreshold5 = log(QThreshold5);

    // const double QThreshold6 = 1e-4;
    // const double logQThreshold6 = log(QThreshold6);

    // const double QThreshold7 = 1e-3;
    // const double logQThreshold7 = log(QThreshold7);

    // const double QThreshold8 = 1e-2;
    // const double logQThreshold8 = log(QThreshold8);


    const double QThresholdForBreaks = 1e-6;
    const double logQThresholdForBreaks = log(QThresholdForBreaks);

    //
    // 1. Order alignments in order of increasing Q. 
    //


    // Gather in alignmentTable[alignmentID, Q]
    // alignments in order of increasing Q.
    // Q(n) = (1 + δ/2ε)^n * e-δL
    // ε = 1e-4, δ = 5e-4
    // logQ(n) = αn - δL
    vector< pair<uint64_t, double> > alignmentTable1;
    vector< pair<uint64_t, double> > alignmentTable2;
    vector< pair<uint64_t, double> > alignmentTable3;
    vector< pair<uint64_t, double> > alignmentTable4;
    vector< pair<uint64_t, double> > alignmentTable5;
    vector< pair<uint64_t, double> > alignmentTable6;
    vector< pair<uint64_t, double> > alignmentTable7;
    vector< pair<uint64_t, double> > alignmentTable8;
    vector< pair<uint64_t, double> > alignmentTable9;
    vector< pair<uint64_t, double> > alignmentTable10;
    vector< pair<uint64_t, double> > alignmentTable11;
    vector< pair<uint64_t, double> > alignmentTable12;
    vector< pair<uint64_t, double> > alignmentTable13;
    vector< pair<uint64_t, double> > alignmentTable14;
    vector< pair<uint64_t, double> > alignmentTable15;
    vector< pair<uint64_t, double> > alignmentTable16;
    vector< pair<uint64_t, double> > alignmentTable17;
    vector< pair<uint64_t, double> > alignmentTable18;
    vector< pair<uint64_t, double> > alignmentTable19;
    vector< pair<uint64_t, double> > alignmentTable20;
    vector< pair<uint64_t, double> > alignmentTable21;
    vector< pair<uint64_t, double> > alignmentTable22;
    vector< pair<uint64_t, double> > alignmentTable23;
    vector< pair<uint64_t, double> > alignmentTable24;
    vector< pair<uint64_t, double> > alignmentTable25;
    vector< pair<uint64_t, double> > alignmentTable26;
    vector< pair<uint64_t, double> > alignmentTable27;
    vector< pair<uint64_t, double> > alignmentTable28;
    vector< pair<uint64_t, double> > alignmentTable29;
    vector< pair<uint64_t, double> > alignmentTable30;
    vector< pair<uint64_t, double> > alignmentTable31;
    vector< pair<uint64_t, double> > alignmentTable32;
    vector< pair<uint64_t, double> > alignmentTable33;
    vector< pair<uint64_t, double> > alignmentTable34;
    vector< pair<uint64_t, double> > alignmentTable35;
    vector< pair<uint64_t, double> > alignmentTable36;
    vector< pair<uint64_t, double> > alignmentTable37;
    vector< pair<uint64_t, double> > alignmentTable38;
    vector< pair<uint64_t, double> > alignmentTable39;
    vector< pair<uint64_t, double> > alignmentTable40;
    vector< pair<uint64_t, double> > alignmentTable41;
    vector< pair<uint64_t, double> > alignmentTable42;
    vector< pair<uint64_t, double> > alignmentTable43;
    vector< pair<uint64_t, double> > alignmentTable44;
    vector< pair<uint64_t, double> > alignmentTable45;
    vector< pair<uint64_t, double> > alignmentTable46;
    vector< pair<uint64_t, double> > alignmentTable47;
    vector< pair<uint64_t, double> > alignmentTable48;
    vector< pair<uint64_t, double> > alignmentTable49;
    vector< pair<uint64_t, double> > alignmentTable50;
    vector< pair<uint64_t, double> > alignmentTable51;
    vector< pair<uint64_t, double> > alignmentTable52;
    vector< pair<uint64_t, double> > alignmentTable53;
    vector< pair<uint64_t, double> > alignmentTable54;
    vector< pair<uint64_t, double> > alignmentTable55;
    vector< pair<uint64_t, double> > alignmentTable56;



    vector< pair<uint64_t, double> > alignmentTable;
    vector< pair<uint64_t, double> > alignmentTableNotPassFilter;
    vector< pair<uint64_t, double> > QAlignments;
    const double epsilon = 1e-4;
    const double delta = 5e-4;
    const double alpha = log(1 + delta/(2*epsilon));

    // Get stats about the reads
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    
    // Keep track of which readIds were used in alignments
    vector<bool> readUsed(readCount, false);

    // Flag alignments to be kept for breaks detections.
    vector<bool> keepAlignmentsForBreaks(alignmentCount, false);


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
        const double errorRateRle = thisAlignmentData.info.errorRateRle;
        const double markerCount = thisAlignmentData.info.markerCount;

        // logQ(n) = αn - δL
        const double logQ = alpha * double(n) - delta * L;

        // if(errorRateRle<= maxErrorRateRle) {
        //     QAlignments.push_back(make_pair(alignmentId, errorRateRle));
        //     alignmentTable.push_back(make_pair(alignmentId, logQ));
        // }

        // if (logQ <= logQThreshold and errorRateRle <= maxErrorRateRle) {

        // if (logQ <= logQThreshold and errorRateRle <= maxErrorRateRle) {
        //     alignmentTable.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThresholdForBreaks){
        //     alignmentTableNotPassFilter.push_back(make_pair(alignmentId, logQ));
        //     keepAlignmentsForBreaks[alignmentId] = true;
        // }




        if (logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 50000) {
            alignmentTable1.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 50000){
            alignmentTable2.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 50000){
            alignmentTable3.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 30000) {
            alignmentTable4.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 30000){
            alignmentTable5.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 30000){
            alignmentTable6.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 10000) {
            alignmentTable7.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 10000){
            alignmentTable8.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 10000){
            alignmentTable9.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 8000) {
            alignmentTable10.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 8000){
            alignmentTable11.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 8000){
            alignmentTable12.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 6000) {
            alignmentTable13.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 6000){
            alignmentTable14.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 6000){
            alignmentTable15.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 4000) {
            alignmentTable16.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 4000){
            alignmentTable17.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 4000){
            alignmentTable18.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if (logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 50000) {
            alignmentTable19.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 50000){
            alignmentTable20.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 50000){
            alignmentTable21.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 30000) {
            alignmentTable22.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 30000){
            alignmentTable23.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 30000){
            alignmentTable24.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 10000) {
            alignmentTable25.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 10000){
            alignmentTable26.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 10000){
            alignmentTable27.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 8000) {
            alignmentTable28.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 8000){
            alignmentTable29.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 8000){
            alignmentTable30.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 6000) {
            alignmentTable31.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 6000){
            alignmentTable32.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 6000){
            alignmentTable33.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 4000) {
            alignmentTable34.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 4000){
            alignmentTable35.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 4000){
            alignmentTable36.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60 and markerCount >= 1500) {
            alignmentTable37.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50 and markerCount >= 1500){
            alignmentTable38.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40 and markerCount >= 1500){
            alignmentTable39.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60 and markerCount >= 1500) {
            alignmentTable40.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50 and markerCount >= 1500){
            alignmentTable41.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40 and markerCount >= 1500){
            alignmentTable42.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle60) {
            alignmentTable43.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle50){
            alignmentTable44.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle40){
            alignmentTable45.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle30){
        //     alignmentTable46.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle60) {
            alignmentTable47.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle50){
            alignmentTable48.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle40){
            alignmentTable49.push_back(make_pair(alignmentId, logQ));
            readUsed[readId0] = true;
            readUsed[readId1] = true;
        }
        // } else if(logQ <= logQThreshold6 and errorRateRle <= maxErrorRateRle30){
        //     alignmentTable50.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold3 and errorRateRle <= maxErrorRateRle30){
        //     alignmentTable51.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // }





        // if (logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 50000) {
        //     alignmentTable1.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 50000){
        //     alignmentTable2.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 50000){
        //     alignmentTable3.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 50000){
        //     alignmentTable4.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 30000) {
        //     alignmentTable5.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 30000){
        //     alignmentTable6.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 30000){
        //     alignmentTable7.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 30000){
        //     alignmentTable8.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 10000) {
        //     alignmentTable9.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 10000){
        //     alignmentTable10.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 10000){
        //     alignmentTable11.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 10000){
        //     alignmentTable12.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 8000) {
        //     alignmentTable13.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 8000){
        //     alignmentTable14.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 8000){
        //     alignmentTable15.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 8000){
        //     alignmentTable16.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 6000) {
        //     alignmentTable17.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 6000){
        //     alignmentTable18.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 6000){
        //     alignmentTable19.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 6000){
        //     alignmentTable20.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 4000) {
        //     alignmentTable21.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 4000){
        //     alignmentTable22.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 4000){
        //     alignmentTable23.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 4000){
        //     alignmentTable24.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60 and markerCount >= 2000) {
        //     alignmentTable25.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50 and markerCount >= 2000){
        //     alignmentTable26.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40 and markerCount >= 2000){
        //     alignmentTable27.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30 and markerCount >= 2000){
        //     alignmentTable28.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle60) {
        //     alignmentTable29.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle50){
        //     alignmentTable30.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle40){
        //     alignmentTable31.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // } else if(logQ <= logQThreshold and errorRateRle <= maxErrorRateRle30){
        //     alignmentTable32.push_back(make_pair(alignmentId, logQ));
        //     readUsed[readId0] = true;
        //     readUsed[readId1] = true;
        // }

    }


    // Sort by increasing logQ
    sort(alignmentTable1.begin(), alignmentTable1.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable2.begin(), alignmentTable2.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable3.begin(), alignmentTable3.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable4.begin(), alignmentTable4.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable5.begin(), alignmentTable5.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable6.begin(), alignmentTable6.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable7.begin(), alignmentTable7.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable8.begin(), alignmentTable8.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable9.begin(), alignmentTable9.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable10.begin(), alignmentTable10.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable11.begin(), alignmentTable11.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable12.begin(), alignmentTable12.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable13.begin(), alignmentTable13.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable14.begin(), alignmentTable14.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable15.begin(), alignmentTable15.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable16.begin(), alignmentTable16.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable17.begin(), alignmentTable17.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable18.begin(), alignmentTable18.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable19.begin(), alignmentTable19.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable20.begin(), alignmentTable20.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable21.begin(), alignmentTable21.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable22.begin(), alignmentTable22.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable23.begin(), alignmentTable23.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable24.begin(), alignmentTable24.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable25.begin(), alignmentTable25.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable26.begin(), alignmentTable26.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable27.begin(), alignmentTable27.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable28.begin(), alignmentTable28.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable29.begin(), alignmentTable29.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable30.begin(), alignmentTable30.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable31.begin(), alignmentTable31.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable32.begin(), alignmentTable32.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable33.begin(), alignmentTable33.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable34.begin(), alignmentTable34.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable35.begin(), alignmentTable35.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable36.begin(), alignmentTable36.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable37.begin(), alignmentTable37.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable38.begin(), alignmentTable38.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable39.begin(), alignmentTable39.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable40.begin(), alignmentTable40.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable41.begin(), alignmentTable41.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable42.begin(), alignmentTable42.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable43.begin(), alignmentTable43.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable44.begin(), alignmentTable44.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable45.begin(), alignmentTable45.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable46.begin(), alignmentTable46.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable47.begin(), alignmentTable47.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable48.begin(), alignmentTable48.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable49.begin(), alignmentTable49.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable50.begin(), alignmentTable50.end(), OrderPairsBySecondOnly<uint64_t, double>());
    sort(alignmentTable51.begin(), alignmentTable51.end(), OrderPairsBySecondOnly<uint64_t, double>());
    // sort(alignmentTable52.begin(), alignmentTable52.end(), OrderPairsBySecondOnly<uint64_t, double>());
    // sort(alignmentTable53.begin(), alignmentTable53.end(), OrderPairsBySecondOnly<uint64_t, double>());
    // sort(alignmentTable54.begin(), alignmentTable54.end(), OrderPairsBySecondOnly<uint64_t, double>());
    // sort(alignmentTable55.begin(), alignmentTable55.end(), OrderPairsBySecondOnly<uint64_t, double>());
    // sort(alignmentTable56.begin(), alignmentTable56.end(), OrderPairsBySecondOnly<uint64_t, double>());

    

    //sort(alignmentTable.begin(), alignmentTable.end(), OrderPairsBySecondOnly<uint64_t, double>());
    //sort(alignmentTableNotPassFilter.begin(), alignmentTableNotPassFilter.end(), OrderPairsBySecondOnly<uint64_t, double>());
    // cout << "The alignmentTable has " << alignmentTable.size() << " entries." << endl; // 123863
    // cout << "The alignmentTableNotPassFilter has " << alignmentTableNotPassFilter.size() << " entries." << endl;

    

    
    


    // Create the read graph using all the alignments.
    // createReadGraphUsingAllAlignments(keepAlignmentsForBreaks);



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
    vector<uint64_t> verticesDegreeGoodSupport(orientedReadCount, 0);
    vector<uint64_t> verticesDegree30to40(orientedReadCount, 0);


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


    // Flag all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentCount, false);



    // Process alignments in order of increasing Q
    //vector alignmentTablesToProcess({alignmentTable1, alignmentTable2, alignmentTable3, alignmentTable4, alignmentTable5, alignmentTable6, alignmentTable7, alignmentTable8});
    vector alignmentTablesToProcess({
        alignmentTable1, alignmentTable2, alignmentTable3, alignmentTable4, alignmentTable5, alignmentTable6, alignmentTable7, alignmentTable8,
        alignmentTable9, alignmentTable10, alignmentTable11, alignmentTable12, alignmentTable13, alignmentTable14, alignmentTable15, alignmentTable16, 
        alignmentTable17, alignmentTable18, alignmentTable19, alignmentTable20, alignmentTable21, alignmentTable22, alignmentTable23, alignmentTable24,
        alignmentTable25, alignmentTable26, alignmentTable27, alignmentTable28, alignmentTable29, alignmentTable30, alignmentTable31, alignmentTable32,
        alignmentTable33, alignmentTable34, alignmentTable35, alignmentTable36, alignmentTable37, alignmentTable38, alignmentTable39, alignmentTable40,
        alignmentTable41, alignmentTable42, alignmentTable43, alignmentTable44, alignmentTable45, alignmentTable46, alignmentTable47, alignmentTable48,
        alignmentTable49, alignmentTable50, alignmentTable51, alignmentTable52, alignmentTable53, alignmentTable54, alignmentTable55, alignmentTable56
    });
    
    uint64_t crossStrandEdgeCount = 0;

    for (auto alignmentTableToProcess : alignmentTablesToProcess) {
        for(auto it=alignmentTableToProcess.begin(); it!=alignmentTableToProcess.end(); ++it) {
            const pair<uint64_t, double>& p = *it;
            const uint64_t alignmentId = p.first;
            const double logQ = p.second;

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
            const uint64_t degreeA0GoodSupport = verticesDegreeGoodSupport[A0.getValue()];
            const uint64_t degreeB0GoodSupport = verticesDegreeGoodSupport[B0.getValue()];
            const uint64_t degreeA030to40 = verticesDegree30to40[A0.getValue()];
            const uint64_t degreeB030to40 = verticesDegree30to40[B0.getValue()];


            // if(degreeA030to40 >= 1 and degreeB030to40 >= 1) {
            //     // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
            //     continue;
            // }

            if(degreeA0GoodSupport >= 20 and degreeB0GoodSupport >= 20) {
                // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
                continue;
            }


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
            
            if(alignment.info.markerCount >= 1500) {
                verticesDegreeGoodSupport[A0.getValue()]++;
                verticesDegreeGoodSupport[B0.getValue()]++;
                verticesDegreeGoodSupport[A1.getValue()]++;
                verticesDegreeGoodSupport[B1.getValue()]++;
            }

            if(alignment.info.errorRateRle <= maxErrorRateRle30 and alignment.info.errorRateRle >= maxErrorRateRle40) {
                verticesDegree30to40[A0.getValue()]++;
                verticesDegree30to40[B0.getValue()]++;
                verticesDegree30to40[A1.getValue()]++;
                verticesDegree30to40[B1.getValue()]++;
            }
            

            // Update disjoint sets
            disjointSets.union_set(a0, b0);
            disjointSets.union_set(a1, b1);

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
    }

    // Track the size of each set in the disjoint sets
    vector<uint64_t> setSizes(orientedReadCount, 0);
    for (int i = 0; i < orientedReadCount; ++i) {
        setSizes[disjointSets.find_set(i)]++;
    }


    // Print how many alignments were kept in this step
    const size_t keepCountR1 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Finding strict disjointSets step: Keeping " << keepCountR1 << " alignments of " << keepAlignment.size() << endl;















    //*
    //
    // Create the dynamically adjustable boost readGraph using the alignments we selected.
    //
    //*
    
    //createReadGraphUsingSelectedAlignments(keepAlignment);

    // Create the dynamically adjusted strictly fitlered readGraph.
    using boost::add_vertex;
    using boost::add_edge;

    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph4 readGraph(orientedReadCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // Record whether this alignment is used in the read graph.
        const bool keepThisAlignment = keepAlignment[alignmentId];
        const AlignmentData& alignment = alignmentData[alignmentId];

        // If this alignment is not used in the read graph, we are done.
        if(not keepThisAlignment) {
            continue;
        }

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(alignment.readIds[0], 0);
        OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4Edge(alignmentId), readGraph);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    }


    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph4AllAlignments readGraphAllAlignments(orientedReadCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // Record whether this alignment is used in the read graph.
        const bool keepThisAlignment = keepAlignmentsForBreaks[alignmentId];
        const AlignmentData& alignment = alignmentData[alignmentId];

        // If this alignment is not used in the read graph, we are done.
        if(not keepThisAlignment) {
            continue;
        }

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(alignment.readIds[0], 0);
        OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);
    }




    // vector<pair<uint64_t, size_t>> edgesToRemove;

    // // Keep only the top 10 neighbors with the highest markerCounts
    // BGL_FORALL_VERTICES(v, readGraph, ReadGraph4) {
    //     // Create a vector to store for each vertex the all of it's neighrbors markerCounts in distance 1
    //     vector<pair<uint64_t, double>> neighborsMarkerCount;
    //     vector<pair<uint64_t, double>> neighborsMarkerCount2;

    //     const ReadGraph4Vertex& vertex = readGraph[v];
    //     const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
        
    //     //cout << "Vertex " << orientedReadId << " has degree " << out_degree(v, readGraph) << endl;


    //     if (verticesDegree[orientedReadId.getValue()] <= 10) {
    //         continue;
    //     }

    //     uint64_t numberOfAlignmentsThatCanBeRemoved = verticesDegree[orientedReadId.getValue()] - 10;

    //     cout << "Vertex " << orientedReadId << " has degree " << out_degree(v, readGraph) << endl;
        

    //     BGL_FORALL_OUTEDGES(v, e, readGraph, ReadGraph4) {
    //         const uint32_t alignmentId = readGraph[e].alignmentId;
    //         const AlignmentData& alignment = alignmentData[alignmentId];
    //         const double errorRateRle = alignment.info.errorRateRle;
    //         const double markerCount = alignment.info.markerCount;

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         SHASTA_ASSERT(readId0 < readId1);
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         if (verticesDegree[B0.getValue()] <= 10) {
    //             continue;
    //         }

    //         neighborsMarkerCount.push_back(make_pair(target(e, readGraph), markerCount));
    //         neighborsMarkerCount2.push_back(make_pair(B0.getValue(), markerCount));
    //     }

    //     //cout << "Vertex " << orientedReadId << " has out degree " << out_degree(v, readGraph) << endl;
    //     //cout << "Vertex " << orientedReadId << " has in degree " << in_degree(v, readGraph) << endl;

    //     // Sort the neighbors by markerCount
    //     sort(neighborsMarkerCount.begin(), neighborsMarkerCount.end(), OrderPairsBySecondOnly<uint64_t, double>());
    //     sort(neighborsMarkerCount2.begin(), neighborsMarkerCount2.end(), OrderPairsBySecondOnly<uint64_t, double>());

    //     //cout << "Vertex " << orientedReadId << " has out degree " << out_degree(v, readGraph) << endl;
    //     //cout << "Vertex " << orientedReadId << " has in degree " << in_degree(v, readGraph) << endl;


    //     // Remove edges in the edgesToRemove vector
    //     for (uint32_t i = 0; i < numberOfAlignmentsThatCanBeRemoved; ++i) {
    //         edgesToRemove.push_back(make_pair(orientedReadId.getValue(), neighborsMarkerCount[i].first));
    //         //remove_edge(orientedReadId.getValue(), neighborsMarkerCount[i].first, readGraph);
    //         verticesDegree[orientedReadId.getValue()]--;
    //         verticesDegree[neighborsMarkerCount2[i].first]--;
    //         // Also decrease degree for reverse complement vertices
    //         OrientedReadId reverseOrientedReadId = orientedReadId;
    //         reverseOrientedReadId.flipStrand();
    //         OrientedReadId reverseNeighborId(OrientedReadId::fromValue(neighborsMarkerCount2[i].first));
    //         reverseNeighborId.flipStrand();
    //         verticesDegree[reverseOrientedReadId.getValue()]--;
    //         verticesDegree[reverseNeighborId.getValue()]--;
    //     }


    //     // cout << "Vertex " << orientedReadId << " has out degree " << out_degree(v, readGraph) << endl;
    //     //cout << "Vertex " << orientedReadId << " has in degree " << in_degree(v, readGraph) << endl;


    // }

    // cout << "Number of edges to remove: " << edgesToRemove.size() << endl;
    // // Remove them
    // for (uint32_t i = 0; i < edgesToRemove.size(); ++i) {
    //     //remove_edge(edgesToRemove[i], readGraph);
    //     remove_edge(edgesToRemove[i].first, edgesToRemove[i].second, readGraph);
    // }


    // // Initially, each alignment generates two edges.
    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

    //     // Record whether this alignment is used in the read graph.
    //     const bool keepThisAlignment = keepAlignment[alignmentId];
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // If this alignment is not used in the read graph, we are done.
    //     if(not keepThisAlignment) {
    //         continue;
    //     }

    //     // Get the OrientedReadIds.
    //     OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //     OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //     SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //     // Get the alignment data
    //     const ReadId readId0 = alignment.readIds[0];
    //     const ReadId readId1 = alignment.readIds[1];
    //     const bool isSameStrand = alignment.isSameStrand;
    //     SHASTA_ASSERT(readId0 < readId1);
    //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //     if (verticesDegree[A0.getValue()] < 10 or verticesDegree[B0.getValue()] < 10) {
    //         continue;
    //     }

        
        
    
    // }

    
    
    // // iterate over all alignments in readGraph and check if the alignment are in QAlignments
    
    // // Print the number of alignments in QAlignments
    // cout << "There are " << QAlignments.size() << " alignments in QAlignments." << endl;

    // // Now iterate over alignments in readGraph and check if they exist in QAlignments
    // for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        
    //     // Look for this alignmentId in QAlignments
    //     auto it = std::find_if(QAlignments.begin(), QAlignments.end(),
    //         [alignmentId](const pair<uint64_t, double>& p) { 
    //             return p.first == alignmentId; 
    //         });

    //     // If this alignment is not used in readGraph, and the alignment is in QAlignments, add the edge
    //     if(!alignmentData[alignmentId].info.isInReadGraph) {
    //         if(it != QAlignments.end()) {
    //             AlignmentData& alignment = alignmentData[alignmentId];
    //             // Get the OrientedReadIds.
    //             OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //             OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //             SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //             // Remove the edge.
    //             add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

    //             // Also remove the reverse complemented edge.
    //             orientedReadId0.flipStrand();
    //             orientedReadId1.flipStrand();
    //             add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph4AllAlignmentsEdge(alignmentId), readGraphAllAlignments);

    //             alignment.info.isInReadGraph = 1;
    //             keepAlignment[alignmentId] = true;

    //             continue;
    //         }
    //     }
        
        
    //     // // If this alignment is used in readGraph, and the alignment is NOT in QAlignments, remove the edge
    //     // if(it == QAlignments.end()) {
    //     //     AlignmentData& alignment = alignmentData[alignmentId];
    //     //     // Get the OrientedReadIds.
    //     //     OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    //     //     OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
    //     //     SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

    //     //     // Remove the edge.
    //     //     remove_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), readGraph);

    //     //     // Also remove the reverse complemented edge.
    //     //     orientedReadId0.flipStrand();
    //     //     orientedReadId1.flipStrand();
    //     //     remove_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), readGraph); 

    //     //     alignment.info.isInReadGraph = 0;
    //     //     keepAlignment[alignmentId] = false;
    //     // }
    // }


    
    







    // //*
    // //
    // // Now it is time to find the breaks in the readGraph
    // // In this case we add the alignments that did not pass the initial strict logQ filter
    // // involving reads that HAVE ALREADY BEEN USED (they are not isolated reads) AND are IN A DIFFERENT COMPONENT.
    // //
    // //*
    // vector<bool> endNodes(orientedReadCount, false);
    
    
    // bool weFoundAlignmentsInBreaksToAdd = true;

    // vector<bool> alignmentTableNotPassFilterAlreadyConsidered(alignmentTableNotPassFilter.size(), false);
    // const double QThresholdForAlignmentsInBreaks = 1e-10;
    // const double logQThresholdForAlignmentsInBreaks = log(QThresholdForAlignmentsInBreaks);

    // uint64_t allignmentsAdded = 0;
    // vector<vector<uint32_t>> betweenComponentsAlignments;
    // vector<pair<uint32_t, uint32_t>> componentsToUnite;


    // for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //     if((i % 100000) == 0) {
    //         cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //     }

    //     const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //     const double logQ = alignmentTableNotPassFilter[i].second;
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     const ReadId readId0 = alignment.readIds[0];
    //     const ReadId readId1 = alignment.readIds[1];
    //     const bool isSameStrand = alignment.isSameStrand;
    //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //     // Get the connected components that these oriented reads are in.
    //     const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //     const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //     const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //     const uint32_t b1 = disjointSets.find_set(B1.getValue());


    //     // Case 1: There is an alignment between two added reads
    //     // that passed the initial strict logQ filter and these reads are 
    //     // IN A DIFFERENT COMPONENT.
    //     // This case checks for direct connections between to nodes in a break
    //     // using alignments that did not pass the initial strict logQ filter.
    //     // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
    //     bool isInDifferentComponent = (a0 != b0) && (a0 != b1);
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && isInDifferentComponent && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
    //         // Find neighbors of the orientedReadId in the filtered readGraph.
    //         // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
    //         // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
    //         uint64_t maxDistance = 3;
    //         vector<OrientedReadId> readGraphNeighbors;
    //         readGraph.findNeighbors(A0, maxDistance, readGraphNeighbors);

    //         // Check if B0 is NOT in the neighbors
    //         if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //             cout << "Connecting components: b0 " << b0 << " and a0 " << a0 << endl;
    //             cout << "Components sizes are: b0 " << setSizes[b0] << " a0 " << setSizes[a0] << " b1 " << setSizes[b1] << " a1 " << setSizes[a1] << endl;
    //             cout << "Adding alignment of Read1: " << readId0 << " with the Read2: " << readId1 << endl;
    //             endNodes[A0.getValue()] = true;
    //             endNodes[B0.getValue()] = true;
    //             add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //             keepAlignment[alignmentId] = true;
    //             allignmentsAdded++;
    //             alignmentTableNotPassFilterAlreadyConsidered[i] = true;
    //             setSizes[a0] += setSizes[b0];
    //             setSizes[b0] = setSizes[a0];
    //             setSizes[a1] += setSizes[b1];
    //             setSizes[b1] = setSizes[a1];
    //             cout << "Components sizes are: b0 " << setSizes[b0] << " a0 " << setSizes[a0] << " b1 " << setSizes[b1] << " a1 " << setSizes[a1] << endl;
    //             componentsToUnite.push_back(make_pair(a0, b0));
    //             componentsToUnite.push_back(make_pair(a1, b1));
                
    //         }
    //         continue;
    //     }
    // }
    // for (auto& p : componentsToUnite) {
    //     disjointSets.union_set(p.first, p.second);
    // }

    // cout << "We added " << allignmentsAdded << " Alignments" << endl;







    // //*
    // //
    // // Now it is time to find the breaks in the readGraph and add the
    // // alignments between two added reads that passed the initial strict 
    // // logQ filter and these reads are IN A DIFFERENT COMPONENT.
    // //
    // //*
    // weFoundAlignmentsInBreaksToAdd = true;
    // while (weFoundAlignmentsInBreaksToAdd) {
    //     weFoundAlignmentsInBreaksToAdd = false;
    //     cout << "We added " << allignmentsAdded << " Alignments" << endl;
    //     for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //         if((i % 100000) == 0) {
    //             cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //         }

    //         const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //         const double logQ = alignmentTableNotPassFilter[i].second;
    //         const AlignmentData& alignment = alignmentData[alignmentId];

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // Case 1a: There is an alignment between two added reads
    //         // that passed the initial strict logQ filter and these reads are 
    //         // IN THE SAME CONNECTED COMPONENT.
    //         // This case checks for direct connections between to nodes in a break
    //         // using alignments that did not pass the initial strict logQ filter.
    //         // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
    //         if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && a0 == b0 && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
    //             // Find neighbors of the orientedReadId in the filtered readGraph.
    //             // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
    //             // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
                
    //             uint64_t maxDistance = 10;
    //             vector<OrientedReadId> readGraphNeighbors;
    //             readGraphAllAlignments.findNeighborsEarlyStopWhenReachSameComponentNode(A0, disjointSets, maxDistance, readGraphNeighbors);

    //             // Get the last node in the neighbors, which contains the first node in the same component as A0
    //             OrientedReadId lastNode = readGraphNeighbors.back();

    //             // Verify it is in the same component as A0
    //             if(disjointSets.find_set(lastNode.getValue()) == a0) {
    //                 SHASTA_ASSERT(a1 == b0);
    //             }

    //             // Check if it is in the same component as A0
    //             if(disjointSets.find_set(lastNode.getValue()) != a0) {
    //                 // We did not found a node in the same component as A0
    //                 // Tha means that this is probably not a break situation
    //                 // Or that we did not search further enough
    //                 continue;

    //             }

    //             // Check if it is in the same component as A0
                
    //             endNodes[A0.getValue()] = true;
    //             endNodes[B0.getValue()] = true;
    //             add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //             keepAlignment[alignmentId] = true;
    //             weFoundAlignmentsInBreaksToAdd = true;
    //             allignmentsAdded++;
    //             alignmentTableNotPassFilterAlreadyConsidered[i] = true;
                

    //             // 


    //             // Check if B0 is NOT in the neighbors
    //             if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //                 endNodes[A0.getValue()] = true;
    //                 endNodes[B0.getValue()] = true;
    //                 add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
    //                 keepAlignment[alignmentId] = true;
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //                 allignmentsAdded++;
    //                 alignmentTableNotPassFilterAlreadyConsidered[i] = true;
    //             }
    //             continue;
    //         }
    //     }
    // }






    //*
    //
    // Now it is time to find the breaks in the readGraph and add the
    // alignments between two added reads that passed the initial strict 
    // logQ filter and these reads are IN A DIFFERENT COMPONENT.
    //
    //*
    //vector<bool> endNodes(orientedReadCount, false);
    //vector<bool> alignmentTableNotPassFilterAlreadyConsidered(alignmentTableNotPassFilter.size(), false);
    //bool weFoundAlignmentsInBreaksToAdd = true;

    //const double QThresholdForAlignmentsInBreaks = 1e-10;
    //const double logQThresholdForAlignmentsInBreaks = log(QThresholdForAlignmentsInBreaks);

    //uint64_t allignmentsAdded = 0;

    // while (weFoundAlignmentsInBreaksToAdd) {
    //     weFoundAlignmentsInBreaksToAdd = false;
    //     cout << "We added " << allignmentsAdded << " Alignments" << endl;
    //     for(size_t i=0; i<alignmentTableNotPassFilter.size(); i++) {

    //         if((i % 100000) == 0) {
    //             cout << timestamp << i << "/" << alignmentTableNotPassFilter.size() << endl;
    //         }

    //         const uint64_t alignmentId = alignmentTableNotPassFilter[i].first;
    //         const double logQ = alignmentTableNotPassFilter[i].second;
    //         const AlignmentData& alignment = alignmentData[alignmentId];

    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

            // // Case 1a: There is an alignment between two added reads
            // // that passed the initial strict logQ filter and these reads are 
            // // IN THE SAME CONNECTED COMPONENT.
            // // This case checks for direct connections between to nodes in a break
            // // using alignments that did not pass the initial strict logQ filter.
            // // TODO: check if these alignments introduce bad connections (logQThresholdForAlignmentsInBreaks)
            // if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && a0 == b0 && logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
            //     // Find neighbors of the orientedReadId in the filtered readGraph.
            //     // We know there is a distance 1 connection between the 2 reads in the readGraphAllAlignments.
            //     // If there is no connection in the filtered readGraph in maxDistance apart, we can add the alignment.
            //     uint64_t maxDistance = 5;
            //     vector<OrientedReadId> readGraphNeighbors;
            //     readGraph.findNeighbors(A0, maxDistance, readGraphNeighbors);

            //     // Check if B0 is NOT in the neighbors
            //     if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
            //         endNodes[A0.getValue()] = true;
            //         endNodes[B0.getValue()] = true;
            //         add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
            //         keepAlignment[alignmentId] = true;
            //         weFoundAlignmentsInBreaksToAdd = true;
            //         allignmentsAdded++;
            //         alignmentTableNotPassFilterAlreadyConsidered[i] = true;
            //     }
            //     continue;
            // }

            

            // // Case 2: There is an alignment between an originaly used read that was in an alignment
            // // that passed the initial strict logQ filter, and an isolated read that was not used.
            // // if(readUsed[alignment.readIds[0]] && !readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks) {
            // if(readUsed[alignment.readIds[0]] && !readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
            //     // Find neighbors of the orientedReadId B0 (the isolated one) in the allAlignments readGraph.
            //     // Next, check which of these neighbors are in the same connected component as A0.
            //     // Next, check if these neighbors of B0 in the same connected component of A0 are reachable 
            //     // in a set distance from A0 in the filtered readGraph.
            //     // If some of them are not, we can add the alignment because we are in a break.
            //     // TODO: The maxDistances need to be adjusted probably.
            //     uint64_t maxDistance = 8;
            //     vector<OrientedReadId> readGraphAllAlignmentsNeighbors;
            //     readGraphAllAlignments.findNeighbors(B0, maxDistance, readGraphAllAlignmentsNeighbors);
                
            //     maxDistance = 12;
            //     vector<OrientedReadId> readGraphNeighbors;
            //     readGraph.findNeighbors(A0, maxDistance, readGraphNeighbors);

            //     // Check which of these neighbors of B0 are in the same connected component as A0.
            //     vector<OrientedReadId> neighborsInSameComponent;
            //     for(const OrientedReadId& neighbor: readGraphAllAlignmentsNeighbors) {
            //         if(disjointSets.find_set(neighbor.getValue()) == a0) {
            //             neighborsInSameComponent.push_back(neighbor);
            //         }
            //     }

            //     // Check if these neighbors of B0 in the same component of A0 are reachable from A0 in the filtered readGraph.
            //     // If any of them are not, we can add the alignment.
            //     bool weCanAddAlignment = false;
            //     for(const OrientedReadId& neighbor: neighborsInSameComponent) {
            //         if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), neighbor) == readGraphNeighbors.end()) {
            //             weCanAddAlignment = true;
            //             break;
            //         }
            //     }

            //     if(weCanAddAlignment) {
            //         endNodes[A0.getValue()] = true;
            //         endNodes[B0.getValue()] = true;
            //         add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
            //         keepAlignment[alignmentId] = true;
            //         readUsed[alignment.readIds[1]] = true;
            //         // Update disjoint sets
            //         disjointSets.union_set(a0, b0);
            //         disjointSets.union_set(a1, b1);
            //         weFoundAlignmentsInBreaksToAdd = true;
            //         allignmentsAdded++;
            //         alignmentTableNotPassFilterAlreadyConsidered[i] = true;
            //     }

            //     continue;
            // }
            
            // // Case 3: There is an alignment between an isolated read that was not used and
            // // an originaly used read that was in an alignment that passed the initial strict logQ filter.
            // // if(!readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks) {
            // if(!readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]] and logQ <= logQThresholdForAlignmentsInBreaks && alignmentTableNotPassFilterAlreadyConsidered[i] == false) {
            //     // Find neighbors of the orientedReadId A0 (the isolated one) in the allAlignments readGraph.
            //     // Next, check which of these neighbors are in the same connected component as B0.
            //     // Next, check if these neighbors of A0 in the same connected component of B0 are reachable 
            //     // in a set distance from B0 in the filtered readGraph.
            //     // If some of them are not, we can add the alignment because we are in a break.

            //     uint64_t maxDistance = 8;
            //     vector<OrientedReadId> readGraphAllAlignmentsNeighbors;
            //     readGraphAllAlignments.findNeighbors(A0, maxDistance, readGraphAllAlignmentsNeighbors);

            //     maxDistance = 12;
            //     vector<OrientedReadId> readGraphNeighbors;
            //     readGraph.findNeighbors(B0, maxDistance, readGraphNeighbors);

            //     // Check which of these neighbors of A0 are in the same connected component as B0.
            //     vector<OrientedReadId> neighborsInSameComponent;
            //     for(const OrientedReadId& neighbor: readGraphAllAlignmentsNeighbors) {
            //         if(disjointSets.find_set(neighbor.getValue()) == b0) {
            //             neighborsInSameComponent.push_back(neighbor);
            //         }
            //     }

            //     // Check if these neighbors of A0 in the same component of B0 are reachable from B0 in the filtered readGraph.
            //     // If any of them are not, we can add the alignment.
            //     bool weCanAddAlignment = false;
            //     for(const OrientedReadId& neighbor: neighborsInSameComponent) {
            //         if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), neighbor) == readGraphNeighbors.end()) {
            //             weCanAddAlignment = true;
            //             break;
            //         }
            //     }

            //     if(weCanAddAlignment) {
            //         endNodes[A0.getValue()] = true;
            //         endNodes[B0.getValue()] = true;
            //         add_edge(A0.getValue(), B0.getValue(), ReadGraph4Edge(alignmentId), readGraph);
            //         keepAlignment[alignmentId] = true;
            //         readUsed[alignment.readIds[0]] = true;
            //         // Update disjoint sets
            //         disjointSets.union_set(a0, b0);
            //         disjointSets.union_set(a1, b1);
            //         weFoundAlignmentsInBreaksToAdd = true;
            //         allignmentsAdded++;
            //         alignmentTableNotPassFilterAlreadyConsidered[i] = true;
            //     }

            //     continue;
            // }
    //     }
    // }


    // Remove previously created graphs
    //readGraphAllAlignments.remove();












    // // In this step now we will add additional alignments that were filterd by
    // // the Q filter and that involve isolated oriented reads that are added (readAdded) 
    // // and they belong to the same connected component.
    // const double QThresholdForFinalRoundAddedAlignments = 1e-4;
    // const double logQThresholdForFinalRoundAddedAlignments = log(QThresholdForFinalRoundAddedAlignments);
    // bool weFoundAlignmentsInBreaksToAdd = true;
    // // while (weFoundAlignmentsInBreaksToAdd) {
    // //     weFoundAlignmentsInBreaksToAdd = false;
    // for(size_t i=0; i<alignmentTableNotPassFilterNotIncludedAlignments.size(); i++) {
    //     if((i % 10000) == 0) {
    //         cout << timestamp << i << "/" << alignmentTableNotPassFilterNotIncludedAlignments.size() << endl;
    //     }

    //     const uint64_t alignmentId = alignmentTableNotPassFilterNotIncludedAlignments[i].first;
    //     const double logQ = alignmentTableNotPassFilterNotIncludedAlignments[i].second;
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // Case 1: There is an alignment between two added reads
    //     if(readAdded[alignment.readIds[0]] && readAdded[alignment.readIds[1]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 //cout << "Case 1 happened" << endl;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //             }
    //         }
    //         continue;
    //     }

    //     // Case 2a: There is an alignment between an added read and a originaly used read
    //     if(readAdded[alignment.readIds[0]] && readUsed[alignment.readIds[1]] && !readAdded[alignment.readIds[1]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 //cout << "Case 2a happened" << endl;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //             }
    //         }
    //         continue;
    //     }

    //     // Case 2b: There is an alignment between an added read and a originaly used read
    //     if(readAdded[alignment.readIds[1]] && readUsed[alignment.readIds[0]] && !readAdded[alignment.readIds[0]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 //cout << "Case 2b happened" << endl;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //                 weFoundAlignmentsInBreaksToAdd = true;
    //             }
    //         }
    //         continue;
    //     }



    //     // Case 3: There is an alignment between two originaly used read
    //     // We want to add only alignments that connect breaks in the read graph.
    //     // To do that we will check if the two oriented reads are in the same connected component
    //     // and if their distance with BFS is greater than 1 given that we found an alignment that connects them
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 cout << "Case 3a happened" << endl;

    //                 // Find neighbors of the orientedReadId
    //                 uint64_t maxDistance = 3;

    //                 vector<OrientedReadId> readGraphNeighbors;
    //                 readGraphAllAlignments.findNeighbors(A0, maxDistance, readGraphNeighbors);

    //                 // Check if B0 is NOT in the neighbors
    //                 if(find(readGraphNeighbors.begin(), readGraphNeighbors.end(), B0) == readGraphNeighbors.end()) {
    //                     cout << "Case 3b happened" << endl;
    //                     keepAlignment[alignmentId] = true;
    //                     // Update disjoint sets
    //                     disjointSets.union_set(a0, b0);
    //                     disjointSets.union_set(a1, b1);
    //                     weFoundAlignmentsInBreaksToAdd = true;
    //                 }
                    
    //             }
                
    //         }
            
    //     }
    // }



    
    


    // // In this step now we will add additional alignments that were filterd by
    // // the Q filter and that involve isolated oriented reads that are added (readAdded) 
    // // and they belong to the same connected component.
    // const double QThresholdForFinalRoundAddedAlignments = 1e-5;
    // const double logQThresholdForFinalRoundAddedAlignments = log(QThresholdForFinalRoundAddedAlignments);
    // for(size_t i=0; i<alignmentTableNotPassFilterNotIncludedAlignments.size(); i++) {
    //     const uint64_t alignmentId = alignmentTableNotPassFilterNotIncludedAlignments[i].first;
    //     const double logQ = alignmentTableNotPassFilterNotIncludedAlignments[i].second;
    //     const AlignmentData& alignment = alignmentData[alignmentId];

    //     // check if the reads were used in alignments
    //     if(readUsed[alignment.readIds[0]] && readUsed[alignment.readIds[1]]) {
    //         // if(readAdded[alignment.readIds[0]] && readAdded[alignment.readIds[1]]) {
    //         //     const ReadId readId0 = alignment.readIds[0];
    //         //     const ReadId readId1 = alignment.readIds[1];
    //         //     const bool isSameStrand = alignment.isSameStrand;
    //         //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         //     // Get the connected components that these oriented reads are in.
    //         //     const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         //     const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         //     const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         //     const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         //     // check if the reads are in the same connected component
    //         //     if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //         //         keepAlignment[alignmentId] = true;
    //         //         // Update disjoint sets
    //         //         disjointSets.union_set(a0, b0);
    //         //         disjointSets.union_set(a1, b1);
    //         //     }
    //         // }
    //         const ReadId readId0 = alignment.readIds[0];
    //         const ReadId readId1 = alignment.readIds[1];
    //         const bool isSameStrand = alignment.isSameStrand;
    //         const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //         const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //         const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //         const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //         // Get the connected components that these oriented reads are in.
    //         const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //         const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //         const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //         const uint32_t b1 = disjointSets.find_set(B1.getValue());

    //         // check if the reads are in the same connected component
    //         if(disjointSets.find_set(A0.getValue()) == disjointSets.find_set(B0.getValue())) {
    //             if(logQ <= logQThresholdForFinalRoundAddedAlignments) {
    //                 // check the degree of the vertices
    //                 const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //                 const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //                 // if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //                 //     continue;
    //                 // }
    //                 // if(alignment.info.markerCount <= 1500) {
    //                 //     continue;
    //                 // }
    //                 keepAlignment[alignmentId] = true;
    //                 // Update disjoint sets
    //                 disjointSets.union_set(a0, b0);
    //                 disjointSets.union_set(a1, b1);
    //             }
                
    //         }
            
    //     }
    // }

    // Print how many alignments were kept
    const size_t keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Adding isolated reads step: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;
    

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);


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


    // // We iterate over each conected component we have constructed
    // // using strict disjointSets and we try to detect breaks.
    // // If there is a break in the readGraph and there is not a break in the AllAlignments graph
    // // we need to add the alignments that connect the two oriented reads to the readGraph.
    // for(auto it=componentMap.begin(); it!=componentMap.end(); ++it) {
    //     const vector<OrientedReadId>& orientedReadIds = it->second;
    //     for(const OrientedReadId& orientedReadId: orientedReadIds) {
    //         // Find neighbors of the orientedReadId
    //         uint64_t maxDistance = 5;

    //         vector<OrientedReadId> readGraphNeighbors;
    //         readGraph.findNeighbors(orientedReadId, maxDistance, readGraphNeighbors);

    //         vector<OrientedReadId> readGraphAllAlignmentsNeighbors;
    //         readGraphAllAlignments.findNeighbors(orientedReadId, maxDistance, readGraphAllAlignmentsNeighbors);

    //         // Find the neighbors that exist in the AllAlignments graph but not in the readGraph and also
    //         // they are in the disjointSets of the orientedReadId
    //         for(const OrientedReadId& orientedReadId1: readGraphAllAlignmentsNeighbors) {
    //             if(disjointSets.find_set(orientedReadId1.getValue()) == disjointSets.find_set(orientedReadId.getValue())) {

    //                 // Check if the distance between the OrientedReadIs are significant different in the readGraph and AllAlignments graph
    //                 // Find the shortest path between the two oriented reads
    //                 vector<uint32_t> distanceReadGraph(2*readCount, ReadGraph::infiniteDistance);
    //                 vector<OrientedReadId> reachedVerticesReadGraph;
    //                 vector<uint32_t> parentEdgesReadGraph(2*readCount);
    //                 vector<uint32_t> shortestPathReadGraph;
    //                 readGraph.computeShortPath(orientedReadId, orientedReadId1,
    //                     maxDistance, shortestPathReadGraph,
    //                     distanceReadGraph, reachedVerticesReadGraph, parentEdgesReadGraph);

    //                 vector<uint32_t> distanceReadGraphAllAlignments(2*readCount, ReadGraph::infiniteDistance);
    //                 vector<OrientedReadId> reachedVerticesReadGraphAllAlignments;
    //                 vector<uint32_t> parentEdgesReadGraphAllAlignments(2*readCount);
    //                 vector<uint32_t> shortestPathReadGraphAllAlignments;
    //                 readGraphAllAlignments.computeShortPath(orientedReadId, orientedReadId1,
    //                     maxDistance, shortestPathReadGraphAllAlignments,
    //                     distanceReadGraphAllAlignments, reachedVerticesReadGraphAllAlignments, parentEdgesReadGraphAllAlignments);

    //                 // If the distance between the OrientedReadIs are about the same then we can skip
    //                 // It means that there is not really a break in the readGraph. If there was a break
    //                 // we would have a significant difference in the distance between the two OrientedReadIs
    //                 // because a way bigger path would be needed to connect them.
    //                 uint64_t maxDistance = 5;
    //                 if (distanceReadGraph[orientedReadId1.getValue()] != ReadGraph::infiniteDistance && 
    //                     distanceReadGraphAllAlignments[orientedReadId1.getValue()] != ReadGraph::infiniteDistance) {
    //                     if(std::abs(int(distanceReadGraph[orientedReadId1.getValue()] - distanceReadGraphAllAlignments[orientedReadId1.getValue()])) < maxDistance) {
    //                         continue;
    //                     }

    //                 }

    //                 // Else, we need to add the alignments that connect the two oriented reads
    //                 // and add them to the readGraph if they do not break the disjoint sets introducing cross-strand edges
    //                 // Find all alignments that connect orientedReadId and orientedReadId1 even if they do not connect them directly
    //                 vector<uint32_t> alignmentsBetweenReads;

                    


    //                 if(std::find(readGraphNeighbors.begin(), readGraphNeighbors.end(), orientedReadId1) == readGraphNeighbors.end()) {
    //                     // Find all alignments that connect orientedReadId and orientedReadId1
    //                     vector<uint32_t> alignmentsBetweenReads;
    //                     for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //                         const AlignmentData& alignment = alignmentData[alignmentId];
    //                         if(alignment.info.isInReadGraph) {
    //                             continue;  // Skip alignments already in read graph
    //                         }
                            
    //                         const ReadId readId0 = alignment.readIds[0];
    //                         const ReadId readId1 = alignment.readIds[1];
    //                         const bool isSameStrand = alignment.isSameStrand;
                            
    //                         OrientedReadId readId0Oriented(readId0, 0);
    //                         OrientedReadId readId1Oriented(readId1, isSameStrand ? 0 : 1);

    //                         // Check if this alignment connects our two reads
    //                         if((readId0Oriented == orientedReadId && readId1Oriented == orientedReadId1) ||
    //                            (readId0Oriented == orientedReadId1 && readId1Oriented == orientedReadId)) {
    //                             alignmentsBetweenReads.push_back(alignmentId);
    //                         }
    //                     }

    //                     // Add these alignments to the read graph
    //                     for(uint32_t alignmentId : alignmentsBetweenReads) {
    //                         AlignmentData& alignment = alignmentData[alignmentId];
    //                         alignment.info.isInReadGraph = true;
    //                         readGraph.addEdge(alignment.readIds[0], alignment.readIds[1], 
    //                                          alignment.isSameStrand, alignmentId);
    //                     }
    //                 }
    //             }
    //         }
            


    //         // iterate over the neighbors
    //         for(const OrientedReadId& orientedReadId1: readGraphNeighbors) {
    //             // Find the shortest path between the two oriented reads
    //             vector<uint32_t> distanceReadGraph(2*readCount, ReadGraph::infiniteDistance);
    //             vector<OrientedReadId> reachedVerticesReadGraph;
    //             vector<uint32_t> parentEdgesReadGraph(2*readCount);
    //             vector<uint32_t> shortestPathReadGraph;
    //             readGraph.computeShortPath(orientedReadId, orientedReadId1,
    //                 maxDistance, shortestPathReadGraph,
    //                 distanceReadGraph, reachedVerticesReadGraph, parentEdgesReadGraph);

    //             vector<uint32_t> distanceReadGraphAllAlignments(2*readCount, ReadGraph::infiniteDistance);
    //             vector<OrientedReadId> reachedVerticesReadGraphAllAlignments;
    //             vector<uint32_t> parentEdgesReadGraphAllAlignments(2*readCount);
    //             vector<uint32_t> shortestPathReadGraphAllAlignments;
    //             readGraph.computeShortPath(orientedReadId, orientedReadId1,
    //                 maxDistance, shortestPathReadGraphAllAlignments,
    //                 distanceReadGraphAllAlignments, reachedVerticesReadGraphAllAlignments, parentEdgesReadGraphAllAlignments);

    //             // If the sistance between the OrientedReadIs are about the same then we can skip
    //             if (distanceReadGraph[orientedReadId1.getValue()] != ReadGraph::infiniteDistance && 
    //                 distanceReadGraphAllAlignments[orientedReadId1.getValue()] != ReadGraph::infiniteDistance) {
    //                 if(std::abs(int(distanceReadGraph[orientedReadId1.getValue()] - distanceReadGraphAllAlignments[orientedReadId1.getValue()])) < 5) {
    //                     continue;
    //                 }

    //             }

    //             // Iterate over all reach vertices in the AllAlignments graph and check if they
    //             // are in the same set of the disjointSets
    //             for(const OrientedReadId& reachedVertex: reachedVerticesReadGraphAllAlignments) {
    //                 if(disjointSets.find_set(reachedVertex.getValue()) == disjointSets.find_set(orientedReadId.getValue())) {
    //                     // now check if it is reachable from the readGraph
    //                     if(distanceReadGraph[reachedVertex.getValue()] == ReadGraph::infiniteDistance) {
    //                         // It is not reachable from the readGraph
    //                         // This means we need to find all the alignments that connect the two oriented reads
    //                         // and add them to the readGraph.
    //                         // Find all alignments that connect orientedReadId and reachedVertex
    //                         vector<uint32_t> alignmentsBetweenReads;
    //                         for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
    //                             const AlignmentData& alignment = alignmentData[alignmentId];
    //                             if(alignment.info.isInReadGraph) {
    //                                 continue;  // Skip alignments already in read graph
    //                             }
                                
    //                             const ReadId readId0 = alignment.readIds[0];
    //                             const ReadId readId1 = alignment.readIds[1];
    //                             const bool isSameStrand = alignment.isSameStrand;
                                
    //                             OrientedReadId readId0Oriented(readId0, 0);
    //                             OrientedReadId readId1Oriented(readId1, isSameStrand ? 0 : 1);

    //                             // Check if this alignment connects our two reads
    //                             if((readId0Oriented == orientedReadId && readId1Oriented == reachedVertex) ||
    //                                (readId0Oriented == reachedVertex && readId1Oriented == orientedReadId)) {
    //                                 alignmentsBetweenReads.push_back(alignmentId);
    //                             }
    //                         }

    //                         // Add these alignments to the read graph
    //                         for(uint32_t alignmentId : alignmentsBetweenReads) {
    //                             AlignmentData& alignment = alignmentData[alignmentId];
    //                             alignment.info.isInReadGraph = true;
    //                             readGraph.addEdge(alignment.readIds[0], alignment.readIds[1], 
    //                                              alignment.isSameStrand, alignmentId);
    //                         }


    //                     }
                        
    //                 }
    //             }

    //             // Check if the distance is the same
    //             if(distanceReadGraph[orientedReadId1.getValue()] == ReadGraph::infiniteDistance && distanceReadGraphAllAlignments[orientedReadId1.getValue()]) {
    //             }

                
    //         }
            


    //     }

    // }













    // // Process notIncludedReadsToAlignments in order of increasing Q
    // // but now use the already constructed disjointSets as a "good starting point"
    // uint64_t crossStrandEdgeCountR2 = 0;
    // for(auto it=notIncludedReadsToAlignments.begin(); it!=notIncludedReadsToAlignments.end(); ++it) {
    //     const pair<uint64_t, double>& p = *it;
    //     const uint64_t alignmentId = p.first;

    //     // Get the alignment data
    //     AlignmentData& alignment = alignmentData[alignmentId];
    //     const ReadId readId0 = alignment.readIds[0];
    //     const ReadId readId1 = alignment.readIds[1];
    //     const bool isSameStrand = alignment.isSameStrand;
    //     SHASTA_ASSERT(readId0 < readId1);
    //     const OrientedReadId A0 = OrientedReadId(readId0, 0);
    //     const OrientedReadId B0 = OrientedReadId(readId1, isSameStrand ? 0 : 1);
    //     const OrientedReadId A1 = OrientedReadId(readId0, 1);
    //     const OrientedReadId B1 = OrientedReadId(readId1, isSameStrand ? 1 : 0);

    //     SHASTA_ASSERT(A0.getReadId() == A1.getReadId());
    //     SHASTA_ASSERT(B0.getReadId() == B1.getReadId());
    //     SHASTA_ASSERT(A0.getStrand() == 1 - A1.getStrand());
    //     SHASTA_ASSERT(B0.getStrand() == 1 - B1.getStrand());

    //     // Get the connected components that these oriented reads are in.
    //     const uint32_t a0 = disjointSets.find_set(A0.getValue());
    //     const uint32_t b0 = disjointSets.find_set(B0.getValue());
    //     const uint32_t a1 = disjointSets.find_set(A1.getValue());
    //     const uint32_t b1 = disjointSets.find_set(B1.getValue());


    //     // If the alignment breaks strand separation, it is skipped.
    //     // If A0 and B1 are in the same connected component,
    //     // A1 and B0 also must be in the same connected component.
    //     // Adding this pair of edges would create a self-complementary
    //     // connected component containing A0, B0, A1, and B1,
    //     // and to ensure strand separation we don't want to do that.
    //     // So we mark these edges as cross-strand edges
    //     // and don't use them to update the disjoint set data structure.
    //     if(a0 == b1) {
    //         SHASTA_ASSERT(a1 == b0);
    //         crossStrandEdgeCountR2 += 2;
    //         continue;
    //     }

    //     // If both vertices of the potential edge have at least the required minimum number 
    //     // of neighbors, the alignment is also skipped. 
    //     const uint64_t degreeA0 = verticesDegree[A0.getValue()];
    //     const uint64_t degreeB0 = verticesDegree[B0.getValue()];
    //     if(degreeA0 >= maxAlignmentCount && degreeB0 >= maxAlignmentCount) {
    //         // cout << "Skipping alignment " << alignmentId << " because vertex " << A0.getValue() << " has degree " << degreeA0 << " and vertex " << B0.getValue() << " has degree " << degreeB0 << endl;
    //         continue;
    //     }

    //     // Add the alignment to the read graph.
    //     keepAlignment[alignmentId] = true;
    //     alignment.info.isInReadGraph = 1;

    //     // Update vertex degrees
    //     verticesDegree[A0.getValue()]++;
    //     verticesDegree[B0.getValue()]++;
    //     verticesDegree[A1.getValue()]++;
    //     verticesDegree[B1.getValue()]++;
        

    //     // Update disjoint sets
    //     disjointSets.union_set(a0, b0);
    //     disjointSets.union_set(a1, b1);
    // }

    // // Verify that for any read the two oriented reads are in distinct
    // // connected components.
    // for(ReadId readId=0; readId<readCount; readId++) {
    //     const OrientedReadId orientedReadId0(readId, 0);
    //     const OrientedReadId orientedReadId1(readId, 1);
    //     SHASTA_ASSERT(
    //         disjointSets.find_set(orientedReadId0.getValue()) !=
    //         disjointSets.find_set(orientedReadId1.getValue())
    //     );
    // }


    // // Print how many alignments were kept
    // const size_t keepCountR2 = count(keepAlignment.begin(), keepAlignment.end(), true);
    // cout << "Adding isolated reads step: Keeping " << keepCountR2 << " alignments of " << keepAlignment.size() << endl;

    

    cout << timestamp << "Done processing alignments." << endl;


    

    cout << timestamp << "createReadGraph4 with strand separation ends." << endl;

    cout << "Strand separation flagged " << crossStrandEdgeCount <<
        " read graph edges out of " << num_edges(readGraph) << " total in round 1." << endl;
    
    // cout << "Strand separation flagged " << crossStrandEdgeCountR2 <<
    //     " read graph edges out of " << readGraph.edges.size() << " total in round 2." << endl;


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

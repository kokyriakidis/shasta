// Shasta.
#include "mode3a-BubbleCleaner.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "MarkerGraph.hpp"
#include "invalid.hpp"
#include "mode3a.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
using namespace shasta;using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "utility.hpp"



BubbleCleaner::BubbleCleaner(
    const PackedMarkerGraph& packedMarkerGraph) :
    packedMarkerGraph(packedMarkerGraph)
{
    BubbleCleaner& bubbleCleaner = *this;

    // Construct an edge for each segment of the PackedMarkerGraph.
    // Add vertices as we go.
    for(uint64_t segmentId=0; segmentId<packedMarkerGraph.segments.size(); segmentId++) {
        const uint64_t vertexId0 = packedMarkerGraph.getFirstSegmentVertex(segmentId);
        const uint64_t vertexId1 = packedMarkerGraph.getLastSegmentVertex(segmentId);
        const vertex_descriptor v0 = getVertex(vertexId0);
        const vertex_descriptor v1 = getVertex(vertexId1);
        add_edge(v0, v1, BubbleCleanerEdge(segmentId), bubbleCleaner);
    }

}



string BubbleCleanerEdge::representation() const
{
    string s;
    for(uint64_t i=0; i<segments.size(); i++) {
        if(i != 0) {
            s += "-";
        }
        s += to_string(segments[i]);
    }
    return s;
}



// Compute average marker graph edge coverage for an edge.
double BubbleCleaner::averageEdgeCoverage(edge_descriptor e) const
{
    const BubbleCleanerEdge& edge = (*this)[e];

    uint64_t sum = 0;
    uint64_t n = 0;
    for(const uint64_t segmentId: edge.segments) {
        const auto segment = packedMarkerGraph.segments[segmentId];
        for(const uint64_t edgeId: segment) {
            sum += packedMarkerGraph.markerGraph.edgeCoverage(edgeId);
            n++;
        }
    }
    return double(sum) / double(n);
}



void BubbleCleanerEdge::assembledSequence(
    const PackedMarkerGraph& packedMarkerGraph,
    vector<Base>& sequence) const
{
    sequence.clear();
    for(const uint64_t segmentId: segments) {
        const span<const Base> segmentSequence =
            packedMarkerGraph.segmentClippedSequence(segmentId);
        copy(segmentSequence.begin(), segmentSequence.end(), back_inserter(sequence));

    }
}





// Get the vertex corresponding to a given marker graph vertex,
// creating if necessary.
BubbleCleaner::vertex_descriptor BubbleCleaner::getVertex(uint64_t markerGraphVertexId)
{
    BubbleCleaner& bubbleCleaner = *this;

    auto it = vertexMap.find(markerGraphVertexId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(BubbleCleanerVertex(markerGraphVertexId), bubbleCleaner);
        vertexMap.insert(make_pair(markerGraphVertexId, v));
        return v;
    } else {
        return it->second;
    }
}



// Given the initial PackedMarkerGraph, cleanup the bubbles
// due to errors and store the result in a new PackedMarkerGraph.
void BubbleCleaner::cleanup(MarkerGraph& markerGraph)
{
    // The maximum period length that this will cleanup.
    // EXPOSE AFTER CODE STABILIZES?   ************
    const uint64_t maxPeriod = 4;

    const bool debug = false;
    BubbleCleaner& bubbleCleaner = *this;

    Bubbles bubbles;
    findBubbles(bubbles);

    if(debug) {
        cout << "Found " << bubbles.size() << " bubbles." << endl;
    }



    // Cleanup the bubbles we found.
    // For now, this is not recursive.
    // This means it will only deal with one level of "bubble within a bubble"
    // situations.
    while(not bubbles.empty()) {
        const auto it = bubbles.begin();
        vertex_descriptor v0;
        vertex_descriptor v1;
        tie(v0, v1) = it->first;
        const Bubble& bubble = it->second;

        // Sanity check.
        for(const edge_descriptor e: bubble) {
            SHASTA_ASSERT(source(e, bubbleCleaner) == v0);
            SHASTA_ASSERT(target(e, bubbleCleaner) == v1);
        }

        // Get the sequence of the branches.
        vector< vector<Base> > sequences;
        getBubbleSequences(bubble, sequences);

        if(debug) {
            cout << "Working on a bubble with " << bubble.size() << " branches." << endl;
            cout << "MarkerGraph vertices " <<
                bubbleCleaner[v0].markerGraphVertexId << " " <<
                bubbleCleaner[v1].markerGraphVertexId << endl;

            for(uint64_t i=0; i<bubble.size(); i++) {
                const edge_descriptor e = bubble[i];
                copy(sequences[i].begin(), sequences[i].end(), ostream_iterator<Base>(cout));
                cout << " " << bubbleCleaner[e].representation() << endl;
            }
        }

        // See if the branches differ by a repeat count in a short repeat, with period
        // up to maxPeriod.
        const uint64_t period = computeCopyNumberDifferencePeriod(sequences, maxPeriod);
        if(period == 0) {
            bubbles.erase(it);
            continue;
        }
        if(period > maxPeriod) {
            bubbles.erase(it);
            continue;
        }

        if(debug) {
            cout << "This bubble describes copy number changes in a repeat of period " << period << endl;
        }

        // Compute average edge coverage for the branches of this bubble.
        vector<double> coverage(bubbles.size());
        computeBubbleCoverage(bubble, coverage);

        if(debug) {
            for(uint64_t i=0; i<bubble.size(); i++) {
                const edge_descriptor e = bubble[i];
                if(debug) {
                    cout << bubbleCleaner[e].representation() <<
                        " has length " << sequences[i].size() << " and coverage " << coverage[i] << endl;
                }
            }
        }

        // Compute a coverage-weighted average of the lengths.
        double sum = 0.;
        double sumCoverage = 0.;
        for(uint64_t i=0; i<bubble.size(); i++) {
            sum += coverage[i] * double(sequences[i].size());
            sumCoverage += coverage[i];
        }
        const double weightedAverageLength = double(sum) / double(sumCoverage);

        if(debug) {
            cout << "Coverage-weighted average length " << weightedAverageLength << endl;
        }

        // Find the index of the branch that is closest to the weighted average.
        uint64_t iBest = invalid<uint64_t>;
        double bestDelta = 0.;
        for(uint64_t i=0; i<bubble.size(); i++) {
            const double delta = fabs(double(sequences[i].size()) - weightedAverageLength);
            if(i == 0 or delta < bestDelta) {
                iBest = i;
                bestDelta = delta;
            }
        }
        SHASTA_ASSERT(iBest != invalid<uint64_t>);

        if(debug) {
            cout << "Best approximation branch is " << bubbleCleaner[bubble[iBest]].representation() << endl;
        }

        // Flag as removed the marker graph edges of the other branches.
        // Remove the BubbleCleaner edges corresponding to the other branches.
        for(uint64_t i=0; i<bubble.size(); i++) {
            if(i == iBest) {
                continue;
            }
            const edge_descriptor e = bubble[i];
            const BubbleCleanerEdge& edge = bubbleCleaner[e];
            for(const uint64_t segmentId: edge.segments) {
                const auto segment = packedMarkerGraph.segments[segmentId];
                for(const uint64_t edgeId: segment) {
                    markerGraph.edges[edgeId].isSuperBubbleEdge = 1;
                }
            }
            boost::remove_edge(e, bubbleCleaner);
        }

        // Merge the surviving bubble with the previous and/or next edge, if possible.
        edge_descriptor e = bubble[iBest];
        edge_descriptor eNew;
        const bool mergedWithPrevious = mergeWithPreviousIfPossible(e, eNew);
        e = eNew;
        const bool mergedWithNext = mergeWithNextIfPossible(e, eNew);
        e = eNew;



        // If merge occurred, see if we can add the merged edge to an existing
        // unprocessed bubble, or if we can use it to create a new bubble.
        if(mergedWithPrevious or mergedWithNext) {
            const vertex_descriptor v0 = source(e, bubbleCleaner);
            const vertex_descriptor v1 = target(e, bubbleCleaner);
            auto it = bubbles.find(VertexPair(v0, v1));
            if(it != bubbles.end()) {
                // Add it to this bubble.
                it->second.push_back(e);
            } else {
                // We can't add it to an existing unprocessed bubble.
                // See if we can create a new bubble.
                Bubble bubble;
                bool noWay = false;
                BGL_FORALL_OUTEDGES(v0, e, bubbleCleaner, BubbleCleaner) {
                    if(target(e, bubbleCleaner) != v1) {
                        noWay = true;
                        break;
                    }
                    bubble.push_back(e);
                }
                if(not noWay and bubble.size() > 1) {
                    BGL_FORALL_INEDGES(v1, e, bubbleCleaner, BubbleCleaner) {
                        if(source(e, bubbleCleaner) != v0) {
                            noWay = true;
                            break;
                        }
                    }
                    if(not noWay) {
                        bubbles.insert(make_pair(VertexPair(v0, v1), bubble));
                    }
                }
            }
        }

        // Remove this bubble from our list.
        bubbles.erase(it);
    }

}



// Compute average edge coverage for the branches of a bubble.
void BubbleCleaner::computeBubbleCoverage(
    const Bubble& bubble,
    vector<double>& coverage) const
{
    coverage.resize(bubble.size());
    for(uint64_t i=0; i<bubble.size(); i++) {
        const edge_descriptor e = bubble[i];
        coverage[i] = averageEdgeCoverage(e);
    }
}



// Get assembled sequences of the branches of a bubble.
void BubbleCleaner::getBubbleSequences(
    const Bubble& bubble,
    vector< vector<Base> >& sequences) const
{
    const BubbleCleaner& bubbleCleaner = *this;

    sequences.resize(bubble.size());
    for(uint64_t i=0; i<bubble.size(); i++) {
        const edge_descriptor e = bubble[i];
        bubbleCleaner[e].assembledSequence(packedMarkerGraph, sequences[i]);
    }

}



// A bubble is a set of parallel edges between two vertices (v0, v1).
void BubbleCleaner::findBubbles(Bubbles& bubbles) const
{
    const BubbleCleaner& bubbleCleaner = *this;
    bubbles.clear();

    BGL_FORALL_EDGES(e, bubbleCleaner, BubbleCleaner) {
        const vertex_descriptor v0 = source(e, bubbleCleaner);
        const vertex_descriptor v1 = target(e, bubbleCleaner);
        bubbles[VertexPair(v0, v1)].push_back(e);
    }

    // Only keep the ones with at least 2 edges.
    for(auto it=bubbles.begin(); it!=bubbles.end(); /* Increment later */) {
        auto itNext = it;
        ++itNext;
        if(it->second.size() < 2) {
            bubbles.erase(it);
        }
        it = itNext;
    }
}



// Given assembled sequences of the branches of a bubble,
// figure out if this is a bubble caused by copy number
// differences in repeats of period up to maxPeriod.
// If this is the case, returns the shortest period for which this is true.
// Otherwise, returns 0.
uint64_t BubbleCleaner::computeCopyNumberDifferencePeriod(
    const vector< vector<Base> >& sequences,
    uint64_t maxPeriod)
{

    // Check all pairs of branches.
    vector<uint64_t> periods;
    for(uint64_t i=0; i<sequences.size()-1; i++) {
        for(uint64_t j=i+1; j<sequences.size(); j++) {
            const uint64_t pairPeriod =
                shasta::isCopyNumberDifference(sequences[i], sequences[j], maxPeriod);
            if(pairPeriod == 0) {
                return 0;
            }
            periods.push_back(pairPeriod);
        }
    }
    deduplicate(periods);


    if(periods.size() == 1) {
        return periods.front();
    } else {
        return 0;
    }
}



// Merge an edge with its only previous edge, if possible.
// If the merge  was done this returns true and eNew is the
// newly created edge.
// Otherwise, this returns false and eNew is set equal to e.
bool BubbleCleaner::mergeWithPreviousIfPossible(
    edge_descriptor e,
    edge_descriptor& eNew
    )
{
    BubbleCleaner& bubbleCleaner = *this;

    const vertex_descriptor v0 = source(e, bubbleCleaner);
    if(in_degree(v0, bubbleCleaner) == 1 and out_degree(v0, bubbleCleaner) == 1) {

        // We can merge (except for exceptional case below).
        const vertex_descriptor v1 = target(e, bubbleCleaner);

        in_edge_iterator it;
        tie(it, ignore) = in_edges(v0, bubbleCleaner);
        const edge_descriptor ePrevious = *it;
        if(ePrevious == e) {
            // We can't merge.
            eNew = e;
            return false;
        }
        const vertex_descriptor v2 = source(ePrevious, bubbleCleaner);

        // Add the new edge.
        bool edgeWasAdded = false;
        tie(eNew, edgeWasAdded) = add_edge(v2, v1, bubbleCleaner);
        SHASTA_ASSERT(edgeWasAdded);
        BubbleCleanerEdge& newEdge = bubbleCleaner[eNew];

        // Store the segments of the new edge.
        const BubbleCleanerEdge& previousEdge = bubbleCleaner[ePrevious];
        newEdge.segments = previousEdge.segments;
        const BubbleCleanerEdge& edge = bubbleCleaner[e];
        copy(edge.segments.begin(), edge.segments.end(), back_inserter(newEdge.segments));

        // Remove the old edges.
        boost::remove_edge(ePrevious, bubbleCleaner);
        boost::remove_edge(e, bubbleCleaner);

        return true;
    } else {

        // We can't merge.
        eNew = e;
        return false;
    }
}






// Merge an edge with its only next edge, if possible.
// If the merge  was done this returns true and eNew is the
// newly created edge.
// Otherwise, this returns false and eNew is set equal to e.
bool BubbleCleaner::mergeWithNextIfPossible(
    edge_descriptor e,
    edge_descriptor& eNew
    )
{
    BubbleCleaner& bubbleCleaner = *this;

    const vertex_descriptor v1 = target(e, bubbleCleaner);
    if(in_degree(v1, bubbleCleaner) == 1 and out_degree(v1, bubbleCleaner) == 1) {

        // We can merge (except for exceptional case below).
        const vertex_descriptor v0 = source(e, bubbleCleaner);

        out_edge_iterator it;
        tie(it, ignore) = out_edges(v1, bubbleCleaner);
        const edge_descriptor eNext = *it;
        if(eNext == e) {
            // We can't merge.
            eNew = e;
            return false;
        }
        const vertex_descriptor v2 = target(eNext, bubbleCleaner);

        // Add the new edge.
        bool edgeWasAdded = false;
        tie(eNew, edgeWasAdded) = add_edge(v0, v2, bubbleCleaner);
        SHASTA_ASSERT(edgeWasAdded);
        BubbleCleanerEdge& newEdge = bubbleCleaner[eNew];

        // Store the segments of the new edge.
        const BubbleCleanerEdge& edge = bubbleCleaner[e];
        newEdge.segments = edge.segments;
        const BubbleCleanerEdge& nextEdge = bubbleCleaner[eNext];
        copy(nextEdge.segments.begin(), nextEdge.segments.end(), back_inserter(newEdge.segments));

        // Remove the old edges.
        boost::remove_edge(eNext, bubbleCleaner);
        boost::remove_edge(e, bubbleCleaner);

        return true;
    } else {

        // We can't merge.
        eNew = e;
        return false;
    }
}


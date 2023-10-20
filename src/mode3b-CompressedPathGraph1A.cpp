// Shasta.
#include "mode3b-PathGraph1.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "enumeratePaths.hpp"
#include "findLinearChains.hpp"
#include "orderVectors.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3b;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<CompressedPathGraph1A>;



CompressedPathGraph1A::CompressedPathGraph1A(
    const PathGraph1& graph,
    uint64_t componentId,
    const Assembler& assembler) :
    MultithreadedObject<CompressedPathGraph1A>(*this),
    graph(graph),
    componentId(componentId),
    assembler(assembler)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t detangleThresholdLow = 2;
    const uint64_t detangleThresholdHigh = 6;
    const uint64_t pathLengthForChokePoints = 10;
    const uint64_t maxBubbleIndexDelta = 10;

    create();
    writeGfaAndGraphviz("Initial");

    detangle(
        detangleThresholdLow,
        detangleThresholdHigh,
        pathLengthForChokePoints,
        maxBubbleIndexDelta);
    writeGfaAndGraphviz("Final");
}



// Initial creation from the PathGraph1.
void CompressedPathGraph1A::create()
{
    CompressedPathGraph1A& cGraph = *this;

    // Create a filtered version of the PathGraph1, containing only the
    // transitive reduction edges.
    class EdgePredicate {
    public:
        bool operator()(const PathGraph1::edge_descriptor e) const
        {
            return not (*graph)[e].isNonTransitiveReductionEdge;
        }
        EdgePredicate(const PathGraph1& graph) : graph(&graph) {}
        EdgePredicate() : graph(0) {}
    private:
        const PathGraph1* graph;
    };
    using FilteredPathGraph1 = boost::filtered_graph<PathGraph1, EdgePredicate>;
    FilteredPathGraph1 filteredGraph(graph, EdgePredicate(graph));

    // Find linear chains in the PathGraph1 after transitive reduction.
    vector< vector<PathGraph1::edge_descriptor> > chains;
    findLinearChains(filteredGraph, 0, chains);

    // Each chain generates an edge.
    // Vertices are added as needed.
    for(const vector<PathGraph1::edge_descriptor>& chain: chains) {
        const PathGraph1::vertex_descriptor v0 = source(chain.front(), graph);
        const PathGraph1::vertex_descriptor v1 = target(chain.back(), graph);
        const vertex_descriptor cv0 = getCompressedVertex(v0);
        const vertex_descriptor cv1 = getCompressedVertex(v1);

        CompressedPathGraph1AEdge edge;
        edge.id = nextEdgeId++;
        for(const PathGraph1::edge_descriptor e: chain) {
            const PathGraph1::vertex_descriptor v = source(e, graph);
            edge.chain.push_back(graph[v].edgeId);
        }
        const PathGraph1::edge_descriptor eLast = chain.back();
        const PathGraph1::vertex_descriptor vLast = target(eLast, graph);
        edge.chain.push_back(graph[vLast].edgeId);

        add_edge(cv0, cv1, edge, cGraph);
    };

}



// Get the vertex_descriptor corresponding to a PathGraph1::vertex_descriptor,
// adding a vertex if necessary.
CompressedPathGraph1A::vertex_descriptor CompressedPathGraph1A::getCompressedVertex(PathGraph1::vertex_descriptor v)
{
    CompressedPathGraph1A& cGraph = *this;

    auto it = vertexMap.find(v);
    if(it == vertexMap.end()) {
        const vertex_descriptor cv = add_vertex({v}, cGraph);
        vertexMap.insert({v, cv});
        return cv;
    } else {
        return it->second;
    }
}



MarkerGraphEdgeId CompressedPathGraph1A::firstMarkerGraphEdgeId(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;

    SHASTA_ASSERT(chain.size() >= 2);
    return chain.front();
}



MarkerGraphEdgeId CompressedPathGraph1A::lastMarkerGraphEdgeId(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;

    SHASTA_ASSERT(chain.size() >= 2);
    return chain.back();
}



MarkerGraphEdgeId CompressedPathGraph1A::secondMarkerGraphEdgeId(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;

    SHASTA_ASSERT(chain.size() >= 2);
    return chain[1];
}



MarkerGraphEdgeId CompressedPathGraph1A::secondToLastMarkerGraphEdgeId(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const auto& chain = cGraph[ce].chain;

    SHASTA_ASSERT(chain.size() >= 2);
    return chain[chain.size() - 2];
}



void CompressedPathGraph1A::writeGfaAndGraphviz(const string& fileNamePrefix) const
{
    writeGfa(fileNamePrefix);
    writeGraphviz(fileNamePrefix);
}



void CompressedPathGraph1A::writeGfa(const string& fileNamePrefix) const
{
    const CompressedPathGraph1A& cGraph = *this;

    ofstream gfa("CompressedPathGraph1A-" + to_string(componentId) + "-" + fileNamePrefix + ".gfa");

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each edge.
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1A) {

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edgeStringId(ce) << "\t";

        // Sequence.
        gfa << "*\t";

        // Sequence length in bases.
        gfa << "LN:i:" << totalBaseOffset(ce) << "\n";
    }


    // For each vertex, write links between each pair of incoming/outgoing edges.
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1A) {
        BGL_FORALL_INEDGES(cv, ceIn, cGraph, CompressedPathGraph1A) {
            BGL_FORALL_OUTEDGES(cv, ceOut, cGraph, CompressedPathGraph1A) {
                gfa <<
                    "L\t" <<
                    edgeStringId(ceIn) << "\t+\t" <<
                    edgeStringId(ceOut) << "\t+\t*\n";
            }
        }
    }

}



void CompressedPathGraph1A::writeGraphviz(const string& fileNamePrefix) const
{
    const CompressedPathGraph1A& cGraph = *this;

    cout << fileNamePrefix << ": " << num_vertices(cGraph) <<
        " vertices, " << num_edges(cGraph) << " edges." << endl;

    ofstream dot("CompressedPathGraph1A-" + to_string(componentId) + "-" + fileNamePrefix + ".dot");
    dot << "digraph CompressedPathGraph1A_" << componentId << " {\n";


    // Vertices.
    BGL_FORALL_VERTICES(cv, cGraph, CompressedPathGraph1A) {
        const PathGraph1::vertex_descriptor v = cGraph[cv].v;
        dot << graph[v].edgeId << ";\n";
    }



    // Edges.
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1A) {
        const vertex_descriptor cv0 = source(ce, cGraph);
        const vertex_descriptor cv1 = target(ce, cGraph);

        const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v;
        const PathGraph1::vertex_descriptor v1 = cGraph[cv1].v;

        const auto& chain = cGraph[ce].chain;
        SHASTA_ASSERT(chain.size() >= 2);
        dot << graph[v0].edgeId << "->" << graph[v1].edgeId;

        // Label.
        dot << " [label=\"" << edgeStringId(ce) <<
            "\\n" << totalBaseOffset(ce);
        if(chain.size() == 2) {
            // Nothing else
        } else if(chain.size() == 3) {
            dot << "\\n" <<
            secondMarkerGraphEdgeId(ce);
        } else if(chain.size() == 4) {
            dot << "\\n" <<
                secondMarkerGraphEdgeId(ce) << "\\n" <<
                secondToLastMarkerGraphEdgeId(ce);
        } else {
            dot << "\\n" <<
                secondMarkerGraphEdgeId(ce) << "\\n...\\n" <<
                secondToLastMarkerGraphEdgeId(ce);
        }
        dot << "\"];\n";
    }

    dot << "}\n";
}



string CompressedPathGraph1A::edgeStringId(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    return to_string(componentId) + "-" + to_string(cGraph[ce].id);
}



void CompressedPathGraph1A::detangle(
    uint64_t detangleThresholdLow,
    uint64_t detangleThresholdHigh,
    uint64_t pathLengthForChokePoints,
    uint64_t maxBubbleIndexDelta
    )
{
    while(
        detangleEdges(
            detangleThresholdLow,
            detangleThresholdHigh));

    detangleUsingChokePoints(
        pathLengthForChokePoints,
        maxBubbleIndexDelta,
        detangleThresholdLow,
        detangleThresholdHigh);
}



uint64_t CompressedPathGraph1A::detangleEdges(
    uint64_t detangleThresholdLow,
    uint64_t detangleThresholdHigh
    )
{
    const CompressedPathGraph1A& cGraph = *this;

    // The detagling process will create new edges and remove some.
    // Create a set of the edges that exist now.
    // We will erase edges from the set as they get removed,
    // but never add to it.
    std::set<edge_descriptor> initialEdges;
    BGL_FORALL_EDGES(ce, cGraph, CompressedPathGraph1A) {
        initialEdges.insert(ce);
    }
    SHASTA_ASSERT(initialEdges.size() == num_edges(cGraph));

    vector<edge_descriptor> removedEdges;
    uint64_t detangledCount = 0;
    for(auto it=initialEdges.begin(); it!=initialEdges.end(); /* Increment later */) {
        const edge_descriptor ce = *it;
        const bool wasDetangled = detangleEdge(
            ce,
            detangleThresholdLow,
            detangleThresholdHigh,
            removedEdges);
        if(wasDetangled) {
            ++detangledCount;

            // Before removing any edges from the initialEdges set,
            // increment the iterator to point to an edge that will not be removed.
            sort(removedEdges.begin(), removedEdges.end());
            while(it != initialEdges.end()) {
                if(not binary_search(removedEdges.begin(), removedEdges.end(), *it)) {
                    break;
                }
                ++it;
            }

            for(const edge_descriptor ceRemoved: removedEdges) {
                initialEdges.erase(ceRemoved);
            }
        } else {
            ++it;
        }
    }
    cout << detangledCount << " edges were detangled." << endl;
    return detangledCount;
}



bool CompressedPathGraph1A::detangleEdge(
    edge_descriptor ce,
    uint64_t detangleThresholdLow,
    uint64_t detangleThresholdHigh,
    vector<edge_descriptor>& removedEdges)
{
    CompressedPathGraph1A& cGraph = *this;
    removedEdges.clear();

    // The source vertex must have out-degree 1
    // (this is the only edge out of it).
    const vertex_descriptor cv0 = source(ce, cGraph);
    if(out_degree(cv0, cGraph) > 1) {
        return false;
    }

    // The target vertex must have out-degree 1
    // (this is the only edge in from it).
    const vertex_descriptor cv1 = target(ce, cGraph);
    if(in_degree(cv1, cGraph) > 1) {
        return false;
    }

    // Compute the TangleMatrix.
    TangleMatrix tangleMatrix;
    computeTangleMatrix(cv0, cv1, tangleMatrix);

    // Only detangle if in-degree and out-degree are both 2.
    if(tangleMatrix.inDegree() != 2) {
        return false;
    }
    if(tangleMatrix.outDegree() != 2) {
        return false;
    }

#if 0
    cout << "Tangle matrix for " << edgeStringId(ce) << endl;
    for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
        const edge_descriptor ce0 = tangleMatrix.inEdges[i];

        for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
            const edge_descriptor ce1 = tangleMatrix.outEdges[j];

            cout << edgeStringId(ce0) << " ";
            cout << edgeStringId(ce1) << " ";
            cout << tangleMatrix.m[i][j] << endl;
        }
    }
#endif


    // We can only detangle if each column and row of the tangle matrix contains
    // exactly one element >= detangleThresholdHigh
    // and all other elements are <= detangleThresholdLow.
    for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
        uint64_t highCount = 0;
        for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
            const uint64_t m = tangleMatrix.m[i][j];
            if(m >= detangleThresholdHigh) {
                ++highCount;
            } else if(m > detangleThresholdLow) {
                return false;
            }
        }
        if(highCount != 1) {
            return false;
        }
    }

    for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
        uint64_t highCount = 0;
        for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
            const uint64_t m = tangleMatrix.m[i][j];
            if(m >= detangleThresholdHigh) {
                ++highCount;
            } else if(m > detangleThresholdLow) {
                return false;
            }
        }
        if(highCount != 1) {
            return false;
        }
    }
    // cout << "Can detangle." << endl;

    // Create pairs of edges that are going to be merged.
    // This code assumes degree = 2.
    vector< pair<edge_descriptor, edge_descriptor> > mergePairs;
    if(tangleMatrix.m[0][0] >= detangleThresholdHigh) {
        mergePairs.push_back({tangleMatrix.inEdges[0], tangleMatrix.outEdges[0]});
        mergePairs.push_back({tangleMatrix.inEdges[1], tangleMatrix.outEdges[1]});
    } else {
        mergePairs.push_back({tangleMatrix.inEdges[0], tangleMatrix.outEdges[1]});
        mergePairs.push_back({tangleMatrix.inEdges[1], tangleMatrix.outEdges[0]});
    }



    // Create the merged edges.
    for(const auto& p: mergePairs) {
        const edge_descriptor ce0 = p.first;
        const edge_descriptor ce1 = p.second;
        // cout << "Merging edges " << edgeStringId(ce0) << " " << edgeStringId(ce1) << endl;

        const auto& chain0 = cGraph[ce0].chain;
        const auto& chain1 = cGraph[ce1].chain;

        CompressedPathGraph1AEdge newEdge;
        newEdge.id = nextEdgeId++;
        copy(chain0.begin(), chain0.end() - 1, back_inserter(newEdge.chain));
        copy(chain1.begin() + 1, chain1.end(), back_inserter(newEdge.chain));

        add_edge(source(ce0, cGraph), target(ce1, cGraph), newEdge, cGraph);
    }

    // Remove the old edges.
    boost::remove_edge(ce, cGraph);
    removedEdges.push_back(ce);
    for(const edge_descriptor ce: tangleMatrix.inEdges) {
        boost::remove_edge(ce, cGraph);
        removedEdges.push_back(ce);
    }
    for(const edge_descriptor ce: tangleMatrix.outEdges) {
        boost::remove_edge(ce, cGraph);
        removedEdges.push_back(ce);
    }

    // Remove the vertices.
    clear_vertex(cv0, cGraph);
    clear_vertex(cv1, cGraph);
    remove_vertex(cv0, cGraph);
    remove_vertex(cv1, cGraph);

    return true;
}



void CompressedPathGraph1A::writeTangleMatrix(const TangleMatrix& tangleMatrix) const
{
    for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
        const edge_descriptor ce0 = tangleMatrix.inEdges[i];

        for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
            const edge_descriptor ce1 = tangleMatrix.outEdges[j];

            cout << edgeStringId(ce0) << " ";
            cout << edgeStringId(ce1) << " ";
            cout << tangleMatrix.m[i][j] << endl;
        }
    }

}


#if 0
uint64_t CompressedPathGraph1A::detangleBubbleChains(
    uint64_t detangleThresholdLow,
    uint64_t detangleThresholdHigh
    )
{
    // Find the bubble chains.
    vector<BubbleChain> bubbleChains;
    findBubbleChains(bubbleChains);

    return 0;
}



void CompressedPathGraph1A::findBubbleChains(vector<BubbleChain>& bubbleChains) const
{
    const CompressedPathGraph1A& cGraph = *this;

    // Find the bubbles.
    vector<Bubble> bubbles;
    findBubbles(bubbles);

    // Index the bubbles by their source vertex.
    std::map<vertex_descriptor, Bubble*> bubbleMap;
    for(Bubble& bubble: bubbles) {
        SHASTA_ASSERT(not bubbleMap.contains(bubble.source));
        bubbleMap.insert({bubble.source, &bubble});
    }

    // Find the bubble following/preceding each bubble.
    for(Bubble& bubble: bubbles) {

        // See if there is bubble immediately following.
        auto it = bubbleMap.find(bubble.target);
        if(it != bubbleMap.end()) {
            SHASTA_ASSERT(bubble.nextBubble == 0);
            bubble.nextBubble = it->second;
            SHASTA_ASSERT(bubble.nextBubble->previousBubble == 0);
            bubble.nextBubble->previousBubble = &bubble;
            continue;
        }

        // See if there is a bubble following, with an edge in between.
        if(out_degree(bubble.target, cGraph) != 1) {
            continue;
        }
        edge_descriptor nextEdge;
        BGL_FORALL_OUTEDGES(bubble.target, ce, cGraph, CompressedPathGraph1A) {
            nextEdge = ce;
        }
        const vertex_descriptor cv1 = target(nextEdge, cGraph);
        if(in_degree(cv1, cGraph) != 1) {
            continue;
        }
        it = bubbleMap.find(cv1);
        if(it != bubbleMap.end()) {
            SHASTA_ASSERT(bubble.nextBubble == 0);
            bubble.nextBubble = it->second;
            SHASTA_ASSERT(bubble.nextBubble->previousBubble == 0);
            bubble.nextBubble->previousBubble = &bubble;
        }
    }

    // Sanity check.
    uint64_t noNextCount = 0;
    uint64_t noPreviousCount = 0;
    for(const Bubble& bubble: bubbles) {
        if(bubble.nextBubble == 0) {
            noNextCount++;
        }
        if(bubble.previousBubble == 0) {
            noPreviousCount++;
        }
    }
    SHASTA_ASSERT(noNextCount == noPreviousCount);


    // Now we can construct the BubbleChains.
    bubbleChains.clear();
    for(const Bubble& bubble: bubbles) {

        if(bubble.previousBubble) {
            // A BubbleChain does not begin at this bubble.
            continue;
        }

        bubbleChains.resize(bubbleChains.size() + 1);
        BubbleChain& bubbleChain = bubbleChains.back();
        const Bubble* b = &bubble;
        while(b) {
            bubbleChain.bubbles.push_back(b);
            b = b->nextBubble;
        }
    }
    cout << "Found " << bubbleChains.size() << " bubble chains." << endl;

    // Find the InterBubbles.
    for(BubbleChain& bubbleChain: bubbleChains) {
        for(uint64_t i=1; i<bubbleChain.bubbles.size(); i++) {
            const Bubble* bubble0 = bubbleChain.bubbles[i-1];
            const Bubble* bubble1 = bubbleChain.bubbles[i];

            if(bubble0->target == bubble1->source) {
                // There is no edge in between.
                bubbleChain.interBubbles.push_back({edge_descriptor(), false});
            } else {
                // Find the edge in between.
                edge_descriptor ce;
                bool edgeWasFound;
                tie(ce, edgeWasFound) = boost::edge(bubble0->target, bubble1->source, cGraph);
                SHASTA_ASSERT(edgeWasFound);
                SHASTA_ASSERT(source(ce, cGraph) == bubble0->target);
                SHASTA_ASSERT(target(ce, cGraph) == bubble1->source);
                SHASTA_ASSERT(out_degree(bubble0->target, cGraph) == 1);
                SHASTA_ASSERT(in_degree(bubble1->source, cGraph) == 1);
                bubbleChain.interBubbles.push_back({ce, true});
            }
        }
    }

    if(true) {
        for(const BubbleChain& bubbleChain: bubbleChains) {
            cout << "Bubble chain with " << bubbleChain.bubbles.size() << " bubbles:" << flush;
            for(uint64_t i=0; /* Check later */; i++) {
                const Bubble* bubble = bubbleChain.bubbles[i];
                SHASTA_ASSERT(bubble);
                cout << " ";
                writeBubble(*bubble, cout);
                cout << flush;

                if(i == bubbleChain.bubbles.size() - 1) {
                    break;
                }

                const InterBubble& interbubble = bubbleChain.interBubbles[i];
                if(interbubble.edgeExists) {
                    cout << " " << edgeStringId(interbubble.edge);
                }
            }
            cout << endl;


            // Write tangle matrices between pairs of bubbles in this bubble chain.
            for(uint64_t i0=0; i0<bubbleChain.bubbles.size()-1; i0++) {
                const Bubble& bubble0 = *bubbleChain.bubbles[i0];
                const vertex_descriptor cv0 = bubble0.target;
                for(uint64_t i1=i0+1; i1<bubbleChain.bubbles.size(); i1++) {
                    const Bubble& bubble1 = *bubbleChain.bubbles[i1];
                    const vertex_descriptor cv1 = bubble1.source;
                    TangleMatrix tangleMatrix;
                    computeTangleMatrix(cv0, cv1, tangleMatrix);

                    cout << "Tangle matrix between bubbles " << i0 << " ";
                    writeBubble(bubble0, cout);
                    cout << " " << i1 << " ";
                    writeBubble(bubble1, cout);
                    cout << endl;

                    for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
                        const edge_descriptor ce0 = tangleMatrix.inEdges[i];

                        for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
                            const edge_descriptor ce1 = tangleMatrix.outEdges[j];

                            cout << edgeStringId(ce0) << " ";
                            cout << edgeStringId(ce1) << " ";
                            cout << tangleMatrix.m[i][j] << endl;
                        }
                    }
                }
            }
        }
    }

}




void CompressedPathGraph1A::findBubbles(vector<Bubble>& bubbles) const
{
    const CompressedPathGraph1A& cGraph = *this;
    bubbles.clear();

    vector<vertex_descriptor> targets;
    BGL_FORALL_VERTICES(cv0, cGraph, CompressedPathGraph1A) {

        // For cv0 to be the source of a bubble, all its out-edges must have the same target.
        targets.clear();
        BGL_FORALL_OUTEDGES(cv0, ce, cGraph, CompressedPathGraph1A) {
            targets.push_back(target(ce, cGraph));
        }
        if(targets.size() < 2) {
            continue;
        }
        deduplicate(targets);
        if(targets.size() != 1) {
            // No bubble with cv0 as the source.
            continue;
        }
        const vertex_descriptor cv1 = targets.front();

        // We must also check that all in-edges of cv1 have cv0 as their source.
        bool isBubble = true;
        BGL_FORALL_INEDGES(cv1, ce, cGraph, CompressedPathGraph1A) {
            if(source(ce, cGraph) != cv0) {
                isBubble = false;
                break;
            }
        }
        if(not isBubble) {
            continue;
        }

        // Create a bubble.
        Bubble bubble;
        bubble.source = cv0;
        bubble.target = cv1;
        BGL_FORALL_OUTEDGES(cv0, ce, cGraph, CompressedPathGraph1A) {
            bubble.edges.push_back(ce);
        }
        bubbles.push_back(bubble);

    }

    cout << "Found " << bubbles.size() << " bubbles." << endl;

    // Bubble size histogram.
    {
        vector<uint64_t> histogram;
        for(const Bubble& bubble: bubbles) {
            const uint64_t degree = bubble.degree();
            if(degree >= histogram.size()) {
                histogram.resize(degree+1, 0);
            }
            ++histogram[degree];
        }
        cout << "Bubble degree histogram" << endl;
        for(uint64_t degree=0; degree<histogram.size(); degree++) {
            const uint64_t frequency = histogram[degree];
            if(frequency) {
                cout << degree << " " << frequency << endl;
            }
        }
    }

    // Bubble output.
    if(false) {
        for(const Bubble& bubble: bubbles) {
            const PathGraph1::vertex_descriptor v0 = cGraph[bubble.source].v;
            const PathGraph1::vertex_descriptor v1 = cGraph[bubble.target].v;

            cout << "Bubble between vertices " <<
                graph[v0].edgeId << " " <<
                graph[v1].edgeId << ":";
            for(const edge_descriptor ce: bubble.edges) {
                cout << " " << edgeStringId(ce);
            }
            cout << endl;
        }

    }

}
#endif


#if 0
void CompressedPathGraph1A::analyzeChokePoints() const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t pathLength = 10;

    const CompressedPathGraph1A& cGraph = *this;

    using Path = vector<edge_descriptor>;
    class PathInspector {
    public:
        PathInspector(
            const CompressedPathGraph1A& cGraph,
            uint64_t direction  // 0=forward, 1=backward
            ) :
            cGraph(cGraph), direction(direction) {}
        void operator()(const Path& path)
        {
            if(path.size() == pathLength) {
                ++pathCount;
                // cout << "Path: ";
                for(const edge_descriptor ce: path) {
                    const vertex_descriptor cv = (direction == 0) ? target(ce, cGraph) : source(ce, cGraph);
                    ++vertexCountMap[cv];
                    // const PathGraph1::vertex_descriptor v = cGraph[cv].v;
                    // cout << " " << cGraph.graph[v].edgeId;
                }
                // cout << endl;
            }
        }

        const CompressedPathGraph1A& cGraph;
        uint64_t direction;

        // The number of paths with the specified length.
        uint64_t pathCount = 0;

        // The number of paths each vertex appears in.
        std::map<vertex_descriptor, uint64_t> vertexCountMap;
    };

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(cv0, cGraph, CompressedPathGraph1A) {

        // Forward
        {
            PathInspector pathInspector(cGraph, 0);
            enumeratePaths(cGraph, cv0, pathLength, pathInspector);

            for(const auto& p: pathInspector.vertexCountMap) {
                if(p.second == pathInspector.pathCount) {
                    // This vertex appears in all paths.
                    const vertex_descriptor cv1 = p.first;
                    forwardPairs.push_back({cv0, cv1});
                }
            }
        }

        // Backward
        {
            PathInspector pathInspector(cGraph, 1);
            enumeratePathsReverse(cGraph, cv0, pathLength, pathInspector);

            for(const auto& p: pathInspector.vertexCountMap) {
                if(p.second == pathInspector.pathCount) {
                    // This vertex appears in all paths.
                    const vertex_descriptor cv1 = p.first;
                    backwardPairs.push_back({cv1, cv0});
                }
            }
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());
    sort(backwardPairs.begin(), backwardPairs.end());

    // Find the pairs that appear in both directions.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));
    cout << "Found " << forwardPairs.size() << " forward pairs." << endl;
    cout << "Found " << backwardPairs.size() << " backward pairs." << endl;
    cout << "Found " << bidirectionalPairs.size() << " bidirectional pairs." << endl;



    // The pairs that appear in both directions define a directed graph.
    // Its vertices are the "choke points" or "bottlenecks" of the CompressedPathGraph1A.

    ofstream dot("BidirectionalPairs.dot");
    dot << "digraph BidirectionalPairs {\n";
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor cv0 = p.first;
        const vertex_descriptor cv1 = p.second;

        const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v;
        const PathGraph1::vertex_descriptor v1 = cGraph[cv1].v;

        dot << graph[v0].edgeId << "->" << graph[v1].edgeId << ";\n";
    }
    dot << "}\n";

    // Gather the choke points.
    vector<vertex_descriptor> chokePoints;
    for(const auto& p: bidirectionalPairs) {
        chokePoints.push_back(p.first);
        chokePoints.push_back(p.second);
    }
    deduplicate(chokePoints);
    cout << "Found " << chokePoints.size() << " choke points." << endl;

    // Create the choke point graph.
    // Each vertex corresponds to a vertex of the CompressedPathGraph1A.
    using ChokePointGraph = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            vertex_descriptor>;
    ChokePointGraph chokePointGraph;
    for(const vertex_descriptor cv: chokePoints) {
        add_vertex({cv}, chokePointGraph);
    }
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor cv0 = p.first;
        const vertex_descriptor cv1 = p.second;
        const uint64_t iv0 = lower_bound(chokePoints.begin(), chokePoints.end(), cv0) - chokePoints.begin();
        const uint64_t iv1 = lower_bound(chokePoints.begin(), chokePoints.end(), cv1) - chokePoints.begin();
        add_edge(iv0, iv1, chokePointGraph);
    }



    // If the ChokePointGraph happens to have any non-trivial strongly connected components,
    // remove them by disconnecting their vertices. We don't want to remove the vertices
    // because we are using vecS for the ChokePoint to simplify operations.
    std::map<ChokePointGraph::vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        chokePointGraph,
        boost::make_assoc_property_map(componentMap));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<ChokePointGraph::vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    // Disconnect vertices belonging to non-trivial strongly connected components.
    // A non-trivial strong component has at least one internal edge.
    // This means that it either has more than one vertex,
    // or it consists of a single vertex with a self-edge.
    uint64_t stronglyConnectedChokePointsCount = 0;;
    for(const auto& p: componentVertices) {

        // Figure out if it is non-trivial.
        bool isNonTrivial;
        if(p.second.size() > 1) {

            // More than one vertex. Certainly non-trivial.
            isNonTrivial = true;
        } else if (p.second.size() == 1) {

            // Only one vertex. Non-trivial if self-edge present.
            const ChokePointGraph::vertex_descriptor v = p.second.front();
            bool selfEdgeExists = false;
            tie(ignore, selfEdgeExists) = edge(v, v, chokePointGraph);
            isNonTrivial = selfEdgeExists;
        } else {

            // Empty. This should never happen.
            SHASTA_ASSERT(0);
        }

        // If non-trivial, disconnect all of its vertices.
        if(isNonTrivial) {
            for(const ChokePointGraph::vertex_descriptor v: p.second) {
                clear_vertex(v, chokePointGraph);
                ++stronglyConnectedChokePointsCount;
            }
        }
    }
    cout << "Found " << stronglyConnectedChokePointsCount <<
        " choke points in strongly connected components and disconnected them." << endl;

    // Now we can compute the transitive reduction of the ChokePointGraph.
    transitiveReduction(chokePointGraph);
    cout << "After transitive reduction, the choke point graph has " <<
        num_vertices(chokePointGraph) << " vertices and " <<
        num_edges(chokePointGraph) << " edges." << endl;
    {
        ofstream dot("ChokePointGraph.dot");
        dot << "digraph ChokePointGraph{\n";
        BGL_FORALL_EDGES(e, chokePointGraph, ChokePointGraph) {
            const ChokePointGraph::vertex_descriptor cpv0 = source(e, chokePointGraph);
            const ChokePointGraph::vertex_descriptor cpv1 = target(e, chokePointGraph);
            const vertex_descriptor cv0 = chokePoints[cpv0];
            const vertex_descriptor cv1 = chokePoints[cpv1];
            const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v;
            const PathGraph1::vertex_descriptor v1 = cGraph[cv1].v;
            dot << graph[v0].edgeId << "->" << graph[v1].edgeId << ";\n";
        }
        dot << "}\n";
    }

    // Linear chains in the ChokePointGraph will be used
    // to create bubble chains for detangling.
    vector< vector<ChokePointGraph::edge_descriptor> > chains;
    findLinearChains(chokePointGraph, 0, chains);
    cout << "Found " << chains.size() << " linear chains in the choke point graph." << endl;

    // Compute a histogram of chain lengths.
    {
        vector<uint64_t> histogram;
        for(const auto& chain: chains) {
            const uint64_t length = chain.size();
            if(histogram.size() <= length) {
                histogram.resize(length + 1, 0);
            }
            ++histogram[length];
        }
        ofstream csv("ChokePointChainLengthHistogram.csv");
        csv << "Length,Frequency\n";
        for(uint64_t length=0; length<histogram.size(); length++) {
            const uint64_t frequency = histogram[length];
            if(frequency) {
                csv << length << "," << frequency << "\n";
            }
        }
    }

}
#endif



void CompressedPathGraph1A::detangleUsingChokePoints(
    uint64_t pathLengthForChokePoints,
    uint64_t maxBubbleIndexDelta,
    uint64_t detangleThresholdLow,
    uint64_t detangleThresholdHigh)
{
    vector<ChokePointChain> chokePointChains;
    findChokePointChains(pathLengthForChokePoints, chokePointChains);

    for(uint64_t chokePointChainId=0; chokePointChainId<chokePointChains.size(); chokePointChainId++) {
        ChokePointChain& chokePointChain = chokePointChains[chokePointChainId];
        detangleChokePointChain(
            chokePointChain,
            chokePointChainId,
            maxBubbleIndexDelta,
            detangleThresholdLow,
            detangleThresholdHigh);
    }

}



void CompressedPathGraph1A::detangleChokePointChain(
    ChokePointChain& chain,
    uint64_t chokePointChainId,
    uint64_t maxBubbleIndexDelta,
    uint64_t detangleThresholdLow,
    uint64_t detangleThresholdHigh)
{
    // If the chain has less that 2 diploid bubbles, don't do anything.
    if(chain.diploidBubblesIndexes.size() < 2) {
        return;
    }

    cout << "Detangling a choke point chain with " << chain.chokePoints.size() << " choke points." << endl;
    writeChokePointChain(chain);

    // Use tangle matrices for near pairs of bubbles to create a PhasingGraph.
    PhasingGraph phasingGraph(chain.diploidBubblesIndexes.size());
    for(uint64_t i0=0; i0<chain.diploidBubblesIndexes.size()-1; i0++) {
        const uint64_t j0 = chain.diploidBubblesIndexes[i0];
        const Superbubble& diploidBubble0 = chain.superbubbles[j0];
        SHASTA_ASSERT(diploidBubble0.isDiploidBubble);
        SHASTA_ASSERT(diploidBubble0.diploidEdges.size() == 2);
        for(uint64_t i1=i0+1; i1<min(i0+maxBubbleIndexDelta+1, chain.diploidBubblesIndexes.size()); i1++) {
            const uint64_t j1 = chain.diploidBubblesIndexes[i1];
            const Superbubble& diploidBubble1 = chain.superbubbles[j1];
            SHASTA_ASSERT(diploidBubble1.isDiploidBubble);
            SHASTA_ASSERT(diploidBubble1.diploidEdges.size() == 2);

            // Compute the tangle matrix between these two vertices.
            TangleMatrix tangleMatrix;
            computeTangleMatrix(diploidBubble0.diploidEdges, diploidBubble1.diploidEdges, tangleMatrix);

            cout << "Tangle matrix for bubbles " << i0 << " " << i1 << endl;
            for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
                const edge_descriptor ce0 = tangleMatrix.inEdges[i];

                for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
                    const edge_descriptor ce1 = tangleMatrix.outEdges[j];

                    cout << edgeStringId(ce0) << " ";
                    cout << edgeStringId(ce1) << " ";
                    cout << tangleMatrix.m[i][j] << endl;
                }
            }

            // Analyze the tangle matrix.
            int64_t phase;
            uint64_t minConcordant;
            uint64_t maxDiscordant;
            tangleMatrix.analyze(detangleThresholdLow, detangleThresholdHigh,
                phase, minConcordant, maxDiscordant);
            if(phase == 0) {
                cout << "Ambiguous";
            } else {
                if (phase == 1) {
                    cout << "In phase";
                } else {
                    cout << "Out of phase";
                }
                cout << " minConcordant " << minConcordant;
                cout << " maxDiscordant " << maxDiscordant;
            }
            cout << endl;

            // If not ambiguous, add an edge to the PhasingGraph.
            if(phase != 0) {
                phasingGraph.addEdge(i0, i1, {phase, minConcordant, maxDiscordant});
            }
        }
    }

    phasingGraph.writeGraphviz(cout, "PhasingGraph_" + to_string(chokePointChainId));
}



CompressedPathGraph1A::PhasingGraph::PhasingGraph(uint64_t bubbleCount)
{
    PhasingGraph& phasingGraph = *this;

    for(uint64_t bubbleId=0; bubbleId<bubbleCount; bubbleId++) {
        vertexMap.insert(make_pair(bubbleId, add_vertex({bubbleId}, phasingGraph)));
    }
}


void CompressedPathGraph1A::PhasingGraph::writeGraphviz(ostream& s, const string& name) const
{
    const PhasingGraph& phasingGraph = *this;

    s << "graph " << name << " {\n";

    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        s << phasingGraph[v].bubbleId << ";\n";
    }

    BGL_FORALL_EDGES(e, phasingGraph, PhasingGraph) {
        const vertex_descriptor v0 = source(e, phasingGraph);
        const vertex_descriptor v1 = target(e, phasingGraph);
        s << phasingGraph[v0].bubbleId << "--";
        s << phasingGraph[v1].bubbleId << ";\n";
    }

    s << "}\n";
}


void CompressedPathGraph1A::PhasingGraph::addEdge(
    uint64_t bubbleId0,
    uint64_t bubbleId1,
    const PhasingGraphEdge& edge)
{
    PhasingGraph& phasingGraph = *this;

    // Get the vertices corresponding to these bubbles.
    auto it0 = vertexMap.find(bubbleId0);
    auto it1 = vertexMap.find(bubbleId1);
    SHASTA_ASSERT(it0 != vertexMap.end());
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const vertex_descriptor v1 = it1->second;

    // Add the edge.
    add_edge(v0, v1, edge, phasingGraph);
}



void CompressedPathGraph1A::findChokePointChains(
    uint64_t pathLengthForChokePoints,
    vector<ChokePointChain>& chokePointChains) const
{
    const CompressedPathGraph1A& cGraph = *this;

    // To find choke points and their connectivity, enumerate forward
    // and backward paths of the specified length and starting at every vertex.
    // If v1 appears in all forward paths starting at v0
    // and v0 appears in all backward aths starting at v1,
    // v0 and v1 are choke points and v0->v1 generates an edge
    // in the ChokePointGraph.


    using Path = vector<edge_descriptor>;
    class PathInspector {
    public:
        PathInspector(
            const CompressedPathGraph1A& cGraph,
            uint64_t direction,  // 0=forward, 1=backward
            uint64_t pathLengthForChokePoints
            ) :
            cGraph(cGraph),
            direction(direction),
            pathLengthForChokePoints(pathLengthForChokePoints) {}
        void operator()(const Path& path)
        {
            if(path.size() == pathLengthForChokePoints) {
                ++pathCount;
                // cout << "Path: ";
                for(const edge_descriptor ce: path) {
                    const vertex_descriptor cv = (direction == 0) ? target(ce, cGraph) : source(ce, cGraph);
                    ++vertexCountMap[cv];
                    // const PathGraph1::vertex_descriptor v = cGraph[cv].v;
                    // cout << " " << cGraph.graph[v].edgeId;
                }
                // cout << endl;
            }
        }

        const CompressedPathGraph1A& cGraph;
        uint64_t direction;
        uint64_t pathLengthForChokePoints;

        // The number of paths with the specified length.
        uint64_t pathCount = 0;

        // The number of paths each vertex appears in.
        std::map<vertex_descriptor, uint64_t> vertexCountMap;
    };

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(cv0, cGraph, CompressedPathGraph1A) {

        // Forward
        {
            PathInspector pathInspector(cGraph, 0, pathLengthForChokePoints);
            enumeratePaths(cGraph, cv0, pathLengthForChokePoints, pathInspector);

            for(const auto& p: pathInspector.vertexCountMap) {
                if(p.second == pathInspector.pathCount) {
                    // This vertex appears in all paths.
                    const vertex_descriptor cv1 = p.first;
                    forwardPairs.push_back({cv0, cv1});
                }
            }
        }

        // Backward
        {
            PathInspector pathInspector(cGraph, 1, pathLengthForChokePoints);
            enumeratePathsReverse(cGraph, cv0, pathLengthForChokePoints, pathInspector);

            for(const auto& p: pathInspector.vertexCountMap) {
                if(p.second == pathInspector.pathCount) {
                    // This vertex appears in all paths.
                    const vertex_descriptor cv1 = p.first;
                    backwardPairs.push_back({cv1, cv0});
                }
            }
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());
    sort(backwardPairs.begin(), backwardPairs.end());

    // Find the pairs that appear in both directions.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));


    // The pairs that appear in both directions define a directed graph.
    // Its vertices are the "choke points" or "bottlenecks" of the CompressedPathGraph1A.

    if(false) {
        ofstream dot("BidirectionalPairs.dot");
        dot << "digraph BidirectionalPairs {\n";
        for(const auto& p: bidirectionalPairs) {
            const vertex_descriptor cv0 = p.first;
            const vertex_descriptor cv1 = p.second;

            const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v;
            const PathGraph1::vertex_descriptor v1 = cGraph[cv1].v;

            dot << graph[v0].edgeId << "->" << graph[v1].edgeId << ";\n";
        }
        dot << "}\n";
    }

    // Gather the choke points.
    vector<vertex_descriptor> chokePoints;
    for(const auto& p: bidirectionalPairs) {
        chokePoints.push_back(p.first);
        chokePoints.push_back(p.second);
    }
    deduplicate(chokePoints);
    cout << "Found " << chokePoints.size() << " choke points." << endl;

    // Create the choke point graph.
    // Each vertex corresponds to a vertex of the CompressedPathGraph1A.
    using ChokePointGraph = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            vertex_descriptor>;
    ChokePointGraph chokePointGraph;
    for(const vertex_descriptor cv: chokePoints) {
        add_vertex({cv}, chokePointGraph);
    }
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor cv0 = p.first;
        const vertex_descriptor cv1 = p.second;
        const uint64_t iv0 = lower_bound(chokePoints.begin(), chokePoints.end(), cv0) - chokePoints.begin();
        const uint64_t iv1 = lower_bound(chokePoints.begin(), chokePoints.end(), cv1) - chokePoints.begin();
        add_edge(iv0, iv1, chokePointGraph);
    }



    // If the ChokePointGraph happens to have any non-trivial strongly connected components,
    // remove them by disconnecting their vertices. We don't want to remove the vertices
    // because we are using vecS for the ChokePoint to simplify operations.
    std::map<ChokePointGraph::vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        chokePointGraph,
        boost::make_assoc_property_map(componentMap));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<ChokePointGraph::vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    // Disconnect vertices belonging to non-trivial strongly connected components.
    // A non-trivial strong component has at least one internal edge.
    // This means that it either has more than one vertex,
    // or it consists of a single vertex with a self-edge.
    uint64_t stronglyConnectedChokePointsCount = 0;;
    for(const auto& p: componentVertices) {

        // Figure out if it is non-trivial.
        bool isNonTrivial;
        if(p.second.size() > 1) {

            // More than one vertex. Certainly non-trivial.
            isNonTrivial = true;
        } else if (p.second.size() == 1) {

            // Only one vertex. Non-trivial if self-edge present.
            const ChokePointGraph::vertex_descriptor v = p.second.front();
            bool selfEdgeExists = false;
            tie(ignore, selfEdgeExists) = edge(v, v, chokePointGraph);
            isNonTrivial = selfEdgeExists;
        } else {

            // Empty. This should never happen.
            SHASTA_ASSERT(0);
        }

        // If non-trivial, disconnect all of its vertices.
        if(isNonTrivial) {
            for(const ChokePointGraph::vertex_descriptor v: p.second) {
                clear_vertex(v, chokePointGraph);
                ++stronglyConnectedChokePointsCount;
            }
        }
    }
    // cout << "Found " << stronglyConnectedChokePointsCount <<
    //     " choke points in strongly connected components and disconnected them." << endl;

    // Now we can compute the transitive reduction of the ChokePointGraph.
    transitiveReduction(chokePointGraph);
    cout << "After transitive reduction, the choke point graph has " <<
        num_vertices(chokePointGraph) << " vertices and " <<
        num_edges(chokePointGraph) << " edges." << endl;
    {
        ofstream dot("ChokePointGraph.dot");
        dot << "digraph ChokePointGraph{\n";
        BGL_FORALL_EDGES(e, chokePointGraph, ChokePointGraph) {
            const ChokePointGraph::vertex_descriptor cpv0 = source(e, chokePointGraph);
            const ChokePointGraph::vertex_descriptor cpv1 = target(e, chokePointGraph);
            const vertex_descriptor cv0 = chokePoints[cpv0];
            const vertex_descriptor cv1 = chokePoints[cpv1];
            const PathGraph1::vertex_descriptor v0 = cGraph[cv0].v;
            const PathGraph1::vertex_descriptor v1 = cGraph[cv1].v;
            dot << graph[v0].edgeId << "->" << graph[v1].edgeId << ";\n";
        }
        dot << "}\n";
    }

    // Linear chains in the ChokePointGraph will be used
    // to create ChokePointChains for detangling.
    vector< vector<ChokePointGraph::edge_descriptor> > chains;
    findLinearChains(chokePointGraph, 0, chains);
    sort(chains.begin(), chains.end(), OrderVectorsByDecreasingSize<ChokePointGraph::edge_descriptor>());
    cout << "Found " << chains.size() << " choke point chains." << endl;

    // Compute a histogram of chain lengths.
    if(false) {
        vector<uint64_t> histogram;
        for(const auto& chain: chains) {
            const uint64_t length = chain.size();
            if(histogram.size() <= length) {
                histogram.resize(length + 1, 0);
            }
            ++histogram[length];
        }
        ofstream csv("ChokePointChainLengthHistogram.csv");
        csv << "Length,Frequency\n";
        for(uint64_t length=0; length<histogram.size(); length++) {
            const uint64_t frequency = histogram[length];
            if(frequency) {
                csv << length << "," << frequency << "\n";
            }
        }
    }



    // Now we can create the chokePointChains.
    chokePointChains.clear();
    for(const vector<ChokePointGraph::edge_descriptor>& chain: chains) {
        chokePointChains.emplace_back();
        ChokePointChain& chokePointChain = chokePointChains .back();

        // Fill in the chokepoints.
        SHASTA_ASSERT(not chain.empty());
        chokePointChain.chokePoints.push_back(chokePoints[source(chain.front(), cGraph)]);
        for(const ChokePointGraph::edge_descriptor cv: chain) {
            chokePointChain.chokePoints.push_back(chokePoints[target(cv, cGraph)]);
        }
        SHASTA_ASSERT(chokePointChain.chokePoints.size() == chain.size() + 1);

        // Fill in the superbubbles.
        for(uint64_t i=1; i<chokePointChain.chokePoints.size(); i++) {
            const vertex_descriptor cv0 = chokePointChain.chokePoints[i-1];
            const vertex_descriptor cv1 = chokePointChain.chokePoints[i];
            chokePointChain.superbubbles.emplace_back();
            Superbubble& superbubble = chokePointChain.superbubbles.back();
            findVerticesAndEdgesBetweenChokePoints(cv0, cv1,
                superbubble.internalVertices, superbubble.internalEdges);

            // Figure out if this superbubble is a diploid bubble.
            if(out_degree(cv0, cGraph) != 2) {
                continue;
            }
            if(in_degree(cv1, cGraph) != 2) {
                continue;
            }
            if(not superbubble.internalVertices.empty()) {
                continue;
            }
            // If getting here, this must be a diploid bubble.
            superbubble.isDiploidBubble = true;
            chokePointChain.diploidBubblesIndexes.push_back(chokePointChain.superbubbles.size() - 1);
            BGL_FORALL_OUTEDGES(cv0, e, cGraph, CompressedPathGraph1A) {
                superbubble.diploidEdges.push_back(e);
            }
        }
    }

}



void CompressedPathGraph1A::findVerticesAndEdgesBetweenChokePoints(
    vertex_descriptor cv0,
    vertex_descriptor cv1,
    vector<vertex_descriptor>& interveningVertices,
    vector<edge_descriptor>& interveningEdges) const
{
    const CompressedPathGraph1A& cGraph = *this;

    // Do a BFS starting at cv0, stopping it at cv1.
    std::queue<vertex_descriptor> q;
    q.push(cv0);
    std::set<vertex_descriptor> visitedVertices;

    while(not q.empty()) {
        const vertex_descriptor cvA = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(cvA, e, cGraph, CompressedPathGraph1A) {
            const vertex_descriptor cvB = target(e, cGraph);
            if((cvB != cv0) and (cvB != cv1) and (not visitedVertices.contains(cvB))) {
                q.push(cvB);
                visitedVertices.insert(cvB);
            }
        }
    }
    visitedVertices.erase(cv0);
    interveningVertices.clear();
    copy(visitedVertices.begin(), visitedVertices.end(), back_inserter(interveningVertices));


    // The intervening edges are the out-edges of cv0 plus the out-edges of the internal vertices.
    interveningEdges.clear();
    BGL_FORALL_OUTEDGES(cv0, ce, cGraph, CompressedPathGraph1A) {
        interveningEdges.push_back(ce);
    }
    for(const vertex_descriptor cv: interveningVertices) {
        BGL_FORALL_OUTEDGES(cv, ce, cGraph, CompressedPathGraph1A) {
            interveningEdges.push_back(ce);
        }
    }
}



void CompressedPathGraph1A::writeChokePointChain(const ChokePointChain& chain) const
{
    const CompressedPathGraph1A& cGraph = *this;

    for(uint64_t i=0; /* Check later */; i++) {

        const vertex_descriptor cv = chain.chokePoints[i];
        const PathGraph1::vertex_descriptor v = cGraph[cv].v;
        cout << "Choke point at position " << i << ": " << graph[v].edgeId << "\n";

        if(i == chain.chokePoints.size() - 1) {
            break;
        }

        const Superbubble& superbubble = chain.superbubbles[i];
        cout << "Superbubble at position " << i << ":";
        if(superbubble.internalVertices.empty()) {
            cout << " no vertices";
        } else {
            cout << " vertices ";
            for(const vertex_descriptor cv: superbubble.internalVertices) {
                const PathGraph1::vertex_descriptor v = cGraph[cv].v;
                cout << " " << graph[v].edgeId;
            }
        }
        if(superbubble.internalEdges.empty()) {
            cout << ", no edges";
        } else {
            cout << ", edges";
            for(const edge_descriptor ce: superbubble.internalEdges) {
                cout << " " << edgeStringId(ce);
            }
        }
        if(superbubble.isDiploidBubble) {
            SHASTA_ASSERT(superbubble.diploidEdges.size() == 2);
            cout << ", diploid bubble, edges ";
            cout << edgeStringId(superbubble.diploidEdges[0]) << " " <<
                edgeStringId(superbubble.diploidEdges[1]);
        }
        cout << "\n";
    }

    // Also list the diploid bubbles again.
    cout << "Diploid bubbles in this choke point chain:\n";
    for(uint64_t i=0; i<chain.diploidBubblesIndexes.size(); i++) {
        const uint64_t j = chain.diploidBubblesIndexes[i];
        const Superbubble& superbubble = chain.superbubbles[j];
        SHASTA_ASSERT(superbubble.isDiploidBubble);
        SHASTA_ASSERT(superbubble.diploidEdges.size() == 2);
        cout << "Diploid bubble " << i << " is superbubble at position " << j <<
            " with edges " <<
            edgeStringId(superbubble.diploidEdges[0]) << " " <<
            edgeStringId(superbubble.diploidEdges[1]) << "\n";
    }


    cout << flush;
}



void CompressedPathGraph1A::computeTangleMatrix(
    vertex_descriptor cv0,
    vertex_descriptor cv1,
    TangleMatrix& tangleMatrix
    ) const
{
    const CompressedPathGraph1A& cGraph = *this;

    tangleMatrix.inEdges.clear();
    BGL_FORALL_INEDGES(cv0, ce, cGraph, CompressedPathGraph1A) {
        tangleMatrix.inEdges.push_back(ce);
    }

    tangleMatrix.outEdges.clear();
    BGL_FORALL_OUTEDGES(cv1, ce, cGraph, CompressedPathGraph1A) {
        tangleMatrix.outEdges.push_back(ce);
    }

    computeTangleMatrix(tangleMatrix);
}



void CompressedPathGraph1A::computeTangleMatrix(
    const vector<edge_descriptor>& inEdges,
    const vector<edge_descriptor>& outEdges,
    TangleMatrix& tangleMatrix
    ) const
{
    tangleMatrix.inEdges = inEdges;
    tangleMatrix.outEdges = outEdges;
    computeTangleMatrix(tangleMatrix);
}



// This version only fills in m. The inEdges and out_edgesmust have already been filled in.
void CompressedPathGraph1A::computeTangleMatrix(TangleMatrix& tangleMatrix) const
{
    tangleMatrix.m.resize(tangleMatrix.inEdges.size(), vector<uint64_t>(tangleMatrix.outEdges.size()));
    MarkerGraphEdgePairInfo info;
    for(uint64_t i=0; i<tangleMatrix.inEdges.size(); i++) {
        const edge_descriptor ce0 = tangleMatrix.inEdges[i];
        const MarkerGraphEdgeId markerGraphEdgeId0 = secondToLastMarkerGraphEdgeId(ce0);

        for(uint64_t j=0; j<tangleMatrix.outEdges.size(); j++) {
            const edge_descriptor ce1 = tangleMatrix.outEdges[j];
            const MarkerGraphEdgeId markerGraphEdgeId1 = secondToLastMarkerGraphEdgeId(ce1);

            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(markerGraphEdgeId0, markerGraphEdgeId1, info));
            tangleMatrix.m[i][j] = info.common;
        }
    }

}



uint64_t CompressedPathGraph1A::TangleMatrix::inDegree() const
{
    return inEdges.size();
}
uint64_t CompressedPathGraph1A::TangleMatrix::outDegree() const
{
    return outEdges.size();
}



// Analyze the TangleMatrix.
// A matrix element is considered negigible and treated as zero if it is <= lowThreshold.
// A matrix element is considered significantly meaningful if it is >= highThreshold.
// Otherwise a matrix element is considered ambiguous.
void CompressedPathGraph1A::TangleMatrix::analyze(
    uint64_t lowThreshold,
    uint64_t highThreshold,
    int64_t& phase,
    uint64_t& minConcordant,
    uint64_t& maxDiscordant) const
{
    SHASTA_ASSERT(inDegree() == 2);
    SHASTA_ASSERT(inDegree() == 2);

    // Classify matrix elements:
    // 0 = low (<=lowThreshold)
    // 1 = ambiguous (>lowThreshold, <highThreshold
    // 2 = high (>=highThreshold)
    array< array<uint64_t, 2>, 2> c;
    for(uint64_t i=0; i<2; i++) {
        for(uint64_t j=0; j<2; j++) {
            const uint64_t matrixElement = m[i][j];
            uint64_t& classification = c[i][j];
            if(matrixElement <= lowThreshold) {
                classification = 0;
            } else if(matrixElement >= highThreshold) {
                classification = 2;
            } else {
                classification = 1;
            }
        }
    }

    // Check if this tangle matrix is unambiguously in phase.
    if(c[0][0]==2 and c[1][1]==2 and c[0][1]==0 and c[1][0]==0) {
        phase = +1;
        minConcordant = min(m[0][0], m[1][1]);
        maxDiscordant = max(m[0][1], m[1][0]);
    }

    // Check if this tangle matrix is unambiguously out of phase.
    else if(c[0][1]==2 and c[1][0]==2 and c[0][0]==0 and c[1][1]==0) {
        phase = -1;
        minConcordant = min(m[0][1], m[1][0]);
        maxDiscordant = max(m[0][0], m[0][0]);
    }

    // Otherwise, it is ambiguous.
    else {
        phase = 0;
        minConcordant = 0;
        maxDiscordant = 0;
    }
}



uint64_t CompressedPathGraph1A::totalBaseOffset(edge_descriptor ce) const
{
    const CompressedPathGraph1A& cGraph = *this;
    const CompressedPathGraph1AEdge& edge = cGraph[ce];
    const vector<MarkerGraphEdgeId>& chain = edge.chain;

    uint64_t totalOffset = 0;
    MarkerGraphEdgePairInfo info;
    for(uint64_t i=0; i<chain.size(); i++) {
        const MarkerGraphEdgeId edgeId0 = chain[i-1];
        const MarkerGraphEdgeId edgeId1 = chain[i];
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));

        if(info.common > 0 and info.offsetInBases > 0) {
            totalOffset += info.offsetInBases;
        }
    }
    return totalOffset;
}



#if 0
void CompressedPathGraph1A::writeBubble(const Bubble& bubble, ostream& s) const
{
    s << "{";
    for(uint64_t i=0; i<bubble.edges.size(); i++) {
        if(i != 0) {
            s << ",";
        }
        s << edgeStringId(bubble.edges[i]);
    }
    s << "}";
}
#endif

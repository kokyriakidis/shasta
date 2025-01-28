// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "AssemblerOptions.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;



void AssemblyGraph::run4(
    uint64_t threadCount,
    bool /* assembleSequence */,
    bool debug)
{
    cout << "AssemblyGraph::run4 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    write("A");

    // Bubble cleanup.
    compress();
    for(uint64_t iteration=0; ; iteration ++) {
        performanceLog << timestamp << "Iteration " << iteration <<
            " of bubble cleanup begins." << endl;
        const uint64_t cleanedUpBubbleCount = cleanupBubbles(
            debug,
            options.assemblyGraphOptions.bubbleCleanupMaxOffset,
            options.assemblyGraphOptions.chainTerminalCommonThreshold,
            threadCount);
        if(cleanedUpBubbleCount == 0) {
            break;
        }
        cout << "Cleaned up " << cleanedUpBubbleCount << " bubbles." << endl;
        compressBubbleChains();
        compress();
    }

    // Clean up short superbubbles.
    cleanupSuperbubbles(false,
        options.assemblyGraphOptions.superbubbleLengthThreshold1,
        options.assemblyGraphOptions.chainTerminalCommonThreshold);
    compress();

    // Edge detangling.
    expand();
    write("B");

    return;

    for(uint64_t iteration=0; ; iteration++) {
        write("Before-" + to_string(iteration));
        const uint64_t detangledCount = detangleEdges4(true,
            options.assemblyGraphOptions.detangleToleranceHigh);
        write("After-" + to_string(iteration));
        while(compressSequentialEdges());
        compressBubbleChains();
        if(detangledCount == 0) {
            break;
        }
    }
    write("C");

    // Expand, so we can do read following on the AssemblyGraph in the http server.
    expand();
    write("D");

    // Assemble sequence.
    assembleAllChainsMultithreaded(
        options.assemblyGraphOptions.chainTerminalCommonThreshold,
        threadCount);
    writeAssemblyDetails();
    write("Final", true);
}



// This requires all BubbleChains to consist of a single Chain.
// It returns the number of edges that were detangled.
uint64_t AssemblyGraph::detangleEdges4(
    bool debug,
    uint64_t detangleToleranceHigh)
{
    if(debug) {
        cout << "Detangling edges." << endl;
    }

    AssemblyGraph& assemblyGraph = *this;

    // To safely iterate over edges while removing edges we must use edge ids
    // as unique identifiers, because edge descriptors can be reused as edges are
    // deleted and new edges are created.
    std::map<uint64_t, edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgeMap.insert({assemblyGraph[e].id, e});
    }

    uint64_t detangleCount = 0;;
    for(auto it=edgeMap.begin(); it!=edgeMap.end(); /* Incremented safely by detangleEdge4 */) {
        if(detangleEdge4(debug, edgeMap, it, detangleToleranceHigh)) {
            ++detangleCount;
        }
    }

    if(true) {
        cout << "Detangled " << detangleCount << " edges." << endl;
    }

    return detangleCount;

}



// This requires all BubbleChains to consist of a single Chain.
// It returns the true if the edge was detangled and false
// if no changes were made to the AssemblyGraph.
bool AssemblyGraph::detangleEdge4(
    bool debug,
    std::map<uint64_t, edge_descriptor>& edgeMap,
    std::map<uint64_t, edge_descriptor>::iterator& it,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& assemblyGraph = *this;

    // Get the edge to be detangled, then tentatively increment the iterator.
    // However, we may to increment it again below if it ends up pointing to an edge
    // that gets removed.
    const edge_descriptor e = it->second;
    ++it;

    // Sanity check.
    SHASTA_ASSERT(assemblyGraph[e].isSimpleChain());

    // Check that the edge connectivity permits detangling.
    const vertex_descriptor v0 = source(e, assemblyGraph);
    const vertex_descriptor v1 = target(e, assemblyGraph);
    if(out_degree(v0, assemblyGraph) != 1) {
        return false;
    }
    if(in_degree(v1, assemblyGraph) != 1) {
        return false;
    }

    if(debug) {
        cout << "Attempting to detangle edge " << bubbleChainStringId(e) << endl;
    }

    // Gather the in-edges
    // Ignore in-edges coming from cv1 (back-edges).
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> backEdges;
    BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
        if(source(e, assemblyGraph) != v1) {
            inEdges.push_back(e);
        }else {
            backEdges.push_back(e);
        }
    }
    const uint64_t inDegree = inEdges.size();

    // Gather the out-edges.
    // Ignore out-edges going to cv0 (back-edges).
    vector<edge_descriptor> outEdges;
    BGL_FORALL_OUTEDGES(v1, e, assemblyGraph, AssemblyGraph) {
        if(target(e, assemblyGraph) != v0) {
            outEdges.push_back(e);
        }
    }
    const uint64_t outDegree = outEdges.size();

    if((inDegree < 2) or (outDegree < 2)) {
        if(debug) {
            cout << "Not detangling due to degree." << endl;
        }
        return false;
    }

    // If there are common edges between the in-edges and out-edges, skip.
    // The code below does not work for this case.
    for(const edge_descriptor e0: inEdges) {
        if(find(outEdges.begin(), outEdges.end(), e0) != outEdges.end()) {
            if(true) {
                cout << "Not detangling due to cycle." << endl;
            }
            return false;
        }
    }

    // Gather the second to last AnchorId of each inEdge.
    vector<AnchorId> inAnchors;
    for(const edge_descriptor e: inEdges) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        inAnchors.push_back(chain.secondToLast());
    }

    // Gather the second AnchorId of each outEdge.
    vector<AnchorId> outAnchors;
    for(const edge_descriptor e: outEdges) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        outAnchors.push_back(chain.second());
    }

    // If an AnchorId appears both in the inEdges and in the outEdges,
    // detangling could generate a chain with two consecutive copies of the same
    // AnchorId. Don't detangle.
    for(const AnchorId anchorId: inAnchors) {
        if(find(outAnchors.begin(), outAnchors.end(), anchorId) != outAnchors.end()) {
            if(debug) {
                cout << "Not detangling due to cycle." << endl;
            }
            return false;
        }
    }

    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix(inDegree, vector<uint64_t>(outDegree));
    for(uint64_t i0=0; i0<inDegree; i0++) {
        const AnchorId anchorId0 = inAnchors[i0];
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const AnchorId anchorId1 = outAnchors[i1];
            tangleMatrix[i0][i1] = anchors.countCommon(anchorId0, anchorId1);
        }
    }



    if(debug) {
        cout << "Tangle matrix for edge " << bubbleChainStringId(e) << " is " <<
            inDegree << " by " << outDegree << endl;

        cout << "In-edges: ";
        for(const edge_descriptor e: inEdges) {
            cout << " " << bubbleChainStringId(e);
        }
        cout << endl;

        cout << "Out-edges: ";
        for(const edge_descriptor e: outEdges) {
            cout << " " << bubbleChainStringId(e);
        }
        cout << endl;

        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor e0 = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor e1 = outEdges[i1];
                cout <<
                    bubbleChainStringId(e0) << " " <<
                    bubbleChainStringId(e1) << " " <<
                    tangleMatrix[i0][i1] << endl;
            }
        }
    }

    // Gather significant entries of the tangle matrix.
    // These are the ones that are equal to or greater than
    // detangleToleranceHigh.
    vector< pair<uint64_t, uint64_t> > significantPairs;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                significantPairs.push_back(make_pair(i0, i1));
            }
        }
    }

    // If there are no significant pairs, we can't detangle.
    if(significantPairs.empty()) {
        if(debug) {
            cout << "Not detangling because the tangle matrix "
                "contains no significant entries." << endl;
        }
        return false;
    }

    // If all tangle matrix entries are significant, we can't detangle.
    if(significantPairs.size() == inDegree * outDegree) {
        if(debug) {
            cout << "Not detangling because all tangle matrix "
                "entries are significant." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "This edge will be detangled." << endl;
    }

    // Create truncated versions of the inEdges and outEdges.
    vector<vertex_descriptor> inVertices;
    for(const edge_descriptor ce: inEdges) {
        inVertices.push_back(cloneAndTruncateAtEnd(debug, ce));
    }
    vector<vertex_descriptor> outVertices;
    for(const edge_descriptor ce: outEdges) {
        outVertices.push_back(cloneAndTruncateAtBeginning(debug, ce));
    }


    // Each significant element of the tangle matrix generates a new edge.
    for(const auto& p: significantPairs) {
        const uint64_t i0 = p.first;
        const uint64_t i1 = p.second;
        const edge_descriptor eNew = connect(inVertices[i0], outVertices[i1]);
        if(debug) {
            const edge_descriptor e0 = inEdges[i0];
            const edge_descriptor e1 = outEdges[i1];
            cout << "Created " << bubbleChainStringId(eNew) <<
                " by connecting " << bubbleChainStringId(e0) <<
                " with " << bubbleChainStringId(e1) << endl;
        }
    }


    // Now we can remove v0, v1, e, and all of the in-edges and out-edges.
    // We have to do this while safely incrementing the edge iterator to point to the
    // next edge that was not removed.
    // We already incremented the iterator to point past e.
    boost::remove_edge(e, assemblyGraph);
    for(const edge_descriptor e: inEdges) {
        if(it != edgeMap.end() and assemblyGraph[e].id == it->first) {
            ++it;
        }
        edgeMap.erase(assemblyGraph[e].id);
        boost::remove_edge(e, assemblyGraph);
    }
    for(const edge_descriptor e: outEdges) {
        if(it != edgeMap.end() and assemblyGraph[e].id == it->first) {
            ++it;
        }
        edgeMap.erase(assemblyGraph[e].id);
        boost::remove_edge(e, assemblyGraph);
    }
    for(const edge_descriptor e: backEdges) {
        if(it != edgeMap.end() and assemblyGraph[e].id == it->first) {
            ++it;
        }
        edgeMap.erase(assemblyGraph[e].id);
        boost::remove_edge(e, assemblyGraph);
    }
    removeVertex(v0);
    removeVertex(v1);

    return true;
}


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
    AssemblyGraph& assemblyGraph = *this;

    cout << "AssemblyGraph::run4 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    // write("A");

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

    // write("B");

    // Clean up short superbubbles.
    cleanupSuperbubbles(false,
        options.assemblyGraphOptions.superbubbleLengthThreshold1,
        options.assemblyGraphOptions.chainTerminalCommonThreshold);
    compressBubbleChains();
    compress();

    // write("C");

    // For detangling the Assemblyraph needs to be in expanded form.
    expand();

    // Edge detangling.
    Superbubbles superbubbles(assemblyGraph, Superbubbles::FromEdges{});
    write("D");
    detangleShortSuperbubbles4(false, superbubbles);
    write("E");

#if 0

    // Detangling.
    expand();
    write("D");

    detangleEdges(false,
        2, // options.assemblyGraphOptions.detangleToleranceLow,
        8, //options.assemblyGraphOptions.detangleToleranceHigh,
        false, //useBayesianModel,
        options.assemblyGraphOptions.epsilon,
        options.assemblyGraphOptions.minLogP);
    while(compressSequentialEdges());
    compressBubbleChains();

    detangleVertices(false,
        2, // options.assemblyGraphOptions.detangleToleranceLow,
        8, // options.assemblyGraphOptions.detangleToleranceHigh,
        false, //useBayesianModel,
        options.assemblyGraphOptions.epsilon,
        options.assemblyGraphOptions.minLogP);
    while(compressSequentialEdges());
    compressBubbleChains();

    write("E");

    compress();
    write("F");

    // Assemble sequence.
    assembleAllChainsMultithreaded(
        options.assemblyGraphOptions.chainTerminalCommonThreshold,
        threadCount);
    writeAssemblyDetails();
    write("Final", true);
#endif
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

    if(debug) {
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
            if(debug) {
                cout << "Created " << bubbleChainStringId(eNew) <<
                    " by connecting " << bubbleChainStringId(e0) <<
                    " with " << bubbleChainStringId(e1) << endl;
            }
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



void AssemblyGraph::detangleVertices4()
{
    AssemblyGraph& assemblyGraph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if(detangleVertex4(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }
}



// This only detangles pathological vertices with high in-degree
// and out-degree.
bool AssemblyGraph::detangleVertex4(vertex_descriptor v)
{
    AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minDegree = 5;
    const uint64_t minCommon = 3;

    const bool debug = true;

    if(debug) {
        cout << "detangleVertex4 for " << anchorIdToString(assemblyGraph[v].anchorId) <<
            " begins." << endl;

    }

    if(in_degree(v, assemblyGraph) < minDegree) {
        if(debug) {
            cout << "Skipped because in-degree is " << in_degree(v, assemblyGraph) << endl;
        }
        return false;
    }
    if(out_degree(v, assemblyGraph) < minDegree) {
        if(debug) {
            cout << "Skipped because out-degree is " << out_degree(v, assemblyGraph) << endl;
        }
        return false;
    }

    // Gather the inEdges with their source vertices.
    vector< pair<edge_descriptor, vertex_descriptor> > inEdges;
    BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        inEdges.push_back(make_pair(e, v0));
    }

    // Gather the outEdges with their target vertices.
    vector< pair<edge_descriptor, vertex_descriptor> > outEdges;
    BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = target(e, assemblyGraph);
        outEdges.push_back(make_pair(e, v0));
    }

    // If there are common edges between inEdges and outEdges, don't do anything.
    for(const auto& p0: inEdges) {
        const edge_descriptor e0 = p0.first;
        for(const auto& p1: outEdges) {
            const edge_descriptor e1 = p1.first;
            if(e0 == e1) {
                if(debug) {
                    cout << "detangleVertex4 for " << anchorIdToString(assemblyGraph[v].anchorId) <<
                        ": skipped due to cycle." << endl;

                }
                return false;
            }
        }
    }


    // Loop over pairs of inEdges and outEdges.
    // If we find common oriented reads, generate a new edge.
    uint64_t newEdgeCount = 0;
    for(const auto& p0: inEdges) {
        const vertex_descriptor v0 = p0.second;
        const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
        for(const auto& p1: outEdges) {
            const vertex_descriptor v1 = p1.second;
            const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

            if(anchors.countCommon(anchorId0, anchorId1) >= minCommon) {

                // Generate the new edge.
                // Make it of a single Chain connecting v0 to v1.
                edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(v0, v1, assemblyGraph);
                ++newEdgeCount;
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = nextEdgeId++;
                BubbleChain& bubbleChain = newEdge;
                bubbleChain.resize(1);
                Bubble& bubble = bubbleChain.front();
                bubble.resize(1);
                Chain& chain = bubble.front();
                chain.push_back(anchorId0);
                chain.push_back(anchorId1);
            }
        }
    }

    const int64_t deltaEdgeCount = int64_t(newEdgeCount) - int64_t(inEdges.size()) - int64_t(outEdges.size());

    if(debug) {
        cout << "detangleVertex4 for " << anchorIdToString(assemblyGraph[v].anchorId) <<
            ": in-degree " << inEdges.size() <<
            ", out-degree " << outEdges.size() <<
            ", edges created " << newEdgeCount <<
            ", delta edge count " << deltaEdgeCount << endl;
    }


    // Now we can remove the inEdges and outEdges.
    // The vertex we detangled will be removed later by detangleVertices4.
    for(const auto& p: inEdges) {
        boost::remove_edge(p.first, assemblyGraph);
    }
    for(const auto& p: outEdges) {
        boost::remove_edge(p.first, assemblyGraph);
    }

    return true;


}




uint64_t AssemblyGraph::detangleEdges4Strict(bool debug)
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
        if(detangleEdge4Strict(debug, edgeMap, it)) {
            ++detangleCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangleCount << " edges." << endl;
    }

    return detangleCount;

}



bool AssemblyGraph::detangleEdge4Strict(
    bool debug,
    std::map<uint64_t, edge_descriptor>& edgeMap,
    std::map<uint64_t, edge_descriptor>::iterator& it)
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

    if((inDegree != 2) or (outDegree != 2)) {
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


    const bool inPhase =
        (tangleMatrix[0][0] > 0 ) and (tangleMatrix[1][1] >  0) and
        (tangleMatrix[0][1] == 0) and (tangleMatrix[1][0] == 0);
    const bool outOfPhase =
        (tangleMatrix[0][1] > 0 ) and (tangleMatrix[1][0] >  0) and
        (tangleMatrix[0][0] == 0) and (tangleMatrix[1][1] == 0);


    // If neither in phase not out of phase, we can't detangle.
    if(not(inPhase or outOfPhase)) {
        if(debug) {
            cout << "Not detangling because not in phase nor out of phase." << endl;
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


    // Each non-zero element of the tangle matrix generates a new edge.
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            if(tangleMatrix[i0][i1] > 0) {
                const edge_descriptor eNew = connect(inVertices[i0], outVertices[i1]);
                if(debug) {
                    const edge_descriptor e0 = inEdges[i0];
                    const edge_descriptor e1 = outEdges[i1];
                    if(debug) {
                        cout << "Created " << bubbleChainStringId(eNew) <<
                            " by connecting " << bubbleChainStringId(e0) <<
                            " with " << bubbleChainStringId(e1) << endl;
                    }
                }
            }
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



uint64_t AssemblyGraph::detangleShortSuperbubbles4(
    bool debug,
    const Superbubbles& superbubbles)
{
    uint64_t n = 0;

    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(detangleShortSuperbubble4(debug, superbubbles, superbubbleId)) {
            ++n;
        }
    }

    return n;
}



// This only detangles 2 by 2 superbubbles.
// It uses a chi-squared test for phasing.
bool AssemblyGraph::detangleShortSuperbubble4(
    bool debug,
    const Superbubbles& superbubbles,
    uint64_t superbubbleId)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double epsilon = 0.05;
    const double chiSquareThreshold = 10.;

    AssemblyGraph& assemblyGraph = *this;
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

    if(debug) {
        cout << "Found a superbubble with " << superbubble.size() <<
            " vertices:";
        for(const vertex_descriptor cv: superbubble) {
            cout << " " << anchorIdToString(assemblyGraph[cv].getAnchorId());
        }
        cout << endl;
    }

    // Fill in the in-edges and out-edges.
    // These cannot be computed while constructing the superbubbles
    // as they can change when other superbubbles are detangled.
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> outEdges;
    for(const vertex_descriptor v0: superbubble) {
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = source(e, assemblyGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, v1)) {
                 inEdges.push_back(e);
            }
        }
    }
    for(const vertex_descriptor v0: superbubble) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, v1)) {
                 outEdges.push_back(e);
            }
        }
    }
    const uint64_t inDegree = inEdges.size();
    const uint64_t outDegree = outEdges.size();

    if(debug) {
        cout << inDegree << " in-edges:";
        for(const edge_descriptor e: inEdges) {
            cout << " " << bubbleChainStringId(e);
        }
        cout << endl;
        cout << outDegree << " out-edges:";
        for(const edge_descriptor e: outEdges) {
            cout << " " << bubbleChainStringId(e);
        }
        cout << endl;
    }

    // If an inEdge is also an outEdge, don't do anything.
    for(const edge_descriptor e: inEdges) {
        if(find(outEdges.begin(), outEdges.end(), e) != outEdges.end()) {
            if(debug) {
                cout << "Not detangling because " << bubbleChainStringId(e) <<
                    " is both an in-edge and out-edge." << endl;
            }
            return false;
        }
    }

    // This only detangles superbubbles with 2 entrances and 2 exits.
    if(not ((inDegree == 2) and (outDegree == 2))) {
        if(debug) {
            cout << "Not detangling because it is not a 2 by 2 superbubble." << endl;
        }
        return false;
    }

    // Gather the second to last AnchorId of each inEdge and the second AnchorId
    // of each outEdge.
    vector<AnchorId> inAnchors;
    for(const edge_descriptor e: inEdges) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        inAnchors.push_back(chain.secondToLast());
    }
    vector<AnchorId> outAnchors;
    for(const edge_descriptor e: outEdges) {
        const Chain& chain = assemblyGraph[e].getOnlyChain();
        outAnchors.push_back(chain.second());
    }

    if(debug) {
        cout << inDegree << " in-anchors:";
        for(const AnchorId anchorId: inAnchors) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;
        cout << outDegree << " out-anchors:";
        for(const AnchorId anchorId: outAnchors) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;
    }

    // If an AnchorId appears both in the inAnchors and in the outAnchors,
    // detangling could generate a chain with two consecutive copies of the same
    // AnchorId. Don't detangle.
    for(const AnchorId anchorId: inAnchors) {
        if(find(outAnchors.begin(), outAnchors.end(), anchorId) != outAnchors.end()) {
            if(debug) {
                cout << "Not detangling because " << anchorIdToString(anchorId) <<
                    " is both an in-anchor and out-anchor." << endl;
            }
            return false;
        }
    }


    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix(2, vector<uint64_t>(2));
    uint64_t N = 0;
    vector<uint64_t> inCoverage(inDegree, 0);
    vector<uint64_t> outCoverage(inDegree, 0);
    for(uint64_t i0=0; i0<inDegree; i0++) {
        const AnchorId anchorId0 = inAnchors[i0];
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const AnchorId anchorId1 = outAnchors[i1];
            const uint64_t n = anchors.countCommon(anchorId0, anchorId1);
            tangleMatrix[i0][i1] = n;
            N += n;
            inCoverage[i0] += n;
            outCoverage[i1] += n;
        }
    }

    if(debug) {
        cout << "Tangle matrix with total coverage " << N << ":" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrix[i0][i1];

                cout << endl;
            }
        }
        cout << "In-coverage: " << inCoverage[0] << " " << inCoverage[1] << endl;
        cout << "Out-coverage: " << outCoverage[0] << " " << outCoverage[1] << endl;
        SHASTA_ASSERT(inCoverage[0] + inCoverage[1] == N);
        SHASTA_ASSERT(outCoverage[0] + outCoverage[1] == N);
    }

    // If the inCoverage or outCoverage have zero entries, do nothing.
    for(uint64_t i=0; i<inDegree; i++) {
        if(inCoverage[i] == 0) {
            if(debug) {
                cout << "Not detangling because of zero common coverage on in-edge " << bubbleChainStringId(inEdges[i]) << endl;
            }
            return false;
        }
    }
    for(uint64_t i=0; i<outDegree; i++) {
        if(outCoverage[i] == 0) {
            if(debug) {
                cout << "Not detangling because of zero common coverage on out-edge " << bubbleChainStringId(outEdges[i]) << endl;
            }
            return false;
        }
    }


    // Create expected values for the tangle matrix under various assumptions.

    // Random.
    vector< vector<double> > tangleMatrixRandom(2, vector<double>(2, 0.));
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixRandom[i0][i1] = double(inCoverage[i0]) * double(outCoverage[i1]) / double(N);
        }
    }

    // In phase, ideal.
    vector< vector<double> > tangleMatrixInPhaseIdeal(2, vector<double>(2, 0.));
    tangleMatrixInPhaseIdeal[0][0]    = sqrt(double(inCoverage[0]) * double(outCoverage[0]));
    tangleMatrixInPhaseIdeal[1][1]    = sqrt(double(inCoverage[1]) * double(outCoverage[1]));
    // Normalize.
    double inPhaseSum = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            inPhaseSum += tangleMatrixInPhaseIdeal[i0][i1];
        }
    }
    double inPhaseFactor = double(N) / inPhaseSum;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixInPhaseIdeal[i0][i1] *= inPhaseFactor;
        }
    }

    // Out of phase, ideal.
    vector< vector<double> > tangleMatrixOutOfPhaseIdeal(2, vector<double>(2, 0.));
    tangleMatrixOutOfPhaseIdeal[0][1] = sqrt(double(inCoverage[0]) * double(outCoverage[1]));
    tangleMatrixOutOfPhaseIdeal[1][0] = sqrt(double(inCoverage[1]) * double(outCoverage[0]));
    // Normalize.
    double outOfPhaseSum = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            outOfPhaseSum += tangleMatrixOutOfPhaseIdeal[i0][i1];
        }
    }
    double outOfPhaseFactor = double(N) / outOfPhaseSum;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixOutOfPhaseIdeal[i0][i1] *= outOfPhaseFactor;
        }
    }

    // In phase and out of phase, non-ideal.
    vector< vector<double> > tangleMatrixInPhase(2, vector<double>(2, 0.));
    vector< vector<double> > tangleMatrixOutOfPhase(2, vector<double>(2, 0.));
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            tangleMatrixInPhase[i0][i1] = epsilon * tangleMatrixRandom[i0][i1] + (1. - epsilon) * tangleMatrixInPhaseIdeal[i0][i1];
            tangleMatrixOutOfPhase[i0][i1] = epsilon * tangleMatrixRandom[i0][i1] + (1. - epsilon) * tangleMatrixOutOfPhaseIdeal[i0][i1];
        }
    }


    if(false) {
        cout << "Random tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixRandom[i0][i1];

                cout << endl;
            }
        }
        cout << "Ideal in phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixInPhaseIdeal[i0][i1];

                cout << endl;
            }
        }
        cout << "Ideal out of phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixOutOfPhaseIdeal[i0][i1];

                cout << endl;
            }
        }
        cout << "Non-ideal in phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixInPhase[i0][i1];

                cout << endl;
            }
        }
        cout << "Non-ideal out of phase tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrixOutOfPhase[i0][i1];

                cout << endl;
            }
        }
    }


    // Do a chi-square test.
    double chi2InPhase = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const double expected = tangleMatrixInPhase[i0][i1];
            const double delta = double(tangleMatrix[i0][i1]) - expected;
            chi2InPhase += delta * delta / expected;
        }
    }
    double chi2OutOfPhase = 0.;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const double expected = tangleMatrixOutOfPhase[i0][i1];
            const double delta = double(tangleMatrix[i0][i1]) - expected;
            chi2OutOfPhase += delta * delta / expected;
        }
    }

    if(debug) {
        cout << "Chi square test: in phase " << chi2InPhase <<
            ", out of phase " << chi2OutOfPhase << endl;
    }

    const bool isInPhase = (chi2InPhase < chiSquareThreshold) and (chi2OutOfPhase > chiSquareThreshold);
    const bool isOutOfPhase = (chi2InPhase > chiSquareThreshold) and (chi2OutOfPhase < chiSquareThreshold);

    if(not (isInPhase or isOutOfPhase)) {
        if(debug) {
            cout << "Not detangling bacause phasing is not reliable." << endl;
        }
        return false;
    }

    if(debug) {
        if(isInPhase) {
            cout << "In-phase." << endl;
        }
        if(isOutOfPhase) {
            cout << "Out-of-phase." << endl;
        }
    }



    // We are in-phase or out-of-phase. Generate two new edges.
    for(uint64_t i=0; i<2; i++) {

        // Get the two edges to be connected.
        const edge_descriptor e0 = inEdges[i];
        const edge_descriptor e1 = (isInPhase ? outEdges[i] : outEdges[1 - i]);

        // Get the corresponding Chains.
        const Chain& chain0 = assemblyGraph[e0].getOnlyChain();
        const Chain& chain1 = assemblyGraph[e1].getOnlyChain();

        // Get the two vertices for the new edge.
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        const vertex_descriptor v1 = target(e1, assemblyGraph);

        // Create the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = boost::add_edge(v0, v1, assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.id = nextEdgeId++;
        BubbleChain& newBubbleChain = edgeNew;
        newBubbleChain.resize(1);   // The new BubbleChain has a single Bubble
        Bubble& newBubble = newBubbleChain.front();
        newBubble.resize(1);        // The new Bubble is haploid.
        Chain& newChain = newBubble.front();

        // Build the new chain.
        copy(chain0.begin(), chain0.end() - 1, back_inserter(newChain));
        copy(chain1.begin() + 1, chain1.end(), back_inserter(newChain));
    }

    // Now we can remove  all the vertices inside the superbubble
    // and their edges. This includes the inEdges and outEdges.
    for(const vertex_descriptor v: superbubble) {
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }

    return false;
}


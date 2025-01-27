// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-Detangle3.hpp"
#include "findLinearChains.hpp"
#include "orderPairs.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include <queue>



void AssemblyGraph::run3(
    uint64_t /* threadCount */,
    bool /* assembleSequence */,
    bool /* debug */)
{
    cout << "AssemblyGraph::run3 begins for component " << componentId << endl;

    SHASTA_ASSERT(std::is_sorted(orientedReadIds.begin(), orientedReadIds.end()));
    SHASTA_ASSERT(std::is_sorted(anchorIds.begin(), anchorIds.end()));

    write("A");
    detangle3();
    write("B");
}



void AssemblyGraph::detangle3()
{
    AssemblyGraph& assemblyGraph = *this;

    Detangle3Graph detangle3Graph(assemblyGraph);
    cout << "The Detangle3Graph has " << num_vertices(detangle3Graph) <<
        " vertices and " << num_edges(detangle3Graph) << " edges." << endl;
}



Detangle3Graph::Detangle3Graph(AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 0;
    const uint64_t maxCoverage = 1000000;
    const uint64_t maxPruneLength = 10000;

    // Create the graph.
    createVertices(minCoverage, maxCoverage);
    createEdges();
    removeIsolatedVertices();
    prune(maxPruneLength);
    setHasMaximumCommonFlags();
    writeGraphviz("Detangle3Graph-A.dot");

    // Find strong chains that can be used to stitch AssemblyGraph Chains together.
    findStrongChains(maxPruneLength);
    writeGraphviz("Detangle3Graph-C.dot");

    // Update the AssemblyGraph with the strong chains we found.
    updateAssemblyGraph();

#if 0
    // Further cleanup of the Detangle3Graph.
    // Remove edges between vertices of the same strong chain,
    // except for the edges which form the strong chain itself.
    removeInternalStrongChainEdges();
    removeEdgesIncidentInsideStrongChains();
    setHasMaximumCommonFlags();
    writeGraphviz("Detangle3Graph-D.dot");
    transitiveReduction();
    writeGraphviz("Detangle3Graph-E.dot");
#endif

#if 0
    // Debugging only.
    prune(maxPruneLength);
    removeStrongComponents();
    transitiveReduction();
    setHasMaximumCommonFlags();
    removeVeryWeakEdges();
    transitiveReduction();
    writeGraphviz("Detangle3Graph-Z.dot");
#endif
}



void Detangle3Graph::findStrongChains(uint64_t maxPruneLength)
{
    Detangle3Graph& detangle3Graph = *this;

    // Make a copy that will be used to find the strong chains.
    Detangle3Graph simplifiedDetangle3Graph(*this);

    // Clean it up.
    simplifiedDetangle3Graph.removeStrongComponents();
    simplifiedDetangle3Graph.removeWeakEdges();
    simplifiedDetangle3Graph.removeIsolatedVertices();
    simplifiedDetangle3Graph.transitiveReduction();
    simplifiedDetangle3Graph.prune(maxPruneLength);

    simplifiedDetangle3Graph.writeGraphviz("Detangle3Graph-B.dot");

    // Find linear chains on the simplified graph.
    // The strong chains will be subsets of these, to make sure
    // no vertex ends up in more than one strong chain.
    vector< vector<vertex_descriptor> > linearChainsVertices;
    simplifiedDetangle3Graph.findLinearChains(linearChainsVertices);


    // Clear any preexisting strong chains.
    strongChains.clear();
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        detangle3Graph[v].resetStrongChain();
    }

    // Create the strong chains.
    for(const vector<vertex_descriptor>& linearChainVertices: linearChainsVertices) {
        const vertex_descriptor v0 = linearChainVertices.front();
        const vertex_descriptor v1 = linearChainVertices.back();

        // Figure out if v0 and v1 should be included in the strong chain.
        const bool include0 = (out_degree(v0, simplifiedDetangle3Graph) < 2);
        const bool include1 = (in_degree(v1, simplifiedDetangle3Graph) < 2);

        // Compute the length of our strong chain.
        uint64_t strongChainLength = linearChainVertices.size();
        if(not include0) {
            strongChainLength -= 1;
        }
        if(not include1) {
            strongChainLength -= 1;
        }

        // If too short, skip it.
        if(strongChainLength < 2) {
            continue;
        }

        // Generate this strong chain.
        const uint64_t strongChainId = strongChains.size();
        strongChains.resize(strongChains.size() + 1);
        vector<vertex_descriptor>& strongChain = strongChains.back();
        for(uint64_t i=0; i<linearChainVertices.size(); i++) {

            // Skip it if necessary.
            if((i == 0) and not include0) {
                continue;
            }
            if((i == linearChainVertices.size() - 1) and not include1) {
                continue;
            }

            // Find the corresponding vertex in the original Detangle3Graph.
            const vertex_descriptor u = linearChainVertices[i]; // Vertex descriptor in the simplifiedDetangle3Graph
            const AssemblyGraph::edge_descriptor e = simplifiedDetangle3Graph[u].e;
            auto it = vertexMap.find(e);
            SHASTA_ASSERT(it != vertexMap.end());
            const vertex_descriptor v = it->second;     // // Vertex descriptor in the detagle3Graph

            // Store strong chain information in the vertex.
            Detangle3GraphVertex& vertex = detangle3Graph[v];
            vertex.strongChainId = strongChainId;
            vertex.positionInStrongChain = strongChain.size();

            // Add this vertex to the strong chain.
            strongChain.push_back(v);
        }

        cout << "Strong chain:" << endl;
        for(const vertex_descriptor v: strongChain) {
            cout << vertexStringId(v) << " ";
        }
        cout << endl;
    }
}



// Copy constructor.
Detangle3Graph::Detangle3Graph(const Detangle3Graph& that) :
    Detangle3GraphBaseClass(that),
    assemblyGraph(that.assemblyGraph)
{
    Detangle3Graph& detangle3Graph = *this;

    // Create the vertexMap.
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        vertexMap.insert(make_pair(detangle3Graph[v].e, v));
    }

}



void Detangle3Graph::removeWeakEdges()
{
    Detangle3Graph& detangle3Graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphEdge& edge = detangle3Graph[e];
        if(not edge.isStrong()) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle3Graph);
    }
}



void Detangle3Graph::removeVeryWeakEdges()
{
    Detangle3Graph& detangle3Graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphEdge& edge = detangle3Graph[e];
        if(edge.isVeryWeak()) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle3Graph);
    }
}


void Detangle3Graph::removeIsolatedVertices()
{
    Detangle3Graph& detangle3Graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        if(
            (in_degree(v, detangle3Graph) == 0) and
            (out_degree(v, detangle3Graph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        const AssemblyGraph::edge_descriptor e = detangle3Graph[v].e;
        vertexMap.erase(e);
        boost::remove_vertex(v, detangle3Graph);
    }
}



void Detangle3Graph::removeStrongComponents()
{
    Detangle3Graph& detangle3Graph = *this;

    // Map the vertices to integers.
    // This is needed for the computation of strong components below.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        vertexIndexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        detangle3Graph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));

    // Gather the vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents(vertexIndexMap.size());
    for(const auto& p: componentMap) {
        const vertex_descriptor v = p.first;
        const uint64_t componentId = p.second;
        strongComponents[componentId].push_back(v);
    }

    // Remove the vertices in the non-trivial strong components.
    uint64_t strongComponentCount = 0;
    uint64_t removedCount = 0;
    for(const vector<vertex_descriptor>& strongComponent: strongComponents) {
        if(strongComponent.size() > 1) {
            ++strongComponentCount;
            for(const vertex_descriptor v: strongComponent) {
                ++removedCount;
                const AssemblyGraph::edge_descriptor e = detangle3Graph[v].e;
                vertexMap.erase(e);
                boost::remove_vertex(v, detangle3Graph);
            }
        }
    }
    cout << "Removed " << removedCount <<
        " vertices in " << strongComponentCount <<
        " strongly connected components." << endl;
}



// Each assembly graph edge, which must consist of a single chain, can generate
// a vertex of the Detangle3Graph.
void Detangle3Graph::createVertices(uint64_t minCoverage, uint64_t maxCoverage)
{
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        createVertex(e, minCoverage, maxCoverage);
    }
}



void Detangle3Graph::createVertex(
    AssemblyGraph::edge_descriptor e,
    uint64_t minCoverage,
    uint64_t maxCoverage)
{

    Detangle3Graph& detangle3Graph = *this;

    const BubbleChain& bubbleChain = assemblyGraph[e];
    SHASTA_ASSERT(bubbleChain.isSimpleChain());
    const Chain& chain = bubbleChain.getOnlyChain();

    // A Chain must have at least one internal anchor to generate
    // a vertex of the Detangle3Graph.
    if(chain.size() <= 2) {
        return;
    }

    // Compute average coverage of the internal anchors.
    double sum = 0.;
    for(uint64_t i=1; i<chain.size()-1; i++) {
        const AnchorId anchorId = chain[i];
        sum += double(assemblyGraph.anchors[anchorId].size());
    }
    sum /= double(chain.size() - 2);
    const uint64_t coverage = uint64_t(std::round(sum));

    // If coverage is not in the requested range, don't generate a vertex.
    if((coverage < minCoverage) or (coverage>maxCoverage)) {
        return;
    }

    // Compute the total offset between the first and last internal anchors.
    uint64_t offset = 0;
    for(uint64_t i=1; i<chain.size()-2; i++) {
        const AnchorId anchorId0 = chain[i];
        const AnchorId anchorId1 = chain[i + 1];
        AnchorPairInfo info;
        assemblyGraph.anchors.analyzeAnchorPair(anchorId0, anchorId1, info);
        offset += info.offsetInBases;
    }

    // Add the vertex.
    const AnchorId firstInternalAncorId = chain[1];
    const AnchorId lastInternalAncorId = chain[chain.size() - 2];
    const vertex_descriptor v = add_vertex(
        Detangle3GraphVertex(e, firstInternalAncorId, lastInternalAncorId, coverage, offset),
        detangle3Graph);

    // Store it in the vertexMap.
    vertexMap.insert(make_pair(e, v));
}



void Detangle3Graph::createEdges()
{
    Detangle3Graph& detangle3Graph = *this;

    BGL_FORALL_VERTICES(v0, detangle3Graph, Detangle3Graph) {
        createEdges(v0, 0);
        createEdges(v0, 1);
    }

}



// Set the hasMaximumCommon on all edges.
void Detangle3Graph::setHasMaximumCommonFlags()
{
    Detangle3Graph& detangle3Graph = *this;

    array<vector<edge_descriptor>, 2> adjacentEdges;
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        detangle3Graph[e].resetHasMinimumOffsetFlags();
    }
    vector< pair<edge_descriptor, uint64_t> > edgesWithCommonCount;
    BGL_FORALL_VERTICES(v0, detangle3Graph, Detangle3Graph) {

        // Gather the out-edges (direction=0) and in-edges (direction=1).
        adjacentEdges[0].clear();
        BGL_FORALL_OUTEDGES(v0, e, detangle3Graph, Detangle3Graph) {
            adjacentEdges[0].push_back(e);
        }
        adjacentEdges[1].clear();
        BGL_FORALL_INEDGES(v0, e, detangle3Graph, Detangle3Graph) {
            adjacentEdges[1].push_back(e);
        }

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            if(adjacentEdges[direction].empty()) {
                continue;
            }
            edgesWithCommonCount.clear();
            for(const edge_descriptor e: adjacentEdges[direction]) {
                const Detangle3GraphEdge& edge = detangle3Graph[e];
                const uint64_t common = edge.info.common;
                edgesWithCommonCount.push_back(make_pair(e, common));
            }

            sort(edgesWithCommonCount.begin(), edgesWithCommonCount.end(),
                OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
            const uint64_t maxCommon = edgesWithCommonCount.front().second;
            for(const auto& p: edgesWithCommonCount) {
                if(p.second < maxCommon) {
                    break;
                }
                detangle3Graph[p.first].hasMaximumCommon[direction] = true;

            }
        }
    }
}



// Create the edges starting at v0:
// If direction is 0, move forward in the AssemblyGraph.
// If direction is 1, move backward in the AssemblyGraph.
// The code is similar to Mode3Assembler::exploreReadFollowingAssemblyGraph.
void Detangle3Graph::createEdges(vertex_descriptor v0, uint64_t direction)
{
    // EXPOSE WHEN CODE STABILIZES.
    uint64_t minCommon = 4;
    double minJaccard = 0.;
    double minCorrectedJaccard = 0.8;

    Detangle3Graph& detangle3Graph = *this;
    const Detangle3GraphVertex& vertex0 = detangle3Graph[v0];
    const AssemblyGraph::edge_descriptor e0 = vertex0.e;
    const AnchorId anchorId0 = (direction == 0) ? vertex0.lastInternalAncorId : vertex0.firstInternalAncorId;

    // Do a BFS in the AssemblyGraph, starting at e0.
    // It is a forward BFS if direction is 0 and a backward BFS if direction is 1.
    std::queue<AssemblyGraph::edge_descriptor> q;
    q.push(e0);
    std::set<AssemblyGraph::edge_descriptor> s;
    s.insert(e0);



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a Chain.
        const AssemblyGraph::edge_descriptor e1 = q.front();
        q.pop();
        const BubbleChain& bubbleChain1 = assemblyGraph[e1];
        const Chain& chain1 = bubbleChain1.getOnlyChain();
        bool hasInternalAnchors = (chain1.size() > 2);


        // If the Chain has internal anchors, see if we can
        // generate an edge v0->v1 (if direction is 0) or v1->v0 (if direction is 1).
        if((e1 != e0) and hasInternalAnchors) {
            const auto it1 = vertexMap.find(e1);

            if(it1 != vertexMap.end()) {
                const vertex_descriptor v1 = it1->second;
                const Detangle3GraphVertex& vertex1 = detangle3Graph[v1];

                const AnchorId anchorId1 = (direction == 0) ? vertex1.firstInternalAncorId : vertex1.lastInternalAncorId;

                // Analyze this pair of anchors.
                AnchorPairInfo info;
                if(direction == 0) {
                    assemblyGraph.anchors.analyzeAnchorPair(anchorId0, anchorId1, info);
                } else {
                    assemblyGraph.anchors.analyzeAnchorPair(anchorId1, anchorId0, info);
                }



                // If good enough, generate an edge if we don't already have it.
                if(
                    (info.common >= minCommon) and
                    (info.jaccard() >= minJaccard) and
                    (info.correctedJaccard() >= minCorrectedJaccard)
                    ) {

                    if(direction == 0) {
                        // Generate an edge v0->v1 if we don't already have it.
                        edge_descriptor e;
                        bool edgeExists = false;
                        tie(e, edgeExists) = boost::edge(v0, v1, detangle3Graph);
                        if(not edgeExists) {
                            add_edge(v0, v1, Detangle3GraphEdge(info), detangle3Graph);
                        }
                    } else {
                        // Generate an edge v1->v0 if we don't already have it.
                        edge_descriptor e;
                        bool edgeExists = false;
                        tie(e, edgeExists) = boost::edge(v1, v0, detangle3Graph);
                        if(not edgeExists) {
                            add_edge(v1, v0, Detangle3GraphEdge(info), detangle3Graph);
                        }
                    }

                }

                // If there are no common reads with anchorId0, don't continue the BFS past this point
                // (but AssemblyGraph edges still in the queue will continue to be processed.
                if(info.common == 0) {
                    continue;
                }
            }
        }


        // Enqueue the children (if direction=0) or parents (if direction=1).
        if(direction == 0) {
            const AssemblyGraph::vertex_descriptor v2 = target(e1, assemblyGraph);
            BGL_FORALL_OUTEDGES(v2, e2, assemblyGraph, AssemblyGraph) {
                if(not s.contains(e2)) {
                    q.push(e2);
                    s.insert(e2);
                }
            }
        } else {
            const AssemblyGraph::vertex_descriptor v2 = source(e1, assemblyGraph);
            BGL_FORALL_INEDGES(v2, e2, assemblyGraph, AssemblyGraph) {
                if(not s.contains(e2)) {
                    q.push(e2);
                    s.insert(e2);
                }
            }

        }
    }

}



void Detangle3Graph::transitiveReduction()
{
    transitiveReductionAny(*this);
}

string Detangle3Graph::vertexStringId(vertex_descriptor v) const
{
    const Detangle3Graph& detangle3Graph = *this;
    const Detangle3GraphVertex& vertex = detangle3Graph[v];
    const AssemblyGraph::edge_descriptor e = vertex.e;
    return assemblyGraph.bubbleChainStringId(e);

}


void Detangle3Graph::writeGraphviz(const string& fileName) const
{
    const Detangle3Graph& detangle3Graph = *this;

    ofstream dot(fileName);
    dot << "digraph detangle3Graph {\n";



    // Write the vertices.
    BGL_FORALL_VERTICES(v, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphVertex& vertex = detangle3Graph[v];
        const AssemblyGraph::edge_descriptor e = vertex.e;
        const Chain& chain = assemblyGraph[e].getOnlyChain();

        dot << "\"" << vertexStringId(v) << "\"";

        // Begin attributes.
        dot << "[";

        // Label
        dot << "label=\"" <<
            vertexStringId(v) << "\\n" <<
            "n = " << chain.size() - 2 << "\\n" <<
            "c = " << vertex.coverage << "\\n" <<
            "o = " << vertex.offset;
        if(vertex.strongChainId != invalid<uint64_t>) {
            dot << "\\n" << vertex.strongChainId << ":" << vertex.positionInStrongChain;
        }
        dot << "\"";

        // Color.
        if(vertex.strongChainId != invalid<uint64_t>) {
            dot << " style=filled fillcolor=pink";
        }

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot <<";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {
        const Detangle3GraphEdge& edge = detangle3Graph[e];
        const vertex_descriptor v0 = source(e, detangle3Graph);
        const vertex_descriptor v1 = target(e, detangle3Graph);

        // Write the source and target vertices.
        dot << "\"" << vertexStringId(v0) <<
            "\"->\"" << vertexStringId(v1) << "\"";

        // Begin attributes.
        dot << "[";

        // Label.
        dot << "label=\"" <<
            edge.info.common << "\\n" <<
            std::fixed << std::setprecision(2)  << edge.info.correctedJaccard() << "\\n" <<
            edge.info.offsetInBases <<
            "\"";

        // Style.
        if(edge.hasMaximumCommon[0]) {
            if(edge.hasMaximumCommon[1]) {
                // Leave if black.
            } else {
                dot << " color=green";
            }
        } else {
            if(edge.hasMaximumCommon[1]) {
                dot << " color=blue";
            } else {
                dot << " color=red";
            }
        }

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";
    }

    dot << "}\n";
}



void Detangle3Graph::prune(uint64_t maxLength)
{
    for(uint64_t iteration=0; ; ++iteration) {
        writeGraphviz("Detangle3Graph-BeforePruneIteration-" + to_string(iteration) + ".dot");
        cout << "Prune iteration " << iteration << " begins." << endl;
        if(not pruneIteration(maxLength)) {
            break;
        }
    }
}



bool Detangle3Graph::pruneIteration(uint64_t maxLength)
{
    Detangle3Graph& detangle3Graph = *this;
    const bool debug = false;

    // Find linear chains.
    vector< vector<vertex_descriptor> > linearChainsVertices;
    vector< vector<edge_descriptor> > linearChainsEdges;
    findLinearChains(linearChainsVertices, linearChainsEdges);

    std::set<vertex_descriptor> verticesToBeRemoved;
    vector<edge_descriptor> edgesToBeRemoved;

    // Loop over all linear chains.
    for(uint64_t i=0; i<linearChainsVertices.size(); i++) {
        const vector<vertex_descriptor>& linearChainVertices = linearChainsVertices[i];
        const vector<edge_descriptor>& linearChainEdges = linearChainsEdges[i];
        const vertex_descriptor v0 = linearChainVertices.front();
        const vertex_descriptor v1 = linearChainVertices.back();

        if(debug) {
            cout << "Found a linear chain that begins at " << vertexStringId(v0) <<
                " and ends at " << vertexStringId(v1) << endl;
        }

        // Check if it is hanging.
        const bool isHangingBackward = (in_degree(v0, detangle3Graph) == 0) and (out_degree(v0, detangle3Graph) == 1);
        const bool isHangingForward = (out_degree(v1, detangle3Graph) == 0) and (in_degree(v1, detangle3Graph) == 1);

        if(debug) {
            cout << "Hanging flags: " << int(isHangingBackward) << int(isHangingForward) << endl;
        }

        // If it is not hanging, skip it.
        if(not (isHangingBackward or isHangingForward)) {
            continue;
        }

        // Compute the length of the hanging portion. This is the sum of the offsets of its vertices.
        // But the first vertex is excluded if the chain is hanging forward and not backward,
        // and the last vertex is excluded if the chain is hanging backward but not forward.
        uint64_t hangingLength = 0;
        /*
        for(const edge_descriptor e: linearChainEdges) {
            hangingLength += detangle3Graph[e].info.offsetInBases;
        }
        */
        for(uint64_t i=0; i<linearChainVertices.size(); i++) {
            if(isHangingForward and not isHangingBackward and (i == 0)) {
                continue;
            }
            if(isHangingBackward and not isHangingForward and (i == linearChainVertices.size() - 1)) {
                continue;
            }
            const vertex_descriptor v = linearChainVertices[i];
            hangingLength += detangle3Graph[v].offset;
        }
        if(debug) {
            cout << "Hanging length is " << hangingLength << endl;
        }

        // If the hanging length is more than maxLength, don't prune it.
        if(hangingLength > maxLength) {
            continue;
        }

        if(debug) {
            cout << "Pruning linear chain " << vertexStringId(v0) <<
                "..." << vertexStringId(v1) <<
                " with hanging length " << hangingLength << endl;
        }

        // Flag vertices and edges to be removed.
        copy(linearChainEdges.begin(), linearChainEdges.end(), back_inserter(edgesToBeRemoved));
        for(uint64_t i=0; i<linearChainVertices.size(); i++) {
            if(isHangingForward and not isHangingBackward and (i == 0)) {
                continue;
            }
            if(isHangingBackward and not isHangingForward and (i == linearChainVertices.size() - 1)) {
                continue;
            }
            const vertex_descriptor v = linearChainVertices[i];
            if(verticesToBeRemoved.contains(v)) {
                cout << "Assertion failing for " << vertexStringId(v) << endl;
            }
            SHASTA_ASSERT(not verticesToBeRemoved.contains(v));
            verticesToBeRemoved.insert(v);
        }
    }

    // Remove all the vertices and edges we flagged to be removed.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle3Graph);
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        if(debug) {
            cout << "Removing " << vertexStringId(v) << endl;
        }
        boost::clear_vertex(v, detangle3Graph);
        const AssemblyGraph::edge_descriptor e = detangle3Graph[v].e;
        vertexMap.erase(e);
        boost::remove_vertex(v, detangle3Graph);
    }

    return not (verticesToBeRemoved.empty() and edgesToBeRemoved.empty());
}



void Detangle3Graph::findLinearChains(
    vector< vector<vertex_descriptor> >& linearChainsVertices) const
{
    vector< vector<edge_descriptor> >linearChainsEdges;
    findLinearChains(linearChainsVertices, linearChainsEdges);
}



void Detangle3Graph::findLinearChains(
    vector< vector<edge_descriptor> >& linearChainsEdges) const
{
    shasta::findLinearChains(*this, 0, linearChainsEdges);
}



void Detangle3Graph::findLinearChains(
    vector< vector<vertex_descriptor> >& linearChainsVertices,
    vector< vector<edge_descriptor> >& linearChainsEdges) const
{
    const Detangle3Graph& detangle3Graph = *this;

    // Find the edges of the linear chains.
    findLinearChains(linearChainsEdges);

    // Gather the vertices of each linear chain.
    linearChainsVertices.clear();
    for(const vector<edge_descriptor>& linearChainEdges: linearChainsEdges) {

        linearChainsVertices.resize(linearChainsVertices.size() + 1);
        vector<vertex_descriptor>& linearChainVertices = linearChainsVertices.back();
        linearChainVertices.push_back(source(linearChainEdges.front(), detangle3Graph));
        for(const edge_descriptor e: linearChainEdges) {
            linearChainVertices.push_back(target(e, detangle3Graph));
        }
    }
}



// Remove edges between vertices of the same strong chain,
// except for the edges which form the strong chain itself.
void Detangle3Graph::removeInternalStrongChainEdges()
{
    Detangle3Graph& detangle3Graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {

        const vertex_descriptor v0 = source(e, detangle3Graph);
        const vertex_descriptor v1 = target(e, detangle3Graph);

        const Detangle3GraphVertex& vertex0 = detangle3Graph[v0];
        const Detangle3GraphVertex& vertex1 = detangle3Graph[v1];

        // If v0 and v1 are not in the same strong chain, skip this edge.
        if(vertex0.strongChainId == invalid<uint64_t>) {
            continue;
        }
        if(vertex1.strongChainId == invalid<uint64_t>) {
            continue;
        }
        if(vertex0.strongChainId  != vertex1.strongChainId) {
            continue;
        }

        // If getting here, v0 and v1 are in the same strong chain.
        // If they are on adjacent positions, skip this edge.
        if(vertex1.positionInStrongChain == vertex0.positionInStrongChain + 1) {
            continue;
        }

        // If getting here, we will remove this edge.
        edgesToBeRemoved.push_back(e);
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle3Graph);
    }

}



// Remove edges incident on vertices that are internal to a strong chain,
// except for the edges of strong chains themselves.
void Detangle3Graph::removeEdgesIncidentInsideStrongChains()
{
    Detangle3Graph& detangle3Graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, detangle3Graph, Detangle3Graph) {

        const vertex_descriptor v0 = source(e, detangle3Graph);
        const vertex_descriptor v1 = target(e, detangle3Graph);

        const Detangle3GraphVertex& vertex0 = detangle3Graph[v0];
        const Detangle3GraphVertex& vertex1 = detangle3Graph[v1];

        const uint64_t chainId0 = vertex0.strongChainId;
        const uint64_t chainId1 = vertex1.strongChainId;

        // Are v0 and v1 in strong chains?
        const bool isInStrongChain0 = chainId0 != invalid<uint64_t>;
        const bool isInStrongChain1 = chainId1 != invalid<uint64_t>;

        // Are v0 and v1 internal to strong chains?
        bool isInternal0 = false;
        if(isInStrongChain0) {
            const vector<vertex_descriptor>& strongChain = strongChains[chainId0];
            const uint64_t position = vertex0.positionInStrongChain;
            if((position != 0) and (position != strongChain.size() - 1)) {
                isInternal0 = true;
            }
        }
        bool isInternal1 = false;
        if(isInStrongChain1) {
            const vector<vertex_descriptor>& strongChain = strongChains[chainId1];
            const uint64_t position = vertex1.positionInStrongChain;
            if((position != 0) and (position != strongChain.size() - 1)) {
                isInternal1 = true;
            }
        }

        // If neither v0 nor v1 is internal to a strong chain, skip this edge so we will keep it.
        if(not (isInternal0 or isInternal1)) {
            continue;
        }

        // If v0 and v1 are internal to the same internal chain and at adjacent position,
        // skip this edge so we will keep it.
        if(isInternal0 and isInternal1 and (chainId0 == chainId1) and
            (vertex1.positionInStrongChain == vertex0.positionInStrongChain + 1)) {
            continue;
        }

        // If getting here, we know that at least one of v0 and v1 is internal to a strong chain,
        // and they are not adjacent in the same strong chain.
        // So this edge will be removed.
        edgesToBeRemoved.push_back(e);
    }



    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, detangle3Graph);
    }
}



// Use the strong chains to update the AssemblyGraph.
void Detangle3Graph::updateAssemblyGraph()
{
    for(const vector<vertex_descriptor>& strongChain: strongChains) {
        updateAssemblyGraph(strongChain);
    }
}



void Detangle3Graph::updateAssemblyGraph(const vector<vertex_descriptor>& strongChain)
{
    Detangle3Graph& detangle3Graph = *this;

    // Create the new Chain to be added to the assembly graph.
    Chain newChain;
    for(uint64_t i=0; i<strongChain.size(); i++) {

        // Get the AssemblyGraph Chain corresponding to this vertex.
        const vertex_descriptor v = strongChain[i];
        const Detangle3GraphVertex& vertex = detangle3Graph[v];
        const AssemblyGraph::edge_descriptor e = vertex.e;
        const Chain& chain = assemblyGraph[e].getOnlyChain();

        // If this is the first vertex of the strong chain, add the first AnchorId.
        if(i == 0) {
            newChain.push_back(chain.front());
        }

        // Add the internal AnchorIds.
        copy(chain.begin() + 1, chain.end() -1, back_inserter(newChain));

        // If this is the last vertex of the strong chain, add the last AnchorId.
        if(i == strongChain.size() - 1) {
            newChain.push_back(chain.back());
        }
    }



    // Add the new Chain to the assemblyGraph.

    // Find the source and target vertices for the new AssemblyGraph edge.
    const AssemblyGraph::edge_descriptor e0 = detangle3Graph[strongChain.front()].e;
    const AssemblyGraph::edge_descriptor e1 = detangle3Graph[strongChain.back()].e;
    const AssemblyGraph::vertex_descriptor v0 = source(e0, assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = target(e1, assemblyGraph);

    // Create the new AssemblyGraph edge (BubbleChain consisting of a single chain.
    AssemblyGraph::edge_descriptor e;
    tie(e,ignore) = add_edge(v0, v1, assemblyGraph);
    AssemblyGraphEdge& newEdge = assemblyGraph[e];
    newEdge.id = assemblyGraph.nextEdgeId++;

    // Store out new Chain as the only Chain of this BubbleChain.
    BubbleChain& newBubbleChain = newEdge;
    newBubbleChain.resize(1);
    Bubble& newBubble = newBubbleChain.front();
    newBubble.resize(1);
    newBubble.front().swap(newChain);

    cout << "Created " << assemblyGraph.bubbleChainStringId(e) <<
        " by stitching " << assemblyGraph.bubbleChainStringId(e0) <<
        " ... " << assemblyGraph.bubbleChainStringId(e1) <<
        " (" << strongChain.size() << " segments)." << endl;

    // Now remove the AssemblyGraph edges containing the old Chains.
    for(const vertex_descriptor v: strongChain) {
        boost::remove_edge(detangle3Graph[v].e, assemblyGraph);
    }
}


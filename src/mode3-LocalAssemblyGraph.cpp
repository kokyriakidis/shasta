// Shasta.
#include "mode3-LocalAssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>



LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraphPostprocessor& assemblyGraph,
    const vector<ChainIdentifier>& startingChains,
    uint64_t maxDistance) :
    assemblyGraph(assemblyGraph)
{
    addVertices(startingChains, maxDistance);
    addEdges();
}



void LocalAssemblyGraph::addVertices(
    const vector<ChainIdentifier>& startingChains,
    uint64_t maxDistance)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    // We will use a BFS to create the vertices.
    std::queue<vertex_descriptor> q;

    // Initialize the queue by adding vertices corresponding to the initial chains.
    /// These vertices have distance 0.
    for(const ChainIdentifier& chainIdentifier: startingChains) {

        /*
        cout << "Creating starting vertices for Chain " <<
            assemblyGraph.chainStringId(chainIdentifier.e, chainIdentifier.positionInBubbleChain, chainIdentifier.indexInBubble)
            << endl;
        */

        // Add a LocalAssemblyGraphVertex for the source of this Chain.
        {
            const LocalAssemblyGraphVertex vertex = createVertexAtChainSource(chainIdentifier);
            if(not vertexMap.contains(vertex)) {
                const vertex_descriptor v = add_vertex(vertex, localAssemblyGraph);
                vertexMap.insert({vertex, v});
                q.push(v);
                /*
                cout << "Added starting vertex: ";
                writeVertex(vertex, cout);
                cout << endl;
                */
            }
        }

        // Add a LocalAssemblyGraphVertex for the target of this Chain.
        {
            const LocalAssemblyGraphVertex vertex = createVertexAtChainTarget(chainIdentifier);
            if(not vertexMap.contains(vertex)) {
                const vertex_descriptor v = add_vertex(vertex, localAssemblyGraph);
                vertexMap.insert({vertex, v});
                q.push(v);
                /*
                cout << "Added starting vertex: ";
                writeVertex(vertex, cout);
                cout << endl;
                */
            }
        }
    }



    // BFS loop.
    vector<LocalAssemblyGraphVertex> neighbors;
    while(not q.empty()) {

        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];

        /*
        cout << "Dequeued ";
        writeVertex(vertex0, cout);
        cout << endl;
        */

        neighbors.clear();
        getNeighbors(vertex0, neighbors);

        for(const LocalAssemblyGraphVertex& vertex1: neighbors) {
            if(vertexMap.contains(vertex1)) {
                continue;
            }

            /*
            cout << "Found ";
            writeVertex(vertex1, cout);
            cout << endl;
            */

            const vertex_descriptor v1 = add_vertex(vertex1, localAssemblyGraph);
            LocalAssemblyGraphVertex& addedVertex1 = localAssemblyGraph[v1];
            addedVertex1.distance = vertex0.distance + 1;
            vertexMap.insert({vertex1, v1});
            if(addedVertex1.distance < maxDistance) {
                q.push(v1);
            }
        }
    }

    uint64_t id = 0;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        localAssemblyGraph[v].id = id++;
    }


    /*
    cout << "Vertices: " << endl;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        writeVertex(vertex, cout);
        cout << " distance " << vertex.distance << endl;
    }
    */
}



void LocalAssemblyGraph::addEdges()
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    BGL_FORALL_VERTICES(v0, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];

        if(vertex0.isTypeA()) {

            // Vertex0 corresponds to a vertex of the assembly graph.
            // Loop over all outgoing chains.
            const AssemblyGraph::vertex_descriptor av0 = vertex0.v;
            BGL_FORALL_OUTEDGES(av0, ae, assemblyGraph, AssemblyGraph) {
                const AssemblyGraphEdge& aEdge = assemblyGraph[ae];
                const BubbleChain& bubbleChain = aEdge;
                const Bubble& firstBubble = bubbleChain.front();
                for(uint64_t indexInBubble=0; indexInBubble<firstBubble.size(); indexInBubble++) {
                    const LocalAssemblyGraphEdge edge(ChainIdentifier(ae, 0, indexInBubble));
                    if(bubbleChain.size() == 1) {
                        const LocalAssemblyGraphVertex& vertex1(target(ae, assemblyGraph));
                        auto it = vertexMap.find(vertex1);
                        if(it != vertexMap.end()) {
                            const vertex_descriptor v1 = it->second;
                            add_edge(v0, v1, edge, localAssemblyGraph);
                        }
                    } else {
                        const LocalAssemblyGraphVertex vertex1(ae, 0);
                        auto it = vertexMap.find(vertex1);
                        if(it != vertexMap.end()) {
                            const vertex_descriptor v1 = it->second;
                            add_edge(v0, v1, edge, localAssemblyGraph);
                        }

                    }
                }
            }

        } else {

            cout << "Adding edges with source vertex ";
            writeVertex(vertex0, cout);
            cout << endl;

            // Vertex0 does not correspond to a vertex of the assembly graph.
            // Loop over all outgoing chains.
            const edge_descriptor ae = vertex0.e;
            const AssemblyGraphEdge& aEdge = assemblyGraph[ae];
            const BubbleChain& bubbleChain = aEdge;
            const Bubble& nextBubble = bubbleChain[vertex0.positionInBubbleChain + 1];
            for(uint64_t indexInBubble=0; indexInBubble<nextBubble.size(); indexInBubble++) {
                const LocalAssemblyGraphEdge edge(ChainIdentifier(ae, vertex0.positionInBubbleChain + 1, indexInBubble));
                if(vertex0.positionInBubbleChain == bubbleChain.size() - 2) {
                    cout << "Case 1" << endl;
                    const LocalAssemblyGraphVertex& vertex1(target(ae, assemblyGraph));
                    auto it = vertexMap.find(vertex1);
                    if(it != vertexMap.end()) {
                        cout << "Case 1 added" << endl;
                        const vertex_descriptor v1 = it->second;
                        add_edge(v0, v1, edge, localAssemblyGraph);
                    }
                } else {
                    cout << "Case 2" << endl;
                    const LocalAssemblyGraphVertex vertex1(ae, vertex0.positionInBubbleChain + 1);
                    auto it = vertexMap.find(vertex1);
                    if(it != vertexMap.end()) {
                        cout << "Case 2 added." << endl;
                        const vertex_descriptor v1 = it->second;
                        add_edge(v0, v1, edge, localAssemblyGraph);
                    }
                }
            }
        }


    }


    cout << " Edges:" << endl;
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        writeEdge(e, cout);
        cout << endl;
    }
}


LocalAssemblyGraphVertex LocalAssemblyGraph::createVertexAtChainSource(
    const ChainIdentifier& chainIdentifier) const
{
    if(chainIdentifier.positionInBubbleChain == 0) {

        // Type A vertex.
        const AssemblyGraph::vertex_descriptor av = source(chainIdentifier.e, assemblyGraph);
        return LocalAssemblyGraphVertex(av);

    } else {

        // Type B vertex.
        return LocalAssemblyGraphVertex(chainIdentifier.e, chainIdentifier.positionInBubbleChain - 1);

    }
}



LocalAssemblyGraphVertex LocalAssemblyGraph::createVertexAtChainTarget(
    const ChainIdentifier& chainIdentifier) const
{
    const BubbleChain& bubbleChain = assemblyGraph[chainIdentifier.e];

    if(chainIdentifier.positionInBubbleChain == bubbleChain.size() - 1) {

        // Type A vertex.
        const AssemblyGraph::vertex_descriptor av = target(chainIdentifier.e, assemblyGraph);
        return LocalAssemblyGraphVertex(av);

    } else {

        // Type B vertex.
        return LocalAssemblyGraphVertex(chainIdentifier.e, chainIdentifier.positionInBubbleChain);
    }
}



void LocalAssemblyGraph::getNeighbors(
    const LocalAssemblyGraphVertex& vertex,
    vector<LocalAssemblyGraphVertex>& neighbors) const
{
    getChildren(vertex, neighbors);
    getParents(vertex, neighbors);
}



void LocalAssemblyGraph::getChildren(
    const LocalAssemblyGraphVertex& vertex0,
    vector<LocalAssemblyGraphVertex>& children) const
{

    if(vertex0.isTypeA()) {
        BGL_FORALL_OUTEDGES(vertex0.v, ae, assemblyGraph, AssemblyGraph) {
            ChainIdentifier chainIdentifier;
            chainIdentifier.e = ae;
            chainIdentifier.positionInBubbleChain = 0;
            chainIdentifier.indexInBubble = 0;
            children.push_back(createVertexAtChainTarget(chainIdentifier));
        }
    } else {
        ChainIdentifier chainIdentifier;
        chainIdentifier.e = vertex0.e;
        chainIdentifier.positionInBubbleChain = vertex0.positionInBubbleChain + 1;
        chainIdentifier.indexInBubble = 0;
        children.push_back(createVertexAtChainTarget(chainIdentifier));
    }
}



void LocalAssemblyGraph::getParents(
    const LocalAssemblyGraphVertex& vertex0,
    vector<LocalAssemblyGraphVertex>& parents) const
{
    /*
    cout << "LocalAssemblyGraph::getParents called for ";
    writeVertex(vertex0, cout);
    cout << endl;
    */

    if(vertex0.isTypeA()) {
        BGL_FORALL_INEDGES(vertex0.v, ae, assemblyGraph, AssemblyGraph) {
            const BubbleChain& bubbleChain = assemblyGraph[ae];
            ChainIdentifier chainIdentifier;
            chainIdentifier.e = ae;
            chainIdentifier.positionInBubbleChain = bubbleChain.size() - 1;
            chainIdentifier.indexInBubble = 0;
            parents.push_back(createVertexAtChainSource(chainIdentifier));

            /*
            cout << "LocalAssemblyGraph::getParents added ";
            writeVertex(parents.back(), cout);
            cout << endl;
            */
        }
    } else {
        ChainIdentifier chainIdentifier;
        chainIdentifier.e = vertex0.e;
        chainIdentifier.positionInBubbleChain = vertex0.positionInBubbleChain;
        chainIdentifier.indexInBubble = 0;
        parents.push_back(createVertexAtChainSource(chainIdentifier));

        /*
        cout << "LocalAssemblyGraph::getParents added ";
        writeVertex(parents.back(), cout);
        cout << endl;
        */
    }
    cout << "LocalAssemblyGraph::getParents ends" << endl;;
}



void LocalAssemblyGraph::writeVertex(
    const LocalAssemblyGraphVertex& vertex,
    ostream& s) const
{
    if(vertex.isTypeA()) {
        const AssemblyGraphVertex aVertex = assemblyGraph[vertex.v];
        s << "Type A vertex with AnchorId " << anchorIdToString(aVertex.anchorId);
    } else {
        const AssemblyGraphEdge& aEdge = assemblyGraph[vertex.e];
        s << "Type B vertex with BubbleChain " << aEdge.id <<
            " position in BubbleChain " << vertex.positionInBubbleChain;
    }
}


void LocalAssemblyGraph::writeEdge(edge_descriptor e, ostream& s) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;
    const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
    s << assemblyGraph.getChainStringId(edge);
}



void LocalAssemblyGraph::writeGraphviz(ostream& s) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    s << "digraph LocalAssemblyGraph {\n";

    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);

        s << localAssemblyGraph[v0].id << "->" << localAssemblyGraph[v1].id << " [";

        // Label.
        s << "label=\"";
        writeEdge(e, s);
        s << "\"";

        // Tooltip.
        s << " tooltip=\"";
        writeEdge(e, s);
        s << "\"";

        s << "];\n";
    }

    s << "}\n";
}



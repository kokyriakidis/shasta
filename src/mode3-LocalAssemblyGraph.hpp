#pragma once

// Shasta.
#include "mode3-AssemblyGraphPostprocessor.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "tuple.hpp"



// The LocalAssemblyGraph is a local version of the AssemblyGraph for a component.
// Each edge corresponds to a Chain in the AssemblyGraph.
// Each can correspond to either:
// - A vertex of the AssemblyGraph.
// - An anchor at the end of a Bubble that is not the last Bubble of its BubbleChain.

namespace shasta {
    namespace mode3 {
        class LocalAssemblyGraph;
        class LocalAssemblyGraphVertex;
        class LocalAssemblyGraphEdge;

        using LocalAssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            LocalAssemblyGraphVertex,
            LocalAssemblyGraphEdge>;

        class LocalAssemblyGraphDisplayOptions;
    }
}



class shasta::mode3::LocalAssemblyGraphDisplayOptions {
public:
    uint64_t sizePixels;
    string layoutMethod;

    // Vertices.
    double vertexSize;
    bool vertexSizeByCoverage;
    bool vertexLabels;

    // Edges.
    string edgeColoring;
    double edgeThickness;
    double minimumEdgeLength;
    double additionalEdgeLengthPerKb;
    double arrowSize;
    bool edgeLabels;

    // Construct from an html request.
    LocalAssemblyGraphDisplayOptions(const vector<string>& request);

    // Write the form.
    void writeForm(ostream& html) const;
};



// A LocalAssemblyGraphVertex can be of two types, type A and type B.
// - A type A LocalAssemblyGraphVertex corresponds to a vertex of the Assembly graph.
//   In this case only v is set.
// - A type B LocalAssemblyGraphVertex corresponds to an anchor
//   at the end of a Bubble which is not the last Bubble of its BubbleChain.
//   In this case only e and positionInBubbleChain are set.
class shasta::mode3::LocalAssemblyGraphVertex {
public:
    AssemblyGraph::vertex_descriptor v = AssemblyGraph::null_vertex();

    AssemblyGraph::edge_descriptor e;
    uint64_t positionInBubbleChain = invalid<uint64_t>;

    uint64_t distance = 0;
    array<double, 2> xy;

    bool isTypeA() const
    {
        return v != AssemblyGraph::null_vertex();
    }
    bool isTypeB() const
    {
        return v == AssemblyGraph::null_vertex();
    }

    // Constructor for a type A LocalAssemblyGraphVertex.
    LocalAssemblyGraphVertex(AssemblyGraph::vertex_descriptor v) :
        v(v)
    {
    }

    // Default constructor.
    LocalAssemblyGraphVertex() {}

    // Constructor for a type B LocalAssemblyGraphVertex.
    LocalAssemblyGraphVertex(
        AssemblyGraph::edge_descriptor e,
        uint64_t positionInBubbleChain) :
        e(e),
        positionInBubbleChain(positionInBubbleChain)
    {
    }

    bool operator<(const LocalAssemblyGraphVertex& that) const {
        return std::tie(v, e, positionInBubbleChain) < std::tie(that.v, that.e, that.positionInBubbleChain);
    }
};



// A LocalAssemblyGraphEdge always corresponds to a Chain in the AssemblyGraph.
class shasta::mode3::LocalAssemblyGraphEdge : public ChainIdentifier {
public:
    LocalAssemblyGraphEdge(const ChainIdentifier& chainIdentifier) :
        ChainIdentifier(chainIdentifier) {}
    array<double, 2> xy = {0., 0.};
};



class shasta::mode3::LocalAssemblyGraph : public LocalAssemblyGraphBaseClass {
public:
    LocalAssemblyGraph(
        const AssemblyGraphPostprocessor&,
        const vector<ChainIdentifier>& startingChains,
        uint64_t maxDistance,
        const string& assemblyStage);

    void writeHtml(
        ostream& html,
        const LocalAssemblyGraphDisplayOptions&,
        const string& assemblyStage);

private:


    const AssemblyGraphPostprocessor& assemblyGraph;
    uint64_t maxDistance;
    string assemblyStage;

    std::map<LocalAssemblyGraphVertex, vertex_descriptor> vertexMap;

    void addVertices(
        const vector<ChainIdentifier>& startingChains,
        uint64_t maxDistance);
    void addEdges();

    LocalAssemblyGraphVertex createVertexAtChainSource(const ChainIdentifier&) const;
    LocalAssemblyGraphVertex createVertexAtChainTarget(const ChainIdentifier&) const;

    // Use the assemblyGraph to gather neighbors of a given vertex.
    void getNeighbors(const LocalAssemblyGraphVertex&, vector<LocalAssemblyGraphVertex>&) const;
    void getChildren(const LocalAssemblyGraphVertex&, vector<LocalAssemblyGraphVertex>&) const;
    void getParents(const LocalAssemblyGraphVertex&, vector<LocalAssemblyGraphVertex>&) const;

    void writeVertex(const LocalAssemblyGraphVertex&, ostream&) const;
    void writeEdge(edge_descriptor, ostream&) const;

    // Get the AnchorId corresponding to a vertex.
    AnchorId getAnchorId(vertex_descriptor) const;

    // Return exact chain length, if available, or estimated value otherwise.
    uint64_t getChainLength(edge_descriptor) const;

    // Html/svg output using svg output created by Graphviz.
    void writeHtml1(
        ostream& html,
        const LocalAssemblyGraphDisplayOptions&) const;
    void writeGraphviz(
        const string& fileName,
        const LocalAssemblyGraphDisplayOptions&
        ) const;
    void writeGraphviz(
        ostream&,
        const LocalAssemblyGraphDisplayOptions&
        ) const;

    // Html/svg output without using svg output created by Graphviz.
    void writeHtml2(
        ostream& html,
        const LocalAssemblyGraphDisplayOptions&,
        const string& assemblyStage);
    void computeLayout(const LocalAssemblyGraphDisplayOptions&);
    // The bounding box of the computed layout.
    class Box {
    public:
        double xMin;
        double xMax;
        double yMin;
        double yMax;
        double xSize() {return xMax - xMin;}
        double ySize() {return yMax - yMin;}
        void makeSquare();
        void extend(double factor);
    };
    Box boundingBox;
    void computeLayoutBoundingBox();

    void writeVertices(
        ostream& html,
        const LocalAssemblyGraphDisplayOptions&) const;
    double vertexRadius(vertex_descriptor, const LocalAssemblyGraphDisplayOptions&) const;

    void writeEdges(
        ostream& html,
        const LocalAssemblyGraphDisplayOptions&,
        const string& assemblyStage) const;
    void writeSvgControls(
        ostream& html,
        const LocalAssemblyGraphDisplayOptions&) const;
};

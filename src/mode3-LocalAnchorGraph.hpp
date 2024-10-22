#pragma once

// Shasta.
#include "mode3-Anchor.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library;
#include <map>
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class LocalAnchorGraph;
        class LocalAnchorGraphVertex;
        class LocalAnchorGraphEdge;

        using LocalAnchorGraphAssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            LocalAnchorGraphVertex,
            LocalAnchorGraphEdge>;

        class LocalAnchorGraphDisplayOptions;
    }
}



class shasta::mode3::LocalAnchorGraphDisplayOptions {
public:
    uint64_t sizePixels;
    string layoutMethod;

    // Vertices.
    double vertexSize;
    bool vertexSizeByCoverage;
    bool vertexLabels;
    string vertexColoring;
    string similarityMeasure;
    string referenceAnchorIdString;

    // Edges.
    string edgeColoring;
    double edgeThickness;
    bool edgeThicknessByCoverage;
    double minimumEdgeLength;
    double additionalEdgeLengthPerKb;
    double arrowSize;
    bool edgeLabels;

    // Construct from an html request.
    LocalAnchorGraphDisplayOptions(const vector<string>& request);

    // Write the form.
    void writeForm(ostream& html) const;
};



class shasta::mode3::LocalAnchorGraphVertex {
public:
    AnchorId anchorId;
    uint64_t distance;
    LocalAnchorGraphVertex(
        AnchorId anchorId,
        uint64_t distance) :
        anchorId(anchorId),
        distance(distance)
    {}
};


class shasta::mode3::LocalAnchorGraphEdge {
public:
    AnchorPairInfo info;
    uint64_t coverage;

    double coverageLoss() const
    {
        return double(info.common - coverage) / double(info.common);
    }
};



class shasta::mode3::LocalAnchorGraph : public LocalAnchorGraphAssemblyGraphBaseClass {
public:
    LocalAnchorGraph(
        const Anchors&,
        const vector<AnchorId>&,
        uint64_t maxDistance,
        bool filterEdgesByCoverageLoss,
        double maxCoverageLoss);

    const Anchors& anchors;
    uint64_t maxDistance;
    std::map<AnchorId, vertex_descriptor> vertexMap;

    void writeGraphviz(
        const string& fileName,
        const LocalAnchorGraphDisplayOptions&
        ) const;
    void writeGraphviz(
        ostream&,
        const LocalAnchorGraphDisplayOptions&
        ) const;
};

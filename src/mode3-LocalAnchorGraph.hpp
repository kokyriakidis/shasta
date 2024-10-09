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
    }
}



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
};



class shasta::mode3::LocalAnchorGraph : public LocalAnchorGraphAssemblyGraphBaseClass {
public:
    LocalAnchorGraph(
        const Anchors&,
        const vector<AnchorId>&,
        uint64_t maxDistance);

    const Anchors& anchors;
    uint64_t maxDistance;
    std::map<AnchorId, vertex_descriptor> vertexMap;

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
};

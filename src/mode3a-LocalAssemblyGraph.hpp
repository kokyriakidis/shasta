#ifndef SHASTA_MODE3A_LOCAL_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_LOCAL_ASSEMBLY_GRAPH_HPP

// Shasta.
#include "Base.hpp"

// Boost libraries.
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "span.hpp"
#include "string.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3a {

        class LocalAssemblyGraph;
        class LocalAssemblyGraphEdge;
        class LocalAssemblyGraphVertex;

        class AssemblyGraphSnapshot;

        using Point = boost::geometry::model::d2::point_xy<double>;
    }

}


// Classes used to display in the http server a local portion of the AssemblyGraph.
class shasta::mode3a::LocalAssemblyGraphVertex {
public:
    uint64_t vertexId;  // In the AssemblyGraphSnapshot
    uint64_t distance;  // From the start vertex.
    LocalAssemblyGraphVertex(
        uint64_t vertexId,
        uint64_t distance);
    LocalAssemblyGraphVertex();

    // The positions of the auxiliary graph vertices corresponding
    // to this segment.
    vector<Point> position;

    // Unit vectors for the outward pointing tangents at the two ends of the segment.
    // The are computed as averages of the directions of the
    // incoming/outgoing links.
    // They are used to display the links as a cubic spline.
    Point t1;
    Point t2;
};



class shasta::mode3a::LocalAssemblyGraphEdge {
public:
    uint64_t edgeId; // In the AssemblyGraphSnapshot.
    LocalAssemblyGraphEdge(uint64_t edgeId=0) :
        edgeId(edgeId)
        {}
};



class shasta::mode3a::LocalAssemblyGraph :
    public boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
    LocalAssemblyGraphVertex, LocalAssemblyGraphEdge> {
public:

    // If dontDisconnect is true, minLinkCoverage is locally reduced
    // where necessary to avoid the creation of vertices with
    // no incoming or outgoing edges.
    LocalAssemblyGraph(
        const AssemblyGraphSnapshot&,
        uint64_t startVertexId,
        uint64_t maxDistance,
        uint64_t minLinkCoverage,
        bool dontDisconnect);

    // This constructor creates the limited detailed representation
    // of the LocalAssemblyGraph.
    // It follows (only) the reads present in the start vertex.
    // It has vertices onl - no edges.
    LocalAssemblyGraph(
        const AssemblyGraphSnapshot&,
        uint64_t startVertexId);

    const AssemblyGraphSnapshot& assemblyGraphSnapshot;
    uint64_t maxDistance;

    vertex_descriptor addVertex(
        uint64_t segmentId,
        uint64_t distance);

    span<const shasta::Base> vertexCompleteSequence(vertex_descriptor) const;


    class SvgOptions {
    public:

        double sizePixels = 600.;
        string layoutMethod;



        // Segment length and thickness.

        // The display length of a segment is computed as
        // minimumSegmentLength + (n-1) * additionalSegmentLengthPerMarker
        // where n is the path length of the segment, in markers.
        double minimumSegmentLength = 1.;
        double additionalSegmentLengthPerBase = 0.001;

        // The thickness of a segment is computed as
        // minimumSegmentThickness + coverage * additionalSegmentThicknessPerUnitCoverage
        // where coverage is average marker graph edge coverage on the segment path.
        double minimumSegmentThickness = 0.3;
        double additionalSegmentThicknessPerUnitCoverage = 0.000;
        double additionalSegmentThicknessPerJourneyEntry = 0.000;
        uint64_t auxiliaryVertexCountPerSegment = 2;

        // Segment coloring
        string segmentColoring = "random";
        uint64_t referenceSegmentId = 0;
        uint64_t referenceSegmentReplicaIndex = 0;

        // Link length and thickness.
        double linkLength = 1;

        // The display thickness of a link is computed as
        // minimumLinkThickness + (linkCoverage-1) * additionalLinkThicknessPerRead
        double minimumLinkThickness = 0.05;
        double additionalLinkThicknessPerRead = 0.005;

        // Construct the options from an html request.
        SvgOptions(const vector<string>& request);

        // Add rows to the html request form.
        void addFormRows(ostream& html);

    };
    void writeHtml(ostream& html, const SvgOptions&, uint64_t snapshotIndex) const;
    void writeSvg(
        const string& fileName,
        const SvgOptions&,
        uint64_t snapshotIndex) const;
    void writeSvg(
        ostream&,
        const SvgOptions&,
        uint64_t snapshotIndex) const;
    void computeLayout(const SvgOptions&, double timeout);
    void computeSegmentTangents();
    void computeSegmentTangents(vertex_descriptor);

    // The detailed representation shows the details of path entries.
    // Each vertex of the LocalAssemblyGraph is displayed as a graphviz cluster.
    // Each cluster contains one graphviz vertex for each journey entry in the
    // corresponding LocalAssemblyGraph vertex.
    // Journey entries with positions differing by one are joined by
    // directed edges.
    void writeDetailedHtml(ostream& html, double timeout, bool limitedRepresentation) const;
    void writeDetailedDot(ostream& html, bool limitedRepresentation) const;

    // Return the svg color for a segment.
    // All copies of a segment are alway displayed in the same color.
    static string randomSegmentColor(uint64_t segmentId);



    bool haveConsecutivePaths(
        vertex_descriptor v1,
        vertex_descriptor v2) const;

    // Return the average link separation for the Link
    // described by an edge.
    int32_t linkSeparation(edge_descriptor) const;

    // Write the local assembly graph in gfa format.
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
};
#endif


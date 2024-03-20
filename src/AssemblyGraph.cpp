#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode0;

#include "fstream.hpp"
#include "iterator.hpp"



void AssemblyGraph::createMarkerToAssemblyTable(uint64_t markerGrapEdgeCount)
{
    markerToAssemblyTable.beginPass1(markerGrapEdgeCount);
    for(EdgeId assemblyGraphEdgeId=0; assemblyGraphEdgeId<edgeLists.size(); assemblyGraphEdgeId++) {
        const span<EdgeId> chain = edgeLists[assemblyGraphEdgeId];
        for(uint32_t position=0; position!=chain.size(); position++) {
            const EdgeId markerGraphEdgeId = chain[position];
            markerToAssemblyTable.incrementCount(markerGraphEdgeId);
        }
    }
    markerToAssemblyTable.beginPass2();
    for(EdgeId assemblyGraphEdgeId=0; assemblyGraphEdgeId<edgeLists.size(); assemblyGraphEdgeId++) {
        const span<EdgeId> chain = edgeLists[assemblyGraphEdgeId];
        for(uint32_t position=0; position!=chain.size(); position++) {
            const EdgeId markerGraphEdgeId = chain[position];
            markerToAssemblyTable.store(
                markerGraphEdgeId, make_pair(assemblyGraphEdgeId, position));
        }
    }
    markerToAssemblyTable.endPass2();

}



// Close all open data.
void AssemblyGraph::close()
{
    if(vertices.isOpen) {
        vertices.close();
    }

    if(reverseComplementVertex.isOpen) {
        reverseComplementVertex.close();
    }

    if(edges.isOpen) {
        edges.close();
    }

    if(reverseComplementEdge.isOpen) {
        reverseComplementEdge.close();
    }

    if(edgesBySource.isOpen()) {
        edgesBySource.close();
    }

    if(edgesByTarget.isOpen()) {
        edgesByTarget.close();
    }

    if(edgeLists.isOpen()) {
        edgeLists.close();
    }

    if(markerToAssemblyTable.isOpen()) {
        markerToAssemblyTable.close();
    }

    if(sequences.isOpen()) {
        sequences.close();
    }

    if(repeatCounts.isOpen()) {
        repeatCounts.close();
    }

    if(orientedReadsByEdge.isOpen()) {
        orientedReadsByEdge.close();
    }
}


// Close and remove all open data.
void AssemblyGraph::remove()
{
    if(vertices.isOpen) {
        vertices.remove();
    }

    if(reverseComplementVertex.isOpen) {
    	reverseComplementVertex.remove();
    }

    if(edges.isOpen) {
        edges.remove();
    }

    if(reverseComplementEdge.isOpen) {
    	reverseComplementEdge.remove();
    }

    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }

    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }

    if(edgeLists.isOpen()) {
        edgeLists.remove();
    }

    if(markerToAssemblyTable.isOpen()) {
        markerToAssemblyTable.remove();
    }

    if(sequences.isOpen()) {
        sequences.remove();
    }

    if(repeatCounts.isOpen()) {
        repeatCounts.remove();
    }

    if(orientedReadsByEdge.isOpen()) {
        orientedReadsByEdge.remove();
    }
}


// Basic Graphviz output of the global assembly graph.
void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream graphOut(fileName);
    graphOut << "digraph AssemblyGraph {\n";

    // Write the vertices.
    // The label contains the corresponding marker graph vertex id.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        graphOut <<
            vertexId <<
            " [label=\"" <<
            vertexId << "\\n" << vertices[vertexId] <<
             "\"];\n";
    }

    // Write the edges.
    // The label contains the edge id and the number of maker graph edges
    // that correspond to this assembly graph edge.
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        const Edge& edge = edges[edgeId];
        graphOut <<
            edge.source << "->" << edge.target <<
            " [label=\"" << edgeId << "\\n" <<
            edgeLists.size(edgeId) << "\\n" <<
            edge.averageEdgeCoverage <<
            "\"];\n";
    }

    graphOut << "}\n";
}



// Create a csv file that can be loaded in Bandage to color assembled segments
// by similarity (number of common oriented reads) with a given assembled segment.
void AssemblyGraph::colorGfaBySimilarityToSegment(
    EdgeId edgeId0,
    uint64_t minVertexCount,
    uint64_t minEdgeCount)
{
    // Compute the number of common oriented reads with edgeId0.
    vector<uint64_t> commonCount(edges.size(), 0);
    uint64_t maximumValue = 0;
    for(EdgeId edgeId1=0; edgeId1<edges.size(); edgeId1++) {
        commonCount[edgeId1] = commonOrientedReadCount(edgeId0, edgeId1, minVertexCount, minEdgeCount);
        if(edgeId1 != edgeId0) {
            maximumValue = max(maximumValue, commonCount[edgeId1]);
        }
    }

    ofstream csv("Assembly-BothStrands-Color.csv");
    csv << "Id,Number of common oriented reads,Color\n";
    for(EdgeId edgeId1=0; edgeId1<edges.size(); edgeId1++) {
        const uint64_t n = commonCount[edgeId1];

        std::ostringstream color;
        if(edgeId1 == edgeId0) {
            color << "blue";
        } else if(n == 0) {
            color << "grey";
        } else {
            const double ratio = double(n) /double(maximumValue);
#if 0
            const double angle = M_PI_2 * ratio;
            const int red = int(255. * std::cos(angle));
            const int green = int(255. * std::sin(angle));
#endif
            int red, green;
            if(ratio < 0.5) {
                red = 255;
                green = int(510. * ratio);
            } else {
                red = int(510. * (1.-ratio));
                green = 255;
            }
            const int blue = 0;
            color.fill('0');
            color << "#" << hex << std::setw(2) << red;
            color << hex << std::setw(2) << green;
            color << hex << std::setw(2) << blue;
        }

        csv << edgeId1 << ",";
        if(n) {
            csv << n;
        }
        csv << "," << color.str() << "\n";
    }
}



// Compute the number of oriented reads in common between two segments.
uint64_t AssemblyGraph::commonOrientedReadCount(
    EdgeId edgeId0,
    EdgeId edgeId1,
    uint64_t minVertexCount,
    uint64_t minEdgeCount) const
{
    const span<const OrientedReadInfo> info0 = orientedReadsByEdge[edgeId0];
    const span<const OrientedReadInfo> info1 = orientedReadsByEdge[edgeId1];
    uint64_t n = 0;
    auto it0 = info0.begin();
    auto it1 = info1.begin();
    while(it0 != info0.end() and it1 != info1.end()){
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
        } else if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
        } else {
            // Only count it if they have common edges.
            if(
                it0->vertexCount >= minVertexCount and
                it1->vertexCount >= minVertexCount and
                it0->edgeCount   >= minEdgeCount and
                it1->edgeCount   >= minEdgeCount) {
                ++n;
            }
            ++it0;
            ++it1;
        }
    }
    return n;
}



// Find the out-degree or in-degree of a vertex.
// This is not simply the same as counting edgesBySource
// and edgesByTarget, because we have to skip edges
// that were removed.
AssemblyGraph::VertexId AssemblyGraph::inDegree(VertexId vertexId) const
{
    const auto e = edgesByTarget[vertexId];
    VertexId inDegree = 0;
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            ++inDegree;
        }
    }
    return inDegree;
}
AssemblyGraph::VertexId AssemblyGraph::outDegree(VertexId vertexId) const
{
    const auto e = edgesBySource[vertexId];
    VertexId outDegree = 0;
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            ++outDegree;
        }
    }
    return outDegree;
}



// Fill in edgesBySource and edgesByTarget.
void AssemblyGraph::computeConnectivity()
{
    edgesBySource.beginPass1(vertices.size());
    edgesByTarget.beginPass1(vertices.size());
    for(const Edge& edge: edges) {
        edgesBySource.incrementCount(edge.source);
        edgesByTarget.incrementCount(edge.target);
    }
    edgesBySource.beginPass2();
    edgesByTarget.beginPass2();
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        const Edge& edge = edges[edgeId];
        edgesBySource.store(edge.source, edgeId);
        edgesByTarget.store(edge.target, edgeId);
    }
    edgesBySource.endPass2();
    edgesByTarget.endPass2();

    // Make sure edges by source and by target are sorted.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        const auto es = edgesBySource[vertexId];
        const auto et = edgesByTarget[vertexId];
        sort(es.begin(), es.end());
        sort(et.begin(), et.end());
    }

}



// Find incoming/outgoing edges of a vertex
// that were not removed.
// They are returned sorted by edge id.
void AssemblyGraph::findInEdges(VertexId vertexId, vector<EdgeId>& edgeIds) const
{
    const auto e = edgesByTarget[vertexId];
    edgeIds.clear();
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            edgeIds.push_back(edgeId);
        }
    }
}
void AssemblyGraph::findOutEdges(VertexId vertexId, vector<EdgeId>& edgeIds) const
{
    const auto e = edgesBySource[vertexId];
    edgeIds.clear();
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            edgeIds.push_back(edgeId);
        }
    }
}



// Write the AssemblyGraph to GFA without including sequence.
// The sequence length of each edge is written as the number of
// marker graph edges.
// Equivalent functions including output of assembled sequence
// are in class Assembler and should be moved here.
// The GFA 1.0 format is described here:
// https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
void AssemblyGraph::writeGfa1BothStrandsNoSequence(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa1BothStrandsNoSequence(gfa);
}
void AssemblyGraph::writeGfa1BothStrandsNoSequence(ostream& gfa) const
{

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment record for each edge.
    for(EdgeId edgeId=0; edgeId<edgeLists.size(); edgeId++) {

        // If this edge was removed, skip it.
        if(edges[edgeId].wasRemoved()) {
            continue;
        }

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edgeId << "\t";

        // Sequence.
        gfa << "*\t";

        // Sequence length is written out expressed in markers.
        gfa << "LN:i:" << edgeLists.size(edgeId) << "\n";
    }


    // Write GFA links.
    // For each vertex in the assembly graph there is a link for
    // each combination of in-edges and out-edges.
    // Therefore each assembly graph vertex generates a number of
    // links equal to the product of its in-degree and out-degree.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {

        // In-edges.
        const span<const EdgeId> edges0 = edgesByTarget[vertexId];

        // Out-edges.
        const span<const EdgeId> edges1 = edgesBySource[vertexId];

        // Loop over in-edges.
        for(const EdgeId edge0: edges0) {
            if(edges[edge0].wasRemoved()) {
                continue;
            }

            // Loop over out-edges.
            for(const EdgeId edge1: edges1) {
                if(edges[edge1].wasRemoved()) {
                    continue;
                }

                // Write out the link record for this edge,
                // with the CIGAR string left unspecified.
                // Note all links are written with orientation ++.
                gfa << "L\t" <<
                    edge0 << "\t+\t" <<
                    edge1 << "\t+\t*\n";
            }
        }
    }


}



// Return the number of edges that were not removed.
AssemblyGraph::EdgeId AssemblyGraph::edgeCount() const
{
    EdgeId count = 0;
    for(const Edge& edge: edges) {
        if(!edge.wasRemoved()) {
            ++count;
        }
    }
    return count;
}


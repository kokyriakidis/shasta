// Shasta.
#include "Assembler.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"

// Experiments with a directed version of the read graph.

namespace shasta {
    class ReadGraph3;
    class ReadGraph3Vertex;
    class ReadGraph3Edge;

    // The ReadGraph3 has a vertex for each OrientedReadId.
    // The vertex descriptor is OrientedReadId::getValue().
    using ReadGraph3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ReadGraph3Vertex,
        ReadGraph3Vertex>;
}



class shasta::ReadGraph3Vertex {
public:

};



class shasta::ReadGraph3Edge {
public:

};



class shasta::ReadGraph3: public ReadGraph3BaseClass {
public:
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
private:
    static void writeQuotedVertexDescriptor(ostream&, vertex_descriptor);

};



void Assembler::createReadGraph3(uint64_t maxAlignmentCount)
{
    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);

    // Vector to keep the alignments for each read,
    // with their number of markers.
    // Contains pairs(marker count, alignment id).
    vector< pair<uint32_t, uint32_t> > readAlignments;

    const bool debug = true;
    if(debug) {
        cout << "createReadGraph3 begins, maxAlignmentCount " << maxAlignmentCount << endl;
    }


    // Loop over reads.
    for(ReadId readId=0; readId<readCount; readId++) {

        // Gather the alignments for this read, each with its number of markers.
        readAlignments.clear();
        for(const uint32_t alignmentId: alignmentTable[OrientedReadId(readId, 0).getValue()]) {
            const AlignmentData& alignment = alignmentData[alignmentId];
            readAlignments.push_back(make_pair(alignment.info.markerCount, alignmentId));
        }

        // Keep the best maxAlignmentCount.
        if(readAlignments.size() > maxAlignmentCount) {
            std::nth_element(
                readAlignments.begin(),
                readAlignments.begin() + maxAlignmentCount,
                readAlignments.end(),
                std::greater< pair<uint32_t, uint32_t> >());
            readAlignments.resize(maxAlignmentCount);
        }

        // Mark the surviving alignments as to be kept.
        for(const auto& p: readAlignments) {
            const uint32_t alignmentId = p.second;
            keepAlignment[alignmentId] = true;
        }
    }
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);

    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;



    // Create the directed read graph.
    // It has a vertex for each OrientedReadId.
    // The vertex descriptor is OrientedReadId::getValue().
    ReadGraph3 readGraph(2 * readCount);
    for(size_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const bool keepThisAlignment = keepAlignment[alignmentId];
        const AlignmentData& alignment = alignmentData[alignmentId];

        /*
        // Record whether this alignment is used in the read graph.
        alignment.info.isInReadGraph = uint8_t(keepThisAlignment);
        */

        // If this alignment is not used in the read graph, we are done.
        if(not keepThisAlignment) {
            continue;
        }

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(alignment.readIds[0], 0);
        OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

        // Swap them if necessary, depending on the average alignment offset at center.
        double offsetAtCenter = alignment.info.offsetAtCenter();
        if(offsetAtCenter == 0.) {
            offsetAtCenter = 0.01;
        }
        if(offsetAtCenter < 0.) {
            swap(orientedReadId0, orientedReadId1);
            offsetAtCenter = -offsetAtCenter;
        }

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), readGraph);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();

        // Create the edge.
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), readGraph);
    }

    if(debug) {
        cout << "The directed read graph has " << num_vertices(readGraph) <<
            " vertices and " << num_edges(readGraph) << " edges." << endl;
    }
    readGraph.writeGraphviz("ReadGraph3.dot");


    if(debug) {
        cout << "createReadGraph3 ends." << endl;
    }

}



void ReadGraph3::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void ReadGraph3::writeGraphviz(ostream& s) const
{
    s << "digraph ReadGraph3 {\n";

    BGL_FORALL_VERTICES(v, *this, ReadGraph3) {
        writeQuotedVertexDescriptor(s, v);
        s << ";\n";
    }

    BGL_FORALL_EDGES(e, *this, ReadGraph3) {
        const vertex_descriptor v0 = source(e, *this);
        const vertex_descriptor v1 = target(e, *this);
        writeQuotedVertexDescriptor(s, v0);
        s << "->";
        writeQuotedVertexDescriptor(s, v1);
        s << ";\n";
    }

    s << "}\n";
}



void ReadGraph3::writeQuotedVertexDescriptor(ostream& s, vertex_descriptor v)
{
    const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
    s << "\"" << orientedReadId << "\"";
}

// Shasta.
#include "Assembler.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include <random>




namespace shasta {
    class ReadGraph3;
    class ReadGraph3Vertex;
    class ReadGraph3Edge;

    using ReadGraph3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ReadGraph3Vertex,
        ReadGraph3Edge>;
}



class shasta::ReadGraph3Vertex {
public:

    // The strong component this vertex belongs to.
    uint64_t strongComponentId = invalid<uint64_t>;

    // The distances from a starting vertex and its reverse complement.
    uint64_t distance0 = invalid<uint64_t>;
    uint64_t distance1 = invalid<uint64_t>;
};



class shasta::ReadGraph3Edge {
public:
    uint64_t alignmentId;
    ReadGraph3Edge(uint64_t alignmentId = invalid<uint64_t>) : alignmentId(alignmentId) {}
};



// The vertex_descriptor is OrientedReadId::getValue().
class shasta::ReadGraph3: public ReadGraph3BaseClass {
public:

    ReadGraph3(uint64_t n) : ReadGraph3BaseClass(n) {}

    void computeStrongComponents(bool debug);
    void processSelfComplementaryStrongComponents(bool debug);

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    // The vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents;

    // The ids of the self-complementary strong components.
    vector<uint64_t> selfComplementaryStrongComponentIds;

};



void ReadGraph3::computeStrongComponents(bool debug)
{
    ReadGraph3& readGraph = *this;

    // Compute strongly connected components.
    const uint64_t strongComponentCount =
        boost::strong_components(readGraph, boost::get(&ReadGraph3Vertex::strongComponentId, readGraph));


    // Store the vertices in each component.
    strongComponents.clear();
    strongComponents.resize(strongComponentCount);
    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        const uint64_t strongComponentId = readGraph[v].strongComponentId;
        SHASTA_ASSERT(strongComponentId < strongComponentCount);
        strongComponents[strongComponentId].push_back(v);
    }

    // Sort them by vertex descriptor, equivalent to sorting by OrientedReadId.
    for(auto& strongComponent: strongComponents) {
        sort(strongComponent.begin(), strongComponent.end());
    }

    // Histogram the strong component sizes.
    vector<uint64_t> histogram;
    for(const auto& strongComponent: strongComponents) {
        const uint64_t size = strongComponent.size();
        if(size >= histogram.size()) {
            histogram.resize(size + 1, 0);
        }
        ++histogram[size];
    }

    if(debug) {
        cout << "Histogram of strong component sizes:" << endl;
        for(uint64_t size=0; size<histogram.size(); size++) {
            const uint64_t frequency = histogram[size];
            if(frequency) {
                cout << size << "," << frequency << endl;
            }
        }
    }

    // Find the strong components that are self-complementary.
    selfComplementaryStrongComponentIds.clear();
    for(uint64_t strongComponentId=0; strongComponentId<strongComponentCount; strongComponentId++) {
        const auto& strongComponent = strongComponents[strongComponentId];
        SHASTA_ASSERT(not strongComponent.empty());

        // A self-complementary component must have more than one vertex.
        if(strongComponent.size() == 1) {
            continue;
        }

        // The first two ReadIds must be the same, on opposite strands.
        const vertex_descriptor v0 = strongComponent[0];
        const vertex_descriptor v1 = strongComponent[1];
        const OrientedReadId orientedReadId0 = OrientedReadId::fromValue(ReadId(v0));
        const OrientedReadId orientedReadId1 = OrientedReadId::fromValue(ReadId(v1));
        if(orientedReadId0.getReadId() != orientedReadId1.getReadId()) {
            continue;
        }

        if(debug) {
            cout << "Strong component " << strongComponentId << " is self-complementary "
                "and has " << strongComponent.size() << " vertices." << endl;
        }
        selfComplementaryStrongComponentIds.push_back(strongComponentId);

        // Sanity checks.
        SHASTA_ASSERT((strongComponent.size() %2) == 0);
        for(uint64_t i0=0; i0<strongComponent.size(); i0+=2) {
            const uint64_t i1 = i0 + 1;
            const vertex_descriptor v0 = strongComponent[i0];
            const vertex_descriptor v1 = strongComponent[i1];
            const OrientedReadId orientedReadId0 = OrientedReadId::fromValue(ReadId(v0));
            const OrientedReadId orientedReadId1 = OrientedReadId::fromValue(ReadId(v1));
            SHASTA_ASSERT(orientedReadId0.getReadId() == orientedReadId1.getReadId());
            SHASTA_ASSERT(orientedReadId0.getStrand() == 0);
            SHASTA_ASSERT(orientedReadId1.getStrand() == 1);
        }

    }

}


void ReadGraph3::processSelfComplementaryStrongComponents(bool debug)
{
    ReadGraph3& readGraph = *this;

    for(const uint64_t strongComponentId: selfComplementaryStrongComponentIds) {
        const vector<vertex_descriptor>& strongComponent = strongComponents[strongComponentId];
        if(debug) {
            cout << "Processing self-complementary strong component " << strongComponentId <<
                " with " << strongComponent.size() << " vertices." << endl;
        }

        // We will compute distances from the two vertices corresponding
        // to the lowest numbered ReadId.
        const vertex_descriptor v0 = strongComponent[0];
        const vertex_descriptor v1 = strongComponent[1];
        const OrientedReadId orientedReadId0 = OrientedReadId::fromValue(ReadId(v1));
        const OrientedReadId orientedReadId1 = OrientedReadId::fromValue(ReadId(v0));

        /* Use this to start from a manually selected ReadId.
        const ReadId readId = 3471;
        const OrientedReadId orientedReadId0(readId, 1);
        const OrientedReadId orientedReadId1(readId, 0);
        const vertex_descriptor v0 = orientedReadId0.getValue();
        const vertex_descriptor v1 = orientedReadId1.getValue();
        */

        BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
            ReadGraph3Vertex& vertex = readGraph[v];
            vertex.distance0 = invalid<uint64_t>;
            vertex.distance1 = invalid<uint64_t>;
        }

        if(debug) {
            cout << "Computing distances from " << orientedReadId0 << endl;
        }
        const auto visitor0 =  boost::make_bfs_visitor(
            boost::record_distances(get(&ReadGraph3Vertex::distance0, readGraph), boost::on_tree_edge()));
        boost::breadth_first_search(readGraph, v0, visitor(visitor0));

        if(debug) {
            cout << "Computing distances from " << orientedReadId1 << endl;
        }
        const auto visitor1 =  boost::make_bfs_visitor(
            boost::record_distances(get(&ReadGraph3Vertex::distance1, readGraph), boost::on_tree_edge()));
        boost::breadth_first_search(boost::reverse_graph<ReadGraph3>(readGraph), v1, visitor(visitor1));
    }

    // Use the distance information to remove edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, readGraph, ReadGraph3) {
        const vertex_descriptor v0 = source(e, readGraph);
        const vertex_descriptor v1 = target(e, readGraph);
        // const bool debug = (v0 == OrientedReadId(5521,1).getValue()) and (v1 == OrientedReadId(5529,1).getValue());
        const ReadGraph3Vertex& vertex0 = readGraph[v0];
        const ReadGraph3Vertex& vertex1 = readGraph[v1];

        // Remove all edges involving vertices at equal distance.
        if( (vertex0.distance0 != invalid<uint64_t>) and
            (vertex0.distance1 != invalid<uint64_t>) and
            (vertex0.distance0 == vertex0.distance1)) {
            edgesToBeRemoved.push_back(e);
            continue;
        }
        if( (vertex1.distance0 != invalid<uint64_t>) and
            (vertex1.distance1 != invalid<uint64_t>) and
            (vertex1.distance0 == vertex1.distance1)) {
            edgesToBeRemoved.push_back(e);
            continue;
        }

        if(
            (vertex0.distance0 == invalid<uint64_t>) or
            (vertex0.distance1 == invalid<uint64_t>) or
            (vertex1.distance0 == invalid<uint64_t>) or
            (vertex1.distance1 == invalid<uint64_t>)) {
            continue;
        }

        if(
            (vertex0.distance0 < vertex0.distance1) and
            (vertex1.distance0 < vertex1.distance1)) {
            continue;
        }

        if(
            (vertex0.distance0 > vertex0.distance1) and
            (vertex1.distance0 > vertex1.distance1)) {
            continue;
        }


        edgesToBeRemoved.push_back(e);
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, readGraph);
    }

    // Reset.
    strongComponents.clear();
    selfComplementaryStrongComponentIds.clear();
    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        ReadGraph3Vertex& vertex = readGraph[v];
        vertex.strongComponentId = invalid<uint64_t>;
        vertex.distance0 = invalid<uint64_t>;
        vertex.distance1 = invalid<uint64_t>;
    }
}



void ReadGraph3::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void ReadGraph3::writeGraphviz(ostream& s) const
{
    const ReadGraph3& readGraph = *this;

    s << "digraph ReadGraph3 {\n";

    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        const ReadGraph3Vertex& vertex = readGraph[v];
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
        s << "\"" << orientedReadId << "\"";
        s << " [tooltip=\"" << orientedReadId << " " << vertex.strongComponentId << " ";
        if(vertex.distance0 == invalid<uint64_t>) {
            s << "-";
        } else {
            s << vertex.distance0;
        }
        s << " ";
        if(vertex.distance1 == invalid<uint64_t>) {
            s << "-";
        } else {
            s << vertex.distance1;
        }
        s << "\"";


        // Vertex color.
        if(vertex.distance0 == invalid<uint64_t>) {
            if(vertex.distance1 == invalid<uint64_t>) {
                s << " color=black"; // Inaccessible by both.
            } else {
                s << " color=red";    // Accessible by 1 only
            }
        } else {
            if(vertex.distance1 == invalid<uint64_t>) {
                s << " color=green"; // Inaccessible by 0 only.
            } else {
                // Accessible by both
                if(vertex.distance0 < vertex.distance1) {
                    s << " color=LightGreen"; // Closer to 0
                }
                if(vertex.distance0 < vertex.distance1) {
                    s << " color=LightGreen"; // Closer to 0
                } else if(vertex.distance1 < vertex.distance0) {
                    s << " color=pink"; // Closer to 1
                } else {
                    s << " color=purple"; // Equal distance
                }
            }

        }

        s << "]";
        s << ";\n";
    }

    BGL_FORALL_EDGES(e, *this, ReadGraph3) {
        const vertex_descriptor v0 = source(e, *this);
        const vertex_descriptor v1 = target(e, *this);
        const OrientedReadId orientedReadId0 = OrientedReadId::fromValue(ReadId(v0));
        const OrientedReadId orientedReadId1 = OrientedReadId::fromValue(ReadId(v1));
        s << "\"" << orientedReadId0 << "\"->\"" << orientedReadId1 << "\"";
        //
        if(readGraph[v0].strongComponentId == readGraph[v1].strongComponentId) {
            // s << " [color=red]";
        }
        //
        s << ";\n";
    }

    s << "}\n";
}



void Assembler::createReadGraph3(uint64_t /* maxAlignmentCount */)
{
    const bool debug = true;
    using boost::add_vertex;
    using boost::add_edge;

    // Find the number of reads.
    const ReadId readCount = getReads().readCount();

    // The vertex_descriptor is OrientedReadId::getValue().
    ReadGraph3 readGraph(2 * readCount);

    // Initially, each alignment generates two edges.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const AlignmentData& alignment = alignmentData[alignmentId];

        // Get the OrientedReadIds.
        OrientedReadId orientedReadId0(alignment.readIds[0], 0);
        OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(orientedReadId0 < orientedReadId1);

        // Swap them if necessary, depending on the average alignment offset at center.
        if(alignment.info.offsetAtCenter() < 0.) {
            swap(orientedReadId0, orientedReadId1);
        }

        // Create the edge.
        add_edge(orientedReadId0.getValue(), orientedReadId1.getValue(), ReadGraph3Edge(alignmentId), readGraph);

        // Also create the reverse complemented edge.
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        add_edge(orientedReadId1.getValue(), orientedReadId0.getValue(), ReadGraph3Edge(alignmentId), readGraph);
    }

    if(debug) {
        cout << "Initially, the read graph has " << num_vertices(readGraph) <<
            " vertices and " << num_edges(readGraph) << " edges." << endl;
    }

    readGraph.computeStrongComponents(debug);
    if(debug) {
        cout << "The initial read graph has " << readGraph.selfComplementaryStrongComponentIds.size() <<
            " self-complementary strongly connected components." << endl;
        readGraph.writeGraphviz("ReadGraph3-A.dot");
    }

    readGraph.processSelfComplementaryStrongComponents(debug);
    if(debug) {
        cout << "After processing self-complementary strongly connected components, "
            "the read graph has " << num_vertices(readGraph) <<
            " vertices and " << num_edges(readGraph) << " edges." << endl;
    }

    if(debug) {
        readGraph.writeGraphviz("ReadGraph3-B.dot");
    }

    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);
    BGL_FORALL_EDGES(e, readGraph, ReadGraph3) {
        keepAlignment[readGraph[e].alignmentId] = true;
    }

    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);
}


#if 0
// Read graph with one vertex per ReadId (not one per OrientedReadId).
// This can represent the entire read graph or just one of its connected
// components.
namespace shasta {
    class ReadGraph3;
    class ReadGraph3Vertex;
    class ReadGraph3Edge;
    class ReadGraph3ComponentPredicate;

    using ReadGraph3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::undirectedS,
        ReadGraph3Vertex,
        ReadGraph3Edge>;
}



#if 0
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
        ReadGraph3Edge>;
}



class shasta::ReadGraph3Vertex {
public:
    uint64_t strongComponent;
    int64_t distance0 = invalid<int64_t>;
    int64_t distance1 = invalid<int64_t>;
};



class shasta::ReadGraph3Edge {
public:

};



class shasta::ReadGraph3: public ReadGraph3BaseClass {
public:

    void computeStrongComponents();
    vector<uint64_t> strongComponentSize;

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
private:
    static void writeQuotedVertexDescriptor(ostream&, vertex_descriptor);

};
#endif



// In the global read graph the ReadId is the same as the vertex descriptor,
// but for individual connected components this is not the case.
class shasta::ReadGraph3Vertex {
public:
    ReadId readId;
    Strand strand = invalid<Strand>;
    ReadGraph3Vertex(ReadId readId = invalid<ReadId>) : readId(readId) {}
};



class shasta::ReadGraph3Edge {
public:
    uint64_t alignmentId;
    bool isSpanningTreeEdge = false;
    bool isForbidden = false;
    bool isInconsistent = false;
    ReadGraph3Edge(uint64_t alignmentId) : alignmentId(alignmentId) {}
};



class shasta::ReadGraph3: public ReadGraph3BaseClass {
public:
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
private:
};



void Assembler::createReadGraph3(uint64_t maxAlignmentCount)
{
    using boost::add_vertex;
    using boost::add_edge;
    using vertex_descriptor = ReadGraph3::vertex_descriptor;
    using edge_descriptor = ReadGraph3::edge_descriptor;

    // Find the number of reads.
    const ReadId readCount = getReads().readCount();

    // Create a ReadGraph3 to represent the global read graph with one vertex per ReadId
    // (not per OrientedReadId). In the global ReadGraph3, for every vertex descriptor v,
    // readGraph[v].readId == v. That is, the vertex descriptor equals the ReadId.
    ReadGraph3 readGraph;
    for(ReadId readId=0; readId<readCount; readId++) {
        add_vertex(ReadGraph3Vertex(readId), readGraph);
    }

    // Initially, each alignment generates an edge.
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const AlignmentData& alignment = alignmentData[alignmentId];
        add_edge(alignment.readIds[0], alignment.readIds[1], ReadGraph3Edge(alignmentId), readGraph);
    }

    cout << "Initially, the read graph with one vertex per ReadId has " << num_vertices(readGraph) <<
        " vertices and " << num_edges(readGraph) << " edges." << endl;



    // By construction parallel edges can only happen if there are two alignments
    // between the same ReadIds and opposite isSameStrand.
    // These can generate strand-strand contact.
    // Remove them all now.
    class NeighborInfo {
    public:
        edge_descriptor e;
        vertex_descriptor u;
        bool operator<(const NeighborInfo& that) const {
            return u < that.u;
        }
    };
    vector<NeighborInfo> neighbors;
    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        neighbors.clear();
        BGL_FORALL_OUTEDGES(v, e, readGraph, ReadGraph3) {
            const vertex_descriptor u = target(e, readGraph);
            neighbors.push_back({e, u});
        }
        sort(neighbors.begin(), neighbors.end());

        // Look for forbidden pairs.
        for(uint64_t i1=1; i1<neighbors.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const NeighborInfo& neighbor0 = neighbors[i0];
            const NeighborInfo& neighbor1 = neighbors[i1];
            const vertex_descriptor u0 = neighbor0.u;
            const vertex_descriptor u1 = neighbor1.u;
            if(u0 == u1) {
                const edge_descriptor e0 = neighbor0.e;
                const edge_descriptor e1 = neighbor1.e;
                const uint64_t alignmentId0 = readGraph[e0].alignmentId;
                const uint64_t alignmentId1 = readGraph[e1].alignmentId;
                const AlignmentData& alignment0 = alignmentData[alignmentId0];
                const AlignmentData& alignment1 = alignmentData[alignmentId1];
                SHASTA_ASSERT(alignment0.isSameStrand == not (alignment1.isSameStrand));
                readGraph[e0].isForbidden = true;
                readGraph[e1].isForbidden = true;
            }
        }
    }

    // Remove the edges marked as forbidden.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, readGraph, ReadGraph3) {
        if(readGraph[e].isForbidden) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, readGraph);
    }
    cout << "Removed " << edgesToBeRemoved.size() << " forbidden edges." << endl;
    cout << "The read graph with one vertex per ReadId now has " << num_vertices(readGraph) <<
        " vertices and " << num_edges(readGraph) << " edges." << endl;

    readGraph.writeGraphviz("ReadGraph3.dot");





    // Compute connected components.

    // If v=readId is a vertex_descriptor of the global read graph, component[v]
    // will contain the componentId of the connected componenthe vertex belongs to.
    vector<uint64_t> component(readCount);
    const uint64_t componentCount = boost::connected_components(readGraph,
        boost::make_iterator_property_map(component.begin(), get(boost::vertex_index, readGraph)));

    // Each connected component is a ReadGraph3.
    vector<ReadGraph3> components(componentCount);

    // If v=readId is a vertex_descriptor of the global read graph, vertexTable[v]
    // will contain the vertex_descriptor of the corresponding vertex in its
    // connected component given by component[v].
    vector<vertex_descriptor> vertexTable(readCount);

    // Construct the vertices of the connected components.
    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        const uint64_t componentId = component[v];
        const vertex_descriptor vComponent =
            add_vertex(ReadGraph3Vertex(readGraph[v].readId), components[componentId]);
        vertexTable[v] = vComponent;
    }

    // Construct the edges of the connected components.
    BGL_FORALL_EDGES(e, readGraph, ReadGraph3) {
        const vertex_descriptor v0 = source(e, readGraph);
        const vertex_descriptor v1 = target(e, readGraph);
        const uint64_t componentId = component[v0];
        SHASTA_ASSERT(component[v1] == componentId);
        add_edge(vertexTable[v0], vertexTable[v1], readGraph[e], components[componentId]);
    }

    // For each non-trivial connected component, compute a random spanning tree.
    const uint32_t seed = 233;
    std::mt19937 randomSource(seed);
    std::minstd_rand generator(randomSource());
    vector<vertex_descriptor> predecessor(readCount);
    for(uint64_t componentId=0; componentId<componentCount; componentId++) {
        ReadGraph3& component = components[componentId];
        const uint64_t vertexCount = num_vertices(component);
        SHASTA_ASSERT(vertexCount > 0);
        if(vertexCount == 1) {
            continue;
        }
        component.writeGraphviz("ReadGraph3-" + to_string(componentId) + ".dot");

#if 0
        for(uint64_t iteration=0; iteration<100; iteration++) {

            // Create a random spanning tree on this component.
            vector<vertex_descriptor> predecessor(vertexCount, ReadGraph3::null_vertex());
            boost::random_spanning_tree(
                component, generator,
                predecessor_map(make_iterator_property_map(predecessor.begin(), get(boost::vertex_index, component)))
                );

            // Flag the spanning tree edges.
            uint64_t spanningTreeEdgeCount = 0;
            BGL_FORALL_EDGES(e, component, ReadGraph3) {
                const ReadGraph3::vertex_descriptor v0 = source(e, component);
                const ReadGraph3::vertex_descriptor v1 = target(e, component);
                const bool isSpanningTreeEdge = (predecessor[v0] == v1) or (predecessor[v1] == v0);
                component[e].isSpanningTreeEdge = isSpanningTreeEdge;
                if(isSpanningTreeEdge) {
                    ++spanningTreeEdgeCount;
                }
            }
            cout << "Found a connected component with " << num_vertices(component) << " vertices, " <<
                num_edges(component) << " edges and " << spanningTreeEdgeCount <<
                " spanning tree edges." << endl;

            // Sanity check on the spanning tree.
            SHASTA_ASSERT(spanningTreeEdgeCount == num_vertices(component) - 1);


            // Do a BFS on the spanning tree to assign strands to vertices.
            BGL_FORALL_VERTICES(v, component, ReadGraph3) {
                component[v].strand = invalid<Strand>;
            }
            component[0].strand = 0;
            std::queue<vertex_descriptor> q;
            q.push(0);
            while(not q.empty()) {
                const vertex_descriptor v0 = q.front();
                q.pop();
                const Strand strand0 = component[v0].strand;
                SHASTA_ASSERT(strand0 != invalid<Strand>);
                BGL_FORALL_OUTEDGES(v0, e, component, ReadGraph3) {
                    const ReadGraph3Edge& edge = component[e];
                    if(not edge.isSpanningTreeEdge) {
                        continue;
                    }
                    const vertex_descriptor v1 = target(e, component);
                    if(component[v1].strand != invalid<Strand>) {
                        continue;
                    }
                    if(alignmentData[component[e].alignmentId].isSameStrand) {
                        component[v1].strand = strand0;
                    } else {
                        component[v1].strand = 1 - strand0;
                    }
                    q.push(v1);
                }
            }
            BGL_FORALL_VERTICES(v, component, ReadGraph3) {
                SHASTA_ASSERT(component[v].strand != invalid<Strand>);
            }

            // Flag and count the inconsistent edges.
            uint64_t inconsistentEdgeCount = 0;
            BGL_FORALL_EDGES(e, component, ReadGraph3) {
                const vertex_descriptor v0 = source(e, component);
                const vertex_descriptor v1 = target(e, component);
                const bool isSameStrand01 = (component[v0].strand == component[v1].strand);
                const bool isGood = alignmentData[component[e].alignmentId].isSameStrand == isSameStrand01;
                if(not isGood) {
                    ++inconsistentEdgeCount;
                }
                component[e].isInconsistent = not isGood;
            }
            cout << "Iteration " << iteration << " found " << inconsistentEdgeCount << " inconsistent edges." << endl;
            // component.writeGraphviz("ReadGraph3-" + to_string(componentId) + "-" + to_string(iteration) + ".dot");
        }
#endif
        SHASTA_ASSERT(0);
    }

#if 0
    // Create a random spanning tree.
    // const uint32_t seed = 231;
    // std::mt19937 randomGenerator(seed);
    // std::uniform_real_distribution<double> uniformDistribution;
    std::random_device rand;
    std::minstd_rand generator(rand());
    vector<ReadGraph3::vertex_descriptor> predecessor;
    random_spanning_tree(
        readGraph, generator,
        // predecessor_map(make_iterator_property_map(predecessor.begin(), get(boost::vertex_index, readGraph)))
        predecessor_map(boost::get(&ReadGraph3Vertex::predecessor, readGraph))
        );

    // Flag the spanning tree edges.
    BGL_FORALL_EDGES(e, readGraph, ReadGraph3) {
        const ReadGraph3::vertex_descriptor v0 = source(e, readGraph);
        const ReadGraph3::vertex_descriptor v1 = target(e, readGraph);
        const ReadGraph3Vertex& vertex0 = readGraph[v0];
        const ReadGraph3Vertex& vertex1 = readGraph[v1];
        readGraph[e].isSpanningTreeEdge = (vertex0.predecessor == v1) or (vertex1.predecessor == v0);
    }

    readGraph.writeGraphviz("ReadGraph3.dot");
#endif

#if 0
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


    // Conmpute distances from the two vertices corresponding ot a given ReadId.
    {
        const ReadId readId = 0;
        const OrientedReadId orientedReadId0(readId, 1);
        const OrientedReadId orientedReadId1(readId, 0);

        auto distance0Map = boost::get(&ReadGraph3Vertex::distance0, readGraph);
        auto distance1Map = boost::get(&ReadGraph3Vertex::distance1, readGraph);

        boost::breadth_first_search(readGraph, orientedReadId0.getValue(),
            boost::visitor(
            boost::make_bfs_visitor(
            boost::record_distances(distance0Map, boost::on_tree_edge()))));
        boost::reverse_graph<ReadGraph3> reverseGraph(readGraph);
        boost::breadth_first_search(reverseGraph, orientedReadId1.getValue(),
            boost::visitor(
            boost::make_bfs_visitor(
            boost::record_distances(distance1Map, boost::on_tree_edge()))));
    }



    readGraph.computeStrongComponents();
    readGraph.writeGraphviz("ReadGraph3.dot");
#endif


}


#if 0
void ReadGraph3::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void ReadGraph3::writeGraphviz(ostream& s) const
{
    const ReadGraph3& readGraph = *this;

    s << "digraph ReadGraph3 {\n";

    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        const ReadGraph3Vertex& vertex = readGraph[v];
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
        writeQuotedVertexDescriptor(s, v);
        /*
        const uint64_t strongComponent = readGraph[v].strongComponent;
        if(strongComponentSize[strongComponent] > 1) {
            s << " [color=red]";
        }
        */

        s << "[";
        s << "label=\"" <<  orientedReadId << " " << vertex.distance0 << " " << vertex.distance1 << "\"";

        if((vertex.distance0 != invalid<int64_t>) and (vertex.distance1 != invalid<int64_t>)) {
            if(vertex.distance0 - vertex.distance1 >= 3) {
                s << " color=green";
            } else if(vertex.distance1 - vertex.distance0 >= 3) {
                s << " color=red";
            }
        }
        s << "];\n";
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



void ReadGraph3::computeStrongComponents()
{
    ReadGraph3& readGraph = *this;

    // Compute strongly connected components.
    const uint64_t strongComponentCount =
        boost::strong_components(readGraph, boost::get(&ReadGraph3Vertex::strongComponent, readGraph));

    // Store the size of each component.
    strongComponentSize.clear();
    strongComponentSize.resize(strongComponentCount, 0);
    BGL_FORALL_VERTICES(v, readGraph, ReadGraph3) {
        ++strongComponentSize[readGraph[v].strongComponent];
    }

    // Histogram the strongly component sizes.
    vector<uint64_t> histogram;
    for(const uint64_t size: strongComponentSize) {
        if(size >= histogram.size()) {
            histogram.resize(size + 1, 0);
        }
        ++histogram[size];
    }
    cout << "Histogram of strong component sizes:" << endl;
    for(uint64_t size=0; size<histogram.size(); size++) {
        const uint64_t frequency = histogram[size];
        if(frequency) {
            cout << size << "," << frequency << endl;
        }
    }

}
#endif
#endif


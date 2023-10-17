#include "enumeratePaths.hpp"
using namespace shasta;

#include <boost/graph/adjacency_list.hpp>

#include "iostream.hpp"


void shasta::testEnumeratePaths()
{
    using Graph = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::bidirectionalS>;
    Graph graph(10);
    using edge_descriptor = Graph::edge_descriptor;
    using Path = vector<edge_descriptor>;

    class PathInspector {
    public:
        const Graph& graph;
        uint64_t length;
        PathInspector(const Graph& graph, uint64_t length) : graph(graph), length(length) {}
        void operator()(const Path& path)
        {
            if(path.size() == length) {
                for(const edge_descriptor e: path) {
                    cout << source(e, graph) << "->" << target(e, graph) << " ";
                }
                cout << endl;
            }
        }
    };
    const uint64_t length = 4;
    PathInspector pathInspector(graph, length);

    add_edge(0, 1, graph);
    add_edge(1, 2, graph);
    add_edge(2, 3, graph);
    add_edge(3, 4, graph);
    add_edge(1, 5, graph);
    add_edge(5, 6, graph);
    add_edge(6, 7, graph);
    add_edge(6, 3, graph);
    add_edge(3, 8, graph);
    add_edge(7, 9, graph);

    enumeratePaths(graph, 0, length, pathInspector);

}

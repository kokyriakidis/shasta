#include "longestPath.hpp"
#include "iostream.hpp"
using namespace shasta;

void shasta::testLongestPath()
{
    using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS>;
    Graph graph(7);
    add_edge(0, 1, graph);
    add_edge(1, 2, graph);
    add_edge(2, 3, graph);
    add_edge(4, 1, graph);
    add_edge(2, 5, graph);
    add_edge(6, 4, graph);

    vector<Graph::vertex_descriptor> longestPath;
    shasta::longestPath(graph, longestPath);

    for(const auto v: longestPath) {
        cout << v << " ";
    }
    cout << endl;
}

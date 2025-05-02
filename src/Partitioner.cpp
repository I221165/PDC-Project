#include "Partitioner.h"
#include "utils.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <algorithm>

using namespace boost;
typedef adjacency_list<vecS, vecS, directedS> Graph;

void computePartition(
    const CSRGraph &G,
    std::vector<int> &scc,
    std::vector<int> &levels,
    std::vector<int> &cac_id)
{
    int n = G.n;
    Graph g(n);
    for (auto &e : G.edges)
        add_edge(e.src, e.dst, g);

    scc.assign(n, 0);
    int num_scc = strong_components(
        g,
        make_iterator_property_map(scc.begin(), get(vertex_index, g))
    );

    // build component graph
    Graph cg(num_scc);
    for (auto &e : G.edges) {
        int u = scc[e.src], v = scc[e.dst];
        if (u != v) add_edge(u, v, cg);
    }

    // topo sort
    std::vector<int> topo;
    topological_sort(cg, std::back_inserter(topo));

    std::vector<int> clevel(num_scc);
    for (int i = 0; i < (int)topo.size(); ++i)
        clevel[topo[i]] = i;

    levels.assign(n, 0);
    for (int i = 0; i < n; ++i)
        levels[i] = clevel[scc[i]];

    // for now, CAC = SCC
    cac_id = scc;
}

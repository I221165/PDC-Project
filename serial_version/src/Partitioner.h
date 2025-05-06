#pragma once
#include "GraphLoader.h"
#include <vector>

namespace psaim {

/**
 * Partition G into SCC and CAC communities, build the CAC-DAG,
 * assign a level to each node's community, and report timing.
 *
 * @param G         Input CSRGraph
 * @param scc       Output: scc[u] = component ID if node u is in an SCC; -1 otherwise
 * @param cac_id    Output: cac_id[u] = component ID if node u is in a CAC; -1 otherwise
 * @param levels    Output: levels[u] = topological level of u's community in the CAC-DAG
 * @param cac_dag   Output: adjacency list of the CAC-DAG (size = #components)
 */
void computePartition(
    const CSRGraph &G,
    std::vector<int> &scc,
    std::vector<int> &cac_id,
    std::vector<int> &levels,
    std::vector<std::vector<int>> &cac_dag
);

} // namespace psaim 
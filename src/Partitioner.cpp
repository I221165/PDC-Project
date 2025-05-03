// src/Partitioner.cpp

#include "Partitioner.h"
#include "utils.h"
#include <stack>
#include <vector>
#include <functional>
#include <queue>
#include <mpi.h>

void computePartition(
    const CSRGraph &G,
    std::vector<int> &scc,
    std::vector<int> &levels,
    std::vector<int> &cac_id)
{
    int n = G.n;

    // 1) Prepare arrays for Tarjan’s SCC/CAC
    scc.assign(n, -1);
    cac_id.assign(n, -1);

    std::vector<int> index(n, -1), lowlink(n, 0), onStack(n, 0), comp_id(n, -1);
    std::stack<int> st;
    int idx = 0, compCount = 0;

    // 2) Recursive DFS lambda
    std::function<void(int)> dfs = [&](int u) {
        index[u] = lowlink[u] = idx++;
        st.push(u);
        onStack[u] = 1;

        // Explore all neighbors via CSR
        for (int ei = G.xadj[u]; ei < G.xadj[u+1]; ++ei) {
            int v = G.edges[ei].dst;
            if (index[v] < 0) {
                dfs(v);
                lowlink[u] = std::min(lowlink[u], lowlink[v]);
            }
            else if (onStack[v]) {
                lowlink[u] = std::min(lowlink[u], index[v]);
            }
        }

        // If u is root of SCC/CAC
        if (lowlink[u] == index[u]) {
            std::vector<int> comp;
            int w;
            do {
                w = st.top(); st.pop();
                onStack[w] = 0;
                comp_id[w] = compCount;
                comp.push_back(w);
            } while (w != u);

            // If size>1 ⇒ true SCC
            if (comp.size() > 1) {
                for (int x : comp) scc[x] = compCount;
            }
            // Every comp element is in a CAC
            for (int x : comp) cac_id[x] = compCount;
            ++compCount;
        }
    };

    // Launch DFS on all unvisited nodes
    for (int u = 0; u < n; ++u) {
        if (index[u] < 0) dfs(u);
    }

    // 3) Build the inter-component DAG
    //    Nodes = 0..compCount-1
    std::vector<std::vector<int>> dag(compCount);
    std::vector<int> indegree(compCount, 0);

    for (int u = 0; u < n; ++u) {
        int cu = comp_id[u];
        for (int ei = G.xadj[u]; ei < G.xadj[u+1]; ++ei) {
            int v = G.edges[ei].dst;
            int cv = comp_id[v];
            if (cu != cv) {
                dag[cu].push_back(cv);
                ++indegree[cv];
            }
        }
    }

    // 4) Topological BFS (Kahn) to compute comp-levels
    std::queue<int> q;
    levels.assign(compCount, 0);
    for (int c = 0; c < compCount; ++c) {
        if (indegree[c] == 0) q.push(c);
    }
    while (!q.empty()) {
        int c = q.front(); q.pop();
        for (int nx : dag[c]) {
            levels[nx] = std::max(levels[nx], levels[c] + 1);
            if (--indegree[nx] == 0) q.push(nx);
        }
    }

    // 5) Map component-levels back to each node
    std::vector<int> nodeLevel(n);
    for (int u = 0; u < n; ++u) {
        nodeLevel[u] = levels[comp_id[u]];
    }
    levels.swap(nodeLevel);

    // Log summary
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        log("Partitioning complete: "
            + std::to_string(compCount)
            + " components, max level "
            + std::to_string(*std::max_element(levels.begin(), levels.end())));
    }
}

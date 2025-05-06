#include "Partitioner.h"
#include <stack>
#include <vector>
#include <functional>
#include <queue>
#include <algorithm>
#include <string>
#include <iostream>

namespace psaim {

void log(const std::string& msg) {
    std::cout << msg << std::endl;
}

void computePartition(
    const CSRGraph &G,
    std::vector<int> &scc,
    std::vector<int> &cac_id,
    std::vector<int> &levels,
    std::vector<std::vector<int>> &cac_dag
) {
    int n = G.n;
    scc   .assign(n, -1);
    cac_id.assign(n, -1);

    // Tarjan state
    std::vector<int> index(n, -1), lowlink(n), onStack(n), comp_id(n, -1);
    std::stack<int> S;
    int idx = 0, compCount = 0;
    std::vector<int> compLevel;  // will be filled after DAG

    // Recursive DFS for SCC/CAC
    std::function<void(int)> dfs = [&](int u) {
        index[u]   = lowlink[u] = idx++;
        onStack[u] = 1;
        S.push(u);

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

        if (lowlink[u] == index[u]) {
            std::vector<int> comp;
            int w;
            do {
                w = S.top(); S.pop();
                onStack[w] = 0;
                comp_id[w] = compCount;
                comp.push_back(w);
            } while (w != u);

            bool isSCC = (comp.size() > 1);
            for (int x : comp) {
                if (isSCC)         scc[x]    = compCount;
                else /* singleton */ cac_id[x] = compCount;
            }
            ++compCount;
        }
    };

    // Run DFS to find components
    for (int u = 0; u < n; ++u) {
        if (index[u] < 0) dfs(u);
    }
    compLevel.assign(compCount, 0);

    // Build CAC-DAG
    cac_dag.assign(compCount, {});
    for (int u = 0; u < n; ++u) {
        int cu = comp_id[u];
        for (int ei = G.xadj[u]; ei < G.xadj[u+1]; ++ei) {
            int cv = comp_id[G.edges[ei].dst];
            if (cu != cv)
                cac_dag[cu].push_back(cv);
        }
    }

    // Remove duplicates from CAC-DAG edges
    for (auto& edges : cac_dag) {
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    }

    // Topological leveling (Kahn's algorithm)
    std::vector<int> indegree(compCount, 0);
    for (int c = 0; c < compCount; ++c)
        for (int v : cac_dag[c])
            ++indegree[v];

    std::queue<int> Q;
    for (int c = 0; c < compCount; ++c)
        if (indegree[c] == 0) Q.push(c);

    while (!Q.empty()) {
        int c = Q.front(); Q.pop();
        for (int v : cac_dag[c]) {
            compLevel[v] = std::max(compLevel[v], compLevel[c] + 1);
            if (--indegree[v] == 0) Q.push(v);
        }
    }

    // Assign per-node levels
    levels.assign(n, 0);
    for (int u = 0; u < n; ++u) {
        levels[u] = compLevel[ comp_id[u] ];
    }

    log("Partitioning complete: "
        + std::to_string(compCount) + " communities, "
        + std::to_string(n)         + " nodes");
}

} // namespace psaim 
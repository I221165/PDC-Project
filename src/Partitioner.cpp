// ──────────────────────────────────────────────────────────────────────────────
// File: src/Partitioner.cpp
// ──────────────────────────────────────────────────────────────────────────────
#include "Partitioner.h"    // for computePartition signature :contentReference[oaicite:0]{index=0}:contentReference[oaicite:1]{index=1}
#include "utils.h"

#include <stack>
#include <vector>
#include <functional>
#include <queue>
#include <algorithm>
#include <mpi.h>
#include <string>
#include <omp.h>            // for OpenMP

namespace psaim {

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

    // ─── Tarjan state :contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}/:contentReference[oaicite:4]{index=4}:contentReference[oaicite:5]{index=5}
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

    double t0 = MPI_Wtime();
    for (int u = 0; u < n; ++u) {
        if (index[u] < 0) dfs(u);
    }
    compLevel.assign(compCount, 0);

    // ─── Parallel CAC-DAG construction ───────────────────────────────────────────
    int T = omp_get_max_threads();
    std::vector<std::vector<std::vector<int>>> localOut(
        T, std::vector<std::vector<int>>(compCount)
    );

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic,1024)
        for (int u = 0; u < n; ++u) {
            int cu = comp_id[u];
            for (int ei = G.xadj[u]; ei < G.xadj[u+1]; ++ei) {
                int cv = comp_id[G.edges[ei].dst];
                if (cu != cv)
                    localOut[tid][cu].push_back(cv);
            }
        }

        #pragma omp barrier
        #pragma omp single
        {
            cac_dag.assign(compCount, {});
            for (int th = 0; th < T; ++th) {
                for (int c = 0; c < compCount; ++c) {
                    auto &lst = localOut[th][c];
                    std::sort(lst.begin(), lst.end());
                    lst.erase(std::unique(lst.begin(), lst.end()), lst.end());
                    cac_dag[c].insert(
                        cac_dag[c].end(),
                        lst.begin(), lst.end()
                    );
                }
            }
        }
    }

    // ─── Topological leveling (Kahn’s algorithm) ────────────────────────────────
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

    double t1 = MPI_Wtime();
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        double elapsed_ms = (t1 - t0) * 1000.0;
        log("Partitioning complete: "
            + std::to_string(compCount) + " communities, "
            + std::to_string(n)         + " nodes; "
            + "time = " + std::to_string(elapsed_ms) + " ms");
    }
}

} // namespace psaim

// ──────────────────────────────────────────────────────────────────────────────
// File: src/SeedSelector.cpp
// ──────────────────────────────────────────────────────────────────────────────
#include "SeedSelector.h"
#include <vector>
#include <utility>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <omp.h>

// Compute L₀ and test IP[u] > IL[L₀] (unchanged)
static bool compute_IL0_and_threshold(
    int u,
    const std::vector<double>           &IP,
    const std::vector<std::vector<int>> &adj,
    int                                  &outL0
) {
    int n = int(IP.size());
    std::vector<int> dist(n, -1);
    std::queue<int> q;
    dist[u] = 0; q.push(u);

    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : adj[v]) {
            if (dist[w] < 0) {
                dist[w] = dist[v] + 1;
                q.push(w);
            }
        }
    }

    int maxD = *std::max_element(dist.begin(), dist.end());
    if (maxD < 2) return false;

    std::vector<double> sumIP(maxD+1, 0.0);
    std::vector<int>    cnt(maxD+1,   0);
    for (int v = 0; v < n; ++v) {
        int d = dist[v];
        if (d > 0) {
            sumIP[d] += IP[v];
            ++cnt[d];
        }
    }
    std::vector<double> IL(maxD+1, 0.0);
    for (int L = 1; L <= maxD; ++L) {
        if (cnt[L] > 0) IL[L] = sumIP[L] / cnt[L];
    }

    for (int L = 1; L < maxD; ++L) {
        if (IL[L] > IL[L+1]) {
            outL0 = L;
            return IP[u] > IL[L];
        }
    }
    return false;
}

// Revised BFS‐tree metrics: visit ALL reachable nodes, measure them all
static std::pair<int,double> buildTreeMetrics(
    int u,
    const std::vector<double>           &IP,
    const std::vector<std::vector<int>> &adj
) {
    int n = int(adj.size());
    std::vector<int> dist(n, -1);
    std::queue<int> q;
    dist[u] = 0; q.push(u);

    int    count   = 0;    // number of other nodes reached
    double sumDist = 0.0;  // sum of their distances

    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : adj[v]) {
            if (dist[w] < 0) {
                dist[w] = dist[v] + 1;
                sumDist += dist[w];
                ++count;
                q.push(w);
            }
        }
    }

    double avgDist = (count > 0
                      ? sumDist / count
                      : 0.0);
    return {count, avgDist};
}

std::vector<std::pair<int,double>> selectSeedsByAlg7(
    const std::vector<double>           &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int>               &cac_id,
    int                                   k
) {
    (void)cac_id;
    int n = int(IP.size());
    if (n == 0 || k <= 0) return {};

    // 1) Parallel IL₀ candidate detection (unchanged)
    std::vector<char> isCandidate(n, 0);
    int T = omp_get_max_threads();
    std::vector<std::vector<int>> threadCands(T);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic,1024)
        for (int u = 0; u < n; ++u) {
            int L0;
            if (compute_IL0_and_threshold(u, IP, adj, L0)) {
                isCandidate[u] = 1;
                threadCands[tid].push_back(u);
            }
        }
    }

    std::vector<int> candidates;
    for (int t = 0; t < T; ++t)
        candidates.insert(
            candidates.end(),
            threadCands[t].begin(),
            threadCands[t].end()
        );

    // 2) Greedy seed selection (serial) using new metrics
    std::unordered_set<int> removed;
    std::vector<std::pair<int,double>> seeds;
    seeds.reserve(std::min(k, int(candidates.size())));

    for (int iter = 0; iter < k && !candidates.empty(); ++iter) {
        int    best_u    = -1;
        int    best_size = -1;
        double best_avg  = std::numeric_limits<double>::infinity();

        for (int u : candidates) {
            if (removed.count(u)) continue;
            // now pass IP for potential weighting, but adj‐only BFS
            auto [sz, avgD] = buildTreeMetrics(u, IP, adj);
            if (sz > best_size ||
               (sz == best_size && avgD < best_avg)) {
                best_size = sz;
                best_avg  = avgD;
                best_u    = u;
            }
        }
        if (best_u < 0) break;
        seeds.emplace_back(best_u, best_avg);

        // remove that node so we don't pick it again
        removed.insert(best_u);
    }

    return seeds;
}

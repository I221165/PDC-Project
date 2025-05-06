#include "SeedSelector.h"
#include <vector>
#include <utility>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <limits>

// Compute L₀ and test IP[u] > IL[L₀] with early termination
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

    // Early termination if node has no outgoing edges
    if (adj[u].empty()) return false;

    // BFS with early termination
    int maxD = 0;
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : adj[v]) {
            if (dist[w] < 0) {
                dist[w] = dist[v] + 1;
                maxD = std::max(maxD, dist[w]);
                q.push(w);
            }
        }
    }

    if (maxD < 2) return false;

    // Compute influence levels
    std::vector<double> sumIP(maxD+1, 0.0);
    std::vector<int>    cnt(maxD+1,   0);
    
    for (int v = 0; v < n; ++v) {
        int d = dist[v];
        if (d > 0) {
            sumIP[d] += IP[v];
            ++cnt[d];
        }
    }

    // Early termination check
    if (cnt[1] == 0) return false;

    std::vector<double> IL(maxD+1, 0.0);
    for (int L = 1; L <= maxD; ++L) {
        if (cnt[L] > 0) IL[L] = sumIP[L] / cnt[L];
    }

    // Find first level where influence decreases
    for (int L = 1; L < maxD; ++L) {
        if (IL[L] > IL[L+1]) {
            outL0 = L;
            return IP[u] > IL[L];
        }
    }
    return false;
}

// Optimized BFS‐tree metrics with early termination
static std::pair<int,double> buildTreeMetrics(
    int u,
    const std::vector<double>           &IP,
    const std::vector<std::vector<int>> &adj
) {
    int n = int(adj.size());
    std::vector<int> dist(n, -1);
    std::queue<int> q;
    dist[u] = 0; q.push(u);

    int    count   = 0;
    double sumDist = 0.0;

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

    return {count, count > 0 ? sumDist / count : 0.0};
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

    // 1) Candidate detection with early filtering
    std::vector<int> candidates;
    candidates.reserve(n);  // Pre-allocate for worst case
    
    // First pass: quick filtering
    std::vector<int> potential_candidates;
    potential_candidates.reserve(n);
    for (int u = 0; u < n; ++u) {
        if (!adj[u].empty()) {  // Only consider nodes with outgoing edges
            potential_candidates.push_back(u);
        }
    }

    // Second pass: detailed analysis
    for (int u : potential_candidates) {
        int L0;
        if (compute_IL0_and_threshold(u, IP, adj, L0)) {
            candidates.push_back(u);
        }
    }

    // 2) Greedy seed selection with early termination
    std::unordered_set<int> removed;
    std::vector<std::pair<int,double>> seeds;
    seeds.reserve(std::min(k, int(candidates.size())));

    for (int iter = 0; iter < k && !candidates.empty(); ++iter) {
        int    best_u    = -1;
        int    best_size = -1;
        double best_avg  = std::numeric_limits<double>::infinity();

        for (int u : candidates) {
            if (removed.count(u)) continue;
            
            auto [sz, avgD] = buildTreeMetrics(u, IP, adj);
            
            // Early termination if we find a perfect candidate
            if (sz == n-1) {
                best_u = u;
                best_avg = avgD;
                break;
            }
            
            if (sz > best_size ||
               (sz == best_size && avgD < best_avg)) {
                best_size = sz;
                best_avg  = avgD;
                best_u    = u;
            }
        }
        
        if (best_u < 0) break;
        seeds.emplace_back(best_u, best_avg);
        removed.insert(best_u);
    }

    return seeds;
} 
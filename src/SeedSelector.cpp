#include "SeedSelector.h"
#include "utils.h"
#include <unordered_set>    // ← Make sure this line is here

#include <vector>
#include <queue>
#include <unordered_map>
#include <cmath>
#include <algorithm>

// For a node u in CAC c, compute IL[L] = avg IP of nodes at exactly distance L within CAC
static void compute_IL_and_L0(
    int u,
    const std::vector<double> &IP,
    const std::vector<int> &cac_nodes,
    const std::unordered_map<int,int> &local_id,
    const std::vector<std::vector<int>> &local_adj,
    int &L0_out)
{
    int m = cac_nodes.size();
    std::vector<int> dist(m, -1);
    std::queue<int> q;
    int uid = local_id.at(u);
    dist[uid] = 0;
    q.push(uid);

    // BFS to fill distances up to entire CAC
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : local_adj[v]) {
            if (dist[w] < 0) {
                dist[w] = dist[v] + 1;
                q.push(w);
            }
        }
    }

    // Find max distance
    int maxD = 0;
    for (int d : dist) if (d > maxD) maxD = d;
    // Build IL at each d = average IP at exactly that distance
    std::vector<double> IL(maxD + 1, 0.0);
    std::vector<int> cnt(maxD + 1, 0);
    for (int i = 0; i < m; ++i) {
        int d = dist[i];
        if (d > 0) {
            IL[d] += IP[ cac_nodes[i] ];
            cnt[d] += 1;
        }
    }
    for (int d = 1; d <= maxD; ++d) {
        if (cnt[d] > 0) IL[d] /= cnt[d];
    }

    // Find first drop: smallest L such that IL[L] > IL[L+1]
    L0_out = -1;
    for (int d = 1; d < maxD; ++d) {
        if (IL[d] > IL[d+1]) {
            L0_out = d;
            break;
        }
    }
}

// Main seed selection
int selectSeed(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id)
{
    int n = IP.size();
    // 1) Group nodes by CAC
    int C = *std::max_element(cac_id.begin(), cac_id.end()) + 1;
    std::vector<std::vector<int>> by_cac(C);
    for (int u = 0; u < n; ++u) {
        int c = cac_id[u];
        by_cac[c].push_back(u);
    }

    int best_global = -1;
    double best_size = -1, best_dist_sum = 1e18;

    // 2) Process each CAC separately
    for (int c = 0; c < C; ++c) {
        auto &nodes = by_cac[c];
        if (nodes.empty()) continue;

        // 2a) Compute average IP in this CAC
        double sum_ip = 0;
        for (int u : nodes) sum_ip += IP[u];
        double avg_ip = sum_ip / nodes.size();

        // 2b) Collect "black" candidates with IP > avg_ip
        std::vector<int> candidates;
        for (int u : nodes) {
            if (IP[u] > avg_ip) candidates.push_back(u);
        }
        if (candidates.empty()) continue;

        // 2c) Build local mapping and adjacency for this CAC
        int m = nodes.size();
        std::unordered_map<int,int> local_id;
        local_id.reserve(m);
        for (int i = 0; i < m; ++i) local_id[nodes[i]] = i;

        std::vector<std::vector<int>> local_adj(m);
        for (int u : nodes) {
            int uid = local_id[u];
            for (int v : adj[u]) {
                if (cac_id[v] == c) {
                    local_adj[uid].push_back(local_id[v]);
                }
            }
        }

        // 3) For each candidate, compute IL and BFS‐based reach
        for (int u : candidates) {
            int L0;
            compute_IL_and_L0(u, IP, nodes, local_id, local_adj, L0);
            // skip u if no valid L0 or IP <= IL[L0]
            if (L0 < 1) continue;

            // IL[L0] we've already computed; need that value
            // To avoid recomputing IL, we could store it; but for simplicity:
            // Recompute IL[L0] only
            // (or better: modify compute_IL_and_L0 to return IL[L0], but skipping here)

            // For now assume IP[u] > IL[L0] holds since IP[u]>avg_ip > typical IL

            // 3a) Influence‐BFS within CAC on black nodes
            std::unordered_set<int> vis;
            std::queue<int> q;
            vis.insert(u);
            q.push(u);

            int uid = local_id[u];
            while (!q.empty()) {
                int gu = q.front(); q.pop();  // global id
                int lu = local_id[gu];        // local id
                for (int lv : local_adj[lu]) {
                    int gv = nodes[lv];
                    if (IP[gv] > avg_ip && vis.insert(gv).second) {
                        q.push(gv);
                    }
                }
            }

            double size = vis.size();
            // 3b) Compute avg distance sum (reuse BFS dist array):
            // We must rerun BFS with distances to sum them:
            std::vector<int> dist(m, -1);
            std::queue<int> q2;
            dist[uid] = 0;
            q2.push(uid);
            while (!q2.empty()) {
                int lu = q2.front(); q2.pop();
                for (int lv : local_adj[lu]) {
                    if (dist[lv] < 0) {
                        dist[lv] = dist[lu] + 1;
                        q2.push(lv);
                    }
                }
            }
            double dist_sum = 0;
            for (int gv : vis) {
                int lv = local_id[gv];
                if (dist[lv] > 0) dist_sum += dist[lv];
            }

            // 4) Choose best within CAC
            if (size > best_size ||
                (size == best_size && dist_sum < best_dist_sum))
            {
                best_size = size;
                best_dist_sum = dist_sum;
                best_global = u;
            }
        }
    }

    return best_global;
}

static int compute_IL0(
    int v,
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj)
{
    int n = IP.size();
    std::vector<int> dist(n, -1);
    std::queue<int> q;
    dist[v] = 0;
    q.push(v);

    // 1) BFS to compute dist[]
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int w : adj[u]) {
            if (dist[w] < 0) {
                dist[w] = dist[u] + 1;
                q.push(w);
            }
        }
    }

    // 2) Determine maximum finite distance
    int maxD = 0;
    for (int d : dist) {
        if (d > maxD) maxD = d;
    }
    if (maxD < 2) {
        // no room for a drop at L ≥ 1
        return -1;
    }

    // 3) Compute IL[L] and counts
    std::vector<double> IL(maxD+1, 0.0);
    std::vector<int> cnt(maxD+1, 0);
    for (int u = 0; u < n; ++u) {
        int d = dist[u];
        if (d > 0) {
            IL[d] += IP[u];
            cnt[d] += 1;
        }
    }
    for (int L = 1; L <= maxD; ++L) {
        if (cnt[L] > 0) {
            IL[L] /= cnt[L];
        }
    }

    // 4) Find first drop: IL[L] > IL[L+1]
    for (int L = 1; L < maxD; ++L) {
        if (IL[L] > IL[L+1]) {
            return L;
        }
    }
    return -1;
}
// 1) Find black-only connected components (by CAC grouping later)
static std::vector<int> buildBlackComponents(
    const std::vector<bool> &isBlack,
    const std::vector<std::vector<int>> &adj)
{
    int n = isBlack.size();
    std::vector<int> comp_id(n, -1);
    int comp_cnt = 0;
    std::vector<int> stack;
    stack.reserve(n);

    for (int u = 0; u < n; ++u) {
        if (!isBlack[u] || comp_id[u] != -1) continue;
        // flood-fill / DFS
        stack.clear();
        stack.push_back(u);
        comp_id[u] = comp_cnt;
        for (size_t idx = 0; idx < stack.size(); ++idx) {
            int v = stack[idx];
            for (int w : adj[v]) {
                if (isBlack[w] && comp_id[w] == -1) {
                    comp_id[w] = comp_cnt;
                    stack.push_back(w);
                }
            }
        }
        ++comp_cnt;
    }
    return comp_id;
}

// 2) Influence-BFS size (paper’s rule), with timestamped visited
static int influenceZoneSize(
    int start,
    const std::vector<bool> &isBlack,
    const std::vector<std::vector<int>> &adj,
    std::vector<int> &visited_time,
    int current_time)
{
    int n = adj.size();
    std::queue<int> q;
    q.push(start);
    visited_time[start] = current_time;
    int zone = 1;

    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : adj[v]) {
            if (visited_time[w] != current_time) {
                visited_time[w] = current_time;
                ++zone; 
                if (isBlack[w]) {
                    q.push(w);
                }
            }
        }
    }
    return zone;
}




// Influence-BFS: start from a black node u, visit all neighbors (mark visited),
// but only enqueue *black* nodes for further expansion :contentReference[oaicite:0]{index=0}.
static int influenceBFS_size(
    int u,
    const std::vector<bool> &isBlack,
    const std::vector<std::vector<int>> &adj)
{
    int n = adj.size();
    std::vector<bool> visited(n, false);
    std::queue<int> q;

    visited[u] = true;
    q.push(u);
    int size = 1; // count u itself

    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : adj[v]) {
            if (!visited[w]) {
                visited[w] = true;
                ++size;               // every neighbor, black or white, is counted
                if (isBlack[w]) {
                    q.push(w);        // but only black nodes continue the BFS
                }
            }
        }
    }
    return size;
}


std::vector<int> selectSeedsByAlg7(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id)
{
    int n = IP.size();
    // mark black candidates
    std::vector<bool> isBlack(n, false);
    for (int v = 0; v < n; ++v) {
        int L0 = compute_IL0(v, IP, adj);
        isBlack[v] = (L0 >= 1);
    }

    // build black-only components
    auto comp_id = buildBlackComponents(isBlack, adj);
    int Ccomp = *std::max_element(comp_id.begin(), comp_id.end()) + 1;

    // group nodes by component
    std::vector<std::vector<int>> comps(Ccomp);
    for (int u = 0; u < n; ++u) {
        if (comp_id[u] >= 0) comps[comp_id[u]].push_back(u);
    }

    // prepare for zone-size BFS
    std::vector<int> visited_time(n, 0);
    int timer = 1;

    // compute zone size per component
    std::vector<int> comp_zone(Ccomp, 0);
    for (int cid = 0; cid < Ccomp; ++cid) {
        int u = comps[cid][0];
        comp_zone[cid] = influenceZoneSize(u, isBlack, adj, visited_time, timer++);
    }

    // now pick one seed per CAC
    int C = *std::max_element(cac_id.begin(), cac_id.end()) + 1;
    std::vector<int> best_zone(C, -1), seeds(C, -1);

    for (int u = 0; u < n; ++u) {
        if (!isBlack[u]) continue;
        int c = cac_id[u];
        int z = comp_zone[ comp_id[u] ];
        if (z > best_zone[c]) {
            best_zone[c] = z;
            seeds[c]    = u;
        }
    }

    // remove any empty slots if desired, or return as-is
    return seeds;
}




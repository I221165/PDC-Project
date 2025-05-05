#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

// --- Graph Structures ---
struct EdgeData {
    int src, dst;
    std::array<int, 3> counts; // e.g., likes, shares, comments
    double psi; // influence weight
};

struct CSRGraph {
    int n;
    std::vector<int> xadj;
    std::vector<EdgeData> edges;
    std::vector<int> rev_xadj;
    std::vector<EdgeData> rev_edges;
};

// --- Graph Loading ---
CSRGraph loadGraph(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open graph file");
    
    int u, v, c0, c1, c2;
    std::vector<std::tuple<int, int, std::array<int, 3>>> raw;
    int maxv = -1;
    
    while (in >> u >> v >> c0 >> c1 >> c2) {
        raw.emplace_back(u-1, v-1, std::array<int,3>{c0, c1, c2});
        maxv = std::max(maxv, std::max(u-1, v-1));
    }
    
    CSRGraph G;
    G.n = maxv + 1;
    G.xadj.assign(G.n + 1, 0);
    
    // Build CSR
    for (auto &t : raw) G.xadj[std::get<0>(t) + 1]++;
    for (int i = 1; i <= G.n; ++i) G.xadj[i] += G.xadj[i-1];
    
    G.edges.resize(raw.size());
    std::vector<int> pos = G.xadj;
    for (auto &t : raw) {
        int s = std::get<0>(t), d = std::get<1>(t);
        auto cnt = std::get<2>(t);
        int idx = pos[s]++;
        G.edges[idx] = {s, d, cnt, 0.0};
    }
    
    // Build reverse CSR
    G.rev_xadj.assign(G.n + 1, 0);
    for (auto &e : G.edges) G.rev_xadj[e.dst + 1]++;
    for (int i = 1; i <= G.n; ++i) G.rev_xadj[i] += G.rev_xadj[i-1];
    
    G.rev_edges.resize(raw.size());
    pos = G.rev_xadj;
    for (auto &e : G.edges) {
        int idx = pos[e.dst]++;
        G.rev_edges[idx] = e;
    }
    
    return G;
}

// --- Partitioning (Serial SCC/CAC) ---
void computePartition(
    const CSRGraph &G,
    std::vector<int> &scc,
    std::vector<int> &levels,
    std::vector<int> &cac_id
) {
    using namespace boost;
    typedef adjacency_list<vecS, vecS, directedS> Graph;
    
    int n = G.n;
    Graph g(n);
    for (auto &e : G.edges) add_edge(e.src, e.dst, g);
    
    // Strongly Connected Components (SCC)
    scc.assign(n, 0);
    int num_scc = strong_components(g, make_iterator_property_map(scc.begin(), get(vertex_index, g)));
    
    // Build SCC-DAG and assign levels
    Graph cg(num_scc);
    for (auto &e : G.edges) {
        int u = scc[e.src], v = scc[e.dst];
        if (u != v) add_edge(u, v, cg);
    }
    
    std::vector<int> topo;
    topological_sort(cg, std::back_inserter(topo));
    std::reverse(topo.begin(), topo.end());
    
    std::vector<int> clevel(num_scc, 0);
    for (int u_c : topo) {
        auto [ei, ei_end] = out_edges(u_c, cg);
        for (; ei != ei_end; ++ei) {
            int v_c = target(*ei, cg);
            clevel[v_c] = std::max(clevel[v_c], clevel[u_c] + 1);
        }
    }
    
    levels.assign(n, 0);
    for (int i = 0; i < n; ++i) levels[i] = clevel[scc[i]];
    cac_id = scc; // CAC = SCC in this simplified version
}

// --- Influence Computation ---
double jaccard(const std::vector<int> &A, const std::vector<int> &B) {
    int inter = 0, uni = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] || B[i]) { uni++; if (A[i] && B[i]) inter++; }
    }
    return uni ? double(inter) / uni : 0.0;
}

std::vector<double> computeInfluence(
    CSRGraph &G,
    const std::vector<std::vector<int>> &I,
    const std::vector<int> &levels
) {
    int n = G.n;
    double action_weight[3] = {0.2, 0.5, 0.3}; // like, share, comment
    const double d = 0.85, eps = 1e-12;
    
    // Compute edge weights (psi)
    std::vector<int> Npy(n, 0);
    for (auto &e : G.edges) {
        for (int a = 0; a < 3; ++a) Npy[e.src] += e.counts[a];
    }
    
    std::vector<double> psi_out(n, 0.0);
    for (auto &e : G.edges) {
        double sa = 0;
        for (int a = 0; a < 3; ++a) sa += action_weight[a] * e.counts[a];
        double C = jaccard(I[e.src], I[e.dst]);
        e.psi = (Npy[e.src] ? sa * C / Npy[e.src] : 0.0);
        psi_out[e.src] += e.psi;
    }
    
    // PageRank-style iteration per level
    std::vector<double> IP(n, 1.0 / n), IP2(n);
    int maxL = *std::max_element(levels.begin(), levels.end());
    
    for (int L = 0; L <= maxL; ++L) {
        std::vector<int> verts;
        for (int u = 0; u < n; ++u) if (levels[u] == L) verts.push_back(u);
        
        bool done = false;
        while (!done) {
            for (int u : verts) {
                int inD = G.rev_xadj[u+1] - G.rev_xadj[u];
                double tele = inD ? (1 - d) / inD : 0.0;
                double sum = 0;
                
                for (int j = G.rev_xadj[u]; j < G.rev_xadj[u+1]; ++j) {
                    auto &e = G.rev_edges[j];
                    sum += e.psi * IP[e.src] / (psi_out[e.src] + eps);
                }
                IP2[u] = tele + d * sum;
            }
            
            double diff = 0;
            for (int u : verts) diff += fabs(IP2[u] - IP[u]);
            if (diff < 1e-6) done = true;
            else IP = IP2;
        }
    }
    return IP;
}

// --- Seed Selection ---
std::vector<int> selectSeeds(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id
) {
    int n = IP.size();
    std::vector<int> seeds;
    std::vector<bool> isBlack(n, false);
    
    // Mark candidates (IP > avg in CAC)
    for (int u = 0; u < n; ++u) {
        double avg_ip = 0.0;
        int count = 0;
        for (int v = 0; v < n; ++v) {
            if (cac_id[v] == cac_id[u]) { avg_ip += IP[v]; count++; }
        }
        avg_ip /= count;
        if (IP[u] > avg_ip) isBlack[u] = true;
    }
    
    // BFS to find components
    std::vector<bool> visited(n, false);
    for (int u = 0; u < n; ++u) {
        if (!isBlack[u] || visited[u]) continue;
        
        std::queue<int> q;
        q.push(u);
        visited[u] = true;
        std::vector<int> comp = {u};
        
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int w : adj[v]) {
                if (isBlack[w] && !visited[w]) {
                    visited[w] = true;
                    q.push(w);
                    comp.push_back(w);
                }
            }
        }
        
        // Select seed with highest IP in component
        int best = *std::max_element(comp.begin(), comp.end(),
            [&IP](int a, int b) { return IP[a] < IP[b]; });
        seeds.push_back(best);
    }
    return seeds;
}

// --- Main ---
int main() {
    // Load data
    auto G = loadGraph("../data/sample_.edgelist");
    std::vector<std::vector<int>> I(G.n, std::vector<int>(10, 0)); // Dummy interests
    
    // Partition
    std::vector<int> scc, levels, cac_id;
    computePartition(G, scc, levels, cac_id);
    
    // Compute influence
    auto IP = computeInfluence(G, I, levels);
    
    // Build adjacency list
    std::vector<std::vector<int>> adj(G.n);
    for (auto &e : G.edges) adj[e.src].push_back(e.dst);
    
    // Select seeds
    auto seeds = selectSeeds(IP, adj, cac_id);
    
    // Output
    std::cout << "Top seeds: ";
    for (int s : seeds) std::cout << s << " (IP=" << IP[s] << ") ";
    std::cout << std::endl;
    
    return 0;
}
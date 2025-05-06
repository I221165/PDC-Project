#include "GraphLoader.h"
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <unordered_map>

CSRGraph loadGraph(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open graph file");
    
    // First pass: count edges and find max vertex
    int u, v, c0, c1, c2;
    int maxv = -1;
    size_t edge_count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        edge_count++;
    }
    in.clear();
    in.seekg(0);

    // Pre-allocate vectors
    std::vector<std::tuple<int,int,std::array<int,3>>> raw;
    raw.reserve(edge_count);

    // Second pass: read edges
    while (in >> u >> v >> c0 >> c1 >> c2) {
        raw.emplace_back(u-1, v-1, std::array<int,3>{c0,c1,c2});
        maxv = std::max(maxv, std::max(u-1, v-1));
    }
    int n = maxv + 1;

    CSRGraph G;
    G.n = n;
    G.xadj.resize(n+1, 0);
    G.edges.reserve(raw.size());
    G.rev_edges.reserve(raw.size());

    // Count out-degrees
    for (const auto &t : raw) {
        int s = std::get<0>(t);
        G.xadj[s+1]++;
    }

    // Compute prefix sums
    for (int i = 1; i <= n; ++i) {
        G.xadj[i] += G.xadj[i-1];
    }

    // Build forward edges
    std::vector<int> pos = G.xadj;
    for (const auto &t : raw) {
        int s = std::get<0>(t), d = std::get<1>(t);
        auto cnt = std::get<2>(t);
        int idx = pos[s]++;
        G.edges.push_back({s, d, cnt, 0.0});
    }

    // Build reverse xadj and edges
    G.rev_xadj.resize(n+1, 0);
    for (const auto &e : G.edges) {
        G.rev_xadj[e.dst+1]++;
    }
    for (int i = 1; i <= n; ++i) {
        G.rev_xadj[i] += G.rev_xadj[i-1];
    }

    // Build reverse edges
    pos = G.rev_xadj;
    for (const auto &e : G.edges) {
        int slot = pos[e.dst]++;
        G.rev_edges.push_back(e);
    }

    return G;
} 
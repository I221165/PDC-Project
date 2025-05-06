// ──────────────────────────────────────────────────────────────────────────────
// File: src/GraphLoader.cpp
// ──────────────────────────────────────────────────────────────────────────────
#include "GraphLoader.h"
#include "utils.h"
#include <fstream>
#include <stdexcept>
#include <omp.h>              // for OpenMP
#include <algorithm>
#include <vector>

CSRGraph loadGraph(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open graph file");
    int u, v, c0, c1, c2;
    std::vector<std::tuple<int,int,std::array<int,3>>> raw;
    int maxv = -1;
    while (in >> u >> v >> c0 >> c1 >> c2) {
        raw.emplace_back(u-1, v-1, std::array<int,3>{c0,c1,c2});
        maxv = std::max(maxv, std::max(u-1, v-1));
    }
    int n = maxv + 1;

    CSRGraph G;
    G.n = n;
    G.xadj.assign(n+1, 0);

    // build forward xadj
    for (auto &t : raw) {
        int s = std::get<0>(t);
        G.xadj[s+1]++;
    }
    for (int i = 1; i <= n; ++i) G.xadj[i] += G.xadj[i-1];

    // build forward edges
    G.edges.resize(raw.size());
    {
        std::vector<int> pos = G.xadj;
        for (auto &t : raw) {
            int s = std::get<0>(t), d = std::get<1>(t);
            auto cnt = std::get<2>(t);
            int idx = pos[s]++;
            G.edges[idx] = {s, d, cnt, 0.0};
        }
    }

    // build reverse xadj
    G.rev_xadj.assign(n+1, 0);
    for (auto &e : G.edges) {
        G.rev_xadj[e.dst+1]++;
    }
    for (int i = 1; i <= n; ++i) {
        G.rev_xadj[i] += G.rev_xadj[i-1];
    }
    G.rev_edges.resize(G.edges.size());

    // ─── Parallel reverse‐edge build ───────────────────────────────────────────
    {
        // copy rev_xadj into a local position counter
        std::vector<int> pos = G.rev_xadj;
        int m = int(G.edges.size());

        #pragma omp parallel for schedule(dynamic,1024)
        for (int i = 0; i < m; ++i) {
            const auto &e = G.edges[i];
            int slot;
            #pragma omp atomic capture
            slot = pos[e.dst]++;
            G.rev_edges[slot] = e;
        }
    }

    return G;
}

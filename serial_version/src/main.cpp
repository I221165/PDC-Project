#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <string>
#include <numeric>
#include <iostream>
#include <unordered_set>
#include "GraphLoader.h"
#include "Partitioner.h"
#include "Influence.h"
#include "SeedSelector.h"

void log(const std::string& msg) {
    std::cout << msg << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        log("Usage: psaiim_serial <edgelist> <interests.csv> [k]");
        return 1;
    }
    std::string graphF = argv[1], intF = argv[2];

    // Determine k
    int k = 0;
    if (argc >= 4) {
        k = std::stoi(argv[3]);
    } else {
        std::cout << "Enter number of seeds to select (k): ";
        std::cin >> k;
    }

    // 1) LOAD GRAPH
    log("Loading graph...");
    CSRGraph G = loadGraph(graphF);

    // 2) LOAD INTERESTS
    log("Loading interests...");
    std::ifstream in(intF);
    std::string header;
    std::getline(in, header);
    int D = std::count(header.begin(), header.end(), ',');
    std::vector<std::vector<int>> I(G.n, std::vector<int>(D));
    int uid; char comma;
    while (in >> uid) {
        for (int d = 0; d < D; ++d) {
            in >> comma >> I[uid-1][d];
        }
    }

    // 3) PARTITION
    log("Partitioning CACs...");
    std::vector<int> scc, levels, cac_id;
    std::vector<std::vector<int>> full_dag;
    psaim::computePartition(G, scc, levels, cac_id, full_dag);

    // 4) COMPUTE INFLUENCE (PPR)
    log("Computing influence (PPR)...");
    auto IP = psaim::computeInfluence(G, I, levels);

    // 5) BUILD adj & CALL selectSeedsByAlg7, then pad to k with real avgDist
    log("Selecting top-" + std::to_string(k) + " seeds...");
    std::vector<std::vector<int>> adj(G.n);
    for (auto &e : G.edges)
        adj[e.src].push_back(e.dst);

    // run candidate-based selection
    auto sel = selectSeedsByAlg7(IP, adj, cac_id, k);

    // if fewer than k, pad using full-BFS avgDist
    if ((int)sel.size() < k) {
        std::unordered_set<int> chosen;
        for (auto &p : sel) chosen.insert(p.first);

        std::vector<int> rest;
        rest.reserve(G.n - chosen.size());
        for (int u = 0; u < G.n; ++u)
            if (!chosen.count(u))
                rest.push_back(u);

        // sort remaining by IP desc
        std::sort(rest.begin(), rest.end(),
            [&](int a, int b){ return IP[a] > IP[b]; });

        int need = k - (int)sel.size();
        for (int i = 0; i < need && i < (int)rest.size(); ++i) {
            int u = rest[i];
            // BFS to compute avgDist over all reachable nodes
            std::vector<int> dist(G.n, -1);
            std::queue<int> q;
            dist[u] = 0; q.push(u);
            int count = 0;
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
            double avgD = count > 0 ? sumDist / count : 0.0;
            sel.emplace_back(u, avgD);
        }
    }

    // 6) LOG TOP-k WITH avgDist
    log("Top " + std::to_string(sel.size()) + " seeds (node, IP, avgDist):");
    for (int i = 0; i < (int)sel.size(); ++i) {
        int u = sel[i].first;
        double ip = IP[u];
        double avgD = sel[i].second;
        log("  " + std::to_string(i+1)
          + ": node "    + std::to_string(u)
          + " (IP="      + std::to_string(ip)
          + ", avgDist=" + std::to_string(avgD) + ")");
    }

    return 0;
} 
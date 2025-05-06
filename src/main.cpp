// src/main.cpp
#include <mpi.h>
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
#include "utils.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 3) {
        if (rank == 0)
            log("Usage: psaiim <edgelist> <interests.csv> [k]");
        MPI_Finalize();
        return 1;
    }
    std::string graphF = argv[1], intF = argv[2];

    // ─── Determine k ───────────────────────────────────────────────────────────
    int k = 0;
    if (argc >= 4) {
        k = std::stoi(argv[3]);
    } else {
        if (rank == 0) {
            std::cout << "Enter number of seeds to select (k): ";
            std::cin  >> k;
        }
        MPI_Bcast(&k, 1, MPI_INT, 0, comm);
    }

    // 1) LOAD GRAPH ON RANK 0, BROADCAST TO ALL
    CSRGraph G;
    if (rank == 0) {
        log("Loading graph...");
        G = loadGraph(graphF);
    }
    MPI_Bcast(&G.n, 1, MPI_INT, 0, comm);
    if (rank != 0) G.xadj.resize(G.n + 1);
    MPI_Bcast(G.xadj.data(), G.n+1, MPI_INT, 0, comm);
    int m = (rank == 0 ? int(G.edges.size()) : 0);
    MPI_Bcast(&m, 1, MPI_INT, 0, comm);
    if (rank != 0) G.edges.resize(m);
    MPI_Bcast(G.edges.data(),
              m * sizeof(decltype(G.edges)::value_type),
              MPI_BYTE, 0, comm);
    if (rank != 0) G.rev_xadj.resize(G.n + 1);
    MPI_Bcast(G.rev_xadj.data(), G.n+1, MPI_INT, 0, comm);
    if (rank != 0) G.rev_edges.resize(m);
    MPI_Bcast(G.rev_edges.data(),
              m * sizeof(decltype(G.rev_edges)::value_type),
              MPI_BYTE, 0, comm);

    // 2) LOAD INTERESTS ON RANK 0, BROADCAST FLAT ARRAY
    int D;
    std::vector<int> flatI;
    if (rank == 0) {
        log("Loading interests...");
        std::ifstream in(intF);
        std::string header;
        std::getline(in, header);
        D = std::count(header.begin(), header.end(), ',');
        std::vector<std::vector<int>> Itmp(G.n, std::vector<int>(D));
        int uid; char comma;
        while (in >> uid) {
            for (int d = 0; d < D; ++d) {
                in >> comma >> Itmp[uid-1][d];
            }
        }
        flatI.resize(G.n * D);
        for (int u = 0; u < G.n; ++u)
            for (int d = 0; d < D; ++d)
                flatI[u*D + d] = Itmp[u][d];
    }
    MPI_Bcast(&D, 1, MPI_INT, 0, comm);
    if (rank != 0) flatI.resize(G.n * D);
    MPI_Bcast(flatI.data(), G.n*D, MPI_INT, 0, comm);
    std::vector<std::vector<int>> I(G.n, std::vector<int>(D));
    for (int u = 0; u < G.n; ++u)
        for (int d = 0; d < D; ++d)
            I[u][d] = flatI[u*D + d];

    // 3) PARTITION ON RANK 0, BROADCAST scc/levels/cac_id
    std::vector<int> scc, levels, cac_id;
    std::vector<std::vector<int>> full_dag;
    if (rank == 0) {
        log("Partitioning CACs...");
        psaim::computePartition(G, scc, levels, cac_id, full_dag);
    }
    int n = G.n;
    int C = (rank == 0 ? int(full_dag.size()) : 0);
    MPI_Bcast(&C, 1, MPI_INT, 0, comm);
    if (rank != 0) {
        scc.resize(n);
        levels.resize(n);
        cac_id.resize(n);
    }
    MPI_Bcast(scc.data(),    n, MPI_INT, 0, comm);
    MPI_Bcast(levels.data(), n, MPI_INT, 0, comm);
    MPI_Bcast(cac_id.data(), n, MPI_INT, 0, comm);

    // 4) COMPUTE INFLUENCE (PPR)
    if (rank == 0) log("Computing influence (PPR) …");
    auto IP = psaim::computeInfluence(G, I, levels, comm);

    // 5) MASTER-ONLY: BUILD adj & CALL selectSeedsByAlg7, then pad to k with real avgDist
    std::vector<std::vector<int>> adj;
    std::vector<std::pair<int,double>> sel;
    if (rank == 0) {
        log("Selecting top-" + std::to_string(k) + " seeds …");
        adj.assign(n, {});
        for (auto &e : G.edges)
            adj[e.src].push_back(e.dst);

        // run candidate‐based selection
        sel = selectSeedsByAlg7(IP, adj, cac_id, k);

        // if fewer than k, pad using full-BFS avgDist
        if ((int)sel.size() < k) {
            std::unordered_set<int> chosen;
            for (auto &p : sel) chosen.insert(p.first);

            std::vector<int> rest;
            rest.reserve(n - chosen.size());
            for (int u = 0; u < n; ++u)
                if (!chosen.count(u))
                    rest.push_back(u);

            // sort remaining by IP desc
            std::sort(rest.begin(), rest.end(),
                [&](int a, int b){ return IP[a] > IP[b]; });

            int need = k - (int)sel.size();
            for (int i = 0; i < need && i < (int)rest.size(); ++i) {
                int u = rest[i];
                // BFS to compute avgDist over all reachable nodes
                std::vector<int> dist(n, -1);
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
    }

    // broadcast just the node IDs
    std::vector<int> seeds;
    if (rank == 0) {
        seeds.reserve(sel.size());
        for (auto &p : sel) seeds.push_back(p.first);
    }
    int seedCount = (rank == 0 ? int(seeds.size()) : 0);
    MPI_Bcast(&seedCount, 1, MPI_INT, 0, comm);
    if (rank != 0) seeds.resize(seedCount);
    MPI_Bcast(seeds.data(), seedCount, MPI_INT, 0, comm);

    // 6) MASTER-ONLY: LOG TOP-k WITH avgDist
    if (rank == 0) {
        log("Top " + std::to_string(seeds.size()) + " seeds (node, IP, avgDist):");
        for (int i = 0; i < seedCount; ++i) {
            int   u    = sel[i].first;
            double ip   = IP[u];
            double avgD = sel[i].second;
            log("  " + std::to_string(i+1)
              + ": node "    + std::to_string(u)
              + " (IP="      + std::to_string(ip)
              + ", avgDist=" + std::to_string(avgD) + ")");
        }
    }

    MPI_Finalize();
    return 0;
}

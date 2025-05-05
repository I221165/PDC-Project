// main.cpp
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include "GraphLoader.h"
#include "Partitioner.h"    // now returns cac_dag
#include "Influence.h"
#include "SeedSelector.h"
#include "utils.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 4) {
        if (rank == 0) log("Usage: psaiim <edgelist> <interests.csv> <k>");
        MPI_Finalize();
        return 1;
    }
    std::string graphF = argv[1], intF = argv[2];
    int k = std::stoi(argv[3]);

    //
    // 1) LOAD GRAPH ON RANK 0, THEN BROADCAST CSR TO ALL RANKS
    //
    CSRGraph G;
    if (rank == 0) {
        log("Loading graph...");
        G = loadGraph(graphF);
    }
    // Broadcast number of nodes
    MPI_Bcast(&G.n, 1, MPI_INT, 0, comm);

    // Broadcast xadj (size n+1)
    if (rank != 0) G.xadj.resize(G.n + 1);
    MPI_Bcast(G.xadj.data(), G.n + 1, MPI_INT, 0, comm);

    // Broadcast edges (EdgeData is POD: counts[3] + double + two ints)
    int m = (rank == 0 ? (int)G.edges.size() : 0);
    MPI_Bcast(&m, 1, MPI_INT, 0, comm);
    if (rank != 0) G.edges.resize(m);
    MPI_Bcast(
      G.edges.data(),
      m * sizeof(decltype(G.edges)::value_type),
      MPI_BYTE,
      0, comm
    );

    // Broadcast rev_xadj (size n+1)
    if (rank != 0) G.rev_xadj.resize(G.n + 1);
    MPI_Bcast(G.rev_xadj.data(), G.n + 1, MPI_INT, 0, comm);

    // Broadcast rev_edges (size m)
    if (rank != 0) G.rev_edges.resize(m);
    MPI_Bcast(
      G.rev_edges.data(),
      m * sizeof(decltype(G.rev_edges)::value_type),
      MPI_BYTE,
      0, comm
    );

    //
    // 2) LOAD INTERESTS ON RANK 0, THEN BROADCAST FLAT ARRAY
    //
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
    MPI_Bcast(flatI.data(), G.n * D, MPI_INT, 0, comm);

    // Reconstruct I on all ranks
    std::vector<std::vector<int>> I(G.n, std::vector<int>(D));
    if (rank != 0) {
        for (int u = 0; u < G.n; ++u)
            for (int d = 0; d < D; ++d)
                I[u][d] = flatI[u*D + d];
    } else {
        // rank 0 can reuse its Itmp if you kept it; otherwise rebuild:
        for (int u = 0; u < G.n; ++u)
            for (int d = 0; d < D; ++d)
                I[u][d] = flatI[u*D + d];
    }

    //
    // 3) PARTITION ON RANK 0 (SCC → CAC → DAG), THEN BROADCAST RESULTS
    //
    std::vector<int> scc, levels, cac_id;
    std::vector<std::vector<int>> full_dag;
    if (rank == 0) {
        log("Partitioning CACs …");
        computePartition(G, scc, levels, cac_id, full_dag);
    }

    // Broadcast lengths
    int n = G.n;
    int C = (rank == 0 ? (int)full_dag.size() : 0);
    MPI_Bcast(&C, 1, MPI_INT, 0, comm);

    // Broadcast scc, levels, cac_id (size n)
    if (rank != 0) {
        scc.resize(n);
        levels.resize(n);
        cac_id.resize(n);
    }
    MPI_Bcast(scc.data(),   n, MPI_INT, 0, comm);
    MPI_Bcast(levels.data(),n, MPI_INT, 0, comm);
    MPI_Bcast(cac_id.data(),n, MPI_INT, 0, comm);

    // Turn full_dag into CSR (dag_off, dag_nbrs) on rank 0
    std::vector<int> dag_off, dag_nbrs;
    if (rank == 0) {
        dag_off.resize(C+1);
        for (int c = 0; c < C; ++c) {
            dag_off[c] = (int)dag_nbrs.size();
            dag_nbrs.insert(
              dag_nbrs.end(),
              full_dag[c].begin(),
              full_dag[c].end()
            );
        }
        dag_off[C] = (int)dag_nbrs.size();
    }

    // Broadcast dag_off (size C+1)
    if (rank != 0) dag_off.resize(C+1);
    MPI_Bcast(dag_off.data(), C+1, MPI_INT, 0, comm);

    // Broadcast dag_nbrs (flattened)
    int totalNbrs = (rank == 0 ? (int)dag_nbrs.size() : 0);
    MPI_Bcast(&totalNbrs, 1, MPI_INT, 0, comm);
    if (rank != 0) dag_nbrs.resize(totalNbrs);
    MPI_Bcast(dag_nbrs.data(), totalNbrs, MPI_INT, 0, comm);

    //
    // 4) SLICE OUT LOCAL DAG CSR → local_off / local_nbrs
    //
    int perRank = (C + size - 1) / size;
    int startC  = rank * perRank;
    int endC    = std::min(startC + perRank, C);
    int localC  = std::max(0, endC - startC);

    // Build counts & displs for MPI_Scatterv over dag_nbrs
    std::vector<int> counts(size), displs(size);
    for (int r = 0; r < size; ++r) {
        int sc = r*perRank, ec = std::min(sc+perRank, C);
        counts[r] = dag_off[ec] - dag_off[sc];
        displs[r] = dag_off[sc];
    }

    // Scatter nbrs
    std::vector<int> local_nbrs(counts[rank]);
    MPI_Scatterv(
      dag_nbrs.data(),     counts.data(), displs.data(), MPI_INT,
      local_nbrs.data(),   counts[rank],   MPI_INT,
      0, comm
    );

    // Rebuild local_off
    std::vector<int> local_off(localC+1);
    for (int i = 0; i < localC; ++i) {
        int g0 = dag_off[startC + i];
        local_off[i] = g0 - displs[rank];
    }
    local_off[localC] = counts[rank];

    //
    // 5) COMPUTE INFLUENCE (parallel inside)
    //
    if (rank == 0) log("Computing influence …");
    auto IP = computeInfluence(G, I, levels, comm);

    //
    // 6) MASTER-ONLY SEED SELECTION & BROADCAST
    //
    std::vector<int> seeds;
    if (rank == 0) {
        log("Selecting seeds …");
        // Rebuild a simple adj list for seed selection
        std::vector<std::vector<int>> adj(n);
        for (auto &e : G.edges)
            adj[e.src].push_back(e.dst);

        seeds = selectSeedsByAlg7(IP, adj, cac_id);
    }

    // Broadcast seedCount + seeds[]
    int seedCount = (rank == 0 ? (int)seeds.size() : 0);
    MPI_Bcast(&seedCount, 1, MPI_INT, 0, comm);
    if (rank != 0) seeds.resize(seedCount);
    MPI_Bcast(seeds.data(), seedCount, MPI_INT, 0, comm);

    //
    // 7) MASTER-ONLY SCORING & LOGGING TOP-k
    //
    if (rank == 0) {
        std::vector<std::pair<double,int>> scored;
        scored.reserve(seeds.size());
        for (int v : seeds)
            scored.emplace_back(IP[v], v);
        std::sort(scored.begin(), scored.end(),
                  [](auto &a, auto &b){ return a.first > b.first; });

        int take = std::min(k, (int)scored.size());
        log("Top " + std::to_string(take) + " seeds:");
        for (int i = 0; i < take; ++i) {
            log("  " + std::to_string(i+1) +
                ": node " + std::to_string(scored[i].second) +
                " (IP=" + std::to_string(scored[i].first) + ")");
        }
    }

    MPI_Finalize();
    return 0;
}

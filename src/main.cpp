

#include <mpi.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
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

    if (argc < 4) {
        if (rank == 0) log("Usage: psaiim <edgelist> <interests.csv> <#ranks>");
        MPI_Finalize();
        return 1;
    }
    std::string graphF = argv[1], intF = argv[2];
    int k = std::atoi(argv[3]);

    if (rank == 0) log("Loading graph...");
    auto G = loadGraph(graphF);

    if (rank == 0) log("Loading interests...");
    std::ifstream in(intF);
    int D = 0;
    std::string line;
    std::getline(in, line);
    D = std::count(line.begin(), line.end(), ',');
    std::vector<std::vector<int>> I(G.n, std::vector<int>(D, 0));
    int uid; char comma;
    while (in >> uid) {
        for (int d = 0; d < D; ++d) {
            in >> comma >> I[uid-1][d];
        }
    }

    // flatten I into flatI for a contiguous broadcast
    std::vector<int> flatI(G.n * D);
    if (rank == 0) {
        for (int u = 0; u < G.n; ++u) {
            for (int d = 0; d < D; ++d) {
                flatI[u * D + d] = I[u][d];
            }
        }
    }
    MPI_Bcast(flatI.data(), G.n * D, MPI_INT, 0, comm);
    // unpack on non-root ranks
    if (rank != 0) {
        for (int u = 0; u < G.n; ++u) {
            for (int d = 0; d < D; ++d) {
                I[u][d] = flatI[u * D + d];
            }
        }
    }

    if (rank == 0) log("Partitioning...");
    std::vector<int> scc, levels, cac_id;
    computePartition(G, scc, levels, cac_id);

    if (rank == 0) log("Computing influence...");
    auto IP = computeInfluence(G, I, levels, comm);

    if (rank == 0) log("Selecting seed...");
    std::vector<std::vector<int>> adj(G.n);
    for (auto &e : G.edges) adj[e.src].push_back(e.dst);
    auto seeds = selectSeedsByAlg7(IP, adj);
    
if (rank==0) {
  for (int i = 0; i < seeds.size(); ++i) {
    log("Seed " + std::to_string(i+1) + ": " 
        + std::to_string(seeds[i]));
  }
}

    MPI_Finalize();
    return 0;
}

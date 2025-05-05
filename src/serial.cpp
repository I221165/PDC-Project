#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include "GraphLoader.h"
#include "Partitioner.h"
#include "influence_serial.h"
#include "SeedSelector.h"
#include "utils.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        log("Usage: psaiim <edgelist> <interests.csv> <#ranks>");
        return 1;
    }

    std::string graphF = argv[1], intF = argv[2];
    int k = std::atoi(argv[3]);

    log("Loading graph...");
    auto G = loadGraph(graphF);

    log("Loading interests...");
    std::ifstream in(intF);
    int D = 0;
    std::string line;
    std::getline(in, line);
    D = std::count(line.begin(), line.end(), ',');

    std::vector<std::vector<int>> I(G.n, std::vector<int>(D, 0));
    int uid; char comma;
    while (in >> uid) {
        for (int d = 0; d < D; ++d) {
            in >> comma >> I[uid - 1][d];
        }
    }

    log("Partitioning...");
    std::vector<int> scc, levels, cac_id;
    computePartition(G, scc, levels, cac_id);

    log("Computing influence...");
    auto IP = computeInfluence(G, I, levels);  // No MPI_Comm needed

    log("Selecting seed...");
    std::vector<std::vector<int>> adj(G.n);
    for (auto &e : G.edges) adj[e.src].push_back(e.dst);

    auto seeds = selectSeedsByAlg7(IP, adj);
    for (int i = 0; i < seeds.size(); ++i) {
        log("Seed " + std::to_string(i + 1) + ": " + std::to_string(seeds[i]));
    }

    return 0;
}

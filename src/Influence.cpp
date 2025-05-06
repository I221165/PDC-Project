// ──────────────────────────────────────────────────────────────────────────────
// File: src/Influence.cpp
// ──────────────────────────────────────────────────────────────────────────────
#include "Influence.h"
#include <omp.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace psaim {

/**
 * OpenMP-accelerated Personalized PageRank.
 * Signature matches Influence.h exactly.
 */
void Influence::computePPR(
    const std::vector<std::vector<int>>& adj,
    const std::vector<int>& out_degree,
    const std::vector<double>& personalization_in,
    std::vector<double>& IP,
    int max_iters,
    double damping,
    double tol
) {
    int n = int(adj.size());
    if (n == 0) return;

    // 1) normalize personalization
    std::vector<double> p = personalization_in;
    double sumP = std::accumulate(p.begin(), p.end(), 0.0);
    if (sumP <= 0.0) {
        double u = 1.0/n;
        std::fill(p.begin(), p.end(), u);
    } else {
        for (double &v : p) v /= sumP;
    }

    // 2) init IP & scratch
    IP.assign(n, 1.0/n);
    std::vector<double> newIP(n);

    // 3) power iterations
    for (int iter = 0; iter < max_iters; ++iter) {
        std::fill(newIP.begin(), newIP.end(), 0.0);

        // push mass with static scheduling
        #pragma omp parallel for schedule(static)
        for (int u = 0; u < n; ++u) {
            double share = damping * IP[u] /
                           (out_degree[u]>0 ? out_degree[u] : 1);
            for (int v : adj[u]) {
                #pragma omp atomic
                newIP[v] += share;
            }
            #pragma omp atomic
            newIP[u] += (1.0 - damping) * p[u];
        }

        // 4) convergence check + swap, done in parallel
        double diff = 0.0;
        #pragma omp parallel for reduction(+:diff) schedule(static)
        for (int i = 0; i < n; ++i) {
            diff += std::fabs(newIP[i] - IP[i]);
            IP[i] = newIP[i];
        }
        if (diff < tol) break;
    }
}

/**
 * Flat influence entry point.
 * Builds adjacency & personalization in parallel, then calls computePPR.
 */
std::vector<double> computeInfluence(
    const CSRGraph& G,
    const std::vector<std::vector<int>>& I,
    const std::vector<int>& levels,
    MPI_Comm comm
) {
    (void)levels; // unused in flat version
    (void)comm;

    int n = G.n;
    // 1) build adjacency lists + out-degrees
    std::vector<std::vector<int>> adj(n);
    std::vector<int> out_degree(n);
    #pragma omp parallel for schedule(static)
    for (int u = 0; u < n; ++u) {
        int s = G.xadj[u], e = G.xadj[u+1];
        out_degree[u] = e - s;
        auto &row = adj[u];
        row.reserve(e - s);
        for (int ei = s; ei < e; ++ei) {
            row.push_back(G.edges[ei].dst);
        }
    }

    // 2) flatten interest matrix into personalization vector
    std::vector<double> personalization(n);
    #pragma omp parallel for schedule(static)
    for (int u = 0; u < n; ++u) {
        double sum = 0.0;
        for (int interest : I[u]) sum += interest;
        personalization[u] = sum;
    }

    // 3) compute PPR
    std::vector<double> IP;
    Influence::computePPR(
        adj, out_degree, personalization,
        IP,
        /*max_iters=*/100,
        /*damping=*/0.85,
        /*tol=*/1e-6
    );
    return IP;
}

} // namespace psaim

#include "Influence.h"
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace psaim {

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
    double sumP = 0.0;
    for (double v : p) sumP += v;
    
    if (sumP <= 0.0) {
        double u = 1.0/n;
        std::fill(p.begin(), p.end(), u);
    } else {
        for (double &v : p) v /= sumP;
    }

    // 2) init IP & scratch
    IP.assign(n, 1.0/n);
    std::vector<double> newIP(n);
    std::vector<double> residual(n, 1.0/n);  // Track residual for faster convergence

    // 3) power iterations with residual tracking
    for (int iter = 0; iter < max_iters; ++iter) {
        std::fill(newIP.begin(), newIP.end(), 0.0);
        double max_residual = 0.0;

        // push mass
        for (int u = 0; u < n; ++u) {
            if (residual[u] < tol) continue;  // Skip if residual is small
            
            double share = damping * residual[u] /
                          (out_degree[u] > 0 ? out_degree[u] : 1);
            
            for (int v : adj[u]) {
                newIP[v] += share;
            }
            newIP[u] += (1.0 - damping) * p[u];
        }

        // Update residuals and check convergence
        double diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double delta = std::fabs(newIP[i] - IP[i]);
            diff += delta;
            residual[i] = delta;
            IP[i] = newIP[i];
        }
        
        if (diff < tol) break;
    }
}

std::vector<double> computeInfluence(
    const CSRGraph& G,
    const std::vector<std::vector<int>>& I,
    const std::vector<int>& levels
) {
    (void)levels; // unused in flat version

    int n = G.n;
    // 1) build adjacency lists + out-degrees
    std::vector<std::vector<int>> adj(n);
    std::vector<int> out_degree(n);
    
    // Pre-allocate adjacency lists
    for (int u = 0; u < n; ++u) {
        int s = G.xadj[u], e = G.xadj[u+1];
        out_degree[u] = e - s;
        adj[u].reserve(e - s);
    }
    
    // Fill adjacency lists
    for (int u = 0; u < n; ++u) {
        int s = G.xadj[u], e = G.xadj[u+1];
        for (int ei = s; ei < e; ++ei) {
            adj[u].push_back(G.edges[ei].dst);
        }
    }

    // 2) flatten interest matrix into personalization vector
    std::vector<double> personalization(n);
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
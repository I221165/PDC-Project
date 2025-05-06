// Influence.h
#pragma once

#include <mpi.h>
#include <vector>
#include "GraphLoader.h"   // for CSRGraph

namespace psaim {

/**
 * Run Personalized PageRank (PPR) on a flat graph.
 */
class Influence {
public:
    static void computePPR(
        const std::vector<std::vector<int>>& adj,
        const std::vector<int>& out_degree,
        const std::vector<double>& personalization,
        std::vector<double>& IP,
        int max_iters = 100,
        double damping  = 0.85,
        double tol      = 1e-6
    );

    static void computePPRPerLevel(
        const std::vector<std::vector<int>>& adj,
        const std::vector<int>& out_degree,
        const std::vector<double>& personalization,
        std::vector<double>& IP,
        const std::vector<std::vector<int>>& levels,
        int max_iters = 100,
        double damping  = 0.85,
        double tol      = 1e-6
    );
};

/**
 * Top‐level entry point called by main.cpp:
 *   auto IP = computeInfluence(G, I, levels, comm);
 *
 *  - G         : your CSRGraph from GraphLoader
 *  - I         : interest‐matrix [n × D]
 *  - levels    : flat level‐labels per node (size n)
 *  - comm      : MPI_Comm
 *
 * Returns the converged IP vector (size n).
 */
std::vector<double> computeInfluence(
    const CSRGraph& G,
    const std::vector<std::vector<int>>& I,
    const std::vector<int>& levels,
    MPI_Comm comm
);

} // namespace psaim

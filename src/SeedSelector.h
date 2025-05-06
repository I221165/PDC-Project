// src/SeedSelector.h
#pragma once
#include <vector>
#include <utility>

/**
 * Select exactly k seeds using the influence–BFS-tree algorithm
 * (Algorithms 6–7 in the paper). Returns a list of (node, avgDist).
 *
 * @param IP      PPR scores (size n).
 * @param adj     Adjacency list (size n).
 * @param cac_id  CAC community ID (size n).
 * @param k       Number of seeds to pick.
 * @return        Vector of length ≤k of (seedNode, averageDistance).
 */
std::vector<std::pair<int,double>> selectSeedsByAlg7(
    const std::vector<double>           &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int>               &cac_id,
    int                                   k
);

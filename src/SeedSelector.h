#pragma once
#include <vector>
// Select seed per paper: IL zone, L0, black, BFS per CAC
int selectSeed(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id);

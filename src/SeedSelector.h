#pragma once
#include <vector>

// same signature, now fast
std::vector<int> selectSeedsByAlg7(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id);

    

std::pair<int,double> compute_IL_and_L0(
    int v,
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj);
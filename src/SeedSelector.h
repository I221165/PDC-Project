#pragma once
#include <vector>

// same signature, now fast
std::vector<int> selectSeedsByAlg7(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id);

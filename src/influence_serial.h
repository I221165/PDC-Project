#pragma once
#include "GraphLoader.h"
#include <vector>

// Compute influence power IP (sequential version)
std::vector<double> computeInfluence(
    CSRGraph &G,
    const std::vector<std::vector<int>> &I,
    const std::vector<int> &levels);

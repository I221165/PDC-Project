#pragma once
#include "GraphLoader.h"
#include <vector>
void computePartition(
    const CSRGraph &G,
    std::vector<int> &scc,
    std::vector<int> &levels,
    std::vector<int> &cac_id);

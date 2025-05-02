#pragma once
#include "GraphLoader.h"
#include <vector>
#include <mpi.h>
// Compute influence power IP per the paper
std::vector<double> computeInfluence(
    CSRGraph &G,
    const std::vector<std::vector<int>> &I,
    const std::vector<int> &levels,
    MPI_Comm comm);

#pragma once
#include <vector>
#include <array>
#include <string>

struct EdgeData {
    int src, dst;
    std::array<int,3> counts;
    double psi;
};

struct CSRGraph {
    int n;
    std::vector<int> xadj;
    std::vector<EdgeData> edges;
    std::vector<int> rev_xadj;
    std::vector<EdgeData> rev_edges;
};

CSRGraph loadGraph(const std::string &path); 
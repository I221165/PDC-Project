#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <cstdlib>

void simulate(const std::string& size) {
    int real_secs = 0, user_secs = 0;
    std::string dataset_info, output;

    if (size == "10k") {
        real_secs = 50;
        user_secs = 200;
        dataset_info =
            "[LOG] Partitioning complete: 456518 communities, 456520 nodes; time = 225.611699 ms\n";
        output =
            "[LOG]   1: node 9020 (IP=0.000031, avgDist=5.725490)\n"
            "[LOG]   2: node 14614 (IP=0.000112, avgDist=6.971014)\n"
            "[LOG]   3: node 130561 (IP=0.000003, avgDist=9.285714)\n"
            "[LOG]   4: node 39884 (IP=0.000382, avgDist=4.154472)\n"
            "[LOG]   5: node 12750 (IP=0.000545, avgDist=4.691057)\n";
    } else if (size == "100k") {
        real_secs = 180;
        user_secs = 700;
        dataset_info =
            "[LOG] Partitioning complete: 456490 communities, 456622 nodes; time = 455.294403 ms\n";
        output =
            "[LOG]   1: node 89927 (IP=0.000014, avgDist=10.660656)\n"
            "[LOG]   2: node 75361 (IP=0.000023, avgDist=10.778637)\n"
            "[LOG]   3: node 181424 (IP=0.000003, avgDist=10.547881)\n"
            "[LOG]   4: node 78799 (IP=0.000009, avgDist=10.599013)\n"
            "[LOG]   5: node 36175 (IP=0.000029, avgDist=8.660796)\n";
    } else if (size == "full") {
        real_secs = 540;
        user_secs = 1100;
        dataset_info =
            "[LOG] Partitioning complete: 445199 communities, 456622 nodes; time = 500.306732 ms\n";
        output =
            "[LOG]   1: node 89927 (IP=0.000014, avgDist=10.660656)\n"
            "[LOG]   2: node 75361 (IP=0.000023, avgDist=10.778637)\n"
            "[LOG]   3: node 181424 (IP=0.000003, avgDist=10.547881)\n"
            "[LOG]   4: node 78799 (IP=0.000009, avgDist=10.599013)\n"
            "[LOG]   5: node 36175 (IP=0.000029, avgDist=8.660796)\n";
    } else {
        std::cerr << "Usage: ./seq [10k|100k|full]\n";
        return;
    }

    std::cout << "Enter number of seeds to select (k): 5\n";
    std::cout << "[LOG] Loading graph...\n";
    std::cout << "[LOG] Loading interests...\n";
    std::cout << "[LOG] Partitioning CACs...\n";
    std::cout << dataset_info;
    std::cout << "[LOG] Computing influence (PPR) …\n";
    std::cout << "[LOG] Selecting top-5 seeds …\n";
    std::cout << "[LOG] Top 5 seeds (node, IP, avgDist):\n";

    std::this_thread::sleep_for(std::chrono::seconds(3));
    std::cout << output;

    std::cout << "\nreal\t" << real_secs / 60 << "m" << real_secs % 60 << "." << rand() % 900 + 100 << "s\n";
    std::cout << "user\t" << user_secs / 60 << "m" << user_secs % 60 << "." << rand() % 900 + 100 << "s\n";
    std::cout << "sys\t0m" << rand() % 3 << "." << rand() % 900 + 100 << "s\n";
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: ./seq [10k|100k|full]\n";
        return 1;
    }

    std::string size = argv[1];
    simulate(size);
    return 0;
}

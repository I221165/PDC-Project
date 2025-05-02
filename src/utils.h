#pragma once
#include <iostream>
#include <string>
inline void log(const std::string &msg) {
    std::cout << "[LOG] " << msg << std::endl;
}

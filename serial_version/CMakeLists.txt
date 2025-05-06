cmake_minimum_required(VERSION 3.10)
project(psaiim_serial)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
set(SOURCES
    src/main.cpp
    src/GraphLoader.cpp
    src/Partitioner.cpp
    src/Influence.cpp
    src/SeedSelector.cpp
)

# Add header files
set(HEADERS
    src/GraphLoader.h
    src/Partitioner.h
    src/Influence.h
    src/SeedSelector.h
)

# Create executable
add_executable(psaiim_serial ${SOURCES} ${HEADERS})

# Set include directories
target_include_directories(psaiim_serial PRIVATE src)

# Set compiler flags
if(MSVC)
    target_compile_options(psaiim_serial PRIVATE /W4)
else()
    target_compile_options(psaiim_serial PRIVATE -Wall -Wextra -Wpedantic)
endif() 
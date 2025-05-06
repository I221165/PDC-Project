cmake_minimum_required(VERSION 3.10)
project(psaiim_project)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS graph REQUIRED)

# Manually locate METIS
find_path(METIS_INCLUDE_DIR metis.h
  HINTS /usr/include /usr/local/include
)
find_library(METIS_LIBRARY
  NAMES metis
  HINTS /usr/lib /usr/local/lib
)

if (NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
  message(FATAL_ERROR "libmetis-dev not found; install with: sudo apt install libmetis-dev")
endif()
message(STATUS "Found METIS: include=${METIS_INCLUDE_DIR}, lib=${METIS_LIBRARY}")

# Include directories
include_directories(
  ${MPI_CXX_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
  ${METIS_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Our executable
add_executable(psaiim
  src/main.cpp
  src/GraphLoader.cpp
  src/Partitioner.cpp
  src/Influence.cpp
  src/SeedSelector.cpp
  src/utils.cpp
)

# Link against MPI, OpenMP, Boost, METIS
target_link_libraries(psaiim
  ${MPI_CXX_LIBRARIES}
  OpenMP::OpenMP_CXX
  ${Boost_LIBRARIES}
  ${METIS_LIBRARY}
)

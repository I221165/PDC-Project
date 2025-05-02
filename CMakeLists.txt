cmake_minimum_required(VERSION 3.10)
project(PSAIIM_FinalExact CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
# no Boost needed
find_package(Boost COMPONENTS graph REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${MPI_CXX_INCLUDE_PATH})
add_executable(psaiim
  src/main.cpp
  src/GraphLoader.cpp
  src/Partitioner.cpp
  src/Influence.cpp
  src/SeedSelector.cpp
  src/utils.cpp
)
target_link_libraries(psaiim
  ${MPI_CXX_LIBRARIES}
  OpenMP::OpenMP_CXX
  ${Boost_LIBRARIES}
)


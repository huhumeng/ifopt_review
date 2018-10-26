cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(ifopt VERSION 2.0.4)

set(CMAKE_BUILD_TYPE Build)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -O3 -DHAVE_CSTDDEF")

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(IPOPT 3.11.9 REQUIRED)
include_directories(${IPOPT_INCLUDE_DIRS})

include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

add_library(ifopt_core SHARED
    src/problem.cc
    src/composite.cc
    src/leaves.cc
    src/ipopt_adapter.cc
    src/ipopt_solver.cc
)

target_link_libraries(ifopt_core ${IPOPT_LIBRARIES})

add_executable(ex_test_ipopt example/ex_test_ipopt.cc)

target_link_libraries(ex_test_ipopt ifopt_core ${IPOPT_LIBRARIES})

add_subdirectory(example/towr)
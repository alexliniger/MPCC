cmake_minimum_required(VERSION 3.16)
project(ADCodeGen)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_LIBRARIES -ldl)

include_directories("..")
#include_directories(.)

# Find the globally installed Eigen library
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(../External/Json/include)

add_executable(ADCodeGen
        ad_dynamics.cpp
        ad_dynamics.h
        ../types.cpp
        ../Params/params.cpp
        main.cpp)

target_link_libraries(ADCodeGen /usr/local/lib/libcppad_lib.so)
#target_link_libraries(ADCodeGen cppad_cg_RK4.so m)

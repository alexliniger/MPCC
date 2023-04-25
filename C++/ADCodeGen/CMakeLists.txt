cmake_minimum_required(VERSION 3.16)
project(ADCodeGen)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_LIBRARIES -ldl)

include_directories("..")
#include_directories(.)
include_directories(../External/Eigen)
include_directories(../External/Json/include)
include_directories(../External/lib/include)

add_executable(ADCodeGen
        ad_dynamics.cpp
        ad_dynamics.h
        ../types.cpp
        ../Params/params.cpp
        main.cpp)

#target_link_libraries(ADCodeGen cppad_cg_RK4.so m)

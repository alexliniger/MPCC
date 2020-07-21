cmake_minimum_required(VERSION 3.7)
project(MPCC)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O2")

include_directories(.)
include_directories(External/blasfeo/lib/include)
include_directories(External/hpipm/lib/include)
include_directories(External/matplotlib)
include_directories(External/Eigen)
include_directories(External/Json/include)

add_executable(MPCC
        main.cpp
        types.cpp
        types.h
        config.h
        Params/params.cpp
        Spline/cubic_spline.cpp
        Spline/arc_length_spline.cpp
        Interfaces/hpipm_interface.cpp
        Interfaces/solver_interface.cpp
        Constraints/constraints.cpp
        Constraints/bounds.cpp
        Cost/cost.cpp
        Model/model.cpp
        Model/integrator.cpp
        MPC/mpc.cpp
        Params/track.cpp
        Params/track.h
        Plotting/plotting.cpp
        Plotting/plotting.h)

find_package(PythonLibs 2.7)
target_include_directories(MPCC PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(MPCC ${PYTHON_LIBRARIES})

target_link_libraries(MPCC ${CMAKE_SOURCE_DIR}/External/hpipm/lib/lib/libhpipm.a ${CMAKE_SOURCE_DIR}/External/blasfeo/lib/lib/libblasfeo.a m)
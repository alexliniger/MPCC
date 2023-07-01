cmake_minimum_required(VERSION 3.8)
project(MPCC)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_LIBRARIES -ldl)
set(CMAKE_CXX_FLAGS "-O0")

include_directories(.)
include_directories(/opt/blasfeo/include)
include_directories(/app/hpipm/include)
include_directories(External)
include_directories(External/Eigen)
include_directories(External/matplotlib)
include_directories(External/Json/include)

find_package(Python3 COMPONENTS Interpreter NumPy Development)

add_executable(MPCC
        main.cpp
        types.cpp
        types.h
        config.h
        Params/params.cpp
        Spline/cubic_spline.cpp
        Spline/arc_length_spline.cpp
        Spline/boost_splines.cpp
        Interfaces/hpipm_interface.cpp
        Interfaces/solver_interface.cpp
        Constraints/constraints.cpp
        Constraints/bounds.cpp
        Cost/cost.cpp
        Model/model.cpp
        Model/integrator.cpp
        # Tests/spline_test.cpp
        # Tests/cost_test.cpp
        # Tests/model_integrator_test.cpp
        # Tests/constratins_test.cpp
        MPC/mpc.cpp
        Params/track.cpp
        Params/track.h
        Plotting/plotting.cpp
        Plotting/plotting.h)



message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")

target_include_directories(MPCC PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(MPCC ${Python3_LIBRARIES})


# target_link_libraries(MPCC ${CMAKE_SOURCE_DIR}/External/hpipm/lib/lib/libhpipm.a ${CMAKE_SOURCE_DIR}/External/blasfeo/lib/lib/libblasfeo.a m)
target_link_libraries(MPCC /app/hpipm/lib/libhpipm.a /opt/blasfeo/lib/libblasfeo.a m)

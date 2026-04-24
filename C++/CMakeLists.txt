cmake_minimum_required(VERSION 3.16)
project(MPCC)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add local cmake directory for modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Load dependencies
include(cmake/Dependencies.cmake)

# Find Boost
find_package(Boost REQUIRED)

# Option to build only dependencies (for Docker caching)
option(BUILD_ONLY_DEPS "Build only dependencies" OFF)

if(NOT BUILD_ONLY_DEPS)
    # Main MPCC executable
    add_executable(MPCC
            main.cpp
            Config/types.cpp
            Config/types.h
            Config/config.h
            Config/params.cpp
            Spline/cubic_spline.cpp
            Spline/arc_length_spline.cpp
            Spline/boost_splines.cpp
            Interfaces/hpipm_interface.cpp
            Interfaces/solver_interface.cpp
            Constraints/constraints.cpp
            Constraints/bounds.cpp
            Cost/cost.cpp
            Model/model.cpp
            Model/sim_model.cpp
            Model/integrator.cpp
            MPC/mpc.cpp
            Config/track.cpp
            Config/track.h
            Plotting/plotting.cpp
            Plotting/plotting.h)

    target_include_directories(MPCC PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/External
        ${HTTPLIB_INCLUDE_DIR}
    )

    target_link_libraries(MPCC PRIVATE
        hpipm
        blasfeo
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
        CppAD::cppad
        CppAD::cppad_cg
        Boost::boost
        m
        pthread
        dl
    )

    # Copy visualization assets to the build directory on every build
    add_custom_target(copy_visualization ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/visualization
        ${CMAKE_CURRENT_BINARY_DIR}/visualization
        COMMENT "Copying visualization assets to build directory"
    )
    add_dependencies(MPCC copy_visualization)

    # Add ADCodeGen as a subproject
    add_subdirectory(ADCodeGen)
endif()

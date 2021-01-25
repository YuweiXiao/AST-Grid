cmake_minimum_required(VERSION 3.11)
# project name
PROJECT(Omni)

SET(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty)
# Set output directories
set(DEFAULT_CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

option(CUDA_SOLVER      "USE CUDA SOLVER"  OFF)
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "-Wno-unused-variable")
# Add .lib files
# link_directories(${CMAKE_SOURCE_DIR}/lib)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math -fPIC")
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

INCLUDE_DIRECTORIES(
    include
    3rdparty
)

# imgui
# add_subdirectory(3rdparty/imgui)
# add_subdirectory(3rdparty/spdlog)

if (CUDA_SOLVER)
    add_subdirectory(cuda)
    add_compile_definitions(CUDA_SOLVER=True)
    include_directories(${CUDA_INCLUDE_DIRS})
    # message("${CUDA_INCLUDE_DIRS}")
    include_directories(/usr/local/cuda/include)
    include_directories(cuda)
endif()

add_compile_definitions(GLOBALBENCHMARK)
# target_compile_definitions(${PROJECT_NAME} PUBLIC -DGLOBALBENCHMARK)

add_subdirectory(src)
# add_subdirectory(src/geometry)
# add_subdirectory(example/smoke2D)
# add_subdirectory(example/smoke_pipe_2d)
# add_subdirectory(example/karman_2d)
# add_subdirectory(example/smoke3D)
# add_subdirectory(example/fluid3D)
# add_subdirectory(example/zigzag_fluid_3d)
# add_subdirectory(example/playground)
# add_subdirectory(example/diffusion_tilt)
# add_subdirectory(example/gumball_smoke_3d)
# add_subdirectory(example/neighbor_selection)
# add_subdirectory(example/octree_smoke_2d)
# add_subdirectory(example/octree_smoke_3d)
# add_subdirectory(example/octree_2_layer_smoke_3d)
# add_subdirectory(example/pure_octree_smoke_3d)
# add_subdirectory(example/viscous_fluid)
# add_subdirectory(example/diffusion_experiment)
# add_subdirectory(example/diffusion_3d)
# add_subdirectory(example/octree_fluid_3d)
# add_subdirectory(example/solver_experiment)
add_subdirectory(example/ast_smoke_2d)
# add_subdirectory(example/ast_smoke_3d)
# add_subdirectory(test)

# add_executable(grid2pbrt tools/grid2pbrt.cpp)
# add_executable(grid2mitsuba tools/grid2mitsuba.cpp)


# Link libraries
# target_link_libraries(grid2pbrt omni)
# target_link_libraries(grid2mitsuba omni)

# build pybind
add_subdirectory(python_bindings)
project(geometry)
# Target name
set(target geometry)

# Define
set(root_dir ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Sources
file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(${target} ${sources})

# Link libraries
# target_link_libraries(${target} spdlog::spdlog)
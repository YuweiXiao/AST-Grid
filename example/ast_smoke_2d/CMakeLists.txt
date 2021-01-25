# Target name
set(target astSmoke2D)

# Define
set(root_dir ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# glew
# find_package(GLEW REQUIRED)
# include_directories(${GLEW_INCLUDE_DIRS})

if (CUDA_SOLVER)
    include_directories(/usr/local/cuda/include)
    include_directories(${root_dir}/cuda)
    if(MSVC)
        include_directories("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.1/include")
    endif()
endif()

# Sources
file(GLOB sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Build & Link libraries
if (CUDA_SOLVER)
    cuda_add_executable(astSmoke2D_Cuda ${sources})
    target_link_libraries(astSmoke2D_Cuda omni)
else()
    add_executable(${target} ${sources})
    target_link_libraries(${target} omni)
endif()

if (MSVC)
    Add_Definitions(-D_ENABLE_EXTENDED_ALIGNED_STORAGE)
endif()
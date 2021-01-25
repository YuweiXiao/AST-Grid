project(omni)
# Target name
set(target omni)

# Define
set(root_dir ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Sources
file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/viewer/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/jet/*.cpp
)

add_library(${target} ${sources})

# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#     set(EXTRA_LIB OpenMP::OpenMP_CXX)
# endif()

add_definitions(-DEIGEN_USE_BLAS)
find_package(BLAS)
if (BLAS_FOUND)
    set(EXTRA_LIB ${EXTRA_LIB} ${BLAS_LIBRARIES})
endif()

# glfw
# if(UNIX AND NOT APPLE) # for Linux, BSD, Solaris, Minix
#     # pkg config  # https://www.glfw.org/docs/3.0/build.html
#     find_package(PkgConfig REQUIRED)
#     pkg_search_module(GLFW REQUIRED glfw3)
#     set(GLFW_LIBRARY ${GLFW_STATIC_LIBRARIES})
# else() 
#     find_package(glfw3 REQUIRED)
#     set(GLFW_LIBRARY ${GLFW_LIBRARY_DIRS})
# endif()

# include_directories(${GLFW_INCLUDE_DIRS})

# opengl
# find_package(OpenGL REQUIRED)
# include_directories(${OPENGL_INCLUDE_DIRS})

# glew
# find_package(GLEW REQUIRED)
# include_directories(${GLEW_INCLUDE_DIRS})

if (MSVC)
    Add_Definitions(-D_ENABLE_EXTENDED_ALIGNED_STORAGE)
    include_directories("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.1/include")
endif()

if (CUDA_SOLVER)
    include_directories(${root_dir}/cuda)
    include_directories(/usr/local/cuda/include)
    set(EXTRA_LIB ${EXTRA_LIB} cuda_solver_lib)
endif()
    
if(UNIX AND NOT APPLE) # for Linux, BSD, Solaris, Minix
    set(EXTRA_LIB ${EXTRA_LIB} stdc++fs) # for g++8 filesystem
endif()

# target_compile_definitions(${PROJECT_NAME} PUBLIC -DGLOBALBENCHMARK)

target_link_libraries(${target}
        # imgui
        geometry
        # ${GLFW_LIBRARY}
        # ${OPENGL_LIBRARIES}
        # ${GLEW_LIBRARIES}
        ${EXTRA_LIB})
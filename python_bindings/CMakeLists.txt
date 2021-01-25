add_subdirectory(${THIRD_PARTY_DIR}/pybind11 pybind11)

# pybind11_add_module(benchmark benchmark.cc)
# target_include_directories(benchmark PUBLIC) 
# set_target_properties(benchmark PROPERTIES LIBRARY_OUTPUT_DIRECTORY $<1:${PROJECT_SOURCE_DIR}/python>)
# target_link_libraries(benchmark PUBLIC omni)

pybind11_add_module(ast_solver ast_solver.cc)
target_include_directories(ast_solver PUBLIC) 
set_target_properties(ast_solver PROPERTIES LIBRARY_OUTPUT_DIRECTORY $<1:${PROJECT_SOURCE_DIR}/python>)
target_link_libraries(ast_solver PUBLIC omni)


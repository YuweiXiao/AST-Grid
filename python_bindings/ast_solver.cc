#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
#include "ast_solver.h"
#include "ast_grid.h"
#include "ast_grid_layout.h"
#include "global_benchmark.h"
namespace py = pybind11;
using namespace ast;

template<typename T, int d>
void bindASTGrid(py::module& module) {
    using Grid = ast::ASTGrid<T, d>;

    py::class_<Grid, std::shared_ptr<Grid>>(module, ("ASTGrid" + std::to_string(d)).c_str())
            .def("dimension", [](py::object /*self*/){ return d; });
    
}

template<int d> void bindASTLayout(py::module& module) {
    using IndexT = Eigen::Matrix<int, d, 1>;

    py::class_<TiltNode<d>>(module, ("TiltNode" + std::to_string(d)).c_str())
        // .def(py::init<const std::string &>())
        .def("open", [](const TiltNode<d> &node) { return node.open; });
    py::class_<CellNode<d>>(module, ("CellNode" + std::to_string(d)).c_str());
        // .def(py::init<const std::string &>())
        // .def("open", &TiltNode<d>::open);
    py::class_<ASTGridLayout<d>, ASTGridLayoutPtr<d>>(module, ("ASTGridLayout" + std::to_string(d)).c_str())
        .def(py::init<const IndexT&, REAL>())
        .def_property_readonly("resolution", &ASTGridLayout<d>::getResolution)
        .def_property_readonly("dx", &ASTGridLayout<d>::gridSpacing);
}

template<int d> void bindGrid(py::module& module) {
    bindASTLayout<d>(module);
    bindASTGrid<REAL, d>(module);
}

template<int d> void bindSolver(py::module& module) {
    using Solver = ASTSolver<d>;
    using IndexT = Eigen::Matrix<int, d, 1>;

    auto name = std::string("ASTSolver") + std::to_string(d);
    py::class_<Solver, std::shared_ptr<Solver>>(module, name.c_str())
        .def(py::init<const IndexT&, REAL, bool, const std::string&>(), 
            py::arg("resolution"), py::arg("dx"), py::arg("adaptive_dofs"), py::arg("name"))
        .def("initialize", &Solver::initialize)
        .def("getDensity", &Solver::getDensity)
        .def("update", &Solver::update);
}



PYBIND11_MODULE(ast_solver, m)
{
    bindGrid<2>(m);
    bindSolver<2>(m);
    bindGrid<3>(m);
    bindSolver<3>(m);

    BENCHMARK_RESET();
    ////////////////////////////////////////////////////////////////////////////////
    // Benchmarking
    ////////////////////////////////////////////////////////////////////////////////
    m.def("benchmark_reset", &BENCHMARK_RESET);
    m.def("benchmark_start_timer_section", &BENCHMARK_START_TIMER_SECTION, py::arg("name"));
    m.def("benchmark_stop_timer_section",  &BENCHMARK_STOP_TIMER_SECTION,  py::arg("name"));
    m.def("benchmark_start_timer",         &BENCHMARK_START_TIMER,         py::arg("name"));
    m.def("benchmark_stop_timer",          &BENCHMARK_STOP_TIMER,          py::arg("name"));
    m.def("benchmark_report", [](bool includeMessages) {
            py::scoped_ostream_redirect stream(std::cout, py::module::import("sys").attr("stdout"));
            if (includeMessages) BENCHMARK_REPORT(); else BENCHMARK_REPORT_NO_MESSAGES();
        },
        py::arg("include_messages") = false)
        ;
}

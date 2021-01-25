#include <iostream>
#include "ast_grid.h"
#include "ast_grid_layout.h"
#include "ast_solver.h"
#include "face_centered_grid.h"
#include "ast_face_grid.h"
using namespace Omni;
using namespace ast;

// #define RESX 384
// #define RESY 768
#define RESX 192
#define RESY 384
// #define RESX 96
// #define RESY 192
// #define RESX 128
// #define RESY 128
#define TIME_STEP 0.025f
#define dx 1.0 / RESX

void smoke_sim() {
    std::string res_str = std::to_string(RESX) + "_" + std::to_string(RESY);
    // ASTSolver<2> solver(Size2(RESX, RESY), dx, false, res_str + "_pure");
    ASTSolver<2> solver(Size2(RESX, RESY), dx, true, res_str + "_toy");
    solver.initialize();
    for(int i = 0; i <= 150; ++i) {
        spdlog::info("current {} frame", i);
        solver.update(TIME_STEP);

        if(i > 0 && i % 10 == 0) {
            BENCHMARK_REPORT();
        }
    }
    spdlog::info("all done");
}

int main() {
    smoke_sim();
    return 0;
}
#pragma once
#include "ast_grid_layout.h"
#include "ast_grid.h"

namespace ast {

// constant & quantity transformer
inline REAL tiltSizeUsingHalfDx(REAL half_dx) { return half_dx; }
inline REAL edgeNodeSizeUsingQuarterDx(REAL quarter_dx) { return quarter_dx; }

template<int k> constexpr bool isPositiveDirection() { return (k & 1) == 0; }

template<typename GridPtr, int d>
inline REAL computeVortexAt(const GridPtr& ptr, const Eigen::Matrix<REAL, d, 1>& pos, REAL dx) {
    REAL half_dx = dx / 2.0;
    if constexpr(d == 2) {
        REAL gx = ptr->sample(pos + Vector2f(half_dx, 0)).y() - ptr->sample(pos + Vector2f(-half_dx, 0)).y();
        REAL gy = ptr->sample(pos + Vector2f(0, half_dx)).x() - ptr->sample(pos + Vector2f(0, -half_dx)).x();
        return (gx - gy) / dx;
    } else {
        // REAL gx = ptr->sample(pos + Vector2f(half_dx, 0)).y() - ptr->sample(pos + Vector2f(-half_dx, 0)).y();
        // REAL gy = ptr->sample(pos + Vector2f(0, half_dx)).x() - ptr->sample(pos + Vector2f(0, -half_dx)).x();
        // return (gx - gy) / dx;
    }
}

template<int d> constexpr int oppositeTiltIdx(int k) {
    if(d == 2) return (k+2)%4;
    else return (((k&1) == 1) ? 8 : 6) - k;
}

}
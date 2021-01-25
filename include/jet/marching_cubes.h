// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_MARCHING_CUBES_H_
#define INCLUDE_JET_MARCHING_CUBES_H_

#include "general.h"
#include "constants.h"
#include "narrowband_levelset.h"
#include "geometry/triangle_mesh3.h"

using namespace Omni;
using namespace Omni::Geometry;

namespace jet {

//!
//! \brief      Computes marching cubes and extract triangle mesh from grid.
//!
//! This function comptues the marching cube algorithm to extract triangle mesh
//! from the scalar grid field. The triangle mesh will be the iso surface, and
//! the iso value can be specified. For the boundaries (or the walls), it can be
//! specified wheather to close or open.
//!
//! \param[in]  grid     The grid.
//! \param[in]  gridSize The grid size.
//! \param[in]  origin   The origin.
//! \param      mesh     The output triangle mesh.
//! \param[in]  isoValue The iso-surface value.
//! \param[in]  bndFlag  The boundary direction flag.
//!
void marchingCubes(
    const Omni::CellCenteredScalarGrid3Ptr& grid,
    REAL gridSize,
    const Vector3f& origin,
    TriangleMesh3* mesh,
    FLOAT isoValue = 0,
    int bndFlag = Omni::kDirectionAll);

void marchingCubes(
    const Omni::NarrowbandLevelSetPtr& grid,
    REAL gridSize,
    const Vector3f& origin,
    TriangleMesh3* mesh,
    FLOAT isoValue = 0,
    int bndFlag = Omni::kDirectionAll);

}  // namespace jet

#endif  // INCLUDE_JET_MARCHING_CUBES_H_

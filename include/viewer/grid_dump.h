#pragma once
#include <string>
#include "cell_centered_grid.h"
#include "vertex_centered_grid.h"
#include "octree_omni_grid3.h"
#include "octree_grid_layout3.h"
#include "narrowband_levelset.h"
#include "octree_omni_face_centered_grid3.h"

namespace Omni {

inline std::string getGridFilename(const std::string& foldername, int index) {
    return foldername + '/' + "grid" + std::to_string(index) + ".txt";
}

void dumpSdf3D(const std::string &filename, CellCenteredScalarGrid3Ptr grid);

bool readSdf3D(const std::string &filename, CellCenteredScalarGrid3Ptr grid);

void dumpLevelSetGrid3DBinary(const std::string &filename, CellCenteredScalarGrid3 &grid);
bool readLevelSetGrid3DBinary(const std::string &filename, CellCenteredScalarGrid3 &grid);
void dumpLevelSetGrid3DBinary(const std::string &filename, NarrowbandLevelSet &grid);
bool readLevelSetGrid3DBinary(const std::string &filename, NarrowbandLevelSet &grid);

void dumpOctreeOmniGrid3DBinary(const std::string &filename, Size3 res, OctreeOmniScalarGrid3 &grid);
void dumpOctreeGrid3DBinary(const std::string &filename, const Size3& res, OctreeCellCenteredScalarGrid3 &grid);

void dumpOctreeLayout3dBinary(const std::string &filename, OctreeGridLayout3 &grid);
bool readOctreeLayout3dBinary(const std::string &filename, OctreeGridLayout3 &grid);

void dumpOctreeOmniScalar3dBinary(const std::string &filename, OctreeOmniScalarGrid3 &grid);
bool readOctreeOmniScalar3dBinary(const std::string &filename, OctreeOmniScalarGrid3 &grid);

void dumpOctreeOmniTiltE3dBinary(const std::string &filename, OctreeTiltENodeGrid3 &grid);
bool readOctreeOmniTiltE3dBinary(const std::string &filename, OctreeTiltENodeGrid3 &grid);

void dumpOctreeOmniFaceCentered3dBinary(const std::string &filename, OctreeOmniFaceCenteredGrid3 &grid);
bool readOctreeOmniFaceCentered3dBinary(const std::string &filename, OctreeOmniFaceCenteredGrid3 &grid);

} // end of namespace Omni
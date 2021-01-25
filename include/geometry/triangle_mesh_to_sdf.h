#pragma once
#include "geometry/triangle_mesh3.h"
#include "cell_centered_grid3.h"
#include "general.h"

namespace Omni {
namespace Geometry {
	void triangleMeshToSdf(	const TriangleMesh3& mesh, CellCenteredScalarGrid3Ptr sdf, const unsigned int exactBand = 1);
}
}

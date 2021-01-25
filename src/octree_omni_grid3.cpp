#include "octree_omni_grid3.h"

using namespace Omni;


template<>
void OctreeOmniGrid3<int>::updateGhost() {
    vertexPtr->updateGhost();
    // no need to update cell ghost cause there is no direct reference of octagon cell across level
}
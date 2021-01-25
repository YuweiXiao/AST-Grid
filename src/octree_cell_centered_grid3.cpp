#include "octree_cell_centered_grid3.h"

using namespace Omni;

template<>
void OctreeCellCenteredGrid3<int>::updateGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
        // spdlog::info("ghost cell:{}, idx:{}, {}, flag: {}", l, node.xIdx, node.yIdx, node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE);

        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            set(l, node.xIdx, node.yIdx, node.zIdx, get(l+1, node.xIdx>>1, node.yIdx>>1, node.zIdx>>1));
        } else {
            // skip fine to coarse propagation
        }
    });
}
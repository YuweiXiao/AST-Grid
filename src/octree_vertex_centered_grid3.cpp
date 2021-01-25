#include "octree_vertex_centered_grid3.h"

namespace Omni {

template<>
void OctreeVertexCenteredGrid3<TiltENode>::updateGhost() {
	ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
		if (node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
			if ((node.xIdx & 1) || (node.yIdx & 1) || (node.zIdx & 1)) {
				// skip may be half-tilt
				// set(l, node.xIdx, node.yIdx, node.zIdx, TiltENode(0, false, true, true));// closed tilt
			}
			else {
				set(l, node.xIdx, node.yIdx, node.zIdx, get(l+1, node.xIdx>>1, node.yIdx>>1, node.zIdx>>1));
			}
		}
		else {
			set(l, node.xIdx, node.yIdx, node.zIdx, get(l-1, node.xIdx<<1, node.yIdx<<1, node.zIdx<<1));
		}

	});
}

template<>
void OctreeVertexCenteredGrid3<TiltENodeS>::updateGhost() {
	ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
		if (node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
			if ((node.xIdx & 1) || (node.yIdx & 1) || (node.zIdx & 1)) {
				// skip may be half-tilt
				// set(l, node.xIdx, node.yIdx, node.zIdx, TiltENode(0, false, true, true));// closed tilt
			}
			else {
				set(l, node.xIdx, node.yIdx, node.zIdx, get(l+1, node.xIdx>>1, node.yIdx>>1, node.zIdx>>1));
			}
		}
		else {
			set(l, node.xIdx, node.yIdx, node.zIdx, get(l-1, node.xIdx<<1, node.yIdx<<1, node.zIdx<<1));
		}
	});
}


template<>
void OctreeVertexCenteredGrid3<int>::updateGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
		if (node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
			if (node.xIdx & 1 || node.yIdx & 1 || node.zIdx & 1) {
				// skip
			}
			else {
				set(l, node.xIdx, node.yIdx, node.zIdx, get(l+1, node.xIdx>>1, node.yIdx>>1, node.zIdx>>1));
			}
		}
		else {
			set(l, node.xIdx, node.yIdx, node.zIdx, get(l-1, node.xIdx<<1, node.yIdx<<1, node.zIdx<<1));
		}
	});
}


template<>
void OctreeVertexCenteredGrid3<Vector3f>::updateGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
        // spdlog::info("ghost cell:{}, node.zIdx, idx:{}, {}, flag: {}", l, node.xIdx, node.yIdx, node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE);
        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            Vector3f v(0, 0, 0);
            int count = 0;
            for(int dx = 0; dx <= (node.xIdx&1); ++dx) {
                for(int dy = 0; dy <= (node.yIdx&1); ++dy) {
                    for(int dz = 0; dz <= (node.zIdx&1); ++dz) {
                        LCOOR_T c = layout->getUntil(l+1, (node.xIdx>>1)+dx, (node.yIdx>>1)+dy, (node.zIdx>>1)+dz, true);
                        if(c>=0) {
                            v += get(c);
                            count += 1;
                        }
                    }
                }
            }
            set(l, node.xIdx, node.yIdx, node.zIdx, v * (1.0/(count+EPSILON)));
        } else {
            set(l, node.xIdx, node.yIdx, node.zIdx, get(layout->getUntil(l-1, node.xIdx<<1, node.yIdx<<1, node.zIdx<<1, true)));
        }
    });
}


} // end of namespace Omni
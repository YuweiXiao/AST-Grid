#include "octree_face_centered_grid3.h"

using namespace Omni;

OctreeFaceCenteredGrid3::OctreeFaceCenteredGrid3(const Size3& res, REAL spacing, OctreeGridLayout3Ptr _layout, REAL init) 
    : BaseOctreeGrid3(res, spacing, _layout)
{
    ux.resize(layout->level(), std::vector<REAL>());
    for(int i = 0; i < layout->level(); ++i) {
        ux[i].resize(layout->levelSize(i), init);
    }

    uy.resize(layout->level(), std::vector<REAL>());
    for(int i = 0; i < layout->level(); ++i) {
        uy[i].resize(layout->levelSize(i), init);
    }
    
    uz.resize(layout->level(), std::vector<REAL>());
    for(int i = 0; i < layout->level(); ++i) {
        uz[i].resize(layout->levelSize(i), init);
    }
}

Omni::OctreeFaceCenteredGrid3::OctreeFaceCenteredGrid3(OctreeFaceCenteredGrid3 & old)
	: BaseOctreeGrid3(old.resolution(), old.gridSpacing(), old.getLayout())
{
	ux.resize(layout->level(), std::vector<REAL>());
	copy(old.getUx().begin(), old.getUx().end(), ux.begin());
	uy.resize(layout->level(), std::vector<REAL>());
	copy(old.getUy().begin(), old.getUy().end(), uy.begin());
	uz.resize(layout->level(), std::vector<REAL>());
	copy(old.getUz().begin(), old.getUz().end(), uz.begin());
}


Vector3f OctreeFaceCenteredGrid3::sample(REAL x_pos, REAL y_pos, REAL z_pos) const {
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    int level = UNPACK_LEVEL3(levelCoor);
    REAL delta_h = layout->levelGridSpacing(level);
    FLOAT rx = x_pos / delta_h;
    FLOAT ry = y_pos / delta_h;
    FLOAT rz = z_pos / delta_h;

    return Vector3f(sampleX(level, rx, ry, rz), sampleY(level, rx, ry, rz), sampleZ(level, rx, ry, rz));
}

REAL OctreeFaceCenteredGrid3::sampleX(int level, REAL rx, REAL ry, REAL rz) const {
    int xIdx, yIdx, zIdx;
    FLOAT fx, fy, fz;

    getBarycentric(rx, 0, layout->levelResolution(level).x() - 1, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, layout->levelResolution(level).y() - 2, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, layout->levelResolution(level).z() - 2, &zIdx, &fz);

    int xIdx1 = std::min(xIdx + 1, layout->levelResolution(level).x()-1);
    int yIdx1 = std::min(yIdx + 1, layout->levelResolution(level).y()-2);
    int zIdx1 = std::min(zIdx + 1, layout->levelResolution(level).z()-2);

    return trilerp(getX(level, xIdx, yIdx, zIdx),
                    getX(level, xIdx1, yIdx, zIdx),
                    getX(level, xIdx, yIdx1, zIdx),
                    getX(level, xIdx1, yIdx1, zIdx),
                    getX(level, xIdx, yIdx, zIdx1),
                    getX(level, xIdx1, yIdx, zIdx1),
                    getX(level, xIdx, yIdx1, zIdx1),
                    getX(level, xIdx1, yIdx1, zIdx1), fx, fy, fz);
}


REAL OctreeFaceCenteredGrid3::sampleY(int level, REAL rx, REAL ry, REAL rz) const {
    int xIdx, yIdx, zIdx;
    FLOAT fx, fy, fz;

    getBarycentric(rx-0.5, 0, layout->levelResolution(level).x() - 2, &xIdx, &fx);
    getBarycentric(ry, 0, layout->levelResolution(level).y()-1, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, layout->levelResolution(level).z() - 2, &zIdx, &fz);

    int xIdx1 = std::min(xIdx + 1, layout->levelResolution(level).x()-2);
    int yIdx1 = std::min(yIdx + 1, layout->levelResolution(level).y()-1);
    int zIdx1 = std::min(zIdx + 1, layout->levelResolution(level).z()-2);

    return trilerp(getY(level, xIdx, yIdx, zIdx),
                    getY(level, xIdx1, yIdx, zIdx),
                    getY(level, xIdx, yIdx1, zIdx),
                    getY(level, xIdx1, yIdx1, zIdx),
                    getY(level, xIdx, yIdx, zIdx1),
                    getY(level, xIdx1, yIdx, zIdx1),
                    getY(level, xIdx, yIdx1, zIdx1),
                    getY(level, xIdx1, yIdx1, zIdx1), fx, fy, fz);
}


REAL OctreeFaceCenteredGrid3::sampleZ(int level, REAL rx, REAL ry, REAL rz) const {
    int xIdx, yIdx, zIdx;
    FLOAT fx, fy, fz;

    getBarycentric(rx-0.5, 0, layout->levelResolution(level).x() - 2, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, layout->levelResolution(level).y() - 2, &yIdx, &fy);
    getBarycentric(rz, 0, layout->levelResolution(level).z()-1, &zIdx, &fz);

    int xIdx1 = std::min(xIdx + 1, layout->levelResolution(level).x()-2);
    int yIdx1 = std::min(yIdx + 1, layout->levelResolution(level).y()-2);
    int zIdx1 = std::min(zIdx + 1, layout->levelResolution(level).z()-1);

    return trilerp(getZ(level, xIdx, yIdx, zIdx),
                    getZ(level, xIdx1, yIdx, zIdx),
                    getZ(level, xIdx, yIdx1, zIdx),
                    getZ(level, xIdx1, yIdx1, zIdx),
                    getZ(level, xIdx, yIdx, zIdx1),
                    getZ(level, xIdx1, yIdx, zIdx1),
                    getZ(level, xIdx, yIdx1, zIdx1),
                    getZ(level, xIdx1, yIdx1, zIdx1), fx, fy, fz);
}


void OctreeFaceCenteredGrid3::setAllGhost(REAL v) {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        setX(iter.levelCoor(), v);
        setY(iter.levelCoor(), v);
        setZ(iter.levelCoor(), v);
    });
}


void OctreeFaceCenteredGrid3::updateGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
        // spdlog::info("ghost cell:{}, idx:{}, {}, flag: {}", l, node.xIdx, node.yIdx, node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE);
        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            Size3 cIdx(node.xIdx, node.yIdx, node.zIdx);
            for(int axis = 0; axis < 3; ++axis) {
                REAL v = 0; 
                int count = 0;
                for(int i = 0; i <= (cIdx[axis] & 1); ++i) {
                    Size3 idx(node.xIdx>>1, node.yIdx>>1, node.zIdx>>1);
                    idx[axis] += i;
                    LCOOR_T c = layout->getUntil(l+1, idx.x(), idx.y(), idx.z(), true);
                    if(c >= 0) {    // TODO: may have inconsistent issue.
                        v += get(axis, c);
                        count += 1;
                    }
                }
                set(axis, l, node.xIdx, node.yIdx, node.zIdx, v / (REAL)count);
            }
        } else {
            for(int axis = 0; axis < 3; ++axis) {
                REAL v = 0;
                int count = 0;
                for(int i = 0; i <= 1; ++i) {
                    for(int j = 0; j <= 1; ++j) {
                        Size3 idx(node.xIdx<<1, node.yIdx<<1, node.zIdx<<1);
                        idx[(axis+1)%3] += i;
                        idx[(axis+2)%3] += j;
                        LCOOR_T c = layout->getUntil(l-1, idx.x(), idx.y(), idx.z(), true);
                        if(c >= 0) {
                            v += get(axis, c);
                            count += 1;
                        }
                    }
                }
                set(axis, l, node.xIdx, node.yIdx, node.zIdx, v / (REAL)count);
            }
        }
    });
}
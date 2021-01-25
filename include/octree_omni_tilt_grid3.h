#pragma once
#include "base_octree_omni_grid.h"
#include "octree_vertex_centered_grid3.h"

namespace Omni {

class OctreeOmniTiltGrid3 : public BaseOctreeOmniGrid3 {
public:
    OctreeOmniTiltGrid3(const Size3& res, REAL spac, OctreeGridLayout3Ptr layout, OctreeTiltENodeGrid3Ptr tiltE) 
        : BaseOctreeOmniGrid3(res, spac, layout, tiltE)
    {
        gridPtr = std::make_shared<OctreeVertexCenteredVector8fGrid3>(res, spac, layout, Vector8f::Zero());
    }
    ~OctreeOmniTiltGrid3() {};

	OctreeOmniTiltGrid3(OctreeOmniTiltGrid3 & old) 
		: BaseOctreeOmniGrid3(old.resolution(), old.gridSpacing(), old.getLayout(), old.getTiltEGridPtr()) 
	{
		gridPtr = std::make_shared<OctreeVertexCenteredVector8fGrid3>(*old.getGridPtr());
	}

    Vector3f positionCenter(LCOOR_T levelCoor) const {return gridPtr->position(levelCoor);}
    Vector3f positionCenter(int level, int x, int y, int z) const {return gridPtr->position(level, x, y, z);}
    Vector3f positionTilt(int level, int x, int y, int z, int dir) const {
        return positionTilt(positionCenter(level, x, y, z), 
                dir, getTiltEGridPtr()->get(level, x, y, z).e);
    }
    Vector3f positionTilt(LCOOR_T levelCoor, int dir) const {
        return positionTilt(levelCoor, dir, getTiltEGridPtr()->get(levelCoor).e);
    }
    Vector3f positionTilt(LCOOR_T levelCoor, int dir, REAL e) const {
        auto center = positionCenter(levelCoor);
        return positionTilt(center, dir, e);
    }
    Vector3f positionTilt(const Vector3f& center, int dir, REAL e) const {
        int dirm4 = dir & 0x3;
        REAL tz = (dir < 4 ? 1 : -1) * e / 3.0; // sign * INV_SQRT3 * INV_SQRT3 * e;
        REAL ty = (dirm4 < 2 ? 1 : -1) * e / 3.0;
        REAL tx = ((dirm4 == 0 || dirm4 == 3) ? 1 : -1) * e / 3.0;

        return center + Vector3f(tx, ty, tz);
    }

    // bool isValid(int level, int xIdx, int yIdx) const { return Valid(Size3(xIdx, yIdx), layout->levelResolution(level));}
    const Vector8f& get(LCOOR_T levelCoor) const {return gridPtr->get(levelCoor);}
    const Vector8f& get(int level, int x, int y, int z) const {return gridPtr->get(level, x, y, z);}
    const Vector8f& get(int level, const Size3& idx) const {return gridPtr->get(level, idx.x(), idx.y(), idx.z());}
    void set(LCOOR_T levelCoor, const Vector8f& v) {gridPtr->set(levelCoor, v);}
    void set(int level, int x, int y, int z, const Vector8f& v) {gridPtr->set(level, x, y, z, v);}
    void set(LCOOR_T levelCoor, int k, REAL v) {
        auto t = gridPtr->get(levelCoor); 
        t(k) = v;
        gridPtr->set(levelCoor, t);
    }

    virtual void updateGhost() override {
        // seems same with tiltE
        ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
            if(getTiltEGridPtr()->get(iter.levelCoor()).is_half_tilt) {
                return;
            }
			if (iter.node().flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
				if (iter.xIdx() & 1 || iter.yIdx() & 1 || iter.zIdx() & 1) {
					set(iter.levelCoor(), Vector8f::Zero());   // closed tilt
				}
				else {
					set(iter.levelCoor(), get(iter.level + 1, iter.xIdx() >> 1, iter.yIdx() >> 1, iter.zIdx() >> 1));
				}
			}
			else {
				set(iter.levelCoor(), get(iter.level - 1, iter.xIdx() << 1, iter.yIdx() << 1, iter.zIdx() << 1));
			}
		});
    }

	void setAll(const Vector8f& v) { gridPtr->setAll(v); }

	OctreeVertexCenteredVector8fGrid3Ptr & getGridPtr() { return gridPtr; }

private:
    OctreeVertexCenteredVector8fGrid3Ptr gridPtr;
};

typedef std::shared_ptr<OctreeOmniTiltGrid3> OctreeOmniTiltGrid3Ptr;

} // end of namespace Omni
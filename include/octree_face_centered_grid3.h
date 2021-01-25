#pragma once
#include "base_octree_grid.h"

namespace Omni {

class OctreeFaceCenteredGrid3 : public BaseOctreeGrid3 {
public:
    OctreeFaceCenteredGrid3(const Size3& res, REAL spacing, OctreeGridLayout3Ptr layout, REAL init=0);
    ~OctreeFaceCenteredGrid3() {}

	OctreeFaceCenteredGrid3(OctreeFaceCenteredGrid3 &old);

    // sample using absolute positon
    Vector3f sample(const Vector3f& pos) const { return sample(pos.x(), pos.y(), pos.z()); }
    Vector3f sample(REAL x_pos, REAL y_pos, REAL z_pos) const;

    Vector3f positionX(LCOOR_T levelCoor) const {return layout->position(levelCoor) + Vector3f(0.0, 0.5, 0.5) * layout->levelGridSpacing(UNPACK_LEVEL3(levelCoor));}
    Vector3f positionX(int level, int xIdx, int yIdx, int zIdx) const {return layout->position(level, xIdx, yIdx, zIdx) + Vector3f(0.0, 0.5, 0.5) * layout->levelGridSpacing(level);}
    Vector3f positionY(LCOOR_T levelCoor) const {return layout->position(levelCoor) + Vector3f(0.5, 0.0, 0.5) * layout->levelGridSpacing(UNPACK_LEVEL3(levelCoor));}
    Vector3f positionY(int level, int xIdx, int yIdx, int zIdx) const {return layout->position(level, xIdx, yIdx, zIdx) + Vector3f(0.5, 0.0, 0.5) * layout->levelGridSpacing(level);}
    Vector3f positionZ(LCOOR_T levelCoor) const {return layout->position(levelCoor) + Vector3f(0.5, 0.5, 0) * layout->levelGridSpacing(UNPACK_LEVEL3(levelCoor));}
    Vector3f positionZ(int level, int xIdx, int yIdx, int zIdx) const {return layout->position(level, xIdx, yIdx, zIdx) + Vector3f(0.5, 0.5, 0.0) * layout->levelGridSpacing(level);}
    Vector3f position(int axis, int level, int xIdx, int yIdx, int zIdx) const {
        if(axis == 0) return positionX(level, xIdx, yIdx, zIdx);
        else if(axis == 1) return positionY(level, xIdx, yIdx, zIdx);
        else return positionZ(level, xIdx, yIdx, zIdx);
    }

    bool isValidX(int level, int xIdx, int yIdx, int zIdx) const { return Valid(Size3(xIdx, yIdx, zIdx), layout->levelResolution(level) - Size3(0, 1, 1));}
    bool isValidY(int level, int xIdx, int yIdx, int zIdx) const { return Valid(Size3(xIdx, yIdx, zIdx), layout->levelResolution(level) - Size3(1, 0, 1));}
    bool isValidZ(int level, int xIdx, int yIdx, int zIdx) const { return Valid(Size3(xIdx, yIdx, zIdx), layout->levelResolution(level) - Size3(1, 1, 0));}
    bool isValid(int axis, int level, const Size3& idx) const { 
        if(axis == 0) return isValidX(level, idx.x(), idx.y(), idx.z());
        else if(axis == 1) return isValidY(level, idx.x(), idx.y(), idx.z());
        else return isValidZ(level, idx.x(), idx.y(), idx.z());
    }
    bool isValid(int axis, int level, int xIdx, int yIdx, int zIdx) const { 
        if(axis == 0) return isValidX(level, xIdx, yIdx, zIdx);
        else if(axis == 1) return isValidY(level, xIdx, yIdx, zIdx);
        else return isValidZ(level, xIdx, yIdx, zIdx);
    }
    // axis 0: x, 1: y, 2: z
    REAL get(int axis, LCOOR_T levelCoor) const {
        if(axis == 0) return getX(levelCoor);
        else if(axis == 1) return getY(levelCoor);
        else return getZ(levelCoor);
    }
    REAL get(int axis, int level, const Size3& idx) const {
        if(axis == 0) return getX(level, idx.x(), idx.y(), idx.z());
        else if(axis == 1) return getY(level, idx.x(), idx.y(), idx.z());
        else return getZ(level, idx.x(), idx.y(), idx.z());
    }
    REAL get(int axis, int level, int x, int y, int z) const {
        if(axis == 0) return getX(level, x, y, z);
        else if(axis == 1) return getY(level, x, y, z);
        else return getZ(level, x, y, z);
    }
    REAL getX(LCOOR_T levelCoor) const {return getXByIndex(UNPACK_LEVEL3(levelCoor), layout->get(levelCoor)); }
    REAL getX(int level, int x, int y, int z) const {return getXByIndex(level, layout->get(level, x, y, z)); }
    REAL getY(LCOOR_T levelCoor) const {return getYByIndex(UNPACK_LEVEL3(levelCoor), layout->get(levelCoor)); }
    REAL getY(int level, int x, int y, int z) const {return getYByIndex(level, layout->get(level, x, y, z)); }
    REAL getZ(LCOOR_T levelCoor) const {return getZByIndex(UNPACK_LEVEL3(levelCoor), layout->get(levelCoor)); }
    REAL getZ(int level, int x, int y, int z) const {return getZByIndex(level, layout->get(level, x, y, z)); }
    // axis 0: x, 1: y, 2: z
    void set(int axis, LCOOR_T levelCoor, REAL v) { 
        if(axis == 0) setX(levelCoor, v);
        else if(axis == 1) setY(levelCoor, v);
        else setZ(levelCoor, v);
    }
    void set(int axis, int level, int x, int y, int z, REAL v) { 
        if(axis == 0) setX(level, x, y, z, v); 
        else if(axis == 1) setY(level, x, y, z, v);
        else setZ(level, x, y, z, v);
    }
    void setX(LCOOR_T levelCoor, REAL v) {setX(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor), v);}
    void setX(int level, int x, int y, int z, REAL v) {
        auto index = layout->get(level, x, y, z);
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), ux[level].size())), "OctreeFaceCenteredGrid3::setX out of bound");
        ux[level][index] = v;
    }
    bool trySetX(int level, int x, int y, int z, REAL v) {
        auto index = layout->get(level, x, y, z);
        if(index >= 0) {
            ux[level][index] = v;
            return true;
        }
        return false;
    }
    void setY(LCOOR_T levelCoor, REAL v) {setY(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor), v);}
    void setY(int level, int x, int y, int z, REAL v) {
        auto index = layout->get(level, x, y, z);
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), uy[level].size())), "OctreeFaceCenteredGrid3::setY out of bound.");
        uy[level][index] = v;
    }
    bool trySetY(int level, int x, int y, int z, REAL v) {
        auto index = layout->get(level, x, y, z);
        if(index >= 0) {
            uy[level][index] = v;
            return true;
        }
        return false;
    }
    void setZ(LCOOR_T levelCoor, REAL v) {setZ(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor), v);}
    void setZ(int level, int x, int y, int z, REAL v) {
        auto index = layout->get(level, x, y, z);
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), uz[level].size())), "OctreeFaceCenteredGrid3::setZ out of bound.");
        uz[level][index] = v;
    }
    bool trySetZ(int level, int x, int y, int z, REAL v) {
        auto index = layout->get(level, x, y, z);
        if(index >= 0) {
            uz[level][index] = v;
            return true;
        }
        return false;
    }

    void setAllGhost(REAL v);
    void updateGhost() override;
	void setAll(REAL v) {
		for (int i = 0; i < layout->level(); ++i) {
			std::fill(ux[i].begin(), ux[i].end(), v);
			std::fill(uy[i].begin(), uy[i].end(), v);
			std::fill(uz[i].begin(), uz[i].end(), v);
		}
	}

	std::vector<std::vector<REAL>> &getU(int axis) { 
        if (axis == 0) return ux;
        else if (axis == 1) return uy;
        else {
            ASSERT(axis == 2, "non standard axis value");
            return uz;
        }
    }
	std::vector<std::vector<REAL>> &getUx() { return ux; }
	std::vector<std::vector<REAL>> &getUy() { return uy; }
	std::vector<std::vector<REAL>> &getUz() { return uz; }

private:

    REAL getXByIndex(int level, int index) const {
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), ux[level].size())), "OctreeFaceCenteredGrid3::getX out of bound.");
        return ux[level][index];
    }
    REAL getYByIndex(int level, int index) const {
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), uy[level].size())), "OctreeFaceCenteredGrid3::getY out of bound.");
        return uy[level][index];
    }
    REAL getZByIndex(int level, int index) const {
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), uz[level].size())), "OctreeFaceCenteredGrid3::getZ out of bound.");
        return uz[level][index];
    }

    REAL sampleX(int level, REAL rx, REAL ry, REAL rz) const;
    REAL sampleY(int level, REAL rx, REAL ry, REAL rz) const;
    REAL sampleZ(int level, REAL rx, REAL ry, REAL rz) const;

    std::vector<std::vector<REAL>> ux;
    std::vector<std::vector<REAL>> uy;
    std::vector<std::vector<REAL>> uz;
};

typedef std::shared_ptr<OctreeFaceCenteredGrid3> OctreeFaceCenteredGrid3Ptr;

} // end of namespace Omni
#pragma once
#include <unordered_map>
#include "flat_hash_map/flat_hash_map.hpp"
#include "octree_layer.h"
#include "constants.h"
#include "macro.h"
#include "cell_centered_grid.h"

namespace Omni {

struct LayerNode3 {
    size_t index;
    u_int16_t xIdx;
    u_int16_t yIdx;
    u_int16_t zIdx;
    uint32_t flag;  // refer to LayerNodeFlagMask

    LayerNode3(size_t idx=0, u_int16_t x=0, u_int16_t y=0, u_int16_t z=0, uint32_t f=0)
        :index(idx), xIdx(x), yIdx(y), zIdx(z), flag(f) {}
};

typedef std::function<void(int l, const LayerNode3&)> LayerNode3LoopFunc;

class OctreeLayer3 {
public:
    OctreeLayer3(const Size3& res, REAL spacing, int l, CellCenteredFlagGrid3Ptr& region, bool singleLayer=false);

    Vector3f position(LCOOR_T coor) { return position(UNPACK_X3(coor), UNPACK_Y3(coor), UNPACK_Z3(coor));}
    Vector3f position(int xIdx, int yIdx, int zIdx) {return Vector3f(xIdx, yIdx, zIdx) * delta_h;}

    // coordinatie -> index, include ghost cell
    int get(LCOOR_T coordinate) const {
        auto iter = mapTable.find(coordinate);
        return iter == mapTable.end() ? -1 : iter->second.index;
    }
    int size() const {return mapTable.size();}  // include ghost node size
    int ghostNodeSize() const {return ghostNodeIdx.size();}
    int nodeSize() const {return nodeIdx.size();}
    LayerNode3& ghostNode(int idx) {return mapTable[ghostNodeIdx[idx]];}
    LayerNode3& nodeByIdx(int idx) {return mapTable[nodeIdx[idx]];}
    LayerNode3& node(LCOOR_T coor) {
        DEBUG_ONLY(if(get(coor) == -1){throw std::runtime_error("layer get node fail");}); 
        return mapTable[coor];
    }
    bool isExist(LCOOR_T coor, bool countGhost) {
        auto iter = mapTable.find(coor);
        if(iter == mapTable.end())
            return false;
        if(!countGhost && (iter->second.flag & LayerNodeFlagMask::IS_GHOST)) {
            return false;
        }
        return true;
    }
    // 0: do not exist, 1: is ghost, 2: normal cell
    int cellType(LCOOR_T coor) {
        auto iter = mapTable.find(coor);
        if(iter == mapTable.end())
            return 0;
        return (iter->second.flag & LayerNodeFlagMask::IS_GHOST) ? 1 : 2;
    }

    Real gridSpacing() const {return delta_h;}
    const Size3& resolution() const {return res;}

	ska::flat_hash_map<LCOOR_T, LayerNode3>& getMapTable() { return mapTable; }
	std::vector<LCOOR_T>& getGhostNodeIdx() { return ghostNodeIdx; }
	std::vector<LCOOR_T>& getNodeIdx() { return nodeIdx; }

private:
    int level;
    // current layer origin & resolution
    Size3 res;
    REAL delta_h;
    bool isSingleLayer;
    // coordinate -> index
    // std::unordered_map<LCOOR_T, LayerNode3> mapTable;
    ska::flat_hash_map<LCOOR_T, LayerNode3> mapTable;
    std::vector<LCOOR_T> ghostNodeIdx;
    std::vector<LCOOR_T> nodeIdx;
};

typedef std::shared_ptr<OctreeLayer3> OctreeLayer3Ptr;

} // end of namespace Omni
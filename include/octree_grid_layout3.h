#pragma once
#include "grid.h"
#include "macro.h"
#include "octree_layer3.h"

namespace Omni {

class OctreeGridLayout3: public Grid3 {
public:
    OctreeGridLayout3(const Size3& resolution, REAL spacing, int totalLevel, CellCenteredFlagGrid3Ptr splitMap, bool singleLayer=false);

    // return level coordinate of valid and non-ghost cell
    // start searching from finest level
    LCOOR_T levelCoordinate(const Vector3f& pos) const {return levelCoordinate(pos.x(), pos.y(), pos.z());}
    LCOOR_T levelCoordinate(REAL pos_x, REAL pos_y, REAL pos_z) const;

    REAL levelGridSpacing(int level) {
        DEBUG_ONLY(if(level < 0 || level >= totalLevel){throw std::runtime_error("levelGridSpacing::level out of range");}); 
        return layers[level]->gridSpacing();
    }
    const Size3& levelResolution(int level) const {return layers[level]->resolution();}
    int levelSize(int level, bool includeGhost=true) const {return includeGhost ? layers[level]->size() : layers[level]->nodeSize();}
    int levelGhostSize(int level) const {return layers[level]->ghostNodeSize();}
    int totalSize(bool includeGhost = true) const {
        int sum = 0;
        for (int i = 0; i < totalLevel; i++) {
            sum += includeGhost ? layers[i]->size() : layers[i]->nodeSize();
        }
        return sum;
    }
    int totalGhostSize() const {
        int sum = 0;
        for (int i = 0; i < totalLevel; i++) {
            sum += levelGhostSize(i);
        }
        return sum;
    }
    int level() {return totalLevel;}
    int size() const;
    LayerNode3& layerNodeByIndex(int level, int index) {return layers[level]->nodeByIdx(index);}
    LayerNode3& layerGhostNodeByIndex(int level, int index) {return layers[level]->ghostNode(index);}
    LayerNode3& node(int level, int xIdx, int yIdx, int zIdx) {return layers[level]->node(PACK_COOR3(xIdx, yIdx, zIdx));}

    // return index, if invalid return -1
    int get(LCOOR_T levelCoor) const {return get(UNPACK_LEVEL3(levelCoor), UNPACK_COOR3(levelCoor));}
    int get(int level, int x, int y, int z) const {return get(level, PACK_COOR3(x, y, z));}
    int get(int level, LCOOR_T coor) const {
        DEBUG_ONLY(if(level<0 || level>=layers.size()) { spdlog::info("level: {}, coor:{}", level, coor); throw std::runtime_error("OctreeLayout::level out of range");});
        return layers[level]->get(coor);
    }
    bool isExist(int level, int x, int y, int z, bool countGhost=false) { return layers[level]->isExist(PACK_COOR3(x, y, z), countGhost);}
    bool isExist(LCOOR_T levelCoor, bool countGhost=false) { return layers[UNPACK_LEVEL3(levelCoor)]->isExist(UNPACK_COOR3(levelCoor), countGhost);}
    // return levelCoor, include ghost cell
    LCOOR_T getUntil(int l, int x, int y, int z, bool toFinest);
    bool isValid(int l, int x, int y, int z) const {return Valid(Size3(x, y, z), levelResolution(l));}

    Vector3f position(LCOOR_T levelCoor) const {return position(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor));}
    Vector3f position(int level, LCOOR_T coor) const {return position(level, UNPACK_X3(coor), UNPACK_Y3(coor), UNPACK_Z3(coor));}
    Vector3f position(int level, int x, int y, int z) const {return layers[level]->position(x, y, z);}

    std::vector<LCOOR_T>& getHalfGhostIdx() {return halfGhostIdx;};
    std::vector<LCOOR_T>& getTJointArr() {return TJointIdxArr;};
    std::vector<LCOOR_T>& getclosedByHalfTiltArr() {return closedByHalfTiltArr;};

    void forEachGhost(const LayerNode3LoopFunc &func) {
        for(int i = 0; i < totalLevel; ++i) {
            for(int j = 0; j < layers[i]->ghostNodeSize(); ++j) {
                DEBUG_ONLY(if(!(layers[i]->ghostNode(j).flag & LayerNodeFlagMask::IS_GHOST)) {
                    spdlog::error("{}, {}, {}", layers[i]->ghostNode(j).xIdx, layers[i]->ghostNode(j).yIdx, layers[i]->ghostNode(j).zIdx);
                    throw std::runtime_error("OctreeLoopGhost exists a non-ghost");}
                )
                func(i, layers[i]->ghostNode(j));
            }
        }
    }

    OctreeLayer3Ptr getLayerPtr(int level) { return layers[level]; }

private:
    int totalLevel;
    bool isSingleLayer;
    std::vector<OctreeLayer3Ptr> layers;
    std::vector<LCOOR_T> halfGhostIdx;
    std::vector<LCOOR_T> TJointIdxArr;
    std::vector<LCOOR_T> closedByHalfTiltArr;
};


typedef std::shared_ptr<OctreeGridLayout3> OctreeGridLayout3Ptr;

} // end of namespace Omni
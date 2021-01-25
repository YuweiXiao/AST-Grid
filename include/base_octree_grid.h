#pragma once
#include "general.h"
#include "grid.h"
#include "octree_grid_layout3.h"
#include <taskflow/taskflow.hpp>

// #define SINGLE_LEVEL

namespace Omni {

class OctreeGridIterator3 {
public:
    OctreeGridIterator3(OctreeGridLayout3Ptr& _layout)
        :layout(_layout), internalIdxUpperBound(-1)
    {
        reset();
    }
    // must followed by isValid check
    void next() {
        idx += 1;
        internalIdx += 1;
        while(level < layout->level() && layout->levelSize(level, false) <= idx) {
            level += 1;
            idx = 0;
        }
    }
    bool isValid() {
        if(doInternalIndexCheck && internalIdx >= internalIdxUpperBound) {
            return false;
        }
        if(level >=0 && level < layout->level() && idx >= 0 && idx < layout->levelSize(level)) {
            const LayerNode3& n = node();
            x = n.xIdx; y = n.yIdx; z = n.zIdx;
            return true;
        }
        return false;
    }
    void reset() {
        level = 0; idx = 0, internalIdx = 0;
        while(layout->levelSize(level) == 0 && level < layout->level()) {
            level +=1;
        }
    }
    int index() {return idx;}
    LCOOR_T levelCoor() const {return PACK_LCOOR3(level, xIdx(), yIdx(), zIdx());}
    LayerNode3& node() const { return layout->layerNodeByIndex(level, idx);}
    int totalSize() const { return layout->totalSize(false); }
    int xIdx() const {return x;}
    int yIdx() const {return y;}
    int zIdx() const {return z;}
    void setUpperBound(int upper) {
        doInternalIndexCheck = true;
        internalIdx = 0;
        internalIdxUpperBound = upper;
    }

    int level;
protected:
    OctreeGridLayout3Ptr layout;
    int idx, x, y, z;
    bool doInternalIndexCheck = false;
    int internalIdx, internalIdxUpperBound;
};

class OctreeGhostCellIterator3 : public GridIterator {
public:
    OctreeGhostCellIterator3(OctreeGridLayout3Ptr& _layout)
        :layout(_layout)
    {
        reset();
    }
    // must followed by isValid check
    void next() override {
        idx += 1;
        while(level < layout->level() && layout->levelGhostSize(level) <= idx) {
            level += 1;
            idx = 0;
        }
    }
    bool isValid() override {
        if(level >=0 && level < layout->level() 
            && idx >= 0 && idx < layout->levelGhostSize(level)) {
            return true;
        }
        return false;
    }
    void reset() override {
        level = 0; idx = 0;
        while(level < layout->level() && layout->levelGhostSize(level) == 0) {
            level +=1;
        }
    }
    int index() override {return idx;}
    LCOOR_T levelCoor() const {return PACK_LCOOR3(level, xIdx(), yIdx(), zIdx());}
    LayerNode3& node() const { return layout->layerGhostNodeByIndex(level, idx);}
    virtual int totalSize() const { return layout->totalGhostSize(); }
    virtual int xIdx() const {return node().xIdx;}
    virtual int yIdx() const {return node().yIdx;}
    virtual int zIdx() const {return node().zIdx;}
    int level;
protected:
    OctreeGridLayout3Ptr layout;
    int idx;
};


class BaseOctreeGrid3: public Grid3 {
public:
    BaseOctreeGrid3(const Size3& res, REAL spac, const OctreeGridLayout3Ptr& _layout)
        : Grid3(res, spac), layout(_layout)
    {}
    virtual ~BaseOctreeGrid3() {}
    
    virtual void updateGhost() {spdlog::warn("Octree3D::updateGhost::not implemented yet");}
    virtual void setLayout(OctreeGridLayout3Ptr _layout) {layout = _layout;}
    OctreeGridLayout3Ptr& getLayout() {return layout;}
    const OctreeGridLayout3Ptr& getLayout() const {return layout;}
    OctreeGridIterator3 getIterator() {return OctreeGridIterator3(layout);}
    std::vector<OctreeGridIterator3> getParallelIteratorVec() { // should marked as const
        if(cachedParallelIteratorVec.size() == 0) {
            cachedParallelIteratorVec = ThreadPool::initParallelIterator(getIterator());
        }
        return cachedParallelIteratorVec;
    }
    OctreeGhostCellIterator3 getGhostIterator() {return OctreeGhostCellIterator3(layout);}

protected:
    OctreeGridLayout3Ptr layout;
    std::vector<OctreeGridIterator3> cachedParallelIteratorVec;
};

struct LCoor3 {
    int level, xIdx, yIdx, zIdx;
    LCoor3(int l=0, int x=0, int y=0, int z=0) :level(l), xIdx(x), yIdx(y), zIdx(z) {}
    void set(int l, int x, int y, int z) {
        level = l; xIdx = x; yIdx = y; zIdx = z;
    }
    LCOOR_T levelCoor() const { 
        ASSERT(level>=0 && xIdx>=0 && yIdx>=0 && zIdx>=0, "LCoor::packed value is negative");
        return PACK_LCOOR3(level, xIdx, yIdx, zIdx);
    }
};

#define iterate_grid(iter, grid) \
for(auto iter = grid->getIterator(); iter.isValid();iter.next())

#define iterate_half_tilt(lCoor, layout)\
for(auto lCoor: layout->getHalfGhostIdx())

template<typename LayoutPtr, typename Func>
inline void parallel_iterate_half_tilt(const LayoutPtr& layout, const Func& func) { ThreadPool::parallelForTF(layout->getHalfGhostIdx(), func); }

} // end of namespace Omni
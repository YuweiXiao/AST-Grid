#pragma once
#include "general.h"
#include "octree_grid_single_layout3.h"
#include "octree_cell_centered_grid3.h"
#include "jet/fmm_levelset_solver3.h"
#include "bitmap.h"
#include "octree_util.h"
#include "thread_pool.h"
#include "advection.h"
#include "bitmap_cell_centered_grid3.h"

namespace Omni {

class NarrowbandLevelSet {
public:
    NarrowbandLevelSet(const Size3& _res, REAL _delta_h, const std::string& _name="narrowband_level_set") 
        :res(_res), delta_h(_delta_h), name(_name)
    {}
    NarrowbandLevelSet(const Size3& _res, REAL _delta_h, BitMapPtr _inOutMarker, BitMapPtr _narrowbandMarker, const std::string& _name="narrowband_level_set") 
        :res(_res), delta_h(_delta_h), name(_name)
    {
        inOutMarker = _inOutMarker;
        narrowbandMarker = _narrowbandMarker;
        levelSet = make_shared<BitMapCellCenteredScalarGrid3>(res, delta_h, narrowbandMarker, 1000);
    }
    NarrowbandLevelSet(const NarrowbandLevelSet& old)
        :res(old.res), delta_h(old.delta_h), name(old.name)
    {
        inOutMarker = make_shared<BitMap>(*old.inOutMarker);
        narrowbandMarker = make_shared<BitMap>(*old.narrowbandMarker);
        levelSet = make_shared<BitMapCellCenteredScalarGrid3>(res, delta_h, narrowbandMarker, *old.levelSet);
    }

    template<typename SDF>
    void initialize(const SDF& sdf, REAL distance) {
        inOutMarker = make_shared<BitMap>(res, delta_h);
        narrowbandMarker = make_shared<BitMap>(res, delta_h);
        narrowbandMarker->fill(0);
        refineSplitMapByObstacleSdf(narrowbandMarker, sdf, Vector2f(-distance, distance), true);
        levelSet = make_shared<BitMapCellCenteredScalarGrid3>(res, delta_h, narrowbandMarker, 1000);

        ThreadPool::parallelForTF(0, res.z(), [&](int zIdx){
            for(int yIdx = 0; yIdx < res.y(); ++yIdx) {
                for(int xIdx = 0; xIdx < res.x(); ++xIdx) {
                    set(xIdx, yIdx, zIdx, sdf->sample(position(xIdx, yIdx, zIdx)));
                }
            }
        });
    }

    void initialize(const NarrowbandLevelSet& old) {
        inOutMarker = make_shared<BitMap>(*old.inOutMarker);
        narrowbandMarker = make_shared<BitMap>(*old.narrowbandMarker);
        levelSet = make_shared<BitMapCellCenteredScalarGrid3>(res, delta_h, narrowbandMarker, *old.levelSet);
    }

    Vector3f position(int xIdx, int yIdx, int zIdx) {return levelSet->position(xIdx, yIdx, zIdx);}
    void set(int xIdx, int yIdx, int zIdx, REAL v) {
        inOutMarker->set(xIdx, yIdx, zIdx, isInsideSdf(v) ? true : false);
        if(narrowbandMarker->get(xIdx, yIdx, zIdx))
            levelSet->set(xIdx, yIdx, zIdx, v);
    }
    REAL get(const Size3& idx) const {return get(idx.x(), idx.y(), idx.z());}
    REAL get(int xIdx, int yIdx, int zIdx) const {
        if(narrowbandMarker->get(xIdx, yIdx, zIdx)) 
            return levelSet->get(xIdx, yIdx, zIdx);
        REAL sign = inOutMarker->get(xIdx, yIdx, zIdx) ? -1: 1;
        return sign * 1000;
    }
    // NOTE: if pos is out of narrow band, return sign(pos) * 1000
    REAL sample(const Vector3f& pos) const {return sample(pos.x(), pos.y(), pos.z());}
    REAL sample(REAL pos_x, REAL pos_y, REAL pos_z) const;
    bool isExist(int x, int y, int z) {return narrowbandMarker->get(x, y, z);}
    
    void parallelForEach(const LoopFunc3& func) {
        ThreadPool::parallelForTF(0, resolution().z(), [&](int zIdx) {
        for(int y = 0; y < resolution().y(); ++y) {
                for(int x = 0; x < resolution().x(); ++x) {
                    func(x, y, zIdx);
                }   
            }
        });
    }

    int getSign(const Size3& idx) const {return getSign(idx.x(), idx.y(), idx.z());}
    int getSign(int xIdx, int yIdx, int zIdx) const {return inOutMarker->get(xIdx, yIdx, zIdx) ? -1:1;}

    void reinitialize(REAL maxDistance);

    template<typename VT, typename BT>
    void semiLagrangian(const VT& flow, const BT& bSdf, REAL dt);

    const Size3& resolution() const {return res;}
    REAL gridSpacing() const {return delta_h;}
    BitMapPtr getInOutMarker() {return inOutMarker;}
    void setInOutMarker(BitMapPtr v) {inOutMarker = v;}
    BitMapPtr getNarrowbandMarker() {return narrowbandMarker;}
    void setNarrowbandMarker(BitMapPtr v) {narrowbandMarker = v;}
    void setLevelSet(BitMapCellCenteredScalarGrid3Ptr v) {levelSet = v;}

    // SdfFunc accepts a position and return its sdf value
    template<typename SdfFunc>
    static std::shared_ptr<NarrowbandLevelSet> buildLevelsetFromSDF(const Size3& res, double gridSpacing, const SdfFunc& sdfFunc);

private:
    BitMapPtr inOutMarker = nullptr, narrowbandMarker = nullptr;
    BitMapCellCenteredScalarGrid3Ptr levelSet = nullptr;
    std::vector<std::pair<Size3, REAL>> container;
    Size3 res;
    REAL delta_h;
    std::string name;
};

typedef std::shared_ptr<NarrowbandLevelSet> NarrowbandLevelSetPtr;

template<typename VT, typename BT>
void NarrowbandLevelSet::semiLagrangian(const VT& flow, const BT& bSdf, REAL dt) {
    auto tInOutMarker = make_shared<BitMap>(*inOutMarker);
    auto tLevelSet = make_shared<BitMapCellCenteredScalarGrid3>(res, delta_h, narrowbandMarker);
    ThreadPool::parallelForTF(0, res.z(), [&](int z){
        for(int y = 0; y < res.y(); ++y) {
            for(int x = 0; x < res.x(); ++x) {
                if(narrowbandMarker->get(x, y, z)) {
                    Vector3f pos = position(x, y, z);
                    if(bSdf && isInsideSdf(bSdf->sample(pos))) {
                        tLevelSet->set(x, y, z, get(x, y, z));
                        continue;
                    }
                    Vector3f oldPos = BackTraceHighPrecision<OCTREE_LEVEL_RES_FACTOR>(bSdf, flow, pos, dt, gridSpacing());
                    REAL l = sample(oldPos);
                    tLevelSet->set(x, y, z, l);
                    tInOutMarker->set(x, y, z, isInsideSdf(l) ? true : false);
                }
            }
        }
    });
    levelSet = tLevelSet;
    inOutMarker = tInOutMarker;
}


template<typename SdfFunc>
std::shared_ptr<NarrowbandLevelSet> NarrowbandLevelSet::buildLevelsetFromSDF(const Size3& res, double gridSpacing, const SdfFunc& sdfFunc) {
    auto narrowBand = make_shared<BitMap>(res, gridSpacing);
    ThreadPool::parallelForTF(0, res.z(), [&](int z){
        for(int y = 0; y < res.y(); ++y) {
        for(int x = 0; x < res.x(); ++x) {
            auto pos = narrowBand->position(x, y, z);
            if(fabs(sdfFunc(pos)) <= 10 * narrowBand->gridSpacing()) {
                narrowBand->set(x, y, z, 1);
            }
        }}
    });

    auto inOut = make_shared<BitMap>(res, gridSpacing);
    auto narrowbandSdf = make_shared<NarrowbandLevelSet>(res, gridSpacing, inOut, narrowBand);
    ThreadPool::parallelForTF(0, res.z(), [&](int z){
        for(int y = 0; y < res.y(); ++y) {
        for(int x = 0; x < res.x(); ++x) {
            narrowbandSdf->set(x, y, z, sdfFunc(narrowbandSdf->position(x, y, z)));
        }}
    });

    return narrowbandSdf;
}

}   // end of namespace Omni
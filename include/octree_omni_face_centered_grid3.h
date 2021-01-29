#pragma once
#include "base_octree_omni_grid.h"
#include "octree_face_centered_grid3.h"
#include "octree_omni_tilt_grid3.h"
#include "octree_omni_grid3.h"
#include "global_benchmark.h"

namespace Omni {

class OctreeOmniFaceCenteredGrid3 : public BaseOctreeOmniGrid3 {
public:
    OctreeOmniFaceCenteredGrid3(const Size3& res, REAL spac, OctreeGridLayout3Ptr layout, OctreeTiltENodeGrid3Ptr tiltE, bool constructDual = true)
        : BaseOctreeOmniGrid3(res, spac, layout, tiltE) 
    {
        MACGridPtr = std::make_shared<OctreeFaceCenteredGrid3>(res, spac, layout);
        tiltGridPtr = std::make_shared<OctreeOmniTiltGrid3>(res, spac, layout, tiltE);
        if(constructDual) {
            dualGridPtr = std::make_shared<OctreeOmniVector3Grid3>(res, spac, layout, tiltE, Vector3f::Zero());
        } else {
            dualGridPtr = nullptr;
        }
    }
    ~OctreeOmniFaceCenteredGrid3() {}

	OctreeOmniFaceCenteredGrid3(OctreeOmniFaceCenteredGrid3 & old) 
		: BaseOctreeOmniGrid3(old.resolution(), old.gridSpacing(), old.getLayout(), old.getTiltEGridPtr()) 
	{
		tiltE = make_shared<OctreeTiltENodeGrid3>(*old.getTiltEGridPtr());
		MACGridPtr = std::make_shared<OctreeFaceCenteredGrid3>(*old.getMACGrid());
		tiltGridPtr = std::make_shared<OctreeOmniTiltGrid3>(*old.getTiltGrid());
		tiltGridPtr->setTiltEGridPtr(tiltE);
		dualGridPtr = std::make_shared<OctreeOmniVector3Grid3>(*old.getDualGrid());
		dualGridPtr->setTiltEGridPtr(tiltE);
	}

    void setAll(REAL v);
    Vector3f sample(const Vector3f &pos) const {
        // return dualGridPtr->sample(pos);
#ifdef SINGLE_LEVEL
        Vector3f rxyz = pos / gridSpacing() - Vector3f::Ones()*0.5;
        Vector3f f; 
        Size3 idx;
        getBaryCentric(rxyz, resolution()-Size3::Ones(), idx, f);
        Vector3f v = MACGridPtr->sample(pos);
        if(getTiltEGridPtr()->get(0, idx.x()+1, idx.y()+1, idx.z()+1).is_closed)
            return v;
        Vector3f diff = dualGridPtr->getTilt(0, idx.x()+1, idx.y()+1, idx.z()+1) -
            MACGridPtr->sample(getTiltEGridPtr()->position(0, idx.x()+1, idx.y()+1, idx.z()+1));
        return v + diff * 2 * std::min(f.minCoeff(), (Vector3f::Ones()-f).minCoeff());
#elif defined(OCTREE_CORRECTION_BASED)
        LCOOR_T levelCoor = layout->levelCoordinate(pos.x(), pos.y(), pos.z());
        int level = UNPACK_LEVEL3(levelCoor);
        Real delta_h = layout->levelGridSpacing(level);
        const Size3& cRes = layout->levelResolution(level);

        Vector3f rxyz = pos / delta_h - Vector3f::Ones()*0.5;
        Vector3f f; 
        Size3 idx;
        getBaryCentric(rxyz, cRes-Size3::Ones(), idx, f);
        const auto& tNode = getTiltEGridPtr()->get(level, idx.x()+1, idx.y()+1, idx.z()+1);
        if((getLayout()->node(level, idx.x(), idx.y(), idx.z()).flag & LayerNodeFlagMask::IS_GHOST) || 
            (level == 0 && (getLayout()->get(level+1, idx.x()>>1, idx.y()>>1, idx.z()>>1) >= 0)) ||
            (level == 1 && (getLayout()->get(level-1, (idx.x()<<1)+1, (idx.y()<<1)+1, (idx.z()<<1)+1) >= 0)) ||
            (level == 1 && (getLayout()->get(level-1, (idx.x()<<1), (idx.y()<<1), (idx.z()<<1)) >= 0)) ) {
            // return dualGridPtr->sample(pos);
            return dualGridPtr->sampleSub(pos.x(), pos.y(), pos.z());
        }

        Vector3f v = MACGridPtr->sample(pos);
        if(!tNode.is_half_tilt && tNode.is_closed) {
            return v;
        } else {
            Vector3f diff = dualGridPtr->getTilt(level, idx.x()+1, idx.y()+1, idx.z()+1) -
                MACGridPtr->sample(getTiltEGridPtr()->position(level, idx.x()+1, idx.y()+1, idx.z()+1));
            return v + diff * 2 * std::min(f.minCoeff(), (Vector3f::Ones()-f).minCoeff());
        }
#else
        return dualGridPtr->sample(pos);
#endif
    }
    Vector3f sample(REAL x, REAL y, REAL z) const { return sample(Vector3f(x, y, z)); }

    void setTiltEGridPtr(OctreeTiltENodeGrid3Ptr p) override {
        BaseOctreeOmniGrid3::setTiltEGridPtr(p);
        tiltGridPtr->setTiltEGridPtr(p);
        dualGridPtr->setTiltEGridPtr(p);
    } 

    OctreeFaceCenteredGrid3Ptr& getMACGrid() { return MACGridPtr; }
    OctreeOmniTiltGrid3Ptr& getTiltGrid() { return tiltGridPtr; }
    OctreeOmniVector3Grid3Ptr& getDualGrid() { return dualGridPtr; }
    void setDualGrid(OctreeOmniVector3Grid3Ptr ptr) {dualGridPtr = ptr;}
    void interpolateFromDualGraph(OctreeOmniVector3Grid3Ptr ptr);
    void updateDualGrid();
    // template<typename SDF>
    // void updateDualGrid(const SDF& sdf, REAL threshold = 0.1);
    virtual void updateGhost() override {
        BENCHMARK_SCOPED_TIMER_SECTION t("face, update ghost cell");
        
        MACGridPtr->setAllGhost(0);
        ThreadPool::parallelIterateGrid(MACGridPtr->getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
            if(!(iter.node().flag & LayerNodeFlagMask::GHOST_FROM_COARSE)) {
                for(int axis = 0; axis < 3; ++axis) {
                    auto pos = MACGridPtr->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
                    MACGridPtr->set(axis, iter.levelCoor(), MACGridPtr->sample(pos)[axis]);
                }
            }
        });
        updateGhostFromHalfTilt();  // override face value with half-tilted cell
        ThreadPool::parallelIterateGrid(MACGridPtr->getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
            if(iter.node().flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
                for(int axis = 0; axis < 3; ++axis) {
                    auto pos = MACGridPtr->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
                    MACGridPtr->set(axis, iter.levelCoor(), MACGridPtr->sample(pos)[axis]);
                }
            }
        });

        update_ghost_exclude_half_tilt(Vector8f, tiltGridPtr);
    }
    // NOTE: Some non-ghost MAC faces may also be updated
    void updateGhostFromHalfTilt();

private:
    OctreeFaceCenteredGrid3Ptr MACGridPtr;
    OctreeOmniTiltGrid3Ptr tiltGridPtr;
    OctreeOmniVector3Grid3Ptr dualGridPtr;
};

typedef std::shared_ptr<OctreeOmniFaceCenteredGrid3> OctreeOmniFaceCenteredGrid3Ptr;

} // end of namespace Omni
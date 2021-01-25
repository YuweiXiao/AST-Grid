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

// template<typename SDF>
// void OctreeOmniFaceCenteredGrid3::updateDualGrid(const SDF& sdf, REAL threshold) {
//     auto& dualTiltGrid = dualGridPtr->getTiltGrid();
//     auto& dualOctagonGrid = dualGridPtr->getOctagonGrid();
// 	auto& tiltE = getTiltEGridPtr();
    
// 	ThreadPool::parallelIterateGrid(dualTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
//         if(tiltE->get(iter.levelCoor()).is_closed) {	// skip half-tilt for now
//             dualTiltGrid->set(iter.levelCoor(), Vector3f::Zero());
//         } else {
//             auto pos = dualTiltGrid->position(iter.levelCoor());
// 			if(sdf->sample(pos) > threshold)
// 				return;
//             const Vector8f& gVal = tiltGridPtr->get(iter.levelCoor());
//             Vector4f v((gVal(0) + gVal(6)) * 0.5, (gVal(1) +gVal(7)) * 0.5,
//                         (gVal(2) + gVal(4)) * 0.5, (gVal(3) + gVal(5)) * 0.5);
//             Vector3f nv(0, 0, 0);
//             nv(0) = INV_SQRT3 * v(0) - INV_SQRT3 * v(1) - INV_SQRT3 * v(2) + INV_SQRT3 * v(3);
//             nv(1) = INV_SQRT3 * v(0) + INV_SQRT3 * v(1) - INV_SQRT3 * v(2) - INV_SQRT3 * v(3);
//             nv(2) = INV_SQRT3 * v(0) + INV_SQRT3 * v(1) + INV_SQRT3 * v(2) + INV_SQRT3 * v(3);
//             nv *= 3.0/4.0;
// 			ASSERT(!std::isnan(nv.x()) && !std::isnan(nv.y()) && !std::isnan(nv.z()), "dual velocity is nan");
//             dualTiltGrid->set(iter.levelCoor(), nv);
//         }
// 	});


// 	ThreadPool::parallelIterateGrid(dualOctagonGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
//         if(!dualOctagonGrid->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
//             dualOctagonGrid->set(iter.levelCoor(), Vector3f::Zero());
//             return;
//         }

// 		auto pos = dualOctagonGrid->position(iter.levelCoor());
// 		if(sdf->sample(pos) > threshold)
// 			return;


// 		REAL delta_h = layout->levelGridSpacing(iter.level);
// 		REAL sqrt3_delta_h = SQRT3 * delta_h;
// 		REAL half_sqrt3_delta_h = sqrt3_delta_h * 0.5;

// 		int n_joint[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
// 		for (int i = 0; i < 8; i++) {
// 			const auto& node = tiltE->get(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[i](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[i](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[i](2));
// 			n_joint[i] = (node.is_T_joint && fabs(node.e - delta_h) < EPSILON) || node.is_half_tilt;
// 		}

// 		Vector7f v = Vector7f::Zero();
// 		int cnt[7] = { 0, 0, 0, 0, 0, 0, 0 };

// 		// left (x-)
// 		if (!((n_joint[1] && n_joint[6]) || (n_joint[2] && n_joint[5]))) {
// 			v[0] += MACGridPtr->getX(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
// 			cnt[0]++;
// 		}
// 		// right (x+)
// 		if (!((n_joint[0] && n_joint[7]) || (n_joint[3] && n_joint[4]))) {
// 			v[0] += MACGridPtr->getX(iter.level, iter.xIdx() + 1, iter.yIdx(), iter.zIdx());
// 			cnt[0]++;
// 		}
// 		// front (y-)
// 		if (!((n_joint[2] && n_joint[7]) || (n_joint[3] && n_joint[6]))) {
// 			v[1] += MACGridPtr->getY(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
// 			cnt[1]++;
// 		}
// 		// back (y+)
// 		if (!((n_joint[0] && n_joint[5]) || (n_joint[1] && n_joint[4]))) {
// 			v[1] += MACGridPtr->getY(iter.level, iter.xIdx(), iter.yIdx() + 1, iter.zIdx());
// 			cnt[1]++;
// 		}
// 		// down (z-)
// 		if (!((n_joint[4] && n_joint[6]) || (n_joint[5] && n_joint[7]))) {
// 			v[2] += MACGridPtr->getZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
// 			cnt[2]++;
// 		}
// 		// up (z+)
// 		if (!((n_joint[0] && n_joint[2]) || (n_joint[1] && n_joint[3]))) {
// 			v[2] += MACGridPtr->getZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx() + 1);
// 			cnt[2]++;
// 		}

// 		for(int i = 0; i < 4; ++i) {
// 			int slot = i;
// 			int reSlot = ((i&1)?8:6) - i;	// opposite index of slot
// 			int vSlot = i+3;
// 			const auto& node = tiltE->get(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[i](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[i](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[i](2));
// 			const auto& reNode = tiltE->get(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[reSlot](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[reSlot](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[reSlot](2));
// 			if (!node.is_closed || node.is_half_tilt) {
// 				v[vSlot] += tiltGridPtr->get(iter.level, 
// 						iter.xIdx() + OCT_NB_TILT_DXDYDZ[slot].x(), 
// 						iter.yIdx() + OCT_NB_TILT_DXDYDZ[slot].y(), 
// 						iter.zIdx() + OCT_NB_TILT_DXDYDZ[slot].z())(reSlot);
// 				cnt[vSlot]++;
// 			}
// 			if (!reNode.is_closed || reNode.is_half_tilt) {
// 				REAL tv = tiltGridPtr->get(iter.level, 
// 						iter.xIdx() + OCT_NB_TILT_DXDYDZ[reSlot].x(), 
// 						iter.yIdx() + OCT_NB_TILT_DXDYDZ[reSlot].y(), 
// 						iter.zIdx() + OCT_NB_TILT_DXDYDZ[reSlot].z())(slot);
// 				v[vSlot] = cnt[vSlot] == 0? tv : (((v[vSlot] * (half_sqrt3_delta_h - reNode.e * INV_SQRT3)
// 					+ tv * (half_sqrt3_delta_h - node.e * INV_SQRT3))
// 					/ (sqrt3_delta_h - node.e * INV_SQRT3 - reNode.e * INV_SQRT3)));
// 				cnt[vSlot]++;
// 			}
// 		}

// 		for (int i = 0; i < 3; i++) {
// 			v[i] = v[i] / (cnt[i] + EPSILON);
// 		}
// 		Vector3f finalV = Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[cnt[0] > 0][cnt[1] > 0][cnt[2] > 0]
// 						[cnt[3] > 1 || n_joint[0] || n_joint[6]] [cnt[4] > 1 || n_joint[1] || n_joint[7]]
// 						[cnt[5] > 1 || n_joint[2] || n_joint[4]] [cnt[6] > 1 || n_joint[3] || n_joint[5]] * v;

// 		ASSERT(!std::isnan(finalV.x()) && !std::isnan(finalV.y()) && !std::isnan(finalV.z()), "dual velocity is nan");
// 		dualOctagonGrid->set(iter.levelCoor(), finalV);
//     });

// 	iterate_half_tilt(lCoor, layout) {
// 		const auto& tNode = tiltE->get(lCoor);
// 		const auto& gVal = tiltGridPtr->get(lCoor);
// 		Vector7f v = Vector7f::Zero();
// 		int axis = tNode.half_direction % 3;	// x:0, y:1, z:2

// 		v[axis] = gVal(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]);
// 		for(int i = 0; i < 4; ++i) {
// 			int slot = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
// 			int reSlot = ((slot&1)?8:6) - slot;
// 			int idx = slot < reSlot ? slot : reSlot;
// 			v[idx+3] = gVal(slot);
// 		}

// 		Vector3f finalV = Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[axis == 0][axis == 1][axis == 2][1][1][1][1] * v;
// 		ASSERT(!std::isnan(finalV.x()) && !std::isnan(finalV.y()) && !std::isnan(finalV.z()), "dual velocity is nan");
// 		dualTiltGrid->set(lCoor, finalV);
// 	}


// }

typedef std::shared_ptr<OctreeOmniFaceCenteredGrid3> OctreeOmniFaceCenteredGrid3Ptr;

} // end of namespace Omni
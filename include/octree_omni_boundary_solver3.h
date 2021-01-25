#pragma once
#include "general.h"
#include "constants.h"
#include "geometry/collider3.h"
#include "octree_omni_grid3.h"
#include "octree_omni_face_centered_grid3.h"
#include "octree_cell_centered_grid3.h"
#include "octree_vertex_centered_grid3.h"

namespace Omni {

using namespace Geometry;

class OctreeOmniBoundarySolver3 {
public:
    OctreeOmniBoundarySolver3(const Size3& res, REAL spacing, OctreeGridLayout3Ptr& layout, OctreeTiltENodeGrid3Ptr& tiltE);
    ~OctreeOmniBoundarySolver3() {}

    const Collider3Ptr& collider() const {return _collider;}

    void setClosedDomainBoundaryFlag(int f) {_closedDomainBoundaryFlag = f;}
    int closedDomainBoundaryFlag() const {return _closedDomainBoundaryFlag;}
    bool isClosedDomain(int axis = 0);
    // input index(level, x, y, z) is a octagon index
    bool isClosedDomain(int level, int x, int y, int z, OctreeGridLayout3Ptr& layout);
    // Applies new collider and build the internals.
    void updateCollider(const Collider3Ptr& newCollider);
    void updateSdf(CellCenteredScalarGrid3Ptr& sdf);

    // only domain velocity
    void constrainDomainVelocity(OctreeOmniFaceCenteredGrid3Ptr& velocity);
    // dual graph velocity
    void constrainFaceVelocity(OctreeOmniFaceCenteredGrid3Ptr& velocity);
    void constrainDualVelocity(OctreeOmniVector3Grid3Ptr& dualVelocity);
    void constrainVelocity(OctreeOmniFaceCenteredGrid3Ptr& velocity) {
        auto& tMACGrid = velocity->getMACGrid();
        auto& tTiltGrid = velocity->getTiltGrid();
        auto& tiltE = velocity->getTiltEGridPtr();
        auto& layout = velocity->getLayout();

        constrainDomainVelocity(velocity);

        if(noObstacle)
            return;

        ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
            for(int axis = 0; axis < 3; ++axis) {
                auto pos = tMACGrid->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
                if(_colliderSdf->sample(pos) < 0) {
                    tMACGrid->set(axis, iter.levelCoor(), 0);
                }
            }
        });
        ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
            if(!tiltE->get(iter.levelCoor()).is_closed) {
                auto nV = tTiltGrid->get(iter.levelCoor());
                bool flag = false;
                for(int k = 0; k < 8; ++k) {
                    auto pos = tTiltGrid->positionTilt(iter.levelCoor(), k);
                    if(_colliderSdf->sample(pos) < 0) {
                        nV[k] = 0;
                    }
                }
                if(flag)
                    tTiltGrid->set(iter.levelCoor(), nV);
            }
        });
        parallel_iterate_half_tilt(layout, [&](auto lCoor){
            Vector8f newF = tTiltGrid->get(lCoor);
            int level = UNPACK_LEVEL3(lCoor);
            bool flag = false;
            const auto& tNode = tiltE->get(lCoor);
            for(int i = 0; i < 5; ++i) {
                int k = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
                auto pos = (i == 4) ? tTiltGrid->positionCenter(lCoor) : tTiltGrid->positionTilt(lCoor, k);
                if(_colliderSdf->sample(pos) < 0) {
                    newF(k) = 0;
                }
            }
            if(flag)
                tTiltGrid->set(lCoor, newF);
        });
    }

    OctreeOmniScalarGrid3Ptr colliderSdf() const {return _colliderSdf;}
    REAL sample(const Vector3f& pos) const {
        if(noObstacle)  return 1;
        else return colliderSdf()->sample(pos);
    }

    bool noObstacle = true;
    
private:

    Collider3Ptr _collider = nullptr;
	OctreeOmniScalarGrid3Ptr _colliderSdf = nullptr;
    OctreeOmniVector3Grid3Ptr _colliderGradient = nullptr;
    Size3 resolution;
    REAL delta_h;
    int _closedDomainBoundaryFlag = kDirectionNone;
};

typedef std::shared_ptr<OctreeOmniBoundarySolver3> OctreeOmniBoundarySolver3Ptr;

}  // end of namespace Omni

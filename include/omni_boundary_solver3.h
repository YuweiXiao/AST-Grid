#pragma once
#include "general.h"
#include "geometry/collider3.h"
#include "octree_face_centered_grid3.h"
#include "octree_omni_face_centered_grid3.h"
#include "global_benchmark.h"

namespace Omni {

using namespace Geometry;

class OmniBoundarySolver3 final {
public:
    OmniBoundarySolver3(const Size3& res, REAL spacing);
    ~OmniBoundarySolver3() {}

    const Collider3Ptr& collider() const {return _collider;}

    void setClosedDomainBoundaryFlag(int f) {_closedDomainBoundaryFlag = f;}
    int closedDomainBoundaryFlag() const {return _closedDomainBoundaryFlag;}
    // Applies new collider and build the internals.
    void updateCollider(const Collider3Ptr& newCollider);

    // void constrainVelocity(OmniFaceCenteredGrid3Ptr& velocity, const OmniFaceCenteredGrid3Ptr& weightsCC = nullptr);

    void constrainVelocity(OctreeOmniTiltGrid3Ptr& velocity, OctreeFaceCenteredGrid3Ptr& macV) {
        Vector3f sCenter(0.5, 0.5, 0.45);
        ThreadPool::parallelIterateGrid(velocity->getParallelIteratorVec(), [&](const auto& iter){
            if(velocity->getTiltEGridPtr()->get(iter.levelCoor()).is_closed)
                return;
            auto posCenter = velocity->positionCenter(iter.levelCoor());
            Real dx = velocity->getLayout()->levelGridSpacing(iter.level);

            if(isInsideSdf(_colliderSdf->sample(posCenter) - dx)) {
                auto v = velocity->get(iter.levelCoor());
                REAL e = dx * 0.5;
                for(int k = 0; k < 8; ++k) {
                    auto pos1 = velocity->positionTilt(posCenter, k, e);
                    if(isInsideSdf(_colliderSdf->sample(pos1))) {
                        Vector3f n = pos1 - sCenter;
                        if(n.norm() < EPSILON)
                            v[k] = 0;
                        else {
                            auto v = macV->sample(pos1);
                            n.normalize();
                            v(k) = (v - v.dot(n) * n).dot(TILT3_UNIT_DIR_UP[k]);
                        }
                    }
                }
                velocity->set(iter.levelCoor(), v);
            }
        });
    }

    void constrainVelocity(OctreeOmniVector3Grid3Ptr& velocity) {
        auto& tiltGrid = velocity->getTiltGrid();
        Vector3f sCenter(0.5, 0.5, 0.45);
        ThreadPool::parallelIterateGrid(tiltGrid->getParallelIteratorVec(), [&](const auto& iter){
            auto pos = tiltGrid->position(iter.levelCoor());
            if(isInsideSdf(_colliderSdf->sample(pos))) {
                Vector3f n = pos - sCenter;
                if(n.norm() < EPSILON)
                    tiltGrid->set(iter.levelCoor(), Vector3f::Zero());
                else {
                    auto v = tiltGrid->get(iter.levelCoor());
                    n.normalize();
                    tiltGrid->set(iter.levelCoor(), (v - v.dot(n) * n));
                }
            }
        });
    }

    void constrainVelocity(OctreeFaceCenteredGrid3Ptr& velocity) {
        BENCHMARK_SCOPED_TIMER_SECTION t_boundary("constrain velocity");
        auto& layout = velocity->getLayout();

        Vector3f sCenter(0.5, 0.5, 0.45);

        ThreadPool::parallelIterateGrid(velocity->getParallelIteratorVec(), [&](const auto& iter) {
            for(int axis = 0; axis < 3; ++axis) {
                if(velocity->isValid(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
                    auto pos = velocity->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
                    if(isInsideSdf(_colliderSdf->sample(pos))) {
                        Vector3f n = pos - sCenter;
                        if(n.norm() < EPSILON)
                            velocity->set(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx(), 0);
                        else {
                            auto v = velocity->sample(pos);
                            n.normalize();
                            velocity->set(axis, iter.levelCoor(), (v - v.dot(n) * n)[axis]);
                        }
                    }
                }
            }
        });

        if(closedDomainBoundaryFlag() & kDirectionXPositive) {
            for (int i = 0; i < layout->level(); ++i) {
                const auto& res = layout->levelResolution(i);
                for (int zIdx = 0; zIdx < res.z(); ++zIdx) {
                    for (int yIdx = 0; yIdx < res.y(); ++yIdx) {
                        velocity->trySetX(i, res.x() - 1, yIdx, zIdx, 0);
                    }
                }
            }
        }
        
        if(closedDomainBoundaryFlag() & kDirectionXNegative) {
            for (int i = 0; i < layout->level(); ++i) {
                const auto& res = layout->levelResolution(i);
                for (int zIdx = 0; zIdx < res.z(); ++zIdx) {
                    for (int yIdx = 0; yIdx < res.y(); ++yIdx) {
                        velocity->trySetX(i, 0, yIdx, zIdx, 0);
                    }
                }
            }
        }

        if(closedDomainBoundaryFlag() & kDirectionYPositive) {
            for (int i = 0; i < layout->level(); ++i) {
                const auto& res = layout->levelResolution(i);
                for (int zIdx = 0; zIdx < res.z(); ++zIdx) {
                    for (int xIdx = 0; xIdx < res.x(); ++xIdx) {
                        velocity->trySetY(i, xIdx, res.y()-1, zIdx, 0);
                    }
                }
            }
        }

        if(closedDomainBoundaryFlag() & kDirectionYNegative) {
            for (int i = 0; i < layout->level(); ++i) {
                const auto& res = layout->levelResolution(i);
                for (int zIdx = 0; zIdx < res.z(); ++zIdx) {
                    for (int xIdx = 0; xIdx < res.x(); ++xIdx) {
                        velocity->trySetY(i, xIdx, 0, zIdx, 0);
                    }
                }
            }
        }

        if(closedDomainBoundaryFlag() & kDirectionZPositive) {
            for (int i = 0; i < layout->level(); ++i) {
                const auto& res = layout->levelResolution(i);
                for (int yIdx = 0; yIdx < res.y(); ++yIdx) {
                    for (int xIdx = 0; xIdx < res.x(); ++xIdx) {
                        velocity->trySetZ(i, xIdx, yIdx, res.z() - 1, 0);
                    }
                }
            }
        }

        if(closedDomainBoundaryFlag() & kDirectionZNegative) {
            for (int i = 0; i < layout->level(); ++i) {
                const auto& res = layout->levelResolution(i);
                for (int yIdx = 0; yIdx < res.y(); ++yIdx) {
                    for (int xIdx = 0; xIdx < res.x(); ++xIdx) {
                        velocity->trySetZ(i, xIdx, yIdx, 0, 0);
                    }
                }
            }
        }
    }

    bool isClosedDomain(int axis = 0) {
        switch (axis) {
            case 0:
                return closedDomainBoundaryFlag() & kDirectionXPositive;
            case 1:
                return closedDomainBoundaryFlag() & kDirectionXNegative;
            case 2:
                return closedDomainBoundaryFlag() & kDirectionYPositive;
            case 3:
                return closedDomainBoundaryFlag() & kDirectionYNegative;
            case 4:
                return closedDomainBoundaryFlag() & kDirectionZPositive;
            case 5:
                return closedDomainBoundaryFlag() & kDirectionZNegative;

            default:
                throw std::runtime_error("unknown axis");
        }
    }

    const CellCenteredScalarGrid3Ptr& colliderSdf() const {return _colliderSdf;}
    void updateSdf(CellCenteredScalarGrid3Ptr& sdf) {_colliderSdf = sdf;}

protected:    
    // Returns the size/gridSpacing/gridOrigin of the velocity grid to be constrained.
    const Size3& resolution() const {return _resolution;}
    REAL gridSpacing() const {return _gridSpacing;}

private:
    Collider3Ptr _collider = nullptr;
    CellCenteredScalarGrid3Ptr _colliderSdf = nullptr;
    Size3 _resolution;
    REAL _gridSpacing;
    int _closedDomainBoundaryFlag = kDirectionNone;
};

typedef std::shared_ptr<OmniBoundarySolver3> OmniBoundarySolver3Ptr;

}  // end of namespace Omni

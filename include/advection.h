#pragma once
#include "general.h"
#include "cell_centered_grid.h"
#include "vertex_centered_grid.h"
#include "octree_omni_grid3.h"
#include "octree_omni_face_centered_grid3.h"
#include "octree_util.h"
#include "global_benchmark.h"

using std::shared_ptr;

namespace Omni {

template<class BType, class VType, class PType>
inline PType BackTrace(const BType& boundarySdf, const VType& flow, const PType& p, REAL dt, REAL delta_h=-1) {
#ifdef USE_HIGH_BACKTRACE
    ASSERT(delta_h >= 0, "high precision backtrace need delta_h parameter");
    REAL remainingT = dt;
    PType pt0 = p;
    PType pt1 = p;

    while (remainingT > EPSILON) {
        // Adaptive time-stepping
        PType vel0 = flow->sample(pt0);
        REAL numSubSteps = std::max(std::ceil(vel0.norm() * remainingT / delta_h), 1.0);
        dt = remainingT / numSubSteps;

        // Mid-point rule
        PType midPt = pt0 - 0.5 * dt * vel0;
        PType midVel = flow->sample(midPt);
        pt1 = pt0 - dt * midVel;

        // Boundary handling
        REAL phi0 = boundarySdf->sample(pt0);
        REAL phi1 = boundarySdf->sample(pt1);

        if (phi0 * phi1 < 0.0) {
            REAL w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
            pt1 = w * pt0 + (1.0 - w) * pt1;
            break;
        }

        remainingT -= dt;
        pt0 = pt1;
    }

    return pt1;
#else
    PType midPt = p - flow->sample(p) * dt * 0.5;
    PType ret = p - flow->sample(midPt) * dt;

    if(boundarySdf) {
        REAL phi0 = boundarySdf->sample(p);
        REAL phi1 = boundarySdf->sample(ret);

        if (phi0 * phi1 < 0.0) {
            REAL w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
            return w * p + (1.0 - w) * ret;
        }
    }

    return ret;
#endif
}

template<class BType, class VType, class PType>
inline PType BackTrace(const BType& boundarySdf, const VType& flow, const PType& v0, const PType& p, REAL dt) {
    PType midPt = p - v0 * dt * 0.5;
    PType ret = p - flow->sample(midPt) * dt;

    if(boundarySdf) {
        // Boundary handling
        REAL phi0 = boundarySdf->sample(p);
        REAL phi1 = boundarySdf->sample(ret);

        if (phi0 * phi1 < 0.0) {
            REAL w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
            return w * p + (1.0 - w) * ret;
        }
    }

    return ret;
}


template<class VType, class PType>
inline PType BackTrace(const VType& flow, const PType& p, Real dt) {
    PType midPt = p - flow->sample(p) * dt * 0.5;
    PType ret = p - flow->sample(midPt) * dt;
    return ret;
}

template<class VType, class PType>
inline PType BackTrace(const VType& flow, const PType& v0, const PType& p, Real dt) {
    PType midPt = p - v0 * dt * 0.5;
    PType ret = p - flow->sample(midPt) * dt;
    return ret;
}


// NOTE: maximum backtrace 5 times
template<int maxFactor = 5, class BType, class VType, class PType>
inline PType BackTraceHighPrecision(const BType& boundarySdf, const VType& flow, const PType& p, REAL dt, REAL delta_h=-1) {
    ASSERT(delta_h >= 0, "high precision backtrace need delta_h parameter");
    PType pt0 = p;
    PType pt1 = p;

    PType vel0 = flow->sample(pt0);
    int numSubSteps = std::min(maxFactor, std::max(static_cast<int>(std::ceil(vel0.norm() * dt / delta_h)), 1));
    REAL subDt = dt / numSubSteps;
    for(int i = 0; i < numSubSteps; ++i) {
        // Mid-point rule
        PType vel0 = flow->sample(pt0);
        PType midPt = pt0 - 0.5 * subDt * vel0;
        PType midVel = flow->sample(midPt);
        pt1 = pt0 - subDt * midVel;

        if(boundarySdf) {
            // Boundary handling
            REAL phi0 = boundarySdf->sample(pt0);
            REAL phi1 = boundarySdf->sample(pt1);

            if (phi0 * phi1 < 0.0) {
                REAL w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
                pt1 = w * pt0 + (1.0 - w) * pt1;
                break;
            }
        }

        pt0 = pt1;
    }

    return pt1;
}


template<class BType>
void SemiLagrangianLevelSet3(const BType& boundary, const OctreeOmniFaceCenteredGrid3Ptr& vel,
    CellCenteredScalarGrid3Ptr& dst, const CellCenteredScalarGrid3Ptr& src, REAL dt, REAL thresholdMin=-0.3, REAL thresholdMax=0.3) {
    
    const auto& res = dst->resolution();
    ThreadPool::taskParallelFor(0, res.z(), [&](int zIdx){
        for(int yIdx = 0; yIdx < res.y(); ++yIdx) {
            for(int xIdx = 0; xIdx < res.x(); ++xIdx) {
                Vector3f pos = dst->position(xIdx, yIdx, zIdx);
                REAL cL = src->get(xIdx, yIdx, zIdx);
                if(cL >= thresholdMin && cL <= thresholdMax) {
                    Vector3f oldPos = BackTraceHighPrecision(boundary, vel, pos, dt, dst->gridSpacing());
                    dst->set(xIdx, yIdx, zIdx, src->sample(oldPos));
                } else {
                    dst->set(xIdx, yIdx, zIdx, sign(src->get(xIdx, yIdx, zIdx)*1000));// src->sample(oldPos));
                }
            }
        }
    }, 8);

    // dst->parallelForEach([&](int xIdx, int yIdx, int zIdx){
    //     Vector3f pos = dst->position(xIdx, yIdx, zIdx);
    //     if(src->get(xIdx, yIdx, zIdx) < threshold) {
    //         Vector3f oldPos = BackTraceHighPrecision(boundary, vel, pos, dt, dst->gridSpacing());
    //         dst->set(xIdx, yIdx, zIdx, src->sample(oldPos));
    //     }
    // });
}

template<typename T>
void SemiLagrangianOctreeOmni3(const OctreeOmniScalarGrid3Ptr& boundary, const OctreeOmniFaceCenteredGrid3Ptr& vel,
    shared_ptr<OctreeOmniGrid3<T>>& dst, const shared_ptr<OctreeOmniGrid3<T>>& src, REAL dt) {

    auto& tCellGrid = dst->getOctagonGrid();
    auto& tVertexGrid = dst->getTiltGrid();
    auto& layout = dst->getLayout();
    ThreadPool::parallelIterateGrid(tCellGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
        auto pos = tCellGrid->position(iter.levelCoor());
        auto oldPos = BackTrace(boundary, vel, pos, dt, layout->levelGridSpacing(iter.level));
        tCellGrid->set(iter.levelCoor(), src->sample(oldPos));
    });

    ThreadPool::parallelIterateGrid(tVertexGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
        if(!dst->getTiltEGridPtr()->get(iter.levelCoor()).is_closed) {
            auto pos = tVertexGrid->position(iter.levelCoor());
            auto oldPos = BackTrace(boundary, vel, pos, dt, layout->levelGridSpacing(iter.level));
            tVertexGrid->set(iter.levelCoor(), src->sample(oldPos));
        }
    });

    iterate_half_tilt(lCoor, dst->getLayout()) {
        auto pos = tVertexGrid->position(lCoor);
        auto oldPos = BackTrace(boundary, vel, pos, dt, layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)));
        tVertexGrid->set(lCoor, src->sample(oldPos));
    }
}


template<typename BType>
void SemiLagrangianOctreeOmniFace3(const BType& boundary, const OctreeOmniFaceCenteredGrid3Ptr& vel,
    OctreeOmniFaceCenteredGrid3Ptr& dst, const OctreeOmniFaceCenteredGrid3Ptr& src, REAL dt) {

    auto& tMACGrid = dst->getMACGrid();
    auto& tTiltGrid = dst->getTiltGrid();
    auto& layout = dst->getLayout();
    {BENCHMARK_SCOPED_TIMER_SECTION t2("mac velocity advection");
    for(int axis:{0, 1, 2}) {
        ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
            auto pos = tMACGrid->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
            auto oldPos = BackTrace(boundary, vel, pos, dt);
            tMACGrid->set(axis, iter.levelCoor(), src->sample(oldPos)[axis]);
        });
    }
    }

    {BENCHMARK_SCOPED_TIMER_SECTION t2("tily velocity advection");
    ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        if (!dst->getTiltEGridPtr()->get(iter.levelCoor()).is_closed) {
            Vector8f newValue;
            for (int k = 0; k < 8; ++k) {
                auto pos = tTiltGrid->positionTilt(iter.levelCoor(), k);
                auto oldPos = BackTrace(boundary, vel, pos, dt);
                newValue(k) = src->sample(oldPos).dot(TILT3_UNIT_DIR_UP[k]);
            }
            tTiltGrid->set(iter.levelCoor(), newValue);
        }
    });
    }

    {BENCHMARK_SCOPED_TIMER_SECTION t2("half tilt velocity advection");
    parallel_iterate_half_tilt(dst->getLayout(), [&](LCOOR_T lCoor) {
        Vector8f newValue = Vector8f::Zero();
        const auto& tNode = dst->getTiltEGridPtr()->get(lCoor);
        for (int i = 0; i < 5; ++i) {
            int k = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
            auto pos = (i == 4) ? tTiltGrid->positionCenter(lCoor) : tTiltGrid->positionTilt(lCoor, k);
            auto oldPos = BackTrace(boundary, vel, pos, dt);
            newValue(k) = src->sample(oldPos).dot((i == 4) ? HALF_TILT_XYZ_DIR_POSITIVE[k] : TILT3_UNIT_DIR_UP[k]);
        }
        tTiltGrid->set(lCoor, newValue);
    });
    }
}


template<typename BType, typename LV>
void SemiLagrangianOctreeOmniFace3WithLevelSet(const BType& boundary, const OctreeOmniFaceCenteredGrid3Ptr& vel,
    OctreeOmniFaceCenteredGrid3Ptr& dst, const OctreeOmniFaceCenteredGrid3Ptr& src, const LV& levelSet, REAL dt, REAL thresholdMax=0.3) {
    auto& tMACGrid = dst->getMACGrid();
    auto& tTiltGrid = dst->getTiltGrid();
    auto& layout = dst->getLayout();
    // int groupSize = dst->resolution().x() * dst->resolution().y() * 4;
    
    ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
    // ThreadPool::parallelIterateGrid(tMACGrid->getIterator(), groupSize, [&](const OctreeGridIterator3& iter) {
        auto pos = dst->getDualGrid()->positionOctagon(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
        REAL cL = levelSet->sample(pos);
        if(cL > thresholdMax) {
            return;
        }

        for(int axis : {0, 1, 2}) {
            pos = tMACGrid->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
            auto oldPos = BackTrace(boundary, vel, pos, dt, layout->levelGridSpacing(iter.level));
            tMACGrid->set(axis, iter.levelCoor(), src->sample(oldPos)[axis]);
        }
    });

    ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
    // ThreadPool::parallelIterateGrid(tTiltGrid->getIterator(), groupSize, [&](const OctreeGridIterator3& iter) {
        if(dst->getTiltEGridPtr()->get(iter.levelCoor()).is_closed)
            return;
        auto pos = tTiltGrid->positionCenter(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
        REAL cL = levelSet->sample(pos);
        if(cL > thresholdMax) {
            return;
        }

        Vector8f newValue;
        for (int k = 0; k < 8; ++k) {
            auto pos = tTiltGrid->positionTilt(iter.levelCoor(), k);
            auto oldPos = BackTrace(boundary, vel, pos, dt, layout->levelGridSpacing(iter.level));
            newValue(k) = src->sample(oldPos).dot(TILT3_UNIT_DIR_UP[k]);
        }
        tTiltGrid->set(iter.levelCoor(), newValue);
    });

    iterate_half_tilt(lCoor, dst->getLayout()) {
        auto pos = tTiltGrid->positionCenter(lCoor);
        REAL cL = levelSet->sample(pos);
        if(cL > thresholdMax) {
            continue;
        }
        
        Vector8f newValue = Vector8f::Zero();
        const auto& tNode = dst->getTiltEGridPtr()->get(lCoor);
        for (int i = 0; i < 5; ++i) {
            int k = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
            auto pos = (i == 4) ? tTiltGrid->positionCenter(lCoor) : tTiltGrid->positionTilt(lCoor, k);
            auto oldPos = BackTrace(boundary, vel, pos, dt, layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)));
            newValue(k) = src->sample(oldPos).dot((i == 4) ? HALF_TILT_XYZ_DIR_POSITIVE[k] : TILT3_UNIT_DIR_UP[k]);
        }
        tTiltGrid->set(lCoor, newValue);
    }
}

}
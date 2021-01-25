#include "octree_omnigrid_solver3.h"
#include "octree_util.h"
#include "advection.h"
#include "viewer/image_dump.h"
#include "viewer/grid_dump.h"
#include "geometry/sphere3.h"
#include "geometry/rectangle3.h"
#include "geometry/cylinder3.h"
#include "geometry/rigid_body_collider3.h"
#include "geometry/implicit_surface3.h"
#include "geometry/implicit_surface3set.h"
#include "geometry/shape_to_implicit3.h"
#include "global_benchmark.h"

using namespace Omni;
// #define EXTRAPOLATE_DENSITY

template<typename F> 
inline REAL phi_curvature(const F& field, const Vector3f& pos, REAL delta_x) {
	REAL center = field->sample(pos);
	REAL v = 0;
    REAL dx = delta_x;
	for(int k = 0; k < 6; ++k) {
		Vector3f nP = pos + OCT_NB_DXDYDZ[k].template cast<REAL>() * dx;
		v += field->sample(nP);
	}
	return (v - center * 6) / dx / dx;
}

// pos0: inside fluid, pos1: outside fluid
template<typename F>
inline REAL surface_tension_helper(const F& field, REAL theta, const Vector3f& pos0, const Vector3f pos1, REAL delta_x) {
    REAL c0 = phi_curvature(field, pos0, delta_x);
    REAL c1 = phi_curvature(field, pos1, delta_x);
    return c0 * (1/theta - 1) + c1;
}

// return theta
template<typename F>
inline REAL ghost_fluid_ratio(const F& field, REAL phi0, const Vector3f& pos) {
    REAL phi1 = field->sample(pos);
    return std::max(fractionInside(phi0, phi1), 1e-3);
}

OctreeOmniGridSolver3::OctreeOmniGridSolver3(const Size3& res, REAL spac, const std::string& name) 
    : Solver(name), resolution(res), sceneSetup(nullptr), delta_h(spac), 
        splitMapSetup(nullptr), forceModifyCB(nullptr), constrainVelCB(nullptr)
{}


void OctreeOmniGridSolver3::initialize() {
    spdlog::info("OctreeOmniGridSolver3::initialize");
    is_init = true;
    initData();
    initObstacle();
	initCuttingCellWeights();
}


void OctreeOmniGridSolver3::initData() {
    // solver init
#ifndef CUDA_SOLVER
    pcgSolver.set_solver_parameters(SOLVER_ACCURACY, 3000);
#endif
    splitMap = make_shared<CellCenteredFlagGrid3>(resolution, delta_h, 0);
    ASSERT(splitMapSetup, "no splitMap initialization function");
    TIMING(timer, splitMapSetup(this), "\tsplit map set up: {} s"); timer.reset();
    layout = make_shared<OctreeGridLayout3>(resolution, delta_h, numOctreeLayer, splitMap);
    spdlog::info("\tlayout init:{}, layout DOF: {}", timer.durationInSeconds(), layout->totalSize(false));
    tiltE = make_shared<OctreeTiltENodeGrid3>(resolution, delta_h, layout, TiltENode(0, false, true, false));
    
    tilt_layout = make_shared<OctreeOmniGridLayout<3>>(resolution, delta_h, layout);
    // tilt_layout_tmp = make_shared<OctreeOmniGridLayout<3>>(resolution, delta_h, layout);

    velocity = make_shared<OctreeOmniFaceCenteredGrid3>(resolution, delta_h, layout, tiltE);
    // force = make_shared<OctreeOmniFaceCenteredGrid3>(resolution, delta_h, layout, tiltE);
    weightsCC = make_shared<OctreeOmniFaceCenteredGrid3>(resolution, delta_h, layout, tiltE, false);
    weightsCC->setAll(1);
    heat = make_shared<OctreeOmniScalarGrid3>(resolution, delta_h, layout, tiltE, 0);
    density = make_shared<OctreeOmniScalarGrid3>(resolution, delta_h, layout, tiltE, 0);
    pressure = make_shared<OctreeOmniScalarGrid3>(resolution, delta_h, layout, tiltE, 0);
    pIndex = make_shared<OctreeOmniFlagGrid3>(resolution, delta_h, layout, tiltE, -1);
    
    tTiltE = make_shared<OctreeTiltENodeGrid3>(resolution, delta_h, layout, TiltENode(0, false, true, false));
    tHeat = make_shared<OctreeOmniScalarGrid3>(resolution, delta_h, layout, tTiltE, 0);
    tVelocity = make_shared<OctreeOmniFaceCenteredGrid3>(resolution, delta_h, layout, tTiltE);
    tDensity = make_shared<OctreeOmniScalarGrid3>(resolution, delta_h, layout, tTiltE, 0);

	boundarySolver = std::make_shared<OmniBoundarySolver3>(resolution, delta_h);

    timer.reset();
    constrainOctreeTiltE3(tiltE, true, true, true);
    tiltE->updateGhost();  // record half tilt once
    constrainOctreeTiltE3(tTiltE);
    tTiltE->updateGhost();
    spdlog::info("\tconstrain layout: {} s", timer.durationInSeconds());

// #ifdef DEBUG
//     ThreadPool::parallelIterateGrid(tiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
//         const auto& node = tilt_layout->getTiltNode(iter.levelCoor());
//         const auto& t_node = tiltE->get(iter.levelCoor());
//         if(!(node.is_open == !t_node.is_closed)) { throw std::runtime_error("open state error"); }
//         if(!(node.is_TJoint == t_node.is_T_joint)) { throw std::runtime_error("t joint state error"); }
//         if(!(node.is_fixed == t_node.is_fixed)) { 
//             spdlog::error("{}:{}, {}, {}: {}, {}", iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx(), node.is_fixed, t_node.is_fixed);
//             throw std::runtime_error("is fixed state error"); 
//         }
//         if(!(node.is_closed_by_TJoint == t_node.is_closed_by_T_joint)) { 
//             spdlog::error("is closed by t joint state error");
//             throw std::runtime_error("is closed by t joint state error");
//         }
//         if(!(node.is_half_tilt == t_node.is_half_tilt)) { 
//             spdlog::error("is half tilt state error");
//             throw std::runtime_error("is half tilt state error");
//         }
//         if(!(node.half_tilt_direction == t_node.half_direction)) { 
//             spdlog::error("half tilt direction state error: {}, {}", node.half_tilt_direction, t_node.half_direction);
//             throw std::runtime_error("half tilt direction state error"); 
//         }
//         Real delta_h = layout->levelGridSpacing(iter.level);
//         if(t_node.e < delta_h && node.is_full_size) { 
//             spdlog::error("tilt size state error");
//             throw std::runtime_error("tilt size state error");
//         }
//     });
// #endif
    
    ThreadPool::parallelIterateGrid(tiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        auto node = tiltE->get(iter.levelCoor());
        if(node.is_T_joint || node.is_fixed) {
            tilt_layout->getTiltNode(iter.levelCoor()).is_open = !node.is_closed;
            tilt_layout->getTiltNode(iter.levelCoor()).is_fixed = true;
            return;
        }
        auto res = layout->levelResolution(iter.level);
        if(iter.xIdx() == 0 || iter.yIdx() == 0 || iter.zIdx() == 0 
            || iter.xIdx() == res.x()-1 || iter.yIdx() == res.y()-1 || iter.zIdx() == res.z()-1) {
            node.setE(0);
            node.is_fixed = true;
            tilt_layout->getTiltNode(iter.levelCoor()).is_open = !node.is_closed;
            tilt_layout->getTiltNode(iter.levelCoor()).is_fixed = true;
            tiltE->set(iter.levelCoor(), node);
            tTiltE->set(iter.levelCoor(), node);
        }
    });

	if(!isRestart)
		dumpOctreeLayout3dBinary("./" + name + "/layout.bin", *layout);
}


void OctreeOmniGridSolver3::initObstacle() {
    spdlog::info("OmniGridSolver3::initObstacle");
#if 1
    boundarySolver->setClosedDomainBoundaryFlag( kDirectionXNegative | kDirectionXPositive 
                                            | kDirectionYNegative | kDirectionYPositive);
    Sphere3Ptr sphere = std::make_shared<Sphere3>(0.1, Transform3(Vector3f(0.5, 0.5, 0.45), Quaternion()));
	boundarySolver->updateCollider(std::make_shared<RigidBodyCollider3>(sphere));
#else   // TODO smoke net setup, need test
    auto shapeSet = make_shared<ImplicitSurface3Set>();
	Cylinder3Ptr cylinders[34];
	for (int i = 0; i < 17; i++) {
		REAL offset = 0.05 + i * 0.025;
		cylinders[i] = std::make_shared<Cylinder3>(0.0125, 2, Transform3(Vector3f(offset, offset, 0.3), Quaternion(Vector3f(0, 0, 1), PI / 4.0)));
		shapeSet->addShape(cylinders[i]);
	}
	for (int i = 0; i < 17; i++) {
		REAL offsetx = 0.05 + i * 0.025;
		REAL offsety = 0.45 - i * 0.025;
		cylinders[i + 17] = std::make_shared<Cylinder3>(0.0125, 2, Transform3(Vector3f(offsetx, offsety, 0.3), Quaternion(Vector3f(0, 0, 1), -PI / 4.0)));
		shapeSet->addShape(cylinders[i + 17]);
	}
    // auto rect0 = make_shared<Rectangle3>(Vector3f(Vector3f(0.03, 10, 10)));
    // auto rect1 = make_shared<Rectangle3>(Vector3f(Vector3f(10, 0.03, 10)));
    // auto rect2 = make_shared<Rectangle3>(Vector3f(Vector3f(0.03, 10, 10)), Transform3(Vector3f(0.5, 0, 0)));
    // auto rect3 = make_shared<Rectangle3>(Vector3f(Vector3f(10, 0.03, 10)), Transform3(Vector3f(0, 0.5, 0)));
    // shapeSet->addShape(rect0);
    // shapeSet->addShape(rect1);
    // shapeSet->addShape(rect2);
    // shapeSet->addShape(rect3);
	
	auto objsdf = make_shared<CellCenteredScalarGrid3>(splitMap->resolution(), splitMap->gridSpacing(), 10000);
	objsdf->parallelForEach([&](int xIdx, int yIdx, int zIdx) {
		Vector3f pos = objsdf->position(xIdx, yIdx, zIdx);
		//REAL cl = std::min(shapeSet->signedDistance(pos), -cylinder_implicit->signedDistance(pos));
		//REAL cl = -cylinder_implicit->signedDistance(pos);
		REAL cl = shapeSet->signedDistance(pos);
		objsdf->set(xIdx, yIdx, zIdx, cl);
	});
	boundarySolver->updateSdf(objsdf);
#endif


	boundarySdf = boundarySolver->colliderSdf();
	Viewer::grid2PNGImage("./"+name+"/boundarySdf3D.png", Size2(boundarySdf->resolution().x(), boundarySdf->resolution().z()),
		[&](int xIdx, int zIdx) {
		return boundarySdf->sample(boundarySdf->gridSpacing() * xIdx, 0.5, boundarySdf->gridSpacing() * zIdx) > 0 ? 1 : 0;
	});

    // ThreadPool::parallelIterateGrid(tiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
    //     auto node = tiltE->get(iter.levelCoor());
    //     if(node.is_T_joint || node.is_fixed) {
    //         return;
    //     }
        
    //     auto pos = tiltE->position(iter.levelCoor());
    //     if(boundarySdf->sample(pos) < 0) {
    //         node.setE(0);
    //         node.is_fixed = true;
    //         tilt_layout->getTiltNode(iter.levelCoor()).is_open = !node.is_closed;
    //         tilt_layout->getTiltNode(iter.levelCoor()).is_fixed = true;
    //         tiltE->set(iter.levelCoor(), node);
    //         tTiltE->set(iter.levelCoor(), node);
    //     }
    // });

	// if (!isRestart)
	// 	dumpOctreeOmniScalar3dBinary("./" + name + "/boundary.bin", *boundarySdf);
}


void OctreeOmniGridSolver3::onUpdate(REAL timeStep) {
    Timer t;
    if(sceneSetup)
        sceneSetup(this);
    wipeBoundary();
    TIMING(t, addForce(timeStep), "\t add force finished. Time cost: {} s");

    // extrapolateVelocity(velocity->getMACGrid());
    boundarySolver->constrainVelocity(velocity->getMACGrid());
    boundarySolver->constrainVelocity(velocity->getTiltGrid(), velocity->getMACGrid());
	velocity->updateGhost();
    
    TIMING(t, project(timeStep), "\t project finished. Time cost: {} s");
    
    // extrapolateVelocity(velocity->getMACGrid());
    boundarySolver->constrainVelocity(velocity->getMACGrid());
    boundarySolver->constrainVelocity(velocity->getTiltGrid(), velocity->getMACGrid());
    velocity->updateGhost();
    velocity->updateDualGrid();
    boundarySolver->constrainVelocity(velocity->getDualGrid());
	velocity->getDualGrid()->updateGhost();
    // static bool first = true;
    // if(first) {
    //     saveResult(0);
    //     first = false;
    // }
    TIMING(t, advect(timeStep), "\t advection finished. Time cost: {} s");
}


// void OctreeOmniGridSolver3::addForce(REAL dt) {
//     BENCHMARK_SCOPED_TIMER_SECTION t("add force");
//     auto& tMACGrid = velocity->getMACGrid();
//     auto& tTiltGrid = velocity->getTiltGrid();
//     auto& forceMACGrid = force->getMACGrid();
//     auto& forceTiltGrid = force->getTiltGrid();
  
//     for(int axis = 0; axis < 3; ++axis) {
//         ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
//             REAL newV = tMACGrid->get(axis, iter.levelCoor()) + dt * forceMACGrid->get(axis, iter.levelCoor());
//             tMACGrid->set(axis, iter.levelCoor(), newV);
//         });
//     }

//     ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
//         if(!tiltE->get(iter.levelCoor()).is_closed) {
//             Vector8f nV = tTiltGrid->get(iter.levelCoor()) + dt * forceTiltGrid->get(iter.levelCoor());
//             tTiltGrid->set(iter.levelCoor(), nV);
//         }
//     });

//     // just set all face, no need to care half-tilt direction

//     parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor) {
//         Vector8f nV = tTiltGrid->get(lCoor) + dt * forceTiltGrid->get(lCoor);
//         tTiltGrid->set(lCoor, nV);
//     });
// }


int OctreeOmniGridSolver3::updateIndex() {
    int cIdx = 0;
    static REAL threshold = 3;
    iterate_grid(iter, pIndex->getOctagonGrid()) {
        auto pos = pIndex->positionOctagon(iter.levelCoor());
        if(pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())
            && boundarySdf->sample(pos) > -layout->levelGridSpacing(iter.level) * threshold) {
            pIndex->setOctagon(iter.levelCoor(), cIdx++);
        } else {
            pIndex->setOctagon(iter.levelCoor(), -1);
        }
    }
    int octagonIdx = cIdx;

    iterate_grid(iter, pIndex->getTiltGrid()) {
        auto pos = pIndex->positionTilt(iter.levelCoor());
        if(!tiltE->get(iter.levelCoor()).is_closed 
            && boundarySdf->sample(pos) > -layout->levelGridSpacing(iter.level) * threshold) {
            pIndex->setTilt(iter.levelCoor(), cIdx++);
        } else {
            pIndex->setTilt(iter.levelCoor(), -1);
        }
    }
    int tiltIdx = cIdx;
    // Because some of half-tilt are stored in ghost cell, 
    // so here we do setAll first, then set index for half-tilt
    pIndex->setAllGhost(-1);    
    iterate_half_tilt(lCoor, layout) {
        auto pos = pIndex->positionTilt(lCoor);
        if(boundarySdf->sample(pos) > -layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)) * threshold)
            pIndex->setTilt(lCoor, cIdx++);
        else 
            pIndex->setTilt(lCoor, -1);
    }
    spdlog::info("\t\t\tOctagon index: {}, tilt index: {}, half tilt index: {}", octagonIdx, tiltIdx, cIdx);
    pIndex->updateGhost();
    return cIdx;
}


int OctreeOmniGridSolver3::updateIndexParallel() {
    BENCHMARK_SCOPED_TIMER_SECTION t("update index parallel");
    static REAL threshold = 3;
    
    int coreNum = TaskManager::tf.num_workers();
    std::vector<int> tCount(coreNum, 0);

    ThreadPool::parallelIterateGridWithCoreIdx(pIndex->getOctagonGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter, int coreID){
        auto pos = pIndex->positionOctagon(iter.levelCoor());
        if(pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())
            && boundarySdf->sample(pos) > -layout->levelGridSpacing(iter.level) * threshold) {
            pIndex->setOctagon(iter.levelCoor(), tCount[coreID]++);
        } else {
            pIndex->setOctagon(iter.levelCoor(), -1);
        }

        pos = pIndex->positionTilt(iter.levelCoor());
        // if(!tiltE->get(iter.levelCoor()).is_closed
        if(tilt_layout->isTiltOpen(iter.levelCoor())
            && boundarySdf->sample(pos) > -layout->levelGridSpacing(iter.level) * threshold) {
            pIndex->setTilt(iter.levelCoor(), tCount[coreID]++);
        } else {
            pIndex->setTilt(iter.levelCoor(), -1);
        }
    });
    for(int i = 1; i <coreNum; ++i) {
        tCount[i] += tCount[i-1];
    }

    ThreadPool::parallelIterateGridWithCoreIdx(pIndex->getOctagonGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter, int coreID){
        if(coreID == 0)
            return;
        int cIdx = pIndex->getOctagon(iter.levelCoor());
        if( cIdx >= 0)
            pIndex->setOctagon(iter.levelCoor(), cIdx+tCount[coreID-1]);

        cIdx = pIndex->getTilt(iter.levelCoor());
        if( cIdx >= 0)
            pIndex->setTilt(iter.levelCoor(), cIdx+tCount[coreID-1]);
    });

    int totalIdx = tCount[coreNum-1];
    // Because some of half-tilt are stored in ghost cell, 
    // so here we do setAll first, then set index for half-tilt
    pIndex->setAllGhost(-1);    
    iterate_half_tilt(lCoor, layout) {
        auto pos = pIndex->positionTilt(lCoor);
        if(boundarySdf->sample(pos) > -layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)) * threshold)
            pIndex->setTilt(lCoor, totalIdx++);
        else 
            pIndex->setTilt(lCoor, -1);
    }
    pIndex->updateGhost();
    spdlog::info("\t\t\tTotal count: {}, half_tilt: {}", totalIdx, totalIdx-tCount[coreNum-1]);
    return totalIdx;
}


void OctreeOmniGridSolver3::project(Real dt) {
    BENCHMARK_SCOPED_TIMER_SECTION t("project");
    // TIMING(timer, int totalSize = updateIndex(), "\t\t[update index] timer: {} s");
    int totalSize = updateIndexParallel();
    A.resize(totalSize); A.zero();
#ifdef CUDA_SOLVER
    if(b.size() != totalSize) {
#else
	if (b.rows() != totalSize) {
#endif
		b.resize(totalSize);// = VectorXd(totalSize);
		vRes.resize(totalSize);// = VectorXd(totalSize);
	}
    timer.reset();
    
    BENCHMARK_START_TIMER_SECTION("A setup");

	ThreadPool::parallelIterateGrid(pIndex->getOctagonGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        int cIdx = pIndex->getOctagon(iter.levelCoor()), tIdx = 0;
        if(cIdx < 0 || !pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx()))
            return;
        REAL delta_h = layout->levelGridSpacing(iter.level);
        REAL center = 0, tb = 0, tCenter, tCoeff;
        Size3 idx(iter.xIdx(), iter.yIdx(), iter.zIdx());
        REAL phi0 = getLevelSet() ? getLevelSet()->sample(pIndex->positionOctagon(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) : 0;
        Vector3f cPos = pIndex->positionOctagon(iter.levelCoor());
        // loop over all adjacent octagon
        for(int i = 0; i < 6; ++i) {
            int axis = i / 2;
            Size3 nbCoord = idx + OCT_NB_DXDYDZ[i];
            
            // Though it is ghost, we still need to deal with half-tilt.
            // The half-tilt will be processed when we loop over half-tilt array
            if(pIndex->getOctagonGrid()->isValid(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())
                && (layout->node(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z()).flag & LayerNodeFlagMask::IS_GHOST)) {
                continue;
            }

            // loop over all adjacent tilt cell to get area
            REAL area = delta_h * delta_h;
            for(int m = 0; m <= 1; ++m) {
                for(int n = 0; n <= 1; ++n) {
                    Size3 tIdx = idx; tIdx(axis) += ((i & 1) ? 0 : 1); 
                    tIdx((axis+1)%3) += m; tIdx((axis+2)%3) += n;
                    const auto& tNode = tiltE->get(iter.level, tIdx.x(), tIdx.y(), tIdx.z());
                    area -= tNode.e * tNode.e * 0.5;
                }
            }
            if(pIndex->getOctagonGrid()->isValid(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())) {
                REAL term = weightsCC->getMACGrid()->get(axis, iter.level, iter.xIdx()+OCT_NB_FACE_DXDYDZ[i].x(), iter.yIdx()+OCT_NB_FACE_DXDYDZ[i].y(), iter.zIdx()+OCT_NB_FACE_DXDYDZ[i].z());
                if((tIdx = pIndex->getOctagon(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())) >= 0) {
                    tCenter = term * area / delta_h; 
                    tCoeff = area * term;
                    A.set_element(cIdx, tIdx, -tCenter);
                } else if(tIdx == -2) { // free surface
                    tCoeff = area * term;
                    tCenter = area * term / delta_h;
#ifdef ENABLE_GHOST_FLUID
                    REAL sdfRatio = ghost_fluid_ratio(getLevelSet(), phi0, pIndex->positionOctagon(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z()));
                    tCenter = tCenter / sdfRatio;
#endif 
#ifdef ENABLE_SURFACE_TENSION
                    tb += area * (dt / rho) / delta_h * gamma * surface_tension_helper(getLevelSet(), sdfRatio, cPos, cPos + OCT_NB_DXDYDZ[i].cast<REAL>()*delta_h, getLevelSet()->gridSpacing());
#endif
                } else if(tIdx == -1) { // Neumann boundary
                    // NOTE: Neumann boundary is already handled by cutting-cell, so actually they are inside solid boundary and weights==0
                    ASSERT(fabs(term) < EPSILON, "cutting-cell weights is not zero");
                    tCenter = 0; tCoeff = 0;    
                    // tCoeff = area * weightsCC->getMACGrid()->get(axis, iter.level, iter.xIdx()+OCT_NB_FACE_DXDYDZ[i].x(),
                    //             iter.yIdx()+OCT_NB_FACE_DXDYDZ[i].y(), iter.zIdx()+OCT_NB_FACE_DXDYDZ[i].z()); 
                } else {
                    ASSERT(false, "octagon neighbor octagon:should not be here");
                }
            } else if(boundarySolver->isClosedDomain(i)) {
				tCenter = 0;
				tCoeff = 0; // tCoeff = area;
            } else {    // Dirichlet boundary
                // tCenter = area / delta_h; tCoeff = area;
                tCoeff = area;   // here we don't need to consider tilt, out-most tilt is always closed
                tCenter = tCoeff / delta_h;
            }
            center += tCenter;
            // spdlog::info(((i&1) ? 1 : -1) * tCoeff * velocity->getMACGrid()->get(axis, iter.level, iter.xIdx()+OCT_NB_FACE_DXDYDZ[i].x(),
            //             iter.yIdx()+OCT_NB_FACE_DXDYDZ[i].y(), iter.zIdx()+OCT_NB_FACE_DXDYDZ[i].z()));
            tb += ((i&1) ? 1 : -1) * tCoeff * velocity->getMACGrid()->get(axis, iter.level, iter.xIdx()+OCT_NB_FACE_DXDYDZ[i].x(),
                        iter.yIdx()+OCT_NB_FACE_DXDYDZ[i].y(), iter.zIdx()+OCT_NB_FACE_DXDYDZ[i].z()); 
        }

        for(int i = 0; i < 8; ++i) {
            Size3 nbCoord = idx + OCT_NB_TILT_DXDYDZ[i];
            const auto& e = tiltE->get(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z());
            int reSlot = ((i&1)?8:6) - i;
            if(!e.is_closed || e.is_half_tilt) {
                REAL area = e.e * e.e * SQRT3 * 0.5;
                REAL term = weightsCC->getTiltGrid()->get(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())(reSlot);
                if((tIdx = pIndex->getTilt(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())) >= 0) {
                    tCenter = term * area / (delta_h * SQRT3 * 0.5);
                    tCoeff = area * term;
                    A.add_to_element(cIdx, tIdx, -tCenter);
                } else if(tIdx == -2) { // free surface
                    tCoeff = area * term;
                    tCenter = area * term / (delta_h * SQRT3 * 0.5);
#ifdef ENABLE_GHOST_FLUID
                    REAL ratio = ghost_fluid_ratio(getLevelSet(), phi0, pIndex->positionTilt(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z()));
                    tCenter = tCenter / ratio;
#endif
#ifdef ENABLE_SURFACE_TENSION
                    tb += area * (dt / rho) / (delta_h * 0.5 * SQRT3) * gamma 
                        * surface_tension_helper(getLevelSet(), ratio, cPos, pIndex->positionTilt(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z()), getLevelSet()->gridSpacing());
#endif
                } else if(tIdx == -1) { // Neumann
                    // NOTE: Neumann boundary is already handled by cutting-cell, so actually they are inside solid boundary and weights==0
                    ASSERT(fabs(term) < EPSILON, "cutting-cell weights is not zero");
                    tCenter = 0; tCoeff = 0;
                    // tCoeff = area * weightsCC->getTiltGrid()->get(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())(reSlot);
                } else {
                    ASSERT(false, "should not be here");
                }

                center += tCenter;
                tb += (i<4?-1:1) * tCoeff * velocity->getTiltGrid()->get(iter.level, nbCoord.x(), nbCoord.y(), nbCoord.z())(reSlot);
            }
        }

        if(fabs(center) < EPSILON) {center = 1; tb = 0;}
        A.set_element(cIdx, cIdx, center);
        b[cIdx] = tb;
	});

	ThreadPool::parallelIterateGrid(pIndex->getTiltGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        int cIdx = pIndex->getTilt(iter.levelCoor()), tIdx;
        const auto& e = tiltE->get(iter.levelCoor());
        if(!tilt_layout->isTiltOpen(iter.levelCoor()) || cIdx < 0)
            return;
        // if(e.is_closed || cIdx < 0)
        ASSERT(!e.is_half_tilt, "half tilt should not be processed here");
        REAL center = 0, tb = 0, tCenter, tCoeff;
        const auto& u = velocity->getTiltGrid()->get(iter.levelCoor());
        const auto& tw = weightsCC->getTiltGrid()->get(iter.levelCoor());
        auto neighborLCoor = tiltNeighborOctagon3(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx(), layout);
        REAL area = e.e * e.e * SQRT3 * 0.5;
        Vector3f cPos = pIndex->positionTilt(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
        REAL phi0 = getLevelSet() ? getLevelSet()->sample(cPos) : 0;
        for(int i = 0; i < 8; ++i) {
            if(neighborLCoor[i] == -1) {
                if(!pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx()+TILT_NB_OCT_DXDYDZ[i].x(),
                    iter.yIdx()+TILT_NB_OCT_DXDYDZ[i].y(), iter.zIdx()+TILT_NB_OCT_DXDYDZ[i].z())) {
                    // if(boundarySolver->isClosedDomain(iter.level, iter.xIdx()+TILT_NB_OCT_DXDYDZ[i].x(), iter.yIdx()+TILT_NB_OCT_DXDYDZ[i].y(), iter.zIdx()+TILT_NB_OCT_DXDYDZ[i].z(), layout)){
                    if(false){  // TODO: remove me
                        tCenter = 0;
                        tCoeff = 0; //tCoeff = area;
                    } else {
                        tCenter = area / (layout->levelGridSpacing(iter.level) * SQRT3 / 3.0); // delta_h equal to center to boundary
                        tCoeff = area;
                    }
                    center += tCenter;
                    tb += (i<4?-1:1) * tCoeff * u(i);
                }
                continue;
            }

            LCOOR_T lCoor = neighborLCoor[i];
            Size3 idx(UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor));
            if(pIndex->getOctagonGrid()->isValid(UNPACK_LEVEL3(lCoor), idx.x(), idx.y(), idx.z())) {
                auto delta_h = layout->levelGridSpacing(UNPACK_LEVEL3(lCoor));
                if((tIdx = pIndex->getOctagon(lCoor)) >= 0) {
                    tCenter = tw(i) * area / (delta_h * SQRT3 * 0.5); 
                    tCoeff = area * tw(i);
                    A.set_element(cIdx, tIdx, -tCenter);
                } else if(tIdx == -2) { // free surface
                    tCoeff = area * tw(i);
                    tCenter = area * tw(i) / (delta_h * SQRT3 * 0.5);
#ifdef ENABLE_GHOST_FLUID
                    REAL ratio = ghost_fluid_ratio(getLevelSet(), phi0, pIndex->positionOctagon(lCoor));
                    tCenter = tCenter / ratio;
#endif
#ifdef ENABLE_SURFACE_TENSION
                    tb += area * (dt / rho) / (delta_h * 0.5 * SQRT3) * gamma * surface_tension_helper(getLevelSet(), ratio, cPos, pIndex->positionOctagon(lCoor), getLevelSet()->gridSpacing());
#endif
                } else if(tIdx == -1) { // Neumann
                    // NOTE: Neumann boundary is already handled by cutting-cell, so actually they are inside solid boundary and weights==0
                    ASSERT(fabs(tw(i)) < EPSILON, "cutting-cell weights is not zero");
                    tCenter = 0; tCoeff = 0;
                    // tCoeff = area * tw(i);
                } else {
                    throw std::runtime_error("unknown index type");
                }
            // } else if(boundarySolver->isClosedDomain(UNPACK_LEVEL3(lCoor), idx.x(), idx.y(), idx.z(), layout)){
            } else if(false){ // TODO remove me
                tCenter = 0;
                tCoeff = 0; //tCoeff = area;
            } else {
                tCenter = area / (layout->levelGridSpacing(iter.level) * SQRT3 / 3.0); // delta_h equal to center to boundary
                tCoeff = area;
            }
            center += tCenter;
            tb += (i<4?-1:1) * tCoeff * u(i);
        }

        if(fabs(center) < EPSILON) { center = 1; }
        A.set_element(cIdx, cIdx, center);
        b[cIdx] = tb;
	});

    // NOTE: cannot parallize this
    // parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor) {
    iterate_half_tilt(lCoor, layout) {
        // TODO FIX ME: cutting cell boundary wrong, FIXED, NEED CHECK
        // throw std::runtime_error("FIXME: half tilt cutting-cell handle is wrong");
        int level = UNPACK_LEVEL3(lCoor);
        int cIdx = pIndex->getTilt(lCoor), tIdx;    // note: can not skip cIdx < 0, since it have to handle across level
        const auto& tNode = tiltE->get(lCoor);
        REAL area = tNode.e * tNode.e * 2;
        REAL delta_h = layout->levelGridSpacing(level);
        REAL center = 0, tb = 0, tCenter, tCoeff;
        REAL phi0 = getLevelSet() ? getLevelSet()->sample(pIndex->positionTilt(lCoor)) : 0;
        const auto& u = velocity->getTiltGrid()->get(lCoor);
        const auto& tw = weightsCC->getTiltGrid()->get(lCoor);
        // get big octagon neighbor index
        Size3 idx = Size3(UNPACK_X3(lCoor)>>1, UNPACK_Y3(lCoor)>>1, UNPACK_Z3(lCoor)>>1);
        idx((tNode.half_direction)%3) -= tNode.half_direction < 3 ? 1 : 0;
        int octIdx = HALF_TILT_SLOT_INDEX[tNode.half_direction][4];
        if(pIndex->getOctagonGrid()->isValid(level+1, idx.x(), idx.y(), idx.z())) {
            if((tIdx = pIndex->getOctagon(level+1, idx.x(), idx.y(), idx.z())) >= 0) {
                // XYZ face part
                if(cIdx >= 0) {
                    A.add_to_element(tIdx, tIdx, tw(octIdx) * area / delta_h);
                    center += tw(octIdx) * area / delta_h;
                    A.add_to_element(tIdx, cIdx, - tw(octIdx) * area / delta_h);
                    A.add_to_element(cIdx, tIdx, - tw(octIdx) * area / delta_h);
                    b[tIdx] += (tNode.half_direction < 3 ? -1 : 1) * u(octIdx) * tw(octIdx) * area;
                    tb += (tNode.half_direction < 3 ? 1 : -1) * u(octIdx) * tw(octIdx) * area;
                } else if(cIdx == -2) { // free surface
                    // spdlog::info("hi cIdx:{} {} {} {} {} phi {} - {} {} {} tIdx:{}", cIdx, level, UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor), phi0, idx.x(), idx.y(), idx.z(), tIdx);
                    REAL phi1 = getLevelSet()->sample(pIndex->positionOctagon(level+1, idx.x(), idx.y(), idx.z()));
                    A.add_to_element(tIdx, tIdx, area / delta_h * std::min((phi1-phi0)/(phi1-EPSILON), 3000.0));
                    b[tIdx] += (tNode.half_direction < 3 ? -1 : 1) * u(octIdx) * area;
                } else if(cIdx == -1) { // Neumann boundary
                    ASSERT(fabs(tw(octIdx)) < EPSILON, "cutting-cell weights is not zero");
                    // spdlog::info("hi-2 cIdx:{} {} {} {} {} phi {} - {} {} {} tIdx:{}", cIdx, level, UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor), phi0, idx.x(), idx.y(), idx.z(), tIdx);
                    // b[tIdx] += (tNode.half_direction < 3 ? -1 : 1) * u(octIdx) * tw(octIdx) * area;
                }
            } else if(tIdx == -2 && cIdx >= 0) { // Dirichlet
                REAL phi1 = getLevelSet()->sample(pIndex->positionOctagon(level+1, idx.x(), idx.y(), idx.z()));
                A.add_to_element(cIdx, cIdx, area / delta_h * std::min((phi0-phi1)/(phi0-EPSILON), 3000.0));
                tb += (tNode.half_direction < 3 ? 1 : -1) * u(octIdx) * area;
            } else if(tIdx == -1  && cIdx >= 0) { // Neumann
                ASSERT(fabs(tw(octIdx)) < EPSILON, "cutting-cell weights is not zero");
                // tb += (tNode.half_direction < 3 ? 1 : -1) * u(octIdx) * tw(octIdx) * area;
            } else {
                ASSERT(cIdx < 0, "tilt neighbor octagon:should not be here");
            }
        } else if(cIdx >= 0) {
            int axis = (tNode.half_direction % 3) * 2 + (tNode.half_direction >= 3);
            if(boundarySolver->isClosedDomain(axis)) {
                tCenter = 0; tCoeff = area;
            } else {
                tCenter = 1000; tCoeff = area; // TODO check correctness
            }
            center += tCenter;
            tb += (tNode.half_direction < 3 ? 1 : -1) * u(octIdx) * area;
        }
        if(cIdx < 0)
            continue;

        // other neighbor
        area = tNode.e * tNode.e * SQRT3 * 0.5;
        auto neighborLCoor = halfTiltNeighborOctagon3(level, UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor), tNode.half_direction);
        for(int i = 0; i < 8; ++i) {
            if(neighborLCoor[i] == -1) 
                continue;
            ASSERT(neighborLCoor[i] >= 0, "should be positive");
            if(pIndex->getOctagonGrid()->isValid(UNPACK_LEVEL3(neighborLCoor[i]), 
            UNPACK_X3(neighborLCoor[i]), UNPACK_Y3(neighborLCoor[i]), UNPACK_Z3(neighborLCoor[i]))) {
                auto delta_h = layout->levelGridSpacing(UNPACK_LEVEL3(neighborLCoor[i]));
                if((tIdx = pIndex->getOctagon(neighborLCoor[i])) >= 0) {
                    tCenter = area * tw(i) / (delta_h * SQRT3 * 0.5); tCoeff = area * tw(i);
                    A.set_element(cIdx, tIdx, -tCenter);
                } else if(tIdx == -2) { // free surface
                    tCoeff = area; 
                    REAL phi1 = getLevelSet()->sample(pIndex->positionOctagon(neighborLCoor[i]));
                    tCenter = tCoeff * std::min((phi0-phi1)/(phi0-EPSILON), 3000.0) / (delta_h * SQRT3 * 0.5);
                } else if(tIdx == -1) { // Neumann
                    ASSERT(fabs(tw(octIdx)) < EPSILON, "cutting-cell weights is not zero");
                    tCenter = 0; tCoeff = 0;
                    // tCoeff = area * tw(i);
                } else {
                    ASSERT(false, "tilt neighbor:should not be here");
                }
            } else {
                    // normally it won't be here. The out most tilt are all closed
                throw std::runtime_error("it should not be here. All out-most tilt should be closed");
            }
            center += tCenter;
            tb += (i<4?-1:1) * tCoeff * u(i);
        }
        A.add_to_element(cIdx, cIdx, center);
        ASSERT(fabs(A(cIdx, cIdx)) > EPSILON, "should not be 0");
        b[cIdx] = tb;
    }//);

    BLAS::scale(rho / dt, b);
    BENCHMARK_STOP_TIMER_SECTION("A setup");
    
    ASSERT(isSymmetrical(A), "func symmetric!");
#ifdef CUDA_SOLVER
    DEBUG_ONLY(for(int i = 0; i < b.size(); ++i) { if(std::isnan(b[i])) throw std::runtime_error("b is nan");});
#else
    DEBUG_ONLY(for(int i = 0; i < b.rows(); ++i) { if(std::isnan(b[i])) throw std::runtime_error("b is nan");});
#endif

    {
        REAL err; int iter;
        timer.reset();
#ifdef CUDA_SOLVER
        BENCHMARK_SCOPED_TIMER_SECTION tpcg("cuda solve");
        fixedA.construct_from_matrix(A);
        cudaSolver.solve(totalSize, fixedA.rowstart, fixedA.colindex, fixedA.value, b, vRes);
        err = cudaSolver.solver_error; 
        iter = cudaSolver.solver_iters;
#else
        BENCHMARK_SCOPED_TIMER_SECTION tpcg("pcg solve");
        pcgSolver.solve(A, b, vRes, err, iter);
#endif
        spdlog::info("\t\t[project] MICPCG error:{}, iter:{}, timing: {} s", err, iter, timer.durationInSeconds()); timer.reset();
        DEBUG_ONLY(if (iter > 1000) {spdlog::error("[project] abnormal iteration number: {}, error:{}", iter, err);throw std::runtime_error("[project] abnormal iteration number");});
    }

    updateVelocityByPressure(vRes, dt);
}

#ifdef CUDA_SOLVER
void OctreeOmniGridSolver3::updateVelocityByPressure(const std::vector<REAL>& p, REAL dt) {
#else
void OctreeOmniGridSolver3::updateVelocityByPressure(const VectorXd& p, REAL dt) {
#endif
    BENCHMARK_SCOPED_TIMER_SECTION t("update velocity by pressure");
    auto& tMACGrid = velocity->getMACGrid();
    auto& tTiltGrid = velocity->getTiltGrid();

    // REAL maxP = -1;
    // iterate_grid(iter, pressure->getOctagonGrid()) {
    //     int idx = pIndex->getOctagon(iter.levelCoor());
    //     if(idx >= 0) {
    //         pressure->setOctagon(iter.levelCoor(), p[idx]);
    //         maxP = std::max(maxP, p[idx]);
    //     }
    // }
    // spdlog::info("current maxP: {}", maxP);
    // pressure->updateGhost();
    // static int index = 0;
    // Viewer::grid2PNGImage("./"+name+"/pressure/p_sample_"+paddingStr(std::to_string(index), '0', 4)+".png", 
    //     Size2(resolution.x(), resolution.z()) * 2, [&](int xIdx, int zIdx)->REAL{
    //         Vector3f pos = Vector3f(xIdx * delta_h / 2.0, 0.5, zIdx * delta_h / 2.0);
    //         return std::min(std::max(pressure->sample(pos) * 10, 0.0), 1.0);
    // });
    // ++index;

    // for(int i = 0; i < getLevelSet()->resolution().z(); ++i) {
    //     Viewer::grid2PNGImage("./"+name+"/levelset_debug/l_"+paddingStr(std::to_string(i), '0', 4)+".png", 
    //         Size2(resolution.x(), resolution.y()) * 2, [&](int xIdx, int yIdx)->REAL{
    //             Vector3f pos = Vector3f(xIdx * delta_h / 2.0, yIdx * delta_h / 2.0, i * delta_h / 2.0);
    //             return std::min(std::max(fabs(getLevelSet()->sample(pos)) * 10, 0.0), 1.0);
    //     });
    // }
    // spdlog::info(getLevelSet()->sample(Vector3f(0, 0, 38 * delta_h / 2.0), true));
    // spdlog::info(getLevelSet()->sample(Vector3f(0.3, 0.3, 38 * delta_h / 2.0), true));

    ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
    // ThreadPool::parallelIterateGrid(tMACGrid->getIterator(), TaskManager::tf.num_workers() * 4, [&](const OctreeGridIterator3& iter) {
        auto delta_h = layout->levelGridSpacing(iter.level);
        int cIdx = -1, tIdx = -1;
        REAL p0 = 0, p1 = 0;
        if(pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
            if((cIdx = pIndex->getOctagon(iter.levelCoor())) >= 0 )
                p0 = p[cIdx];
        }
        Vector3f cPos = pIndex->positionOctagon(iter.levelCoor());
        REAL l0 = getLevelSet() ? getLevelSet()->sample(cPos) : 0;
        
        for(int axis = 0; axis < 3; ++axis) {
            if(!tMACGrid->isValid(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx()))
                continue;
            if (weightsCC->getMACGrid()->get(axis, iter.levelCoor()) <= 0) {
                tMACGrid->set(axis, iter.levelCoor(), 0);
            } else {
                Size3 idx(iter.xIdx(), iter.yIdx(), iter.zIdx());
                idx(axis) -= 1;
                p1 = 0; // NOTE: assume Dirichlet here. neumann boundary will be enforced after 
                if (pIndex->getOctagonGrid()->isValid(iter.level, idx.x(), idx.y(), idx.z())) {
                    if ((tIdx = pIndex->getOctagon(iter.level, idx.x(), idx.y(), idx.z())) >= 0)
                        p1 = p[tIdx];
                    if(cIdx == -2 && tIdx >= 0) {
                        p0 = 0;
#ifdef ENABLE_GHOST_FLUID
                        REAL ratio = ghost_fluid_ratio(getLevelSet(), l0, pIndex->positionOctagon(iter.level, idx.x(), idx.y(), idx.z()));
                        p0 = (ratio - 1) / ratio * p1;
#endif 
#ifdef ENABLE_SURFACE_TENSION
                        p0 += gamma * surface_tension_helper(getLevelSet(), ratio, cPos+OCT_NB_DXDYDZ[axis*2+1].cast<REAL>()*delta_h, cPos, getLevelSet()->gridSpacing());
#endif
                    } else if(cIdx >= 0 && tIdx == -2) {
                        p1 = 0;
#ifdef ENABLE_GHOST_FLUID
                        REAL ratio = ghost_fluid_ratio(getLevelSet(), l0, pIndex->positionOctagon(iter.level, idx.x(), idx.y(), idx.z()));
                        p1 = (ratio - 1) / ratio * p0;
#endif 
#ifdef ENABLE_SURFACE_TENSION
                        p1 += gamma * surface_tension_helper(getLevelSet(), ratio, cPos, cPos+OCT_NB_DXDYDZ[axis*2+1].cast<REAL>()*delta_h, getLevelSet()->gridSpacing());
#endif
                    }
                }
                tMACGrid->set(axis, iter.levelCoor(), tMACGrid->get(axis, iter.levelCoor()) - dt / rho * (p0 - p1) / delta_h);
            }
        }
    });

    ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        // if(tNode.is_closed)
        if(!tilt_layout->isTiltOpen(iter.levelCoor()))
            return;
        
        const auto& tNode = tiltE->get(iter.levelCoor());
        int cIdx = pIndex->getTilt(iter.levelCoor()); //, tIdx;
        REAL p0 = cIdx >= 0 ? p[cIdx] : 0;
        Vector3f cPos = tTiltGrid->positionCenter(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
        REAL l0 = getLevelSet() ? getLevelSet()->sample(cPos) : 0;
        Vector8f nV = tTiltGrid->get(iter.levelCoor());
        const Vector8f& tw = weightsCC->getTiltGrid()->get(iter.levelCoor());
        auto neighborLCoor = tiltNeighborOctagon3(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx(), layout);
        for(int k = 0; k < 8; ++k) {
            LCOOR_T nLCoor = neighborLCoor(k);
            if(tw(k) > 0 && nLCoor != -1){
                LCoor3 dp(UNPACK_LEVEL3(nLCoor),  UNPACK_X3(nLCoor), UNPACK_Y3(nLCoor), UNPACK_Z3(nLCoor));
                auto delta_h = layout->levelGridSpacing(dp.level);
                REAL p1 = 0;
                int tIdx = -1;
                if (pIndex->getOctagonGrid()->isValid(dp.level, dp.xIdx, dp.yIdx, dp.zIdx)) {
                    if((tIdx = pIndex->getOctagon(nLCoor)) >= 0)
                        p1 = p[tIdx];
                    if(cIdx == -2 && tIdx >= 0) {
#ifdef ENABLE_GHOST_FLUID
                        REAL ratio = ghost_fluid_ratio(getLevelSet(), l0, pIndex->positionOctagon(dp.level, dp.xIdx, dp.yIdx, dp.zIdx));
                        p0 = (ratio - 1) / ratio * p1;
#endif
#ifdef ENABLE_SURFACE_TENSION
                        p0 += gamma * surface_tension_helper(getLevelSet(), ratio, pIndex->positionOctagon(nLCoor), cPos, getLevelSet()->gridSpacing());
#endif
                    } else if(cIdx >=0 && tIdx == -2) {
#ifdef ENABLE_GHOST_FLUID
                        REAL ratio = ghost_fluid_ratio(getLevelSet(), l0, pIndex->positionOctagon(dp.level, dp.xIdx, dp.yIdx, dp.zIdx));
                        p1 = (ratio - 1) / ratio * p0;
#endif
#ifdef ENABLE_SURFACE_TENSION
                        p1 += gamma * surface_tension_helper(getLevelSet(), ratio, cPos, pIndex->positionOctagon(nLCoor), getLevelSet()->gridSpacing());
#endif
                    }
                }   
                nV(k) = nV(k) - dt / rho * (k < 4 ? 1 : -1) / (delta_h * SQRT3 * 0.5) * (p1 - p0);
            } else {
                nV(k) = 0;
            }
        }
        tTiltGrid->set(iter.levelCoor(), nV);
    });

    parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor) {
        const auto& tNode = tiltE->get(lCoor);
        int cIdx = pIndex->getTilt(lCoor), tIdx;
        int level = UNPACK_LEVEL3(lCoor);

        REAL p0 = cIdx >= 0 ? p[cIdx] : 0;
        REAL l0 = getLevelSet() ? getLevelSet()->sample(pIndex->positionTilt(lCoor)) : 0;
        Vector8f nV = tTiltGrid->get(lCoor);
        const Vector8f& tw = weightsCC->getTiltGrid()->get(lCoor);
        Size3 idx = Size3(UNPACK_X3(lCoor)>>1, UNPACK_Y3(lCoor)>>1, UNPACK_Z3(lCoor)>>1);
        idx((tNode.half_direction)%3) -= tNode.half_direction < 3 ? 1 : 0;
        int faceSlot = HALF_TILT_SLOT_INDEX[tNode.half_direction][4];
        if(tw(faceSlot) > 0) {
            int tIdx = pIndex->getOctagon(level+1, idx.x(), idx.y(), idx.z());
            REAL p1 = tIdx >= 0 ? p[tIdx] : 0;
            if(cIdx == -2 && tIdx >= 0) {
                REAL l_air = l0;
                REAL l_fluid = getLevelSet()->sample(pIndex->positionOctagon(level+1, idx.x(), idx.y(), idx.z()));
                p0 = fmax(l_air / (l_fluid-EPSILON), -1000.0) * p1;
            } else if(cIdx >= 0 && tIdx == -2) {
                REAL l_air = getLevelSet()->sample(pIndex->positionOctagon(level+1, idx.x(), idx.y(), idx.z()));
                REAL l_fluid = l0;
                p1 = fmax(l_air / (l_fluid-EPSILON), -1000.0) * p0;
            }
            nV(faceSlot) -= dt / rho * (tNode.half_direction < 3 ? 1 : -1) / layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)) * (p0 - p1);
            // if(UNPACK_LEVEL3(lCoor) == 0 && UNPACK_X3(lCoor) == 65 && UNPACK_Y3(lCoor) == 65 && UNPACK_Z3(lCoor) == 106) {
            //     spdlog::error("update tilt idx:{} {} p:{} {} l0:{}", cIdx, tIdx, p0, p1, l0);
            //     spdlog::info("dir:{}, v:{}, old:{}", tNode.half_direction, nV(faceSlot), tTiltGrid->get(lCoor)(faceSlot));
            // }
        } else {
            nV(faceSlot) = 0;
        }

        auto neighborLCoor = halfTiltNeighborOctagon3(level, UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor), tNode.half_direction);
        for(int k = 0; k < 8; ++k) {
            LCOOR_T nLCoor = neighborLCoor(k);
            if(tw(k) > 0 && nLCoor != -1){
                auto delta_h = layout->levelGridSpacing(UNPACK_LEVEL3(nLCoor));
                REAL p1 = 0;
                if (pIndex->getOctagonGrid()->isValid(UNPACK_LEVEL3(nLCoor),  UNPACK_X3(nLCoor), UNPACK_Y3(nLCoor), UNPACK_Z3(nLCoor))) {
                    if((tIdx = pIndex->getOctagon(nLCoor)) >= 0)
                        p1 = p[tIdx];
                    if(cIdx == -2 && tIdx >= 0) {
                        REAL l_air = l0;
                        REAL l_fluid = getLevelSet()->sample(pIndex->positionOctagon(nLCoor));
                        p0 = fmax(l_air / (l_fluid-EPSILON), -1000.0) * p1;
                    } else if(cIdx >= 0 && tIdx == -2) {
                        REAL l_air = getLevelSet()->sample(pIndex->positionOctagon(nLCoor));
                        REAL l_fluid = l0;
                        p1 = fmax(l_air / (l_fluid-EPSILON), -1000.0) * p0;
                    }
                }
                nV(k) = nV(k) - dt / rho * (k < 4 ? 1 : -1) / (delta_h * SQRT3 * 0.5) * (p1 - p0);
            } else if(k != faceSlot){
                nV(k) = 0;
            }
        }
        tTiltGrid->set(lCoor, nV);
    });
}


void OctreeOmniGridSolver3::advect(REAL dt) {
    BENCHMARK_SCOPED_TIMER_SECTION t("advect");
    if(enableTiltAdaptivity)
        adjustE(0.02, dt, EPSILON, 1000);  // for large-scele demo, it is threshold is 0.2

    // SemiLagrangianOctreeOmni3(boundarySdf, velocity, tDensity, density, dt);
    // density.swap(tDensity);
    // density->updateGhost();

    // SemiLagrangianOctreeOmni3(boundarySdf, velocity, tHeat, heat, dt);
    // heat.swap(tHeat);
    // heat->updateGhost();
    auto& dOctagonGrid = tDensity->getOctagonGrid();
    auto& dTiltGrid = tDensity->getTiltGrid();
    auto& hOctagonGrid = tHeat->getOctagonGrid();
    auto& hTiltGrid = tHeat->getTiltGrid();

    {
        BENCHMARK_SCOPED_TIMER_SECTION t1("heat&density");
        ThreadPool::parallelIterateGrid(dOctagonGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
            auto pos = dOctagonGrid->position(iter.levelCoor());
            auto oldPos = BackTrace(boundarySdf, velocity, pos, dt, layout->levelGridSpacing(iter.level));
            
            dOctagonGrid->set(iter.levelCoor(), density->sample(oldPos));
            hOctagonGrid->set(iter.levelCoor(), heat->sample(oldPos));
        });
    }

    {
        BENCHMARK_SCOPED_TIMER_SECTION t1("heat&density tilt");
        ThreadPool::parallelIterateGrid(dTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
            // if(!tTiltE->get(iter.levelCoor()).is_closed) {
            if(tilt_layout->isTiltOpen(iter.levelCoor())) {
                auto pos = dTiltGrid->position(iter.levelCoor());
                auto oldPos = BackTrace(boundarySdf, velocity, pos, dt, layout->levelGridSpacing(iter.level));
                dTiltGrid->set(iter.levelCoor(), density->sample(oldPos));
                hTiltGrid->set(iter.levelCoor(), heat->sample(oldPos));
            }
        });
    }

    parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor) {
        auto pos = dTiltGrid->position(lCoor);
        auto oldPos = BackTrace(boundarySdf, velocity, pos, dt, layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)));
        dTiltGrid->set(lCoor, density->sample(oldPos));
        hTiltGrid->set(lCoor, heat->sample(oldPos));
    });

    density.swap(tDensity);
    heat.swap(tHeat);
    density->updateGhost();
    heat->updateGhost();

    {
        BENCHMARK_SCOPED_TIMER_SECTION t1("velocity");
        SemiLagrangianOctreeOmniFace3(boundarySdf, velocity, tVelocity, velocity, dt);
    }
    velocity.swap(tVelocity);

    tiltE.swap(tTiltE);

    // boundarySolver->constrainVelocity(velocity);
    density->scale(_smokeDecayFactor);
    heat->scale(_temperatureDecayFactor);
}


void OctreeOmniGridSolver3::addForce(Real dt) {
    BENCHMARK_SCOPED_TIMER_SECTION t("add force");
    auto& tMACGrid = velocity->getMACGrid();
    auto& tTiltGrid = velocity->getTiltGrid();

    Real ambT = 0.0;
    int count = 0;
    for(int i = 0; i < heat->getOctagonGrid()->getData().size(); ++i) {
        Real tSum = BLAS::sum(heat->getOctagonGrid()->getData()[i]);
        ambT += tSum * ( i == 0 ? 1 : 8);
        count += heat->getOctagonGrid()->getData()[i].size() * ( i == 0 ? 1 : 8);
    }
    ambT /= Real(count);
    spdlog::info("ambT: {}, {}, {}, {}", ambT, count, _buoyancyTemperatureFactor, _buoyancySmokeDensityFactor);

	ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        if(tMACGrid->isValidZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
			Vector3f pos = tMACGrid->positionZ(iter.levelCoor());
			REAL newF = _buoyancySmokeDensityFactor * density->sample(pos) + _buoyancyTemperatureFactor * (heat->sample(pos) - ambT);
            ASSERT(!std::isnan(newF), "heat sample is nan");
            Real newV = tMACGrid->getZ(iter.levelCoor()) + dt * newF;
            tMACGrid->setZ(iter.levelCoor(), newV);
        }
	});

	ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        if(tilt_layout->isTiltOpen(iter.levelCoor())) {
        // if(!tiltE->get(iter.levelCoor()).is_closed) {
            Vector8f newF = Vector8f::Zero();
            for(int k = 0; k < 8; ++k) {
				Vector3f pos = tTiltGrid->positionTilt(iter.levelCoor(), k);
				REAL tf = _buoyancySmokeDensityFactor * density->sample(pos) + _buoyancyTemperatureFactor * (heat->sample(pos) - ambT);

                ASSERT(!std::isnan(tf), "heat sample is nan");
                newF(k) = TILT3_UNIT_DIR_UP[k].z() * tf; // dir.dot(vector3f(0, 0, tf));
            }
            Vector8f newV = tTiltGrid->get(iter.levelCoor()) + dt * newF;
            tTiltGrid->set(iter.levelCoor(), newV);
        }
	});

    parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor) {
		Vector8f newF = Vector8f::Zero();
        const auto& tNode = tiltE->get(lCoor);
        for(int i = 0; i < 5; ++i) {
            int k = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
            auto pos = (i == 4) ? tTiltGrid->positionCenter(lCoor) : tTiltGrid->positionTilt(lCoor, k);
			REAL tf = _buoyancySmokeDensityFactor * density->sample(pos) + _buoyancyTemperatureFactor * (heat->sample(pos) - ambT);
            ASSERT(!std::isnan(tf), "heat sample is nan");

            newF(k) = Vector3f(0, 0, tf).dot((i == 4) ? HALF_TILT_XYZ_DIR_POSITIVE[k] : TILT3_UNIT_DIR_UP[k]);
        }
        Vector8f newV = tTiltGrid->get(lCoor) + dt * newF;
        tTiltGrid->set(lCoor, newV);
    });
}

// #define WRITE_RESULT

void OctreeOmniGridSolver3::saveResult(int index) {
    BENCHMARK_SCOPED_TIMER_SECTION t("save result");
    int yIdx = resolution.y() / 2;
    int scale = 2;
	auto copy_density = new OctreeOmniScalarGrid3(*density);
	auto copy_tilt = new OctreeTiltENodeGrid3(*tiltE);
#ifdef WRITE_RESULT
	auto copy_velocity = new OctreeOmniFaceCenteredGrid3(*velocity);
	auto copy_heat = new OctreeOmniScalarGrid3(*heat);
#endif

    aCounter.increase();
    copy_density->updateGhost();    // It is parallel now. Do it now to avoid nested parallel problem
	ResultManager::dispatch([=]() {
		Viewer::grid2PNGImage("./" + name + "/preview_density/d_sample_" + paddingStr(std::to_string(index), '0', 4) + ".png",
			Size2(resolution.x(), resolution.z()) * scale, [&](int xIdx, int zIdx)->REAL {
			Vector3f pos = Vector3f(xIdx, yIdx * scale, zIdx) * delta_h / (REAL)scale;
			return std::min(std::max(copy_density->sample(pos), 0.0), 1.0);
		});
        if(isDumpDensity) {
            dumpOctreeOmniGrid3DBinary("./" + name + "/sampled_density/d_" + paddingStr(std::to_string(index), '0', 4) + ".bin",
                resolution, *copy_density);
        }
        
#ifdef WRITE_RESULT
		dumpOctreeOmniScalar3dBinary("./" + name + "/density/d_" + paddingStr(std::to_string(index), '0', 4) + ".bin", *copy_density);
		dumpOctreeOmniScalar3dBinary("./" + name + "/heat/h_" + paddingStr(std::to_string(index), '0', 4) + ".bin", *copy_heat);
		dumpOctreeOmniFaceCentered3dBinary("./" + name + "/velocity/v_" + paddingStr(std::to_string(index), '0', 4) + ".bin", *copy_velocity);
		delete copy_heat;
		delete copy_velocity;
#endif
		dumpOctreeOmniTiltE3dBinary("./" + name + "/tilt/t_" + paddingStr(std::to_string(index), '0', 4) + ".bin", *copy_tilt);

		delete copy_tilt;
		delete copy_density;
        aCounter.decrease();
	});

    // pressure->updateGhost();
    // Viewer::grid2PNGImage("./"+name+"/pressure_sample_"+paddingStr(std::to_string(index), '0', 4)+".png", 
    //     Size2(resolution.x(), resolution.z()) * scale, [&](int xIdx, int zIdx)->REAL{
    //         Vector3f pos = Vector3f(xIdx, yIdx * scale, zIdx) * delta_h / (REAL)scale;
    //         return std::min(std::max(pressure->sample(pos) * 10, 0.0), 1.0);
    // });
    // Viewer::grid2PNGImage("./"+name+"/density_"+paddingStr(std::to_string(index), '0', 4)+".png", 
    //     Size2(resolution.x(), resolution.z()*2), [&](int xIdx, int zIdx)->REAL{
    //         int c0 = layout->getUntil(0, xIdx, yIdx, zIdx>>1, false);
    //         REAL v = 0;
    //         if(c0 >= 0) {
    //             if(zIdx&1)
    //                 v = density->getOctagon(c0);
    //             else
    //                 v = density->getTilt(c0);
    //         } 
    //         return std::min(std::max(v, 0.0), 1.0);
    // });
    // velocity->updateGhost();
    // velocity->updateDualGrid();
    // velocity->getDualGrid()->updateGhost();
    // Viewer::grid2PNGImage("./"+name+"/velocity_"+paddingStr(std::to_string(index), '0', 4)+".png", 
    //     Size2(resolution.x(), resolution.z()*2), [&](int xIdx, int zIdx)->REAL{
    //         int c0 = layout->getUntil(0, xIdx, yIdx, zIdx>>1, false);
    //         REAL v = 0;
    //         if(c0 >= 0) {
    //             if(zIdx&1)
    //                 v = velocity->getDualGrid()->getOctagon(c0).norm();
    //             else
    //                 v = velocity->getDualGrid()->getTilt(c0).norm();
    //         } 
    //         return std::min(std::max(v*4, 0.0), 1.0);
    // });
    // Viewer::grid2PNGImage("./"+name+"/velocity_sample_"+paddingStr(std::to_string(index), '0', 4)+".png", 
    //    Size2(resolution.x(), resolution.z()) * scale, [&](int xIdx, int zIdx)->REAL{
    //        Vector3f pos = Vector3f(xIdx, yIdx * scale, zIdx) * delta_h / (REAL)scale;
    //        REAL v = velocity->sample(pos).norm();
    //        return std::min(std::max(v * 4, 0.0), 1.0);
    // });
}

bool Omni::OctreeOmniGridSolver3::restart(int frame, const string & layoutfile, const string & boundaryfile, const string & tiltfile, const string & densityfile, const string & heatfile, const string & velocityfile)
{
	isRestart = true;
	spdlog::info("restart");
	cFrame = frame;
	if (!readOctreeLayout3dBinary(layoutfile, *layout))
		return false;
    // TODO fix me
	// if (!readOctreeOmniScalar3dBinary(boundaryfile, *boundarySdf))
	// 	return false;
	if (!readOctreeOmniTiltE3dBinary(tiltfile, *tiltE))
		return false;
	if (!readOctreeOmniTiltE3dBinary(tiltfile, *tTiltE))
		return false;
	if (!readOctreeOmniScalar3dBinary(densityfile, *density))
		return false;
	if (!readOctreeOmniScalar3dBinary(densityfile, *tDensity))
		return false;
	if (!readOctreeOmniScalar3dBinary(heatfile, *heat))
		return false;
	if (!readOctreeOmniScalar3dBinary(heatfile, *tHeat))
		return false;
	if (!readOctreeOmniFaceCentered3dBinary(velocityfile, *velocity))
		return false;
	if (!readOctreeOmniFaceCentered3dBinary(velocityfile, *tVelocity))
		return false;
    velocity->updateGhost();
	velocity->updateDualGrid();
	velocity->getDualGrid()->updateGhost();
    tVelocity->updateGhost();
	tVelocity->updateDualGrid();
	tVelocity->getDualGrid()->updateGhost();
	initCuttingCellWeights();

	int yIdx = resolution.y() / 2;
	int scale = 4;
	Viewer::grid2PNGImage("./" + name + "/restart.png",
		Size2(resolution.x(), resolution.z()) * scale, [&](int xIdx, int zIdx)->REAL {
		Vector3f pos = Vector3f(xIdx, yIdx * scale, zIdx) * delta_h / (REAL)scale;
		return std::min(std::max(density->sample(pos), 0.0), 1.0);
	});
	return true;
}


int OctreeOmniGridSolver3::numberOfSubTimeSteps(REAL dt) const {
    int step = std::ceil(cfl(dt) / maxCfl);
    return step > 1 ? step : 1;
}


int OctreeOmniGridSolver3::cfl(REAL dt) const {
    BENCHMARK_SCOPED_TIMER_SECTION t("cfl computation");

    auto& MACGrid = velocity->getMACGrid();
    Real maxV = 0;

    for(int axis = 0; axis < 3; ++axis) {
        for(auto& data: MACGrid->getU(axis)) {
            if(data.size() > 0)
                maxV = std::max(maxV, ThreadPool::reduceTransformMax(data, [](Real v) -> Real { return v; }));
        }
    }
    // for(auto& data: octagonGrid->getData()) {
    //     if(data.size() > 0)
    //         maxV = std::max(maxV, ThreadPool::reduceTransformMax(data, [](const Vector3f& v) -> Real {
    //             return v.norm();
    //         }));
    // }
    // auto& octagonGrid = velocity->getDualGrid()->getOctagonGrid();
    // Real maxV = 0;
    // for(auto& data: octagonGrid->getData()) {
    //     if(data.size() > 0)
    //         maxV = std::max(maxV, ThreadPool::reduceTransformMax(data, [](const Vector3f& v) -> Real {
    //             return v.norm();
    //         }));
    // }

#ifdef DEBUG
    REAL tMaxV = 0, minV = 1000;
    iterate_grid(iter, octagonGrid) {
        if(octagonGrid->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
            const auto& v = octagonGrid->get(iter.levelCoor());
            tMaxV = tMaxV > v.norm() ? tMaxV : v.norm();
            minV = minV > v.norm() ? v.norm() : minV;
        }
    }
    if(fabs(tMaxV - maxV) > 1e-7) { throw std::runtime_error("inconsistent maxV"); }
#endif

    // iterate_grid(iter, tiltGrid) {
    //     if(!tiltE->get(iter.levelCoor()).is_closed) {
    //         const auto& v = tiltGrid->get(iter.levelCoor());
    //         maxV = maxV > v.norm() ? maxV : v.norm();
    //         minV = minV > v.norm() ? v.norm() : minV;
    //     }
    // }

    // iterate_half_tilt(lCoor, layout) {
    //     const auto& v = tiltGrid->get(lCoor);
    //     maxV = maxV > v.norm() ? maxV : v.norm();
    //     minV = minV > v.norm() ? v.norm() : minV;
    // }
    
    // DEBUG_ONLY(spdlog::info("max V:{}, minV:{}", maxV, minV));
    DEBUG_ONLY(spdlog::info("max V:{}", maxV));
    return dt * maxV / delta_h;
}


void OctreeOmniGridSolver3::adjustE(REAL dThreshold, REAL dt, REAL vThreshold, REAL maxThreshold) {
    BENCHMARK_SCOPED_TIMER_SECTION t("adjust E");
    {
    BENCHMARK_SCOPED_TIMER_SECTION t_layout("adjust E layout");
	ThreadPool::parallelIterateGrid(tTiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
		if(tilt_layout->getTiltNode(iter.levelCoor()).is_fixed) return;
        if(tilt_layout->getTiltNode(iter.levelCoor()).is_open) {
            Real d = density->getTilt(iter.levelCoor());
            if (d < dThreshold || d > maxThreshold) {
                tilt_layout->getTiltNode(iter.levelCoor()).is_open = false;
            }
        } else {
            Vector3f pos = tTiltE->position(iter.levelCoor());
            Vector3f v = velocity->sample(pos);
            // Vector3f v = velocity->getDualGrid()->getOctagon(iter.levelCoor());
            if (v.norm() > vThreshold) {
                Vector3f oldPos = pos - v * dt;
                Real d = density->sample(oldPos);
            // Real d = density->getOctagon(iter.levelCoor());
                if (d > dThreshold && d < maxThreshold) {
                    tilt_layout->getTiltNode(iter.levelCoor()).is_open = true;
                }
            }
        }
	});}
    tilt_layout->updateGhost();

    {
    BENCHMARK_SCOPED_TIMER_SECTION t_copy("adjust E copy states");
	ThreadPool::parallelIterateGrid(tTiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
        if(tilt_layout->getTiltNode(iter.levelCoor()).is_fixed)
            return;
        auto& node = tTiltE->get(iter.levelCoor());
        if(tilt_layout->isTiltOpen(iter.levelCoor())) {
            Real maxE = layout->levelGridSpacing(iter.level) * 0.5;
            node.setE(maxE, maxE);
        } else {
            node.setE(0);
        }
	});}
	tTiltE->updateGhost();
}


void OctreeOmniGridSolver3::initCuttingCellWeights() {
	auto& tDualGrid = velocity->getDualGrid();
	auto& tMACGrid = weightsCC->getMACGrid();
	auto& tTiltGrid = weightsCC->getTiltGrid();
	ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
		REAL current = boundarySdf->sample(tDualGrid->positionOctagon(iter.levelCoor()));
		if (tMACGrid->isValidX(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
			if (pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx() - 1, iter.yIdx(), iter.zIdx())) {
				REAL left = boundarySdf->sample(tDualGrid->positionOctagon(iter.level, iter.xIdx() - 1, 
					iter.yIdx(), iter.zIdx()));
				tMACGrid->setX(iter.levelCoor(), 1 - fractionInside(left, current));
			}
		}
		if (tMACGrid->isValidY(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
			if (pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx(), iter.yIdx() - 1, iter.zIdx())) {
				REAL back = boundarySdf->sample(tDualGrid->positionOctagon(iter.level, iter.xIdx(), iter.yIdx() - 1, iter.zIdx()));
				tMACGrid->setY(iter.levelCoor(), 1 - fractionInside(back, current));
			}
		}
		if (tMACGrid->isValidZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
			if (pIndex->getOctagonGrid()->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx() - 1)) {
				REAL bottom = boundarySdf->sample(tDualGrid->positionOctagon(iter.level, iter.xIdx(), iter.yIdx(), 
					iter.zIdx() - 1));
				tMACGrid->setZ(iter.levelCoor(), 1 - fractionInside(bottom, current));
			}
		}
	});

	ThreadPool::parallelIterateGrid(tTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter) {
		Vector8f w = tTiltGrid->get(iter.levelCoor());
		Vector3f pos = tTiltGrid->positionCenter(iter.levelCoor());
		REAL current = boundarySdf->sample(pos);
		for (int i = 0; i < 8; i++) {
			Size3 nbIdx = Size3(iter.xIdx(), iter.yIdx(), iter.zIdx()) + TILT_NB_OCT_DXDYDZ[i];
			if (pIndex->getOctagonGrid()->isValid(iter.level, nbIdx.x(), nbIdx.y(), nbIdx.z())) {
				pos = tDualGrid->positionOctagon(iter.level, nbIdx.x(), nbIdx.y(), nbIdx.z());
				REAL nb = boundarySdf->sample(pos);
				w(i) = 1 - fractionInside(current, nb);
			}
		}
		tTiltGrid->set(iter.levelCoor(), w);
	});

	parallel_iterate_half_tilt(weightsCC->getLayout(), [&](LCOOR_T lCoor) {
		Vector8f w = Vector8f::Zero();
		Vector3f pos = tTiltGrid->positionCenter(lCoor);
		REAL current = boundarySdf->sample(pos);
		const auto& tNode = tiltE->get(lCoor);
		int level = UNPACK_LEVEL3(lCoor);

		// big octagon
		Size3 idx = Size3(UNPACK_X3(lCoor) >> 1, UNPACK_Y3(lCoor) >> 1, UNPACK_Z3(lCoor) >> 1);
		idx((tNode.half_direction) % 3) -= tNode.half_direction < 3 ? 1 : 0;
		if (pIndex->getOctagonGrid()->isValid(level + 1, idx.x(), idx.y(), idx.z())) {
			pos = tDualGrid->positionOctagon(level, idx.x(), idx.y(), idx.z());
			REAL nb = boundarySdf->sample(pos);
			w(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]) = 1 - fractionInside(current, nb);
		}
		// same level octagon
		for (int i = 0; i < 4; ++i) {
			int k = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
			Size3 nbIdx = Size3(UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor)) + TILT_NB_OCT_DXDYDZ[k];
			if (pIndex->getOctagonGrid()->isValid(level, nbIdx.x(), nbIdx.y(), nbIdx.z())) {
				REAL nb = boundarySdf->sample(pos);
				w(k) = 1 - fractionInside(current, nb);
			}
		}
		tTiltGrid->set(lCoor, w);
	});
}

void Omni::OctreeOmniGridSolver3::wipeBoundary() {
    BENCHMARK_SCOPED_TIMER_SECTION t("wipe boundary");
    auto& dOctagonGrid = density->getOctagonGrid();
    auto& dTiltGrid = density->getTiltGrid();
	ThreadPool::parallelIterateGrid(dOctagonGrid->getParallelIteratorVec(), [&](const auto& iter) {
		auto pos = dOctagonGrid->position(iter.levelCoor());
		if (boundarySdf->sample(pos) < 0) {
            Real t_d = 0, t_h = 0;
            Real count = EPSILON;
#ifdef EXTRAPOLATE_DENSITY
            for(int k = 0; k < 6; ++k) {
                int x = OCT_NB_DXDYDZ[k].x();
                int y = OCT_NB_DXDYDZ[k].y();
                int z = OCT_NB_DXDYDZ[k].z();
                Vector3f nPos = pos + OCT_NB_DXDYDZ[k].template cast<Real>() * delta_h;
                if(boundarySdf->sample(nPos) >= 0) {
                    t_d += dOctagonGrid->get(iter.level, iter.xIdx()+x, iter.yIdx()+y, iter.zIdx()+z);
                    t_h += heat->getOctagonGrid()->get(iter.level, iter.xIdx()+x, iter.yIdx()+y, iter.zIdx()+z);
                    count += 1;
                }
            }
#endif
			dOctagonGrid->set(iter.levelCoor(), t_d / count);
			heat->getOctagonGrid()->set(iter.levelCoor(), t_h / count);
		}
	});
	ThreadPool::parallelIterateGrid(dTiltGrid->getParallelIteratorVec(), [&](const auto& iter) {
		auto pos = dTiltGrid->position(iter.levelCoor());
		if (boundarySdf->sample(pos) < 0) {
            Real t_d = 0, t_h = 0;
            Real count = EPSILON;
#ifdef EXTRAPOLATE_DENSITY
            for(int k = 0; k < 8; ++k) {
                Vector3f nPos = pos + TILT_NB_OCT_DXDYDZ[k].template cast<Real>() * delta_h;
                if(boundarySdf->sample(nPos) >= 0) {
                    t_d += dOctagonGrid->get(iter.level, iter.xIdx()+TILT_NB_OCT_DXDYDZ[k].x(), iter.yIdx()+TILT_NB_OCT_DXDYDZ[k].y(), iter.zIdx()+TILT_NB_OCT_DXDYDZ[k].z());
                    t_h += heat->getOctagonGrid()->get(iter.level, iter.xIdx()+TILT_NB_OCT_DXDYDZ[k].x(), iter.yIdx()+TILT_NB_OCT_DXDYDZ[k].y(), iter.zIdx()+TILT_NB_OCT_DXDYDZ[k].z());
                    count += 1;
                }
            }
#endif
			dTiltGrid->set(iter.levelCoor(), t_d / count);
			heat->getTiltGrid()->set(iter.levelCoor(), t_h / count);
		}
	});

	parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor) {
		Vector3f pos = density->positionTilt(lCoor);
		if (boundarySdf->sample(pos) < 0) {
			density->setTilt(lCoor, 0);
			heat->setTilt(lCoor, 0);
		}
	});
}

void Omni::OctreeOmniGridSolver3::extrapolateVelocity(OctreeFaceCenteredGrid3Ptr& velocity) {
    BENCHMARK_SCOPED_TIMER_SECTION t("extrapolate velocity");

    for(int axis = 0; axis < 3; ++axis) {

        if(!marker0) {
            marker0 = std::make_shared<OctreeCellCenteredCharGrid3>(resolution, delta_h, layout, 0);
            marker1 = std::make_shared<OctreeCellCenteredCharGrid3>(resolution, delta_h, layout, 0);
        } else {
            marker0->setAll(0);
        }

        ThreadPool::parallelIterateGrid(velocity->getParallelIteratorVec(), [&](const auto& iter) {
            if(velocity->isValid(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
                auto pos = velocity->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
                if(weightsCC->getMACGrid()->get(axis, iter.levelCoor()) > EPSILON) {
                    marker0->set(iter.levelCoor(), 1);
                }
            }
        });
        marker1->fill(marker0);

        for(int i = 0; i < 4; ++i) {
            ThreadPool::parallelIterateGrid(velocity->getParallelIteratorVec(), [&](const auto& iter) {
                if(velocity->isValid(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())
                    && !marker0->get(iter.levelCoor())) {
                    Real t_v = 0;
                    Real count = 0;
                    for(int k = 0; k < 6; ++k) {
                        int x = OCT_NB_DXDYDZ[k].x();
                        int y = OCT_NB_DXDYDZ[k].y();
                        int z = OCT_NB_DXDYDZ[k].z();
                        if(velocity->isValid(axis, iter.level, iter.xIdx()+x, iter.yIdx()+y, iter.zIdx()+z)
                            && marker0->get(iter.level, iter.xIdx()+x, iter.yIdx()+y, iter.zIdx()+z)){
                            t_v += velocity->get(axis, iter.level, iter.xIdx()+x, iter.yIdx()+y, iter.zIdx()+z);
                            count += 1;
                        }
                    }
                    if(count > 0) {
                        velocity->set(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx(), t_v / count);
                        marker1->set(iter.levelCoor(), 1);
                    }
                } else {
                    marker1->set(iter.levelCoor(), 1);
                }
            });
            marker0.swap(marker1);
        }
    }
}



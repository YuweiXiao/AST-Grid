#include "general.h"
#include "octree_omnigrid_solver3.h"
#include "viewer/image_dump.h"
#include "geometry/sphere3.h"
#include "geometry/cylinder3.h"
#include "geometry/rectangle3.h"
#include "geometry/shape_to_implicit3.h"
#include "octree_util.h"
#include "args_parsing.h"

#define RESX 128
#define RESY 128
#define RESZ 128

#define TIME_STEP 1.0 / 60

using namespace Omni;
using namespace Omni::Geometry;

Size3 res(RESX, RESY, RESZ);
REAL delta_h = 1.0 / RESX;

void sphereSource(OctreeOmniGridSolver3* solver) {
	Vector3f sphereCenter(0.5, 0.5, 0.07);
    Real smoke_radius = 0.06;
    Real initV = 1.0;
	auto& density = solver->getDensity();
	auto& heat = solver->getHeat();
	auto& dOctagonGrid = density->getOctagonGrid();
    auto& dTiltGrid = density->getTiltGrid();
    auto& hOctagonGrid = heat->getOctagonGrid();
    auto& hTiltGrid = heat->getTiltGrid();
	
	ThreadPool::parallelIterateGrid(dOctagonGrid->getParallelIteratorVec(), [&](const auto& iter){
		if ((dOctagonGrid->position(iter.levelCoor()) - sphereCenter).norm() < smoke_radius) {
			dOctagonGrid->set(iter.levelCoor(), 1);
			hOctagonGrid->set(iter.levelCoor(), 1);
		}
	});
	ThreadPool::parallelIterateGrid(dTiltGrid->getParallelIteratorVec(), [&](const auto& iter){
		if ((dTiltGrid->position(iter.levelCoor()) - sphereCenter).norm() < smoke_radius) {
			dTiltGrid->set(iter.levelCoor(), 1);
			hTiltGrid->set(iter.levelCoor(), 1);
		}
	});
	parallel_iterate_half_tilt(solver->getLayout(), [&](LCOOR_T lCoor) {
		if ((dTiltGrid->position(lCoor) - sphereCenter).norm() < smoke_radius) {
			dTiltGrid->set(lCoor, 1);
			hTiltGrid->set(lCoor, 1);
		}
	});

	auto& velocity = solver->getVelocity();
	auto& tMACGrid = velocity->getMACGrid();
	auto& tTiltGrid = velocity->getTiltGrid();
	ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const auto& iter){
		if ((tMACGrid->positionZ(iter.levelCoor()) - sphereCenter).norm() < smoke_radius) {
			tMACGrid->setZ(iter.levelCoor(), initV);
		}

		if ((tTiltGrid->positionCenter(iter.levelCoor()) - sphereCenter).norm() < smoke_radius) {
			tTiltGrid->set(iter.levelCoor(), Vector8f::Ones() * initV * INV_SQRT3);
		}
	});
	parallel_iterate_half_tilt(solver->getLayout(), [&](LCOOR_T lCoor) {
		auto pos = tTiltGrid->positionCenter(lCoor);
		if ((pos - sphereCenter).norm() < smoke_radius) {
			Vector8f nv = Vector8f::Zero();
			auto tNode = solver->getTiltEGrid()->get(lCoor);
			for (int i = 0; i < 4; ++i) {
				nv(HALF_TILT_SLOT_INDEX[tNode.half_direction][i]) = initV * INV_SQRT3;
			}
			if (tNode.half_direction % 3 == 2) {
				nv(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]) = initV;
			}
			tTiltGrid->set(lCoor, nv);
		}
	});

	bool static first = true;
	if (first) {
		density->updateGhost();
		heat->updateGhost();
		first = false;
	}
}

void cylinderSource(OctreeOmniGridSolver3* solver) {
	Vector2f sphereCenter(0.5, 0.5);
	Vector2f height(0.02, 0.15);
	REAL smoke_radius = 0.08;
	auto& density = solver->getDensity();
	auto& heat = solver->getHeat();
	auto& dOctagonGrid = density->getOctagonGrid();
    auto& dTiltGrid = density->getTiltGrid();
    auto& hOctagonGrid = heat->getOctagonGrid();
    auto& hTiltGrid = heat->getTiltGrid();
	
	ThreadPool::parallelIterateGrid(dOctagonGrid->getParallelIteratorVec(), [&](const auto& iter){
		auto pos = dOctagonGrid->position(iter.levelCoor());
		if (pos.z() > height[0] && pos.z() < height[1] && (Vector2f(pos.x(), pos.y()) - sphereCenter).norm() < smoke_radius) {
			dOctagonGrid->set(iter.levelCoor(), 1);
			hOctagonGrid->set(iter.levelCoor(), 1);
		}
	});
	ThreadPool::parallelIterateGrid(dTiltGrid->getParallelIteratorVec(), [&](const auto& iter){
		auto pos = dTiltGrid->position(iter.levelCoor());
		if (pos.z() > height[0] && pos.z() < height[1] && (Vector2f(pos.x(), pos.y()) - sphereCenter).norm() < smoke_radius) {
			dTiltGrid->set(iter.levelCoor(), 1);
			hTiltGrid->set(iter.levelCoor(), 1);
		}
	});
	parallel_iterate_half_tilt(solver->getLayout(), [&](LCOOR_T lCoor) {
		auto pos = dTiltGrid->position(lCoor);
		if (pos.z() > height[0] && pos.z() < height[1] && (Vector2f(pos.x(), pos.y()) - sphereCenter).norm() < smoke_radius) {
			dTiltGrid->set(lCoor, 1);
			hTiltGrid->set(lCoor, 1);
		}
	});

	REAL initV = 0.2;
	auto& velocity = solver->getVelocity();
	auto& tMACGrid = velocity->getMACGrid();
	auto& tTiltGrid = velocity->getTiltGrid();
	ThreadPool::parallelIterateGrid(tMACGrid->getParallelIteratorVec(), [&](const auto& iter){
		auto pos = tMACGrid->positionZ(iter.levelCoor());
		if (pos.z() > height[0] && pos.z() < height[1] && (Vector2f(pos.x(), pos.y()) - sphereCenter).norm() < smoke_radius) {
			tMACGrid->setZ(iter.levelCoor(), initV);
		}

		pos = tTiltGrid->positionCenter(iter.levelCoor());
		if (pos.z() > height[0] && pos.z() < height[1] && (Vector2f(pos.x(), pos.y()) - sphereCenter).norm() < smoke_radius) {
			tTiltGrid->set(iter.levelCoor(), Vector8f::Ones() * initV * INV_SQRT3);
		}
	});
	parallel_iterate_half_tilt(solver->getLayout(), [&](LCOOR_T lCoor) {
		auto pos = tTiltGrid->positionCenter(lCoor);
		if (pos.z() > height[0] && pos.z() < height[1] && (Vector2f(pos.x(), pos.y()) - sphereCenter).norm() < smoke_radius) {
			Vector8f nv = Vector8f::Zero();
			auto tNode = solver->getTiltEGrid()->get(lCoor);
			for (int i = 0; i < 4; ++i) {
				nv(HALF_TILT_SLOT_INDEX[tNode.half_direction][i]) = initV * INV_SQRT3;
			}
			if (tNode.half_direction % 3 == 2) {
				nv(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]) = initV;
			}
			tTiltGrid->set(lCoor, nv);
		}
	});

	bool static first = true;
	if (first) {
		density->updateGhost();
		heat->updateGhost();
		first = false;
	}
}

void smoke_sim(const sim_args& args) {
	std::string res_str = std::to_string(RESX) + "_" + std::to_string(RESY) + "_" + std::to_string(RESZ);
	std::string solver_name = "octree_omni_smoke_sphere_" + res_str;
	OctreeOmniGridSolver3 solver(res, delta_h, solver_name);
	solver.sceneSetup = sphereSource;
	// solver.sceneSetup = cylinderSource;
	solver.isDumpDensity = args.dumpDensity;
	solver.isDumpDensity = true;
	solver.setMaxCfl(2);
	solver.numOctreeLayer = 2;
	solver.splitMapSetup = [&](OctreeOmniGridSolver3* solver) {
		auto& splitMap = solver->getSplitMap();
#if 1
        splitMap->fill(1);  // default all cell is in level 1
        auto boundarySdf = std::make_shared<CellCenteredScalarGrid3>(splitMap->resolution(), splitMap->gridSpacing(), 10000);
        auto cylinder0 = std::make_shared<Cylinder3>(0.08, 0.5, Transform3(Vector3f(0.5, 0.5, 0.2), Quaternion(Vector3f(1, 0, 0), PI / 2)));
        // // cylinder1 = make_shared<Cylinder3>(0.15, 0.35, Transform3(Vector3f(0.5, 0.5, 0.55), Quaternion(Vector3f(1, 0,
        // // 0), PI/2)));
        auto sphere = std::make_shared<Sphere3>(0.2, Transform3(Vector3f(0.5, 0.5, 0.50), Quaternion()));
        // auto sphere = std::make_shared<Sphere3>(0.18, Transform3(Vector3f(0.5, 0.5, 0.50), Quaternion()));
        auto cylinder0_implicit = std::make_shared<ShapeToImplicit3>(cylinder0);
        auto sphere_implicit = std::make_shared<ShapeToImplicit3>(sphere);
        // // cylinder1_implicit = std::make_shared<ShapeToImplicit3>(cylinder1);
        boundarySdf->fill(10000);
        boundarySdf->parallelForEach([&](const Size3& idx) {
            auto pos = boundarySdf->position(idx);
        //     // REAL cl = std::min(sphere_implicit->signedDistance(pos), rect0_implicit->signedDistance(pos));
        //     // REAL cl = std::min(cylinder0_implicit->signedDistance(pos), cylinder1_implicit->signedDistance(pos));
            REAL cl = sphere_implicit->signedDistance(pos);
            // REAL cl = std::min(cylinder0_implicit->signedDistance(pos), sphere_implicit->signedDistance(pos));
            boundarySdf->set(idx, cl);
        });

        refineSplitMapByObstacleSdf(splitMap, 1, boundarySdf, Vector2f(-1000, 0.02));  // level 1 range
#elif 0
		splitMap->fill(3);		// default all cell is in level 3
		//Rectangle3Ptr rectangle_top = std::make_shared<Rectangle3>(Vector3f(0.5, 0.5, 0.05),
			//Transform3(Vector3f(0.25, 0.25, 0.3)));
		//Rectangle3Ptr rectangle_top_coarse = std::make_shared<Rectangle3>(Vector3f(0.5, 0.5, 0.08),
			//Transform3(Vector3f(0.25, 0.25, 0.25)));

		Sphere3Ptr sphere_top = std::make_shared<Sphere3>(0.10, Transform3(Vector3f(0.25, 0.25, 0.36)));
		Cylinder3Ptr cylinder_top = std::make_shared<Cylinder3>(0.25, 0.02, Transform3(Vector3f(0.25, 0.25, 0.27), Quaternion(Vector3f(1, 0, 0), PI / 2)));
		Cylinder3Ptr cylinder_top_coarse = std::make_shared<Cylinder3>(0.25, 0.1, Transform3(Vector3f(0.25, 0.25, 0.22), Quaternion(Vector3f(1, 0, 0), PI / 2)));
		Cylinder3Ptr cylinder_smoke = std::make_shared<Cylinder3>(0.055, 0.18, Transform3(Vector3f(0.25, 0.25, 0.25), Quaternion(Vector3f(1, 0, 0), PI / 2)));

		auto implicit0 = std::make_shared<ShapeToImplicit3>(sphere_top);
		auto implicit1 = std::make_shared<ShapeToImplicit3>(cylinder_smoke);
		auto implicit2 = std::make_shared<ShapeToImplicit3>(cylinder_top);
		auto implicit3 = std::make_shared<ShapeToImplicit3>(cylinder_top_coarse);
		//auto cylinder_implicit2 = std::make_shared<ShapeToImplicit3>(cylinder_boundary);
		auto boundarySdf = make_shared<CellCenteredScalarGrid3>(splitMap->resolution(), splitMap->gridSpacing(), 10000);
		boundarySdf->parallelForEach([&](int xIdx, int yIdx, int zIdx) {
			auto pos = boundarySdf->position(xIdx, yIdx, zIdx);
			REAL cl = std::min(implicit0->signedDistance(pos), implicit3->signedDistance(pos));
			cl = std::min(implicit1->signedDistance(pos), cl);
			boundarySdf->set(xIdx, yIdx, zIdx, cl);
		});
		
		refineSplitMapByObstacleSdf(splitMap, 3, boundarySdf, Vector2f(-1000, 0.4)); // level 2 range
		refineSplitMapByObstacleSdf(splitMap, 2, boundarySdf, Vector2f(-1000, 0.18)); // level 1 range

		boundarySdf->fill(10000);
		boundarySdf->parallelForEach([&](int xIdx, int yIdx, int zIdx) {
			auto pos = boundarySdf->position(xIdx, yIdx, zIdx);
			REAL cl = std::min(implicit0->signedDistance(pos), implicit2->signedDistance(pos));
			cl = std::min(implicit1->signedDistance(pos), cl);
			boundarySdf->set(xIdx, yIdx, zIdx, cl);
		});
		refineSplitMapByObstacleSdf(splitMap, 1, boundarySdf, Vector2f(-1000, 0.04)); // level 0 range
#else
		// return;
		splitMap->fill(3);		// default all cell is in level 3
		auto cylinder0 = make_shared<Cylinder3>(0.1, 0.4, Transform3(Vector3f(0.5, 0.5, 0.2), Quaternion(Vector3f(1, 0, 0), PI/2)));
		auto cylinder1 = make_shared<Cylinder3>(0.15, 0.40, Transform3(Vector3f(0.5, 0.5, 0.45), Quaternion(Vector3f(1, 0, 0), PI/2)));
		auto cylinder2 = make_shared<Cylinder3>(0.1, 0.4, Transform3(Vector3f(0.5, 0.5, 0.6), Quaternion(Vector3f(1, 0, 0), PI/2)));
		auto cylinder0_implicit = std::make_shared<ShapeToImplicit3>(cylinder0);
		auto cylinder1_implicit = std::make_shared<ShapeToImplicit3>(cylinder1);
		auto cylinder2_implicit = std::make_shared<ShapeToImplicit3>(cylinder2);
		auto boundarySdf = make_shared<CellCenteredScalarGrid3>(splitMap->resolution(), splitMap->gridSpacing(), 10000);
		boundarySdf->parallelForEach([&](int xIdx, int yIdx, int zIdx) {
			auto pos = boundarySdf->position(xIdx, yIdx, zIdx);
			// REAL cl = std::min(sphere_implicit->signedDistance(pos), rect0_implicit->signedDistance(pos));
			REAL cl = std::min(cylinder0_implicit->signedDistance(pos), cylinder1_implicit->signedDistance(pos));
			cl = std::min(cl, cylinder2_implicit->signedDistance(pos));
			boundarySdf->set(xIdx, yIdx, zIdx, cl);
		});


		// refineSplitMapByObstacleSdf(splitMap, 4, boundarySdf, Vector2f(-1000, 0.30)); // level 2 range
		refineSplitMapByObstacleSdf(splitMap, 3, boundarySdf, Vector2f(-1000, 0.25)); // level 2 range
		
		cylinder0 = make_shared<Cylinder3>(0.1, 0.4, Transform3(Vector3f(0.5, 0.5, 0.2), Quaternion(Vector3f(1, 0, 0), PI/2)));
		// cylinder1 = make_shared<Cylinder3>(0.15, 0.35, Transform3(Vector3f(0.5, 0.5, 0.55), Quaternion(Vector3f(1, 0, 0), PI/2)));
		auto sphere = std::make_shared<Sphere3>(0.2, Transform3(Vector3f(0.5, 0.5, 0.50), Quaternion()));
		cylinder0_implicit = std::make_shared<ShapeToImplicit3>(cylinder0);
		auto sphere_implicit = std::make_shared<ShapeToImplicit3>(sphere);
		// cylinder1_implicit = std::make_shared<ShapeToImplicit3>(cylinder1);
		boundarySdf->fill(10000);
		boundarySdf->parallelForEach([&](int xIdx, int yIdx, int zIdx) {
			auto pos = boundarySdf->position(xIdx, yIdx, zIdx);
			// REAL cl = std::min(sphere_implicit->signedDistance(pos), rect0_implicit->signedDistance(pos));
			// REAL cl = std::min(cylinder0_implicit->signedDistance(pos), cylinder1_implicit->signedDistance(pos));
			REAL cl = std::min(cylinder0_implicit->signedDistance(pos), sphere_implicit->signedDistance(pos));
			boundarySdf->set(xIdx, yIdx, zIdx, cl);
		});
		
		refineSplitMapByObstacleSdf(splitMap, 2, boundarySdf, Vector2f(-1000, 0.1)); // level 1 range

		cylinder0 = make_shared<Cylinder3>(0.06, 0.3, Transform3(Vector3f(0.5, 0.5, 0.17), Quaternion(Vector3f(1, 0, 0), PI/2)));
		cylinder1 = make_shared<Cylinder3>(0.08, 0.28, Transform3(Vector3f(0.5, 0.5, 0.3), Quaternion(Vector3f(1, 0, 0), PI/2)));
		sphere = std::make_shared<Sphere3>(0.15, Transform3(Vector3f(0.5, 0.5, 0.46), Quaternion()));
		// cylinder2 = make_shared<Cylinder3>(0.18, 0.35, Transform3(Vector3f(0.5, 0.5, 0.4), Quaternion(Vector3f(1, 0, 0), PI/2)));
		cylinder0_implicit = std::make_shared<ShapeToImplicit3>(cylinder0);
		cylinder1_implicit = std::make_shared<ShapeToImplicit3>(cylinder1);
		sphere_implicit = std::make_shared<ShapeToImplicit3>(sphere);
		boundarySdf->fill(10000);
		boundarySdf->parallelForEach([&](int xIdx, int yIdx, int zIdx) {
			auto pos = boundarySdf->position(xIdx, yIdx, zIdx);
			// REAL cl = std::min(sphere_implicit->signedDistance(pos), rect0_implicit->signedDistance(pos));
			// REAL cl = std::min(cylinder0_implicit->signedDistance(pos), cylinder1_implicit->signedDistance(pos));
			REAL cl = std::min(sphere_implicit->signedDistance(pos), cylinder1_implicit->signedDistance(pos));
			// cl = std::min(cl, sphere_implicit->signedDistance(pos));
			boundarySdf->set(xIdx, yIdx, zIdx, cl);
		});

		refineSplitMapByObstacleSdf(splitMap, 1, boundarySdf, Vector2f(-1000, 0.05)); // level 0 range
#endif
		Viewer::grid2PNGImage("./" + solver_name + "/b_fine_scheme.png",
			Size2(splitMap->resolution().x(), splitMap->resolution().z()), [&](int xIdx, int zIdx)->REAL {
			Vector3f pos = Vector3f(xIdx, splitMap->resolution().y() / 2.0, zIdx) * splitMap->gridSpacing();
			return boundarySdf->sample(pos) > 0 ? 0.0 : 1.0;
			// return (boundarySdf->sample(pos) + 1) * 0.5;
		});

		// for split map debug 
		Viewer::grid2PNGImage("./" + solver_name + "/splitMap.png",
			Size2(splitMap->resolution().x(), splitMap->resolution().z()), [&](int xIdx, int zIdx)->REAL {
			// Vector3f pos = Vector3f(xIdx, splitMap->resolution().y()/2.0, zIdx) * splitMap->gridSpacing();
			int p = splitMap->get(xIdx, splitMap->resolution().y()/2.0, zIdx);
			return p / 3.0;
		});
	};

	//solver.isRestart = true;
	solver.initialize();
	// if (!solver.restart(21, "./octree_omni_smoke_128_128_128/layout.bin", "./octree_omni_smoke_128_128_128/boundary.bin", "./octree_omni_smoke_128_128_128/tilt/t_0021.bin",
	// 	"./octree_omni_smoke_128_128_128/density/d_0021.bin", "./octree_omni_smoke_128_128_128/heat/h_0021.bin", "./octree_omni_smoke_128_128_128/velocity/v_0021.bin")) {
	// 	spdlog::info("restart failed");
	// }
	try {
		for (int i = 0; i < args.num_frame; ++i) {
			if(args.showGUI) { /* update viewer */ }
			spdlog::info("current {} frame", i);
			if(i % 10 == 0) { BENCHMARK_REPORT(); }
			solver.update(TIME_STEP);
		}
	} catch(std::runtime_error& e) {
		spdlog::error("main caught exception. {}", e.what());
	}
	
	BENCHMARK_REPORT();
	int c = 0;
    while(solver.getACounter().count != 0) {
        if(c != solver.getACounter().count) {
            c = solver.getACounter().count;
            spdlog::warn("Waiting result saving to complete. Remaining: {}", solver.getACounter().count);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
	spdlog::info("All done");
}

int main(int argc, char *argv[]) {
	auto args = getSimArgs(argc, argv);
	smoke_sim(args);
    return 0;
}
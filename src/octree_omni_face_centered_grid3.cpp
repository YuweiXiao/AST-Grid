#include "octree_omni_face_centered_grid3.h"
#include "util.h"

namespace Omni {

void OctreeOmniFaceCenteredGrid3::setAll(REAL v) {
	MACGridPtr->setAll(v);
	tiltGridPtr->setAll(Vector8f::Ones() * v);
}

void OctreeOmniFaceCenteredGrid3::interpolateFromDualGraph(OctreeOmniVector3Grid3Ptr ptr) {
	ThreadPool::parallelIterateGrid(MACGridPtr->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
		for(int axis = 0; axis < 3; ++axis) {
			auto pos = MACGridPtr->position(axis, iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
			MACGridPtr->set(axis, iter.levelCoor(), ptr->sample(pos)[axis]);
		}
	});
	ThreadPool::parallelIterateGrid(tiltGridPtr->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
		if(tiltE->get(iter.levelCoor()).is_closed)
			return;
		Vector8f v = tiltGridPtr->get(iter.levelCoor());
		for(int k = 0; k < 8; ++k) {
			auto pos = tiltGridPtr->positionTilt(iter.levelCoor(), k);
			v(k) = ptr->sample(pos).dot(TILT3_UNIT_DIR_UP[k]);
		}
		tiltGridPtr->set(iter.levelCoor(), v);
	});
	// TODO half-tilt
}


void OctreeOmniFaceCenteredGrid3::updateDualGrid() {
	BENCHMARK_SCOPED_TIMER_SECTION t("update dual grid");
    auto& dualTiltGrid = dualGridPtr->getTiltGrid();
    auto& dualOctagonGrid = dualGridPtr->getOctagonGrid();
	auto& tiltE = getTiltEGridPtr();

	{
	BENCHMARK_SCOPED_TIMER_SECTION t_tilt("update dual grid - tilt");
	ThreadPool::parallelIterateGrid(dualTiltGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
        if(tiltE->get(iter.levelCoor()).is_closed) {	// skip half-tilt for now
            dualTiltGrid->set(iter.levelCoor(), Vector3f::Zero());
        } else {
            const Vector8f& gVal = tiltGridPtr->get(iter.levelCoor());
			Vector4f v = Vector4f(gVal(0) + gVal(6), gVal(1) +gVal(7), gVal(2) + gVal(4), gVal(3) + gVal(5)) * 0.5;
			// ASSERT(!std::isnan(nv.x()) && !std::isnan(nv.y()) && !std::isnan(nv.z()), "dual velocity is nan");
            dualTiltGrid->set(iter.levelCoor(), Constants::OMNI_DUAL_GRAPH_TILT_EIGEN_MATRIX3 * v);
        }
	});}

	ThreadPool::parallelIterateGrid(dualOctagonGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
        if(!dualOctagonGrid->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
            dualOctagonGrid->set(iter.levelCoor(), Vector3f::Zero());
            return;
        }

		REAL delta_h = layout->levelGridSpacing(iter.level);
		REAL sqrt3_delta_h = SQRT3 * delta_h;
		REAL half_sqrt3_delta_h = sqrt3_delta_h * 0.5;

		int n_joint[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
		for (int i = 0; i < 8; i++) {
			const auto& node = tiltE->get(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[i](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[i](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[i](2));
			n_joint[i] = (node.is_T_joint && fabs(node.e - delta_h) < EPSILON) || node.is_half_tilt;
		}

		Vector7f v = Vector7f::Zero();
		int cnt[7] = { 0, 0, 0, 0, 0, 0, 0 };

		// left (x-)
		if (!((n_joint[1] && n_joint[6]) || (n_joint[2] && n_joint[5]))) {
			v[0] += MACGridPtr->getX(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
			cnt[0]++;
		}
		// right (x+)
		if (!((n_joint[0] && n_joint[7]) || (n_joint[3] && n_joint[4]))) {
			v[0] += MACGridPtr->getX(iter.level, iter.xIdx() + 1, iter.yIdx(), iter.zIdx());
			cnt[0]++;
		}
		// front (y-)
		if (!((n_joint[2] && n_joint[7]) || (n_joint[3] && n_joint[6]))) {
			v[1] += MACGridPtr->getY(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
			cnt[1]++;
		}
		// back (y+)
		if (!((n_joint[0] && n_joint[5]) || (n_joint[1] && n_joint[4]))) {
			v[1] += MACGridPtr->getY(iter.level, iter.xIdx(), iter.yIdx() + 1, iter.zIdx());
			cnt[1]++;
		}
		// down (z-)
		if (!((n_joint[4] && n_joint[6]) || (n_joint[5] && n_joint[7]))) {
			v[2] += MACGridPtr->getZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
			cnt[2]++;
		}
		// up (z+)
		if (!((n_joint[0] && n_joint[2]) || (n_joint[1] && n_joint[3]))) {
			v[2] += MACGridPtr->getZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx() + 1);
			cnt[2]++;
		}

		for(int i = 0; i < 4; ++i) {
			int slot = i;
			int reSlot = ((i&1)?8:6) - i;	// opposite index of slot
			int vSlot = i+3;
			const auto& node = tiltE->get(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[i](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[i](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[i](2));
			const auto& reNode = tiltE->get(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[reSlot](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[reSlot](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[reSlot](2));
			if (!node.is_closed || node.is_half_tilt) {
				v[vSlot] += tiltGridPtr->get(iter.level, 
						iter.xIdx() + OCT_NB_TILT_DXDYDZ[slot].x(), 
						iter.yIdx() + OCT_NB_TILT_DXDYDZ[slot].y(), 
						iter.zIdx() + OCT_NB_TILT_DXDYDZ[slot].z())(reSlot);
				cnt[vSlot]++;
			}
			if (!reNode.is_closed || reNode.is_half_tilt) {
				REAL tv = tiltGridPtr->get(iter.level, 
						iter.xIdx() + OCT_NB_TILT_DXDYDZ[reSlot].x(), 
						iter.yIdx() + OCT_NB_TILT_DXDYDZ[reSlot].y(), 
						iter.zIdx() + OCT_NB_TILT_DXDYDZ[reSlot].z())(slot);
				v[vSlot] = cnt[vSlot] == 0? tv : (((v[vSlot] * (half_sqrt3_delta_h - reNode.e * INV_SQRT3)
					+ tv * (half_sqrt3_delta_h - node.e * INV_SQRT3))
					/ (sqrt3_delta_h - node.e * INV_SQRT3 - reNode.e * INV_SQRT3)));
				cnt[vSlot]++;
			}
		}

		for (int i = 0; i < 3; i++) {
			v[i] = v[i] / (cnt[i] + EPSILON);
		}
		Vector3f finalV = Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[cnt[0] > 0][cnt[1] > 0][cnt[2] > 0]
						[cnt[3] > 1 || n_joint[0] || n_joint[6]] [cnt[4] > 1 || n_joint[1] || n_joint[7]]
						[cnt[5] > 1 || n_joint[2] || n_joint[4]] [cnt[6] > 1 || n_joint[3] || n_joint[5]] * v;

		// for (int k = 0; k < 3; k++) {
		// 	for (int i = 0; i < 7; i++) {
		// 		finalV(k) += v[i] *
		// 			Constants::OCTREE_OMNI_DUAL_GRAPH_MATRIX3[cnt[0] > 0][cnt[1] > 0][cnt[2] > 0]
				// [cnt[3] > 1 || tilts[0].is_T_joint || tilts[6].is_T_joint || tilts[0].is_half_tilt || tilts[6].is_half_tilt]
				// [cnt[4] > 1 || tilts[1].is_T_joint || tilts[7].is_T_joint || tilts[1].is_half_tilt || tilts[7].is_half_tilt]
				// [cnt[5] > 1 || tilts[2].is_T_joint || tilts[4].is_T_joint || tilts[2].is_half_tilt || tilts[4].is_half_tilt]
				// [cnt[6] > 1 || tilts[3].is_T_joint || tilts[5].is_T_joint || tilts[3].is_half_tilt || tilts[5].is_half_tilt]
		// 		[cnt[3] > 1 || n_joint[0] || n_joint[6]] [cnt[4] > 1 || n_joint[1] || n_joint[7]]
		// 		[cnt[5] > 1 || n_joint[2] || n_joint[4]] [cnt[6] > 1 || n_joint[3] || n_joint[5]]
		// 		[k][i];
		// 	}
		// }

		// if(iter.level == 0 && iter.xIdx() == 62 && iter.yIdx() == 31 && iter.zIdx() == 1) {
		// 	spdlog::info("here: {}", finalV.transpose());
		// 	spdlog::info("x: {} {}", MACGridPtr->getX(0, 62,31,1), MACGridPtr->getX(0, 63,31,1));
		// 	spdlog::info("v:{} {} {} {} {} {} {}", v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
		// 	spdlog::info("cnt:{} {} {} {} {} {} {}", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], cnt[5], cnt[6]);
		// }

		ASSERT(!std::isnan(finalV.x()) && !std::isnan(finalV.y()) && !std::isnan(finalV.z()), "dual velocity is nan");
		dualOctagonGrid->set(iter.levelCoor(), finalV);
    });

	{
	BENCHMARK_SCOPED_TIMER_SECTION t_half_tilt("update dual grid - half-tilt");
	parallel_iterate_half_tilt(layout, [&](LCOOR_T lCoor){
	// iterate_half_tilt(lCoor, layout) {
		const auto& tNode = tiltE->get(lCoor);
		const auto& gVal = tiltGridPtr->get(lCoor);
		Vector7f v = Vector7f::Zero();
		int axis = tNode.half_direction % 3;	// x:0, y:1, z:2

		v[axis] = gVal(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]);
		for(int i = 0; i < 4; ++i) {
			int slot = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
			int reSlot = ((slot&1)?8:6) - slot;
			int idx = slot < reSlot ? slot : reSlot;
			v[idx+3] = gVal(slot);
		}

		Vector3f finalV = Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[axis == 0][axis == 1][axis == 2][1][1][1][1] * v;

		// for (int k = 0; k < 3; k++) {
		// 	for (int i = 0; i < 7; i++) {
		// 		finalV(k) += v[i] * 
		// 			Constants::OCTREE_OMNI_DUAL_GRAPH_MATRIX3[axis == 0][axis == 1][axis == 2][1][1][1][1][k][i];
		// 	}
		// }
		ASSERT(!std::isnan(finalV.x()) && !std::isnan(finalV.y()) && !std::isnan(finalV.z()), "dual velocity is nan");
		dualTiltGrid->set(lCoor, finalV);
	});}
}


// void OctreeOmniFaceCenteredGrid3::updateDualGrid(const CellCenteredScalarGrid3Ptr& levelSet) {
//     auto& dualTiltGrid = dualGridPtr->getTiltGrid();
//     auto& dualOctagonGrid = dualGridPtr->getOctagonGrid();
//     REAL levelSetOffset = 0.2;
	
// 	ThreadPool::parallelIterateGrid(dualTiltGrid->getIterator(), [&](const OctreeGridIterator3& iter) {
//         if(getTiltEGridPtr()->get(iter.levelCoor()).is_closed	// skip half-tilt for now
// 			|| !isInsideSdf(layout->levelGridSpacing(iter.level) * levelSetOffset + levelSet->sample(dualTiltGrid->position(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())))) {	
//             dualTiltGrid->set(iter.levelCoor(), Vector3f::Zero());
//         } else {
//             const Vector8f& gVal = tiltGridPtr->get(iter.levelCoor());
// 			Vector4f v;
// 			for(int i = 0; i < 4; ++i) {
// 				int revSlot = (i&1?8:6) - i;
// 				REAL valid_n0 = isInsideSdf(levelSet->sample(tiltGridPtr->positionTilt(iter.levelCoor(), i))) ? 1.0 : EPSILON;
// 				REAL valid_n1 = isInsideSdf(levelSet->sample(tiltGridPtr->positionTilt(iter.levelCoor(), revSlot))) ? 1.0 : EPSILON;
// 				v(i) = (gVal(i) * valid_n0 + gVal(revSlot) * valid_n1) / (valid_n1 + valid_n0);
// 			}
//             Vector3f nv(0, 0, 0);
//             nv(0) = INV_SQRT3 * v(0) - INV_SQRT3 * v(1) - INV_SQRT3 * v(2) + INV_SQRT3 * v(3);
//             nv(1) = INV_SQRT3 * v(0) + INV_SQRT3 * v(1) - INV_SQRT3 * v(2) - INV_SQRT3 * v(3);
//             nv(2) = INV_SQRT3 * v(0) + INV_SQRT3 * v(1) + INV_SQRT3 * v(2) + INV_SQRT3 * v(3);
//             nv *= 3.0/4.0;
// 			// if(iter.level == 1 && iter.xIdx() == 32 && iter.yIdx() == 32 && iter.zIdx() == 53) {
// 			// 	spdlog::info("hihihi {} {} {} {}", v(0), v(1), v(2), v(3));
// 			// 	spdlog::info("\t nV: {}", nv.transpose());
// 			// 	for(int k = 0; k < 8; ++k) {
// 			// 		std::cout<<isInsideSdf(levelSet->sample(tiltGridPtr->positionTilt(iter.levelCoor(), k)))<<' ';
// 			// 	}
// 			// 	std::cout<<std::endl;
// 			// }
// 			ASSERT(!std::isnan(nv.x()) && !std::isnan(nv.y()) && !std::isnan(nv.z()), "tilt cell: dual velocity is nan");
//             dualTiltGrid->set(iter.levelCoor(), nv);
//         }
//     });

// 	ThreadPool::parallelIterateGrid(dualOctagonGrid->getIterator(), [&](const OctreeGridIterator3& iter) {
//         if(!dualOctagonGrid->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())
// 			|| !isInsideSdf(layout->levelGridSpacing(iter.level) * levelSetOffset + levelSet->sample(dualOctagonGrid->position(iter.levelCoor())))) {
//             dualOctagonGrid->set(iter.levelCoor(), Vector3f::Zero());
//             return;
//         }

// 		REAL delta_h = layout->levelGridSpacing(iter.level);
// 		REAL sqrt3_delta_h = SQRT3 * delta_h;
// 		REAL half_sqrt3_delta_h = sqrt3_delta_h * 0.5;

// 		TiltENode tilts[8];
// 		int n_joint[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
// 		bool isInsideLevel[8];
// 		for (int i = 0; i < 8; i++) {
// 			LCoor3 dp(iter.level, iter.xIdx() + OCT_NB_TILT_DXDYDZ[i](0), iter.yIdx() + OCT_NB_TILT_DXDYDZ[i](1), iter.zIdx() + OCT_NB_TILT_DXDYDZ[i](2));
// 			tilts[i] = getTiltEGridPtr()->get(dp.levelCoor());
// 			isInsideLevel[i] = isInsideSdf(levelSet->sample(tiltGridPtr->positionTilt(dp.levelCoor(), reverseTiltIdx(i))));
// 			n_joint[i] = (tilts[i].is_T_joint && fabs(tilts[i].e - delta_h) < EPSILON) || tilts[i].is_half_tilt;
// 		}

// 		REAL v[7] = { 0, 0, 0, 0, 0, 0, 0 };
// 		int cnt[7] = { 0, 0, 0, 0, 0, 0, 0 };
// 		static int OCT_NB_IDX[][4] = {{0, 7, 3, 4}, {1, 6, 2, 5},	// +x, -x
// 									{0, 5, 1, 4}, {2, 7, 3, 6},		// +y, -y
// 									{0, 2, 1, 3}, {4, 6, 5, 7}};	// +z, -z
// 		for(int axis = 0; axis < 3; ++axis) {
// 			Vector3f deltaPos(0, 0, 0); 
// 			deltaPos[axis] = delta_h * 0.5;
// 			int idx0 = axis << 1; 
// 			int idx1 = (axis << 1) + 1;
// 			auto pos = dualOctagonGrid->position(iter.levelCoor());
// 			if (!((n_joint[OCT_NB_IDX[idx0][0]] && n_joint[OCT_NB_IDX[idx0][1]]) || (n_joint[OCT_NB_IDX[idx0][2]] && n_joint[OCT_NB_IDX[idx0][3]]))
// 				&& isInsideSdf(levelSet->sample(pos + deltaPos))) {
// 				v[axis] += MACGridPtr->get(axis, iter.level, iter.xIdx()+OCT_NB_FACE_DXDYDZ[idx0].x(), 
// 								iter.yIdx()+OCT_NB_FACE_DXDYDZ[idx0].y(), iter.zIdx()+OCT_NB_FACE_DXDYDZ[idx0].z());
// 				cnt[axis]++;
// 			}
// 			if (!((n_joint[OCT_NB_IDX[idx1][0]] && n_joint[OCT_NB_IDX[idx1][1]]) || (n_joint[OCT_NB_IDX[idx1][2]] && n_joint[OCT_NB_IDX[idx1][3]]))
// 				&& isInsideSdf(levelSet->sample(pos - deltaPos))) {
// 				v[axis] += MACGridPtr->get(axis, iter.level, iter.xIdx()+OCT_NB_FACE_DXDYDZ[idx1].x(), 
// 								iter.yIdx()+OCT_NB_FACE_DXDYDZ[idx1].y(), iter.zIdx()+OCT_NB_FACE_DXDYDZ[idx1].z());
// 				cnt[axis]++;
// 			}
// 			v[axis] /= (cnt[axis] + EPSILON);
// 		}

// 		for(int i = 0; i < 4; ++i) {
// 			int slot = i;
// 			int reSlot = reverseTiltIdx(i);	// opposite index of slot
// 			int vSlot = i+3;
// 			if ((!tilts[i].is_closed || tilts[i].is_half_tilt) &&  isInsideLevel[i]) {
// 				v[vSlot] += tiltGridPtr->get(iter.level, 
// 						iter.xIdx() + OCT_NB_TILT_DXDYDZ[slot].x(), 
// 						iter.yIdx() + OCT_NB_TILT_DXDYDZ[slot].y(), 
// 						iter.zIdx() + OCT_NB_TILT_DXDYDZ[slot].z())(reSlot);
// 				cnt[vSlot]++;
// 			}
// 			if ((!tilts[reSlot].is_closed || tilts[reSlot].is_half_tilt) && isInsideLevel[reSlot]) {
// 				REAL tv = tiltGridPtr->get(iter.level, 
// 						iter.xIdx() + OCT_NB_TILT_DXDYDZ[reSlot].x(), 
// 						iter.yIdx() + OCT_NB_TILT_DXDYDZ[reSlot].y(), 
// 						iter.zIdx() + OCT_NB_TILT_DXDYDZ[reSlot].z())(slot);
// 				v[vSlot] = cnt[vSlot] == 0? tv : (((v[vSlot] * (half_sqrt3_delta_h - tilts[reSlot].e * INV_SQRT3)
// 					+ tv * (half_sqrt3_delta_h - tilts[slot].e * INV_SQRT3))
// 					/ (sqrt3_delta_h - tilts[slot].e * INV_SQRT3 - tilts[reSlot].e * INV_SQRT3)));
// 				cnt[vSlot]++;
// 			}	
// 		}

// 		Vector3f finalV(0, 0, 0);
// 		for (int k = 0; k < 3; k++) {
// 			for (int i = 0; i < 7; i++) {
// 				finalV(k) += v[i] *
// 					Constants::OCTREE_OMNI_DUAL_GRAPH_MATRIX3[cnt[0] > 0][cnt[1] > 0][cnt[2] > 0]
// 						[cnt[3] > 1 || n_joint[0] || n_joint[6] || (cnt[3]==1 && (isInsideLevel[0]^isInsideLevel[6]))]
// 						[cnt[4] > 1 || n_joint[1] || n_joint[7] || (cnt[4]==1 && (isInsideLevel[1]^isInsideLevel[7]))]
// 						[cnt[5] > 1 || n_joint[2] || n_joint[4] || (cnt[5]==1 && (isInsideLevel[2]^isInsideLevel[4]))]
// 						[cnt[6] > 1 || n_joint[3] || n_joint[5] || (cnt[6]==1 && (isInsideLevel[3]^isInsideLevel[5]))]
// 						[k][i];
// 			}
// 		}
// 		// if(iter.level == 1 && iter.level == 1 && iter.xIdx() == 32 && iter.yIdx() == 32 && iter.zIdx() == 53) {
// 		// 	spdlog::error("octagon!: idx:{} {} {} {}", iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
// 		// 	spdlog::info("\t{} {} {} {} {} {} {}", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], cnt[5], cnt[6]);
// 		// 	spdlog::info("\t {} {} {} {} {} {} {}", v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
// 		// 	spdlog::info("{} {}", MACGridPtr->getZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx()), MACGridPtr->getZ(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx()+1));
// 		// }
// 		ASSERT(!std::isnan(finalV.x()) && !std::isnan(finalV.y()) && !std::isnan(finalV.z()), "octagon cell: dual velocity is nan");
// 		dualOctagonGrid->set(iter.levelCoor(), finalV);
//     });

// 	iterate_half_tilt(lCoor, layout) {
// 		if(!isInsideSdf(layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)) * 2 * levelSetOffset + levelSet->sample(tiltGridPtr->positionCenter(lCoor)))) {
// 			dualTiltGrid->set(lCoor, Vector3f::Zero());
// 			continue;
// 		}
// 		const auto& tNode = getTiltEGridPtr()->get(lCoor);
// 		const auto& gVal = tiltGridPtr->get(lCoor);
// 		REAL v[7] = {0, 0, 0, 0, 0, 0, 0};
// 		bool isInsideLevel[4];	// for 4 tilt axis
// 		int axis = tNode.half_direction % 3;	// x:0, y:1, z:2

// 		v[axis] = gVal(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]);
// 		for(int i = 0; i < 4; ++i) {
// 			int slot = HALF_TILT_SLOT_INDEX[tNode.half_direction][i];
// 			int reSlot = ((slot&1)?8:6) - slot;
// 			int idx = slot < reSlot ? slot : reSlot;
// 			v[idx+3] = gVal(slot);
// 			isInsideLevel[idx] = isInsideSdf(levelSet->sample(tiltGridPtr->positionTilt(lCoor, slot)));
// 		}

// 		Vector3f finalV(0, 0, 0);
// 		for (int k = 0; k < 3; k++) {
// 			for (int i = 0; i < 7; i++) {
// 				finalV(k) += v[i] * 
// 					Constants::OCTREE_OMNI_DUAL_GRAPH_MATRIX3[axis == 0][axis == 1][axis == 2]
// 						[isInsideLevel[0]][isInsideLevel[1]][isInsideLevel[2]][isInsideLevel[3]][k][i];
// 			}
// 		}
// 		// if(UNPACK_LEVEL3(lCoor) == 0 && UNPACK_X3(lCoor) == 65 && UNPACK_Y3(lCoor) == 65 && UNPACK_Z3(lCoor) == 106) {
// 		// 	spdlog::info("tilt! : {} {} {} {} {} {} {}", v[0], v[1],v[2],v[3],v[4],v[5],v[6]);
// 		// 	spdlog::info("finalV: {}", finalV.transpose());
// 		// }
// 		ASSERT(!std::isnan(finalV.x()) && !std::isnan(finalV.y()) && !std::isnan(finalV.z()), "half-tilt, dual velocity is nan");
// 		dualTiltGrid->set(lCoor, finalV);
// 	}
// }


void OctreeOmniFaceCenteredGrid3::updateGhostFromHalfTilt() {
	iterate_half_tilt(lCoor, layout) {
		const auto& tNode = getTiltEGridPtr()->get(lCoor);
		const auto& gVal = tiltGridPtr->get(lCoor);
		int axis = tNode.half_direction % 3;
		// if(UNPACK_LEVEL3(lCoor) == 0 && UNPACK_X3(lCoor) == 65 && UNPACK_Y3(lCoor) == 65 && UNPACK_Z3(lCoor) == 106) {
		// 	spdlog::error("hello ??? {}", gVal(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]));
		// }
		MACGridPtr->set(axis, UNPACK_LEVEL3(lCoor)+1,
			UNPACK_X3(lCoor)>>1, 
			UNPACK_Y3(lCoor)>>1, 
			UNPACK_Z3(lCoor)>>1, 
			gVal(HALF_TILT_SLOT_INDEX[tNode.half_direction][4]));
	}
}

} // end of namespace Omni
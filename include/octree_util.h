#pragma once
#include "general.h"
#include "octree_vertex_centered_grid3.h"
#include "octree_cell_centered_grid3.h"
#include "bitmap.h"
#include <vector>

namespace Omni {

void constrainOctreeTiltE3(OctreeTiltENodeGrid3Ptr& tiltEGrid, bool recordHalfTilt=false, bool recordTJoint=false, bool recordClosedByHalfTilt=false);
Vector8LCoor tiltNeighborOctagon3(int level, int xIdx, int yIdx, int zIdx, const OctreeGridLayout3Ptr& layout);

bool checkLayoutValidity(const OctreeGridLayout3Ptr& layout);

// sdf, inside obstacle is negative
// splitMap should in finest resolution
template<typename SDF>
void refineSplitMapByObstacleSdf(CellCenteredFlagGrid3Ptr& splitMap, int l, const SDF& sdf, const Vector2f splitRange) {
	auto& res = splitMap->resolution();
    {
        int k = (1 << l);
        if(res.x() % k != 0 || res.y() % k != 0 || res.z() % k != 0) {
            throw std::runtime_error("resolution is not dividable by level");
        }
    }
    REAL spacing = splitMap->gridSpacing(); 
    Size3 cRes(res.x() >> l, res.y() >> l, res.z()>>l);
    REAL delta_h = spacing * (1 << l);
    int lCell = 1<<l;

    ThreadPool::parallelForTF(0, cRes.z(), [&](int tz){
        for(int ty = 0; ty < cRes.y(); ++ty) {
        for(int tx = 0; tx < cRes.x(); ++tx) {
            if(splitMap->get(Size3(tx<<l, ty<<l, tz<<l)) != l)
                continue;
            Vector3f pos = (Vector3f(tx, ty, tz) + Vector3f::Ones() * 0.5) * delta_h;
            REAL v = sdf->sample(pos);
            if(v >= splitRange[0] && v <= splitRange[1]) {   /*v in range*/
                for(int zi = (tz << l); zi < (tz << l) + lCell; ++zi) {
                for(int yi = (ty << l); yi < (ty << l) + lCell; ++yi) {
                for(int xi = (tx << l); xi < (tx << l) + lCell; ++xi) {
                    splitMap->set(Size3(xi, yi, zi), l-1);
                }}}
            }
        }}
	});
}


// NOTE: only for single level layout
template<typename SDF>
void refineSplitMapByObstacleSdf(BitMapPtr& splitMap, const SDF& sdf, const Vector2f splitRange, bool defV=false) {
	auto& res = splitMap->resolution();
    REAL delta_h = splitMap->gridSpacing(); 
    Size3 cRes(res.x() >> 1, res.y() >> 1, res.z()>>1);
    
    ThreadPool::taskParallelFor(0, cRes.z(), [&](int tz){
        for(int ty = 0; ty < cRes.y(); ++ty) {
        for(int tx = 0; tx < cRes.x(); ++tx) {
			bool flag = false;
			for(int zz = 0; zz < 2 && !flag; ++ zz) {
			for(int yy = 0; yy < 2 && !flag; ++ yy) {
			for(int xx = 0; xx < 2 && !flag; ++ xx) {
				Size3 idx((tx<<1)+xx, (ty<<1)+yy, (tz<<1)+zz);
            	Vector3f pos = (Vector3f(idx.x(), idx.y(), idx.z()) + Vector3f::Ones() * 0.5) * delta_h;
            	REAL v = sdf->sample(pos);
            	if(v >= splitRange[0] && v <= splitRange[1]) {   /*v in range*/
					flag = true;
				}
			}}}
			if(flag) {
				for(int zz = 0; zz < 2; ++ zz) {
				for(int yy = 0; yy < 2; ++ yy) {
				for(int xx = 0; xx < 2; ++ xx) {
					Size3 idx((tx<<1)+xx, (ty<<1)+yy, (tz<<1)+zz);
                    splitMap->set(idx, defV);
				}}}
			}
        }}
    }, 1);
}

inline Vector8LCoor halfTiltNeighborOctagon3(int level, int xIdx, int yIdx, int zIdx, int halfDirection) {
	Vector8LCoor ret = Vector8LCoor::Ones() * -1;
    Size3 idx(xIdx, yIdx, zIdx);
    for(int i = 0; i < 4; ++i) {
        int validSlot = HALF_TILT_SLOT_INDEX[halfDirection][i];
        Size3 tIdx = idx + TILT_NB_OCT_DXDYDZ[validSlot];
        ret(validSlot) = PACK_LCOOR3(level, tIdx.x(), tIdx.y(), tIdx.z());
    }
    return ret;
}

template<typename T>
inline T getTiltFromOctagonNeighbor3(int level, int xIdx, int yIdx, int zIdx, const std::shared_ptr<OctreeCellCenteredGrid3<T>>& cellPtr) {
	T v = cellPtr->get(level, xIdx, yIdx, zIdx);
	REAL count = 1;
	for(int m = -1; m <= 0; ++m) {  // find all adjacent cell of the vertex
		for (int n = -1; n <= 0; ++n) {
			for (int k = -1; k <= 0; ++k) {
				if (m == 0 && n == 0 && k == 0)
					continue;
				LCoor3 tp; tp.set(level, xIdx+m, yIdx+n, zIdx+k);
				if (cellPtr->isValid(tp.level, tp.xIdx, tp.yIdx, tp.zIdx) && cellPtr->getLayout()->get(tp.level, tp.xIdx, tp.yIdx, tp.zIdx) >= 0) {
					v += cellPtr->get(tp.level, tp.xIdx, tp.yIdx, tp.zIdx);
					count += 1;
				}
			}
		}
	}
	return v / count;
}


// Input: x,y,z \in [0, 1]
// Return: dir: 0-7
inline int getOctreeRegion3D(REAL x, REAL y, REAL z) {
	int dir = z >= 0.5 ? 0 : 4;
	dir += (x>=0.5) ? (y>=0.5?0:3) : (y>=0.5?1:2);
	return dir;
}


template<typename T>
void getNeighborCubeValue(int level, int xIdx, int yIdx, int zIdx, int dir, T* v, 
	const std::shared_ptr<OctreeVertexCenteredGrid3<T>>& vertexPtr,
	const std::shared_ptr<OctreeCellCenteredGrid3<T>>& cellPtr, 
	const OctreeTiltENodeGrid3Ptr& tiltE, bool T_joint=false) {
	int reDir = ((dir & 1) ? 8 : 6) - dir;
	T centerV = cellPtr->get(level, xIdx, yIdx, zIdx);
	v[reDir] = centerV;
	
	for(int i = 0; i < 8; ++i) {
		int marker = OCTREE_SAMPLE_MARKER[dir][i];
		if(marker == XYZ_MARKER::VOIDAXIS)
			continue;
		T tv = centerV;
		REAL count = 1;
		int markerXYZ[] = {(marker >> 8) & 0xf, (marker >> 4) & 0xf, (marker >> 0) & 0xf};

		if(T_joint && level > 0) {
			bool flag = false;
			for(int k = 0; k < 3; ++k) {
				if(markerXYZ[k] != 0 && markerXYZ[(k+1)%3] == 0 && markerXYZ[(k+2)%3] == 0) {	// face node
					Size3 idx(xIdx, yIdx, zIdx);
					idx[k] += (markerXYZ[k] == 0x1);
					idx *= 2;
					idx((k+1)%3) += 1; idx((k+2)%3) += 1;
					if(vertexPtr->getLayout()->get(level-1, idx.x(), idx.y(), idx.z()) >= 0 &&
						tiltE->get(level-1, idx.x(), idx.y(), idx.z()).is_half_tilt) {
						v[i] = vertexPtr->get(level-1, idx.x(), idx.y(), idx.z());
						flag = true;
						// DEBUG_ONLY(if(debug) {
						// 	spdlog::info("\t{}: face node get directly from: {}:{} {} {}, v:{}", i, level-1, idx.x(), idx.y(), idx.z(), v[i]);
						// });
					}
				}
				if(markerXYZ[k] == 0 && markerXYZ[(k+1)%3] != 0 && markerXYZ[(k+2)%3] != 0) {	// edge node
					Size3 idx((xIdx<<1)+(dir==0||dir==3||dir==4||dir==7), (yIdx<<1)+((dir%4)<2), (zIdx<<1)+(dir<4));
					Size3 tIdx = idx + OCT_NB_TILT_DXDYDZ[i];
					if(vertexPtr->getLayout()->get(level-1, tIdx.x(), tIdx.y(), tIdx.z()) >= 0
						&& !(vertexPtr->getLayout()->node(level-1, tIdx.x(), tIdx.y(), tIdx.z()).flag & LayerNodeFlagMask::IS_GHOST)) {
						auto tNode = tiltE->get(level-1, tIdx.x(), tIdx.y(), tIdx.z());
						if(tNode.is_closed_by_T_joint) {
							v[i] = vertexPtr->get(level-1, tIdx.x(), tIdx.y(), tIdx.z());
							flag = true;
							// DEBUG_ONLY(if(debug) {
							// 	spdlog::info("\t{}: edge node get directly from: {}:{} {} {}, v:{}", i, level-1, tIdx.x(), tIdx.y(), tIdx.z(), v[i]);
							// });
						}
					}
				}
			}
			if(flag) {
				continue;
			}
		}
		

		for(int x = -(markerXYZ[0]==0x2); x <= +(markerXYZ[0]==0x1); ++x) {
		for(int y = -(markerXYZ[1]==0x2); y <= +(markerXYZ[1]==0x1); ++y) {
		for(int z = -(markerXYZ[2]==0x2); z <= +(markerXYZ[2]==0x1); ++z) {
			if(x == 0 && y == 0 && z == 0)
				continue;
			LCOOR_T c = cellPtr->clamp(level, xIdx+x, yIdx+y, zIdx+z);
			// DEBUG_ONLY(if(debug) {
			// 	spdlog::info("\t{}: {}, {}, {}, clamped:{}, {}, {}, {}", i, xIdx+x, yIdx+y, zIdx+z, UNPACK_LEVEL3(c), UNPACK_X3(c), UNPACK_Y3(c), UNPACK_Z3(c));
			// });
			tv += cellPtr->get(c);
			count += 1;
		}}}
		v[i] = tv / count;
	}
}

}   // end of namespace Omni
#pragma once
#include "base_octree_omni_grid.h"
#include "global_benchmark.h"
#include "octree_cell_centered_grid3.h"
#include "octree_util.h"
#include "octree_vertex_centered_grid3.h"
#include "timer.h"
#include "util.h"

using std::shared_ptr;
using std::make_shared;

namespace Omni {

template <class T>
class OctreeOmniGrid3 : public BaseOctreeOmniGrid3 {
  public:
    OctreeOmniGrid3(const Size3& res, REAL spac, OctreeGridLayout3Ptr layout, OctreeTiltENodeGrid3Ptr tiltE, const T& init = zero<T>());
    ~OctreeOmniGrid3() {}

    OctreeOmniGrid3(OctreeOmniGrid3<T>& old);

    T sample(const Vector3f& pos) const { 
#ifdef SINGLE_LEVEL
        Vector3f rxyz = pos / gridSpacing() - Vector3f::Ones()*0.5;
        Vector3f f;
        Size3 idx;

        getBaryCentric(rxyz, resolution() - Size3::Ones(), idx, f);
        T v = cellPtr->sample(pos);
        if(getTiltEGridPtr()->get(0, idx.x()+1, idx.y()+1, idx.z()+1).is_closed)
            return v;
        T diff = vertexPtr->get(0, idx.x()+1, idx.y()+1, idx.z()+1) - 
            cellPtr->sample(getTiltEGridPtr()->position(0, idx.x()+1, idx.y()+1, idx.z()+1));
        return v + diff * 2 * std::min(f.minCoeff(), (Vector3f::Ones()-f).minCoeff());
#else
        return sample(pos.x(), pos.y(), pos.z()); 
#endif
    }
    T sample(REAL x, REAL y, REAL z) const;
    T sampleSub(REAL x, REAL y, REAL z) const;

    T get(LCOOR_T levelCoor) const {return get(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor));}
    T get(int level, int x, int y, int z) const {
        if(z & 1)
            return cellPtr->get(level, x, y, z >> 1);
        else
            return vertexPtr->get(level, x, y, z >> 1);
    }
    void set(LCOOR_T levelCoor, const T& v) {set(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor), v);}
    void set(int level, int x, int y, int z, const T& v) {
        if(z & 1)
            cellPtr->set(level, x, y, z >> 1, v);
        else
            vertexPtr->set(level, x, y, z >> 1, v);
    }

    T getOctagon(LCOOR_T levelCoor) const {return cellPtr->get(levelCoor);}
    T getOctagon(int level, int xIdx, int yIdx, int zIdx) const {return cellPtr->get(level, xIdx, yIdx, zIdx);}
    void setOctagon(LCOOR_T levelCoor, const T& v) {cellPtr->set(levelCoor, v);}
    void setOctagon(int level, int x, int y, int z, const T& v) {cellPtr->set(level, x, y, z, v);}
    T getTilt(LCOOR_T levelCoor) const {return vertexPtr->get(levelCoor);}
    T getTilt(int level, int xIdx, int yIdx, int zIdx) const {return vertexPtr->get(level, xIdx, yIdx, zIdx);}
    void setTilt(LCOOR_T levelCoor, const T& v) {vertexPtr->set(levelCoor, v);}
    void setTilt(int level, int x, int y, int z, const T& v) {vertexPtr->set(level, x, y, z, v);}

    Vector3f position(LCOOR_T levelCoor) const { return position(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor));}
    Vector3f position(int level, int x, int y, int z) const {
        if(z & 1)
            return cellPtr->position(level, x, y, z >> 1);
        else
            return vertexPtr->position(level, x, y, z >> 1);
    }

    Vector3f positionOctagon(LCOOR_T levelCoor) const {return cellPtr->position(levelCoor);}
    Vector3f positionOctagon(int level, int x, int y, int z) const {return cellPtr->position(level, x, y, z);}
    Vector3f positionTilt(LCOOR_T levelCoor) const {return vertexPtr->position(levelCoor);}
    Vector3f positionTilt(int level, int x, int y, int z) const {return vertexPtr->position(level, x, y, z);}

    // OctreeOmniGridIterator3 getIterator() {return OctreeOmniGridIterator3(layout);}
    bool isValid(int level, int xIdx, int yIdx, int zIdx) const {
        if(zIdx&1)
            return cellPtr->isValid(level, xIdx, yIdx, zIdx >> 1);
        else 
            return vertexPtr->isValid(level, xIdx, yIdx, zIdx >> 1);
    }

    shared_ptr<OctreeVertexCenteredGrid3<T>>& getTiltGrid() { return vertexPtr; }
    shared_ptr<OctreeCellCenteredGrid3<T>>& getOctagonGrid() { return cellPtr; }

    void updateGhost() override;
    void updateCellGhost();
    void interpolateClosedJoint();

    void setAllGhost(const T& v) {
        vertexPtr->setAllGhost(v);
        cellPtr->setAllGhost(v);
    }

    void setAll(const T& v) {
        vertexPtr->setAll(v);
        cellPtr->setAll(v);
    }

    void scale(Real s) {
        vertexPtr->scale(s);
        cellPtr->scale(s);
    }

private:
    shared_ptr<OctreeVertexCenteredGrid3<T>> vertexPtr;
    shared_ptr<OctreeCellCenteredGrid3<T>> cellPtr;
};

template<class T>
OctreeOmniGrid3<T>::OctreeOmniGrid3(const Size3& res, REAL spac, OctreeGridLayout3Ptr layout, OctreeTiltENodeGrid3Ptr tiltE, const T& init)
    : BaseOctreeOmniGrid3(res, spac, layout, tiltE)
{
    vertexPtr = make_shared<OctreeVertexCenteredGrid3<T>>(res, spac, layout, init);
    cellPtr = make_shared<OctreeCellCenteredGrid3<T>>(res, spac, layout, init);
}

template<class T>
OctreeOmniGrid3<T>::OctreeOmniGrid3(OctreeOmniGrid3<T>& old)
	: BaseOctreeOmniGrid3(old.resolution(), old.gridSpacing(), old.getLayout(), old.getTiltEGridPtr())
{
	tiltE = make_shared<OctreeTiltENodeGrid3>(*old.getTiltEGridPtr());
	vertexPtr = make_shared<OctreeVertexCenteredGrid3<T>>(*old.getTiltGrid());
	cellPtr = make_shared<OctreeCellCenteredGrid3<T>>(*old.getOctagonGrid());
}

// override for flag omnigrid. mainly used for index record before project.
template<>
void OctreeOmniGrid3<int>::updateGhost();

template<class T>
void OctreeOmniGrid3<T>::updateGhost() {
    BENCHMARK_SCOPED_TIMER_SECTION t("omni grid,update ghost cell");
    cellPtr->updateGhost(); // half correct, need updateCellGhost to correct it.
    update_ghost_exclude_half_tilt(T, vertexPtr);
    updateCellGhost();
    interpolateClosedJoint();
}

// always construct a cube to do trilinear interpolation
template <class T>
T OctreeOmniGrid3<T>::sample(Real x_pos, Real y_pos, Real z_pos) const {
#ifdef SINGLE_LEVEL
    return sample(Vector3f(x_pos, y_pos, z_pos));
#elif defined(OCTREE_CORRECTION_BASED)  //correction based
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    int level = UNPACK_LEVEL3(levelCoor);
    Real delta_h = layout->levelGridSpacing(level);
    const Size3& cRes = layout->levelResolution(level);
    Real fx, fy, fz;
    int xIdx, yIdx, zIdx;
    Vector3f rxyz = Vector3f(x_pos, y_pos, z_pos) / delta_h - Vector3f::Ones() * 0.5;

    getBarycentric(rxyz.x(), 0, cRes.x()-2, &xIdx, &fx);
    getBarycentric(rxyz.y(), 0, cRes.y()-2, &yIdx, &fy);
    getBarycentric(rxyz.z(), 0, cRes.z()-2, &zIdx, &fz);
    
    int xIdx1 = std::min(xIdx+1, cRes.x()-2);
    int yIdx1 = std::min(yIdx+1, cRes.y()-2);
    int zIdx1 = std::min(zIdx+1, cRes.z()-2);

    T ret = trilerp(getOctagon(level, xIdx, yIdx, zIdx), getOctagon(level, xIdx1, yIdx, zIdx), getOctagon(level, xIdx, yIdx1, zIdx),
                    getOctagon(level, xIdx1, yIdx1, zIdx), getOctagon(level, xIdx, yIdx, zIdx1), getOctagon(level, xIdx1, yIdx, zIdx1),
                    getOctagon(level, xIdx, yIdx1, zIdx1), getOctagon(level, xIdx1, yIdx1, zIdx1), fx, fy, fz);
    const auto& tNode = getTiltEGridPtr()->get(level, xIdx+1, yIdx+1, zIdx+1);

    if((getLayout()->node(level, xIdx, yIdx, zIdx).flag & LayerNodeFlagMask::IS_GHOST) || 
        (level == 0 && (getLayout()->get(level+1, xIdx>>1, yIdx>>1, zIdx>>1) >= 0)) ||
        (level == 1 && (getLayout()->get(level-1, (xIdx<<1)+1, (yIdx<<1)+1, (zIdx<<1)+1) >= 0)) ||
        (level == 1 && (getLayout()->get(level-1, (xIdx<<1), (yIdx<<1), (zIdx<<1)) >= 0)) ) {
        return sampleSub(x_pos, y_pos, z_pos);
    }

    if(tNode.is_closed && !tNode.is_half_tilt) {
        return ret;
    } else {
        T diff = getTilt(level, xIdx+1, yIdx+1, zIdx+1)
            - (1/8.0) * (getOctagon(level, xIdx, yIdx, zIdx) + getOctagon(level, xIdx1, yIdx, zIdx) + getOctagon(level, xIdx, yIdx1, zIdx)
            + getOctagon(level, xIdx1, yIdx1, zIdx) + getOctagon(level, xIdx, yIdx, zIdx1) + getOctagon(level, xIdx1, yIdx, zIdx1)
            + getOctagon(level, xIdx, yIdx1, zIdx1) + getOctagon(level, xIdx1, yIdx1, zIdx1));
        return ret + diff * 2.0 * std::min(Vector3f(fx, fy, fz).minCoeff(), (Vector3f::Ones() - Vector3f(fx, fy, fz)).minCoeff());
    }
#else
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    int level = UNPACK_LEVEL3(levelCoor);
    Real delta_h = layout->levelGridSpacing(level);
    const Size3& cRes = layout->levelResolution(level);
    Real rx = x_pos / delta_h, ry = y_pos / delta_h, rz = z_pos / delta_h;
    Real fx, fy, fz;
    int xIdx, yIdx, zIdx;

    getBarycentric(rx, 0, cRes.x()-1, &xIdx, &fx);
    getBarycentric(ry, 0, cRes.y()-1, &yIdx, &fy);
    getBarycentric(rz, 0, cRes.z()-1, &zIdx, &fz);


    int region3D = getOctreeRegion3D(fx, fy, fz);
    Size3 idx(xIdx, yIdx, zIdx);
    Size3 tiltIdx = idx + OCT_NB_TILT_DXDYDZ[region3D];
    const auto& tNode = getTiltEGridPtr()->get(level, tiltIdx.x(), tiltIdx.y(), tiltIdx.z());
    // DEBUG_ONLY(if(debug) {
    //     spdlog::info("level:{}, idx:{}, region:{}, fx,fy,fz:{}, {}, {}, tiltIdx:{}", level, idx.transpose(), region3D, fx, fy, fz, tiltIdx.transpose());
    //     spdlog::info("T_joint:{}, closed:{}, half_tilt:{}, is_closed_by_T_joint:{}, dir:{}", tNode.is_T_joint, tNode.is_closed, tNode.is_half_tilt, tNode.is_closed_by_T_joint, int(tNode.half_direction));
    // });
    if(!tNode.is_half_tilt && !tNode.is_closed_by_T_joint) {
        if(tNode.is_closed) {
            return cellPtr->sampleGivenRelativePosition(level, rx, ry, rz);
        }
    }

    T v[8];
    v[region3D] = vertexPtr->get(level, tiltIdx.x(), tiltIdx.y(), tiltIdx.z());
    getNeighborCubeValue(level, xIdx, yIdx, zIdx, region3D, v, vertexPtr, cellPtr, getTiltEGridPtr(), tNode.is_T_joint);

    switch(region3D) {
        case 0: fx-=0.5, fy-=0.5, fz-=0.5; break;
        case 1: fy-=0.5, fz-=0.5; break;
        case 2: fz-=0.5; break;
        case 3: fx-=0.5, fz-=0.5; break;
        case 4: fx-=0.5, fy-=0.5; break;
        case 5: fy-=0.5; break;
        case 6: break;
        case 7: fx-=0.5; break;
    }
    fx *= 2; fy *= 2; fz *= 2;
    ASSERT(fx<=1 && fy <= 1 && fz <= 1, "error fx,fy,fz");
    // DEBUG_ONLY(if(debug) {
    //     spdlog::info("idx:{}, region:{}, fx,fy,fz:{}, {}, {}, tiltIdx:{}", idx.transpose(), region3D, fx, fy, fz, tiltIdx.transpose());
    //     spdlog::info("{}, {}, {}, {}, {}, {}, {}, {}", v[0], v[1],v[2],v[3],v[4],v[5],v[6],v[7]);
    // });

    return trilerp(v[6], v[7], v[5], v[4], v[2], v[3], v[1], v[0], fx, fy, fz);
#endif
}

// always construct a cube to do trilinear interpolation
template <class T>
T OctreeOmniGrid3<T>::sampleSub(Real x_pos, Real y_pos, Real z_pos) const {
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    int level = UNPACK_LEVEL3(levelCoor);
    Real delta_h = layout->levelGridSpacing(level);
    const Size3& cRes = layout->levelResolution(level);
    Real rx = x_pos / delta_h, ry = y_pos / delta_h, rz = z_pos / delta_h;
    Real fx, fy, fz;
    int xIdx, yIdx, zIdx;

    getBarycentric(rx, 0, cRes.x()-1, &xIdx, &fx);
    getBarycentric(ry, 0, cRes.y()-1, &yIdx, &fy);
    getBarycentric(rz, 0, cRes.z()-1, &zIdx, &fz);


    int region3D = getOctreeRegion3D(fx, fy, fz);
    Size3 idx(xIdx, yIdx, zIdx);
    Size3 tiltIdx = idx + OCT_NB_TILT_DXDYDZ[region3D];
    const auto& tNode = getTiltEGridPtr()->get(level, tiltIdx.x(), tiltIdx.y(), tiltIdx.z());

    if(!tNode.is_half_tilt && !tNode.is_closed_by_T_joint) {
        if(tNode.is_closed) {
            return cellPtr->sampleGivenRelativePosition(level, rx, ry, rz);
        }
    }

    T v[8];
    v[region3D] = vertexPtr->get(level, tiltIdx.x(), tiltIdx.y(), tiltIdx.z());
    getNeighborCubeValue(level, xIdx, yIdx, zIdx, region3D, v, vertexPtr, cellPtr, getTiltEGridPtr(), tNode.is_T_joint);

    switch(region3D) {
        case 0: fx-=0.5, fy-=0.5, fz-=0.5; break;
        case 1: fy-=0.5, fz-=0.5; break;
        case 2: fz-=0.5; break;
        case 3: fx-=0.5, fz-=0.5; break;
        case 4: fx-=0.5, fy-=0.5; break;
        case 5: fy-=0.5; break;
        case 6: break;
        case 7: fx-=0.5; break;
    }
    fx *= 2; fy *= 2; fz *= 2;
    ASSERT(fx<=1 && fy <= 1 && fz <= 1, "error fx,fy,fz");

    return trilerp(v[6], v[7], v[5], v[4], v[2], v[3], v[1], v[0], fx, fy, fz);
}


template<class T>
void OctreeOmniGrid3<T>::updateCellGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        const auto& node = iter.node();
        if(!cellPtr->isValid(iter.level, node.xIdx, node.yIdx, node.zIdx))
            return;
        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            // spdlog::info("{}, {}, {}", l, node.xIdx, node.yIdx);
            LCoor3 dp; dp.set(iter.level, node.xIdx+(node.xIdx&1), node.yIdx+(node.yIdx&1), node.zIdx+(node.zIdx&1));
            ASSERT(getTiltEGridPtr()->get(dp.levelCoor()).is_T_joint, "OmniOctree3::UpdateGhost:: should have T_joint across level");
            cellPtr->set(iter.level, node.xIdx, node.yIdx, node.zIdx, 
                0.5 * (vertexPtr->get(dp.levelCoor()) + cellPtr->get(iter.level+1, node.xIdx>>1, node.yIdx>>1, node.zIdx>>1)));
        }
    });
}


template<class T>
void OctreeOmniGrid3<T>::interpolateClosedJoint() {
    ThreadPool::parallelForTF(0, (int)layout->getTJointArr().size(), [&](int idx){
        LCOOR_T lCoor = layout->getTJointArr()[idx];
        const auto& tiltENode = getTiltEGridPtr()->get(lCoor);
        int level = UNPACK_LEVEL3(lCoor);
        LCoor3 dp; dp.set(level, UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor));

        if(tiltENode.e < layout->levelGridSpacing(level)) {    // coarse
            dp.set(level-1, dp.xIdx<<1, dp.yIdx<<1, dp.zIdx<<1);
        }
        
        // loop over tilt's vertex
        for(int i = 0; i < 6; ++i) {
            LCoor3 tDp(dp.level, dp.xIdx+OCT_NB_DXDYDZ[i].x(), dp.yIdx+OCT_NB_DXDYDZ[i].y(), dp.zIdx+OCT_NB_DXDYDZ[i].z());
            if(layout->get(tDp.level, tDp.xIdx, tDp.yIdx, tDp.zIdx) >= 0) {    // vertex may be invalid
                vertexPtr->set(tDp.levelCoor(), getTiltFromOctagonNeighbor3(tDp.level, tDp.xIdx, tDp.yIdx, tDp.zIdx, cellPtr));
            }
        }
    });

    // interplocate tilt closed by half-tilt
    ThreadPool::parallelForTF(0, (int)layout->getclosedByHalfTiltArr().size(), [&](int idx){
        LCOOR_T lCoor = layout->getclosedByHalfTiltArr()[idx];
        const auto& tiltENode = getTiltEGridPtr()->get(lCoor);
        if(tiltENode.is_closed && tiltENode.half_direction == -1) {
            vertexPtr->set(lCoor, getTiltFromOctagonNeighbor3(UNPACK_LEVEL3(lCoor), UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor), cellPtr));
        } else {
            ASSERT(false, "what???");
        }
    });
}


typedef OctreeOmniGrid3<int> OctreeOmniFlagGrid3;
typedef shared_ptr<OctreeOmniFlagGrid3> OctreeOmniFlagGrid3Ptr;
typedef OctreeOmniGrid3<char> OctreeOmniCharGrid3;
typedef shared_ptr<OctreeOmniCharGrid3> OctreeOmniCharGrid3Ptr;
typedef OctreeOmniGrid3<REAL> OctreeOmniScalarGrid3;
typedef shared_ptr<OctreeOmniScalarGrid3> OctreeOmniScalarGrid3Ptr;
typedef OctreeOmniGrid3<Vector3f> OctreeOmniVector3Grid3;
typedef shared_ptr<OctreeOmniVector3Grid3> OctreeOmniVector3Grid3Ptr;

} // end of namespace Omni
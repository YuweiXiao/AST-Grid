#pragma once
#include "base_octree_grid.h"
#include "tilt_e_node.h"

namespace Omni {

template<class T>
class OctreeVertexCenteredGrid3 : public BaseOctreeGrid3 {
public:
    OctreeVertexCenteredGrid3(const Size3& res, REAL spacing, const OctreeGridLayout3Ptr& layout, const T& init=T());
    virtual ~OctreeVertexCenteredGrid3() {}

	OctreeVertexCenteredGrid3(OctreeVertexCenteredGrid3<T> & old);
    
    // sample using absolute positon
    T sample(const Vector3f& pos) const { return sample(pos.x(), pos.y(), pos.z()); }
    T sample(REAL x_pos, REAL y_pos, REAL z_pos) const;

    Vector3f position(LCOOR_T levelCoor) const{return layout->position(levelCoor);}
    Vector3f position(int level, int xIdx, int yIdx, int zIdx) const {return layout->position(level, xIdx, yIdx, zIdx);} 

    const T& get(LCOOR_T levelCoor) const {return getByIndex(UNPACK_LEVEL3(levelCoor), layout->get(levelCoor)); }
    const T& get(int level, int x, int y, int z) const {return getByIndex(level, layout->get(level, x, y, z)); }
    const T& get(int level, const Size3& idx) const {return getByIndex(level, layout->get(level, idx.x(), idx.y(), idx.z())); }
    T& get(LCOOR_T levelCoor) {return getByIndex(UNPACK_LEVEL3(levelCoor), layout->get(levelCoor)); }
    T& get(int level, int x, int y, int z) {return getByIndex(level, layout->get(level, x, y, z)); }
    T& get(int level, const Size3& idx) {return getByIndex(level, layout->get(level, idx.x(), idx.y(), idx.z())); }
    void set(LCOOR_T levelCoor, const T& v) {set(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor), v);}
    void set(int level, int x, int y, int z, const T& v) {
        auto idx = layout->get(level, x, y, z);
        ASSERT(Valid(Size2(level, idx), Size2(layout->level(), data[level].size())), "OctreeVertexCentered3::set::out of range"); 
        data[level][idx] = v;
    }
    void setAll(const T& v) {
        for(int l = 0; l < layout->level(); ++l) {
            std::fill(data[l].begin(), data[l].end(), v);
        }
    }
    void fill(std::shared_ptr<OctreeVertexCenteredGrid3<T> >& old) {
        data.resize(old->getData().size());
        copy(old->getData().begin(), old->getData().end(), data.begin());
    }

    bool isValid(int level, int x, int y, int z) const { return Valid(Size3(x, y, z), layout->levelResolution(level));}
    virtual void updateGhost() override;
    void setAllGhost(const T& v) {
        ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
            set(iter.levelCoor(), v);
        });
    }

    std::vector<std::vector<T>>& getData() { return data; }
    void scale(Real s) {
        for (auto& d : data) {
            std::transform(d.begin(), d.end(), d.begin(), std::bind(std::multiplies<T>(), std::placeholders::_1, s));
        }
    }

private:

    const T& getByIndex(int level, int index) const {
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), data[level].size())), "OctreeVertexCentered3::set::out of range:"+std::to_string(index)); 
        return data[level][index];
    }
    T& getByIndex(int level, int index) {
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), data[level].size())), "OctreeVertexCentered3::set::out of range:"+std::to_string(index)); 
        return data[level][index];
    }

    std::vector<std::vector<T>> data;   // store each layers' data
};

template<class T>
OctreeVertexCenteredGrid3<T>::OctreeVertexCenteredGrid3(const Size3& res, REAL spacing, const OctreeGridLayout3Ptr& _layout, const T& init) 
    : BaseOctreeGrid3(res, spacing, _layout)
{
    data.resize(layout->level(), std::vector<T>());
    for(int i = 0; i < layout->level(); ++i) {
        data[i].resize(layout->levelSize(i), init);
    }
}

template<class T>
OctreeVertexCenteredGrid3<T>::OctreeVertexCenteredGrid3(OctreeVertexCenteredGrid3<T>& old)
	:BaseOctreeGrid3(old.resolution(), old.gridSpacing(), old.getLayout())
{
	data.resize(layout->level(), std::vector<T>());
	copy(old.getData().begin(), old.getData().end(), data.begin());
}

template<class T>
T OctreeVertexCenteredGrid3<T>::sample(REAL x_pos, REAL y_pos, REAL z_pos) const {
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    int level = UNPACK_LEVEL3(levelCoor);
    REAL delta_h = layout->levelGridSpacing(level);
    Size3 cRes = layout->levelResolution(level);
    FLOAT rx = x_pos / delta_h;
    FLOAT ry = y_pos / delta_h;
    FLOAT rz = z_pos / delta_h;
    int xIdx, yIdx, zIdx;
    REAL fx, fy, fz;

    getBarycentric(rx, 0, cRes.x()-1, &xIdx, &fx);
    getBarycentric(ry, 0, cRes.y()-1, &yIdx, &fy);
    getBarycentric(rz, 0, cRes.z()-1, &zIdx, &fz);

    int xIdx1 = std::min(xIdx+1, cRes.x()-1);
    int yIdx1 = std::min(yIdx+1, cRes.y()-1);
    int zIdx1 = std::min(zIdx+1, cRes.z()-1);

    return trilerp(get(level, xIdx, yIdx, zIdx),
                    get(level, xIdx1, yIdx, zIdx),
                    get(level, xIdx, yIdx1, zIdx),
                    get(level, xIdx1, yIdx1, zIdx),
                    get(level, xIdx, yIdx, zIdx1),
                    get(level, xIdx1, yIdx, zIdx1),
                    get(level, xIdx, yIdx1, zIdx1),
                    get(level, xIdx1, yIdx1, zIdx1), fx, fy, fz);
}


template<>void OctreeVertexCenteredGrid3<TiltENode>::updateGhost();
template<>void OctreeVertexCenteredGrid3<TiltENodeS>::updateGhost();
template<>void OctreeVertexCenteredGrid3<int>::updateGhost();
template<>void OctreeVertexCenteredGrid3<Vector3f>::updateGhost();


template<class T>
void OctreeVertexCenteredGrid3<T>::updateGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        int l = iter.level;
        const auto& node = iter.node();
        // spdlog::info("ghost cell:{}, node.zIdx, idx:{}, {}, flag: {}", l, node.xIdx, node.yIdx, node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE);
        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            T v = zero<T>();
            int count = 0;
            for(int dx = 0; dx <= (node.xIdx&1); ++dx) {
                for(int dy = 0; dy <= (node.yIdx&1); ++dy) {
                    for(int dz = 0; dz <= (node.zIdx&1); ++dz) {
                        LCOOR_T c = layout->getUntil(l+1, (node.xIdx>>1)+dx, (node.yIdx>>1)+dy, (node.zIdx>>1)+dz, true);
                        if(c>=0) {
                            v += get(c);
                            count += 1;
                        }
                    }
                }
            }
            set(l, node.xIdx, node.yIdx, node.zIdx, v * (1.0/(count+EPSILON)));
        } else {
            set(l, node.xIdx, node.yIdx, node.zIdx, get(layout->getUntil(l-1, node.xIdx<<1, node.yIdx<<1, node.zIdx<<1, true)));
        }
    });
}

typedef OctreeVertexCenteredGrid3<REAL> OctreeVertexCenteredScalarGrid3;
typedef std::shared_ptr<OctreeVertexCenteredScalarGrid3> OctreeVertexCenteredScalarGrid3Ptr;
typedef OctreeVertexCenteredGrid3<TiltENode> OctreeTiltENodeGrid3;
typedef std::shared_ptr<OctreeTiltENodeGrid3> OctreeTiltENodeGrid3Ptr;
typedef OctreeVertexCenteredGrid3<TiltENodeS> OctreeTiltENodeSGrid3;
typedef std::shared_ptr<OctreeTiltENodeSGrid3> OctreeTiltENodeSGrid3Ptr;
typedef OctreeVertexCenteredGrid3<Vector8f> OctreeVertexCenteredVector8fGrid3;
typedef std::shared_ptr<OctreeVertexCenteredVector8fGrid3> OctreeVertexCenteredVector8fGrid3Ptr;


} // end of namespace Omni
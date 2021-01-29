#pragma once
#include "base_octree_grid.h"
#include "thread_pool.h"

namespace Omni {

template<class T>
class OctreeCellCenteredGrid3 : public BaseOctreeGrid3 {
public:
    OctreeCellCenteredGrid3(const Size3& res, REAL spacing, OctreeGridLayout3Ptr layout, const T& init=zero<T>());
    ~OctreeCellCenteredGrid3() {}

	OctreeCellCenteredGrid3(OctreeCellCenteredGrid3<T> &old);

    // sample using absolute positon
    T sample(const Vector3f& pos) const { return sample(pos.x(), pos.y(), pos.z()); }
    T sample(REAL x_pos, REAL y_pos, REAL z_pos) const;
    T sampleSingleLayer(const Vector3f& pos, T def=zero<T>()) const { return sampleSingleLayer(pos.x(), pos.y(), pos.z(), def); }
    T sampleSingleLayer(REAL x_pos, REAL y_pos, REAL z_pos, T def=zero<T>()) const;
    T sampleGivenRelativePosition(int level, REAL rx, REAL ry, REAL rz, T def=zero<T>()) const;

    Vector3f position(LCOOR_T levelCoor) const { return layout->position(levelCoor) + Vector3f::Ones()*0.5*layout->levelGridSpacing(UNPACK_LEVEL3(levelCoor)); }
    Vector3f position(int level, int xIdx, int yIdx, int zIdx) const { return layout->position(level, xIdx, yIdx, zIdx) + Vector3f::Ones()*0.5*layout->levelGridSpacing(level); }

    T get(LCOOR_T levelCoor, T def=zero<T>()) const {return getByIndex(UNPACK_LEVEL3(levelCoor), layout->get(levelCoor), def); }
    T get(int level, int x, int y, int z, T def=zero<T>()) const {return getByIndex(level, layout->get(level, x, y, z), def); }
    T get(int level, const Size3& idx) const {return getByIndex(level, layout->get(level, idx.x(), idx.y(), idx.z()));}
    void set(LCOOR_T levelCoor, const T& v) {set(UNPACK_LEVEL3(levelCoor), UNPACK_X3(levelCoor), UNPACK_Y3(levelCoor), UNPACK_Z3(levelCoor), v);}
    void set(int level, int x, int y, int z, const T& v) {
        auto idx = layout->get(level, x, y, z);
        ASSERT(Valid(Size2(level, idx), Size2(layout->level(), data[level].size())), "OctreeCellCentered::set::out of range: " + std::to_string(idx)); 
        data[level][idx] = v;
    }

    bool isValid(int level, int xIdx, int yIdx, int zIdx) const {return Valid(Size3(xIdx, yIdx, zIdx), layout->levelResolution(level) - Size3::Ones());}
    LCOOR_T clamp(int level, int xIdx, int yIdx, int zIdx) const {
        auto res = layout->levelResolution(level);
        xIdx = std::min(std::max(0, xIdx), res.x()-2);
        yIdx = std::min(std::max(0, yIdx), res.y()-2);
        zIdx = std::min(std::max(0, zIdx), res.z()-2);
        return PACK_LCOOR3(level, xIdx, yIdx, zIdx);
    }
    void updateGhost() override;
    void setAllGhost(const T& v) {  // TODO could speed up
        ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
            set(iter.levelCoor(), v);
        });
    }
    void setAll(const T& v) {
        for(int l = 0; l < layout->level(); ++l) {
            std::fill(data[l].begin(), data[l].end(), v);
        }
    }

    void fill(std::shared_ptr<OctreeCellCenteredGrid3<T> >& old) {
        data.resize(old->getData().size());
        copy(old->getData().begin(), old->getData().end(), data.begin());
    }
	std::vector<std::vector<T>>& getData() { return data; }
    void scale(Real s) {
        for (auto& d : data) {
            std::transform(d.begin(), d.end(), d.begin(), std::bind(std::multiplies<T>(), std::placeholders::_1, s));
        }
    }
    bool allowNonValidLCoor = false;
private:

    T getByIndex(int level, int index, T def=zero<T>()) const {
        if(allowNonValidLCoor && index == -1)
            return def;
        ASSERT(Valid(Size2(level, index), Size2(layout->level(), data[level].size())), "OctreeCellCentered::get::out of range"); 
        return data[level][index];
    }

    std::vector<std::vector<T>> data;   // store each layers' data
};

template<class T>
OctreeCellCenteredGrid3<T>::OctreeCellCenteredGrid3(const Size3& res, REAL spacing, OctreeGridLayout3Ptr _layout, const T& init) 
    : BaseOctreeGrid3(res, spacing, _layout)
{
    data.resize(layout->level(), std::vector<T>());
    for(int i = 0; i < layout->level(); ++i) {
        data[i].resize(layout->levelSize(i), init);
    }
}

template<class T>
OctreeCellCenteredGrid3<T>::OctreeCellCenteredGrid3(OctreeCellCenteredGrid3<T>& old)
	:BaseOctreeGrid3(old.resolution(), old.gridSpacing(), old.getLayout())
{
	data.resize(layout->level(), std::vector<T>());
	copy(old.getData().begin(), old.getData().end(), data.begin());
}

template<class T>
T OctreeCellCenteredGrid3<T>::sample(REAL x_pos, REAL y_pos, REAL z_pos) const {
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    ASSERT(levelCoor >=0, "[OctreeCellCenteredGrid] level coordinate fail");
    int level = UNPACK_LEVEL3(levelCoor);
    REAL delta_h = layout->levelGridSpacing(level);
    Size3 cRes = layout->levelResolution(level);
    REAL rx = x_pos / delta_h, ry = y_pos / delta_h, rz = z_pos / delta_h;
    int xIdx, yIdx, zIdx;
    REAL fx, fy, fz;

    getBarycentric(rx-0.5, 0, cRes.x()-2, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, cRes.y()-2, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, cRes.z()-2, &zIdx, &fz);

    int xIdx1 = std::min(xIdx+1, cRes.x()-2);
    int yIdx1 = std::min(yIdx+1, cRes.y()-2);
    int zIdx1 = std::min(zIdx+1, cRes.z()-2);

    return trilerp(get(level, xIdx, yIdx, zIdx),
                    get(level, xIdx1, yIdx, zIdx),
                    get(level, xIdx, yIdx1, zIdx),
                    get(level, xIdx1, yIdx1, zIdx),
                    get(level, xIdx, yIdx, zIdx1),
                    get(level, xIdx1, yIdx, zIdx1),
                    get(level, xIdx, yIdx1, zIdx1),
                    get(level, xIdx1, yIdx1, zIdx1), fx, fy, fz);
}

template<class T>
T OctreeCellCenteredGrid3<T>::sampleSingleLayer(REAL x_pos, REAL y_pos, REAL z_pos, T def) const {
    LCOOR_T levelCoor = layout->levelCoordinate(x_pos, y_pos, z_pos);
    if(allowNonValidLCoor && levelCoor == -1)
        return def;
    ASSERT(levelCoor >=0, "[OctreeCellCenteredGrid] level coordinate fail");
    int level = UNPACK_LEVEL3(levelCoor);
    REAL delta_h = layout->levelGridSpacing(level);
    Size3 cRes = layout->levelResolution(level);
    REAL rx = x_pos / delta_h, ry = y_pos / delta_h, rz = z_pos / delta_h;
    int xIdx, yIdx, zIdx;
    REAL fx, fy, fz;

    getBarycentric(rx-0.5, 0, cRes.x()-2, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, cRes.y()-2, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, cRes.z()-2, &zIdx, &fz);

    int xIdx1 = std::min(xIdx+1, cRes.x()-2);
    int yIdx1 = std::min(yIdx+1, cRes.y()-2);
    int zIdx1 = std::min(zIdx+1, cRes.z()-2);

    T c = get(level, xIdx, yIdx, zIdx);
    return trilerp(get(level, xIdx, yIdx, zIdx),
                    get(level, xIdx1, yIdx, zIdx, c),
                    get(level, xIdx, yIdx1, zIdx, c),
                    get(level, xIdx1, yIdx1, zIdx, c),
                    get(level, xIdx, yIdx, zIdx1, c),
                    get(level, xIdx1, yIdx, zIdx1, c),
                    get(level, xIdx, yIdx1, zIdx1, c),
                    get(level, xIdx1, yIdx1, zIdx1, c), fx, fy, fz);
}

template<class T>
T OctreeCellCenteredGrid3<T>::sampleGivenRelativePosition(int level, REAL rx, REAL ry, REAL rz, T def) const {
    Size3 cRes = layout->levelResolution(level);
    int xIdx, yIdx, zIdx;
    REAL fx, fy, fz;

    getBarycentric(rx-0.5, 0, cRes.x()-2, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, cRes.y()-2, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, cRes.z()-2, &zIdx, &fz);

    int xIdx1 = std::min(xIdx+1, cRes.x()-2);
    int yIdx1 = std::min(yIdx+1, cRes.y()-2);
    int zIdx1 = std::min(zIdx+1, cRes.z()-2);

    return trilerp(get(level, xIdx, yIdx, zIdx, def),
                    get(level, xIdx1, yIdx, zIdx, def),
                    get(level, xIdx, yIdx1, zIdx, def),
                    get(level, xIdx1, yIdx1, zIdx, def),
                    get(level, xIdx, yIdx, zIdx1, def),
                    get(level, xIdx1, yIdx, zIdx1, def),
                    get(level, xIdx, yIdx1, zIdx1, def),
                    get(level, xIdx1, yIdx1, zIdx1, def), fx, fy, fz);
}

template<> void OctreeCellCenteredGrid3<int>::updateGhost();

template<class T>
void OctreeCellCenteredGrid3<T>::updateGhost() {
    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter){
        const LayerNode3& node = iter.node();

        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            // set(iter.level, node.xIdx, node.yIdx, node.zIdx, get(iter.level+1, node.xIdx>>1, node.yIdx>>1, node.zIdx>>1));
        } else {
            T v = zero<T>();
            int x = node.xIdx << 1, y = node.yIdx << 1, z = node.zIdx << 1;
            for(int i = 0; i <= 1; ++i) 
                for(int j = 0; j <= 1; ++j) 
                    for(int k = 0; k <= 1; ++k) 
                        v += get(iter.level-1, x+i, y+j, z+k);
            set(iter.level, node.xIdx, node.yIdx, node.zIdx, v * (1/8.0));
        }
    });

    ThreadPool::parallelIterateGrid(getGhostIterator(), [&](const OctreeGhostCellIterator3& iter){
        const LayerNode3& node = iter.node();

        if(node.flag & LayerNodeFlagMask::GHOST_FROM_COARSE) {
            set(iter.level, node.xIdx, node.yIdx, node.zIdx, sample(position(iter.levelCoor())));
        }
    });
}

typedef OctreeCellCenteredGrid3<REAL> OctreeCellCenteredScalarGrid3;
typedef std::shared_ptr<OctreeCellCenteredScalarGrid3> OctreeCellCenteredScalarGrid3Ptr;
typedef OctreeCellCenteredGrid3<int> OctreeCellCenteredFlagGrid3;
typedef std::shared_ptr<OctreeCellCenteredFlagGrid3> OctreeCellCenteredFlagGrid3Ptr;
typedef OctreeCellCenteredGrid3<char> OctreeCellCenteredCharGrid3;
typedef std::shared_ptr<OctreeCellCenteredCharGrid3> OctreeCellCenteredCharGrid3Ptr;

} // end of namespace Omni
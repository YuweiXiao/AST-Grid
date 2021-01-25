#pragma once
#include <sys/mman.h>
#include "grid.h"
#include "thread_pool.h"

namespace Omni {

template<class T>
class BitMapCellCenteredGrid3 : public Grid3 {
public:
    BitMapCellCenteredGrid3(const Size3& res, REAL spacing, BitMapPtr bitmap);
    BitMapCellCenteredGrid3(const Size3& res, REAL spacing, BitMapPtr bitmap, const T& init);
    BitMapCellCenteredGrid3(const Size3& res, REAL spacing, BitMapPtr bitmap, const BitMapCellCenteredGrid3<T>& old);
    ~BitMapCellCenteredGrid3() {int ret = munmap(data, totalSize); if(ret != 0) { spdlog::info("[BitMapCellCenteredGrid3]::unmap error: {}", ret);}}

    // sample using absolute positon
    T sample(const Vector3f& pos, const T& def=zero<T>()) const { return sample(pos.x(), pos.y(), pos.z(), def); }
    T sample(REAL x_pos, REAL y_pos, REAL z_pos, const T& def=zero<T>()) const;
    T sampleGivenRelativePos(int xIdx, int yIdx, int zIdx, REAL fx, REAL fy, REAL fz) const;

    Vector3f position(int xIdx, int yIdx, int zIdx) const {return (Vector3f(xIdx, yIdx, zIdx) + Vector3f::Ones() * 0.5) * gridSpacing();}
    Vector3f position(const Size3& idx) const {return (idx.cast<REAL>() + Vector3f::Ones() * 0.5) * gridSpacing();}
    
    T get(const Size3& idx) const {return get(idx.x(), idx.y(), idx.z());}
    T get(const Size3& idx, const T& def) const {return get(idx.x(), idx.y(), idx.z(), def);}
    T get(int xIdx, int yIdx, int zIdx) const {return data[index(xIdx, yIdx, zIdx)]; }
    T get(int xIdx, int yIdx, int zIdx, const T& def) const {
        if(isExist(xIdx, yIdx, zIdx))
            return data[index(xIdx, yIdx, zIdx)];
        return def;
    }
    void set(const Size3& idx, const T& value) {data[index(idx.x(), idx.y(), idx.z())] = value;} 
    void set(int xIdx, int yIdx, int zIdx, const T& value) {data[index(xIdx, yIdx, zIdx)] = value;}
    bool isExist(int xIdx, int yIdx, int zIdx) const { return bitmap->get(xIdx, yIdx, zIdx) == true; }

    void setBitMap(BitMapPtr _bitmap) {bitmap = _bitmap;}
    void setAll(const T& v) {
        // TODO
        throw std::runtime_error("implement me!");
    }

private:
    ULL index(int xIdx, int yIdx, int zIdx) const {
        ULL ret = static_cast<ULL>(zIdx) * resXY + static_cast<ULL>(yIdx) * resX + static_cast<ULL>(xIdx);
        DEBUG_ONLY(if(xIdx < 0 || yIdx < 0 || zIdx < 0 || xIdx >= resolution().x() || yIdx >= resolution().y() || zIdx >= resolution().z()) {
            spdlog::error("[spgrid_octree_cell_centered_grid3] access out of bound. xIdx:{}, yIdx:{}, zIdx:{}", xIdx, yIdx, zIdx);
            throw std::runtime_error("[spgrid_octree_cell_centered_grid3] access out of bound.");
        });
        return ret;
    }

    T* data;
    ULL resXY, resX;
    ULL totalSize;
    BitMapPtr bitmap;
};

template<class T>
BitMapCellCenteredGrid3<T>::BitMapCellCenteredGrid3(const Size3& res, REAL spacing, BitMapPtr _bitmap) 
    : Grid3(res, spacing), bitmap(_bitmap)
{
    resXY = static_cast<ULL>(res.x()) * static_cast<ULL>(res.y());
    resX = static_cast<ULL>(res.x());

    // void *mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset);
    totalSize = resXY * static_cast<ULL>(res.z()) * sizeof(T);
    data = static_cast<T*>(mmap(0, totalSize, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0));
}

template<class T>
BitMapCellCenteredGrid3<T>::BitMapCellCenteredGrid3(const Size3& res, REAL spacing, BitMapPtr _bitmap, const T& init) 
    : Grid3(res, spacing), bitmap(_bitmap)
{
    resXY = static_cast<ULL>(res.x()) * static_cast<ULL>(res.y());
    resX = static_cast<ULL>(res.x());

    // void *mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset);
    totalSize = resXY * static_cast<ULL>(res.z()) * sizeof(T);
    data = static_cast<T*>(mmap(0, totalSize, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0));

    ThreadPool::parallelForTF(0, res.z(), [&](int z){
        for(int y = 0; y < res.y(); y++) {
            for(int x = 0; x < res.x(); x++) {
                if(isExist(x, y, z)) {
                    set(x, y, z, init);
                }
            }
        }
    });
}

template<class T>
BitMapCellCenteredGrid3<T>::BitMapCellCenteredGrid3(const Size3& res, REAL spacing, BitMapPtr _bitmap, const BitMapCellCenteredGrid3<T>& old)
    : Grid3(res, spacing), bitmap(_bitmap)
{
    resXY = static_cast<ULL>(res.x()) * static_cast<ULL>(res.y());
    resX = static_cast<ULL>(res.x());

    // void *mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset);
    totalSize = resXY * static_cast<ULL>(res.z()) * sizeof(T);
    data = static_cast<T*>(mmap(0, totalSize, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0));

    ThreadPool::parallelForTF(0, res.z(), [&](int z){
        for(int y = 0; y < res.y(); y++) {
            for(int x = 0; x < res.x(); x++) {
                if(isExist(x, y, z)) {
                    set(x, y, z, old.get(x, y ,z));
                }
            }
        }
    });
}

template<class T>
T BitMapCellCenteredGrid3<T>::sample(REAL x_pos, REAL y_pos, REAL z_pos, const T& def) const {
    FLOAT rx = x_pos / gridSpacing();
    FLOAT ry = y_pos / gridSpacing();
    FLOAT rz = z_pos / gridSpacing();
    int xIdx, yIdx, zIdx;
    REAL fx, fy, fz;

    getBarycentric(rx-0.5, 0, resolution().x()-1, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, resolution().y()-1, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, resolution().z()-1, &zIdx, &fz);

    int xIdx1 = std::min(xIdx+1, resolution().x()-1);
    int yIdx1 = std::min(yIdx+1, resolution().y()-1);
    int zIdx1 = std::min(zIdx+1, resolution().z()-1);
    
    if(!isExist(xIdx, yIdx, zIdx)) {
        return def;
    }

    T c = get(xIdx, yIdx, zIdx);
    return trilerp(c,
                    get(xIdx1, yIdx, zIdx, c),
                    get(xIdx, yIdx1, zIdx, c),
                    get(xIdx1, yIdx1, zIdx, c),
                    get(xIdx, yIdx, zIdx1, c),
                    get(xIdx1, yIdx, zIdx1, c),
                    get(xIdx, yIdx1, zIdx1, c),
                    get(xIdx1, yIdx1, zIdx1, c), fx, fy, fz);
}

template<class T>
T BitMapCellCenteredGrid3<T>::sampleGivenRelativePos(int xIdx, int yIdx, int zIdx, REAL fx, REAL fy, REAL fz) const {
    int xIdx1 = std::min(xIdx+1, resolution().x()-1);
    int yIdx1 = std::min(yIdx+1, resolution().y()-1);
    int zIdx1 = std::min(zIdx+1, resolution().z()-1);
    
    T c = get(xIdx, yIdx, zIdx);
    return trilerp(c,
                    get(xIdx1, yIdx, zIdx, c),
                    get(xIdx, yIdx1, zIdx, c),
                    get(xIdx1, yIdx1, zIdx, c),
                    get(xIdx, yIdx, zIdx1, c),
                    get(xIdx1, yIdx, zIdx1, c),
                    get(xIdx, yIdx1, zIdx1, c),
                    get(xIdx1, yIdx1, zIdx1, c), fx, fy, fz);
}

typedef BitMapCellCenteredGrid3<REAL> BitMapCellCenteredScalarGrid3;
typedef std::shared_ptr<BitMapCellCenteredScalarGrid3> BitMapCellCenteredScalarGrid3Ptr;
typedef BitMapCellCenteredGrid3<char> BitMapCellCenteredCharGrid3;
typedef std::shared_ptr<BitMapCellCenteredCharGrid3> BitMapCellCenteredCharGrid3Ptr;

} // end of namespace Omni
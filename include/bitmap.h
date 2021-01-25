#pragma once
#include <vector>
#include "general.h"
#include "grid.h"
#include "util.h"
#include "thread_pool.h"

namespace Omni {

#define ELEMENT_WIDTH (ULL)(sizeof(int)*8)

class BitMap: public Grid3 {
public:
    BitMap(const Size3& res, REAL spacing)
        : Grid3(res, spacing) {
        resXY = static_cast<ULL>(resolution().y()) * static_cast<ULL>(resolution().x());
        resX = static_cast<ULL>(resolution().x());
        totalSize = resXY * static_cast<ULL>(resolution().z());

        while(totalSize % ELEMENT_WIDTH)    // make it divided by 8
            totalSize++;
        ULL len = totalSize / ULL(ELEMENT_WIDTH);
        data.resize(len, 0);
        // spdlog::info("bit map create:res:{}, resXY:{} resX:{} {} {} {}", res.transpose(), resXY, resX, data.size(), totalSize, ELEMENT_WIDTH);
    }
    BitMap(BitMap& old) 
        : Grid3(old.resolution(), old.gridSpacing()), data(old.data.begin(), old.data.end()),
            resXY(old.resXY), resX(old.resX), totalSize(old.totalSize)
    {}

    // void parallelForEach(const LoopFunc3& func);
    // void forEach(const LoopFunc3& func);

    bool get(int xIdx, int yIdx, int zIdx) const {
        ULL idx = index(xIdx, yIdx, zIdx);
        int arrIdx = idx / ELEMENT_WIDTH;
        int pos = (ELEMENT_WIDTH-1 - idx % ELEMENT_WIDTH);
        return data[arrIdx] & (1<<pos); 
    }
    bool get(const Size3& idx) const {return get(idx.x(), idx.y(), idx.z());}
    void set(const Size3& idx, bool value) {set(idx.x(), idx.y(), idx.z(), value);} 
    void set(int xIdx, int yIdx, int zIdx, bool value) {
        ULL idx = index(xIdx, yIdx, zIdx);
        int arrIdx = idx / ELEMENT_WIDTH;
        int pos = (ELEMENT_WIDTH-1 - idx % ELEMENT_WIDTH);
        if(value) {
            data[arrIdx] = data[arrIdx] | (1<<pos);
        } else {
            data[arrIdx] = data[arrIdx] & (~(1<<pos));
        }
    }

    void fill(bool v) {
        int value = 0;
        if(v) value = 0xffffffff;
        std::fill(data.begin(), data.end(), value);
    }
    void fill(const BitMap& b) {
        data.resize(b.data.size());
        copy(b.data.begin(), b.data.end(), data.begin());
        resXY = b.resXY; resX = b.resX; totalSize = b.totalSize;
    }
    Vector3f position(int x, int y ,int z) const {return (Vector3f(x, y, z) + Vector3f::Ones() * 0.5) * gridSpacing();}

    std::vector<int>& getData() {return data;}
    Size3 arrIdxToCoordinate(int arrIdx) {
        ULL idx = static_cast<ULL>(arrIdx) * ELEMENT_WIDTH;
        return Size3((idx % resXY) % resX, (idx % resXY) / resX, idx / resXY);
    }
    
private:
    ULL index(int xIdx, int yIdx, int zIdx) const {
        ULL ret = static_cast<ULL>(zIdx) * resXY + static_cast<ULL>(yIdx) * resX + static_cast<ULL>(xIdx);
        DEBUG_ONLY(if(ret < 0 || xIdx < 0 || yIdx < 0 || zIdx < 0 || xIdx >= resolution().x() || yIdx >= resolution().y() || zIdx >= resolution().z()) {
            spdlog::error("[bitmap] access out of bound. xIdx:{}, yIdx:{}, zIdx:{}", xIdx, yIdx, zIdx);
            throw std::runtime_error("[bitmap] access out of bound.");
        });
        return ret;
    }

    std::vector<int> data;
    ULL resXY, resX, totalSize;
};

// void BitMap::parallelForEach(const LoopFunc3& func){
//     ThreadPool::parallelFor(0, resolution().z(), [&](int zIdx) {
//         for(int y = 0; y < resolution().y(); ++y) {
//             for(int x = 0; x < resolution().x(); ++x) {
//                 func(x, y, zIdx);
//             }   
//         }
//     });
// }


// void BitMap::forEach(const LoopFunc3& func){
//     for(int z = 0; z < resolution().z(); ++z) {
//         for(int y = 0; y < resolution().y(); ++y) {
//             for(int x = 0; x < resolution().x(); ++x) {
//                 func(x, y, z);
//             }   
//         }
//     }
// }

typedef std::shared_ptr<BitMap> BitMapPtr;

} // end namespace Omni
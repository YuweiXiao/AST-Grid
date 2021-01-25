#pragma once
#include "grid.h"
#include "macro.h"
#include "bitmap.h"
#include "thread_pool.h"

namespace Omni {

class OctreeGridSingleLayout3 {
public:
    OctreeGridSingleLayout3(const Size3& resolution, REAL spacing, BitMapPtr _splitMap) 
        : splitMap(_splitMap)
    {
        totalNodeNum = 0;
        auto& data = splitMap->getData();
        std::vector<int> tNum(TaskManager::tf.num_workers());
        ThreadPool::parallelForTFWithCoreID(0, (int)data.size(), [&](int i, int coreId){
            tNum[coreId] += 32 - hammingWeight(data[i]);
        });
        for(int i : tNum) {
            totalNodeNum += i;
        }
    }

    int size() const {return totalNodeNum;}
    int totalSize(bool includeGhost = true) const {return totalNodeNum;}
    bool isExist(int x, int y, int z) { return splitMap->get(x, y, z) == false;}

private:
    ULL totalNodeNum;
    BitMapPtr splitMap;
};


typedef std::shared_ptr<OctreeGridSingleLayout3> OctreeGridSingleLayout3Ptr;

} // end of namespace Omni
#include "octree_layer3.h"
#include <mutex>
#include "timer.h"

using namespace Omni;

OctreeLayer3::OctreeLayer3(const Size3& _res, REAL spacing, int l, CellCenteredFlagGrid3Ptr& region, bool singleLayer)
    : level(l), res(_res.x() >> l, _res.y() >> l, _res.z()>>l), delta_h(spacing * (1 << l)), isSingleLayer(singleLayer)
{
    Timer t;
    mapTable.max_load_factor(0.8);
    std::vector<std::vector<LCOOR_T>> tNode(TaskManager::tf.num_workers());

    ThreadPool::parallelForTFWithCoreID(0, res.z(), [&](int tz, int coreId){
        for(int ty = 0; ty < res.y(); ++ty) {
            for(int tx = 0; tx < res.x(); ++tx) {
                Vector3f pos = position(tx, ty, tz) + Vector3f::Ones() * 0.5 * delta_h;
                // auto v = region->sample(pos);
                auto v = region->get(tx<<l, ty<<l, tz<<l);
                
                if(v == l) {
                    // mapTable[PACK_COOR3(tx, ty, tz)] = LayerNode3(cIdx++, tx, ty, tz, LayerNodeFlagMask::VALID);
                    // nodeIdx.push_back(PACK_COOR3(tx, ty, tz));
                    tNode[coreId].push_back(PACK_COOR3(tx, ty, tz));
                }
            }
        }
    });
    spdlog::info("layer: {}, layer tNode parallel init: {} s", l, t.durationInSeconds()); t.reset();

    int tSize = 0;
    for(int i = 0; i < tNode.size(); ++i) {tSize += tNode[i].size();}
    mapTable.reserve(tSize);
    nodeIdx.reserve(tSize);

    for(int i = 0; i < tNode.size(); ++i) {
        nodeIdx.insert(nodeIdx.end(), tNode[i].begin(), tNode[i].end());
    }
    tNode.clear(); tNode.shrink_to_fit();
    for(int i = 0; i < nodeIdx.size(); ++i) {
        auto lCoor = nodeIdx[i];
        mapTable[lCoor] = LayerNode3(i, UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor), LayerNodeFlagMask::VALID);
    }

    // for(int tz = 0; tz < res.z(); ++tz) {
    //     for(int ty = 0; ty < res.y(); ++ty) {
    //         for(int tx = 0; tx < res.x(); ++tx) {
    //             auto v = region->sample(position(tx, ty, tz) + Vector3f::Ones() * 0.5 * delta_h);
    //             // auto v = region->get(tx<<l, ty<<l, tz<<l);
                
    //             if(v == l) {
    //                 mapTable[PACK_COOR3(tx, ty, tz)] = LayerNode3(cIdx++, tx, ty, tz, LayerNodeFlagMask::VALID);
    //                 nodeIdx.push_back(PACK_COOR3(tx, ty, tz));
    //             }
    //         }
    //     }
    // }
    spdlog::info("layer: {}, layer map init: {} s", l, t.durationInSeconds()); t.reset();
    if(isSingleLayer)
        return;
    // TODO need speed up. ghost creation don't need to check every cell
    int cIdx = nodeIdx.size();
    for(int idx = 0; idx < nodeIdx.size(); ++idx) {
        for(int i = -1; i <= 1; ++i) {
        for(int j = -1; j <= 1; ++j) {
        for(int k = -1; k <= 1; ++k) {
            if(i==0 && j==0 && k==0) 
                continue;
            int nx = UNPACK_X3(nodeIdx[idx]) + i, ny = UNPACK_Y3(nodeIdx[idx]) + j, nz = UNPACK_Z3(nodeIdx[idx]) + k;
            if(Valid(Size3(nx, ny, nz), res)) {
                // auto tv = region->sample(position(nx, ny, nz) + Vector3f::Ones() * 0.5 * delta_h);
                auto tv = region->get(nx<<l, ny<<l, nz<<l);
                LCOOR_T coor = PACK_COOR3(nx, ny, nz);
                if(tv == l-1 && mapTable.find(coor) == mapTable.end()) {
                    mapTable[coor] = LayerNode3(cIdx++, nx, ny, nz,
                                                    LayerNodeFlagMask::VALID | LayerNodeFlagMask::IS_GHOST);
                    ghostNodeIdx.push_back(coor);
                } else if(tv == l+1 && mapTable.find(coor) == mapTable.end()) {
                    mapTable[coor] = LayerNode3(cIdx++, nx, ny, nz,
                                                    LayerNodeFlagMask::VALID | LayerNodeFlagMask::IS_GHOST | LayerNodeFlagMask::GHOST_FROM_COARSE);
                    ghostNodeIdx.push_back(coor);
                }
            }
        }}}
    }

    // spdlog::info("layer: {}, layer ghost init: {} s", l, t.durationInSeconds()); t.reset();
}
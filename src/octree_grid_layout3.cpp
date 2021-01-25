#include "octree_grid_layout3.h"
#include "util.h"
#include "thread_pool.h"

using namespace Omni;

OctreeGridLayout3::OctreeGridLayout3(const Size3& res, REAL spacing, int tLevel, CellCenteredFlagGrid3Ptr splitMap, bool singleLayer) 
    : Grid3(res, spacing), totalLevel(tLevel), isSingleLayer(singleLayer), layers(tLevel), halfGhostIdx(0), TJointIdxArr(0), closedByHalfTiltArr(0)
{   
    // auto source = TaskManager::tf.placeholder();
	// auto target = TaskManager::tf.placeholder();
    for(int i = 0; i < totalLevel; ++i) {
        // auto task = TaskManager::tf.silent_emplace([&, i](){
            layers[i] = std::make_shared<OctreeLayer3>(res, spacing, i, splitMap, singleLayer);
        // });
        // source.precede(task);
		// task.precede(target);
    }
    // TaskManager::tf.wait_for_all();
}


LCOOR_T OctreeGridLayout3::levelCoordinate(REAL pos_x, REAL pos_y, REAL pos_z) const {
    int xIdx, yIdx, zIdx;
    REAL rx = pos_x / gridSpacing(), ry = pos_y / gridSpacing(), rz = pos_z / gridSpacing();
    REAL fx, fy, fz;
    getBarycentric(rx, 0, resolution().x(), &xIdx, &fx);
    getBarycentric(ry, 0, resolution().y(), &yIdx, &fy);
    getBarycentric(rz, 0, resolution().z(), &zIdx, &fz);


    for(int i = 0; i < totalLevel; ++i) {
        if(layers[i]->cellType(PACK_COOR3(xIdx, yIdx, zIdx)) == 2) {
            return PACK_LCOOR3(i, xIdx, yIdx, zIdx);
        }
        // int tmp = layers[i]->get(PACK_COOR3(xIdx, yIdx, zIdx));
        // if(tmp >= 0 && !(layers[i]->node(PACK_COOR3(xIdx, yIdx, zIdx)).flag & LayerNodeFlagMask::IS_GHOST)) {
        //     return PACK_LCOOR3(i, xIdx, yIdx, zIdx);
        // }
        xIdx >>= 1;
        yIdx >>= 1;
        zIdx >>= 1;
    }
    
    if(!isSingleLayer)
        throw std::runtime_error("[OctreeLayout3]:levelCoordinate no valid");
    return -1;
}

int OctreeGridLayout3::size() const {
    int ret = 0;
    for(auto& layer : layers)
        ret += layer->size();
    return ret;
}


LCOOR_T OctreeGridLayout3::getUntil(int l, int x, int y, int z, bool toFinest) {
    const auto& res = levelResolution(l);
    if(x < 0 || y < 0 || z < 0 || x >= res.x() || y >= res.y() || z >= res.z())
        return -1;
    while( get(l, x, y, z) == -1 || (layers[l]->node(PACK_COOR3(x, y, z)).flag & LayerNodeFlagMask::IS_GHOST)) {
        if(toFinest) {
            x <<= 1;y <<= 1;z <<= 1;
            l -= 1;
        } else {
            x >>= 1;y >>= 1;z >>= 1; 
            l += 1;
        }
        if(l < 0 || l >= totalLevel) {
            // spdlog::warn("OctreeGridLayout::getUntil fail");
            return -1;
        }
    }
    return PACK_LCOOR3(l, x, y, z);
}
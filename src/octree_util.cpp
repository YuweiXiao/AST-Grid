#include <mutex>
#include "octree_util.h"

namespace Omni {

bool checkLayoutValidity(const OctreeGridLayout3Ptr& layout) {
    for(int l = layout->level()-1; l > 0; --l) {
        Size3 lRes = layout->levelResolution(l);
        // REAL lDelta_h = layout->levelGridSpacing(l);
        for(int z = 0; z < lRes.z(); ++z) {
        for(int y = 0; y < lRes.y(); ++y) {
        for(int x = 0; x < lRes.x(); ++x) {
            if(layout->get(l, x, y, z) == -1 && layout->get(l-1, x<<1, y<<1, z<<1) >= 0) {
                for(int zz = z<<1; zz < ((z+1)<<1); ++zz) {
                for(int yy = y<<1; yy < ((y+1)<<1); ++yy) {
                for(int xx = x<<1; xx < ((x+1)<<1); ++xx) {
                    if(layout->get(l-1, xx, yy, zz) == -1) {
                        return false;
                    }
                }}}
            }
        }}}
    }

    Size3 res = layout->resolution();
    for(int z = 0; z < res.z(); ++z) {
    for(int y = 0; y < res.y(); ++y) {
    for(int x = 0; x < res.x(); ++x) {
        int c0 = layout->getUntil(0, x, y, z, false);
        if(c0 < 0)
            return false;
    }}}

    spdlog::info("layout validity check pass");
    return true;
}

void constrainOctreeTiltE3(OctreeTiltENodeGrid3Ptr& tiltEGrid, bool recordHalfTilt, bool recordTJoint, bool recordClosedByHalfTilt) {
    auto& layout = tiltEGrid->getLayout();
    std::mutex mm;
    // find level transition region & constrain tilt around T-joint position
    ThreadPool::parallelIterateGrid(tiltEGrid->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){ 
        auto tiltENode = tiltEGrid->get(iter.levelCoor());
        int minL = iter.level;
        int maxL = iter.level;
        for(int i = -1; i <= 0; ++i) {
            for(int j = -1; j <= 0; ++j) {
                for(int k = -1; k <= 0; ++k) {
                    LCOOR_T c = layout->getUntil(iter.level, iter.xIdx()+i, iter.yIdx()+j, iter.zIdx()+k, true);
                    c = c >= 0 ? c : layout->getUntil(iter.level, iter.xIdx()+i, iter.yIdx()+j, iter.zIdx()+k, false);
                    if(c > 0) {
                        minL = std::min(minL, UNPACK_LEVEL3(c));
                        maxL = std::max(maxL, UNPACK_LEVEL3(c));
                    }
                }
            }
        }

        if(minL < iter.level || (maxL > iter.level && !(iter.xIdx() & 1) && !(iter.yIdx() & 1) && !(iter.zIdx() & 1))) {
            tiltENode.is_T_joint = true;
            tiltENode.e = layout->levelGridSpacing(minL);
            tiltENode.is_fixed = true;
            tiltENode.is_closed = false;
            if(recordTJoint) {
                mm.lock();
                DEBUG_ONLY(if(std::find(layout->getTJointArr().begin(), layout->getTJointArr().end(), iter.levelCoor()) != layout->getTJointArr().end()) {throw std::runtime_error("double insert to T_joint Arr");});
                layout->getTJointArr().push_back(iter.levelCoor());
                mm.unlock();
            }
            // spdlog::info("[constrainOctreeTiltE] level:{}, xIdx, yIdx, zIdx: {}, {}, {}", iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
        }
        tiltEGrid->set(iter.levelCoor(), tiltENode);
    });

    // close tilt enforced by T_joint tilt
    iterate_grid(iter, tiltEGrid) {
        auto tiltENode = tiltEGrid->get(iter.levelCoor());
        if(!tiltENode.is_T_joint)
            continue;
        int tLevel = iter.level;
        Size3 idx(iter.xIdx(), iter.yIdx(), iter.zIdx());
        if(tiltENode.e < layout->levelGridSpacing(iter.level)) {
            tLevel -= 1;
            idx *= 2;
        }
        
        for(int i = 0; i < 6; ++i) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[i];
            if(layout->get(tLevel, tIdx.x(), tIdx.y(), tIdx.z()) >= 0) {
                auto tNode = tiltEGrid->get(tLevel, tIdx.x(), tIdx.y(), tIdx.z());
                tNode.setE(0); tNode.is_fixed = true; tNode.is_closed_by_T_joint = true;
                tiltEGrid->set(tLevel, tIdx.x(), tIdx.y(), tIdx.z(), tNode);
            }
        }
    }
    tiltEGrid->updateGhost();

    // half-tilt
    iterate_grid(iter, tiltEGrid) {
        auto tiltENode = tiltEGrid->get(iter.levelCoor());
        if(!tiltENode.is_T_joint)
            continue;
        
        bool is_coarse = tiltENode.e < layout->levelGridSpacing(iter.level);
        Size3 corseIndex = is_coarse ? Size3(iter.xIdx(), iter.yIdx(), iter.zIdx()) : 
            Size3(iter.xIdx()>>1, iter.yIdx()>>1, iter.zIdx()>>1);
        for(int axis = 0; axis < 3; ++axis) {
            if(corseIndex[axis] == 0 || corseIndex[axis] == layout->levelResolution(iter.level-is_coarse)[axis]-1)
                continue;
            bool flag = true;
            for(int i = 0; i <= 1 && flag; ++i) {
                for(int j = 0; j <= 1 && flag; ++j) {
                    if(i == 0 && j == 0) 
                        continue;
                    Size3 idx = corseIndex;
                    idx[(axis+1)%3] += i;idx[(axis+2)%3] += j;
                    LCOOR_T c = layout->getUntil(iter.level+(!is_coarse), idx.x(), idx.y(), idx.z(), true);
                    c = c >= 0 ? c : layout->getUntil(iter.level+(!is_coarse), idx.x(), idx.y(), idx.z(), false);
                    if(c >= 0) {
                        auto tNode = tiltEGrid->get(c);
                        if(tNode.is_T_joint == false)
                            flag = false;
                    } else {
                        flag = false;
                    }
                }
            }
            if(flag) {
                Size3 coarseNeighbor = corseIndex;
                coarseNeighbor(axis) -= 1;
                LCOOR_T c = layout->getUntil(iter.level+(!is_coarse), corseIndex.x(), corseIndex.y(), corseIndex.z(), true);
                LCOOR_T c1 = layout->getUntil(iter.level+(!is_coarse), coarseNeighbor.x(), coarseNeighbor.y(), coarseNeighbor.z(), true);
                ASSERT(c >= 0, "here c should be positive");
                if(c1 < 0 || UNPACK_LEVEL3(c) == UNPACK_LEVEL3(c1)) {
                    flag = false;
                }
            }
            if(flag) {
                Size3 idx = corseIndex * 2;
                idx[(axis+1)%3] += 1;idx[(axis+2)%3] += 1;
                if(layout->get(iter.level-(is_coarse), idx.x(), idx.y(), idx.z()) >= 0) {
                    auto tNode = tiltEGrid->get(iter.level-(is_coarse), idx.x(), idx.y(), idx.z());
                    tNode.is_closed = true;     // NOTE: set it true to avoid normal operation on tilt influencing half-tilt
                    tNode.is_half_tilt = true;  tNode.is_fixed = true;
                    tNode.half_direction = axis + (is_coarse ? 3 : 0);
                    tNode.e = layout->levelGridSpacing(iter.level-is_coarse);
                    // spdlog::info("half_tilt: level:{} idx:{}, {}, {}, dir:{}", iter.level-is_coarse, idx.x(), idx.y(), idx.z(), int(tNode.half_direction));
                    tiltEGrid->set(iter.level-(is_coarse), idx.x(), idx.y(), idx.z(), tNode);
                    if(recordHalfTilt) {
                        LCOOR_T halfLCoor = PACK_LCOOR3(iter.level-(is_coarse), idx.x(), idx.y(), idx.z());
                        DEBUG_ONLY(if(std::find(layout->getHalfGhostIdx().begin(), layout->getHalfGhostIdx().end(), halfLCoor) != layout->getHalfGhostIdx().end()) spdlog::error("insert same half tilt multiple time!, {}, {}, {} , {}", iter.xIdx(), iter.yIdx(), iter.zIdx(), idx.transpose());)
                        layout->getHalfGhostIdx().push_back(halfLCoor);
                    }

                    // close tilt influenced by half-tilt
                    Size3 tiltIdx = idx;   
                    tiltIdx(axis) += tNode.half_direction < 3 ? 1 : -1;
                    LCOOR_T tlCoor = PACK_LCOOR3(iter.level-(is_coarse), tiltIdx.x(), tiltIdx.y(), tiltIdx.z());
                    if(recordClosedByHalfTilt && std::find(layout->getclosedByHalfTiltArr().begin(), layout->getclosedByHalfTiltArr().end(), tlCoor) == layout->getclosedByHalfTiltArr().end()) {
                        layout->getclosedByHalfTiltArr().push_back(tlCoor);
                    }
                    auto tTiltNode = tiltEGrid->get(tlCoor);
                    tTiltNode.setE(0); tTiltNode.is_closed_by_T_joint = true; tTiltNode.half_direction = -1; tTiltNode.is_fixed = true;
                    tiltEGrid->set(tlCoor, tTiltNode);
                }
            }
        }
    }
}


Vector8LCoor tiltNeighborOctagon3(int level, int xIdx, int yIdx, int zIdx, const OctreeGridLayout3Ptr& layout) {
    Vector8LCoor ret = Vector8LCoor::Ones() * -1;
    Size3 idx(xIdx, yIdx, zIdx);
    for(int i = 0; i < 8; ++i) {
        Size3 nbCoord = idx + TILT_NB_OCT_DXDYDZ[i];
        if(!layout->isValid(level, nbCoord.x(), nbCoord.y(), nbCoord.z())) 
            continue;
        
        LCOOR_T c0 = layout->getUntil(level, nbCoord.x(), nbCoord.y(), nbCoord.z(), true);
        if(c0 < 0) 
            c0 = layout->getUntil(level, nbCoord.x(), nbCoord.y(), nbCoord.z(), false);
        if(UNPACK_LEVEL3(c0) < level) {
            int tl = UNPACK_LEVEL3(c0), tx = UNPACK_X3(c0), ty = UNPACK_Y3(c0), tz = UNPACK_Z3(c0);
            int dirm4 = i % 4;
            tz += i >= 4 ? 1 : 0;
            ty += dirm4 >= 2 ? 1: 0;
            tx += (dirm4 == 1 || dirm4 == 2) ? 1 : 0;
            c0 = PACK_LCOOR3(tl, tx, ty, tz);
        }
        ret[i] = c0;
    }
    return ret;
}

} // end of namespace Omni
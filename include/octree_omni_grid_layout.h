#pragma once
#include "general.h"
#include "octree_cell_centered_grid3.h"
#include "octree_vertex_centered_grid3.h"

namespace Omni {

template <int d>
struct CellNode {
    static_assert(d == 3);

    bool tilt_node_status_0 : 1;  // xp_yp
    bool tilt_node_status_1 : 1;  // xn_yp
    bool tilt_node_status_2 : 1;  // xn_yn
    bool tilt_node_status_3 : 1;  // xp_yn
    bool tilt_node_status_4 : 1;  // xp_yp_zn
    bool tilt_node_status_5 : 1;  // xn_yp_zn
    bool tilt_node_status_6 : 1;  // xn_yn_zn
    bool tilt_node_status_7 : 1;  // xp_yn_zn
    template <int k>
    bool getTiltNodeStatus() const {
        if constexpr (k == 0) return tilt_node_status_0;
        if constexpr (k == 1) return tilt_node_status_1;
        if constexpr (k == 2) return tilt_node_status_2;
        if constexpr (k == 3) return tilt_node_status_3;
        if constexpr (k == 4) return tilt_node_status_4;
        if constexpr (k == 5) return tilt_node_status_5;
        if constexpr (k == 6) return tilt_node_status_6;
        if constexpr (k == 7) return tilt_node_status_7;
        throw std::runtime_error("unknown tilt node index k");
    }
};

template <int d>
class OctreeOmniGridLayout {
  public:
    static_assert(d == 3);
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<Real, d, 1>;

    OctreeOmniGridLayout(const IndexT& res, Real spacing, const OctreeGridLayout3Ptr& octree_layout)
        : resolution(res), dx(spacing), layout(octree_layout) {
        tiltE = std::make_shared<OctreeTiltENodeSGrid3>(resolution, dx, layout, TiltENodeS());
    }

    bool isTiltOpen(LCOOR_T lCoor) const { return tiltE->get(lCoor).is_open; }
    bool isTiltOpen(int level, const IndexT& idx) const { return tiltE->get(level, idx).is_open; }
    TiltENodeS& getTiltNode(LCOOR_T lCoor) { return tiltE->get(lCoor); }
    void updateGhost() const { tiltE->updateGhost(); }

  private:
    OctreeGridLayout3Ptr layout;
    OctreeTiltENodeSGrid3Ptr tiltE;
    const IndexT resolution;
    const Real dx;
};

// template<int d>
// void OctreeOmniGridLayout<d>::constrainTiltToOctree(bool recordHalfTilt, bool recordTJoint, bool
// recordClosedByHalfTilt) {
//     std::mutex mm;
//     // find level transition region & constrain tilt around T-joint position
//     ThreadPool::parallelIterateGrid(tiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
//         int minL = iter.level;
//         int maxL = iter.level;
//         for(int i = -1; i <= 0; ++i) {
//             for(int j = -1; j <= 0; ++j) {
//                 for(int k = -1; k <= 0; ++k) {
//                     LCOOR_T c = layout->getUntil(iter.level, iter.xIdx()+i, iter.yIdx()+j, iter.zIdx()+k, true);
//                     c = c >= 0 ? c : layout->getUntil(iter.level, iter.xIdx()+i, iter.yIdx()+j, iter.zIdx()+k,
//                     false); if(c > 0) {
//                         minL = std::min(minL, UNPACK_LEVEL3(c));
//                         maxL = std::max(maxL, UNPACK_LEVEL3(c));
//                     }
//                 }
//             }
//         }

//         if(minL < iter.level || (maxL > iter.level && !(iter.xIdx() & 1) && !(iter.yIdx() & 1) && !(iter.zIdx() &
//         1))) {
//             auto& node = getTiltNode(iter.levelCoor());
//             node.is_open = true;
//             node.is_TJoint = true;
//             node.is_full_size = (iter.level == minL);
//             node.is_fixed = true;
//             if(recordTJoint) {
//                 mm.lock();
//                 DEBUG_ONLY(if(std::find(layout->getTJointArr().begin(), layout->getTJointArr().end(),
//                 iter.levelCoor()) != layout->getTJointArr().end()) {throw std::runtime_error("double insert to
//                 T_joint Arr");}); layout->getTJointArr().push_back(iter.levelCoor()); mm.unlock();
//             }
//             // spdlog::info("[constrainOctreeTiltE] level:{}, xIdx, yIdx, zIdx: {}, {}, {}", iter.level, iter.xIdx(),
//             iter.yIdx(), iter.zIdx());
//         }
//     });

//     // close tilt enforced by T_joint tilt
//     iterate_grid(iter, tiltE) {
//         const auto& node = tiltE->get(iter.levelCoor());
//         if(!node.is_TJoint)
//             continue;
//         int tLevel = iter.level;
//         Size3 idx(iter.xIdx(), iter.yIdx(), iter.zIdx());
//         if(!node.is_full_size) {
//             tLevel -= 1;
//             idx *= 2;
//         }

//         for(int i = 0; i < 6; ++i) {
//             Size3 tIdx = idx + OCT_NB_DXDYDZ[i];
//             if(layout->get(tLevel, tIdx.x(), tIdx.y(), tIdx.z()) >= 0) {
//                 auto& t_node = tiltE->get(tLevel, tIdx.x(), tIdx.y(), tIdx.z());
//                 t_node.is_open = false; t_node.is_fixed = true; t_node.is_closed_by_TJoint = true;
//             }
//         }
//     }
//     tiltE->updateGhost();

//     // half-tilt
//     iterate_grid(iter, tiltE) {
//         const auto& node = tiltE->get(iter.levelCoor());
//         if(!node.is_TJoint)
//             continue;

//         bool is_coarse = !node.is_full_size;
//         Size3 corseIndex = is_coarse ? Size3(iter.xIdx(), iter.yIdx(), iter.zIdx()) :
//             Size3(iter.xIdx()>>1, iter.yIdx()>>1, iter.zIdx()>>1);
//         for(int axis = 0; axis < 3; ++axis) {
//             if(corseIndex[axis] == 0 || corseIndex[axis] == layout->levelResolution(iter.level-is_coarse)[axis]-1)
//                 continue;
//             bool flag = true;
//             for(int i = 0; i <= 1 && flag; ++i) {
//                 for(int j = 0; j <= 1 && flag; ++j) {
//                     if(i == 0 && j == 0)
//                         continue;
//                     Size3 idx = corseIndex;
//                     idx[(axis+1)%3] += i;idx[(axis+2)%3] += j;
//                     LCOOR_T c = layout->getUntil(iter.level+(!is_coarse), idx.x(), idx.y(), idx.z(), true);
//                     c = c >= 0 ? c : layout->getUntil(iter.level+(!is_coarse), idx.x(), idx.y(), idx.z(), false);
//                     if(c >= 0) {
//                         if(tiltE->get(c).is_TJoint == false)
//                             flag = false;
//                     } else {
//                         flag = false;
//                     }
//                 }
//             }
//             if(flag) {
//                 Size3 coarseNeighbor = corseIndex;
//                 coarseNeighbor(axis) -= 1;
//                 LCOOR_T c = layout->getUntil(iter.level+(!is_coarse), corseIndex.x(), corseIndex.y(), corseIndex.z(),
//                 true); LCOOR_T c1 = layout->getUntil(iter.level+(!is_coarse), coarseNeighbor.x(), coarseNeighbor.y(),
//                 coarseNeighbor.z(), true); ASSERT(c >= 0, "here c should be positive"); if(c1 < 0 || UNPACK_LEVEL3(c)
//                 == UNPACK_LEVEL3(c1)) {
//                     flag = false;
//                 }
//             }
//             if(flag) {
//                 Size3 idx = corseIndex * 2;
//                 idx[(axis+1)%3] += 1;idx[(axis+2)%3] += 1;
//                 if(layout->get(iter.level-(is_coarse), idx.x(), idx.y(), idx.z()) >= 0) {
//                     auto& t_node = tiltE->get(iter.level-(is_coarse), idx.x(), idx.y(), idx.z());
//                     t_node.is_open = false;     // NOTE: set it true to avoid normal operation on tilt influencing
//                     half-tilt t_node.is_half_tilt = true;  t_node.is_fixed = true; t_node.half_tilt_direction = axis
//                     + (is_coarse ? 3 : 0); t_node.is_full_size = !is_coarse;
//                     // t_node.e = layout->levelGridSpacing(iter.level-is_coarse);
//                     // spdlog::info("half_tilt: level:{} idx:{}, {}, {}, dir:{}", iter.level-is_coarse, idx.x(),
//                     idx.y(), idx.z(), int(tNode.half_tilt_direction));
//                     // tiltE->set(iter.level-(is_coarse), idx.x(), idx.y(), idx.z(), tNode);
//                     if(recordHalfTilt) {
//                         LCOOR_T halfLCoor = PACK_LCOOR3(iter.level-(is_coarse), idx.x(), idx.y(), idx.z());
//                         DEBUG_ONLY(if(std::find(layout->getHalfGhostIdx().begin(), layout->getHalfGhostIdx().end(),
//                         halfLCoor) != layout->getHalfGhostIdx().end()) spdlog::error("insert same half tilt multiple
//                         time!, {}, {}, {} , {}", iter.xIdx(), iter.yIdx(), iter.zIdx(), idx.transpose());)
//                         layout->getHalfGhostIdx().push_back(halfLCoor);
//                     }

//                     // close tilt influenced by half-tilt
//                     Size3 tiltIdx = idx;
//                     tiltIdx(axis) += t_node.half_tilt_direction < 3 ? 1 : -1;
//                     LCOOR_T tlCoor = PACK_LCOOR3(iter.level-(is_coarse), tiltIdx.x(), tiltIdx.y(), tiltIdx.z());
//                     if(recordClosedByHalfTilt && std::find(layout->getclosedByHalfTiltArr().begin(),
//                     layout->getclosedByHalfTiltArr().end(), tlCoor) == layout->getclosedByHalfTiltArr().end()) {
//                         layout->getclosedByHalfTiltArr().push_back(tlCoor);
//                     }
//                     auto& t_node_2 = tiltE->get(tlCoor);
//                     t_node_2.is_open = false; t_node_2.is_closed_by_TJoint = true; t_node_2.half_tilt_direction = -1;
//                     t_node_2.is_fixed = true;
//                 }
//             }
//         }
//     }
// }

template <int d>
using OctreeOmniGridLayoutPtr = std::shared_ptr<OctreeOmniGridLayout<d>>;

}  // namespace Omni

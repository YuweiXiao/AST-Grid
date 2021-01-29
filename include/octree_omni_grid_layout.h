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

template <int d>
using OctreeOmniGridLayoutPtr = std::shared_ptr<OctreeOmniGridLayout<d>>;

}  // namespace Omni

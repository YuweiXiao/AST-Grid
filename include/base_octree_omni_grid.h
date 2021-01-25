#pragma once
#include "base_octree_grid.h"
#include "octree_vertex_centered_grid3.h"

namespace Omni{

class BaseOctreeOmniGrid3 : public BaseOctreeGrid3 {
  public:
    BaseOctreeOmniGrid3(const Size3& res, REAL spac, OctreeGridLayout3Ptr& layout, OctreeTiltENodeGrid3Ptr& _tiltE)
        : BaseOctreeGrid3(res, spac, layout), tiltE(_tiltE) {}

    virtual void setTiltEGridPtr(OctreeTiltENodeGrid3Ptr p) { tiltE = p; }
    OctreeTiltENodeGrid3Ptr& getTiltEGridPtr() { return tiltE; }
    const OctreeTiltENodeGrid3Ptr& getTiltEGridPtr() const { return tiltE; }

  protected:
    OctreeTiltENodeGrid3Ptr tiltE;
};

#define update_ghost_exclude_half_tilt(type, grid) \
    std::vector<std::pair<LCOOR_T, type>> tmpHalfTilt;\
    iterate_half_tilt(lCoor, layout) { \
        tmpHalfTilt.push_back(std::make_pair(lCoor, grid->get(lCoor)));\
    }\
    grid->updateGhost();\
    for(auto& p : tmpHalfTilt) {\
        grid->set(p.first, p.second);\
    }\

}   // end of namespace Omni
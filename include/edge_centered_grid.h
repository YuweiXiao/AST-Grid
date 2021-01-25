#pragma once
#include "ast_general.h"
#include "face_centered_grid.h"

AST_NAMESPACE_START

// Base case for 2 dimensional
// Reuse face centered grid because face and edge are same in 2D
template<typename T, int d> 
class EdgeCenteredGrid: public FaceCenteredGrid<T, d> { 
public: 
    using FaceCenteredGrid<T, d>::FaceCenteredGrid; // reuse base constructor
};

// TODO 3d edge centered grid
template<typename T> class EdgeCenteredGrid<T, 3> {
    static constexpr int d = 3;
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
    using VecT = Eigen::Matrix<REAL, d, 1>;
public:
    EdgeCenteredGrid(const IndexT& res, REAL delta_h, const T& init_v = Omni::zero<T>()) 
        : resolution(res), dx(delta_h)
    {
        for(int i = 0; i < d; ++i) {
            IndexT t_res = resolution + IndexT::Ones() - IndexT::Unit(i);
            data[i].resize(t_res.prod(), init_v);
        }
    }

    // access interface
    template<int axis> T& get(const IndexT& coor) { return data[axis][vecIndex<axis>(coor)]; }
    template<int axis> const T& get(const IndexT& coor) const { return data[axis][vecIndex<axis>(coor)]; }
    template<int axis> void set(const IndexT& coor, const T& v) { data[axis][vecIndex<axis>(coor)] = v; }

    // position interface
    template<int axis> PosT position(const IndexT& coor) const {
        return (coor.template cast<REAL>() + PosT::Unit(axis) * 0.5) * dx;
    }

private:
    template<int axis> int vecIndex(const IndexT& coor) const {
        static_assert(axis >= 0 && axis < d, "axis should limit in between [0, d]");
        int ret = coor.x() + coor.y() * (resolution.x() + (axis != axis_type::x)) +
            + coor.z() * (resolution.x() + (axis != axis_type::x)) * (resolution.y() + (axis != axis_type::y));
        ASSERT(ret < data[axis].size(), "[face centered grid] exceed range");
        DEBUG_ONLY( for(int i = 0; i < d; ++i) { ASSERT(coor[i] < (resolution[i] + (axis != i)), "[ast vector grid] exceed range"); } );
        return ret;
    }

    const IndexT resolution;
    const REAL dx;
    std::vector<T> data[d];
};

template<typename T, int d> using EdgeCenteredGridPtr = std::shared_ptr<EdgeCenteredGrid<T, d>>;

AST_NAMESPACE_END
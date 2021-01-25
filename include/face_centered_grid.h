#pragma once
#include "general.h"
#include "util.h"
#include "static_for.h"

namespace ast {

template<typename T, int d>
class FaceCenteredGrid {
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
    using VecT = Eigen::Matrix<REAL, d, 1>;
public:
    FaceCenteredGrid(const IndexT& res, REAL delta_h, const T& init_v = Omni::zero<T>()) 
        : resolution(res), dx(delta_h)
    {
        for(int i = 0; i < d; ++i) {
            IndexT t_res = resolution; t_res[i] += 1;
            data[i].resize(t_res.prod(), init_v);
        }
    }

    // access interface
    template<int axis> T& get(const IndexT& coor) { return data[axis][vecIndex<axis>(coor)]; }
    template<int axis> const T& get(const IndexT& coor) const { return data[axis][vecIndex<axis>(coor)]; }
    template<int axis> void set(const IndexT& coor, const T& v) { data[axis][vecIndex<axis>(coor)] = v; }
    template<int axis> std::vector<T>& getData() { return data[axis]; }
    VecT sample(const PosT& pos) const {
        PosT rp = pos / dx; VecT ret;
        static_for<0, d>()([&](auto axis){
            ret[axis] = sampleComponent<axis>(rp);
        });
        return ret;
    }

    // position interface
    template<int axis> PosT position(const IndexT& coor) const {
        PosT offset = (PosT::Ones() - PosT::Unit(axis)) * 0.5;
        return (coor.template cast<Real>() + offset) * dx;
    }

private:
    template<int axis>
    REAL sampleComponent(const PosT& rp) const {
        IndexT idx;
        PosT f;
        IndexT limit_res = resolution - IndexT::Ones() + IndexT::Unit(axis);
        Omni::getBaryCentric(rp - 0.5 * (PosT::Ones()-PosT::Unit(axis)), limit_res, idx, f);

        IndexT idx1 = (idx + IndexT::Ones()).cwiseMin(limit_res);
 
        if constexpr(d == 2) {
            if constexpr (std::is_scalar_v<T>) {
                return Omni::bilerp(get<axis>(idx), get<axis>(IndexT(idx.x(), idx1.y())), 
                    get<axis>(idx1), get<axis>(IndexT(idx1.x(), idx.y())), f.x(), f.y());
            } else {    // edge node
                return Omni::bilerp(get<axis>(idx)[axis*2], get<axis>(IndexT(idx.x(), idx1.y()))[axis*2], 
                    get<axis>(idx1)[axis*2], get<axis>(IndexT(idx1.x(), idx.y()))[axis*2], f.x(), f.y());
            }
        } else {
            if constexpr (std::is_scalar_v<T>) {
                return Omni::trilerp(get<axis>(idx), get<axis>(IndexT(idx1.x(), idx.y(), idx.z())), get<axis>(IndexT(idx.x(), idx1.y(), idx.z())), 
                                get<axis>(IndexT(idx1.x(), idx1.y(), idx.z())), get<axis>(IndexT(idx.x(), idx.y(), idx1.z())), get<axis>(IndexT(idx1.x(), idx.y(), idx1.z())),
                                get<axis>(IndexT(idx.x(), idx1.y(), idx1.z())), get<axis>(idx1), f.x(), f.y(), f.z());
            } else {
                static_assert(always_false<d>());
            }
            return 0.0;
        }
    }

    template<int axis> int vecIndex(const IndexT& coor) const {
        static_assert(axis >= 0 && axis < d, "axis should limit in between [0, d]");
        int ret = coor.x() + coor.y() * (resolution.x() + (axis == 0));
        if constexpr(d == 3)
            ret += coor.z() * (resolution.x() + (axis == 0)) * (resolution.y() + (axis == 1));
        ASSERT(ret < data[axis].size(), "[face centered grid] exceed range");
        DEBUG_ONLY( for(int i = 0; i < d; ++i) { ASSERT(coor[i] < (resolution[i] + (axis == i)), "[ast vector grid] exceed range"); } );
        return ret;
    }

    const IndexT resolution;
    const REAL dx;
    std::vector<T> data[d];
};

template<typename T, int d>
using FaceCenteredGridPtr = std::shared_ptr<FaceCenteredGrid<T, d>>;

}
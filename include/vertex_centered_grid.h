#pragma once
#include "general.h"
#include "tilt_e_node.h"
#include "thread_pool.h"
#include "util.h"

namespace Omni {

template<class T, int d>
class VertexCenteredGrid {
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<Real, d, 1>;
    using  LoopFuncT = std::function<void(const IndexT&)>;
public:
    VertexCenteredGrid(const IndexT& _res, REAL spacing, const T& init=zero<T>())
        : res(_res), dx(spacing), resPlusOne(_res+IndexT::Ones()), data(resPlusOne.prod(), init)
    {
        int v = 1;
        for(int i = 0; i < d; ++i) {
            index_multiplier(i) = v;
            v *= resPlusOne[i];
        }
    }
    ~VertexCenteredGrid() { }

    void forEach(const LoopFuncT& func) {
        for(auto iter=ast::BoxIterator<d>(IndexT::Zero(), res); iter.valid(); iter.next()) 
            func(iter.index());
    }
    void parallelForEach(const LoopFuncT& func) {
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), res, func);
    }

    T sample(const PosT& pos) const { 
        PosT f;
        IndexT idx;
        
        getBaryCentric(pos / dx, res, idx, f);

        // getBarycentric(rx, 0, resolution().x(), &xIdx, &fx);
        // getBarycentric(ry, 0, resolution().y(), &yIdx, &fy);

        // int xIdx1 = std::min(xIdx+1, resolution().x());
        // int yIdx1 = std::min(yIdx+1, resolution().y());
        IndexT idx1 = (idx + IndexT::Ones()).cwiseMin(res);

        if constexpr(d == 2) {
            return bilerp(get(idx), get(IndexT(idx.x(), idx1.y())), get(idx1), get(IndexT(idx1.x(), idx.y())), f.x(), f.y());
        } else {
            return trilerp(get(idx), get(IndexT(idx1.x(), idx.y(), idx.z())), get(IndexT(idx.x(), idx1.y(), idx.z())),
                    get(IndexT(idx1.x(), idx1.y(), idx.z())), get(IndexT(idx.x(), idx.y(), idx1.z())), get(IndexT(idx1.x(), idx.y(), idx1.z())),
                    get(IndexT(idx.x(), idx1.y(), idx1.z())), get(idx1), f.x(), f.y(), f.z());
        }
    }

    // PosT position(int xIdx, int yIdx) const { return PosT(xIdx, yIdx) * gridSpacing(); }
    PosT position(const IndexT& coor) const { return coor.template cast<REAL>() * dx; }
    // T& get(int xIdx, int yIdx) { return data[index(xIdx, yIdx)]; }
    // const T& get(int xIdx, int yIdx) const { return data[index(xIdx, yIdx)]; }
    T& get(const IndexT& coor) { return data[index(coor)];}
    const T& get(const IndexT& coor) const { return data[index(coor)];}
    void set(const IndexT& coor, const T& value) { data[index(coor)] = value; }
    // void set(const IndexT& idx, const T& value) { data[index(idx.x(), idx.y())] = value; }
    // bool isValid(int xIdx, int yIdx) const {return Valid(IndexT(xIdx, yIdx), realRes);}
    bool isValid(const IndexT& idx) const {return Valid(idx, resPlusOne);}
    const IndexT& resolution() const { return res; }
    REAL gridSpacing() const { return dx; }
    void fill(const T& v) { std::fill(data.begin(), data.end(), v); }
    void scale(Real v) { std::transform(data.begin(), data.end(), data.begin(), 
                            std::bind(std::multiplies<T>(), std::placeholders::_1, v)); }

private:
    int index(const IndexT& coor) const {
        DEBUG_ONLY(ASSERT(isValid(coor), "VertexCentered:: access out of range"));
        return coor.dot(index_multiplier);
    }

    const REAL dx;
    const IndexT res, resPlusOne;
    IndexT index_multiplier;
    std::vector<T> data;
};


template<typename T> using VertexCenteredGrid2 = VertexCenteredGrid<T, 2>;
template<typename T, int d> using VertexCenteredGridPtr = std::shared_ptr<VertexCenteredGrid<T, d>>;
typedef VertexCenteredGrid<REAL, 2> VertexCenteredScalarGrid2;
typedef std::shared_ptr<VertexCenteredScalarGrid2> VertexCenteredScalarGrid2Ptr;
typedef VertexCenteredGrid<TiltENode, 2> TiltENodeGrid2;
typedef std::shared_ptr<TiltENodeGrid2> TiltENodeGrid2Ptr;

} // end of namespace Omni
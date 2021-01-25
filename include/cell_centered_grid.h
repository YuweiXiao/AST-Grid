#pragma once
#include "general.h"
#include "util.h"
#include "thread_pool.h"
#include "box_iterator.h"

namespace Omni {

template<typename T, int d>
class CellCenteredGrid {
    static_assert(d==2 || d==3, "only support 2d or 3d");
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
    using LoopFuncT = std::function<void(const IndexT&)>;
public:
    CellCenteredGrid(const IndexT& _res, REAL spacing, const T& init=zero<T>())
        : res(_res), dx(spacing), data(res.prod(), init)
    { 
        int v = 1;
        for(int i = 0; i < d; ++i) {
            coor_multiplier(i) = v;
            v *= res[i];
        }
    }
    ~CellCenteredGrid() { }

    void forEach(const LoopFuncT& func) {
        for(auto iter = ast::BoxIterator<d>(IndexT::Zero(), res - IndexT::Ones()); iter.valid(); iter.next())
            func(iter.index());
    }
    void parallelForEach(const LoopFuncT& func) {
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), res-IndexT::Ones(), func);
    }

    // Sample value based on real position
    T sample(const PosT& pos) const { 
        PosT rp = pos/dx - PosT::Ones()*0.5, f;
        IndexT idx;
    
        getBaryCentric(rp, res-IndexT::Ones(), idx, f);

        IndexT idx1 = (idx + IndexT::Ones()).cwiseMin(res-IndexT::Ones());
        if constexpr(d == 2) {
            return bilerp(get(idx), get(IndexT(idx.x(), idx1.y())), get(idx1), get(IndexT(idx1.x(), idx.y())), f.x(), f.y());
        } else if constexpr(d == 3) {
            return trilerp(get(idx), get(IndexT(idx1.x(), idx.y(), idx.z())), get(IndexT(idx.x(), idx1.y(), idx.z())),
                    get(IndexT(idx1.x(), idx1.y(), idx.z())), get(IndexT(idx.x(), idx.y(), idx1.z())), get(IndexT(idx1.x(), idx.y(), idx1.z())),
                    get(IndexT(idx.x(), idx1.y(), idx1.z())), get(idx1), f.x(), f.y(), f.z());
        }
    }

    // Return real position given the index
    PosT position(const IndexT& coor) const { return (coor.template cast<REAL>() + PosT::Ones() * 0.5) * dx; }

    T& get(const IndexT& coor) {return data[index(coor)];}
    const T& get(const IndexT& coor) const { return data[index(coor)]; }

    // T& get(int x, int y) { return data[index(IndexT(x, y))]; }
    // const T& get(int x, int y) const { return data[index(IndexT(x, y))]; }
    void set(const IndexT& coor, const T& value) { data[index(coor)] = value; }
    void fill(const T& v) { std::fill(data.begin(), data.end(), v); }
    const IndexT& resolution() const { return res; }
    REAL gridSpacing() const { return dx; }
    bool isValid(const IndexT& coor) const { return Valid(coor, res); }
    // bool isValid(int x, int y) const { return Valid(Size2(x, y), res); }
    const std::vector<T>& getData() const { return data; }
    void scale(Real v) {
        std::transform(data.begin(), data.end(), data.begin(),
                       std::bind(std::multiplies<T>(), std::placeholders::_1, v));
    }

private:
    int index(const IndexT& coor) const {
        DEBUG_ONLY(for(int i = 0; i < d; ++i) { ASSERT(coor[i] < res(i), "CellCentered:: access out of range"); });
        return coor.dot(coor_multiplier);
    }

    const IndexT res;
    IndexT coor_multiplier;
    const REAL dx;
    std::vector<T> data;
};

template<typename T> using CellCenteredGrid2 = CellCenteredGrid<T, 2>;
template<typename T, int d> using CellCenteredGridPtr = std::shared_ptr<CellCenteredGrid<T, d>>;
typedef CellCenteredGrid<REAL, 2> CellCenteredScalarGrid2;
typedef std::shared_ptr<CellCenteredScalarGrid2> CellCenteredScalarGrid2Ptr;
typedef CellCenteredGrid<int, 2> CellCenteredFlagGrid2;
typedef std::shared_ptr<CellCenteredFlagGrid2> CellCenteredFlagGrid2Ptr;
typedef CellCenteredGrid<Vector2f, 2> CellCenteredVector2Grid2;
typedef std::shared_ptr<CellCenteredVector2Grid2> CellCenteredVector2Grid2Ptr;

} // end namespace Omni
#pragma once
#include "ast_grid.h"
#include "box_iterator.h"
#include "thread_pool.h"

namespace ast {

template<bool parallel=false, int d, typename Func>
inline void iterateCellCentered(const Eigen::Matrix<int, d, 1>& res, const Func& func) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    IndexT upper = res-IndexT::Ones();
    if constexpr (parallel == false) {
        for(auto iter=BoxIterator<d>(IndexT::Zero(), upper); iter.valid(); iter.next()) {
            func(std::integral_constant<int, node_type::cell_node>{}, iter.index());
        }
    } else {
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), upper, [&](const IndexT& coor){
            func(std::integral_constant<int, node_type::cell_node>{}, coor);
        });
    }
}

template<bool parallel=false, int d, typename Func>
inline void iterateMac(const Eigen::Matrix<int, d, 1>& res, const Func& func) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    IndexT up_x = res - IndexT::Ones();
    up_x[0] += 1;
    if constexpr(parallel == false) {
        for(auto iter=BoxIterator<d>(IndexT::Zero(), up_x); iter.valid(); iter.next())
            func(std::integral_constant<int, 0>{}, iter.index());
    } else {
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), up_x, [&](const IndexT& coor){
            func(std::integral_constant<int, 0>{}, coor);
        });
    }
    IndexT up_y = res - IndexT::Ones();
    up_y[1] += 1;
    if constexpr(parallel == false) {
        for(auto iter=BoxIterator<d>(IndexT::Zero(), up_y); iter.valid(); iter.next())
            func(std::integral_constant<int, 1>{}, iter.index());
    } else {
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), up_y, [&](const IndexT& coor){
            func(std::integral_constant<int, 1>{}, coor);
        });
    }
    if constexpr(d == 3) {
        IndexT up_z = res - IndexT::Ones();
        up_z[2] += 1;
        if constexpr(parallel == false) {
            for(auto iter=BoxIterator<d>(IndexT::Zero(), up_z); iter.valid(); iter.next())
                func(std::integral_constant<int, 2>{}, iter.index());
        } else {
            ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), up_z, [&](const IndexT& coor){
                func(std::integral_constant<int, 2>{}, coor);
            });
        }
    }
}

template<bool parallel=false, int d, typename Func>
inline void iterateASTEdge(const Eigen::Matrix<int, d, 1>& res, const Func& func) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    if constexpr(d == 2) {
        iterateMac<parallel>(res, [&](auto axis, const IndexT& coor){
            func(std::integral_constant<int, node_type::edge_node_x+axis>{}, coor);
        });
    } else {
        IndexT up_x = res - IndexT::Unit(axis_type::x);
        if constexpr(parallel == false) {
            for(auto iter=BoxIterator<d>(IndexT::Zero(), up_x); iter.valid(); iter.next()) {
                func(std::integral_constant<int, node_type::edge_node_x>{}, iter.index());
            }
        } else {
            ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), up_x, [&](const IndexT& coor){
                func(std::integral_constant<int, node_type::edge_node_x>{}, coor);
            });
        }
        IndexT up_y = res - IndexT::Unit(axis_type::y);
        if constexpr(parallel == false) {
            for(auto iter=BoxIterator<d>(IndexT::Zero(), up_y); iter.valid(); iter.next()) {
                func(std::integral_constant<int, node_type::edge_node_y>{}, iter.index());
            }
        } else {
            ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), up_y, [&](const IndexT& coor){
                func(std::integral_constant<int, node_type::edge_node_y>{}, coor);
            });
        }
        IndexT up_z = res - IndexT::Unit(axis_type::z);
        if constexpr(parallel == false) {
            for(auto iter=BoxIterator<d>(IndexT::Zero(), up_z); iter.valid(); iter.next()) {
                func(std::integral_constant<int, node_type::edge_node_z>{}, iter.index());
            }
        } else {
            ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), up_z, [&](const IndexT& coor){
                func(std::integral_constant<int, node_type::edge_node_z>{}, coor);
            });
        }
    }
}

template<int type, int d>
inline bool isSkipClosed(const ASTGridLayoutPtr<d>& layout, const Eigen::Matrix<int, d, 1>& coor) {
    if constexpr(type == node_type::cell_node) { return false; }
    if constexpr(type == node_type::vertex_node) return layout->getTiltNode(coor).open == false;
    if constexpr(is_edge_node<type>)
        return !EDGE_NODE_ADAPTIVITY || layout->template getEdgeNode<edge_node_axis(type)>(coor) == 0;
}

template<bool parallel=false, int d, typename Func>
inline void iterateASTTilt(const Eigen::Matrix<int, d, 1>& res, const Func& func) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    if constexpr(parallel == false) {
        for(auto iter=BoxIterator<d>(IndexT::Zero(), res); iter.valid(); iter.next())
            func(iter.index());
    } else {
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), res, [&](const IndexT& coor){
            func(coor);
        });
    }
}

// iterate over tilted cells and edge cells
template<bool parallel=false, int d, typename Func>
inline void iterateASTAdaptiveDoF(const Eigen::Matrix<int, d, 1>& res, const Func& func) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    iterateASTTilt<parallel>(res, [&](const IndexT& coor){
        func(std::integral_constant<int, node_type::vertex_node>{}, coor);
    });
    iterateASTEdge<parallel>(res, func);
}

// iterate all cell center, including the uniform cell, the tilted cell and the edge cell
template<bool parallel=false, int d, typename Func>
inline void iterateAST(const Eigen::Matrix<int, d, 1>& res, const Func& func) {
    iterateCellCentered<parallel>(res, func);
    iterateASTAdaptiveDoF<parallel>(res, func);
}

}
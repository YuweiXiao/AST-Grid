#pragma once
#include "ast_grid.h"
#include "ast_util.h"
#include "ast_iterator.h"
#include "global_benchmark.h"

namespace ast {
/**
 * edge cell convention:
 *  1. the edge cell is located in edge, which is colocated with the mac cell in 2 dimensional case;
 *  2. In 2D, axis for the edge cell has the same meaning as mac cell (direction of normal vector of the edge). 
 *     In 3D, the axis parameter for the edge cell means the direction of the edge.
 *  3. In 2D, when a edge cell is open, it has same position as the mac cell. So there will be duplication of 
 *     data points. We use edge cell to store the valid value.
 * 
 *  TODO maybe explicit sync edge data to mac data in 2D case, simplify implementation and unify 2D & 3D.
 */
template<typename T, int d>
class ASTFaceGrid {
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
    using EdgeDataT = Eigen::Matrix<T, d==2?4:6, 1>;
    using VectorT = Eigen::Matrix<T, d, 1>;
public:
    using TiltDataT = Eigen::Matrix<T, d==2?4:8, 1>;

    // initialize value to zero
    ASTFaceGrid(const IndexT& res, REAL _dx, const ASTGridLayoutPtr<d>& _layout) 
        : resolution(res), dx(_dx), layout(_layout)
    {
        tilt_data = std::make_shared<VertexCenteredGrid<TiltDataT, d>>(resolution, dx, TiltDataT::Zero());
        mac_data = std::make_shared<FaceCenteredGrid<Real, d>>(resolution, dx, 0);
        edge_data = std::make_shared<EdgeCenteredGrid<EdgeDataT, d>>(resolution, dx, EdgeDataT::Zero());
        dual_grid = std::make_shared<ASTGrid<VectorT, d>>(resolution, dx, layout, VectorT::Zero());
    }

    // access interface
    template<int type, int axis=0> auto& get(const IndexT& coor) {
        if constexpr(type == node_type::cell_node) return mac_data->template get<axis>(coor);
        else if constexpr(type == node_type::vertex_node) return tilt_data->get(coor);
        else if constexpr(is_edge_node<type>) return edge_data->template get<axis>(coor);
        else static_assert(always_false<type>());
    }
    template<int type, int axis=0> const auto& get(const IndexT& coor) const {
        if constexpr(type == node_type::cell_node) return mac_data->template get<axis>(coor);
        else if constexpr(type == node_type::vertex_node) return tilt_data->get(coor);
        else if constexpr(is_edge_node<type>) return edge_data->template get<axis>(coor);
        else static_assert(always_false<type>());
    }
    
    // position interface
    template<int type, int axis = 0> PosT positionCenter(const IndexT& coor) const {
        if constexpr(type == node_type::cell_node) return mac_data->template position<axis>(coor);
        else if constexpr(type == node_type::vertex_node) return tilt_data->position(coor);
        else if constexpr(is_edge_node<type>) return edge_data->template position<axis>(coor);
        else static_assert(always_false<type>());
    }
    // sample
    VectorT sample(const PosT& pos) const {
#ifdef CORRECTION_MAC
        PosT rp = pos / dx - PosT::Ones()*0.5, f;
        IndexT idx;
        getBaryCentric(rp, resolution-IndexT::Ones(), idx, f);
        VectorT ret = mac_data->sample(pos);
        if(!layout->getTiltNode(idx+IndexT::Ones()).open)
            return ret;
        auto delta = dual_grid->template get<node_type::vertex_node>(idx+IndexT::Ones())
            - mac_data->sample(positionCenter<node_type::vertex_node>(idx+IndexT::Ones()));
        return ret + 2 * delta * std::min(f.minCoeff(), (PosT::Ones()-f).minCoeff());
#else
        return dual_grid->sample(pos); 
#endif
    }
    // update dual grid
    void updateDualGrid() {
        BENCHMARK_SCOPED_TIMER_SECTION tt("update dual grid");
#ifndef CORRECTION_MAC
        iterateCellCentered<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            const auto& cell_node = layout->getCellNode(coor);
            // NOTE: the current one works, though better interpolation can be used
            VectorT v;
            static_for<0, d>()([&](auto axis){
                Real x0 = 0.5, x1 = 0.5;
                if constexpr(d==2) {
                    x0 = cell_node.template getEdgeNodeStatus<axis*2>() ? 0.25 : 0.5;
                    x1 = cell_node.template getEdgeNodeStatus<axis*2+1>() ? 0.25 : 0.5;
                }
                REAL v0 = get<node_type::cell_node, axis>(coor+IndexT::Unit(axis));
                REAL v1 = get<node_type::cell_node, axis>(coor);
                if constexpr(d==2) {
                    if(cell_node.template getEdgeNodeStatus<axis*2>()) v0 = get<node_type::edge_node_x+axis, axis>(coor+IndexT::Unit(axis))[axis*2+1];
                    if(cell_node.template getEdgeNodeStatus<axis*2+1>()) v1 = get<node_type::edge_node_x+axis, axis>(coor)[axis*2];
                }
                v[axis] = (v0*x1 + v1*x0) / (x0+x1);
            });
            dual_grid->template set<node_type::cell_node>(coor, v);
        });
#endif
        iterateASTAdaptiveDoF<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            if(isSkipClosed<type, d>(layout, coor)) return;
            VectorT v = VectorT::Zero();
            if constexpr(type == node_type::vertex_node) {
#ifdef CORRECTION_MAC
                const auto& u = tilt_data->get(coor);
                if constexpr(d == 2) {
                    v = TILT_UNIT_DIR[0] * (u(0) + u(2)) + TILT_UNIT_DIR[1] * (u(1) + u(3));
                    // VectorT mac_v(get<node_type::cell_node, 0>(coor)[mac_direction::x_positive] + 
                    //             get<node_type::cell_node, 0>(coor+IndexT(0, -1))[mac_direction::x_positive],
                    //             get<node_type::cell_node, 1>(coor)[mac_direction::y_positive] + 
                    //             get<node_type::cell_node, 1>(coor+IndexT(-1, 0))[mac_direction::y_positive]);
                    // v = 0.5 * (v - mac_v);
                    v *= 0.5;
                } else {    // TODO need code review
                    Vector4f tv = Vector4f(u(0) + u(6), u(1) + u(7), u(2) + u(4), u(3) + u(5)) * 0.5;
                    v(0) = INV_SQRT3 * v(0) - INV_SQRT3 * v(1) - INV_SQRT3 * v(2) + INV_SQRT3 * v(3);
                    v(1) = INV_SQRT3 * v(0) + INV_SQRT3 * v(1) - INV_SQRT3 * v(2) - INV_SQRT3 * v(3);
                    v(2) = INV_SQRT3 * v(0) + INV_SQRT3 * v(1) + INV_SQRT3 * v(2) + INV_SQRT3 * v(3);
                    // VectorT mac_v = mac_data->sample(positionCenter<type>(coor));
                    v = 3.0/4.0 * v;// - mac_v;
                }
#else
                bool all_nb_edge_open = false;
                if constexpr(d==2) all_nb_edge_open = layout->getTiltNode(coor).allNeighborEdgeOpen();
                if(all_nb_edge_open) {
                    v[0] = get<node_type::edge_node_x+1, 1>(coor)[mac_direction::x_negative] 
                        + get<node_type::edge_node_x+1, 1>(coor-IndexT::Unit(0))[mac_direction::x_positive];
                    v[1] = get<node_type::edge_node_x, 0>(coor)[mac_direction::y_negative]
                        + get<node_type::edge_node_x, 0>(coor-IndexT::Unit(1))[mac_direction::y_positive];
                    v *= 0.5;
                } else {
                    const auto& u = tilt_data->get(coor);
                    // TODO make it better using matrix multiplication, unify 2d & 3d
                    if constexpr(d==2) {
                        v = TILT_UNIT_DIR[0] * (u(0) + u(2)) + TILT_UNIT_DIR[1] * (u(1) + u(3));
                        v *= 0.5;
                    } else if constexpr(d==3) {
                        Vector4f nv;
                        static_for<0, 4>()([&](auto i){
                            nv[i] = (u[i] + u[oppositeTiltIdx<d>(i)]) * 0.5;
                        });
                        v(0) = INV_SQRT3 * nv(0) - INV_SQRT3 * nv(1) - INV_SQRT3 * nv(2) + INV_SQRT3 * nv(3);
                        v(1) = INV_SQRT3 * nv(0) + INV_SQRT3 * nv(1) - INV_SQRT3 * nv(2) - INV_SQRT3 * nv(3);
                        v(2) = INV_SQRT3 * nv(0) + INV_SQRT3 * nv(1) + INV_SQRT3 * nv(2) + INV_SQRT3 * nv(3);
                        v *= 3.0/4.0;  
                    }
                }
#endif
            } else if constexpr(d==2 && is_edge_node<type>) {
                const auto& u = get<type, type-node_type::edge_node_x>(coor);
                static_for<0, d*2>()([&](auto k){ v[k/2] += u[k]; });
                v *= 0.5;
            }   // TODO 3D edge node
            dual_grid->template set<type>(coor, v);
        });
    }

    const IndexT& getResolution() const { return resolution; }
    const REAL& getDx() const { return dx; }
    const ASTGridLayoutPtr<d>& getLayout() const { return layout; }
    void setLayout(const ASTGridLayoutPtr<d>& l) { layout = l; }

private:
    std::shared_ptr<VertexCenteredGrid<TiltDataT, d>> tilt_data;
    FaceCenteredGridPtr<Real, d> mac_data;
    EdgeCenteredGridPtr<EdgeDataT, d> edge_data;
    ASTGridPtr<VectorT, d> dual_grid;


    const IndexT resolution;
    const REAL dx;
    ASTGridLayoutPtr<d> layout;
};

template<typename T, int d>
using ASTFaceGridPtr = std::shared_ptr<ASTFaceGrid<T, d>>;

}
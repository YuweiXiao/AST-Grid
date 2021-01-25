#pragma once
#include "general.h"
#include "cell_centered_grid.h"
#include "vertex_centered_grid.h"
#include "edge_centered_grid.h"
#include "ast_grid_layout.h"
#include "util.h"

using namespace Omni;

namespace ast {

enum node_type {
    cell_node = 0,
    vertex_node,
    edge_node_x,
    edge_node_y,
    edge_node_z,
};

// template<int type> struct _is_edge_node : std::false_type {};
// template<> struct _is_edge_node<node_type::edge_node_x> : std::true_type {};
// template<> struct _is_edge_node<node_type::edge_node_y> : std::true_type {};
// template<> struct _is_edge_node<node_type::edge_node_z> : std::true_type {};
template<int type, class = void> struct _is_edge_node : std::false_type {};
template<int type> struct _is_edge_node<type, 
    std::enable_if_t< (type >= node_type::edge_node_x && type <= node_type::edge_node_z), void > > : std::true_type {};
template<int type> constexpr bool is_edge_node = _is_edge_node<type>::value;
constexpr int edge_node_axis(int type) { return type - node_type::edge_node_x; }


template<typename T, int d>
class ASTGrid {
    using Index = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
public:
    ASTGrid(const Index& res, REAL delta_h, const ASTGridLayoutPtr<d>& _layout, const T& init_v = Omni::zero<T>()) 
        : resolution(res), dx(delta_h), layout(_layout), half_dx(delta_h/2.0)
    {
        cell_centered_grid = std::make_shared<CellCenteredGrid<T, d>>(resolution, dx, init_v);
        vertex_centered_grid = std::make_shared<VertexCenteredGrid<T, d>>(resolution, dx, init_v);
        edge_centered_grid = std::make_shared<EdgeCenteredGrid<T, d>>(resolution, dx, init_v);
    }

    T sample3D(const PosT& pos) const {
        PosT f, f_coarse;
        Index idx;
        getBaryCentric(pos - PosT::Ones()*0.5, resolution-Index::Ones(), idx, f);
        
        Size3 idx1 = (idx + Index::Ones()).cwiseMin(resolution-Index::Ones());

        T ret = cell_centered_grid->sample(pos);
        const auto& tilt_node = layout->getTiltNode(idx+Index::Ones());
        if(!tilt_node.open) {
            return ret;
        } else {
            T diff = vertex_centered_grid->get(idx+Index::Ones()) - cell_centered_grid->sample(position<node_type::vertex_node>(idx+Index::Ones()));
            return ret + diff * 2.0 * std::min(f.minCoeff(), (PosT::Ones() - f).minCoeff());
        }
    }

    T sample2D(const PosT& pos) const {
        // get r_pos in the finest grid, which has gridSpcing half_dx and two times of the resolution
        PosT f, f_coarse;
        Index coor;
        getBaryCentric(pos/half_dx, resolution*2, coor, f);
        coor = coor.cwiseMax(Index::Ones());
        Index base_coor = (coor-Index::Ones()) / 2;
        Index tilt_node_coor = base_coor+Index::Ones();
        const auto& tilt_node = layout->getTiltNode(tilt_node_coor);
        if(!tilt_node.open)
            return cell_centered_grid->sample(pos);
        // ASSERT(offset.minCoeff() >= -1 && offset.maxCoeff() <= 0, "!!!");
        if(tilt_node.getEdgeNodeStatus(0, coor.x() & 1) && tilt_node.getEdgeNodeStatus(1, coor.y() & 1)) {
            Index coor1 = (coor + Index::Ones()).cwiseMin(resolution*2);
            // std::cout<<coor.transpose()<<" "<<tilt_node_coor.transpose()<<" tilt node: "<<tilt_node.edge_node_status_0<<' '<<tilt_node.edge_node_status_1<<' '<<tilt_node.edge_node_status_2<<' '<<tilt_node.edge_node_status_3<<std::endl;
            const T& v0 = mixedGet(coor);
            const T& v1 = mixedGet(Index(coor.x(), coor1.y()));
            const T& v2 = mixedGet(coor1);
            const T& v3 = mixedGet(Index(coor1.x(), coor.y()));
            return Omni::bilerp(v0, v1, v2, v3, f.x(), f.y());
        }

        // diff between high_res_base_coor and coor
        Index offset = base_coor * 2 + Index::Ones() - coor;
        offset = offset.cwiseMin(Index::Zero());
        f = f/2.0 - 0.5 * offset.template cast<REAL>();
        ASSERT(f.minCoeff() >= 0 && f.maxCoeff() <= 1, "???");
        Index base_coor1 = tilt_node_coor.cwiseMin(resolution - Index::Ones());
        const T& v0 = cell_centered_grid->get(base_coor);
        const T& v1 = cell_centered_grid->get(Index(base_coor.x(), base_coor1.y()));
        const T& v2 = cell_centered_grid->get(base_coor1);
        const T& v3 = cell_centered_grid->get(Index(base_coor1.x(), base_coor.y()));

        T v = Omni::bilerp(v0, v1, v2, v3, f.x(), f.y());
        T delta = vertex_centered_grid->get(tilt_node_coor) - 0.25 * (v0+v1+v2+v3);
		return v + 2 * delta * std::min(f.minCoeff(), (PosT::Ones()-f).minCoeff());
    }

    T sample(const PosT& pos) const {
        if constexpr(d == 3) return sample3D(pos);
        else return sample2D(pos);
    }

    // position interface
    PosT position(node_type type, const Index& coor) const {
        switch(type) {
            case node_type::cell_node: return position<node_type::cell_node>(coor);
            case node_type::vertex_node: return position<node_type::vertex_node>(coor);
            case node_type::edge_node_x: return position<node_type::edge_node_x>(coor);
            case node_type::edge_node_y: return position<node_type::edge_node_y>(coor);
            case node_type::edge_node_z: return position<node_type::edge_node_z>(coor);
        }
    }
    template<int type> PosT position(const Index& coor) const {
        if constexpr(type == node_type::cell_node) { return cell_centered_grid->position(coor); }
        if constexpr(type == node_type::vertex_node) { return vertex_centered_grid->position(coor); }
        if constexpr(type == node_type::edge_node_x) { return edge_centered_grid->template position<0>(coor); }
        if constexpr(type == node_type::edge_node_y) { return edge_centered_grid->template position<1>(coor); }
        if constexpr(type == node_type::edge_node_z) { return edge_centered_grid->template position<2>(coor); }
        spdlog::error("[position] unknown node type");
    }
    // data access interface
    template<int type> T& get(const Index& coor) {
        if constexpr(type == node_type::cell_node) { return cell_centered_grid->get(coor); }
        if constexpr(type == node_type::vertex_node) { return vertex_centered_grid->get(coor); }
        if constexpr(type == node_type::edge_node_x) { return edge_centered_grid->template get<0>(coor); }
        if constexpr(type == node_type::edge_node_y) { return edge_centered_grid->template get<1>(coor); }
        if constexpr(type == node_type::edge_node_z) { return edge_centered_grid->template get<2>(coor); }
    }
    template<int type> const T& get(const Index& coor) const {
        if constexpr(type == node_type::cell_node) { return cell_centered_grid->get(coor); }
        if constexpr(type == node_type::vertex_node) { return vertex_centered_grid->get(coor); }
        if constexpr(type == node_type::edge_node_x) { return edge_centered_grid->template get<0>(coor); }
        if constexpr(type == node_type::edge_node_y) { return edge_centered_grid->template get<1>(coor); }
        if constexpr(type == node_type::edge_node_z) { return edge_centered_grid->template get<2>(coor); }
    }
    // only for 2d
    const T& mixedGet(const Index& coor) const {
        if(d == 3) { throw std::runtime_error("mixedGet only support 2d"); }
        switch((coor.x()&1)+(coor.y()&1)*2) {
            case 0:
                DEBUG_ONLY(if(!layout->getTiltNode(coor/2).open) { printCallStack(); throw std::runtime_error("tilt node: access to closed DoFs");});
                return get<node_type::vertex_node>(coor/2);
            case 1:
                DEBUG_ONLY(if(!layout->template getEdgeNode<1>(coor/2)) { printCallStack(); throw std::runtime_error("edge_node_y: access to closed DoFs");});
                return get<node_type::edge_node_y>(coor/2);
            case 2:
                DEBUG_ONLY(if(!layout->template getEdgeNode<0>(coor/2)) throw std::runtime_error("edge_node_x: access to closed DoFs"););
                return get<node_type::edge_node_x>(coor/2);
            case 3: return get<node_type::cell_node>(coor/2);
            default:
                throw std::runtime_error("should not be here");
        }
    }
    template<int type> void set(const Index& coor, const T& v) {
        if constexpr(type == node_type::cell_node) { return cell_centered_grid->set(coor, v); }
        if constexpr(type == node_type::vertex_node) { return vertex_centered_grid->set(coor, v); }
        if constexpr(type == node_type::edge_node_x) { return edge_centered_grid->template set<0>(coor, v); }
        if constexpr(type == node_type::edge_node_y) { return edge_centered_grid->template set<1>(coor, v); }
        if constexpr(type == node_type::edge_node_z) { return edge_centered_grid->template set<2>(coor, v); }
    }

    const Index& getResolution() const { return resolution; }
    const REAL& getDx() const { return dx; }
    const ASTGridLayoutPtr<d>& getLayout() const { return layout; }
    void setLayout(const ASTGridLayoutPtr<d>& l) { layout = l; }
    CellCenteredGridPtr<T, d>& getCellCenteredGridPtr() { return cell_centered_grid; }

private:
    CellCenteredGridPtr<T, d> cell_centered_grid;
    VertexCenteredGridPtr<T, d> vertex_centered_grid;
    EdgeCenteredGridPtr<T, d> edge_centered_grid;

    const Index resolution;
    const REAL dx;
    const REAL half_dx;
    ASTGridLayoutPtr<d> layout;
};


template<typename T, int d> using ASTGridPtr = std::shared_ptr<ASTGrid<T, d>>;
};
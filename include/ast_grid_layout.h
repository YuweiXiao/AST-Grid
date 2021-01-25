#pragma once
#include <vector>
#include "general.h"
#include "ast_general.h"
#include "face_centered_grid.h"
#include "static_for.h"
#include "box_iterator.h"
#include "thread_pool.h"
#include "ast_constant.h"

namespace ast {

// TODO suport 3d edge node
template<int d> struct TiltNode {
    static_assert(d == 2 || d == 3, "only support 2D or 3D case");
    bool open : 1;
    bool edge_node_status_0 : 1;    // +x
    bool edge_node_status_1 : 1;    // -x
    bool edge_node_status_2 : 1;    // +y
    bool edge_node_status_3 : 1;    // -y
    // valid for 3D
    bool edge_node_status_4 : 1;    // +z
    bool edge_node_status_5 : 1;    // -z

    TiltNode() : open(false), edge_node_status_0(false), edge_node_status_1(false),
        edge_node_status_2(false), edge_node_status_3(false), edge_node_status_4(false),
        edge_node_status_5(false)
    {}
    template<int k> bool getEdgeNodeStatus() const {
        if constexpr(k == 0) return edge_node_status_0; // +x
        if constexpr(k == 1) return edge_node_status_1; // -x
        if constexpr(k == 2) return edge_node_status_2; // +y
        if constexpr(k == 3) return edge_node_status_3; // -y
        // 3D
        if constexpr(k == 4) return edge_node_status_4;
        if constexpr(k == 5) return edge_node_status_5;
    }

    bool getEdgeNodeStatus(int axis, int offset) const {
        ASSERT(axis >= 0 && axis < d && offset <= 1 && offset >= 0, "get status exceed range");
        switch(axis*2+offset) {
            case 0: return edge_node_status_0;
            case 1: return edge_node_status_1;
            case 2: return edge_node_status_2;
            case 3: return edge_node_status_3;
            case 4: return edge_node_status_4;
            case 5: return edge_node_status_5;
            default: throw std::runtime_error("unknown id");
        }
    }

    bool allNeighborEdgeOpen() const { 
        // bool ret = true;
        // static_for is still capablable of logic short-circuiting
        // static_for<0, 4>()([&](auto k){ret = ret && getEdgeNodeStatus<k>();});
        // return ret;
        return edge_node_status_0 && edge_node_status_1 && edge_node_status_2 && edge_node_status_3;
    }
};

template<int d> struct CellNode {};
template<> struct CellNode<2> {
    bool edge_node_status_0 : 1;    // +x
    bool edge_node_status_1 : 1;    // -x
    bool edge_node_status_2 : 1;    // +y
    bool edge_node_status_3 : 1;    // -y
    bool tilt_node_status_0 : 1;    // xp_yp
    bool tilt_node_status_1 : 1;    // xp_yn
    bool tilt_node_status_2 : 1;    // xn_yn
    bool tilt_node_status_3 : 1;    // xn_yp

    CellNode() : edge_node_status_0(false), edge_node_status_1(false), edge_node_status_2(false), edge_node_status_3(false), 
        tilt_node_status_0(false), tilt_node_status_1(false), tilt_node_status_2(false), tilt_node_status_3(false)
    {}

    template<int k> bool getTiltNodeStatus() const {
        if constexpr(k == 0) return tilt_node_status_0;
        if constexpr(k == 1) return tilt_node_status_1;
        if constexpr(k == 2) return tilt_node_status_2;
        if constexpr(k == 3) return tilt_node_status_3;
        throw std::runtime_error("unknown tilt node index k");
    }
    template<int k> bool getEdgeNodeStatus() const {
        if constexpr(k == 0) return edge_node_status_0; // +x
        if constexpr(k == 1) return edge_node_status_1; // -x
        if constexpr(k == 2) return edge_node_status_2; // +y
        if constexpr(k == 3) return edge_node_status_3; // -y
        throw std::runtime_error("unknown edge node index k");
    }
};

template<> struct CellNode<3> {
    // bool edge_node_status_0 : 1;    // +x
    // bool edge_node_status_1 : 1;    // -x
    // bool edge_node_status_2 : 1;    // +y
    // bool edge_node_status_3 : 1;    // -y
    bool tilt_node_status_0 : 1;    // xp_yp
    bool tilt_node_status_1 : 1;    // xp_yn
    bool tilt_node_status_2 : 1;    // xn_yn
    bool tilt_node_status_3 : 1;    // xn_yp
    bool tilt_node_status_4 : 1;    // xp_yp_zn
    bool tilt_node_status_5 : 1;    // xp_yn_zn
    bool tilt_node_status_6 : 1;    // xn_yn_zn
    bool tilt_node_status_7 : 1;    // xn_yp_zn
    template<int k> bool getTiltNodeStatus() const {
        if constexpr(k == 0) return tilt_node_status_0;
        if constexpr(k == 1) return tilt_node_status_1;
        if constexpr(k == 2) return tilt_node_status_2;
        if constexpr(k == 3) return tilt_node_status_3;
        if constexpr(k == 4) return tilt_node_status_4;
        if constexpr(k == 5) return tilt_node_status_5;
        if constexpr(k == 6) return tilt_node_status_6;
        if constexpr(k == 7) return tilt_node_status_7;
        throw std::runtime_error("unknown tilt node index k");
    }
};

static_assert(sizeof(TiltNode<2>) == 1, "should limit TiltNode<2> size to single byte");
static_assert(sizeof(TiltNode<3>) == 1, "should limit TiltNode<3> size to single byte");
static_assert(sizeof(CellNode<2>) == 1, "should limit CellNode<2> size to single byte");
static_assert(sizeof(CellNode<3>) == 1, "should limit CellNode<3> size to single byte");

template<int d>
class ASTGridLayout {
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
public:
    ASTGridLayout(const IndexT& res, REAL _dx)
        : resolution(res), dx(_dx) {
        edge_node_ptr = std::make_shared<EdgeCenteredGrid<char, d>>(resolution, dx, 0);
        tilt_nodes.resize((resolution+IndexT::Ones()).prod());
        cell_nodes.resize(resolution.prod());
    }

    // [done] access(r/w) tilt_nodes with coordinate
    TiltNode<d>& getTiltNode(const IndexT& coor) { return tilt_nodes[indexTiltNode(coor)]; }
    const TiltNode<d>& getTiltNode(const IndexT& coor) const { return tilt_nodes[indexTiltNode(coor)]; }
    CellNode<d>& getCellNode(const IndexT& coor) { return cell_nodes[indexCellNode(coor)]; }
    const CellNode<d>& getCellNode(const IndexT& coor) const { return cell_nodes[indexCellNode(coor)]; }
    // [done] access(r/w) edge_node with coordinate
    template<int axis> char getEdgeNode(const IndexT& coor) { return edge_node_ptr->template get<axis>(coor); }
    template<int axis> void setEdgeNode(const IndexT& coor, char v) { return edge_node_ptr->template set<axis>(coor, v); }
    void clearEdgeNode() {
        static_for<0, d>()([&](auto k){
            auto& data = edge_node_ptr->template getData<k>();
            std::fill(data.begin(), data.end(), 0);
        });
    }
    // position interface
    PosT positionTilt(const IndexT& coor) const { return coor.template cast<REAL>() * dx; }
    template<int axis> PosT positionEdgeNode(const IndexT& coor) const { return edge_node_ptr->template position<axis>(coor); }
    // fixed node related interface
    const std::vector<int>& getFixedTiltNode() const { return fixed_closed_tilt_node_idx; }
    void addFixedTiltNode(const IndexT& coor) { fixed_closed_tilt_node_idx.push_back(indexTiltNode(coor)); }
    void addFixedTiltNode(const std::vector<int>& fixed_nodes) { 
        fixed_closed_tilt_node_idx.insert(fixed_closed_tilt_node_idx.end(), fixed_nodes.begin(), fixed_nodes.end()); 
    }
    void constrainFixedTiltNode() {
        for(int idx: fixed_closed_tilt_node_idx) tilt_nodes[idx].open = false;
        // static_for<0, d>()([&](auto k){
        //     for(int idx: fixed_closed_edge_node_idx[k]) 
        //         edge_node_ptr->template getData<k>()[idx] = 0;
        // });
    }

    // update cached info stored on tilt and cell nodes.
    void propogateNodeInfo() {
        // cell node: get infos from tilt and edge nodes
        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), resolution-IndexT::Ones(), [&](const IndexT& coor){
        // for(auto iter=BoxIterator<d>(IndexT::Zero(), resolution-IndexT::Ones()); iter.valid(); iter.next()) {
        //     const auto& coor = iter.index();
            auto& node = getCellNode(coor);
            if constexpr(d == 2) {
                node.edge_node_status_0 = edge_node_ptr->template get<0>(coor+IndexT::Unit(0)); // +x
                node.edge_node_status_1 = edge_node_ptr->template get<0>(coor);                 // -x
                node.edge_node_status_2 = edge_node_ptr->template get<1>(coor+IndexT::Unit(1)); // +y
                node.edge_node_status_3 = edge_node_ptr->template get<1>(coor);                 // -y
            }
            node.tilt_node_status_0 = getTiltNode(coor + kOctNbTiltOffset<d, 0>()).open;  // xp_yp
            node.tilt_node_status_1 = getTiltNode(coor + kOctNbTiltOffset<d, 1>()).open;  // xp_yn
            node.tilt_node_status_2 = getTiltNode(coor + kOctNbTiltOffset<d, 2>()).open;  // xn_yn
            node.tilt_node_status_3 = getTiltNode(coor + kOctNbTiltOffset<d, 3>()).open;  // xn_yp
            if constexpr(d == 3) {
                node.tilt_node_status_4 = getTiltNode(coor + kOctNbTiltOffset<d, 4>()).open;  // xp_yp
                node.tilt_node_status_5 = getTiltNode(coor + kOctNbTiltOffset<d, 5>()).open;  // xp_yn
                node.tilt_node_status_6 = getTiltNode(coor + kOctNbTiltOffset<d, 6>()).open;  // xn_yn
                node.tilt_node_status_7 = getTiltNode(coor + kOctNbTiltOffset<d, 7>()).open;  // xn_yp
            }
        }
        );
        // tilt node: get infos from edge nodes

        ThreadPool::parallelForBoxIteratorTF<d>(IndexT::Zero(), resolution, [&](const IndexT& coor){
        // for(auto iter=BoxIterator<d>(IndexT::Zero(), resolution); iter.valid(); iter.next()) {
            // const auto& coor = iter.index();
            auto& node = tilt_nodes[indexTiltNode(coor)];
            if(node.open) {
                if constexpr(d == 2) {
                    node.edge_node_status_0 = edge_node_ptr->template get<1>(coor);                 // +x
                    node.edge_node_status_1 = edge_node_ptr->template get<1>(coor-IndexT::Unit(0)); // -x
                    node.edge_node_status_2 = edge_node_ptr->template get<0>(coor);                 // +y
                    node.edge_node_status_3 = edge_node_ptr->template get<0>(coor-IndexT::Unit(1)); // -y
                    // node.edge_node_status_4 : 1;    // +z
                    // node.edge_node_status_5 : 1;    // -z
                }
                // TODO 3d edge node
                // if(node.edge_node_status_0 || node.edge_node_status_1 || node.edge_node_status_2 || node.edge_node_status_3) {
                //     std::cout<<coor.transpose()<<" tilt status: "<<node.edge_node_status_0<<' '<<node.edge_node_status_1<<' '<<node.edge_node_status_2<<' '<<node.edge_node_status_3<<std::endl;
                // }
            }
        }
        );
    }

    const IndexT& getResolution() const { return resolution; }
    REAL gridSpacing() const { return dx; }

private:
    int indexTiltNode(const IndexT& coor) {
        int ret = coor.x() + coor.y() * (resolution.x()+1);
        if constexpr(d == 3)
            ret += coor.z() * (resolution.x()+1) * (resolution.y()+1);
        ASSERT(ret < tilt_nodes.size(), "layout: tilt node exceeds range");
        DEBUG_ONLY( for(int i = 0; i < d; ++i) {ASSERT(coor[i] <= resolution[i], "[tilt node access] exceed range");} );
        return ret;
    }
    int indexCellNode(const IndexT& coor) {
        int ret = coor.x() + coor.y() * resolution.x();
        if constexpr(d == 3)
            ret += coor.z() * resolution.x() * resolution.y();
        ASSERT(ret < cell_nodes.size(), "layout: cell node exceeds range");
        DEBUG_ONLY( for(int i = 0; i < d; ++i) {ASSERT(coor[i] < resolution[i], "[cell node access] exceed range");} );
        return ret;
    }

    std::vector<int> fixed_closed_tilt_node_idx;
    // std::vector<int> fixed_closed_edge_node_idx[d];
    std::vector<TiltNode<d>> tilt_nodes;
    std::vector<CellNode<d>> cell_nodes;
    // NOTE: cannot use bool type. 
    // https://en.cppreference.com/w/cpp/container/vector_bool
    // https://stackoverflow.com/questions/8399417/why-vectorboolreference-doesnt-return-reference-to-bool
    std::shared_ptr<EdgeCenteredGrid<char, d>> edge_node_ptr;
    const IndexT resolution;
    const REAL dx;
};

template<int d>
using ASTGridLayoutPtr = std::shared_ptr<ASTGridLayout<d>>;

}
#pragma once
#include "ast_grid.h"
#include "ast_face_grid.h"
#include "ast_iterator.h"
#include "ast_util.h"
#include "static_for.h"

namespace ast {


template<class VType, class PType>
inline PType backTrace(const VType& flow, const PType& p, REAL dt) {
    PType midPt = p - flow->sample(p) * dt * 0.5;
    PType ret = p - flow->sample(midPt) * dt;
    return ret;
}

template<typename T, int d>
void semiLagrangianAST(const ASTFaceGridPtr<REAL, d>& vel, ASTGridPtr<T,d>& dst, const ASTGridPtr<T,d>& src, REAL dt) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    iterateAST<ENABLE_PARALLEL>(dst->getResolution(), [&](auto type, const IndexT& coor){
        if(isSkipClosed<type, d>(dst->getLayout(), coor))
            return;
        
        auto pos = dst->template position<type>(coor);
        auto old_pos = backTrace(vel, pos, dt);
        dst->template set<type>(coor, src->sample(old_pos));
    });
}

template<typename T, int d>
void semiLagrangianASTVector(const ASTFaceGridPtr<REAL, d>& vel, ASTFaceGridPtr<T,d>& dst, const ASTFaceGridPtr<T,d>& src, REAL dt) {
    using IndexT = Eigen::Matrix<int, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
    using TiltVelDataT = typename ASTFaceGrid<Real, d>::TiltDataT;
    const REAL half_dx = dst->getDx() * 0.5;
    const REAL quarter_dx = half_dx * 0.5;
    static constexpr int kTiltDimension = TiltVelDataT::RowsAtCompileTime;
    auto& layout = dst->getLayout();

    iterateASTAdaptiveDoF<ENABLE_PARALLEL>(dst->getResolution(), [&](auto type, const IndexT& coor){
        if(isSkipClosed<type>(layout, coor)) return;
        if constexpr(type == node_type::vertex_node) {
            auto& v = dst->template get<type>(coor);
            auto pos_center = dst->template positionCenter<type>(coor);
            static_for<0, kTiltDimension>()([&](auto k){
                PosT pos = pos_center + kTiltPosOffset<d, k>() * half_dx;
                auto old_pos = backTrace(vel, pos, dt);
                v[k] = src->sample(old_pos).dot(kTiltUnitUpDir<d, k>());
            });
        } else if constexpr(EDGE_NODE_ADAPTIVITY && is_edge_node<type>) {
            // TODO 3d edge node
            constexpr auto axis = edge_node_axis(type);
            auto& edge_node = dst->template get<type, axis>(coor);
            auto pos_center = dst->template positionCenter<type, axis>(coor);
            static_for<0, d*2>()([&](auto k){
                PosT pos = pos_center + PosT::Unit(k/2) * quarter_dx * ((k&1)?-1:1);
                auto old_pos = backTrace(vel, pos, dt);
                edge_node[k] = src->sample(old_pos)[k/2];
            });
        }
    });
    
    iterateMac<ENABLE_PARALLEL>(dst->getResolution(), [&](auto axis, const IndexT& coor){
        auto pos_center = dst->template positionCenter<node_type::cell_node, axis>(coor);
        auto old_pos = backTrace(vel, pos_center, dt);
        dst->template get<node_type::cell_node, axis>(coor) = src->sample(old_pos)[axis];
    });
}

    
}
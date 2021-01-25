#pragma once

#include "solver.h"
#include "ast_grid.h"
#include "ast_face_grid.h"
#include "ast_iterator.h"
#include "ast_advection.h"
#include "ast_util.h"
#include "ast_grid_pool.h"
#include "viewer/image_dump.h"
#include "pcgsolver/pcg_solver.h"
#ifdef CUDA_SOLVER
#include "cuda_solver.cuh"
#endif

namespace ast {
    
using namespace Omni;

template<int d>
class ASTSolver : public Solver {
    using IndexT = Eigen::Matrix<int, d, 1>;
    using VelocityT = Eigen::Matrix<REAL, d, 1>;
    using PosT = Eigen::Matrix<REAL, d, 1>;
    using TiltVelDataT = typename ASTFaceGrid<Real, d>::TiltDataT;
    static constexpr int kTiltDimension = TiltVelDataT::RowsAtCompileTime;
public:
    ASTSolver(const IndexT& res, REAL _dx, bool adaptive_dofs = true, const std::string& name="ASTSolver")
        : Solver(name), resolution(res), dx(_dx), half_dx(_dx/2.0), quarter_dx(_dx/4.0), allow_adaptive_dofs(adaptive_dofs)
    {
        if constexpr(d==2) tilt_oct_dx = half_dx * SQRT2;
        else tilt_oct_dx = half_dx * SQRT3;
    }

    void initialize() override {
        is_init = true;
        max_cfl = 2.0;
#ifndef CUDA_SOLVER
        // if(allow_adaptive_dofs)
        pcgSolver.set_solver_parameters(SOLVER_ACCURACY, 5000, 0.97, 0.95);
        // else 
        // pcgSolver.set_solver_parameters(SOLVER_ACCURACY, 5000);
#endif
        initData();
        initObstacle();
        
        updateIndex();
        reuse_octagon_pIndex = true;
        // initCuttingCellWeights();
    }
    int numberOfSubTimeSteps(REAL dt) const override {
        int step = std::ceil(cfl(dt) / max_cfl);
        return step > 1 ? step : 1;
    }
    int cfl(REAL dt) const {
        // TODO parallel reduce
        REAL max_squared_v = 0;
        iterateCellCentered(resolution, [&](auto type, const IndexT& coor){
            VelocityT v = VelocityT::Zero();
            static_for<0, d>()([&](auto axis){ v[axis] = velocity->template get<node_type::cell_node, axis>(coor); });
            if constexpr(d == 2 && EDGE_NODE_ADAPTIVITY) {  // in 2D, the active edge node will replace the original mac face
                static_for<0, d>()([&](auto axis){
                    if(!isSkipClosed<node_type::edge_node_x+axis>(layout, coor)) {
                        v[axis] = velocity->template get<node_type::edge_node_x+axis, axis>(coor)[axis*2];
                    }
                });
            }
            max_squared_v = std::max(max_squared_v, v.squaredNorm());
        });
        spdlog::info("max v: {}", sqrt(max_squared_v));
        return dt * sqrt(max_squared_v) / dx;
    }
    const std::vector<Real>& getDensity() const { return density->getCellCenteredGridPtr()->getData(); }

protected:
    virtual void initData() override {
        data_pool = std::make_shared<ASTGridPool<d>>(resolution, dx);
        layout = std::make_shared<ASTGridLayout<d>>(resolution, dx);
        layout_tmp = std::make_shared<ASTGridLayout<d>>(resolution, dx);
        density = data_pool->template getGrid<grid_type::ast>(layout);
        heat = data_pool->template getGrid<grid_type::ast>(layout);
        velocity = data_pool->template getGrid<grid_type::ast_face>(layout);
        pIndex = std::make_shared<ASTGrid<int, d>>(resolution, dx, layout, 0);

        // constrain outmost tilt cell in closed state
        iterateASTTilt(resolution, [&](const IndexT& coor){
            if(coor.minCoeff() == 0) layout->addFixedTiltNode(coor);
            else {
                static_for<0, d>()([&](auto k) {
                    if(coor[k] == resolution[k]) layout->addFixedTiltNode(coor);
                });
            }
        });
        layout_tmp->addFixedTiltNode(layout->getFixedTiltNode());
    }
    virtual void initObstacle() override {
        // TODO add sphere for 3d
        if constexpr(d == 3) {

        }
    }
    virtual void onUpdate(REAL dt) override {
        if constexpr (d == 2) sphereSource2DSetup();
        else sphereSource3DSetup();
        addForce(dt);
        project(dt);
        if constexpr(d==2) { if(start_of_frame) { saveVortex(cFrame, -30, 30); } }
        advect(dt);
        start_of_frame = false;
    }
    virtual void saveResult(int step) override {
        if constexpr(d==3) { 
            // TODO 3d save result
            return ; 
        } else {
            BENCHMARK_SCOPED_TIMER_SECTION tSolve("save_result");
            Viewer::grid2PNGImage("./"+name+"/density/d_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                    resolution * 2, [&](int xIdx, int yIdx) {
                REAL tmp = density->sample(Vector2f(xIdx+0.5, yIdx+0.5) * half_dx);
                return std::max(0.0, std::min(tmp, 1.0));
            });
            if(allow_adaptive_dofs) {
                Viewer::grid2PNGImage("./"+name+"/tilt_vis/t_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                        resolution + IndexT::Ones(), [&](int xIdx, int yIdx) {
                    return layout->getTiltNode(IndexT(xIdx, yIdx)).open ? 0.5 : 0;
                });
                Viewer::grid2PNGImage("./"+name+"/edge_dof_vis/ex_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                        resolution, [&](int xIdx, int yIdx) {
                    return layout->template getEdgeNode<0>(IndexT(xIdx, yIdx)) ? 1.0 : 0;
                });
                Viewer::grid2PNGImage("./"+name+"/edge_dof_vis/ey_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                        resolution, [&](int xIdx, int yIdx) {
                    return layout->template getEdgeNode<1>(IndexT(xIdx, yIdx)) ? 1.0 : 0;
                });
                Viewer::grid2PNGImage("./"+name+"/dof_vis/dof_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                        resolution, [&](int xIdx, int yIdx) {
                    if(layout->template getEdgeNode<1>(IndexT(xIdx, yIdx)) || layout->template getEdgeNode<0>(IndexT(xIdx, yIdx)))
                        return 1.0;
                    return layout->getTiltNode(IndexT(xIdx, yIdx)).open ? 0.5 : 0;
                });
            }
    #if DEBUG_VELOCITY
            Viewer::grid2PNGImage("./"+name+"/velocity/v_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                    resolution * 2, [&](int xIdx, int yIdx) {
                return std::min(1.0, velocity->sample(Vector2f(xIdx+0.5, yIdx+0.5) * half_dx).norm() / 2);
            });
            // Viewer::grid2PNGImage("./"+name+"/velocity/x_"+paddingStr(std::to_string(step), '0', 4)+".png", 
            //         resolution, [&](int xIdx, int yIdx) {
            //     return std::min(1.0, fabs(velocity->template getMac<0>(Size2(xIdx, yIdx))[mac_direction::x_positive]));
            // });
            // Viewer::grid2PNGImage("./"+name+"/velocity/y_"+paddingStr(std::to_string(step), '0', 4)+".png", 
            //         resolution, [&](int xIdx, int yIdx) {
            //     return std::min(1.0, fabs(velocity->template getMac<1>(Size2(xIdx, yIdx))[mac_direction::y_positive]));
            // });
    #endif
        }

        // TODO dump density, layout
    }

    virtual void endOfFrame() override {
        start_of_frame = true;
    }

private:
    void sphereSource2DSetup() {
        BENCHMARK_SCOPED_TIMER_SECTION tt("scene setup");
        Vector2f sphere_center(0.5, 0.04);
        REAL smoke_radius = 0.03;
        REAL init_v = 1.0;
        Vector4f tilt_init_v = Vector4f(INV_SQRT2, -INV_SQRT2, INV_SQRT2, -INV_SQRT2) * init_v;
        
        iterateAST<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            auto pos = density->template position<type>(coor);
            if((pos-sphere_center).norm() < smoke_radius) {
                density->template set<type>(coor, 1);
                heat->template set<type>(coor, 1);
            }
        });
        // initial velocity setup
        iterateASTAdaptiveDoF<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            if(isSkipClosed<type>(layout, coor)) return;
            if constexpr(type == node_type::vertex_node) {
                auto pos = velocity->template positionCenter<node_type::vertex_node>(coor);
                if((pos-sphere_center).norm() < smoke_radius)
                    velocity->template get<node_type::vertex_node>(coor) = tilt_init_v;
            } else if constexpr(is_edge_node<type>) {
                constexpr int axis = edge_node_axis(type);
                auto pos = velocity->template positionCenter<type, axis>(coor);
                static_for<0, 2>()([&](auto k){
                    PosT n_pos = pos + PosT::Unit(axis_type::y) * (k==0?1:-1) * quarter_dx;
                    if((n_pos - sphere_center).norm() < smoke_radius) {
                        velocity->template get<type, axis>(coor)[mac_direction::y_positive+k] = init_v;
                    }
                });
            }
        });
        iterateMac<ENABLE_PARALLEL>(resolution, [&](auto axis, const IndexT& coor){
            if constexpr(axis == axis_type::y) {   // TODO maybe only iterate through y axis
                auto pos = velocity->template positionCenter<node_type::cell_node, axis>(coor);
                if((pos-sphere_center).norm() < smoke_radius)
                    velocity->template get<node_type::cell_node, axis>(coor) = init_v;
            }
        });
    }

    void sphereSource3DSetup() {    // TODO, check it with omniSmoke3D setup
        Vector3f sphere_center = Vector3f(0.5, 0.5, 0.07);
        Real smoke_radius = 0.06;
        IndexT constrained_res = resolution;
        constrained_res.z() = int((sphere_center.z() + smoke_radius) / dx) + 5;
        Real init_v = 1.0;
        TiltVelDataT tilt_init_v = TiltVelDataT::Ones() * INV_SQRT3 * init_v;

        iterateAST<ENABLE_PARALLEL>(constrained_res, [&](auto type, const IndexT& coor){
            auto pos = density->template position<type>(coor);
            if((pos-sphere_center).norm() < smoke_radius) {
                density->template set<type>(coor, 1);
                heat->template set<type>(coor, 1);
            }
        });
        iterateASTAdaptiveDoF<ENABLE_PARALLEL>(constrained_res, [&](auto type, const IndexT& coor){
            if(isSkipClosed<type>(layout, coor)) return;
            if constexpr(type == node_type::vertex_node) {
                auto pos = velocity->template positionCenter<type>(coor);
                if((pos-sphere_center).norm() < smoke_radius)
                    velocity->template get<node_type::vertex_node>(coor) = tilt_init_v;
            } else if constexpr(is_edge_node<type>) {
                constexpr int axis = edge_node_axis(type);
                auto pos = velocity->template positionCenter<type, axis>(coor);
                static_for<0, 2>()([&](auto k){
                    PosT n_pos = pos + PosT::Unit(axis_type::z) * (k==0?1:-1) * quarter_dx;
                    if((n_pos - sphere_center).norm() < smoke_radius) {
                        velocity->template get<type, axis>(coor)[mac_direction::z_positive+k] = init_v;
                    }
                });
            }
        });
        iterateMac<ENABLE_PARALLEL>(constrained_res, [&](auto axis, const IndexT& coor){
            if constexpr(axis == axis_type::z) {
                auto pos = velocity->template positionCenter<node_type::cell_node, axis>(coor);
                if((pos-sphere_center).norm() < smoke_radius)
                    velocity->template get<node_type::cell_node, axis>(coor) = init_v;
            }
        });
    }

    void addForce(REAL dt) {
        BENCHMARK_SCOPED_TIMER_SECTION tt("add force");
        const REAL beta = 0.5;
        iterateASTAdaptiveDoF<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            if(isSkipClosed<type>(layout, coor)) return;
            if constexpr(type == node_type::vertex_node) {
                TiltVelDataT v;
                auto pos_center = velocity->template positionCenter<node_type::vertex_node>(coor);
                static_for<0, kTiltDimension>()([&](auto k){
                    v[k] = heat->sample(pos_center + kTiltPosOffset<d, k>() * half_dx)
                        * (d == 3 ? INV_SQRT3 : (k&1?-INV_SQRT2:INV_SQRT2));
                });
                velocity->template get<node_type::vertex_node>(coor) += v * beta * dt;
            } else if constexpr(is_edge_node<type>) {   // TODO: check 3D node case
                constexpr int up_axis = d == 2 ? axis_type::y : axis_type::z;
                auto& edge_node = velocity->template get<type, edge_node_axis(type)>(coor);
                auto pos_center = velocity->template positionCenter<type, edge_node_axis(type)>(coor);
                edge_node[up_axis*2] += beta * dt * heat->sample(pos_center+PosT::Unit(up_axis)*quarter_dx);
                edge_node[up_axis*2+1] += beta * dt * heat->sample(pos_center-PosT::Unit(up_axis)*quarter_dx);
            }
        });
        iterateMac<ENABLE_PARALLEL>(resolution, [&](auto axis, const IndexT& coor){
            if constexpr((d == 2 && axis == axis_type::y) || (d == 3 && axis == axis_type::z)) {
                auto pos_center = velocity->template positionCenter<node_type::cell_node, axis>(coor);
                velocity->template get<node_type::cell_node, axis>(coor) += heat->sample(pos_center) * beta * dt;
            }
        });
        
        constrainVelocity();
    }

    void correctVelocityByPressure(REAL dt) {
        BENCHMARK_SCOPED_TIMER_SECTION tt("correct velocity");
        iterateASTAdaptiveDoF<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            if(isSkipClosed<type>(layout, coor)) return;
            if constexpr(type == node_type::vertex_node) {
                int cIdx = pIndex->template get<type>(coor);
                if(cIdx < 0) return;
                auto& u = velocity->template get<node_type::vertex_node>(coor);
                static_for<0, kTiltDimension>()([&](auto k){
                    IndexT nbCoor = coor + kTiltNbOctOffset<d, k>();
                    int tIdx = pIndex->template get<node_type::cell_node>(nbCoor);
                    u[k] += ((k<kTiltDimension/2)?-1:1) * dt / rho * (vRes[tIdx] - vRes[cIdx]) / tilt_oct_dx;
                });
            } else if constexpr(is_edge_node<type>) {
                // TODO edge node 3d case
                constexpr int axis = type - node_type::edge_node_x;
                auto& u = velocity->template get<node_type::edge_node_x+axis, axis>(coor);
                int cIdx = pIndex->template get<type>(coor);
                static_for<0, 2>()([&](auto k){
                    IndexT nbCoor = coor; nbCoor[axis] -= k;
                    u[axis*2+k] -= (k==0?1:-1) * dt / rho * (vRes[pIndex->template get<node_type::cell_node>(nbCoor)] - vRes[cIdx]) / half_dx;
                });
                constexpr int new_axis = axis^1;
                static_for<0, 2>()([&](auto k){
                    IndexT nbCoor = coor; nbCoor[new_axis] += 1 - k;
                    u[new_axis*2+k] -= (k==0?1:-1) * dt / rho * (vRes[pIndex->template get<node_type::vertex_node>(nbCoor)] - vRes[cIdx]) / half_dx;
                });
            } else { static_assert(always_false<type>(), "invalid type"); }
        });
        iterateMac<ENABLE_PARALLEL>(resolution, [&](auto axis, const IndexT& coor){
            REAL p1 = 0, p0 = 0;
            if(Valid(coor, resolution)) {
                int tmp_idx = pIndex->template get<node_type::cell_node>(coor);
                if(tmp_idx >= 0) p1 = vRes[tmp_idx];
            }
            if(Valid(coor-IndexT::Unit(axis), resolution)) {
                int tmp_idx = pIndex->template get<node_type::cell_node>(coor-IndexT::Unit(axis));
                if(tmp_idx >= 0) p0 = vRes[tmp_idx];
            }
            velocity->template get<node_type::cell_node, axis>(coor) -= dt / rho * (p1 - p0) / dx;
        });
    }

    void project(REAL dt) {
        BENCHMARK_SCOPED_TIMER_SECTION tt("project");
        updateIndex();
        // construct A
        A.resize(total_index_size); A.zero();
        b.resize(total_index_size); vRes.resize(total_index_size);
        iterateCellCentered<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            // assume all cell node is active DoFs
            const auto& cell_node = layout->getCellNode(coor);
            REAL center = 0, tb = 0;
            int cIdx = pIndex->template get<node_type::cell_node>(coor);
            REAL face_area[d*2] = {0};
            // neighbor tilt node
            static_for<0, kTiltDimension>()([&](auto k){
                if constexpr(d==2) {    // in 2d, edge can decay because of opened tilted cells
                    if(cell_node.template getEdgeNodeStatus<(TILT_DIRECTION_NEIGHBOR_MAC2[k][0])>()
                        && cell_node.template getEdgeNodeStatus<(TILT_DIRECTION_NEIGHBOR_MAC2[k][1])>()) {
                        return;
                    }
                }
                if(!cell_node.template getTiltNodeStatus<k>())  // closed tilt node
                    return;
                static_for<0, d>()([&](auto i){
                    face_area[kTiltDirectionNbMac<d, k>(i)] -= d==2? half_dx : 0.5*squared(half_dx);
                });
                IndexT nbCoor = coor + kOctNbTiltOffset<d, k>();
                int tIdx = pIndex->template get<node_type::vertex_node>(nbCoor);
                Real area;
                if constexpr(d == 2)
                    area = half_dx * SQRT2 * 
                        (1.0 - 0.5 * (cell_node.template getEdgeNodeStatus<(TILT_DIRECTION_NEIGHBOR_MAC2[k][0])>() ||
                        cell_node.template getEdgeNodeStatus<(TILT_DIRECTION_NEIGHBOR_MAC2[k][1])>()));
                if constexpr(d == 3)    // TODO 3d edge node
                    area = SQRT3 * 0.5 * squared(half_dx);
                REAL coeff = area / tilt_oct_dx;
                DEBUG_ONLY(if(tIdx < 0) { throw std::runtime_error("fuck"); })
                A.set_element(cIdx, tIdx, -coeff);
                center += coeff;
                int sign;
                if constexpr(d==2) sign = k<=tilt_direction::xp_yn?-1:1;
                else sign = k<=tilt_direction::xn_yp?-1:1;
                tb += sign * area * velocity->template get<node_type::vertex_node>(nbCoor)[oppositeTiltIdx<d>(k)];
            });
            // neighbor edge node or cell node
            static_for<0, d*2>()([&](auto k){
                face_area[k] += d == 2 ? dx : squared(dx);
                constexpr int axis = k / 2;
                bool is_neighbor_edge_open = false;
                if constexpr(d == 2) is_neighbor_edge_open = cell_node.template getEdgeNodeStatus<k>();
                // TODO 3d edge node
                if(is_neighbor_edge_open) {
                    IndexT nbCoor = coor; nbCoor[axis] += (k&1) ? 0 : 1;
                    int tIdx = pIndex->template get<node_type::edge_node_x+axis>(nbCoor);
                    DEBUG_ONLY(if(tIdx < 0) { throw std::runtime_error("fuck"); })
                    REAL coeff = half_dx / half_dx;
                    A.set_element(cIdx, tIdx, -coeff);
                    center += coeff;
                    tb += (k&1?1:-1) * half_dx * velocity->template get<node_type::edge_node_x+axis, axis>(nbCoor)[k^1];
                } else if(face_area[k] > EPSILON){
                    IndexT nbCoor = coor; nbCoor[axis] += (k&1) ? -1 : 1;
                    if(Valid(nbCoor, resolution)) {
                        int tIdx = pIndex->template get<node_type::cell_node>(nbCoor);
                        DEBUG_ONLY(if(tIdx < 0) { throw std::runtime_error("fuck"); })
                        REAL coeff = face_area[k] / dx;
                        A.set_element(cIdx, tIdx, -coeff);
                        center += coeff;
                        IndexT face_offset = IndexT::Unit(axis) * ((k&1) ? 0 : 1);
                        tb += (k&1?1:-1) * face_area[k] * velocity->template get<node_type::cell_node, axis>(coor+face_offset);
                    } else if constexpr(axis == 1) {
                        REAL coeff = face_area[k] / dx;
                        ASSERT(coeff == 1, "domain boundary shouldn't open tilt or edge nodes");
                        center += coeff;
                        IndexT face_offset = IndexT::Unit(axis) * ((k&1) ? 0 : 1);
                        tb += (k&1?1:-1) * face_area[k] * velocity->template get<node_type::cell_node, axis>(coor+face_offset);
                    }
                }
            });
            if(center < EPSILON) {
                spdlog::warn("zero diagonal");
                center = 1;
            }
            b[cIdx] = tb;
            A.set_element(cIdx, cIdx, center);
        });
        iterateASTAdaptiveDoF<ENABLE_PARALLEL>(resolution, [&](auto type, const IndexT& coor){
            REAL center = 0, tb = 0;
            int cIdx = pIndex->template get<type>(coor);
            if constexpr(type == node_type::vertex_node) {
                const auto& vertex_node = layout->getTiltNode(coor);
                if(!vertex_node.open) return;
                REAL face_area[d==2?4:8] = {0};
                // edge node neighbor
                if constexpr(d == 2) {
                    // TODO 3d edge node
                    static_for<0, d*2>()([&](auto k){
                        constexpr int axis = k/2;
                        if(vertex_node.template getEdgeNodeStatus<k>()) {
                            face_area[MAC_DIRECTION_NEIGHBOR_TILT2[k][0]] -= half_dx * SQRT2 * 0.5;
                            face_area[MAC_DIRECTION_NEIGHBOR_TILT2[k][1]] -= half_dx * SQRT2 * 0.5;
                            IndexT nbCoor = coor + IndexT::Unit(axis) * ((k&1) ? -1 : 0);
                            int tIdx = pIndex->template get<node_type::edge_node_x+(axis^1)>(nbCoor);
                            DEBUG_ONLY(if(tIdx < 0) { throw std::runtime_error("fuck"); })
                            REAL coeff = half_dx / half_dx;
                            A.set_element(cIdx, tIdx, -coeff);
                            center += coeff;
                            tb += (k&1?1:-1) * half_dx * velocity->template get<node_type::edge_node_x+(axis^1), axis^1>(nbCoor)[k^1];
                        }
                    });
                }
                // octagon node (cell node) neighbor
                const auto& u = velocity->template get<node_type::vertex_node>(coor);
                static_for<0, d==2?4:8>()([&](auto k) {
                    face_area[k] += d == 2 ? half_dx * SQRT2 : SQRT3 * 0.5 * squared(half_dx);
                    if(face_area[k] < EPSILON) return;
                    IndexT nbCoor = coor + kTiltNbOctOffset<d, k>();
                    int tIdx = pIndex->template get<node_type::cell_node>(nbCoor);
                    REAL coeff = face_area[k] / tilt_oct_dx;
                    DEBUG_ONLY(if(tIdx < 0) { throw std::runtime_error("fuck"); })
                    A.set_element(cIdx, tIdx, -coeff);
                    center += coeff;
                    int sign;
                    if constexpr(d==2) sign = (k <= 1 ? -1 : 1);
                    if constexpr(d==3) sign = (k <= tilt_direction::xn_yp ? -1 : 1);
                    tb += sign * face_area[k] * u[k];
                });
                A.set_element(cIdx, cIdx, center);
                b[cIdx] = tb;
            } else if constexpr(d==2 && is_edge_node<type>) {
                constexpr int axis = type - node_type::edge_node_x;
                if(!layout->template getEdgeNode<axis>(coor)) return;
                const auto& u = velocity->template get<node_type::edge_node_x+axis, axis>(coor);
                static_for<0, 2>()([&](auto k){ // cell neighbor
                    IndexT nbCoor = coor - IndexT::Unit(axis) * (k&1);
                    int tIdx = pIndex->template get<node_type::cell_node>(nbCoor);
                    REAL coeff = half_dx / half_dx;
                    A.set_element(cIdx, tIdx, -coeff);
                    center += coeff;
                    tb += (k==0?-1:1) * half_dx * u[axis*2+k];
                });
                static_for<0, 2>()([&](auto k){ // vertex neighbor, k:0 positive, k:1:negative
                    IndexT nbCoor = coor + IndexT::Unit(axis^1) * (1-k);
                    int tIdx = pIndex->template get<node_type::vertex_node>(nbCoor);
                    DEBUG_ONLY(if(tIdx < 0) { throw std::runtime_error("fuck"); })
                    REAL coeff = half_dx / half_dx;
                    A.set_element(cIdx, tIdx, -coeff);
                    center += coeff;
                    tb += (k==0?-1:1) * half_dx * u[(axis^1)*2+k];
                });
                A.set_element(cIdx, cIdx, center);
                b[cIdx] = tb;
            }
        });
        BLAS::scale(rho / dt, b);
        DEBUG_ONLY(isSymmetrical(A));
        // A.printReadable(std::cout);
        // solve
        {   BENCHMARK_SCOPED_TIMER_SECTION tSolve("pcg solve");
            REAL err; int iter;
#ifdef CUDA_SOLVER
            fixedA.construct_from_matrix(A);
            cudaSolver.solve(total_index_size, fixedA.rowstart, fixedA.colindex, fixedA.value, b, vRes);
            err = cudaSolver.solver_error; 
            iter = cudaSolver.solver_iters;
            spdlog::info("\t\t[project] MGPCG error:{}, iter:{}", err, iter);
#else
            pcgSolver.solve(A, b, vRes, err, iter);
            spdlog::info("\t\t[project] MICPCG error:{}, iter:{}", err, iter);
#endif
        }

        correctVelocityByPressure(dt);
        constrainVelocity();

        // DEBUG_ONLY(if(!isDivergenceFree()) {
        //     throw std::runtime_error("not divergence free");
        // });
    }

    void adjustLayout(REAL dt) {
        BENCHMARK_SCOPED_TIMER_SECTION tt("layout adjust");

        auto dst_layout = layout_tmp;
        iterateASTTilt<ENABLE_PARALLEL>(resolution, [&](const IndexT& coor){
            dst_layout->getTiltNode(coor).open = false;
            auto pos = layout->positionTilt(coor);
#ifdef LAYOUT_ADJUST_DENSITY
            auto oldPos = backTrace(velocity, pos, dt);
            REAL tmp = density->sample(oldPos);
            if(tmp > 0.005 && tmp < 0.7) {
                dst_layout->getTiltNode(coor).open = true;
            }
#endif
#ifdef LAYOUT_ADJUST_VORTEX
    #ifdef CORRECTION_MAC
            if(coor.x() == 0 || coor.y() == 0) return;
            REAL gy = velocity->template get<node_type::cell_node, 0>(coor)
                        - velocity->template get<node_type::cell_node, 0>(coor+IndexT(0, -1));
            REAL gx = velocity->template get<node_type::cell_node, 1>(coor)
                        - velocity->template get<node_type::cell_node, 1>(coor+IndexT(-1, 0));
            REAL vortex = (gx-gy)/dx;
    #else
            REAL vortex = computeVortexAt(velocity[data_flag], pos, dx);
    #endif
            if(fabs(vortex) > 10) {
                dst_layout->getTiltNode(coor).open = true;
            }
#endif
        });
        dst_layout->constrainFixedTiltNode();
#ifndef CORRECTION_MAC
        if constexpr(d==2) {
            dst_layout->clearEdgeNode();
            iterateMac<ENABLE_PARALLEL>(resolution, [&](auto axis, const IndexT& coor){
                if(dst_layout->getTiltNode(coor).open
                    && dst_layout->getTiltNode(coor+IndexT::Unit(axis^1)).open) {
                    auto pos = dst_layout->template positionEdgeNode<axis>(coor);
#ifdef LAYOUT_ADJUST_DENSITY
                    auto oldPos = backTrace(velocity, pos, dt); 
                    REAL tmp = density->sample(oldPos);
                    if(tmp > 0.04 && tmp < 0.5) {
                        dst_layout->template setEdgeNode<axis>(coor, 1);
                    }
#endif
#ifdef LAYOUT_ADJUST_VORTEX
                    REAL vortex = computeVortexAt(velocity, pos, dx);
                    if(fabs(vortex) > 15) {
                        dst_layout->template setEdgeNode<axis>(coor, 1);
                    }
#endif
                }
            });
        }
#endif
        dst_layout->propogateNodeInfo();
    }

    void advect(REAL dt) {
        BENCHMARK_SCOPED_TIMER_SECTION tt("advection");

        velocity->updateDualGrid();
        if(allow_adaptive_dofs) 
            adjustLayout(dt);

        auto dst = data_pool->template getGrid<grid_type::ast>(layout_tmp);
        semiLagrangianAST(velocity, dst, density, dt);  
        data_pool->template releaseGrid<grid_type::ast>(density); density = dst;
        
        dst = data_pool->template getGrid<grid_type::ast>(layout_tmp);
        semiLagrangianAST(velocity, dst, heat, dt);
        data_pool->template releaseGrid<grid_type::ast>(heat); heat = dst;

        auto dst_vel = data_pool->template getGrid<grid_type::ast_face>(layout_tmp);
        semiLagrangianASTVector(velocity, dst_vel, velocity, dt);
        data_pool->template releaseGrid<grid_type::ast_face>(velocity); velocity = dst_vel;

        layout.swap(layout_tmp);
    }

    // TODO: doesn't consider obstacles
    void updateIndex() {
        BENCHMARK_SCOPED_TIMER_SECTION tt("index update");
        pIndex->setLayout(layout);
        int cIdx = reuse_octagon_pIndex ? octagon_pIndex_max : 0;
        if(!reuse_octagon_pIndex) {   // re-assign index for ast cell centered DoFs
            iterateCellCentered(resolution, [&](auto type, const IndexT& coor){
                pIndex->template set<type>(coor, cIdx++);
            });
            octagon_pIndex_max = cIdx;
        }
        int vertex_dof_size = 0;
        iterateASTAdaptiveDoF(resolution, [&](auto type, const IndexT& coor){
            if(!isSkipClosed<type, d>(layout, coor)) {
                pIndex->template set<type>(coor, cIdx++);
                if constexpr(type == node_type::vertex_node) vertex_dof_size += 1;
            } else {
                pIndex->template set<type>(coor, -1);
            }
        });
        total_index_size = cIdx;
        spdlog::info("[update index] vertex DoFs: {}, edge Dofs:{}, fixed DoFs: {}", 
            vertex_dof_size, total_index_size - octagon_pIndex_max - vertex_dof_size, octagon_pIndex_max);
    }

    void constrainVelocity() {
        if constexpr(d == 2) {  // closed domain for x direction
            for(int y = 0; y < resolution.y(); ++y) {
                velocity->template get<node_type::cell_node, 0>(IndexT(0, y)) = 0;
                velocity->template get<node_type::cell_node, 0>(IndexT(resolution.x(), y)) = 0;
            }
        } else {    // closed domain for x-positive, x-negative, y-postive, y-negative
            for (int zIdx = 0; zIdx < resolution.z(); ++zIdx) {
                for (int yIdx = 0; yIdx < resolution.y(); ++yIdx) {
                    velocity->template get<node_type::cell_node, axis_type::x>(IndexT(0, yIdx, zIdx)) = 0;
                    velocity->template get<node_type::cell_node, axis_type::x>(IndexT(resolution.x(), yIdx, zIdx)) = 0;
                }
            }
            for (int zIdx = 0; zIdx < resolution.z(); ++zIdx) {
                for (int xIdx = 0; xIdx < resolution.y(); ++xIdx) {
                    velocity->template get<node_type::cell_node, axis_type::y>(IndexT(xIdx, 0, zIdx)) = 0;
                    velocity->template get<node_type::cell_node, axis_type::y>(IndexT(xIdx, resolution.y(), zIdx)) = 0;
                }
            }
        }
    }

    void saveVortex(int step, REAL min_v, REAL max_v) {
        if constexpr(d != 2) return;
        velocity->updateDualGrid();
        // REAL max_vortex = 0;
        Viewer::grid2PNGColorImage("./"+name+"/vortex/vortex_"+paddingStr(std::to_string(step), '0', 4)+".png", 
                resolution, [&](int xIdx, int yIdx) {
            Vector2f pos = Vector2f(xIdx, yIdx)  * dx;
            REAL v = computeVortexAt(velocity, pos, dx);
            // if(fabs(v) > max_vortex)
            //     max_vortex = fabs(v);
            return getHeatColor(v, min_v, max_v);
        });
        // std::cout<<max_vortex<<std::endl;
    }

    // basic
    const IndexT resolution;
    const Real dx, half_dx, quarter_dx;
    Real tilt_oct_dx;   // distance between the tilted cell and the octagon cell
    int max_cfl = 2;
    Real rho = 2.0;
    bool allow_adaptive_dofs = true;
    bool start_of_frame = true;

    // data fields: ast_layout, density, heat, velocity
    ASTGridPoolPtr<d> data_pool;
    ASTGridLayoutPtr<d> layout, layout_tmp;
    ASTGridPtr<REAL, d> density, heat;
    ASTFaceGridPtr<REAL, d> velocity;
    
    // pIndex related
    bool reuse_octagon_pIndex = false;
    int octagon_pIndex_max = 0;
    int total_index_size = 0;
    ASTGridPtr<int, d> pIndex;

    // solver related
    PCG::SparseMatrix<REAL> A;
    std::vector<REAL> b, vRes;
#ifdef CUDA_SOLVER
    FixedSparseMatrix<REAL> fixedA;
    CudaSolver cudaSolver;
#else
    PCGSolver<REAL> pcgSolver;
#endif

};


} // namespace ast

#pragma once
#include "ast_face_grid.h"
#include "ast_grid.h"
#include "ast_grid_layout.h"

namespace ast {

// single thread ast grid pool
template<int d>
class ASTGridPool {
    using IndexT = Eigen::Matrix<int, d, 1>;
public:
    ASTGridPool(const IndexT& res, REAL gridSpacing) 
        : resolution(res), dx(gridSpacing)
    {}

    template<int type>
    auto getGrid(const ASTGridLayoutPtr<d>& layout) {
        if constexpr(type == grid_type::ast) {
            if(ast_pool.size() == 0) {
                auto grid = std::make_shared<ASTGrid<REAL, d>>(resolution, dx, layout, 0);
                ast_pool.push_back(grid);
            }
            auto ret = ast_pool[ast_pool.size()-1];
            ret->setLayout(layout);
            ast_pool.pop_back();
            return ret;
        } else if constexpr(type == grid_type::ast_face) {
            if(ast_face_pool.size() == 0) {
                auto grid = std::make_shared<ASTFaceGrid<REAL, d>>(resolution, dx, layout);
                ast_face_pool.push_back(grid);
            }
            auto ret = ast_face_pool[ast_face_pool.size()-1];
            ret->setLayout(layout);
            ast_face_pool.pop_back();
            return ret;
        } else {
            static_assert(always_false<type>);
        }
    }
    template<int type, typename T>
    void releaseGrid(const T& grid) {
        if constexpr(type == grid_type::ast) 
            ast_pool.push_back(grid);
        else if constexpr(type == grid_type::ast_face)
            ast_face_pool.push_back(grid);
        else 
            static_assert(always_false<type>);
    }

private:
    IndexT resolution;
    REAL dx;
    std::vector<ASTGridPtr<Real, d>> ast_pool;
    std::vector<ASTFaceGridPtr<Real, d>> ast_face_pool;
};

template<int d> using ASTGridPoolPtr = std::shared_ptr<ASTGridPool<d>>;

}
#pragma once

#define AST_NAMESPACE_START namespace ast {
#define AST_NAMESPACE_END }

AST_NAMESPACE_START

#define ENABLE_PARALLEL true
#define DEBUG_VELOCITY false
#define LAYOUT_ADJUST_DENSITY
// #define LAYOUT_ADJUST_VORTEX
#define CORRECTION_MAC
#define EDGE_NODE_ADAPTIVITY false

template <class...> constexpr std::false_type always_false{};

enum axis_type {
    x = 0,
    y, z
};

enum grid_type {
    ast = 0,
    ast_face
};

AST_NAMESPACE_END
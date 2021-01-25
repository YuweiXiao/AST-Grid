#include "general.h"
#include "constants.h"
#include "ast_general.h"

AST_NAMESPACE_START

using Omni::Vector2f;
using Omni::Vector3f;
using Omni::Vector4f;
using Omni::Size3;

enum tilt_direction {
    xp_yp = 0, // short for x positive, y positive
    xp_yn, xn_yn, xn_yp, // 2D case or default z positive
    xp_yp_zn, xp_yn_zn, xn_yn_zn, xn_yp_zn,
};

enum mac_direction {
    x_positive = 0, x_negative,
    y_positive, y_negative,
    z_positive, z_negative,
};

static const Vector2f TILT_POS_OFFSET[4] = {
    Vector2f(INV_SQRT2, INV_SQRT2), Vector2f(INV_SQRT2, -INV_SQRT2), Vector2f(-INV_SQRT2, -INV_SQRT2), Vector2f(-INV_SQRT2, INV_SQRT2)
};
static const Vector3f TILT_POS_OFFSET3[] = {
    Vector3f::Ones()/3.0, Vector3f(1.0/3, -1.0/3, 1.0/3), Vector3f(-1.0/3, -1.0/3, 1.0/3), Vector3f(-1.0/3, 1.0/3, 1.0/3),
    Vector3f(1.0/3, 1.0/3, -1.0/3), Vector3f(1.0/3, -1.0/3, -1.0/3), Vector3f(-1.0/3, -1.0/3, -1.0/3), Vector3f(-1.0/3, 1.0/3, -1.0/3)
};

static const Vector3f TILT_UNIT_DIR3[8] = {
    Vector3f(INV_SQRT3, INV_SQRT3, INV_SQRT3), Vector3f(INV_SQRT3, -INV_SQRT3, INV_SQRT3), 
    Vector3f(-INV_SQRT3, -INV_SQRT3, INV_SQRT3), Vector3f(-INV_SQRT3, INV_SQRT3, INV_SQRT3),
    Vector3f(INV_SQRT3, INV_SQRT3, -INV_SQRT3), Vector3f(INV_SQRT3, -INV_SQRT3, -INV_SQRT3), 
    Vector3f(-INV_SQRT3, -INV_SQRT3, -INV_SQRT3), Vector3f(-INV_SQRT3, INV_SQRT3, -INV_SQRT3)
};

constexpr int TILT_DIRECTION_NEIGHBOR_MAC2[][2] = {
    {mac_direction::x_positive, mac_direction::y_positive},  // xp_yp
    {mac_direction::x_positive, mac_direction::y_negative},  // xp_yn
    {mac_direction::x_negative, mac_direction::y_negative},  // xn_yn
    {mac_direction::x_negative, mac_direction::y_positive}   //xn_yp
};

constexpr int TILT_DIRECTION_NEIGHBOR_MAC3[][3] = {
    {mac_direction::x_positive, mac_direction::y_positive, mac_direction::z_positive},  // xp_yp
    {mac_direction::x_positive, mac_direction::y_negative, mac_direction::z_positive},  // xp_yn
    {mac_direction::x_negative, mac_direction::y_negative, mac_direction::z_positive},  // xn_yn
    {mac_direction::x_negative, mac_direction::y_positive, mac_direction::z_positive},  // xn_yp
    {mac_direction::x_positive, mac_direction::y_positive, mac_direction::z_negative},  // xp_yp_zn
    {mac_direction::x_positive, mac_direction::y_negative, mac_direction::z_negative},  // xp_yn_zn
    {mac_direction::x_negative, mac_direction::y_negative, mac_direction::z_negative},  // xn_yn_zn
    {mac_direction::x_negative, mac_direction::y_positive, mac_direction::z_negative}   // xn_yp_zn
};

const Size3 OCT_NB_TILT_DXDYDZ[8] = {
    Size3(1, 1, 1), Size3(1, 0, 1), Size3(0, 0, 1), Size3(0, 1, 1), 
    Size3(1, 1, 0), Size3(1, 0, 0), Size3(0, 0, 0), Size3(0, 1, 0)
};

const Size3 TILT_NB_OCT_DXDYDZ[8] = {
    Size3(0, 0, 0), Size3(0, -1, 0), Size3(-1, -1, 0), Size3(-1, 0, 0),
    Size3(0, 0, -1), Size3(0, -1, -1), Size3(-1, -1, -1), Size3(-1, 0, -1)
};

template<int d, int k>
constexpr int kTiltDirectionNbMac(int i) {
    if constexpr(d == 2) return TILT_DIRECTION_NEIGHBOR_MAC2[k][i];
    else return TILT_DIRECTION_NEIGHBOR_MAC3[k][i];
}

template<int d, int k>
constexpr Eigen::Matrix<Real, d, 1> kTiltPosOffset() {
    if constexpr(d == 2) {  return TILT_POS_OFFSET[k]; }
    else return TILT_POS_OFFSET3[k];
}

template<int d, int k>
constexpr auto kTiltNbOctOffset() {
    if constexpr(d == 2) return Omni::TILT_NB_OCT_DXDY[k];
    else return TILT_NB_OCT_DXDYDZ[k];
}

template<int d, int k>
constexpr auto kOctNbTiltOffset() {
    if constexpr(d == 2) return Omni::OCT_NB_TILT_DXDY[k];
    else return OCT_NB_TILT_DXDYDZ[k];
}

template<int d, int k>
constexpr auto kTiltUnitDir() {
    if constexpr(d == 2) return Omni::TILT_UNIT_DIR[k];
    else return TILT_UNIT_DIR3[k];
}

template<int d, int k>
constexpr auto kTiltUnitUpDir() {
    if constexpr(d == 2) return Omni::TILT_UNIT_DIR[k%2];
    else return TILT_UNIT_DIR3[k%4];
}


// TODO: make it template and workable for both 2d & 3d
constexpr int MAC_DIRECTION_NEIGHBOR_TILT2[][2] = {
    {tilt_direction::xp_yp, tilt_direction::xp_yn},  // x+
    {tilt_direction::xn_yp, tilt_direction::xn_yn},  // x-
    {tilt_direction::xp_yp, tilt_direction::xn_yp},  // y+
    {tilt_direction::xp_yn, tilt_direction::xn_yn}   // y-
};


AST_NAMESPACE_END
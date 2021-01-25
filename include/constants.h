#pragma once
#include "general.h"

namespace Omni {


//! No direction.
constexpr int kDirectionNone = 0;

//! Left direction.
constexpr int kDirectionXNegative = 1 << 0;

//! RIght direction.
constexpr int kDirectionXPositive = 1 << 1;

//! Down direction.
constexpr int kDirectionZNegative = 1 << 2;

//! Up direction.
constexpr int kDirectionZPositive = 1 << 3;

//! Back direction.
constexpr int kDirectionYNegative = 1 << 4;

//! Front direction.
constexpr int kDirectionYPositive = 1 << 5;

//! All direction.
constexpr int kDirectionAll = kDirectionXNegative | kDirectionXPositive
                            | kDirectionYNegative | kDirectionYPositive
                            | kDirectionZNegative | kDirectionZPositive;

const Size2 UNIT_VECTOR2[2] = {
    Size2(1, 0), Size2(0, 1)
};

const Vector2f TILT_UNIT_DIR[4] = {
    Vector2f(INV_SQRT2, INV_SQRT2), Vector2f(INV_SQRT2, -INV_SQRT2), 
    Vector2f(-INV_SQRT2, -INV_SQRT2), Vector2f(-INV_SQRT2, INV_SQRT2)
};

const Vector3f TILT3_UNIT_DIR[8] = {
    Vector3f(INV_SQRT3, INV_SQRT3, INV_SQRT3), Vector3f(-INV_SQRT3, INV_SQRT3, INV_SQRT3), 
    Vector3f(-INV_SQRT3, -INV_SQRT3, INV_SQRT3), Vector3f(INV_SQRT3, -INV_SQRT3, INV_SQRT3),
    Vector3f(INV_SQRT3, INV_SQRT3, -INV_SQRT3), Vector3f(-INV_SQRT3, INV_SQRT3, -INV_SQRT3), 
    Vector3f(-INV_SQRT3, -INV_SQRT3, -INV_SQRT3), Vector3f(INV_SQRT3, -INV_SQRT3, -INV_SQRT3)
};

const Vector8f TILT3_DIR_UP_Z = [] {
    Vector8f tmp;
    tmp << INV_SQRT3, INV_SQRT3, INV_SQRT3, INV_SQRT3, INV_SQRT3, INV_SQRT3, INV_SQRT3, INV_SQRT3;
    return tmp;
}();

const Vector3f TILT3_UNIT_DIR_UP[8] = {
    Vector3f(INV_SQRT3, INV_SQRT3, INV_SQRT3), Vector3f(-INV_SQRT3, INV_SQRT3, INV_SQRT3), 
    Vector3f(-INV_SQRT3, -INV_SQRT3, INV_SQRT3), Vector3f(INV_SQRT3, -INV_SQRT3, INV_SQRT3),
    Vector3f(-INV_SQRT3, -INV_SQRT3, INV_SQRT3), Vector3f(INV_SQRT3, -INV_SQRT3, INV_SQRT3), 
    Vector3f(INV_SQRT3, INV_SQRT3, INV_SQRT3), Vector3f(-INV_SQRT3, INV_SQRT3, INV_SQRT3)
};

// right, down, left, up
const Size2 OCT_NB_DXDY[4] = {
    Size2(1, 0), Size2(0, -1), Size2(-1, 0), Size2(0, 1)
};
const Size2 OCT_NB_FACE_DXDY[4] = {
    Size2(1, 0), Size2(0, 0), Size2(0, 0), Size2(0, 1)
};

// up-right, lower-right, lower-left, up-left
const Size2 OCT_NB_TILT_DXDY[] = {
    Size2(1, 1), Size2(1, 0), Size2(0, 0), Size2(0, 1), Size2(1, 1)
};

const Size2 TILT_NB_OCT_DXDY[4] = {
    Size2(0, 0), Size2(0, -1), Size2(-1, -1), Size2(-1, 0)
};

const Size2 CELL_NB_DXDY[8] = {
    Size2(1, 0), Size2(1, -1), Size2(0, -1), Size2(-1, -1),
    Size2(-1, 0), Size2(-1, 1), Size2(0, 1), Size2(1, 1)
};

// +x, -x, +y, -y, +z, -z
const Size3 OCT_NB_DXDYDZ[6] = {
    Size3(1, 0, 0), Size3(-1, 0, 0), Size3(0, 1, 0), Size3(0, -1, 0), 
    Size3(0, 0, 1), Size3(0, 0, -1)
};

const Size3 OCT_NB_TILT_DXDYDZ[8] = {
    Size3(1, 1, 1), Size3(0, 1, 1), Size3(0, 0, 1), Size3(1, 0, 1), 
    Size3(1, 1, 0), Size3(0, 1, 0), Size3(0, 0, 0), Size3(1, 0, 0)
};

// follow the order of OCT_NB_DXDYDZ
const Size3 OCT_NB_FACE_DXDYDZ[6] = {
    Size3(1, 0, 0), Size3(0, 0, 0), Size3(0, 1, 0), Size3(0, 0, 0),
    Size3(0, 0, 1), Size3(0, 0, 0)
};

// follow the order of Vector8f
const Size3 TILT_NB_OCT_DXDYDZ[8] = {
    Size3(0, 0, 0), Size3(-1, 0, 0), Size3(-1, -1, 0), Size3(0, -1, 0),
    Size3(0, 0, -1), Size3(-1, 0, -1), Size3(-1, -1, -1), Size3(0, -1, -1)
};

// For each direction:
//  - first 4 indices means the index of valid tilt direction in Vector8 vector
//  - the 5th is index of XYZ surface of half-tilt in Vector8f vector
//  - the 5th index also used to get surface direction combined with HALF_TILT_XYZ_DIR array
const int HALF_TILT_SLOT_INDEX[][5]{
    {0,3,4,7,1}, // x
    {0,1,4,5,2}, // y
    {0,1,2,3,5}, // z
    {1,2,5,6,0}, // -x
    {2,3,6,7,4}, // -y
    {4,5,6,7,3}  // -z
};

// -x, x, y, -z, -y, z
const Vector3f HALF_TILT_XYZ_DIR[6] {
    Vector3f(-1, 0, 0), Vector3f(1, 0, 0), Vector3f(0, 1, 0),
    Vector3f(0, 0, -1), Vector3f(0, -1, 0), Vector3f(0, 0, 1)
};

// -x, x, y, -z, -y, z
const Vector3f HALF_TILT_XYZ_DIR_POSITIVE[6] {
    Vector3f(1, 0, 0), Vector3f(1, 0, 0), Vector3f(0, 1, 0),
    Vector3f(0, 0, 1), Vector3f(0, 1, 0), Vector3f(0, 0, 1)
};

enum XYZ_MARKER {
	NONE,
	X = 0x100, _X = 0x200,   // x, -x
	Y = 0x010, _Y = 0x020,   // y, -y
	Z = 0x001, _Z = 0x002,   // z, -z
    XY = 0x110, XZ = 0x101, YZ = 0x011, XYZ = 0x111,
    _X_Y = 0x220, _X_Z = 0x202, _Y_Z = 0x022, _X_Y_Z = 0x222, 
    _XY = 0x210, _XZ = 0x201, _YZ = 0x021, _XYZ = 0x211,
    X_Y = 0x120, X_YZ = 0x121,
    X_Z = 0x102, Y_Z = 0x012, XY_Z = 0x112,
    _X_YZ = 0x221, X_Y_Z = 0x122, _XY_Z = 0x212,
    VOIDAXIS = 0x1000
};

const XYZ_MARKER OCTREE_SAMPLE_MARKER[8][8] {
    {XYZ_MARKER::VOIDAXIS, XYZ_MARKER::YZ, XYZ_MARKER::Z, XYZ_MARKER::XZ, XYZ_MARKER::XY, XYZ_MARKER::Y, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::X},
    {XYZ_MARKER::YZ, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_XZ, XYZ_MARKER::Z, XYZ_MARKER::Y, XYZ_MARKER::_XY, XYZ_MARKER::_X, XYZ_MARKER::VOIDAXIS},
    {XYZ_MARKER::Z, XYZ_MARKER::_XZ, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_YZ, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_X, XYZ_MARKER::_X_Y, XYZ_MARKER::_Y},
    {XYZ_MARKER::XZ, XYZ_MARKER::Z, XYZ_MARKER::_YZ, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::X, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_Y, XYZ_MARKER::X_Y},
    {XYZ_MARKER::XY, XYZ_MARKER::Y, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::X, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::Y_Z, XYZ_MARKER::_Z, XYZ_MARKER::X_Z},
    {XYZ_MARKER::Y, XYZ_MARKER::_XY, XYZ_MARKER::_X, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::Y_Z, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_X_Z, XYZ_MARKER::_Z},
    {XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_X, XYZ_MARKER::_X_Y, XYZ_MARKER::_Y, XYZ_MARKER::_Z, XYZ_MARKER::_X_Z, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_Y_Z},
    {XYZ_MARKER::X, XYZ_MARKER::VOIDAXIS, XYZ_MARKER::_Y, XYZ_MARKER::X_Y, XYZ_MARKER::X_Z, XYZ_MARKER::_Z, XYZ_MARKER::_Y_Z, XYZ_MARKER::VOIDAXIS}
};

enum INTEGRATION_METHOD {
    EXPLICIT,
    EXPLICIT_IMPLICIT
};

class Constants {
public:
    static void init();

    static REAL OMNI_DUAL_GRAPH_MATRIX2[2][2][2][2][2][4];
    static REAL OMNI_DUAL_GRAPH_MATRIX3[2][2][2][2][3][7];
    // static REAL OCTREE_OMNI_DUAL_GRAPH_MATRIX3[2][2][2][2][2][2][2][3][7];
    static Matrix37f OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[2][2][2][2][2][2][2];
    static Matrix37f OMNI_DUAL_GRAPH_EIGEN_MATRIX3[2][2][2][2];
    static Matrix34d OMNI_DUAL_GRAPH_TILT_EIGEN_MATRIX3;
    static Eigen::Matrix2d ROTATE45;
    static Eigen::Matrix2d TILT_PRINCIPLE_AXIS2;
    static Eigen::Matrix3d TILT_PRINCIPLE_AXIS[8];
};


}
#pragma once
#include "general.h"

namespace Omni {

// simplified
struct TiltENodeS {
    // static int constexpr kdIsOpenIdx = 0;
    // static int constexpr kdIsFullSizeIdx = 1; // is maximized to equal to dx, if not, the size is 0.5dx when it is
    // open static int constexpr kdIsFixedIdx = 2; static int constexpr kdIsTJointIdx = 3; static int constexpr
    // kdIsClosedByTJointIdx = 4; static int constexpr kdIsHalfTiltIdx = 5;             // NOTE: only for 3D octree
    // static int constexpr kdIsHalfTiltDirectionIdx = 6;    // x, y, z, -x, -y, -z; negative means it is tilt closed by
    // half-tilt static int constexpr kdIsHalfTiltDirectionMask = 0x7; // binary: 111

    // int flag;
    // bool isClosed() const { return (flag >> TiltENodeS::kdIsOpenIdx) & 1; }
    // bool isFullSize() const { return (flag >> TiltENodeS::kdIsFullSizeIdx) & 1; }
    // bool isFixed() const { return (flag >> TiltENodeS::kdIsFixedIdx) & 1; }
    // bool isTJoint() const { return (flag >> TiltENodeS::kdIsTJointIdx) & 1; }
    // bool isClosedByTJoint() const { return (flag >> TiltENodeS::kdIsClosedByTJointIdx) & 1; }
    // bool isHalfTilt() const { return (flag >> TiltENodeS::kdIsHalfTiltIdx) & 1; }
    // char getHalfTiltDirection() const { return (flag >> TiltENodeS::kdIsHalfTiltDirectionIdx) &
    // TiltENodeS::kdIsHalfTiltDirectionMask; };

    TiltENodeS() : is_open(false), is_fixed(false) {}
    // TiltENodeS()
    //     : is_open(false),
    //       is_full_size(false),
    //       is_fixed(false),
    //       is_TJoint(false),
    //       is_closed_by_TJoint(false),
    //       is_half_tilt(false),
    //       half_tilt_direction(0) {}

    bool is_open : 1;
    // bool is_full_size : 1;
    bool is_fixed : 1;
    // bool is_TJoint : 1;
    // bool is_closed_by_TJoint : 1;
    // bool is_half_tilt : 1;
    // int half_tilt_direction : 4;

    // Real getSize(Real dx) const {
    //     Real e = is_full_size ? dx : 0.5 * dx;
    //     return (is_open || is_half_tilt) ? e : 0;
    // }
};

struct TiltENode {
    REAL e;
    bool is_largest;
    bool is_closed;
    bool is_fixed;
    bool is_T_joint;
    bool is_closed_by_T_joint;
    bool is_half_tilt;   // NOTE: only for 3D octree
    int half_direction;  // x, y, z, -x, -y, -z; negative means it is tilt closed by half-tilt

    TiltENode(REAL _e = 0, bool _is_largest = false, bool _is_closed = true, bool _is_fixed = false)
        : e(_e), is_largest(_is_largest), is_closed(_is_closed) {
        if (is_closed) e = 0;
        setE(e);

        is_largest = _is_largest;
        is_closed = _is_closed;
        is_fixed = _is_fixed;
        is_T_joint = false;
        is_closed_by_T_joint = false;
        is_half_tilt = false;
        half_direction = 0;
    }

    void setIsFixed(bool p) { is_fixed = p; }

    void setE(REAL _e, REAL maxE = -1) {
        if (is_fixed) return;
        if (_e <= 0) {
            is_closed = true;
            is_largest = false;
            _e = 0;
        } else if (maxE > 0 && _e >= maxE) {
            is_closed = false;
            is_largest = true;
            _e = maxE;
        } else {
            is_closed = false;
            is_largest = false;
        }
        e = _e;
    }

    TiltENode operator+(const TiltENode& rhs) const {
        REAL ne = this->e + rhs.e;
        return TiltENode(ne);
    }

    TiltENode operator*(REAL m) const {
        REAL ne = this->e * m;
        return TiltENode(ne);
    }
};

}  // namespace Omni
#pragma once
#include "general.h"

namespace Omni {

// simplified
struct TiltENodeS {

    TiltENodeS() : is_open(false), is_fixed(false) {}

    bool is_open : 1;
    bool is_fixed : 1;
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
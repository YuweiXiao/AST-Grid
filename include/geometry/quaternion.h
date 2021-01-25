#pragma once
#include "general.h"


namespace Omni {
namespace Geometry {

// Quaternion class defined as q = w + xi + yj + zk.
class Quaternion {
 public:
    FLOAT w;
    FLOAT x;
    FLOAT y;
    FLOAT z;

    Quaternion() {setIdentity();}
    Quaternion(FLOAT _w, FLOAT _x, FLOAT _y, FLOAT _z) { set(_w, _x, _y, _z); }
    //! Constructs a quaternion with given rotation axis and angle.
    Quaternion(const Vector3f& axis, FLOAT angle) {set(axis, angle);}
    Quaternion(const Quaternion& other) {set(other);}


    void set(const Quaternion& other) {set(other.w, other.x, other.y, other.z);}
    
    void set(FLOAT _w, FLOAT _x, FLOAT _y, FLOAT _z) {
        w = _w;
        x = _x;
        y = _y;
        z = _z;
    }
    
    void set(const Vector3f& axis, FLOAT angle) {
        if(axis.squaredNorm() < EPSILON) {
            setIdentity();
        } else {
            FLOAT s = std::sin(angle / 2.0);
            x = axis.normalized().x() * s;
            y = axis.normalized().y() * s;
            z = axis.normalized().z() * s;
            w = std::cos(angle / 2.0);
        }
    }

    // MARK: Basic getters
    Quaternion normalized() const {
        Quaternion q(*this);
        q.normalize();
        return q;
    }

    // MARK: Binary operator methods - new instance = this instance (+) input
    Vector3f mul(const Vector3f& v) const {
        FLOAT _2xx = 2 * x * x;
        FLOAT _2yy = 2 * y * y;
        FLOAT _2zz = 2 * z * z;
        FLOAT _2xy = 2 * x * y;
        FLOAT _2xz = 2 * x * z;
        FLOAT _2xw = 2 * x * w;
        FLOAT _2yz = 2 * y * z;
        FLOAT _2yw = 2 * y * w;
        FLOAT _2zw = 2 * z * w;

        return Vector3f(
            (1 - _2yy - _2zz)*v.x() + (_2xy - _2zw)*v.y() + (_2xz + _2yw)*v.z(),
            (_2xy + _2zw)*v.x() + (1 - _2zz - _2xx)*v.y() + (_2yz - _2xw)*v.z(),
            (_2xz - _2yw)*v.x() + (_2yz + _2xw)*v.y() + (1 - _2yy - _2xx)*v.z());
    }
    Quaternion mul(const Quaternion& other) const {
        return Quaternion(
                w * other.w - x * other.x - y * other.y - z * other.z,
                w * other.x + x * other.w + y * other.z - z * other.y,
                w * other.y - x * other.z + y * other.w + z * other.x,
                w * other.z + x * other.y - y * other.x + z * other.w);
    }
    FLOAT dot(const Quaternion& other) {
        return w * other.w + x * other.x + y * other.y + z * other.z;
    }

    // MARK: Binary operator methods - new instance = input (+) this instance
    Quaternion rmul(const Quaternion& other) const {
        return Quaternion(
                other.w * w - other.x * x - other.y * y - other.z * z,
                other.w * x + other.x * w + other.y * z - other.z * y,
                other.w * y - other.x * z + other.y * w + other.z * x,
                other.w * z + other.x * y - other.y * x + other.z * w);
    }

    // MARK: Augmented operator methods - this instance (+)= input
    void imul(const Quaternion& other) {
        *this = mul(other);
    }

    // MARK: Modifiers
    void setIdentity() {set(1, 0, 0, 0);}
    void rotate(FLOAT angleInRadians) {
        Vector3f axis;
        FLOAT currentAngle;
        getAxisAngle(axis, currentAngle);
        currentAngle += angleInRadians;
        set(axis, currentAngle);
    }
    void normalize() {
        FLOAT norm = l2Norm();
        if(norm > 0) {
            w /= norm;
            x /= norm;
            y /= norm;
            z /= norm;
        }
    }

    // MARK: Complex getters
    Vector3f axis() const {
        Vector3f result(x, y, z);
        result.normalize();

        if (2 * std::acos(w) < PI) {
            return result;
        } else {
            return -result;
        }
    }
    FLOAT angle() const {
        FLOAT result = 2 * std::acos(w);

        if (result < PI) {
            return result;
        } else {
            // Wrap around
            return 2 * PI - result;
        }
    }
    void getAxisAngle(Vector3f& axis, FLOAT& angle) const {
        axis.x() = x; 
        axis.y() = y;
        axis.z() = z;
        axis.normalize();
        angle = 2 * std::acos(w);

        if (angle > PI) {
            // Wrap around
            axis = -axis;
            angle = 2 * PI - angle;
        }
    }
    Quaternion inverse() const {
        FLOAT denom = w * w + x * x + y * y + z * z;
        return Quaternion(w / denom, -x / denom, -y / denom, -z / denom);
    }

    Matrix3f matrix3() const {
        FLOAT _2xx = 2 * x * x;
        FLOAT _2yy = 2 * y * y;
        FLOAT _2zz = 2 * z * z;
        FLOAT _2xy = 2 * x * y;
        FLOAT _2xz = 2 * x * z;
        FLOAT _2xw = 2 * x * w;
        FLOAT _2yz = 2 * y * z;
        FLOAT _2yw = 2 * y * w;
        FLOAT _2zw = 2 * z * w;

        Matrix3f m;
        m << 1 - _2yy - _2zz, _2xy - _2zw, _2xz + _2yw,
            _2xy + _2zw, 1 - _2zz - _2xx, _2yz - _2xw,
            _2xz - _2yw, _2yz + _2xw, 1 - _2yy - _2xx;

        return m;
    }
    Matrix4f matrix4() const {
        FLOAT _2xx = 2 * x * x;
        FLOAT _2yy = 2 * y * y;
        FLOAT _2zz = 2 * z * z;
        FLOAT _2xy = 2 * x * y;
        FLOAT _2xz = 2 * x * z;
        FLOAT _2xw = 2 * x * w;
        FLOAT _2yz = 2 * y * z;
        FLOAT _2yw = 2 * y * w;
        FLOAT _2zw = 2 * z * w;

        Matrix4f m;
        m << 1 - _2yy - _2zz, _2xy - _2zw, _2xz + _2yw, 0,
            _2xy + _2zw, 1 - _2zz - _2xx, _2yz - _2xw, 0,
            _2xz - _2yw, _2yz + _2xw, 1 - _2yy - _2xx, 0,
            0, 0, 0, 1;
        return m;
    }

    //! Returns L2 norm of this quaternion.
    FLOAT l2Norm() const {return std::sqrt(w * w + x * x + y * y + z * z);}
    
    bool operator==(const Quaternion& other) const {
        return w == other.w && x == other.x && y == other.y && z == other.z;
    }
    bool operator!=(const Quaternion& other) const {
        return !((*this) == other);
    }
};

// ! Computes spherical linear interpolation.
// template <typename T>
// Quaternion slerp(const Quaternion& a, const Quaternion& b, T t);
Vector3f operator*(const Quaternion& q, const Vector3f& v);
Quaternion operator*(const Quaternion& a, const Quaternion& b);

} // end of namespace Geometry
} // end of namespace Omni

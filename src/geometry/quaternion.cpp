#include "geometry/quaternion.h"

namespace Omni {
namespace Geometry {

Vector3f operator*(const Quaternion& q, const Vector3f& v) { 
    return q.mul(v); 
}

Quaternion operator*(const Quaternion& a, const Quaternion& b) { 
    return a.mul(b); 
}


} // end of namespace Geometry
} // end of namespace Omni
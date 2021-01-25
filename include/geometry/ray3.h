#pragma once
#include "general.h"

namespace Omni {
namespace Geometry {

// Class for 3-D ray.
class Ray3 {
 public:
    Vector3f origin;
    Vector3f direction;

    //! Constructs an empty ray that points (1, 0, 0) from (0, 0, 0).
    Ray3() 
        : origin(1, 0, 0), direction(0, 0, 0)
    {}

    Ray3(const Vector3f& newOrigin, const Vector3f& newDirection) 
        : origin(newOrigin), direction(newDirection)
    {}

    Ray3(const Ray3& other) 
        :origin(other.origin), direction(other.direction)
    {}

    Vector3f pointAt(FLOAT t) const {
        return origin + direction * t;
    }
};

} // end of namespace geometry
} // end of namespace Omni


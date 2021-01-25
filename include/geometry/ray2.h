#pragma once
#include "general.h"

namespace Omni {
namespace Geometry {

class Ray2 {
public:
    Vector2f origin;

    //! The normalized direction of the ray.
    Vector2f direction;

    //! Constructs an empty ray that points (1, 0) from (0, 0).
    Ray2() 
        :origin(0, 0), direction(1, 0)
    {}

    //! Constructs a ray with given origin and direction.
    Ray2(const Vector2f& _origin, const Vector2f& _dir)
        :origin(_origin), direction(_dir)
    {}

    Ray2(const Ray2& other)
        :origin(other.origin), direction(other.direction)
    {}

    Vector2f pointAt(FLOAT t) const {
        return origin + t * direction;
    }
};

}  // end of namespace Geometry
}  // end of namespace Omni
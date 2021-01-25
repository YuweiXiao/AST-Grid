#include "geometry/sphere2.h"

namespace Omni {
namespace Geometry {

Sphere2::Sphere2(const Vector2f& c, FLOAT r, const Transform2& transform, bool isNormalFlipped)
	: Shape2(transform, isNormalFlipped), center(c), radius(r)
{
	_squaredR = r * r;
}


Vector2f Sphere2::closestPointLocal(const Vector2f& otherPoint) const {
    return radius * closestNormalLocal(otherPoint) + center;
}

FLOAT Sphere2::closestDistanceLocal(const Vector2f& otherPoint) const {
    return std::fabs((center - otherPoint).norm() - radius);
}

Vector2f Sphere2::closestNormalLocal(const Vector2f& otherPoint) const {
    if (center == otherPoint) {
        return Vector2f(1, 0);
    } else { 
        return (otherPoint - center).normalized();
    }
}

bool Sphere2::intersectsLocal(const Ray2& ray) const {
    Vector2f r = ray.origin - center;
    FLOAT b = ray.direction.dot(r);
    FLOAT c = r.squaredNorm() - _squaredR;
    FLOAT d = b * b - c;

    if (d > 0.0) {
        d = std::sqrt(d);
        FLOAT tMin = -b - d;
        FLOAT tMax = -b + d;

        if (tMin < 0.0) {
            tMin = tMax;
        }

        if (tMin >= 0.0) {
            return true;
        }
    }

    return false;
}

ShapeRayIntersection2 Sphere2::closestIntersectionLocal(const Ray2& ray) const {
    ShapeRayIntersection2 intersection;
    Vector2f r = ray.origin - center;
    FLOAT b = ray.direction.dot(r);
    FLOAT c = r.squaredNorm() - _squaredR;
    FLOAT d = b * b - c;

    if (d > 0.0) {
        d = std::sqrt(d);
        FLOAT tMin = -b - d;
        FLOAT tMax = -b + d;

        if (tMin < 0.0) {
            tMin = tMax;
        }

        if (tMin < 0.0) {
            intersection.isIntersecting = false;
        } else {
            intersection.isIntersecting = true;
            intersection.distance = tMin;
            intersection.point = ray.origin + tMin * ray.direction;
            intersection.normal = (intersection.point - center).normalized();
        }
    } else {
        intersection.isIntersecting = false;
    }

    return intersection;
}

} // end of namespace Geometry
} // end of namespace Omni
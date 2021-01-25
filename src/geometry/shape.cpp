#include "geometry/shape.h"

using namespace Omni;
using namespace Omni::Geometry;

bool Shape2::intersects(const Ray2& ray) const {
    return intersectsLocal(transform.toLocal(ray));
}


ShapeRayIntersection2 Shape2::closestIntersection(const Ray2& ray) const {
    auto result = closestIntersectionLocal(transform.toLocal(ray));
    result.point = transform.toWorld(result.point);
    result.normal = transform.toWorldDirection(result.normal);
    result.normal *= (isNormalFlipped) ? -1.0 : 1.0;
    return result;
}

bool Shape2::intersectsLocal(const Ray2& rayLocal) const {
    auto result = closestIntersectionLocal(rayLocal);
    return result.isIntersecting;
}

bool Shape3::intersects(const Ray3& ray) const {
    return intersectsLocal(transform.toLocal(ray));
}

ShapeRayIntersection3 Shape3::closestIntersection(const Ray3& ray) const {
    auto result = closestIntersectionLocal(transform.toLocal(ray));
    result.point = transform.toWorld(result.point);
    result.normal = transform.toWorldDirection(result.normal);
    result.normal *= (isNormalFlipped) ? -1.0 : 1.0;
    return result;
}

bool Shape3::intersectsLocal(const Ray3& ray) const {
    auto result = closestIntersectionLocal(ray);
    return result.isIntersecting;
}

BoundingBox3 Shape3::boundingBox() const {
	return transform.toWorld(boundingBoxLocal());
}

BoundingBox3 Omni::Geometry::Shape3::boundingBoxLocal() const
{
	return BoundingBox3();
}
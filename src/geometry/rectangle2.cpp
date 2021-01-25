#include "geometry/rectangle2.h"

using namespace Omni;
using namespace Omni::Geometry;


Rectangle2::Rectangle2(const Vector2f& _wh, Transform2 transform, bool isNormalFlipped) 
    :Shape2(transform, isNormalFlipped), wh(_wh), halfWH(_wh.x()/2.0, _wh.y()/2.0)
{}

bool Rectangle2::intersectsLocal(const Ray2& ray) const {
    // TODO
    throw std::runtime_error("[Rectangle2::intersectsLocal] not implemented yet");
}

ShapeRayIntersection2 Rectangle2::closestIntersectionLocal(const Ray2& ray) const {
    // TODO
    throw std::runtime_error("[Rectangle2::closestIntersectionLocal] not implemented yet");
}
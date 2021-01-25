#include "geometry/implicit_surface2set.h"

using namespace Omni;
using namespace Omni::Geometry;

Vector2f ImplicitSurface2Set::closestPointLocal(const Vector2f& otherPoint) const {
    ImplicitSurface2Ptr s = nullptr;
    FLOAT min_dis = MAX_FLOAT;
    for(int i = 0; i < _surfaces.size(); ++i) {
        FLOAT dis = _surfaces[i]->closestDistance(otherPoint);
        if(min_dis > dis) {
            min_dis = dis;
            s = _surfaces[i];
        }
    }
    if(s != nullptr) {
        return s->closestPoint(otherPoint);
    } else {
        return Vector2f(MAX_FLOAT, MAX_FLOAT);
    }
}

FLOAT ImplicitSurface2Set::closestDistanceLocal(const Vector2f& otherPoint) const {
    FLOAT min_dis = MAX_FLOAT;
    for(int i = 0; i < _surfaces.size(); ++i) {
        FLOAT dis = _surfaces[i]->closestDistance(otherPoint);
        if(min_dis > dis) {
            min_dis = dis;
        }
    }
    return min_dis;
}

Vector2f ImplicitSurface2Set::closestNormalLocal(const Vector2f& otherPoint) const {
    ImplicitSurface2Ptr s = nullptr;
    FLOAT min_dis = MAX_FLOAT;
    for(int i = 0; i < _surfaces.size(); ++i) {
        FLOAT dis = _surfaces[i]->closestDistance(otherPoint);
        if(min_dis > dis) {
            min_dis = dis;
            s = _surfaces[i];
        }
    }
    if(s != nullptr) {
        return s->closestNormal(otherPoint);
    } else {
        return Vector2f(MAX_FLOAT, MAX_FLOAT);
    }
}

bool ImplicitSurface2Set::intersectsLocal(const Ray2& ray) const {
    throw std::runtime_error("[ImplicitSurface2Set::intersectsLocal] not implemented yet");
    return false;
}

ShapeRayIntersection2 ImplicitSurface2Set::closestIntersectionLocal(const Ray2& ray) const {
    throw std::runtime_error("[ImplicitSurface2Set::closestIntersectionLocal] not implemented yet");
    return ShapeRayIntersection2();
}
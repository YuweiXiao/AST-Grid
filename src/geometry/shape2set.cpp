#include "geometry/shape2set.h"

using namespace Omni;
using namespace Geometry;

Shape2Set::Shape2Set() : Shape2(Transform2()), _shapes(std::vector<Shape2Ptr>()) 
{}

Shape2Set::Shape2Set(const std::vector<Shape2Ptr>& others,
                         const Transform2& transform, bool isNormalFlipped)
    : Shape2(transform, isNormalFlipped), _shapes(others)
{}

Vector2f Shape2Set::closestPointLocal(const Vector2f& otherPoint) const {
    Shape2Ptr s = nullptr;
    FLOAT min_dis = MAX_FLOAT;
    for(int i = 0; i < _shapes.size(); ++i) {
        FLOAT dis = _shapes[i]->closestDistance(otherPoint);
        if(min_dis > dis) {
            min_dis = dis;
            s = _shapes[i];
        }
    }

    if(s != nullptr) {
        return s->closestPoint(otherPoint);
    }
    return Vector2f(MAX_FLOAT, MAX_FLOAT);
}

Vector2f Shape2Set::closestNormalLocal(const Vector2f& otherPoint) const {
    Shape2Ptr s = nullptr;
    FLOAT min_dis = MAX_FLOAT;
    for(int i = 0; i < _shapes.size(); ++i) {
        FLOAT dis = _shapes[i]->closestDistance(otherPoint);
        if(min_dis > dis) {
            min_dis = dis;
            s = _shapes[i];
        }
    }

    if(s != nullptr) {
        return s->closestNormal(otherPoint);
    } 
    return Vector2f(1, 0);
}

FLOAT Shape2Set::closestDistanceLocal(const Vector2f& otherPoint) const {
    FLOAT min_dis = MAX_FLOAT;
    for(int i = 0; i < _shapes.size(); ++i) {
        FLOAT dis = _shapes[i]->closestDistance(otherPoint);
        if(min_dis > dis) {
            min_dis = dis;
        }
    }

    return min_dis;
}

bool Shape2Set::intersectsLocal(const Ray2& ray) const {
    throw std::runtime_error("[Shape2Set::intersectsLocal] not implemented yet");
    
    return true;
}

ShapeRayIntersection2 Shape2Set::closestIntersectionLocal(const Ray2& ray) const {
    throw std::runtime_error("[Shape2Set::closestIntersectionLocal] not implemented yet");
    return ShapeRayIntersection2();
}
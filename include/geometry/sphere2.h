#pragma once
#include "geometry/shape.h"


namespace Omni {
namespace Geometry {

class Sphere2 : public Shape2 {
public:
    Vector2f center;
    FLOAT radius = 1.0;

    Sphere2(const Vector2f& c, FLOAT r, 
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);
            
private:
    Vector2f closestPointLocal(const Vector2f& otherPoint) const override;

    FLOAT closestDistanceLocal(const Vector2f& otherPoint) const override;

    bool intersectsLocal(const Ray2& ray) const override;

    Vector2f closestNormalLocal(const Vector2f& otherPoint) const override;

    ShapeRayIntersection2 closestIntersectionLocal(
        const Ray2& ray) const override;

    FLOAT _squaredR;
};

typedef std::shared_ptr<Sphere2> Sphere2Ptr;

}   // end of namespace Geometry
}   // end of namespace Omni
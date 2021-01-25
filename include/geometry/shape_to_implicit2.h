#pragma once
#include "general.h"
#include "geometry/implicit_surface2.h"

namespace Omni {
namespace Geometry {

class ShapeToImplicit2 : public ImplicitSurface2 {
public:

    ShapeToImplicit2(const Shape2Ptr& shape, const Transform2& transform = Transform2(),
        bool isNormalFlipped = false) 
        :ImplicitSurface2(transform, isNormalFlipped), _shape(shape)
    {}

    //! Returns the raw Shape instance.
    Shape2Ptr shape() const { return _shape;}

protected:
    Vector2f closestPointLocal(const Vector2f& otherPoint) const override {
        return shape()->closestPoint(otherPoint);
    }

    FLOAT closestDistanceLocal(const Vector2f& otherPoint) const override {
        return shape()->closestDistance(otherPoint);
    }

    bool intersectsLocal(const Ray2& ray) const override {
        return shape()->intersects(ray);
    }

    Vector2f closestNormalLocal(const Vector2f& otherPoint) const override {
        return shape()->closestNormal(otherPoint);
    }

    ShapeRayIntersection2 closestIntersectionLocal(const Ray2& ray) const override {
        return shape()->closestIntersection(ray);
    }

    FLOAT signedDistanceLocal(const Vector2f& otherPoint) const override;

private:
    Shape2Ptr _shape;
};

typedef std::shared_ptr<ShapeToImplicit2> ShapeToImplicit2Ptr;

}  // end of namespace Geometry
}  // end of namespace Omni


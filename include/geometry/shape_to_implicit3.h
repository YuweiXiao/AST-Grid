#pragma once
#include "general.h"
#include "geometry/implicit_surface3.h"

namespace Omni {
namespace Geometry {

class ShapeToImplicit3 : public ImplicitSurface3 {
public:

    ShapeToImplicit3(const Shape3Ptr& shape, const Transform3& transform = Transform3(),
        bool isNormalFlipped = false) 
        :ImplicitSurface3(transform, isNormalFlipped), _shape(shape)
    {}

    //! Returns the raw Shape instance.
    const Shape3Ptr& shape() const { return _shape;}

protected:
    Vector3f closestPointLocal(const Vector3f& otherPoint) const override {
        return shape()->closestPoint(otherPoint);
    }

    REAL closestDistanceLocal(const Vector3f& otherPoint) const override {
        return shape()->closestDistance(otherPoint);
    }

    bool intersectsLocal(const Ray3& ray) const override {
        return shape()->intersects(ray);
    }

    Vector3f closestNormalLocal(const Vector3f& otherPoint) const override {
        return shape()->closestNormal(otherPoint);
    }

    ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override {
        return shape()->closestIntersection(ray);
    }

    REAL signedDistanceLocal(const Vector3f& otherPoint) const override {
        Vector3f x = shape()->closestPoint(otherPoint);
        Vector3f n = shape()->closestNormal(otherPoint);
        n = shape()->isNormalFlipped ? -n : n;
        return n.dot(otherPoint - x) < 0.0 ? -(x - otherPoint).norm() : (x - otherPoint).norm();
    }

private:
    Shape3Ptr _shape;
};

typedef std::shared_ptr<ShapeToImplicit3> ShapeToImplicit3Ptr;

}  // end of namespace Geometry
}  // end of namespace Omni

#pragma once
#include "geometry/shape.h"


namespace Omni {
namespace Geometry {

class Sphere3 : public Shape3{
public:
    REAL radius;

    Sphere3(REAL r, const Transform3& _transform = Transform3(),
                bool _isNormalFlipped = false) 
        : Shape3(_transform, _isNormalFlipped),radius(r) {
        _squaredR = r * r;
    }

protected:
    virtual bool intersectsLocal(const Ray3& ray) const override;
    virtual Vector3f closestPointLocal(const Vector3f& otherPoint) const override {
        return closestNormalLocal(otherPoint) * radius;
    }
    virtual ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override;
    virtual Vector3f closestNormalLocal(const Vector3f& otherPoint) const override {
        return otherPoint.norm() < EPSILON ? Vector3f(1, 0, 0) : otherPoint.normalized();
    }

private:
    REAL _squaredR;
};

typedef std::shared_ptr<Sphere3> Sphere3Ptr;

}   // end of namespace Geometry
}   // end of namespace Omni
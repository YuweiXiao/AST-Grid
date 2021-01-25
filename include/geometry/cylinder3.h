#pragma once
#include "geometry/shape.h"
#include "geometry/rectangle2.h"

namespace Omni {
namespace Geometry {

class Cylinder3 : public Shape3{
public:
    FLOAT radius;
	FLOAT height;

	//The cylinder is aligned with the y - axis.
    Cylinder3(FLOAT r, FLOAT h, const Transform3& _transform = Transform3(),
                bool _isNormalFlipped = false) 
        : Shape3(_transform, _isNormalFlipped), radius(r), height(h), rect2(Vector2f(2 * r, h)) {
        _squaredR = r * r;
    }

protected:
    virtual bool intersectsLocal(const Ray3& ray) const override;
    virtual ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override;
    
    virtual Vector3f closestPointLocal(const Vector3f& otherPoint) const override {
        Vector2f otherPoint2d(std::sqrt(otherPoint.x() * otherPoint.x() + otherPoint.z() * otherPoint.z()), otherPoint.y());

		Vector2f cp = rect2.closestPoint(otherPoint2d);
		double angle = std::atan2(otherPoint.z(), otherPoint.x());
		return Vector3f(cp.x() * std::cos(angle), cp.y(), cp.x() * std::sin(angle));
    }
    
    virtual Vector3f closestNormalLocal(const Vector3f& otherPoint) const override {
        Vector2f otherPoint2d(std::sqrt(otherPoint.x() * otherPoint.x() + otherPoint.z() * otherPoint.z()), otherPoint.y());

		Vector2f cn = rect2.closestNormal(otherPoint2d);
		if (cn.y() > 0) {
			return Vector3f(0, 1, 0);
		}
		else if (cn.y() < 0) {
			return Vector3f(0, -1, 0);
		}
		else {
			return Vector3f(otherPoint.x(), 0, otherPoint.z()).normalized();
		}
    }

private:
    FLOAT _squaredR;
    Rectangle2 rect2;
};

typedef std::shared_ptr<Cylinder3> Cylinder3Ptr;

}   // end of namespace Geometry
}   // end of namespace Omni
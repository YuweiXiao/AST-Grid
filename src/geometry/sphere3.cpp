#include "geometry/sphere3.h"

namespace Omni {
namespace Geometry {

ShapeRayIntersection3 Sphere3::closestIntersectionLocal(const Ray3& ray) const {
	ShapeRayIntersection3 intersection;
	Vector3f r = ray.origin;
	REAL b = ray.direction.dot(r);
	REAL c = r.squaredNorm() - _squaredR;
	REAL d = b * b - c;

	if (d > 0.) {
		d = std::sqrt(d);
		REAL tMin = -b - d;
		REAL tMax = -b + d;

		if (tMin < 0.0) {
			tMin = tMax;
		}

		if (tMin < 0.0) {
			intersection.isIntersecting = false;
		} else {
			intersection.isIntersecting = true;
			intersection.distance = tMin;
			intersection.point = ray.origin + tMin * ray.direction;
			intersection.normal = intersection.point.normalized();
		}
	} else {
		intersection.isIntersecting = false;
	}

	return intersection;
}

bool Sphere3::intersectsLocal(const Ray3& ray) const {
	Vector3f r = ray.origin;
	REAL b = ray.direction.dot(r);
	REAL c = r.squaredNorm() - _squaredR;
	REAL d = b * b - c;

	if (d > 0.) {
		d = std::sqrt(d);
		REAL tMin = -b - d;
		REAL tMax = -b + d;

		if (tMin < 0.0) {
			tMin = tMax;
		}

		if (tMin >= 0.0) {
			return true;
		}
	}

	return false;
}

} // end of namespace Geometry  
} // end of namespace Omni  
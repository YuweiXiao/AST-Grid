#include "geometry/shape_to_implicit2.h"

using namespace Omni::Geometry;


Real Omni::Geometry::ShapeToImplicit2::signedDistanceLocal(const Vector2f & otherPoint) const
{
	Vector2f x = shape()->closestPoint(otherPoint);
	Vector2f n = shape()->closestNormal(otherPoint);
	n = shape()->isNormalFlipped ? -n : n;
	if (n.dot(otherPoint - x) < 0.0) {
		return -(x - otherPoint).norm();
	}
	else {
		return (x - otherPoint).norm();
	}
}
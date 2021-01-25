#include "geometry/triangle3.h"

namespace Omni {
namespace Geometry {
	inline Vector3f closestPointOnLine(const Vector3f& v0, const Vector3f& v1, const Vector3f& pt) {
		const double lenSquared = (v1 - v0).squaredNorm();
		if (lenSquared < EPSILON) {
			return v0;
		}
		const double t = (pt - v0).dot(v1 - v0) / lenSquared;
		if (t < 0.0) {
			return v0;
		}
		else if (t > 1.0) {
			return v1;
		}
		return v0 + t * (v1 - v0);
	}

	inline Vector3f closestNormalOnLine(const Vector3f& v0, const Vector3f& v1, const Vector3f& n0, 
									const Vector3f& n1, const Vector3f& pt) {
		const double lenSquared = (v1 - v0).squaredNorm();
		if (lenSquared < EPSILON) {
			return n0;
		}
		const double t = (pt - v0).dot(v1 - v0) / lenSquared;
		if (t < 0.0) {
			return n0;
		}
		else if (t > 1.0) {
			return n1;
		}
		return (n0 + t * (n1 - n0)).normalized();
	}

	Triangle3::Triangle3(const std::array<Vector3f, 3>& newPoints,
		const std::array<Vector3f, 3>& newNormals,
		const std::array<Vector2f, 3>& newUvs,
		const Transform3& transform_, bool isNormalFlipped_) :
		Shape3(transform_, isNormalFlipped_),
		points(newPoints),
		normals(newNormals),
		uvs(newUvs) {}

	double Triangle3::area() const {
		return 0.5 * (points[1] - points[0]).cross(points[2] - points[0]).norm();
	}

	void Triangle3::getBarycentricCoords(const Vector3f& pt, double* b0, double* b1,
		double* b2) const {
		Vector3f q01 = (points[1] - points[0]).cross(pt - points[0]);
		Vector3f q12 = (points[2] - points[1]).cross(pt - points[1]);
		Vector3f q02 = (points[0] - points[2]).cross(pt - points[2]);

		double a = area();
		*b0 = 0.5 * q12.norm() / a;
		*b1 = 0.5 * q02.norm() / a;
		*b2 = 0.5 * q01.norm() / a;
	}

	Vector3f Triangle3::faceNormal() const {
		Vector3f ret = (points[1] - points[0]).cross(points[2] - points[0]);
		return ret.normalized();
	}

	void Triangle3::setNormalsToFaceNormal() {
		normals[0] = normals[1] = normals[2] = faceNormal();
	}

	Vector3f Triangle3::closestPointLocal(const Vector3f& otherPoint) const {
		Vector3f n = faceNormal();
		double nd = n.dot(n);
		double d = n.dot(points[0]);
		double t = (d - n.dot(otherPoint)) / nd;

		Vector3f q = t * n + otherPoint;

		Vector3f q01 = (points[1] - points[0]).cross(q - points[0]);
		if (n.dot(q01) < 0) {
			return closestPointOnLine(points[0], points[1], q);
		}

		Vector3f q12 = (points[2] - points[1]).cross(q - points[1]);
		if (n.dot(q12) < 0) {
			return closestPointOnLine(points[1], points[2], q);
		}

		Vector3f q02 = (points[0] - points[2]).cross(q - points[2]);
		if (n.dot(q02) < 0) {
			return closestPointOnLine(points[0], points[2], q);
		}

		double a = area();
		double b0 = 0.5 * q12.norm() / a;
		double b1 = 0.5 * q02.norm() / a;
		double b2 = 0.5 * q01.norm() / a;

		return b0 * points[0] + b1 * points[1] + b2 * points[2];
	}

	Vector3f Triangle3::closestNormalLocal(const Vector3f& otherPoint) const {
		Vector3f n = faceNormal();
		double nd = n.dot(n);
		double d = n.dot(points[0]);
		double t = (d - n.dot(otherPoint)) / nd;

		Vector3f q = t * n + otherPoint;

		Vector3f q01 = (points[1] - points[0]).cross(q - points[0]);
		if (n.dot(q01) < 0) {
			return closestNormalOnLine(points[0], points[1], normals[0], normals[1],
				q);
		}

		Vector3f q12 = (points[2] - points[1]).cross(q - points[1]);
		if (n.dot(q12) < 0) {
			return closestNormalOnLine(points[1], points[2], normals[1], normals[2],
				q);
		}

		Vector3f q02 = (points[0] - points[2]).cross(q - points[2]);
		if (n.dot(q02) < 0) {
			return closestNormalOnLine(points[0], points[2], normals[0], normals[2],
				q);
		}

		double a = area();
		double b0 = 0.5 * q12.norm() / a;
		double b1 = 0.5 * q02.norm() / a;
		double b2 = 0.5 * q01.norm() / a;

		return (b0 * normals[0] + b1 * normals[1] + b2 * normals[2]).normalized();
	}

	ShapeRayIntersection3 Triangle3::closestIntersectionLocal(const Ray3 & ray) const
	{
		throw std::runtime_error("[Triangle3::closestIntersectionLocal] not implemented yet");
		return ShapeRayIntersection3();
	}

	bool Triangle3::intersectsLocal(const Ray3 & ray) const
	{
		throw std::runtime_error("[Triangle3::intersectionLocal] not implemented yet");
		return false;
	}

	BoundingBox3 Triangle3::boundingBoxLocal() const {
		BoundingBox3 box(points[0], points[1]);
		box.merge(points[2]);
		return box;
	}
}
}
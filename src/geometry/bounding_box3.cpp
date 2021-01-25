#include "geometry/bounding_box3.h"

namespace Omni {
namespace Geometry {


	BoundingBox3::BoundingBox3() {

	}

	BoundingBox3::BoundingBox3(const Vector3f& point1, const Vector3f& point2) {
		lowerCorner.x() = std::min(point1.x(), point2.x());
		lowerCorner.y() = std::min(point1.y(), point2.y());
		lowerCorner.z() = std::min(point1.z(), point2.z());
		upperCorner.x() = std::max(point1.x(), point2.x());
		upperCorner.y() = std::max(point1.y(), point2.y());
		upperCorner.z() = std::max(point1.z(), point2.z());
	}

	FLOAT BoundingBox3::width() const {
		return upperCorner.x() - lowerCorner.x();
	}

	FLOAT BoundingBox3::height() const {
		return upperCorner.y() - lowerCorner.y();
	}

	FLOAT BoundingBox3::depth() const {
		return upperCorner.z() - lowerCorner.z();
	}

	Vector3f BoundingBox3::midPoint() const {
		return (upperCorner + lowerCorner) / static_cast<FLOAT>(2);
	}

	void BoundingBox3::reset() {
		lowerCorner.x() = std::numeric_limits<FLOAT>::max();
		lowerCorner.y() = std::numeric_limits<FLOAT>::max();
		lowerCorner.z() = std::numeric_limits<FLOAT>::max();
		upperCorner.x() = -std::numeric_limits<FLOAT>::max();
		upperCorner.y() = -std::numeric_limits<FLOAT>::max();
		upperCorner.z() = -std::numeric_limits<FLOAT>::max();
	}

	Vector3f BoundingBox3::corner(size_t idx) const {
		static const FLOAT h = 1.0 / 2.0;
		static const Vector3f offset[8] = {
			{-h, -h, -h}, {+h, -h, -h}, {-h, +h, -h}, {+h, +h, -h},
			{-h, -h, +h}, {+h, -h, +h}, {-h, +h, +h}, {+h, +h, +h} };

		return midPoint() + Vector3f(width() * offset[idx][0], height() * offset[idx][1], depth() * offset[idx][2]);
	}

	bool BoundingBox3::overlaps(const BoundingBox3& other) const {
		if (upperCorner.x() < other.lowerCorner.x() ||
			lowerCorner.x() > other.upperCorner.x()) {
			return false;
		}
		if (upperCorner.y() < other.lowerCorner.y() ||
			lowerCorner.y() > other.upperCorner.y()) {
			return false;
		}
		if (upperCorner.z() < other.lowerCorner.z() ||
			lowerCorner.z() > other.upperCorner.z()) {
			return false;
		}
		return true;
	}

	bool BoundingBox3::contains(const Vector3f& point) const {
		if (upperCorner.x() < point.x() || lowerCorner.x() > point.x()) {
			return false;
		}
		if (upperCorner.y() < point.y() || lowerCorner.y() > point.y()) {
			return false;
		}
		if (upperCorner.z() < point.z() || lowerCorner.z() > point.z()) {
			return false;
		}
		return true;
	}

	bool BoundingBox3::intersects(const Ray3& ray) const {
		FLOAT tMin = 0;
		FLOAT tMax = std::numeric_limits<FLOAT>::max();
		const Vector3f& rayInvDir = Vector3f(1.0 / ray.direction.x(), 1.0 / ray.direction.y(), 1.0 / ray.direction.z());

		for (int i = 0; i < 3; ++i) {
			FLOAT tNear = (lowerCorner[i] - ray.origin[i]) * rayInvDir[i];
			FLOAT tFar = (upperCorner[i] - ray.origin[i]) * rayInvDir[i];
			if (tNear > tFar) std::swap(tNear, tFar);
			tMin = tNear > tMin ? tNear : tMin;
			tMax = tFar < tMax ? tFar : tMax;
			if (tMin > tMax) return false;
		}
		return true;
	}

	BoundingBoxRayIntersection3 BoundingBox3::closestIntersection(
		const Ray3& ray) const {
		BoundingBoxRayIntersection3 intersection;

		FLOAT tMin = 0;
		FLOAT tMax = std::numeric_limits<FLOAT>::max();
		const Vector3f& rayInvDir = Vector3f(1.0 / ray.direction.x(), 1.0 / ray.direction.y(), 1.0 / ray.direction.z());

		for (int i = 0; i < 3; ++i) {
			FLOAT tNear = (lowerCorner[i] - ray.origin[i]) * rayInvDir[i];
			FLOAT tFar = (upperCorner[i] - ray.origin[i]) * rayInvDir[i];

			if (tNear > tFar) std::swap(tNear, tFar);
			tMin = tNear > tMin ? tNear : tMin;
			tMax = tFar < tMax ? tFar : tMax;

			if (tMin > tMax) {
				intersection.isIntersecting = false;
				return intersection;
			}
		}

		intersection.isIntersecting = true;

		if (contains(ray.origin)) {
			intersection.tNear = tMax;
			intersection.tFar = std::numeric_limits<FLOAT>::max();
		}
		else {
			intersection.tNear = tMin;
			intersection.tFar = tMax;
		}

		return intersection;
	}


	FLOAT BoundingBox3::diagonalLength() const {
		return (upperCorner - lowerCorner).norm();
	}

	FLOAT BoundingBox3::diagonalLengthSquared() const {
		return (upperCorner - lowerCorner).squaredNorm();
	}

	void BoundingBox3::merge(const Vector3f& point) {
		lowerCorner.x() = std::min(lowerCorner.x(), point.x());
		lowerCorner.y() = std::min(lowerCorner.y(), point.y());
		lowerCorner.z() = std::min(lowerCorner.z(), point.z());
		upperCorner.x() = std::max(upperCorner.x(), point.x());
		upperCorner.y() = std::max(upperCorner.y(), point.y());
		upperCorner.z() = std::max(upperCorner.z(), point.z());
	}

	void BoundingBox3::merge(const BoundingBox3& other) {
		lowerCorner.x() = std::min(lowerCorner.x(), other.lowerCorner.x());
		lowerCorner.y() = std::min(lowerCorner.y(), other.lowerCorner.y());
		lowerCorner.z() = std::min(lowerCorner.z(), other.lowerCorner.z());
		upperCorner.x() = std::max(upperCorner.x(), other.upperCorner.x());
		upperCorner.y() = std::max(upperCorner.y(), other.upperCorner.y());
		upperCorner.z() = std::max(upperCorner.z(), other.upperCorner.z());
	}


	bool BoundingBox3::isEmpty() const {
		return (lowerCorner.x() >= upperCorner.x() || lowerCorner.y() >= upperCorner.y() ||
			lowerCorner.z() >= upperCorner.z());
	}

	Vector3f Omni::Geometry::BoundingBox3::clamp(const Vector3f & point) const
	{
		FLOAT x, y, z;
		x = std::min(point.x(), upperCorner.x());
		x = std::max(point.x(), lowerCorner.x());
		y = std::min(point.y(), upperCorner.y());
		y = std::max(point.y(), lowerCorner.y());
		z = std::min(point.z(), upperCorner.z());
		z = std::max(point.z(), lowerCorner.z());
		return Vector3f(x, y, z);
	}
}
}
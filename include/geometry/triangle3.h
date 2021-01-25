#pragma once
#include "geometry/shape.h"
#include "geometry/bounding_box3.h"
#include <array>

namespace Omni {
namespace Geometry {

	class Triangle3 : public Shape3 {
	public:
		//! Three points.
		std::array<Vector3f, 3> points;
		//! Three normals.
		std::array<Vector3f, 3> normals;
		//! Three UV coordinates.
		std::array<Vector2f, 3> uvs;

		Triangle3() {}

		//! Constructs a triangle with given \p points, \p normals, and \p uvs.
		Triangle3(const std::array<Vector3f, 3>& points,
					const std::array<Vector3f, 3>& normals,
					const std::array<Vector2f, 3>& uvs,
					const Transform3& transform = Transform3(),
					bool isNormalFlipped = false);

		//! Returns the area of this triangle.
		double area() const;

		//! Returns barycentric coordinates for the given point \p pt.
		void getBarycentricCoords(const Vector3f& pt, double* b0, double* b1, double* b2) const;

		//! Returns the face normal of the triangle.
		Vector3f faceNormal() const;

		//! Set Triangle3::normals to the face normal.
		void setNormalsToFaceNormal();


	protected:
		virtual bool intersectsLocal(const Ray3& ray) const override;
		virtual Vector3f closestPointLocal(const Vector3f& otherPoint) const override;
		virtual ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override;
		virtual Vector3f closestNormalLocal(const Vector3f& otherPoint) const override;
		virtual BoundingBox3 boundingBoxLocal() const override;
	};

	//! Shared pointer for the Triangle3 type.
	typedef std::shared_ptr<Triangle3> Triangle3Ptr;
}
}
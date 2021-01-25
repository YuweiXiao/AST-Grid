#pragma once

#include "geometry/shape.h"
#include "geometry/triangle3.h"
#include "geometry/bounding_box3.h"
#include "geometry/bvh3.h"
#include <vector>

namespace Omni {
namespace Geometry {

class TriangleMesh3 : public Shape3 {
public:
	typedef std::vector<Vector2f> Vector2fArray;
	typedef std::vector<Vector3f> Vector3fArray;
	typedef std::vector<Size3> IndexArray;

	typedef Vector3fArray PointArray;
	typedef Vector3fArray NormalArray;
	typedef Vector2fArray UvArray;

	//! Default constructor.
	TriangleMesh3(const Transform3& transform = Transform3(), bool isNormalFlipped = false);

	//! Constructs mesh with points, normals, uvs, and their indices.
	TriangleMesh3(const PointArray& points, const NormalArray& normals,
		const UvArray& uvs, const IndexArray& pointIndices,
		const IndexArray& normalIndices, const IndexArray& uvIndices,
		const Transform3& transform_ = Transform3(),
		bool isNormalFlipped = false);

	//! Returns number of points.
	size_t numberOfPoints() const {return _points.size();}

	//! Returns number of normals.
	size_t numberOfNormals() const {return _normals.size();}

	//! Returns number of UV coordinates.
	size_t numberOfUvs() const {return _uvs.size();}

	//! Returns number of triangles.
	size_t numberOfTriangles() const {return _pointIndices.size();}

	//! Returns constant reference to the i-th point.
	const Vector3f& point(size_t i) const {return _points[i];}

	//! Returns reference to the i-th point.
	Vector3f& point(size_t i) {
		invalidateBvh();
		return _points[i];
	}

	//! Returns constant reference to the point indices of i-th triangle.
	const Size3& pointIndex(size_t i) const { return _pointIndices[i]; }

	//! Returns reference to the point indices of i-th triangle.
	Size3& pointIndex(size_t i) { return _pointIndices[i];}

	//! Returns constant reference to the i-th normal.
	const Vector3f& normal(size_t i) const {return _normals[i];}

	//! Returns reference to the i-th normal.
	Vector3f& normal(size_t i) {return _normals[i];}

	//! Returns true if the mesh has normals.
	bool hasNormals() const {return _normals.size() > 0;}

	//! Returns true if the mesh has UV coordinates.
	bool hasUvs() const {return _uvs.size() > 0;}

	//! Returns i-th triangle.
	Triangle3 triangle(size_t i) const;

	//! Adds a point.
	void addPoint(const Vector3f& pt) {_points.push_back(pt);}

	//! Adds a normal.
	void addNormal(const Vector3f& n) {_normals.push_back(n);}

	//! Adds a UV.
	void addUv(const Vector2f& t) {_uvs.push_back(t);}

	//! Adds a triangle with points.
	void addPointTriangle(const Size3& newPointIndices) {
		_pointIndices.push_back(newPointIndices);
		invalidateBvh();
	}

	//! Adds a triangle with point and normal.
	void addPointNormalTriangle(const Size3& newPointIndices, const Size3& newNormalIndices) {
		_pointIndices.push_back(newPointIndices);
		_normalIndices.push_back(newNormalIndices);
		invalidateBvh();
	}

	//! Adds a triangle with point and UV.
	void addPointUvTriangle(const Size3& newPointIndices, const Size3& newUvIndices) {
		_pointIndices.push_back(newPointIndices);
		_uvIndices.push_back(newUvIndices);
		invalidateBvh();
	}

	//! Adds a triangle with point, normal, and UV.
	void addPointUvNormalTriangle(const Size3& newPointIndices,
								const Size3& newUvIndices,
								const Size3& newNormalIndices) {
		_pointIndices.push_back(newPointIndices);
		_normalIndices.push_back(newNormalIndices);
		_uvIndices.push_back(newUvIndices);
		invalidateBvh();
	}

	void addNormalTriangle(const Size3& newNormalIndices) {
		_normalIndices.push_back(newNormalIndices);
		invalidateBvh();
	}

	void addUvTriangle(const Size3& newUvIndices) {
		_uvIndices.push_back(newUvIndices);
		invalidateBvh();
	}


	//! Reads the mesh in obj format from the input stream.
	bool readObj(std::istream* strm);
	//! Reads the mesh in obj format from the file.
	bool readObj(const std::string& filename);

	//! Writes the mesh in obj format to the file.
	bool writeObj(const std::string& filename) const;

protected:
	bool intersectsLocal(const Ray3& ray) const override {
		throw std::runtime_error("[intersectsLocal] not implemented yet");
		return true;
	}
	Vector3f closestPointLocal(const Vector3f& otherPoint) const override {
		buildBvh();

		const auto distanceFunc = [this](const size_t& triIdx, const Vector3f& pt) {
			Triangle3 tri = triangle(triIdx);
			return tri.closestDistance(pt);
		};
		const auto queryResult = _bvh.nearest(otherPoint, distanceFunc);
		return triangle(*queryResult.item).closestPoint(otherPoint);
	}
	ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override {
		throw std::runtime_error("[closestIntersectionLocal] not implemented yet");
		return ShapeRayIntersection3();
	}
	Vector3f closestNormalLocal(const Vector3f& otherPoint) const override {
		buildBvh();

		const auto distanceFunc = [this](const size_t& triIdx, const Vector3f& pt) {
			Triangle3 tri = triangle(triIdx);
			return tri.closestDistance(pt);
		};

		const auto queryResult = _bvh.nearest(otherPoint, distanceFunc);
		//    printf("%zu\n", *queryResult.item);
		return triangle(*queryResult.item).closestNormal(otherPoint);
	}

	BoundingBox3 boundingBoxLocal() const override {
		buildBvh();
		return _bvh.boundingBox();
	}

private:
	PointArray _points;
	NormalArray _normals;
	UvArray _uvs;
	IndexArray _pointIndices;
	IndexArray _normalIndices;
	IndexArray _uvIndices;

	mutable Bvh3<size_t> _bvh;
	mutable bool _bvhInvalidated = true;

	void invalidateBvh() { _bvhInvalidated = true; }

	void buildBvh() const {
		if (_bvhInvalidated) {
			size_t nTris = numberOfTriangles();
			std::vector<size_t> ids(nTris);
			std::vector<BoundingBox3> bounds(nTris);
			for (size_t i = 0; i < nTris; ++i) {
				ids[i] = i;
				bounds[i] = triangle(i).boundingBox();
			}
			_bvh.build(ids, bounds);
			_bvhInvalidated = false;
		} 
	}
};

//! Shared pointer for the TriangleMesh3 type.
typedef std::shared_ptr<TriangleMesh3> TriangleMesh3Ptr;

} // end of namespace Geometry
}	// end of namespace Omni
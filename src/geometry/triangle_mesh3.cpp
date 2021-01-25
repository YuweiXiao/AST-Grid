#include "geometry/triangle_mesh3.h"
#include <fstream>
#define TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_USE_DOUBLE
#include "tinyobj/tiny_obj_loader.h"

namespace Omni {
namespace Geometry {

inline std::ostream& operator<<(std::ostream& strm, const Vector2f& v) {
    strm << v.x() << ' ' << v.y();
    return strm;
}

inline std::ostream& operator<<(std::ostream& strm, const Vector3f& v) {
    strm << v.x() << ' ' << v.y() << ' ' << v.z();
    return strm;
}

TriangleMesh3::TriangleMesh3(const Transform3& transform_, bool isNormalFlipped_): 
	Shape3(transform_, isNormalFlipped_) {}

TriangleMesh3::TriangleMesh3(const PointArray& points, const NormalArray& normals, const UvArray& uvs,
	const IndexArray& pointIndices, const IndexArray& normalIndices,
	const IndexArray& uvIndices, const Transform3& transform_,
	bool isNormalFlipped_)
	: Shape3(transform_, isNormalFlipped_),
	_points(points),
	_normals(normals),
	_uvs(uvs),
	_pointIndices(pointIndices),
	_normalIndices(normalIndices),
	_uvIndices(uvIndices) {}
	
Triangle3 TriangleMesh3::triangle(size_t i) const {
	Triangle3 tri;
	for (int j = 0; j < 3; j++) {
		tri.points[j] = _points[_pointIndices[i][j]];
		if (hasUvs()) {
			tri.uvs[j] = _uvs[_uvIndices[i][j]];
		}
	}
	Vector3f n = tri.faceNormal();
	for (int j = 0; j < 3; j++) {
		if (hasNormals()) {
			tri.normals[j] = _normals[_normalIndices[i][j]];
		}
		else {
			tri.normals[j] = n;
		}
	}
	return tri;
}

bool TriangleMesh3::readObj(std::istream* strm) {
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	const bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, strm);
	// `err` may contain warning message.
	if (!err.empty()) {
		spdlog::error("{}", err);
		return false;
	}
	// Failed to load obj.
	if (!ret) {
		return false;
	}
	invalidateBvh();
	// Read vertices
	for (size_t idx = 0; idx < attrib.vertices.size() / 3; ++idx) {
		// Access to vertex
		tinyobj::real_t vx = attrib.vertices[3 * idx + 0];
		tinyobj::real_t vy = attrib.vertices[3 * idx + 1];
		tinyobj::real_t vz = attrib.vertices[3 * idx + 2];
		addPoint({ vx, vy, vz });
	}
	// Read normals
	for (size_t idx = 0; idx < attrib.normals.size() / 3; ++idx) {
		// Access to normal
		tinyobj::real_t vx = attrib.normals[3 * idx + 0];
		tinyobj::real_t vy = attrib.normals[3 * idx + 1];
		tinyobj::real_t vz = attrib.normals[3 * idx + 2];
		addNormal({ vx, vy, vz });
	}
	// Read UVs
	for (size_t idx = 0; idx < attrib.texcoords.size() / 2; ++idx) {
		// Access to UV
		tinyobj::real_t tu = attrib.texcoords[2 * idx + 0];
		tinyobj::real_t tv = attrib.texcoords[2 * idx + 1];
		addUv({ tu, tv });
	}
	// Read faces
	for (auto& shape : shapes) {
		size_t idx = 0;
		for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
			const size_t fv = shape.mesh.num_face_vertices[f];
			if (fv == 3) {
				if (!attrib.vertices.empty()) {
					addPointTriangle(
						Size3( shape.mesh.indices[idx].vertex_index,
						 shape.mesh.indices[idx + 1].vertex_index,
						 shape.mesh.indices[idx + 2].vertex_index ));
				}
				if (!attrib.normals.empty()) {
					addNormalTriangle(
						Size3( shape.mesh.indices[idx].normal_index,
						 shape.mesh.indices[idx + 1].normal_index,
						 shape.mesh.indices[idx + 2].normal_index ));
				}
				if (!attrib.texcoords.empty()) {
					addUvTriangle(
						Size3( shape.mesh.indices[idx].texcoord_index,
						 shape.mesh.indices[idx + 1].texcoord_index,
						 shape.mesh.indices[idx + 2].texcoord_index ));
				}
			}
			idx += fv;
		}
	}
	return true;
}

bool TriangleMesh3::readObj(const std::string& filename) {
	std::ifstream file(filename.c_str());
	if (file) {
		bool result = readObj(&file);
		file.close();

		return result;
	}
	else {
		spdlog::error("read file fail:{}", filename);
		return false;
	}
}

//! Writes the mesh in obj format to the file.
bool TriangleMesh3::writeObj(const std::string& filename) const {
	std::ofstream file(filename.c_str());
	if (file) {
		// vertex
		for (const auto& pt : _points) {
			file << "v " << pt << std::endl;
		}

		// uv coords
		for (const auto& uv : _uvs) {
			file << "vt " << uv << std::endl;
		}

		// normals
		for (const auto& n : _normals) {
			file << "vn " << n << std::endl;
		}

		// faces
		bool hasUvs_ = hasUvs();
		bool hasNormals_ = hasNormals();
		for (size_t i = 0; i < numberOfTriangles(); ++i) {
			file << "f ";
			for (int j = 0; j < 3; ++j) {
				file << _pointIndices[i][j] + 1;
				if (hasNormals_ || hasUvs_) {
					file << '/';
				}
				if (hasUvs_) {
					file << _uvIndices[i][j] + 1;
				}
				if (hasNormals_) {
					file << '/' << _normalIndices[i][j] + 1;
				}
				file << ' ';
			}
			file << std::endl;
		}

		file.close();
		return true;
	} else {
		return false;
	}
}


} // end of namespace Geometry
} // end of namespace Omni
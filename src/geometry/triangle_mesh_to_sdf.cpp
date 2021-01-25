#include "geometry/triangle_mesh_to_sdf.h"

#ifdef WIN32
typedef long long ssize_t;
#endif

namespace Omni {
namespace Geometry {
	template <typename T>
	inline T min3(T x, T y, T z) {
		return std::min(std::min(x, y), z);
	}

	template <typename T>
	inline T max3(T x, T y, T z) {
		return std::max(std::max(x, y), z);
	}

	template <typename T>
	inline T clamp(T val, T low, T high) {
		if (val < low) {
			return low;
		}
		else if (val > high) {
			return high;
		}
		else {
			return val;
		}
	}

	static void checkNeighbor(
		const TriangleMesh3& mesh,
		const Vector3f& gx,
		ssize_t i0,
		ssize_t j0,
		ssize_t k0,
		ssize_t i1,
		ssize_t j1,
		ssize_t k1,
		CellCenteredScalarGrid3Ptr sdf,
		CellCenteredFlagGrid3Ptr closestTri) {
		if (closestTri->get(Size3(i1, j1, k1)) != std::numeric_limits<size_t>::max()) {
			size_t t = closestTri->get(Size3(i1, j1, k1));
			Triangle3 tri = mesh.triangle(t);

			double d = tri.closestDistance(gx);

			if (d < sdf->get(Size3(i0, j0, k0))) {
				sdf->get(Size3(i0, j0, k0)) = d;
				closestTri->get(Size3(i0, j0, k0)) = closestTri->get(Size3(i1, j1, k1));
			}
		}
	}

	static void sweep(
		const TriangleMesh3& mesh,
		int di,
		int dj,
		int dk,
		CellCenteredScalarGrid3Ptr sdf,
		CellCenteredFlagGrid3Ptr closestTri) {
		Size3 size = sdf->resolution();
		REAL h = sdf->gridSpacing();
		Vector3f origin = Vector3f(0, 0, 0);

		ssize_t ni = static_cast<ssize_t>(size.x());
		ssize_t nj = static_cast<ssize_t>(size.y());
		ssize_t nk = static_cast<ssize_t>(size.z());

		ssize_t i0, i1;
		if (di > 0) {
			i0 = 1;
			i1 = ni;
		}
		else {
			i0 = ni - 2;
			i1 = -1;
		}

		ssize_t j0, j1;
		if (dj > 0) {
			j0 = 1;
			j1 = nj;
		}
		else {
			j0 = nj - 2;
			j1 = -1;
		}

		ssize_t k0, k1;
		if (dk > 0) {
			k0 = 1;
			k1 = nk;
		}
		else {
			k0 = nk - 2;
			k1 = -1;
		}

		for (ssize_t k = k0; k != k1; k += dk) {
			for (ssize_t j = j0; j != j1; j += dj) {
				for (ssize_t i = i0; i != i1; i += di) {
					Vector3f gx(i * h, j * h, k * h);
					//gx *= h;
					//gx += origin;

					checkNeighbor(
						mesh, gx, i, j, k, i - di, j, k, sdf, closestTri);
					checkNeighbor(
						mesh, gx, i, j, k, i, j - dj, k, sdf, closestTri);
					checkNeighbor(
						mesh, gx, i, j, k, i - di, j - dj, k, sdf, closestTri);
					checkNeighbor(
						mesh, gx, i, j, k, i, j, k - dk, sdf, closestTri);
					checkNeighbor(
						mesh, gx, i, j, k, i - di, j, k - dk, sdf, closestTri);
					checkNeighbor(
						mesh, gx, i, j, k, i, j - dj, k - dk, sdf, closestTri);
					checkNeighbor(
						mesh, gx, i, j, k, i - di, j - dj, k - dk, sdf, closestTri);
				}
			}
		}
	}

	// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
	// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate
	// triangle)
	static int orientation(
		double x1,
		double y1,
		double x2,
		double y2,
		double* twiceSignedArea) {
		(*twiceSignedArea) = y1 * x2 - x1 * y2;

		if ((*twiceSignedArea) > 0) {
			return 1;
		}
		else if ((*twiceSignedArea) < 0) {
			return -1;
		}
		else if (y2 > y1) {
			return 1;
		}
		else if (y2 < y1) {
			return -1;
		}
		else if (x1 > x2) {
			return 1;
		}
		else if (x1 < x2) {
			return -1;
		}
		else {
			// only true when x1==x2 and y1==y2
			return 0;
		}
	}

	// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
	// if true is returned, the barycentric coordinates are set in a,b,c.
	static bool pointInTriangle2D(
		double x0, double y0,
		double x1, double y1, double x2, double y2, double x3, double y3,
		double* a, double* b, double* c) {
		x1 -= x0;
		x2 -= x0;
		x3 -= x0;
		y1 -= y0;
		y2 -= y0;
		y3 -= y0;

		int signa = orientation(x2, y2, x3, y3, a);
		if (signa == 0) {
			return false;
		}

		int signb = orientation(x3, y3, x1, y1, b);
		if (signb != signa) {
			return false;
		}

		int signc = orientation(x1, y1, x2, y2, c);
		if (signc != signa) {
			return false;
		}

		double sum = (*a) + (*b) + (*c);

		// if the SOS signs match and are nonkero, there's no way all of a, b, and c
		// are zero.
		//JET_ASSERT(sum != 0);

		*a /= sum;
		*b /= sum;
		*c /= sum;
		return true;
	}

	void triangleMeshToSdf(
		const TriangleMesh3& mesh,
		CellCenteredScalarGrid3Ptr sdf,
		const unsigned int exactBand) {
		Size3 size = sdf->resolution();

		// Upper bound on distance
		Vector3f domain(sdf->resolution().x() * sdf->gridSpacing(),
			sdf->resolution().y() * sdf->gridSpacing(),
			sdf->resolution().z() * sdf->gridSpacing());
		sdf->fill(sqrt(domain.x() * domain.x() + domain.y() * domain.y() + domain.z() * domain.z()));
		REAL h = sdf->gridSpacing();

		//! origin: data position for the grid point at (0, 0, 0).
		//! Note that this is different from origin() since origin() returns
		//! the lower corner point of the bounding box.
		//!
		Vector3f origin(0, 0, 0);

		//Array3<size_t> closestTri(size, kMaxSize);
		CellCenteredFlagGrid3Ptr closestTri = 
			std::shared_ptr<CellCenteredFlagGrid3>(new CellCenteredFlagGrid3(size, h));
		closestTri->fill(std::numeric_limits<int>::max());
		// Intersection_count(i,j,k) is # of tri intersections in (i-1,i]x{j}x{k}
		//Array3<unsigned int> intersectionCount(size, 0);
		CellCenteredFlagGrid3Ptr intersectionCount = 
			std::shared_ptr<CellCenteredFlagGrid3>(new CellCenteredFlagGrid3(size, h));
		closestTri->fill(0);
		// We begin by initializing distances near the mesh, and figuring out
		// intersection counts

		//auto gridPos = sdf->dataPosition();

		size_t nTri = mesh.numberOfTriangles();
		ssize_t bandwidth = static_cast<ssize_t>(exactBand);
		ssize_t maxSizeX = static_cast<ssize_t>(size.x());
		ssize_t maxSizeY = static_cast<ssize_t>(size.y());
		ssize_t maxSizeZ = static_cast<ssize_t>(size.z());
		for (size_t t = 0; t < nTri; ++t) {
			Size3 indices = mesh.pointIndex(t);

			Triangle3 tri = mesh.triangle(t);

			Vector3f pt1 = mesh.point(indices.x());
			Vector3f pt2 = mesh.point(indices.y());
			Vector3f pt3 = mesh.point(indices.z());

			// Normalize coordinates
			//Vector3f f1 = (pt1 - origin) / h;
			//Vector3f f2 = (pt2 - origin) / h;
			//Vector3f f3 = (pt3 - origin) / h;
			Vector3f f1 = Vector3f(pt1.x() / h, pt1.y() / h, pt1.z() / h);
			Vector3f f2 = Vector3f(pt2.x() / h, pt2.y() / h, pt2.z() / h);
			Vector3f f3 = Vector3f(pt3.x() / h, pt3.y() / h, pt3.z() / h);

			// Do distances nearby
			ssize_t i0 = static_cast<ssize_t>(min3<double>(f1.x(), f2.x(), f3.x()));
			i0 = clamp<ssize_t>(i0 - bandwidth, 0, maxSizeX - 1);
			ssize_t i1 = static_cast<ssize_t>(max3<double>(f1.x(), f2.x(), f3.x()));
			i1 = clamp<ssize_t>(i1 + bandwidth + 1, 0, maxSizeX - 1);

			ssize_t j0 = static_cast<ssize_t>(min3<double>(f1.y(), f2.y(), f3.y()));
			j0 = clamp<ssize_t>(j0 - bandwidth, 0, maxSizeY - 1);
			ssize_t j1 = static_cast<ssize_t>(max3<double>(f1.y(), f2.y(), f3.y()));
			j1 = clamp<ssize_t>(j1 + bandwidth + 1, 0, maxSizeY - 1);

			ssize_t k0 = static_cast<ssize_t>(min3<double>(f1.z(), f2.z(), f3.z()));
			k0 = clamp<ssize_t>(k0 - bandwidth, 0, maxSizeZ - 1);
			ssize_t k1 = static_cast<ssize_t>(max3<double>(f1.z(), f2.z(), f3.z()));
			k1 = clamp<ssize_t>(k1 + bandwidth + 1, 0, maxSizeZ - 1);

			for (ssize_t k = k0; k <= k1; ++k) {
				for (ssize_t j = j0; j <= j1; ++j) {
					for (ssize_t i = i0; i <= i1; ++i) {
						//Vector3f gx = gridPos(i, j, k);
						Vector3f gx = Vector3f(i * h, j * h, k * h);
						double d = tri.closestDistance(gx);
						if (d < sdf->get(Size3(i, j, k))) {
							sdf->get(Size3(i, j, k)) = d;
							closestTri->set(Size3(i, j, k), t);
						}
					}
				}
			}

			// Do intersection counts
			j0 = static_cast<ssize_t>(std::ceil(min3<double>(f1.y(), f2.y(), f3.y())));
			j0 = clamp<ssize_t>(j0 - bandwidth, 0, maxSizeY - 1);
			j1 = static_cast<ssize_t>(std::floor(max3<double>(f1.y(), f2.y(), f3.y())));
			j1 = clamp<ssize_t>(j1 + bandwidth + 1, 0, maxSizeY - 1);
			k0 = static_cast<ssize_t>(std::ceil(min3<double>(f1.z(), f2.z(), f3.z())));
			k0 = clamp<ssize_t>(k0 - bandwidth, 0, maxSizeZ - 1);
			k1 = static_cast<ssize_t>(std::floor(max3<double>(f1.z(), f2.z(), f3.z())));
			k1 = clamp<ssize_t>(k1 + bandwidth + 1, 0, maxSizeZ - 1);

			for (ssize_t k = k0; k <= k1; ++k) {
				for (ssize_t j = j0; j <= j1; ++j) {
					double a, b, c;
					double jD = static_cast<double>(j);
					double kD = static_cast<double>(k);
					if (pointInTriangle2D(
						jD, kD, f1.y(), f1.z(), f2.y(), f2.z(), f3.y(), f3.z(), &a, &b, &c)) {
						// intersection i coordinate
						double fi = a * f1.x() + b * f2.x() + c * f3.x();

						// intersection is in (iInterval - 1, iInterval]
						int iInterval = static_cast<int>(std::ceil(fi));
						if (iInterval < 0) {
							// we enlarge the first interval to include everything
							// to the -x direction
							++intersectionCount ->get(Size3(0, j, k));
						}
						else if (iInterval < static_cast<int>(size.x())) {
							++intersectionCount->get(Size3(iInterval, j, k));
						}
						// we ignore intersections that are beyond the +x side of
						// the grid
					}
				}
			}
		}

		// and now we fill in the rest of the distances with fast sweeping
		for (unsigned int pass = 0; pass < 2; ++pass) {
			std::cout<<"!"<<std::endl;
			sweep(mesh, +1, +1, +1, sdf, closestTri);
			std::cout<<"!!"<<std::endl;
			sweep(mesh, -1, -1, -1, sdf, closestTri);
			std::cout<<"!!!"<<std::endl;
			sweep(mesh, +1, +1, -1, sdf, closestTri);
			std::cout<<"!!!!"<<std::endl;
			sweep(mesh, -1, -1, +1, sdf, closestTri);
			std::cout<<"!!!!!"<<std::endl;
			sweep(mesh, +1, -1, +1, sdf, closestTri);
			std::cout<<"!!!!!!"<<std::endl;
			sweep(mesh, -1, +1, -1, sdf, closestTri);
			std::cout<<"!!!!!!!"<<std::endl;
			sweep(mesh, +1, -1, -1, sdf, closestTri);
			std::cout<<"!!!!!!!!"<<std::endl;
			sweep(mesh, -1, +1, +1, sdf, closestTri);
			std::cout<<"!!!!!!!!!"<<std::endl;
		}

		// then figure out signs (inside/outside) from intersection counts
		for (size_t k = 0; k < size.z(); ++k) {
			for (size_t j = 0; j < size.y(); ++j) {
				unsigned int totalCount = 0U;
				for (size_t i = 0; i < size.x(); ++i) {
					totalCount += intersectionCount->get(Size3(i, j, k));
					// if parity of intersections so far is odd,
					if (totalCount % 2 == 1) {
						// we are inside the mesh
						sdf->get(Size3(i, j, k)) = -sdf->get(Size3(i, j, k));
					}
				}
			}
		}
	}


}
}
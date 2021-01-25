#pragma once

#include "general.h"
#include "geometry/bounding_box3.h"
#include "geometry/ray3.h"
#include <functional>

namespace Omni {
namespace Geometry {

	//! Nearest neighbor query result.
	template <typename T>
	struct NearestNeighborQueryResult3 {
		const T* item = nullptr;
		double distance = std::numeric_limits<double>::max();
	};

	//! Nearest neighbor distance measure function.
	template <typename T>
	using NearestNeighborDistanceFunc3 =
		std::function<double(const T&, const Vector3f&)>;

	//! Abstract base class for 3-D nearest neigbor query engine.
	template <typename T>
	class NearestNeighborQueryEngine3 {
	public:
		//! Returns the nearest neighbor for given point and distance measure
		//! function.
		virtual NearestNeighborQueryResult3<T> nearest(
			const Vector3f& pt,
			const NearestNeighborDistanceFunc3<T>& distanceFunc) const = 0;
	};

}
}
#pragma once
#include "geometry/shape.h"
#include <vector>
#include "geometry/bounding_box3.h"
#include <numeric>
#include "geometry/intersection_query_engine3.h"
#include "geometry/nearest_neighbor_query_engine3.h"

namespace Omni {
namespace Geometry {
	//! \brief Bounding Volume Hierarchy (BVH) in 3D
	//!
	//! This class implements the classic bounding volume hierarchy structure in 3D.
	//! It implements IntersectionQueryEngine3 in order to support box/ray
	//! intersection tests. Also, NearestNeighborQueryEngine3 is implemented to
	//! provide nearest neighbor query.
	//!
	template <typename T>
	class Bvh3 final : public IntersectionQueryEngine3<T>,
		public NearestNeighborQueryEngine3<T> {
	public:
		typedef std::vector<T> ContainerType;
		typedef typename ContainerType::iterator Iterator;
		typedef typename ContainerType::const_iterator ConstIterator;

		//! Default constructor.
		Bvh3();

		//! Builds bounding volume hierarchy.
		void build(const std::vector<T>& items,
			const std::vector<BoundingBox3>& itemsBounds);

		//! Clears all the contents of this instance.
		void clear();

		//! Returns the nearest neighbor for given point and distance measure
		//! function.
		NearestNeighborQueryResult3<T> nearest(
			const Vector3f& pt,
			const NearestNeighborDistanceFunc3<T>& distanceFunc) const override;

		//! Returns true if given \p box intersects with any of the stored items.
		bool intersects(const BoundingBox3& box,
			const BoxIntersectionTestFunc3<T>& testFunc) const override {
			throw std::runtime_error("[intersects] not implemented yet");
			return true;
		}

		//! Returns true if given \p ray intersects with any of the stored items.
		bool intersects(const Ray3& ray,
			const RayIntersectionTestFunc3<T>& testFunc) const override {
			throw std::runtime_error("[intersects] not implemented yet");
			return true;
		}

		//! Invokes \p visitorFunc for every intersecting items.
		void forEachIntersectingItem(
			const BoundingBox3& box, const BoxIntersectionTestFunc3<T>& testFunc,
			const IntersectionVisitorFunc3<T>& visitorFunc) const override {
			throw std::runtime_error("[forEachIntersectingItem] not implemented yet");
		}

		//! Invokes \p visitorFunc for every intersecting items.
		void forEachIntersectingItem(
			const Ray3& ray, const RayIntersectionTestFunc3<T>& testFunc,
			const IntersectionVisitorFunc3<T>& visitorFunc) const override {
			throw std::runtime_error("[forEachIntersectingItem] not implemented yet");
		}

		//! Returns the closest intersection for given \p ray.
		ClosestIntersectionQueryResult3<T> closestIntersection(
			const Ray3& ray,
			const GetRayIntersectionFunc3<T>& testFunc) const override {
			throw std::runtime_error("[closestIntersects] not implemented yet");
			return ClosestIntersectionQueryResult3<T>();
		}

		//! Returns bounding box of every items.
		const BoundingBox3& boundingBox() const;

		//! Returns the begin iterator of the item.
		Iterator begin();

		//! Returns the end iterator of the item.
		Iterator end();

		//! Returns the immutable begin iterator of the item.
		ConstIterator begin() const;

		//! Returns the immutable end iterator of the item.
		ConstIterator end() const;

		//! Returns the number of items.
		size_t numberOfItems() const;

		//! Returns the item at \p i.
		const T& item(size_t i) const;

	private:
		struct Node {
			char flags;
			union {
				size_t child;
				size_t item;
			};
			BoundingBox3 bound;

			Node();
			void initLeaf(size_t it, const BoundingBox3& b);
			void initInternal(uint8_t axis, size_t c, const BoundingBox3& b);
			bool isLeaf() const;
		};

		BoundingBox3 _bound;
		ContainerType _items;
		std::vector<BoundingBox3> _itemBounds;
		std::vector<Node> _nodes;

		size_t build(size_t nodeIndex, size_t* itemIndices, size_t nItems,
			size_t currentDepth);

		size_t qsplit(size_t* itemIndices, size_t numItems, double pivot,
			uint8_t axis);
	};

	template <typename T>
	Bvh3<T>::Node::Node() : flags(0) {
		child = std::numeric_limits<size_t>::max();
	}

	template <typename T>
	void Bvh3<T>::Node::initLeaf(size_t it, const BoundingBox3& b) {
		flags = 3;
		item = it;
		bound = b;
	}

	template <typename T>
	void Bvh3<T>::Node::initInternal(uint8_t axis, size_t c,
		const BoundingBox3& b) {
		flags = axis;
		child = c;
		bound = b;
	}

	template <typename T>
	bool Bvh3<T>::Node::isLeaf() const {
		return flags == 3;
	}

	template <typename T>
	inline NearestNeighborQueryResult3<T> Bvh3<T>::nearest(
		const Vector3f& pt,
		const NearestNeighborDistanceFunc3<T>& distanceFunc) const {
		NearestNeighborQueryResult3<T> best;
		best.distance = std::numeric_limits<double>::max();
		best.item = nullptr;

		// Prepare to traverse BVH
		static const int kMaxTreeDepth = 8 * sizeof(size_t);
		const Node* todo[kMaxTreeDepth];
		size_t todoPos = 0;

		// Traverse BVH nodes
		const Node* node = _nodes.data();
		while (node != nullptr) {
			if (node->isLeaf()) {
				double dist = distanceFunc(_items[node->item], pt);
				if (dist < best.distance) {
					best.distance = dist;
					best.item = &_items[node->item];
				}

				// Grab next node to process from todo stack
				if (todoPos > 0) {
					// Dequeue
					--todoPos;
					node = todo[todoPos];
				}
				else {
					break;
				}
			}
			else {
				const double bestDistSqr = best.distance * best.distance;

				const Node* left = node + 1;
				const Node* right = &_nodes[node->child];

				// If pt is inside the box, then the closestLeft and Right will be
				// identical to pt. This will make distMinLeftSqr and
				// distMinRightSqr zero, meaning that such a box will have higher
				// priority.
				Vector3f closestLeft = left->bound.clamp(pt);
				Vector3f closestRight = right->bound.clamp(pt);

				double distMinLeftSqr = (closestLeft - pt).squaredNorm();
				double distMinRightSqr = (closestRight - pt).squaredNorm();

				bool shouldVisitLeft = distMinLeftSqr < bestDistSqr;
				bool shouldVisitRight = distMinRightSqr < bestDistSqr;

				const Node* firstChild;
				const Node* secondChild;
				if (shouldVisitLeft && shouldVisitRight) {
					if (distMinLeftSqr < distMinRightSqr) {
						firstChild = left;
						secondChild = right;
					}
					else {
						firstChild = right;
						secondChild = left;
					}

					// Enqueue secondChild in todo stack
					todo[todoPos] = secondChild;
					++todoPos;
					node = firstChild;
				}
				else if (shouldVisitLeft) {
					node = left;
				}
				else if (shouldVisitRight) {
					node = right;
				}
				else {
					if (todoPos > 0) {
						// Dequeue
						--todoPos;
						node = todo[todoPos];
					}
					else {
						break;
					}
				}
			}
		}
		return best;
	}

	template <typename T>
	const BoundingBox3& Bvh3<T>::boundingBox() const {
		return _bound;
	}

	template <typename T>
	Bvh3<T>::Bvh3() {}

	template<typename T>
	void Bvh3<T>::build(const std::vector<T>& items, const std::vector<BoundingBox3>& itemsBounds)
	{
		_items = items;
		_itemBounds = itemsBounds;
		if (_items.empty()) {
			return;
		}
		_nodes.clear();
		for (size_t i = 0; i < _items.size(); ++i) {
			_bound.merge(_itemBounds[i]);
		}
		std::vector<size_t> itemIndices(_items.size());
		std::iota(std::begin(itemIndices), std::end(itemIndices), 0);
		build(0, itemIndices.data(), _items.size(), 0);
	}

	template <typename T>
	size_t Bvh3<T>::build(size_t nodeIndex, size_t* itemIndices, size_t nItems,
		size_t currentDepth) {
		// add a node
		_nodes.push_back(Node());

		// initialize leaf node if termination criteria met
		if (nItems == 1) {
			_nodes[nodeIndex].initLeaf(itemIndices[0], _itemBounds[itemIndices[0]]);
			return currentDepth + 1;
		}

		// find the mid-point of the bounding box to use as a qsplit pivot
		BoundingBox3 nodeBound;
		for (size_t i = 0; i < nItems; ++i) {
			nodeBound.merge(_itemBounds[itemIndices[i]]);
		}

		Vector3f d = nodeBound.upperCorner - nodeBound.lowerCorner;

		// choose which axis to split along
		uint8_t axis;
		if (d.x() > d.y() && d.x() > d.z()) {
			axis = 0;
		}
		else {
			axis = (d.y() > d.z()) ? 1 : 2;
		}

		double pivot =
			0.5 * (nodeBound.upperCorner[axis] + nodeBound.lowerCorner[axis]);

		// classify primitives with respect to split
		size_t midPoint = qsplit(itemIndices, nItems, pivot, axis);

		// recursively initialize children _nodes
		size_t d0 = build(nodeIndex + 1, itemIndices, midPoint, currentDepth + 1);
		_nodes[nodeIndex].initInternal(axis, _nodes.size(), nodeBound);
		size_t d1 = build(_nodes[nodeIndex].child, itemIndices + midPoint,
			nItems - midPoint, currentDepth + 1);

		return std::max(d0, d1);
	}

	template <typename T>
	size_t Bvh3<T>::qsplit(size_t* itemIndices, size_t numItems, double pivot,
		uint8_t axis) {
		double centroid;
		size_t ret = 0;
		for (size_t i = 0; i < numItems; ++i) {
			BoundingBox3 b = _itemBounds[itemIndices[i]];
			centroid = 0.5f * (b.lowerCorner[axis] + b.upperCorner[axis]);
			if (centroid < pivot) {
				std::swap(itemIndices[i], itemIndices[ret]);
				ret++;
			}
		}
		if (ret == 0 || ret == numItems) {
			ret = numItems >> 1;
		}
		return ret;
	}

}
}
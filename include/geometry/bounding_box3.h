#pragma once

#include "general.h"
#include "geometry/ray3.h"
#include <limits>

namespace Omni {
namespace Geometry {


struct BoundingBoxRayIntersection3 {
    //! True if the box and ray intersects.
    bool isIntersecting = false;

    //! Distance to the first intersection point.
    FLOAT tNear = std::numeric_limits<FLOAT>::max();

    //! Distance to the second (and the last) intersection point.
    FLOAT tFar = std::numeric_limits<FLOAT>::max();
};


class BoundingBox3{
public:
    //! Lower corner of the bounding box.
    Vector3f lowerCorner;

    //! Upper corner of the bounding box.
    Vector3f upperCorner;

    //! Default constructor.
    BoundingBox3();

    //! Constructs a box that tightly covers two points.
    BoundingBox3(const Vector3f& point1, const Vector3f& point2);

    //! Returns width of the box.
    FLOAT width() const;

    //! Returns height of the box.
    FLOAT height() const;

    //! Returns depth of the box.
    FLOAT depth() const;

    //! Returns length of the box in given axis.
    FLOAT length(size_t axis);

    //! Returns true of this box and other box overlaps.
    bool overlaps(const BoundingBox3& other) const;

    //! Returns true if the input vector is inside of this box.
    bool contains(const Vector3f& point) const;

    //! Returns true if the input ray is intersecting with this box.
    bool intersects(const Ray3& ray) const;

    //! Returns intersection.isIntersecting = true if the input ray is
    //! intersecting with this box. If interesects, intersection.tNear is
    //! assigned with distant to the closest intersecting point, and
    //! intersection.tFar with furthest.
    BoundingBoxRayIntersection3 closestIntersection(
        const Ray3& ray) const;

    //! Returns the mid-point of this box.
    Vector3f midPoint() const;

    //! Returns diagonal length of this box.
    FLOAT diagonalLength() const;

    //! Returns squared diagonal length of this box.
    FLOAT diagonalLengthSquared() const;

    //! Resets this box to initial state (min=infinite, max=-infinite).
    void reset();

    //! Merges this and other point.
    void merge(const Vector3f& point);

    //! Merges this and other box.
    void merge(const BoundingBox3& other);

    //! Expands this box by given delta to all direction.
    //! If the width of the box was x, expand(y) will result a box with
    //! x+y+y width.
    void expand(FLOAT delta);

    //! Returns corner position. Index starts from x-first order.
    Vector3f corner(size_t idx) const;

    //! Returns the clamped point.
    Vector3f clamp(const Vector3f& point) const;

    //! Returns true if the box is empty.
    bool isEmpty() const;
};


}  // namespace Geometry
}  // namespace omni
#pragma once
#include "general.h"
#include "geometry/ray2.h"
#include "geometry/ray3.h"
#include "geometry/transform2.h"
#include "geometry/transform3.h"
#include "geometry/bounding_box3.h"

namespace Omni {
namespace Geometry {


//! Struct that represents ray-surface intersection point.
struct ShapeRayIntersection2 {
    bool isIntersecting = false;
    FLOAT distance = MAX_FLOAT;
    Vector2f point;
    Vector2f normal;
};

struct ShapeRayIntersection3 {
    bool isIntersecting = false;
    FLOAT distance = MAX_FLOAT;
    Vector3f point;
    Vector3f normal;
};

class Shape2 {
public:
    Transform2 transform;
    bool isNormalFlipped = false;

    Shape2(const Transform2& _transform, bool _isNormalFlipped = false) 
        :transform(_transform), isNormalFlipped(_isNormalFlipped) {}
    virtual ~Shape2() {}

    bool intersects(const Ray2& ray) const;
    ShapeRayIntersection2 closestIntersection(const Ray2& ray) const;

    Vector2f closestPoint(const Vector2f& otherPoint) const {
        return transform.toWorld(closestPointLocal(transform.toLocal(otherPoint)));
    }
    
    FLOAT closestDistance(const Vector2f& otherPoint) const {
        return closestDistanceLocal(transform.toLocal(otherPoint));
    }
    
    Vector2f closestNormal(const Vector2f& otherPoint) const {
        auto result = transform.toWorldDirection(
            closestNormalLocal(transform.toLocal(otherPoint)));
        result *= (isNormalFlipped) ? -1.0 : 1.0;
        return result;
    }

protected:
    virtual Vector2f closestPointLocal(const Vector2f& otherPoint) const = 0;
    virtual ShapeRayIntersection2 closestIntersectionLocal(const Ray2& ray) const = 0;
    virtual Vector2f closestNormalLocal(const Vector2f& otherPoint) const = 0;
    virtual bool intersectsLocal(const Ray2& ray) const;
    virtual FLOAT closestDistanceLocal(const Vector2f& otherPoint) const {
        return (otherPoint - (closestPointLocal(otherPoint))).norm();
    }
};


class Shape3 {
public:
    Transform3 transform;
    bool isNormalFlipped = false;

    Shape3(const Transform3& _transform = Transform3(),
            bool _isNormalFlipped = false)
        : transform(_transform), isNormalFlipped(_isNormalFlipped) {}
    virtual ~Shape3() {}

    bool intersects(const Ray3& ray) const;
    ShapeRayIntersection3 closestIntersection(const Ray3& ray) const;
    
    Vector3f closestPoint(const Vector3f& otherPoint) const {
        return transform.toWorld(closestPointLocal(transform.toLocal(otherPoint)));
    }
    
    FLOAT closestDistance(const Vector3f& otherPoint) const {
        return closestDistanceLocal(transform.toLocal(otherPoint));
    }

    Vector3f closestNormal(const Vector3f& otherPoint) const {
        auto result = transform.toWorldDirection(
            closestNormalLocal(transform.toLocal(otherPoint)));
        result *= (isNormalFlipped) ? -1.0 : 1.0;
        return result;
    }

	//! Returns the bounding box of this surface object.
	BoundingBox3 boundingBox() const;
    
protected:
    virtual Vector3f closestPointLocal(const Vector3f& otherPoint) const = 0;
    virtual ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const = 0;
    virtual Vector3f closestNormalLocal(const Vector3f& otherPoint) const = 0;
    virtual bool intersectsLocal(const Ray3& ray) const;
	virtual BoundingBox3 boundingBoxLocal() const;
    
    virtual FLOAT closestDistanceLocal(const Vector3f& otherPoint) const {
        return (otherPoint - (closestPointLocal(otherPoint))).norm();
    }
};

typedef std::shared_ptr<Shape2> Shape2Ptr;
typedef std::shared_ptr<Shape3> Shape3Ptr;

}   // end of namespace Geometry
}   // end of namespace Omni
#pragma once
#include "geometry/shape.h"
#include "util.h"

namespace Omni {
namespace Geometry {

// Rectangle centered at origin
class Rectangle2 : public Shape2 {
public:
    Vector2f wh; // x:width, h:height

    Rectangle2(const Vector2f& wh, 
        Transform2 transform=Transform2(), bool isNormalFlipped=false);
    
    // TODO 加速
    Vector2f closestPointLocal(const Vector2f& otherPoint) const override {
        Vector2f ret;
        ret.x() = fmax(fmin(otherPoint.x(), halfWH.x()), -halfWH.x());
        ret.y() = fmax(fmin(otherPoint.y(), halfWH.y()), -halfWH.y());
        FLOAT signx = sign(ret.x());
        FLOAT signy = sign(ret.y());
        if(fabs(ret.x() - signx * halfWH.x()) < fabs(ret.y() - signy * halfWH.y())) {
            ret.x() = signx * halfWH.x();
        } else {
            ret.y() = signy * halfWH.y();
        }
        return ret;
    }

    FLOAT closestDistanceLocal(const Vector2f& otherPoint) const override {
        return (closestPointLocal(otherPoint) - otherPoint).norm();
    }

    bool intersectsLocal(const Ray2& ray) const override;

    Vector2f closestNormalLocal(const Vector2f& otherPoint) const override {
        Vector2f p = closestPointLocal(otherPoint);
        if(fabs(p.x()) < halfWH.x()) {
            return sign(p.y()) * Vector2f(0, 1);
        } else {
            return sign(p.x()) * Vector2f(1, 0);
        }
    }

    ShapeRayIntersection2 closestIntersectionLocal(
        const Ray2& ray) const override;

public:
    Vector2f halfWH;
};


typedef std::shared_ptr<Rectangle2> Rectangle2Ptr;

} // end of namespace geometry
} // end of namespace Omni
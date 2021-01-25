#pragma once

#include "general.h"
#include "ray2.h"

namespace Omni {
namespace Geometry {

class Transform2 {
public:
    Transform2() 
        :translation(Vector2f(0, 0)), orientation(0)
    {
        cosAngle = std::cos(orientation);
        sinAngle = std::sin(orientation);
    }

    Transform2(const Vector2f& _translation, FLOAT _orientation = 0)
        : translation(_translation), orientation(_orientation) {
        cosAngle = std::cos(orientation);
        sinAngle = std::sin(orientation);
    }

    const Vector2f& getTranslation() const { return translation;}

    void setTranslation(const Vector2f& _translation) {
        translation = _translation;
    }

    FLOAT getOrientation() const {return orientation;}

    void setOrientation(FLOAT _orientation) {
        orientation = _orientation;
        cosAngle = std::cos(orientation);
        sinAngle = std::sin(orientation);
    }

    Vector2f toLocal(const Vector2f& pointInWorld) const {
        Vector2f xmt = pointInWorld - translation;
        return Vector2f( cosAngle * xmt.x() + sinAngle * xmt.y(),
                        -sinAngle * xmt.x() + cosAngle * xmt.y());
    }

    Vector2f toLocalDirection(const Vector2f& dirInWorld) const {
        return Vector2f(
            cosAngle * dirInWorld.x() + sinAngle * dirInWorld.y(),
           -sinAngle * dirInWorld.x() + cosAngle * dirInWorld.y());
    }

    Ray2 toLocal(const Ray2& rayInWorld) const {
        return Ray2(toLocal(rayInWorld.origin),
                        toLocalDirection(rayInWorld.direction));
    }

    Vector2f toWorld(const Vector2f& pointInLocal) const {
        return Vector2f(
            cosAngle * pointInLocal.x() - sinAngle * pointInLocal.y()
                + translation.x(),
            sinAngle * pointInLocal.x() + cosAngle * pointInLocal.y()
                + translation.y());
    }

    Vector2f toWorldDirection(const Vector2f& dirInLocal) const {
        return Vector2f(
            cosAngle * dirInLocal.x() - sinAngle * dirInLocal.y(),
            sinAngle * dirInLocal.x() + cosAngle * dirInLocal.y());
    }

    Ray2 toWorld(const Ray2& rayInLocal) const {
        return Ray2( toWorld(rayInLocal.origin),
                        toWorldDirection(rayInLocal.direction));
    }

private:
    Vector2f translation;
    FLOAT orientation = 0.0;
    FLOAT cosAngle = 1.0;
    FLOAT sinAngle = 0.0;
};


} // end of namespace Geometry
} // end of namespace Omni
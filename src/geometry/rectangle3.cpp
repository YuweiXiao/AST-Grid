#include "geometry/rectangle3.h"
#include "util.h"

using namespace Omni;
using namespace Omni::Geometry;

Rectangle3::Rectangle3(const Vector3f& _lwh, Transform3 transform, bool isNormalFlipped) 
    :Shape3(transform, isNormalFlipped), lwh(_lwh)
{
    half_lwh = _lwh * 0.5;
}

Vector3f Rectangle3::closestPointLocal(const Vector3f& otherPoint) const {
    // inside rectangle
    if(otherPoint.x() > -half_lwh.x() && otherPoint.x() < half_lwh.x() 
        && otherPoint.y() > -half_lwh.y() && otherPoint.y() < half_lwh.y()
        && otherPoint.z() > -half_lwh.z() && otherPoint.z() < half_lwh.z()) {
        FLOAT xDis = fabs(fabs(otherPoint.x()) - half_lwh.x());
        FLOAT yDis = fabs(fabs(otherPoint.y()) - half_lwh.y());
        FLOAT zDis = fabs(fabs(otherPoint.z()) - half_lwh.z());
        if(xDis <= yDis && xDis <= zDis) {
            return Vector3f(half_lwh.x() * sign(otherPoint.x()), otherPoint.y(), otherPoint.z());
        } else if(yDis <= xDis && yDis <= zDis) {
            return Vector3f(otherPoint.x(), half_lwh.y() * sign(otherPoint.y()), otherPoint.z());
        } else {
            return Vector3f(otherPoint.x(), otherPoint.y(), half_lwh.z() * sign(otherPoint.z()));
        }
    } else {
        return Vector3f(fmax(fmin(half_lwh.x(), otherPoint.x()), -half_lwh.x()),
                fmax(fmin(half_lwh.y(), otherPoint.y()), -half_lwh.y()),
                fmax(fmin(half_lwh.z(), otherPoint.z()), -half_lwh.z()));
    }
}

Vector3f Rectangle3::closestNormalLocal(const Vector3f& otherPoint) const {

    if(otherPoint.x() >= -half_lwh.x() && otherPoint.x() <= half_lwh.x() 
        && otherPoint.y() >= -half_lwh.y() && otherPoint.y() <= half_lwh.y()
        && otherPoint.z() >= -half_lwh.z() && otherPoint.z() <= half_lwh.z()) {
        FLOAT xDis = fabs(fabs(otherPoint.x()) - half_lwh.x());
        FLOAT yDis = fabs(fabs(otherPoint.y()) - half_lwh.y());
        FLOAT zDis = fabs(fabs(otherPoint.z()) - half_lwh.z());
        if(xDis <= yDis && xDis <= zDis) {
            return sign(otherPoint.x()) * Vector3f(1, 0, 0);
        } else if(yDis <= xDis && yDis <= zDis) {
            return sign(otherPoint.y()) * Vector3f(0, 1, 0);
        } else {
            return sign(otherPoint.z()) * Vector3f(0, 0, 1);
        }

    } else {
        Vector3f p = closestPointLocal(otherPoint);

        Vector3f closestPointToInputPoint = otherPoint - p;
        Vector3f closestNormal = Vector3f(1, 0, 0);
        FLOAT maxCosineAngle = closestNormal.dot(closestPointToInputPoint);

        for(int i = 0; i < 3; ++i) {
            for(int sign = 0; sign < 2; ++sign) {
                Vector3f tn(0, 0, 0);
                tn(i) = (sign == 0 ? 1 : -1);
                FLOAT cosineAngle = tn.dot(closestPointToInputPoint);

                if (cosineAngle > maxCosineAngle) {
                    closestNormal = tn;
                    maxCosineAngle = cosineAngle;
                }
            }
        }
        return closestNormal;
    }
    

    // if(fabs(p.x()) < half_lwh.x() && fabs(p.y()) < half_lwh.y()) {
    //     return sign(p.z()) * Vector3f(0, 0, 1);
    // } else if(fabs(p.x()) < half_lwh.x() && fabs(p.z()) < half_lwh.z()) {
    //     return sign(p.y()) * Vector3f(0, 1, 0);
    // } else {
    //     return sign(p.x()) * Vector3f(1, 0, 0);
    // }
}

bool Rectangle3::intersectsLocal(const Ray3& ray) const{
    // TODO
    throw std::runtime_error("[Rectangle3::intersectsLocal] not implemented yet");
}

ShapeRayIntersection3 Rectangle3::closestIntersectionLocal(const Ray3& ray) const{
    // TODO
    throw std::runtime_error("[Rectangle3::closestIntersectionLocal] not implemented yet");
}
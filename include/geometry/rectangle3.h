#pragma once
#include "geometry/shape.h"

namespace Omni {
namespace Geometry {

class Rectangle3 : public Shape3 {
public:
    Vector3f lwh; // length:x, width:y, height:z

    Rectangle3(const Vector3f& lwh, 
        Transform3 transform=Transform3(), bool isNormalFlipped=false);
    
protected:
    bool intersectsLocal(const Ray3& ray) const override;
    Vector3f closestPointLocal(const Vector3f& otherPoint) const override;
    ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override;
    Vector3f closestNormalLocal(const Vector3f& otherPoint) const override;

private:

    Vector3f half_lwh;
};

typedef std::shared_ptr<Rectangle3> Rectangle3Ptr;

} // end of namespace Geometry
} // end of namespace Omni
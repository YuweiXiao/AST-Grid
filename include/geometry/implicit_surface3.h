#pragma once
#include "general.h"
#include "geometry/shape.h"

namespace Omni {
namespace Geometry {


class ImplicitSurface3 : public Shape3 {
public:
   ImplicitSurface3(const Transform3& transform = Transform3(), bool isNormalFlipped = false) 
      : Shape3(transform, isNormalFlipped)
   {}
   virtual ~ImplicitSurface3() {}

   FLOAT signedDistance(const Vector3f& otherPoint) const {
      return signedDistanceLocal(transform.toLocal(otherPoint));
   }

protected:
   virtual FLOAT signedDistanceLocal(const Vector3f& otherPoint) const = 0;

private:
   FLOAT closestDistanceLocal(const Vector3f& otherPoint) const override {
      return std::fabs(signedDistanceLocal(otherPoint));
   }
};

typedef std::shared_ptr<ImplicitSurface3> ImplicitSurface3Ptr;


} // end of namespace geometry
} // end of namespace Omni
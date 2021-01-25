#pragma once
#include "general.h"
#include "geometry/shape.h"

namespace Omni {
namespace Geometry {


class ImplicitSurface2 : public Shape2 {
public:
   ImplicitSurface2(const Transform2& transform = Transform2(), bool isNormalFlipped = false) 
      : Shape2(transform, isNormalFlipped)
   {}
   virtual ~ImplicitSurface2() {}

   //! Returns signed distance from the given point \p otherPoint.
   FLOAT signedDistance(const Vector2f& otherPoint) const {
      return signedDistanceLocal(transform.toLocal(otherPoint));
   }

protected:
   //! Returns signed distance from the given point \p otherPoint in local space.
   virtual FLOAT signedDistanceLocal(const Vector2f& otherPoint) const = 0;

private:
   FLOAT closestDistanceLocal(const Vector2f& otherPoint) const override {
      return std::fabs(signedDistanceLocal(otherPoint));
   }
};

typedef std::shared_ptr<ImplicitSurface2> ImplicitSurface2Ptr;


} // end of namespace geometry
} // end of namespace Omni
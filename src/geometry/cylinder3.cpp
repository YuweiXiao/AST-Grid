#include "geometry/cylinder3.h"

namespace Omni {
namespace Geometry {

	bool Cylinder3::intersectsLocal(const Ray3 & ray) const
	{
		throw std::runtime_error("[Cylinder3::intersectionLocal] not implemented yet");
		return false;
	}

	ShapeRayIntersection3 Cylinder3::closestIntersectionLocal(const Ray3 & ray) const
	{
		throw std::runtime_error("[Cylinder3::closestIntersectionLocal] not implemented yet");
		return ShapeRayIntersection3();
	}


} // end of namespace Geometry  
} // end of namespace Omni  
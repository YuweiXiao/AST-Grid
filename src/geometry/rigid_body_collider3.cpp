#include "geometry/rigid_body_collider3.h"

using namespace Omni;
using namespace Omni::Geometry;


RigidBodyCollider3::RigidBodyCollider3(const Shape3Ptr& shape) 
    :linearVelocity(0, 0, 0), angularVelocity(0, 0, 0)
{
    setShape(shape);
}

RigidBodyCollider3::RigidBodyCollider3(
    const Shape3Ptr& shape,
    const Vector3f& _linearVelocity,
    const Vector3f& _angularVelocity)
    : linearVelocity(_linearVelocity), angularVelocity(_angularVelocity) 
{
    setShape(shape);
}

Vector3f RigidBodyCollider3::velocityAt(const Vector3f& point) const {
    Vector3f r = point - shape()->transform.translation();
    return linearVelocity + angularVelocity.cross(r);
}
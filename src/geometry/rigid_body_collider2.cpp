#include "geometry/rigid_body_collider2.h"

using namespace Omni;
using namespace Omni::Geometry;

RigidBodyCollider2::RigidBodyCollider2(const Shape2Ptr& shape) 
    : linearVelocity(0, 0), angularVelocity(0)
{
    setShape(shape);
}

RigidBodyCollider2::RigidBodyCollider2(const Shape2Ptr& shape, const Vector2f& _v, double _aV)
    : linearVelocity(_v), angularVelocity(_aV) 
{
    setShape(shape);
}

Vector2f RigidBodyCollider2::velocityAt(const Vector2f& point) const {
    Vector2f r = point - shape()->transform.getTranslation();
    return linearVelocity + angularVelocity * Vector2f(-r.y(), r.x());
}
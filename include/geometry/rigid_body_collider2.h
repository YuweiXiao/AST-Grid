#pragma once
#include "general.h"
#include "collider2.h"

namespace Omni {
namespace Geometry {


// 2-D rigid body collider.
// Rigid body motion with linear and rotational velocities.
class RigidBodyCollider2 : public Collider2 {
public:
    Vector2f linearVelocity;     // Linear velocity of the rigid body.
    FLOAT angularVelocity = 0.0; // Angular velocity of the rigid body.

    explicit RigidBodyCollider2(const Shape2Ptr& shape);
    RigidBodyCollider2(
        const Shape2Ptr& shape,
        const Vector2f& linearVelocity,
        FLOAT angularVelocity);

    Vector2f velocityAt(const Vector2f& point) const override;
};


typedef std::shared_ptr<RigidBodyCollider2> RigidBodyCollider2Ptr;

} // end of namespace geometry
} // end of namespace Omni
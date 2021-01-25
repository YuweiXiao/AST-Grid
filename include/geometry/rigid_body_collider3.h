#pragma once
#include "general.h"
#include "collider3.h"

namespace Omni {
namespace Geometry {


// 3-D rigid body collider.
// Rigid body motion with linear and rotational velocities.
class RigidBodyCollider3 : public Collider3 {
public:
    Vector3f linearVelocity;     // Linear velocity of the rigid body.
    Vector3f angularVelocity;   // Angular velocity of the rigid body.

    explicit RigidBodyCollider3(const Shape3Ptr& shape);
    RigidBodyCollider3(
        const Shape3Ptr& shape,
        const Vector3f& linearVelocity,
        const Vector3f& angularVelocity);

    Vector3f velocityAt(const Vector3f& point) const override;
};


typedef std::shared_ptr<RigidBodyCollider3> RigidBodyCollider3Ptr;

} // end of namespace geometry
} // end of namespace Omni
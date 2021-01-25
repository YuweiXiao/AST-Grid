#pragma once
#include "general.h"
#include "geometry/shape.h"

namespace Omni {
namespace Geometry {

class Collider3 {
public:
    Collider3() {}
    virtual ~Collider3() {}

    virtual Vector3f velocityAt(const Vector3f& point) const = 0;

    const Shape3Ptr& shape() const {return _shape;}

protected:
    void setShape(const Shape3Ptr& _shape) {
        this->_shape = _shape;
    }

private:
    Shape3Ptr _shape = nullptr;
};

typedef std::shared_ptr<Collider3> Collider3Ptr;

}  // end of namespace geometry
} // end of namespace Omni
#pragma once
#include "general.h"
#include "geometry/shape.h"

namespace Omni {
namespace Geometry {

class Collider2 {
public:
    Collider2() {}
    virtual ~Collider2() {}

    virtual Vector2f velocityAt(const Vector2f& point) const = 0;

    const Shape2Ptr& shape() const {return _shape;}

protected: 
    //! Assigns the shape instance from the subclass.
    void setShape(const Shape2Ptr& _shape) {
        this->_shape = _shape;
    }

private:
    Shape2Ptr _shape = nullptr;
};

typedef std::shared_ptr<Collider2> Collider2Ptr;

}  // end of namespace geometry
} // end of namespace Omni

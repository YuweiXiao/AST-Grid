#include <vector>
#include "geometry/shape.h"

namespace Omni {
namespace Geometry {

//
// 2-D Shape set.
//
class Shape2Set : public Shape2 {
public:
    Shape2Set();
    explicit Shape2Set(const std::vector<Shape2Ptr>& others,
                         const Transform2& transform = Transform2(),
                         bool isNormalFlipped = false);

    size_t numberOfShape() const {return _shapes.size(); }
    const Shape2Ptr& shapeAt(size_t i) const {return _shapes[i]; }
    void addShape(const Shape2Ptr& shape) {_shapes.push_back(shape);}

private:
    std::vector<Shape2Ptr> _shapes;

    Vector2f closestPointLocal(const Vector2f& otherPoint) const override;
    FLOAT closestDistanceLocal(const Vector2f& otherPoint) const override;
    bool intersectsLocal(const Ray2& ray) const override;
    Vector2f closestNormalLocal(const Vector2f& otherPoint) const override;
    ShapeRayIntersection2 closestIntersectionLocal(const Ray2& ray) const override;
};

//! Shared pointer for the Shape2Set type.
typedef std::shared_ptr<Shape2Set> Shape2SetPtr;

} // end of namespace Geometry
} // end of namespace Omni
#pragma once
#include "geometry/collider2.h"
#include "geometry/shape.h"
#include "geometry/shape2set.h"
#include <vector>

namespace Omni {
namespace Geometry {

class Collider2Set : public Collider2 {
public:
    Collider2Set() : Collider2Set(std::vector<Collider2Ptr>()) {}
    explicit Collider2Set(const std::vector<Collider2Ptr>& others) {
        setShape(std::make_shared<Shape2Set>());
        for(auto& collider : others) {
            addCollider(collider);
        }
    }

    Vector2f velocityAt(const Vector2f& point) const override {
        Collider2Ptr c = nullptr;
        FLOAT min_dis = MAX_FLOAT;
        for(int i = 0; i < numberOfColliders(); ++i) {
            FLOAT dis = _colliders[i]->shape()->closestDistance(point);
            if(dis < min_dis) {
                min_dis = dis;
                c = _colliders[i];
            }
        }
        if(c != nullptr) {
            return c->velocityAt(point);
        } else {
            return Vector2f(0, 0);
        }
    }

    void addCollider(const Collider2Ptr& collider) { 
        auto shapeSet = std::dynamic_pointer_cast<Shape2Set>(shape());
        _colliders.push_back(collider); 
        shapeSet->addShape(collider->shape());
    }
    size_t numberOfColliders() const { return _colliders.size();}
    Collider2Ptr collider(size_t i) const { return _colliders[i]; }

private:
    std::vector<Collider2Ptr> _colliders;
};

typedef std::shared_ptr<Collider2Set> Collider2SetPtr;


}  // namespace Geometry
}   // namespace Omni

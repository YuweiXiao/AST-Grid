#pragma once
#ifdef WIN32
#define _ENABLE_EXTENDED_ALIGNED_STORAGE
#endif
#include <vector>
#include "geometry/implicit_surface2.h"
#include "geometry/shape_to_implicit2.h"

namespace Omni {
namespace Geometry {

// 2-D implicit surface set.
class ImplicitSurface2Set final : public ImplicitSurface2 {
public:
    ImplicitSurface2Set() {}
    ImplicitSurface2Set(const std::vector<ImplicitSurface2Ptr>& surfaces,
                        const Transform2& transform = Transform2(),
                        bool isNormalFlipped = false)
        : ImplicitSurface2(transform, isNormalFlipped), _surfaces(surfaces)
    {}
    ImplicitSurface2Set(const std::vector<Shape2Ptr>& surfaces,
                        const Transform2& transform = Transform2(),
                        bool isNormalFlipped = false) 
        : ImplicitSurface2(transform, isNormalFlipped)
    {
        for(const auto& surface : surfaces) {
            addExplicitSurface(surface);
        }
    }

    size_t numberOfSurfaces() const {return _surfaces.size(); }
    const ImplicitSurface2Ptr& surfaceAt(size_t i) const {return _surfaces[i];}

    void addExplicitSurface(const Shape2Ptr& surface) {
        addSurface(std::make_shared<ShapeToImplicit2>(surface));
    }
    void addSurface(const ImplicitSurface2Ptr& surface) {
        _surfaces.push_back(surface);
    }

private:
    std::vector<ImplicitSurface2Ptr> _surfaces;


    Vector2f closestPointLocal(const Vector2f& otherPoint) const override;
    FLOAT closestDistanceLocal(const Vector2f& otherPoint) const override;
    bool intersectsLocal(const Ray2& ray) const override;
    Vector2f closestNormalLocal(const Vector2f& otherPoint) const override;
    ShapeRayIntersection2 closestIntersectionLocal(const Ray2& ray) const override;

    // ImplicitSurface3 implementations TODO
    FLOAT signedDistanceLocal(const Vector2f& otherPoint) const override;
};

} // end of namespace geometry
} // end of namespace Omni
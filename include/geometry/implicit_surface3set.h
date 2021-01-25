#pragma once
#include <vector>
#include "geometry/implicit_surface3.h"
#include "geometry/shape_to_implicit3.h"

namespace Omni {
namespace Geometry {

//!
//! \brief 3-D implicit surface set.
//!
//! This class represents 3-D implicit surface set which extends
//! ImplicitSurface3 by overriding implicit surface-related quries. This is
//! class can hold a collection of other implicit surface instances.
//!
class ImplicitSurface3Set : public ImplicitSurface3 {
public:

    //! Constructs an empty implicit surface set.
    ImplicitSurface3Set() {}

    //! Constructs an implicit surface set using list of other surfaces.
    ImplicitSurface3Set(const std::vector<ImplicitSurface3Ptr>& surfaces, const Transform3& transform = Transform3(), bool isNormalFlipped = false) 
        : ImplicitSurface3(transform, isNormalFlipped), _surfaces(surfaces){}

    //! Constructs an implicit surface set using list of other surfaces.
    ImplicitSurface3Set(const std::vector<Shape3Ptr>& shapes, const Transform3& transform = Transform3(), bool isNormalFlipped = false)
        : ImplicitSurface3(transform, isNormalFlipped)
    {
        for (const auto& shape : shapes) {
            addShape(shape);
        }
    }

    //! Returns the number of implicit surfaces.
    size_t numberOfSurfaces() const {return _surfaces.size();}

    //! Returns the i-th implicit surface.
    const ImplicitSurface3Ptr& surfaceAt(size_t i) const { return _surfaces[i];}

    //! Adds an explicit surface instance.
    void addShape(const Shape3Ptr& shape) { 
        _surfaces.push_back(std::make_shared<ShapeToImplicit3>(shape)); }

    //! Adds an implicit surface instance.
    void addImplicitSurface(const ImplicitSurface3Ptr& surface) {
        _surfaces.push_back(surface);
    }

 private:
    std::vector<ImplicitSurface3Ptr> _surfaces;

    // Surface3 implementations

    Vector3f closestPointLocal(const Vector3f& otherPoint) const override {
        ImplicitSurface3Ptr s = nullptr;
        REAL min_dis = MAX_FLOAT;
        for(int i = 0; i < _surfaces.size(); ++i) {
            FLOAT dis = _surfaces[i]->closestDistance(otherPoint);
            if(min_dis > dis) {
                min_dis = dis;
                s = _surfaces[i];
            }
        }

        if(s != nullptr) {
            return s->closestPoint(otherPoint);
        } else {
            spdlog::warn("[closestPointLocal]::No closed shape");
            return Vector3f::Ones();
        }
    }

    double closestDistanceLocal(const Vector3f& otherPoint) const override {
        REAL min_dis = MAX_FLOAT;
        for(int i = 0; i < _surfaces.size(); ++i) {
            FLOAT dis = _surfaces[i]->closestDistance(otherPoint);
            if(min_dis > dis) {
                min_dis = dis;
            }
        }
        return min_dis;
    }

    bool intersectsLocal(const Ray3& ray) const override {
        throw std::runtime_error("[intersectsLocal] not implementated");
    }

    Vector3f closestNormalLocal(const Vector3f& otherPoint) const override {
        ImplicitSurface3Ptr s = nullptr;
        REAL min_dis = MAX_FLOAT;
        for(int i = 0; i < _surfaces.size(); ++i) {
            FLOAT dis = _surfaces[i]->closestDistance(otherPoint);
            if(min_dis > dis) {
                min_dis = dis;
                s = _surfaces[i];
            }
        }

        if(s != nullptr) {
            return s->closestNormal(otherPoint);
        } else {
            spdlog::warn("[closestPointLocal]::No closed shape");
            return Vector3f::Ones();
        }
    }

    ShapeRayIntersection3 closestIntersectionLocal(const Ray3& ray) const override {
        throw std::runtime_error("[intersectsLocal] not implementated");
    }

    // ImplicitSurface3 implementations

    double signedDistanceLocal(const Vector3f& otherPoint) const override {
        REAL min_dis = MAX_FLOAT;
        for(int i = 0; i < _surfaces.size(); ++i) {
            FLOAT dis = _surfaces[i]->signedDistance(otherPoint);
            if(min_dis > dis) {
                min_dis = dis;
            }
        }
        return min_dis;
    }
};


typedef shared_ptr<ImplicitSurface3Set> ImplicitSurface3SetPtr;

}} // end of namespace Geometry / Omni
#pragma once
#include "general.h"
#include "geometry/quaternion.h"
#include "ray3.h"
#include "geometry/bounding_box3.h"

namespace Omni {
namespace Geometry {


class Transform3 {
 public:
    Transform3() {
        setTranslation(Vector3f(0, 0, 0));
        setOrientation(Quaternion());
    }
    Transform3(const Vector3f& translation, const Quaternion& orientation=Quaternion()) {
        setTranslation(translation);
        setOrientation(orientation);
    }

    const Vector3f& translation() const {return _translation;}
    const Quaternion& orientation() const {return _orientation;}
    
    void setTranslation(const Vector3f& translation) {
        _translation = translation;
    }

    //! Sets the orientation.
    void setOrientation(const Quaternion& orientation) {
        _orientation = orientation;
        _orientationMat3 = orientation.matrix3();
        _inverseOrientationMat3 = orientation.inverse().matrix3();
    }

    Vector3f toLocal(const Vector3f& pointInWorld) const {
        return _inverseOrientationMat3 * (pointInWorld - _translation);
    }
    
    Vector3f toLocalDirection(const Vector3f& dirInWorld) const {
        return _inverseOrientationMat3 * dirInWorld;
    }
    
    Ray3 toLocal(const Ray3& rayInWorld) const {
        return Ray3(toLocal(rayInWorld.origin), 
                    toLocalDirection(rayInWorld.direction));
    }

    Vector3f toWorld(const Vector3f& pointInLocal) const {
        return (_orientationMat3 * pointInLocal) + _translation;
    }

    Vector3f toWorldDirection(const Vector3f& dirInLocal) const {
        return _orientationMat3 * dirInLocal;
    }

    Ray3 toWorld(const Ray3& rayInLocal) const {
        return Ray3(
            toWorld(rayInLocal.origin),
            toWorldDirection(rayInLocal.direction));
    }

	BoundingBox3 toWorld(const BoundingBox3& bboxInLocal) const {
		BoundingBox3 bboxInWorld;
		for (int i = 0; i < 8; ++i) {
			auto cornerInWorld = toWorld(bboxInLocal.corner(i));
			// [WARN] min? max?
			bboxInWorld.lowerCorner
				= Vector3f(std::min(bboxInWorld.lowerCorner.x(), cornerInWorld.x()),
					std::min(bboxInWorld.lowerCorner.y(), cornerInWorld.y()),
					std::min(bboxInWorld.lowerCorner.z(), cornerInWorld.z()));
			bboxInWorld.upperCorner
				= Vector3f(std::max(bboxInWorld.lowerCorner.x(), cornerInWorld.x()),
					std::max(bboxInWorld.lowerCorner.y(), cornerInWorld.y()),
					std::max(bboxInWorld.lowerCorner.z(), cornerInWorld.z()));
		}
		return bboxInWorld;
	}

 private:
    Vector3f _translation;
    Quaternion _orientation;
    Matrix3f _orientationMat3;
    Matrix3f _inverseOrientationMat3;
};


} // end of namespace Geometry
} // end of namespace Omni
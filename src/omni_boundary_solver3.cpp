#include "omni_boundary_solver3.h"
#include "constants.h"
#include "geometry/implicit_surface3.h"
#include "geometry/shape_to_implicit3.h"
#include "geometry/shape.h"
#include "global_benchmark.h"

using namespace Omni;
using namespace Omni::Geometry;

OmniBoundarySolver3::OmniBoundarySolver3(const Size3& res, REAL spacing) 
    : _resolution(res), _gridSpacing(spacing)
{
    _colliderSdf = std::make_shared<CellCenteredScalarGrid3>(res, spacing, 1);
}

void OmniBoundarySolver3::updateCollider(const Collider3Ptr& newCollider) {
    _collider = newCollider;
    if (collider() != nullptr) {
        Shape3Ptr shape = collider()->shape();
        ImplicitSurface3Ptr implicitShape
            = std::dynamic_pointer_cast<ImplicitSurface3>(shape);
        if (implicitShape == nullptr) {
            implicitShape = std::make_shared<ShapeToImplicit3>(shape);
        }

        _colliderSdf->parallelForEach([&](const Size3& idx){
            _colliderSdf->set(idx, implicitShape->signedDistance(_colliderSdf->position(idx)));
        });
    } else {
        _colliderSdf->fill(1);
    }
}
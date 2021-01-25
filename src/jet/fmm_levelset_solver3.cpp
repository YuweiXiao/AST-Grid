// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.


#include "jet/fmm_levelset_solver3.h"
#include <unordered_set>
#include "util.h"
#include "timer.h"

using namespace jet;
using namespace Omni;


// Find geometric solution near the boundary
inline REAL solveQuadNearBoundary(CellCenteredScalarGrid3Ptr& output,
                                    REAL gridSpacing, REAL sign, 
                                    size_t i, size_t j, size_t k) {
    Size3 size = output->resolution();

    bool hasX = false;
    REAL phiX = MAX_FLOAT;

    if (i > 0) {
        if (isInsideSdf(sign * output->get(Size3(i - 1, j, k)))) {
            hasX = true;
            phiX = std::min(phiX, sign * output->get(Size3(i - 1, j, k)));
        }
    }

    if (i + 1 < size.x()) {
        if (isInsideSdf(sign * output->get(Size3(i + 1, j, k)))) {
            hasX = true;
            phiX = std::min(phiX, sign * output->get(Size3(i + 1, j, k)));
        }
    }

    bool hasY = false;
    REAL phiY = MAX_FLOAT;

    if (j > 0) {
        if (isInsideSdf(sign * output->get(Size3(i, j - 1, k)))) {
            hasY = true;
            phiY = std::min(phiY, sign * output->get(Size3(i, j - 1, k)));
        }
    }

    if (j + 1 < size.y()) {
        if (isInsideSdf(sign * output->get(Size3(i, j + 1, k)))) {
            hasY = true;
            phiY = std::min(phiY, sign * output->get(Size3(i, j + 1, k)));
        }
    }

    bool hasZ = false;
    REAL phiZ = MAX_FLOAT;

    if (k > 0) {
        if (isInsideSdf(sign * output->get(Size3(i, j, k - 1)))) {
            hasZ = true;
            phiZ = std::min(phiZ, sign * output->get(Size3(i, j, k - 1)));
        }
    }

    if (k + 1 < size.z()) {
        if (isInsideSdf(sign * output->get(Size3(i, j, k + 1)))) {
            hasZ = true;
            phiZ = std::min(phiZ, sign * output->get(Size3(i, j, k + 1)));
        }
    }

    DEBUG_ONLY(if(!(hasX || hasY || hasZ)) {
        spdlog::info("{}, {}, {}, sign:{}, {}, {}, {}, {}, {}, {}", i, j, k, sign, 
            output->get(i, j, k), output->get(i-1, j, k), output->get(i+1, j, k),
            output->get(i, j-1, k), output->get(i, j+1, k), output->get(i, j, k-1), output->get(i, j, k+1));
        throw std::runtime_error("[fmm: hasX || hasY || hasZ == false]");
    });

    const REAL absCenter = std::fabs(output->get(Size3(i, j, k)));

    REAL distToBndX = gridSpacing * absCenter / (absCenter + std::abs(phiX));

    REAL distToBndY = gridSpacing * absCenter / (absCenter + std::abs(phiY));

    REAL distToBndZ = gridSpacing * absCenter / (absCenter + std::abs(phiZ));

    REAL solution;
    REAL denomSqr = 0.0;

    if (hasX) {
        denomSqr += 1.0 / squared(distToBndX);
    }
    if (hasY) {
        denomSqr += 1.0 / squared(distToBndY);
    }
    if (hasZ) {
        denomSqr += 1.0 / squared(distToBndZ);
    }

    solution = 1.0 / std::sqrt(denomSqr);

    return sign * solution;
}

inline REAL solveQuad(CellCenteredCharGrid3Ptr& markers, CellCenteredScalarGrid3Ptr& output,
                        REAL gridSpacing, REAL sign, const Size3& cell_idx) {
    const Size3& size = output->resolution();
    
    Size3 hasXYZ(0, 0, 0);
    Vector3f phiXYZ(MAX_FLOAT, MAX_FLOAT, MAX_FLOAT);
#if HIGH_ACCURACY
    Vector3f phi2XYZ(MAX_FLOAT, MAX_FLOAT, MAX_FLOAT);
#endif

    for(int i = 0; i < 6; ++i) {
        Size3 tIdx = cell_idx + OCT_NB_DXDYDZ[i];
        if(!Valid(tIdx, size))
            continue;
        if(markers->get(tIdx) == kKnown) {
            int axis = i >> 1;
            hasXYZ[axis] = 1;
            REAL v1 = output->get(tIdx); 
            if(v1 < phiXYZ[axis]) {
                phiXYZ[axis] = v1;
#if HIGH_ACCURACY
                Size3 tIdx2 = cell_idx + 2 * OCT_NB_DXDYDZ[i];
                if(tIdx2.x() >= 0 && tIdx2.y() >= 0 && tIdx2.z() >= 0 && tIdx2.x() < size.x() && tIdx2.y() < size.y() && tIdx2.z() < size.z()
                    && markers->get(tIdx2) == kKnown && output->get(tIdx2) <= v1) {
                    phi2XYZ[axis] = output->get(tIdx2);
                } else {
                    phi2XYZ[axis] = MAX_FLOAT;
                }
#endif
            }
        }
    }

    DEBUG_ONLY(if(hasXYZ.sum() == 0) {
        throw std::runtime_error("[ffm::solveQuad] hasX || hasY || hasZ == false");
    });

    REAL solution = 1000;
    for(int i = 0; i < 3; ++i) {
        if(hasXYZ[i])
            solution = std::min(solution, fabs(phiXYZ[i]));
    }
    solution += gridSpacing;
    solution *= sign;

    // Solve quad
    REAL a = 0.0;
    REAL b = 0.0;
    REAL c = -squared(gridSpacing);

    for(int i = 0; i < 3; ++i) {
        if(hasXYZ[i] == 1) {
#if HIGH_ACCURACY
            if(phi2XYZ[i] < MAX_FLOAT) {
                REAL tp = (1.0 / 3.0) * (4 * phiXYZ[i] - phi2XYZ[i]);
                REAL ta = 9.0 / 4.0;
                a += ta;
                b -= ta * tp;
                c += a * squared(tp);
            } else 
#endif
            {
                a += 1.0;
                b -= phiXYZ[i];
                c += squared(phiXYZ[i]);
            }
        }
    }

    REAL det = b * b - a * c;
    if (det >= 0.0) {
        // spdlog::info("x1:{} x2:{}, sign:{}", (-b + std::sqrt(det)) / a, (-b - std::sqrt(det)) / a, sign);
        solution = (-b + sign * std::sqrt(det)) / a;
    } else {
        // spdlog::error("negative delta in solving quadratic, {}, {}, {}, delta:{}, phi:{}, has:{}, {}", a, b, c, det, phiXYZ.transpose(), hasXYZ.transpose(), sign);
        // throw std::runtime_error();
    }
    
    return solution;
}

inline REAL theta(const REAL phi_1,const REAL phi_2) {
    return phi_1/(phi_1-phi_2);
}

inline Vector3f gradient3(const std::function<REAL(int, int, int)>& data,
    const Size3& ds, REAL gridSpacing, size_t i, size_t j, size_t k) {

    REAL left = data((i > 0) ? i - 1 : i, j, k);
    REAL right = data((i + 1 < ds.x()) ? i + 1 : i, j, k);
    REAL down = data(i, (j > 0) ? j - 1 : j, k);
    REAL up = data(i, (j + 1 < ds.y()) ? j + 1 : j, k);
    REAL back = data(i, j, (k > 0) ? k - 1 : k);
    REAL front = data(i, j, (k + 1 < ds.z()) ? k + 1 : k);

    return 0.5 * Vector3f((right - left) / gridSpacing, 
                        (up - down) / gridSpacing, 
                        (front - back) / gridSpacing);
}


void FmmLevelSetSolver3::reinitialize(CellCenteredScalarGrid3Ptr& inputSdf, REAL maxDistance, 
                    CellCenteredScalarGrid3Ptr& outputSdf) {
    Timer timer;
    const Size3& size = inputSdf->resolution();
    REAL gridSpacing = inputSdf->gridSpacing();
    if(markers == nullptr || markers->resolution() != size) {
        markers = make_shared<CellCenteredCharGrid3>(size, gridSpacing, kUnknown);
    } else {
        markers->fill(kUnknown);
    }
    
    // spdlog::info("\t init marker: {}", timer.durationInSeconds()); timer.reset();   // 10e-4
     auto compare = [](const std::pair<Size3, REAL>& a, const std::pair<Size3, REAL>& b) {
         return a.second > b.second;
     };
    std::priority_queue<std::pair<Size3, REAL>, std::vector<std::pair<Size3, REAL>>, decltype(compare)> trial(compare);

    outputSdf->parallelForEach([&](const Size3& idx){
        outputSdf->set(idx, sign(inputSdf->get(idx))*1000);
    });
    // spdlog::info("\t init output Sdf: {}", timer.durationInSeconds()); timer.reset();   // 10e-3

    std::set<LCOOR_T> intf_cells;
    // std::unordered_set<LCOOR_T> intf_cells;
    outputSdf->forEach([&](const Size3& idx){
        REAL cL = inputSdf->get(idx);
        if(fabs(cL) > maxDistance)
            return;
        int cSign = sign(cL);
        for(int i = 0; i < 6; ++i) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[i];
            if(Valid(tIdx, size) && cSign != sign(inputSdf->get(tIdx))) {
                intf_cells.insert(PACK_COOR3(idx.x(), idx.y(), idx.z()));
                intf_cells.insert(PACK_COOR3(tIdx.x(), tIdx.y(), tIdx.z()));
            }
        }
    });

    // spdlog::info("\t init narrow band: {}", timer.durationInSeconds()); timer.reset();  // 0.1
    for(auto& coor: intf_cells) {
        Size3 cell_idx(UNPACK_X3(coor), UNPACK_Y3(coor), UNPACK_Z3(coor));
        Vector3f correct_phi(MAX_FLOAT, MAX_FLOAT, MAX_FLOAT);
        Size3 correct_axis(0, 0, 0);
        int cSign = sign(inputSdf->get(cell_idx));
        for(int i = 0; i < 6; ++i) {
            Size3 tIdx = cell_idx + OCT_NB_DXDYDZ[i];
            if(!Valid(tIdx, size))
                continue;
            if(intf_cells.find(PACK_COOR3(tIdx.x(), tIdx.y(), tIdx.z())) == intf_cells.end() && markers->get(tIdx) == kUnknown) {
                trial.push(std::make_pair(tIdx, fabs(inputSdf->get(tIdx))));
                //trial.push(fabs(outputSdf->get(tIdx)), tIdx);
                markers->set(tIdx, kTrial);
            } else if( cSign != sign(inputSdf->get(tIdx))) {
                REAL c_phi = jet::theta(inputSdf->get(cell_idx), inputSdf->get(tIdx)) * gridSpacing;
                int axis= i >> 1;
                correct_axis[axis]=1;
                correct_phi[axis]= std::min(correct_phi[axis], c_phi);
            }
        }
        if(correct_axis.sum() > 0){
            REAL hmnc_mean = 0;
            for(int i=0;i<3;i++){
                if(correct_axis[i]==0)
                    continue;
                hmnc_mean += 1.0f/squared(correct_phi[i]);
            }
            hmnc_mean=sqrt(1.0/hmnc_mean);
            outputSdf->set(cell_idx, sign(inputSdf->get(cell_idx)) * hmnc_mean);
        }

        markers->set(cell_idx, kKnown);
    }
    spdlog::info("\t\t\t first ring of narrow band init: {}", timer.durationInSeconds()); timer.reset();    // 0.03

    // Propagate
    int count = 0;
    while (!trial.empty()) {
        // Timer tt;
		Size3 idx = trial.top().first;
        trial.pop();
        if(markers->get(idx) == kKnown) {
            continue;
        }
        ++count;
        markers->set(idx, kKnown);
        outputSdf->get(idx) = solveQuad(markers, outputSdf, gridSpacing, sign(inputSdf->get(idx)), idx);
        DEBUG_ONLY(if(sign(outputSdf->get(idx)) != sign(inputSdf->get(idx)) && inputSdf->get(idx) < 900) {  // sign different & uninitialized
            spdlog::info("old:{} new:{}", inputSdf->get(idx), outputSdf->get(idx));
            throw std::runtime_error("different sign");
        });

        if (fabs(outputSdf->get(idx)) > maxDistance) {
            break;
        }
        // tt.reset();
        for(int m = 0; m < 6; ++m) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[m];
            if(!Valid(tIdx, size))
                continue;
            if(markers->get(tIdx) != kKnown) {
                // tt.reset();
                // a1 += tt.durationInSeconds(); tt.reset();
                REAL p = solveQuad(markers, outputSdf, gridSpacing, sign(inputSdf->get(tIdx)), tIdx);
                if(markers->get(tIdx) == kUnknown//) {
                || fabs(p) < fabs(outputSdf->get(tIdx))) {
                    markers->set(tIdx, kTrial);
                    outputSdf->set(tIdx, p);
                    trial.push(std::make_pair(tIdx, fabs(p)));
                    // trial.push(fabs(p), tIdx);
                }
                // a2 += tt.durationInSeconds();
            }
        }
        // a2 += tt.durationInSeconds();
    }

    spdlog::info("\t\t\t propagation: {}, count: {}", timer.durationInSeconds(), count);// 0.4
}

void FmmLevelSetSolver3::extrapolate(const std::function<REAL(int, int, int)>& input,
    const std::function<REAL(int, int, int)>& sdf,
    REAL gridSpacing,
    const Size3& size,
    REAL maxDistance,
    std::function<REAL&(int, int, int)> output){
    
    Vector3f invGridSpacing = Vector3f(1.0 / gridSpacing, 1.0 / gridSpacing, 1.0 / gridSpacing);

    // Build markers
    CellCenteredFlagGrid3 markers(size, gridSpacing, kUnknown);
    markers.parallelForEach([&](const Size3& idx) {
        if (isInsideSdf(sdf(idx.x(), idx.y(), idx.z()))) {
            markers.set(idx, kKnown);
        }
        output(idx.x(), idx.y(), idx.z()) = input(idx.x(), idx.y(), idx.z());
    });

    auto compare = [&](const Size3& a, const Size3& b) {
        return sdf(a.x(), a.y(), a.z()) > sdf(b.x(), b.y(), b.z());
    };

    // Enqueue initial candidates
    std::priority_queue<Size3, std::vector<Size3>, decltype(compare)> trial(compare);
    
    markers.forEach([&](const Size3& idx) {
        if (markers.get(idx) == kKnown) {
            return;
        }
        for(int k = 0; k < 6; ++k) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[k];
            if(tIdx.x() >= 0 && tIdx.y() >=0 && tIdx.z() >= 0 && tIdx.x() < size.x()
                && tIdx.y() < size.y() && tIdx.z() < size.z()) {
                if(markers.get(tIdx) == kKnown) {
                    trial.push(idx);
                    markers.set(idx, kTrial);
                    return;
                }
            }
        }
    });

    // Propagate
    while (!trial.empty()) {
        Size3 idx = trial.top();
        trial.pop();

        if (sdf(idx.x(), idx.y(), idx.z()) > maxDistance) {
            break;
        }

        Vector3f grad = gradient3(sdf, size, gridSpacing, idx.x(), idx.y(), idx.z()).normalized();

        REAL sum = 0.0;
        REAL count = 0.0;

        for(int k = 0; k < 6; ++k) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[k];
            int axis = int(k / 3);
            if(tIdx.x() >= 0 && tIdx.y() >=0 && tIdx.z() >= 0 && tIdx.x() < size.x()
                && tIdx.y() < size.y() && tIdx.z() < size.z()) {
                if(markers.get(tIdx) == kKnown) {
                    REAL weight = std::max(grad[axis], 0.0) * invGridSpacing[axis];

                    // If gradient is zero, then just assign 1 to weight
                    if (weight < EPSILON) {
                        weight = 1.0;
                    }

                    sum += weight * output(tIdx.x(), tIdx.y(), tIdx.z());
                    count += weight;    
                } else if(markers.get(tIdx) == kUnknown) {
                    markers.set(tIdx, kTrial);
                    trial.push(tIdx);
                }
            }
        }

        if(count == 0.0) {
            spdlog::info("count is zero. {}", idx.transpose());
            throw std::runtime_error("count is zero.");
        }

        output(idx.x(), idx.y(), idx.z()) = sum / count;
        markers.set(idx, kKnown);
    }
}

// void FmmLevelSetSolver3::reinitializeFE(CellCenteredScalarGrid3Ptr& inputSdf, REAL maxDistance, CellCenteredScalarGrid3Ptr& outputSdf) {

//     Size3 size = inputSdf->resolution();
//     REAL gridSpacing = inputSdf->gridSpacing();
//     REAL invGridSpacing = 1.0 / gridSpacing;
//     REAL invGridSpacingSqr = invGridSpacing * invGridSpacing;
//     CellCenteredScalarGrid3 markers(size, gridSpacing, kUnknown);

//     outputSdf->fill(*inputSdf);
//     // outputSdf->forEach([&](size_t i, size_t j, size_t k) {
//     //     outputSdf->set(i, j, k, inputSdf->get(i, j, k));
//     // });

//     // Solve geometrically near the boundary
//     markers.forEach([&](size_t i, size_t j, size_t k) {
//         if (isInsideSdf(outputSdf->get(i, j, k)) &&
//             ((i > 0 && !isInsideSdf(outputSdf->get(i - 1, j, k))) ||
//              (i + 1 < size.x() && !isInsideSdf(outputSdf->get(i + 1, j, k))) ||
//              (j > 0 && !isInsideSdf(outputSdf->get(i, j - 1, k))) ||
//              (j + 1 < size.y() && !isInsideSdf(outputSdf->get(i, j + 1, k))) ||
//              (k > 0 && !isInsideSdf(outputSdf->get(i, j, k - 1))) ||
//              (k + 1 < size.z() && !isInsideSdf(outputSdf->get(i, j, k + 1))))) {
//             outputSdf->get(i, j, k) = solveQuadNearBoundary(
//                 inputSdf, gridSpacing, -1.0, i, j, k);
//         }
//     });
//     markers.forEach([&](size_t i, size_t j, size_t k) {
//         if (!isInsideSdf(outputSdf->get(i, j, k)) &&
//             ((i > 0 && isInsideSdf(outputSdf->get(i - 1, j, k))) ||
//              (i + 1 < size.x() && isInsideSdf(outputSdf->get(i + 1, j, k))) ||
//              (j > 0 && isInsideSdf(outputSdf->get(i, j - 1, k))) ||
//              (j + 1 < size.y() && isInsideSdf(outputSdf->get(i, j + 1, k))) ||
//              (k > 0 && isInsideSdf(outputSdf->get(i, j, k - 1))) ||
//              (k + 1 < size.z() && isInsideSdf(outputSdf->get(i, j, k + 1))))) {
//             outputSdf->get(i, j, k) = solveQuadNearBoundary(
//                 inputSdf, gridSpacing, 1.0, i, j, k);
//         }
//     });

//     for (int sign = 0; sign < 2; ++sign) {
//         // Build markers
//         markers.parallelForEach([&](size_t i, size_t j, size_t k) {
//             if (isInsideSdf(outputSdf->get(i, j, k))) {
//                 markers.get(i, j, k) = kKnown;
//             } else {
//                 markers.get(i, j, k) = kUnknown;
//             }
//         });

//         auto compare = [&](const Size3& a, const Size3& b) {
//             return outputSdf->get(a.x(), a.y(), a.z()) > outputSdf->get(b.x(), b.y(), b.z());
//         };

//         // Enqueue initial candidates
//         std::priority_queue<Size3, std::vector<Size3>, decltype(compare)> trial(compare);
//         markers.forEach([&](size_t i, size_t j, size_t k) {
//             if (markers.get(i, j, k) != kKnown &&
//                 ((i > 0 && markers.get(i - 1, j, k) == kKnown) ||
//                  (i + 1 < size.x() && markers.get(i + 1, j, k) == kKnown) ||
//                  (j > 0 && markers.get(i, j - 1, k) == kKnown) ||
//                  (j + 1 < size.y() && markers.get(i, j + 1, k) == kKnown) ||
//                  (k > 0 && markers.get(i, j, k - 1) == kKnown) ||
//                  (k + 1 < size.z() && markers.get(i, j, k + 1) == kKnown))) {
//                 trial.push(Size3(i, j, k));
//                 markers.set(i, j, k, kTrial);
//             }
//         });

//         // Propagate
//         while (!trial.empty()) {
//             Size3 idx = trial.top();
//             trial.pop();

//             size_t i = idx.x();
//             size_t j = idx.y();
//             size_t k = idx.z();

//             markers.get(i, j, k) = kKnown;
//             // spdlog::info("{}", spdlog::info("{}", outputSdf->get(i, j, k)))
//             outputSdf->get(i, j, k) = solveQuad(markers, outputSdf, gridSpacing,
//                                         invGridSpacingSqr, i, j, k);

//             if (outputSdf->get(i, j, k) > maxDistance) {
//                 break;
//             }

//             if (i > 0) {
//                 if (markers.get(i - 1, j, k) == kUnknown) {
//                     markers.get(i - 1, j, k) = kTrial;
//                     outputSdf->get(i - 1, j, k) =
//                         solveQuad(markers, outputSdf, gridSpacing,
//                                   invGridSpacingSqr, i - 1, j, k);
//                     trial.push(Size3(i - 1, j, k));
//                 }
//             }

//             if (i + 1 < size.x()) {
//                 if (markers.get(i + 1, j, k) == kUnknown) {
//                     markers.get(i + 1, j, k) = kTrial;
//                     outputSdf->get(i + 1, j, k) =
//                         solveQuad(markers, outputSdf, gridSpacing,
//                                   invGridSpacingSqr, i + 1, j, k);
//                     trial.push(Size3(i + 1, j, k));
//                 }
//             }

//             if (j > 0) {
//                 if (markers.get(i, j - 1, k) == kUnknown) {
//                     markers.get(i, j - 1, k) = kTrial;
//                     outputSdf->get(i, j - 1, k) =
//                         solveQuad(markers, outputSdf, gridSpacing,
//                                   invGridSpacingSqr, i, j - 1, k);
//                     trial.push(Size3(i, j - 1, k));
//                 }
//             }

//             if (j + 1 < size.y()) {
//                 if (markers.get(i, j + 1, k) == kUnknown) {
//                     markers.get(i, j + 1, k) = kTrial;
//                     outputSdf->get(i, j + 1, k) =
//                         solveQuad(markers, outputSdf, gridSpacing,
//                                   invGridSpacingSqr, i, j + 1, k);
//                     trial.push(Size3(i, j + 1, k));
//                 }
//             }

//             if (k > 0) {
//                 if (markers.get(i, j, k - 1) == kUnknown) {
//                     markers.get(i, j, k - 1) = kTrial;
//                     outputSdf->get(i, j, k - 1) =
//                         solveQuad(markers, outputSdf, gridSpacing,
//                                   invGridSpacingSqr, i, j, k - 1);
//                     trial.push(Size3(i, j, k - 1));
//                 }
//             }

//             if (k + 1 < size.z()) {
//                 if (markers.get(i, j, k + 1) == kUnknown) {
//                     markers.get(i, j, k + 1) = kTrial;
//                     outputSdf->get(i, j, k + 1) =
//                         solveQuad(markers, outputSdf, gridSpacing,
//                                   invGridSpacingSqr, i, j, k + 1);
//                     trial.push(Size3(i, j, k + 1));
//                 }
//             }
//         }

//         // Flip the sign
//         markers.parallelForEach([&](size_t i, size_t j, size_t k) {
//             outputSdf->get(i, j, k) = -outputSdf->get(i, j, k);
//         });
//     }
// }

// void FmmLevelSetSolver3::extrapolate(const ScalarGrid3& input,
//                                      const ScalarField3& sdf,
//                                      REAL maxDistance, ScalarGrid3* output) {
//     JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));

//     Array3<REAL> sdfGrid(input.dataSize());
//     auto pos = input.dataPosition();
//     sdfGrid.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
//         sdfGrid(i, j, k) = sdf.sample(pos(i, j, k));
//     });

//     extrapolate(input.constDataAccessor(), sdfGrid.constAccessor(),
//                 input.gridSpacing(), maxDistance, output->dataAccessor());
// }

// void FmmLevelSetSolver3::extrapolate(const CollocatedVectorGrid3& input,
//                                      const ScalarField3& sdf,
//                                      REAL maxDistance,
//                                      CollocatedVectorGrid3* output) {
//     JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));

//     Array3<REAL> sdfGrid(input.dataSize());
//     auto pos = input.dataPosition();
//     sdfGrid.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
//         sdfGrid(i, j, k) = sdf.sample(pos(i, j, k));
//     });

//     const Vector3f gridSpacing = input.gridSpacing();

//     Array3<REAL> u(input.dataSize());
//     Array3<REAL> u0(input.dataSize());
//     Array3<REAL> v(input.dataSize());
//     Array3<REAL> v0(input.dataSize());
//     Array3<REAL> w(input.dataSize());
//     Array3<REAL> w0(input.dataSize());

//     input.parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
//         u(i, j, k) = input(i, j, k).x;
//         v(i, j, k) = input(i, j, k).y;
//         w(i, j, k) = input(i, j, k).z;
//     });

//     extrapolate(u, sdfGrid.constAccessor(), gridSpacing, maxDistance, u0);

//     extrapolate(v, sdfGrid.constAccessor(), gridSpacing, maxDistance, v0);

//     extrapolate(w, sdfGrid.constAccessor(), gridSpacing, maxDistance, w0);

//     output->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
//         (*output)(i, j, k).x = u(i, j, k);
//         (*output)(i, j, k).y = v(i, j, k);
//         (*output)(i, j, k).z = w(i, j, k);
//     });
// }
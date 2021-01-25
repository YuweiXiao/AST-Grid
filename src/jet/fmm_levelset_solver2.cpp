// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "jet/fmm_levelset_solver2.h"

#include <algorithm>
#include <vector>
#include <queue>

using namespace jet;

static const char kUnknown = 0;
static const char kKnown = 1;
static const char kTrial = 2;

// Find geometric solution near the boundary
inline double solveQuadNearBoundary( CellCenteredScalarGrid2Ptr& output, REAL gridSpacing,
    FLOAT sign, size_t i, size_t j) {

    Size2 size = output->resolution();

    bool hasX = false;
    double phiX = MAX_FLOAT;

    if (i > 0 && isInsideSdf(sign * output->get(Size2(i - 1, j)))) {
        hasX = true;
        phiX = std::min(phiX, sign * output->get(Size2(i - 1, j)));
    }

    if (i < size.x()-1 && isInsideSdf(sign * output->get(Size2(i + 1, j)))) {
        hasX = true;
        phiX = std::min(phiX, sign * output->get(Size2(i + 1, j)));
    }

    bool hasY = false;
    double phiY = MAX_FLOAT;

    if (j > 0 && isInsideSdf(sign * output->get(Size2(i, j - 1)))) {
        hasY = true;
        phiY = std::min(phiY, sign * output->get(Size2(i, j - 1)));
    }

    if (j < size.y()-1 && isInsideSdf(sign * output->get(Size2(i, j + 1)))) {
        hasY = true;
        phiY = std::min(phiY, sign * output->get(Size2(i, j + 1)));
    }

    if(!(hasX || hasY)) {
        throw std::runtime_error("hasX , hasY is false");
    }

    double distToBndX
        = gridSpacing * std::abs(output->get(Size2(i, j)))
        / (std::abs(output->get(Size2(i, j))) + std::abs(phiX));

    double distToBndY
        = gridSpacing * std::abs(output->get(Size2(i, j)))
        / (std::abs(output->get(Size2(i, j))) + std::abs(phiY));

    double solution;
    double denomSqr = 0.0;

    if (hasX) {
        denomSqr += 1.0 / squared(distToBndX);
    }
    if (hasY) {
        denomSqr += 1.0 / squared(distToBndY);
    }

    solution = 1.0 / std::sqrt(denomSqr);

    return sign * solution;
}

inline double solveQuad(
    CellCenteredFlagGrid2& markers, CellCenteredScalarGrid2Ptr& output,
    REAL gridSpacing, REAL invGridSpacingSqr, size_t i, size_t j) {
    Size2 size = output->resolution();

    bool hasX = false;
    double phiX = MAX_FLOAT;

    if (i > 0 && markers.get(Size2(i - 1, j)) == kKnown) {
        hasX = true;
        phiX = std::min(phiX, output->get(Size2(i - 1, j)));
    }

    if (i < size.x()-1 && markers.get(Size2(i + 1, j)) == kKnown) {
        hasX = true;
        phiX = std::min(phiX, output->get(Size2(i + 1, j)));
    }

    bool hasY = false;
    double phiY = MAX_FLOAT;

    if (j>0 && markers.get(Size2(i, j - 1)) == kKnown) {
        hasY = true;
        phiY = std::min(phiY, output->get(Size2(i, j - 1)));
    }

    if (j < size.y()-1 && markers.get(Size2(i, j + 1)) == kKnown) {
        hasY = true;
        phiY = std::min(phiY, output->get(Size2(i, j + 1)));
    }

    if(!(hasX || hasY)) {
        throw std::runtime_error("hasX || hasY error");
    }

    double solution = 0.0;

    // Initial guess
    if (hasX) {
        solution = phiX + gridSpacing;
    }
    if (hasY) {
        solution = std::max(solution, phiY + gridSpacing);
    }

    // Solve quad
    double a = 0.0;
    double b = 0.0;
    double c = -1.0;

    if (hasX) {
        a += invGridSpacingSqr;
        b -= phiX * invGridSpacingSqr;
        c += squared(phiX) * invGridSpacingSqr;
    }
    if (hasY) {
        a += invGridSpacingSqr;
        b -= phiY * invGridSpacingSqr;
        c += squared(phiY) * invGridSpacingSqr;
    }

    double det = b * b - a * c;

    if (det > 0.0) {
        solution = (-b + std::sqrt(det)) / a;
    }

    return solution;
}

FmmLevelSetSolver2::FmmLevelSetSolver2() {
}

void FmmLevelSetSolver2::reinitialize(
    const CellCenteredScalarGrid2Ptr& inputSdf,
    FLOAT maxDistance,
    CellCenteredScalarGrid2Ptr& outputSdf) {

    Size2 size = inputSdf->resolution();
    REAL gridSpacing = inputSdf->gridSpacing();
    REAL invGridSpacing = 1.0 / gridSpacing;
    REAL invGridSpacingSqr = invGridSpacing * invGridSpacing;
    CellCenteredFlagGrid2 markers(size, gridSpacing, kUnknown);

    throw std::runtime_error("FIX ME");
#if 0
    markers.forEach([&](size_t i, size_t j) {
        outputSdf->get(i, j) = inputSdf->get(i, j);
    });

    // Solve geometrically near the boundary
    markers.forEach([&](size_t i, size_t j) {
        if (!isInsideSdf(outputSdf->get(i, j))
            && ((i > 0 && isInsideSdf(outputSdf->get(i - 1, j)))
             || (i < size.x()-1 && isInsideSdf(outputSdf->get(i + 1, j)))
             || (j > 0 && isInsideSdf(outputSdf->get(i, j - 1)))
             || (j < size.y()-1 && isInsideSdf(outputSdf->get(i, j + 1))))) {
            outputSdf->get(i, j) = solveQuadNearBoundary(outputSdf, gridSpacing, 1.0, i, j);
        } else if (isInsideSdf(outputSdf->get(i, j))
            && ((i > 0 && !isInsideSdf(outputSdf->get(i - 1, j)))
             || (i < size.x()-1 && !isInsideSdf(outputSdf->get(i + 1, j)))
             || (j > 0 && !isInsideSdf(outputSdf->get(i, j - 1)))
             || (j < size.y()-1 && !isInsideSdf(outputSdf->get(i, j + 1))))) {
            outputSdf->get(i, j) = solveQuadNearBoundary(outputSdf, gridSpacing, -1.0, i, j);
        }
    });

    for (int sign = 0; sign < 2; ++sign) {
        // Build markers
        markers.parallelForEach([&](size_t i, size_t j) {
            if (isInsideSdf(outputSdf->get(i, j))) {
                markers.get(i, j) = kKnown;
            } else {
                markers.get(i, j) = kUnknown;
            }
        });

        auto compare = [&](const Size2& a, const Size2& b) {
            return outputSdf->get(a.x(), a.y()) > outputSdf->get(b.x(), b.y());
        };

        // Enqueue initial candidates
        std::priority_queue<Size2, std::vector<Size2>, decltype(compare)> trial(compare);
        markers.forEach([&](size_t i, size_t j) {
            if (markers.get(i, j) != kKnown
                && ((i > 0 && markers.get(i - 1, j) == kKnown) || (i < size.x()-1 && markers.get(i + 1, j) == kKnown)
                 || (j > 0 && markers.get(i, j - 1) == kKnown) || (j < size.y()-1 && markers.get(i, j + 1) == kKnown))) {
                trial.push(Size2(i, j));
                markers.get(i, j) = kTrial;
            }
        });

        // Propagate
        while (!trial.empty()) {
            Size2 idx = trial.top();
            trial.pop();

            size_t i = idx.x();
            size_t j = idx.y();
            if(markers.get(i, j) == kKnown) {
                continue;
            }

            markers.get(i, j) = kKnown;
            outputSdf->get(i, j) = solveQuad(
                markers, outputSdf, gridSpacing, invGridSpacingSqr, i, j);

            if (outputSdf->get(i, j) > maxDistance) {
                break;
            }

            if (i > 0 && markers.get(i - 1, j) == kUnknown) {
                markers.get(i - 1, j) = kTrial;
                outputSdf->get(i - 1, j) = solveQuad(markers, outputSdf, gridSpacing,
                    invGridSpacingSqr, i - 1, j);
                trial.push(Size2(i - 1, j));
            }

            if (i < size.x()-1 && markers.get(i + 1, j) == kUnknown) {
                markers.get(i + 1, j) = kTrial;
                outputSdf->get(i + 1, j) = solveQuad(markers, outputSdf, gridSpacing,
                    invGridSpacingSqr, i + 1, j);
                trial.push(Size2(i + 1, j));
            }

            if (j > 0 && markers.get(i, j - 1) == kUnknown) {
                markers.get(i, j - 1) = kTrial;
                outputSdf->get(i, j - 1) = solveQuad(markers, outputSdf, gridSpacing,
                    invGridSpacingSqr, i, j - 1);
                trial.push(Size2(i, j - 1));
            }

            if (j < size.y()-1 && markers.get(i, j + 1) == kUnknown) {
                markers.get(i, j + 1) = kTrial;
                outputSdf->get(i, j + 1) = solveQuad(markers, outputSdf, gridSpacing,
                    invGridSpacingSqr,i, j + 1);
                trial.push(Size2(i, j + 1));
            }
        }

        // Flip the sign
        markers.parallelForEach([&](size_t i, size_t j) {
            outputSdf->get(i, j) = -outputSdf->get(i, j);
        });
    }
#endif
}
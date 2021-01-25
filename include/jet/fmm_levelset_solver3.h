// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_FMM_LEVEL_SET_SOLVER3_H_
#define INCLUDE_JET_FMM_LEVEL_SET_SOLVER3_H_
#include <algorithm>
#include <queue>
#include <vector>
#include <set>
#include "general.h"
#include "octree_omni_grid3.h"
#include "util.h"
#include "timer.h"
#define BENCHMARK(a) ;

namespace jet {

//!
//! \brief Three-dimensional fast marching method (FMM) implementation.
//!
//! This class implements 3-D FMM. First-order upwind-style differencing is used
//! to solve the PDE.
//!
//! \see https://math.berkeley.edu/~sethian/2006/Explanations/fast_marching_explain.html
//! \see Sethian, James A. "A fast marching level set method for monotonically
//!     advancing fronts." Proceedings of the National Academy of Sciences 93.4
//!     (1996): 1591-1595.
//!

using namespace Omni;
using std::pair;
using std::make_pair;
using std::tuple;
using std::make_tuple;

static const char kUnknown = 0;
static const char kKnown = 1;
static const char kTrial = 2;
static const char kDontTouch = 3;
static const char kGhost = 4;

class FmmLevelSetSolver3 {
public:

    //! Default constructor.
    FmmLevelSetSolver3() 
        : markers(nullptr)
    {}

    //!
    //! Reinitializes given scalar field to signed-distance field.
    //!
    //! \param inputSdf Input signed-distance field which can be distorted.
    //! \param maxDistance Max range of reinitialization.
    //! \param outputSdf Output signed-distance field.
    //!
    void reinitialize(CellCenteredScalarGrid3Ptr& inputSdf, REAL maxDistance, CellCenteredScalarGrid3Ptr& outputSdf);

    //!
    //! Extrapolates given scalar field from negative to positive SDF region.
    //!
    //! \param input Input scalar field to be extrapolated.
    //! \param sdf Reference signed-distance field.
    //! \param maxDistance Max range of extrapolation.
    //! \param output Output scalar field.
    //!
    
   template<typename VT, typename SDF>
   void extrapolate(const std::shared_ptr<OctreeOmniGrid3<VT>>& input, const SDF& levelSet, REAL maxDistanceFactor, std::shared_ptr<OctreeOmniGrid3<VT>>& output,
      OctreeOmniCharGrid3Ptr marker, const std::function<bool(REAL dx, REAL v)>& isValidSdf);

 private:
    void extrapolate(
        const std::function<REAL(int, int, int)>& input, const std::function<REAL(int, int, int)>& sdf,
        REAL gridSpacing, const Size3& size, REAL maxDistance, std::function<REAL&(int, int, int)> output);
    
    CellCenteredCharGrid3Ptr markers;
};

// input & output should have same tiltE
template <typename VT, typename SDF>
void FmmLevelSetSolver3::extrapolate(const std::shared_ptr<OctreeOmniGrid3<VT>>& input, const SDF& Sdf, 
    REAL maxDistanceFactor, std::shared_ptr<OctreeOmniGrid3<VT>>& output, 
    OctreeOmniCharGrid3Ptr cacheMarker, const std::function<bool(REAL dx, REAL v)>& isValidSdf) {
    BENCHMARK_SCOPED_TIMER_SECTION t("velocity extrapolate");
    auto& tiltE = input->getTiltEGridPtr();
    auto& layout = input->getLayout();
    BENCHMARK(Timer tt;)
    OctreeOmniCharGrid3Ptr marker = make_shared<OctreeOmniCharGrid3>(input->resolution(), input->gridSpacing(), layout, tiltE, kUnknown);
    if(cacheMarker == nullptr) {
        auto& markerOctagonGrid = marker->getOctagonGrid();
        ThreadPool::parallelIterateGrid(input->getOctagonGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
            if(markerOctagonGrid->isValid(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx())) {
                Vector3f pos = marker->positionOctagon(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
                if(isValidSdf(layout->levelGridSpacing(iter.level) * SQRT3 * 0.5, Sdf->sample(pos))) {
                    marker->setOctagon(iter.levelCoor(), kKnown);
                }
            } else {
                marker->setOctagon(iter.levelCoor(), kDontTouch);
            }

            // auto pos = marker->positionTilt(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
            // if(tiltE->get(iter.levelCoor()).is_closed) {
            //     marker->setTilt(iter.levelCoor(), kDontTouch);
            // } else if(isValidSdf(tiltE->get(iter.levelCoor()).e, Sdf->sample(pos))) {
            //     marker->setTilt(iter.levelCoor(), kKnown);
            // }
        });
    } else {
        // clone caced marker and invalidate closed tilt cell
        marker->getOctagonGrid()->fill(cacheMarker->getOctagonGrid());
        marker->getTiltGrid()->fill(cacheMarker->getTiltGrid());
        ThreadPool::parallelIterateGrid(tiltE->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
            if(tiltE->get(iter.levelCoor()).is_closed)
                marker->setTilt(iter.levelCoor(), kDontTouch);
        });
    }    
    BENCHMARK(spdlog::info("marker init: {}", tt.durationInSeconds()); tt.reset();)
    ThreadPool::parallelIterateGrid(input->getTiltGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter){
        auto pos = marker->positionTilt(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx());
        if(tiltE->get(iter.levelCoor()).is_closed) {
            marker->setTilt(iter.levelCoor(), kDontTouch);
        } else if(isValidSdf(tiltE->get(iter.levelCoor()).e, Sdf->sample(pos))) {
            marker->setTilt(iter.levelCoor(), kKnown);
        }
    });

    ThreadPool::parallelIterateGrid(marker->getOctagonGrid()->getGhostIterator(), [&](const OctreeGhostCellIterator3& iter) {
        marker->setOctagon(iter.levelCoor(), kGhost);
        marker->setTilt(iter.levelCoor(), kGhost);
    });
    iterate_half_tilt(lCoor, layout) {
        auto pos = marker->positionTilt(lCoor);
        if(isValidSdf(layout->levelGridSpacing(UNPACK_LEVEL3(lCoor)), Sdf->sample(pos))) {
            marker->setTilt(lCoor, kKnown);
        }
    }

    const char TILT_MARK = 0, OCTAGON_MARK = 1, HALF_TILT_MARK = 2;
    // <<char, LCOOR_T>, REAL>, 'char' used to indicate tilt(0), octagon(1), half-tilt(2)
    // NOTE: for extrapolate into boundary, the sdf sort should from positive to negative (large to small)
    //       whereas extraplating from levelset into air should from negative to positive (small to large)
    //       since there is offset when determine whether a point is valid, fabs() isn't suitable to combine these two case.
    REAL isExtrapolateIntoBoundary = (cacheMarker != nullptr) ? -1.0 : 1.0;
    auto compare = [&](const tuple<char, LCOOR_T, REAL>& a, const tuple<char, LCOOR_T, REAL>& b) {
        return std::get<2>(a) * isExtrapolateIntoBoundary > std::get<2>(b) * isExtrapolateIntoBoundary;
    };
    BENCHMARK(spdlog::info("marker init: {}", tt.durationInSeconds()); tt.reset();)

    int num_core = TaskManager::tf.num_workers();
    std::vector< std::vector<tuple<char, LCOOR_T, REAL>>> candidateArr(num_core);
    // initialize first ring of kTrial
    ThreadPool::parallelIterateGridWithCoreIdx(input->getOctagonGrid()->getParallelIteratorVec(), [&](const OctreeGridIterator3& iter, int coreId){
        if(marker->getOctagon(iter.levelCoor()) == kKnown) { // no need to check valid. If invalid, the tag should be kDontTouch
            for(int k = 0; k < 6; ++k) {    // octagon neighbor
                LCoor3 tp(iter.level, iter.xIdx()+OCT_NB_DXDYDZ[k].x(), iter.yIdx()+OCT_NB_DXDYDZ[k].y(), iter.zIdx()+OCT_NB_DXDYDZ[k].z());
                if(marker->getOctagonGrid()->isValid(tp.level, tp.xIdx, tp.yIdx, tp.zIdx) && marker->getOctagon(tp.levelCoor()) == kUnknown) {
                    marker->setOctagon(tp.levelCoor(), kTrial);
                    REAL l = Sdf->sample(marker->positionOctagon(tp.levelCoor()));
                    // heap.push(std::make_pair(std::make_pair(OCTAGON_MARK, tp.levelCoor()), l));
                    candidateArr[coreId].push_back(std::make_tuple(OCTAGON_MARK, tp.levelCoor(), l));
                }
            }
            if(iter.level >= 1) { // check half-tilt neighbor
                Size3 idx(iter.xIdx(), iter.yIdx(), iter.zIdx());
                for(int k = 0; k < 6; ++k) {
                    int axis = k >> 1; 
                    Size3 tIdx = idx * 2 + Size3::Ones();
                    tIdx(axis) += (k&1) ? 1 : -1;
                    LCOOR_T tLCoor = PACK_LCOOR3(iter.level-1, tIdx.x(), tIdx.y(), tIdx.z());
                    if(layout->get(tLCoor) >= 0 && tiltE->get(tLCoor).is_half_tilt) {
                        if(marker->getTilt(tLCoor) == kUnknown) {
                            marker->setTilt(tLCoor, kTrial);
                            REAL l = Sdf->sample(marker->positionTilt(iter.level-1, tIdx.x(), tIdx.y(), tIdx.z()));
                            // heap.push(std::make_pair(std::make_pair(HALF_TILT_MARK, tLCoor), l));
                            candidateArr[coreId].push_back(std::make_tuple(HALF_TILT_MARK, tLCoor, l));
                        }
                    }
                }
            }
            for(int k = 0; k < 8; ++k) {    // tilt neighbor
                LCoor3 tp(iter.level, iter.xIdx()+OCT_NB_TILT_DXDYDZ[k].x(), iter.yIdx()+OCT_NB_TILT_DXDYDZ[k].y(), iter.zIdx()+OCT_NB_TILT_DXDYDZ[k].z());
                if(marker->getTilt(tp.levelCoor()) == kUnknown) {
                    marker->setTilt(tp.levelCoor(), kTrial);
                    REAL l = Sdf->sample(marker->positionTilt(tp.levelCoor()));
                    // heap.push(std::make_pair(std::make_pair(tiltE->get(tp.levelCoor()).is_half_tilt ? HALF_TILT_MARK : TILT_MARK, tp.levelCoor()), l));
                    candidateArr[coreId].push_back(std::make_tuple(tiltE->get(tp.levelCoor()).is_half_tilt ? HALF_TILT_MARK : TILT_MARK, tp.levelCoor(), l));
                } else if(marker->getTilt(tp.levelCoor()) == kGhost) {
                    LCOOR_T c0 = layout->getUntil(tp.level, tp.xIdx, tp.yIdx, tp.zIdx, true);
                    c0 = c0 >= 0 ? c0 : layout->getUntil(tp.level, tp.xIdx, tp.yIdx, tp.zIdx, false);
                    ASSERT(c0 >= 0, "[extrapolation] initialize heap. ghost tilt can not find true tilt");
                    // TODO check this
                    if(UNPACK_LEVEL3(c0) < tp.level || (((UNPACK_X3(c0)<<1)==tp.xIdx) && ((UNPACK_Y3(c0)<<1)==tp.yIdx) && ((UNPACK_Z3(c0)<<1)==tp.zIdx))) {
                        if(marker->getTilt(c0) == kUnknown) {
                            marker->setTilt(c0, kTrial);
                            REAL l = Sdf->sample(marker->positionTilt(c0));
                            // heap.push(std::make_pair(std::make_pair(tiltE->get(c0).is_half_tilt ? HALF_TILT_MARK : TILT_MARK, c0), l));
                            candidateArr[coreId].push_back(std::make_tuple(tiltE->get(c0).is_half_tilt ? HALF_TILT_MARK : TILT_MARK, c0, l));
                        }
                    }
                }
            }
        }
        if(marker->getTilt(iter.levelCoor()) == kKnown && !tiltE->get(iter.levelCoor()).is_half_tilt) {
            auto neighborLCoor = tiltNeighborOctagon3(iter.level, iter.xIdx(), iter.yIdx(), iter.zIdx(), layout);
            for(int k = 0; k < 8; ++k) {
                LCOOR_T tLCoor = neighborLCoor(k);
                if(tLCoor == -1)
                    continue;
                if(marker->getOctagon(tLCoor) == kUnknown) {
                    marker->setOctagon(tLCoor, kTrial);
                    REAL l = Sdf->sample(marker->positionOctagon(tLCoor));
                    // heap.push(std::make_pair(std::make_pair(OCTAGON_MARK, tLCoor), l));
                    candidateArr[coreId].push_back(std::make_tuple(OCTAGON_MARK, tLCoor, l));
                }
            }
        }
    });
    int totalCandidate = 0;
    for(int i = 0; i < candidateArr.size(); ++i)
        totalCandidate += candidateArr[i].size();
    std::vector<tuple<char, LCOOR_T, REAL>> heapArr;
    heapArr.reserve(totalCandidate);
    for(int i = 0; i < candidateArr.size(); ++i) 
        heapArr.insert(heapArr.end(), candidateArr[i].begin(), candidateArr[i].end());

    iterate_half_tilt(lCoor, layout) {
        if(marker->getTilt(lCoor) == kKnown) {
            const auto& tNode = tiltE->get(lCoor);
            LCoor3 tp(UNPACK_LEVEL3(lCoor), UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor));
            // small octagon neighbor
            auto neighboorLCoor = halfTiltNeighborOctagon3(tp.level, tp.xIdx, tp.yIdx, tp.zIdx, tNode.half_direction);
            // big octagon neighbor
            Size3 idx(tp.xIdx>>1, tp.yIdx>>1, tp.zIdx>>1);
            idx((tNode.half_direction)%3) -= tNode.half_direction < 3 ? 1 : 0;
            // combine into single array
            neighboorLCoor[HALF_TILT_SLOT_INDEX[tNode.half_direction][4]] = PACK_LCOOR3(tp.level+1, idx.x(), idx.y(), idx.z());
            
            for(int k = 0; k < 8; ++k) {
                LCOOR_T tLCoor = neighboorLCoor(k);
                if(tLCoor != -1) {
                    if(marker->getOctagon(tLCoor) == kUnknown) {
                        marker->setOctagon(tLCoor, kTrial);
                        REAL l = Sdf->sample(marker->positionOctagon(tLCoor));
                        // heap.push(std::make_pair(std::make_pair(OCTAGON_MARK, tLCoor), l));
                        heapArr.emplace_back(std::make_tuple(OCTAGON_MARK, tLCoor, l));
                    }
                }
            }
        }
    }

    std::priority_queue<tuple<char, LCOOR_T, REAL>, std::vector<tuple<char, LCOOR_T, REAL>>, decltype(compare)> heap{compare, std::move(heapArr)};
    BENCHMARK(spdlog::info("heap init: {}", tt.durationInSeconds()); tt.reset();)

    // Propagate
    DEBUG_ONLY(REAL last = -1;);
    while (!heap.empty()) {
        auto item = heap.top();
        heap.pop();

        int cellMark = std::get<0>(item); //item.first.first;
        LCOOR_T cLCoor = std::get<1>(item);// item.first.second;
        int level = UNPACK_LEVEL3(cLCoor);
        if (fabs(std::get<2>(item)) > layout->levelGridSpacing(0) * maxDistanceFactor) {    // exceed maximum depth
            break;
        }
        // DEBUG_ONLY(
        //     ASSERT(last < item.second, "extrapolation: Heap error");
        //     last = item.second;
        // );
        if((cellMark == OCTAGON_MARK ? marker->getOctagon(cLCoor) : marker->getTilt(cLCoor)) == kKnown)   // skip 'kKnown' cell
            continue;
        int x = UNPACK_X3(cLCoor), y = UNPACK_Y3(cLCoor), z = UNPACK_Z3(cLCoor);
        Size3 cIdx(x, y, z);
        VT sum = zero<VT>();
        REAL count = 0.0;

        if(cellMark == TILT_MARK) {    // current cell is a tilt cell
            auto neighborLCoor = tiltNeighborOctagon3(level, x, y, z, layout);
            for(int k = 0; k < 8; ++k) {
                LCOOR_T tLCoor = neighborLCoor(k);
                if(tLCoor == -1)
                    continue;
                if(marker->getOctagon(tLCoor) == kKnown) {
                    sum += output->getOctagon(tLCoor);
                    count += 1;
                } else if(marker->getOctagon(tLCoor) == kUnknown) {
                    marker->setOctagon(tLCoor, kTrial);
                    heap.push(std::make_tuple(OCTAGON_MARK, tLCoor, Sdf->sample(marker->positionOctagon(tLCoor))));
                } else {
                    ASSERT(marker->getOctagon(tLCoor) != kGhost, "Tilt::neighbor get from util should not be ghost");
                }
            }
        } else if(cellMark == OCTAGON_MARK) {    // current cell is a octagon cell
            for(int k = 0; k < 6; ++k) {    // loop over its octagon neighbors
                Size3 tIdx = cIdx + OCT_NB_DXDYDZ[k];
                LCOOR_T tLCoor = PACK_LCOOR3(level, tIdx.x(), tIdx.y(), tIdx.z());
                if(marker->getOctagonGrid()->isValid(level, tIdx.x(), tIdx.y(), tIdx.z())) {
                    auto tMk = marker->getOctagon(tLCoor);
                    if(tMk == kKnown) {
                        sum += output->getOctagon(tLCoor);
                        count += 1;
                    } else if(tMk == kUnknown) {
                        marker->setOctagon(tLCoor, kTrial);
                        heap.push(std::make_tuple(OCTAGON_MARK, tLCoor, 
                            Sdf->sample(marker->positionOctagon(level, tIdx.x(), tIdx.y(), tIdx.z()))));
                    }
                }
            }
            if(level >= 1) { // check half-tilt neighbor
                for(int k = 0; k < 6; ++k) {
                    int axis = k / 2; 
                    Size3 tIdx = cIdx * 2 + Size3::Ones();
                    tIdx(axis) += (k&1) ? 1 : -1;
                    LCOOR_T tLCoor = PACK_LCOOR3(level-1, tIdx.x(), tIdx.y(), tIdx.z());
                    if(layout->get(tLCoor) >= 0 && tiltE->get(tLCoor).is_half_tilt) {
                        auto tiltMarker = marker->getTilt(tLCoor);
                        if(tiltMarker == kKnown) {
                            sum += output->getTilt(tLCoor);
                            count += 1;
                        } else if(tiltMarker == kUnknown) {
                            marker->setTilt(tLCoor, kTrial);
                            REAL l = Sdf->sample(marker->positionTilt(level-1, tIdx.x(), tIdx.y(), tIdx.z()));
                            heap.push(std::make_tuple(HALF_TILT_MARK, tLCoor, l));
                        }
                    }
                }
            }
            for(int k = 0; k < 8; ++k) {    // tilt neighbor
                LCoor3 tp(level, x+OCT_NB_TILT_DXDYDZ[k].x(), y+OCT_NB_TILT_DXDYDZ[k].y(), z+OCT_NB_TILT_DXDYDZ[k].z());
                auto tiltMarker = marker->getTilt(tp.levelCoor());
                if(tiltMarker == kKnown) {
                    sum += output->getTilt(tp.levelCoor());
                    count += 1;
                } else if(tiltMarker == kUnknown) {
                    marker->setTilt(tp.levelCoor(), kTrial);
                    REAL l = Sdf->sample(marker->positionTilt(tp.levelCoor()));
                    heap.push(std::make_tuple(tiltE->get(tp.levelCoor()).is_half_tilt ? HALF_TILT_MARK : TILT_MARK, tp.levelCoor(), l));
                    // if(tp.level == 1 && tp.xIdx == 22 && tp.yIdx == 15 && tp.zIdx == 20) {
                    //     spdlog::info("3-from:{}, {} {} {}", level, x, y, z);
                    // }
                } else if(tiltMarker == kGhost) {
                    LCOOR_T c0 = layout->getUntil(tp.level, tp.xIdx, tp.yIdx, tp.zIdx, true);
                    c0 = c0 >= 0 ? c0 : layout->getUntil(tp.level, tp.xIdx, tp.yIdx, tp.zIdx, false);
                    ASSERT(c0 >= 0, "[extrapolation] initialize heap. ghost tilt can not find true tilt");
                    if(UNPACK_LEVEL3(c0) < tp.level || ((UNPACK_X3(c0)<<1==tp.xIdx) && (UNPACK_Y3(c0)<<1==tp.yIdx) && (UNPACK_Z3(c0)<<1==tp.zIdx))) {
                        tiltMarker = marker->getTilt(c0);
                        if(tiltMarker == kKnown) {
                            sum += output->getTilt(tp.levelCoor());
                            count += 1;
                        } else if(tiltMarker == kUnknown) {
                            marker->setTilt(c0, kTrial);
                            REAL l = Sdf->sample(marker->positionTilt(c0));
                            heap.push(std::make_tuple(tiltE->get(c0).is_half_tilt ? HALF_TILT_MARK : TILT_MARK, c0, l));
                            // if(UNPACK_LEVEL3(c0) == 1 && UNPACK_X3(c0) == 22 && UNPACK_Y3(c0) == 15 && UNPACK_Z3(c0) == 20) {
                            //     spdlog::info("4-from:{}, {} {} {}", level, x, y, z);
                            // }
                        }
                    }
                }
            }
        } else if(cellMark == HALF_TILT_MARK) {
            const auto& tNode = tiltE->get(cLCoor);
            // small octagon neighbor
            auto neighboorLCoor = halfTiltNeighborOctagon3(level, x, y, z, tNode.half_direction);
            // big octagon neighbor
            Size3 idx(x>>1, y>>1, z>>1);
            idx((tNode.half_direction)%3) -= tNode.half_direction < 3 ? 1 : 0;
            // combine into single array
            neighboorLCoor[HALF_TILT_SLOT_INDEX[tNode.half_direction][4]] = PACK_LCOOR3(level+1, idx.x(), idx.y(), idx.z());

            for(int k = 0; k < 8; ++k) {
                LCOOR_T tLCoor = neighboorLCoor(k);
                if(tLCoor != -1) {
                    auto tMarker = marker->getOctagon(tLCoor);
                    if(tMarker == kKnown) {
                        sum += output->getOctagon(tLCoor);
                        count += 1;
                    } else {
                        ASSERT(tMarker != kGhost, "[halftilt] its neighbor should not be ghost");
                        marker->setOctagon(tLCoor, kTrial);
                        REAL l = Sdf->sample(marker->positionOctagon(tLCoor));
                        heap.push(std::make_tuple(OCTAGON_MARK, tLCoor, l));
                    }
                }
            }
        }

        DEBUG_ONLY(if(count == 0) {
            spdlog::info("mask:{}, level:{}, {} {} {}", int(cellMark), level, x, y, z);
        });
        ASSERT(count > 0, "Count is zero in extrapolation");

        if(cellMark == TILT_MARK || cellMark == HALF_TILT_MARK) {
            output->setTilt(cLCoor, sum/count);
            marker->setTilt(cLCoor, kKnown);
        } else {
            output->setOctagon(cLCoor, sum/count);
            marker->setOctagon(cLCoor, kKnown);
        }
    }

    BENCHMARK(spdlog::info("propagation: {}", tt.durationInSeconds()); tt.reset();)
}

//! Shared pointer type for the FmmLevelSetSolver3.
typedef std::shared_ptr<FmmLevelSetSolver3> FmmLevelSetSolver3Ptr;

}  // namespace jet

#endif  // INCLUDE_JET_FMM_LEVEL_SET_SOLVER3_H_

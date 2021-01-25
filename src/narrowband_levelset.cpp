#include "narrowband_levelset.h"
#include "jet/fmm_levelset_solver3.h"

using namespace Omni;
using namespace std;
using namespace jet;

inline REAL solveQuad(BitMapCellCenteredCharGrid3Ptr& markers, BitMapCellCenteredScalarGrid3Ptr& output,
                        const Size3& size, REAL gridSpacing, REAL sign, const Size3& idx) {
    Size3 hasXYZ = Size3::Zero();
    Vector3f phiXYZ = Vector3f::Zero();

    for(int i = 0; i < 6; ++i) {
        Size3 tIdx = idx + OCT_NB_DXDYDZ[i];
        if(Valid(tIdx, size) && markers->get(tIdx.x(), tIdx.y(), tIdx.z(), kUnknown) == kKnown) {
            ASSERT(fabs(output->get(tIdx)) < 10, "[solve Quad] value should exceed 10");
            int axis = i >> 1;
            REAL v1 = output->get(tIdx); 
            phiXYZ[axis] = hasXYZ[axis] == 0 ? v1 : (v1 < phiXYZ[axis] ? v1 : phiXYZ[axis]);
            hasXYZ[axis] = 1;
        }
    }

    DEBUG_ONLY(if(hasXYZ.sum() == 0) {
        spdlog::info("{} {} {}", idx.x(), idx.y(), idx.z());
        throw std::runtime_error("[ffm::solveQuad] hasX || hasY || hasZ == false");
    });

    REAL solution = 1000.0;
    for(int i = 0; i < 3; ++i) {
        if(hasXYZ[i])
            solution = std::min(solution, fabs(phiXYZ[i]));
    }
    solution += gridSpacing;
    solution *= sign;

    if(hasXYZ.sum() == 1)
        return solution;

    REAL a = hasXYZ.sum();
    REAL b = -phiXYZ.sum();
    REAL c = -squared(gridSpacing) + phiXYZ.squaredNorm();

    REAL det = b * b - a * c;
    solution = det < 0 ? solution : ((-b + sign * std::sqrt(det)) / a);
    
    return solution;
}

REAL NarrowbandLevelSet::sample(REAL pos_x, REAL pos_y, REAL pos_z) const {
    FLOAT rx = pos_x / gridSpacing();
    FLOAT ry = pos_y / gridSpacing();
    FLOAT rz = pos_z / gridSpacing();
    int xIdx, yIdx, zIdx;
    REAL fx, fy, fz;

    getBarycentric(rx-0.5, 0, resolution().x()-1, &xIdx, &fx);
    getBarycentric(ry-0.5, 0, resolution().y()-1, &yIdx, &fy);
    getBarycentric(rz-0.5, 0, resolution().z()-1, &zIdx, &fz);

    if(narrowbandMarker->get(xIdx, yIdx, zIdx) == false) {
        REAL sign = inOutMarker->get(xIdx, yIdx, zIdx) ? -1: 1;
        return sign * 1000;
    }

    return levelSet->sampleGivenRelativePos(xIdx, yIdx, zIdx, fx, fy, fz);
}


void NarrowbandLevelSet::reinitialize(REAL maxDistance) {
    Timer timer;
    auto newNarrowbandMarker = make_shared<BitMap>(res, delta_h);   // default init to 0
    auto marker = make_shared<BitMapCellCenteredCharGrid3>(resolution(), delta_h, narrowbandMarker, kUnknown);

    spdlog::info("\t\t\t init marker: {}, {}, {}", timer.durationInSeconds(), resolution().transpose(), delta_h); timer.reset();   // 10e-4
    auto compare = [](const pair<Size3, REAL>& a, const pair<Size3, REAL>& b) {
        return a.second > b.second;
    };
    
    static int last_count = 0;
    container.resize(0);
    container.reserve(last_count);
    priority_queue<pair<Size3, REAL>, vector<pair<Size3, REAL>>, decltype(compare)> trial(compare, std::move(container));

    std::set<LCOOR_T> intf_cells;
    std::mutex set_mutex;
    // for(int z = 0; z < resolution().z(); ++z) {
    ThreadPool::parallelForTF(0, resolution().z(), [&](int z){
        for(int y = 0; y < resolution().y(); ++y) {
        for(int x = 0; x < resolution().x(); ++x) {
            if(!narrowbandMarker->get(x, y, z)) 
                continue;
            Size3 idx(x, y, z);
            int cSign = getSign(x, y, z);
            for(int i = 0; i < 6; i+=2) {   // skip negative direction
                Size3 tIdx = idx + OCT_NB_DXDYDZ[i];
                if(Valid(tIdx, resolution()) && cSign != getSign(tIdx)) {
                    set_mutex.lock();
                    intf_cells.insert(PACK_COOR3(idx.x(), idx.y(), idx.z()));
                    intf_cells.insert(PACK_COOR3(tIdx.x(), tIdx.y(), tIdx.z()));
                    DEBUG_ONLY(if(fabs(get(tIdx)) > 10) { 
                        spdlog::info("{}: {}, {}: {}", idx.transpose(), get(idx), tIdx.transpose(), get(tIdx)); 
                        throw std::runtime_error("should be small");});
                    set_mutex.unlock();
                }
            }
        }}
    });

    spdlog::info("\t\t\t init narrow band: {} count: {}", timer.durationInSeconds(), intf_cells.size()); timer.reset();  // 0.1
    std::vector<std::pair<Size3, REAL>> intf_result_cache;
    intf_result_cache.resize(intf_cells.size());
    ThreadPool::parallelIterate(intf_cells.begin(), intf_cells.end(), intf_cells.size(), [&](const auto& iter, int arrIdx){
        LCOOR_T lCoor = *iter;
        Size3 idx(UNPACK_X3(lCoor), UNPACK_Y3(lCoor), UNPACK_Z3(lCoor));
        Vector3f correct_phi(MAX_FLOAT, MAX_FLOAT, MAX_FLOAT);
        Size3 correct_axis(0, 0, 0);
        int cSign = getSign(idx);
        for(int i = 0; i < 6; ++i) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[i];
            if(!Valid(tIdx, resolution()))
                continue;
            if(intf_cells.find(PACK_COOR3(tIdx.x(), tIdx.y(), tIdx.z())) == intf_cells.end()) {
                if(marker->get(tIdx, kUnknown) == kUnknown) {
                    set_mutex.lock();
                    trial.push(std::make_pair(tIdx, fabs(get(tIdx))));
                    set_mutex.unlock();
                    marker->set(tIdx, kTrial);
                }
            } else if(cSign != getSign(tIdx)) {
                REAL c_phi = theta(levelSet->get(idx), levelSet->get(tIdx)) * delta_h;
                int axis= i >> 1;
                correct_axis[axis]=1;
                correct_phi[axis]= std::min(correct_phi[axis], c_phi);
            }
        }
        // REAL oldL = levelSet->get(idx);
        ASSERT(correct_axis.sum() > 0, "correct sum zero");
        REAL hmnc_mean = 0;
        for(int i=0;i<3;i++){
            if(correct_axis[i]==0)
                continue;
            hmnc_mean += 1.0f/squared(correct_phi[i]);
        }
        hmnc_mean=sqrt(1.0/hmnc_mean);
        intf_result_cache[arrIdx] = make_pair(idx, getSign(idx) * hmnc_mean);
    });
    for(int i = 0; i < intf_result_cache.size(); ++i) {
        const Size3& idx = intf_result_cache[i].first;
        REAL tL = intf_result_cache[i].second;
        DEBUG_ONLY(if(fabs(tL) > EPSILON && fabs(levelSet->get(idx)) > EPSILON && getSign(idx) != sign(levelSet->get(idx))) {
            spdlog::info("oriL:{}, newL:{}, ori:{} new:{} idx:{}", levelSet->get(idx), tL, getSign(idx), sign(tL), idx.transpose());
            throw std::runtime_error("different sign in narrow band init");
        });
        levelSet->set(idx, tL);
        inOutMarker->set(idx, isInsideSdf(tL) ? true : false);
        marker->set(idx, kKnown);
        newNarrowbandMarker->set(idx, true);
        narrowbandMarker->set(idx, true);
    }
    intf_result_cache.resize(0); intf_result_cache.shrink_to_fit();
    spdlog::info("\t\t\t first ring of narrow band init: {}", timer.durationInSeconds()); timer.reset();    // 0.03
    int count = 0;
    // Propagate
    while (!trial.empty()) {
        // Timer tt;
        Size3 idx = trial.top().first;
        trial.pop();
        ++count;
        if(marker->get(idx) == kKnown)
            continue;
        ASSERT(marker->get(idx) == kTrial, "heap element should be trial");
        REAL tL = solveQuad(marker, levelSet, resolution(), delta_h, getSign(idx), idx);
        marker->set(idx, kKnown);

        DEBUG_ONLY(if(fabs(levelSet->get(idx)) > EPSILON && fabs(tL) > EPSILON && getSign(idx) != sign(tL)) {  // sign different
            spdlog::info("old:{} new:{}, {}", levelSet->get(idx), tL, idx.transpose());
            throw std::runtime_error("different sign");
        });
        levelSet->set(idx, tL);
        newNarrowbandMarker->set(idx, true);

        if (fabs(tL) > maxDistance) {
            break;
        }

        for(int m = 0; m < 6; ++m) {
            Size3 tIdx = idx + OCT_NB_DXDYDZ[m];
            if(Valid(tIdx, resolution()) && marker->get(tIdx.x(), tIdx.y(), tIdx.z(), kUnknown) != kKnown) {
                REAL p = solveQuad(marker, levelSet, resolution(), delta_h, getSign(tIdx), tIdx);
                if(marker->get(tIdx, kUnknown) != kKnown) {
                // || fabs(p) < fabs(get(tIdx))) {
                    marker->set(tIdx, kTrial);
                    narrowbandMarker->set(tIdx, true);
                    trial.push(std::make_pair(tIdx, fabs(p)));
                }
            }
        }
    }

    levelSet->setBitMap(newNarrowbandMarker);
    narrowbandMarker = newNarrowbandMarker;
    last_count = count;
    spdlog::info("\t\t\t propagation: {}, count: {}", timer.durationInSeconds(), count);// 0.4
}
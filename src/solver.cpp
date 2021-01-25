#include <spdlog/spdlog.h>
#include "solver.h"
#include "timer.h"
#include "MemUtilities.h"

using namespace Omni;

void Solver::update(REAL dt) {
    if(!is_init) {
        initialize();
    }
    Timer t;
#ifdef TILT_DOMAIN
    useFixedTimeStep = true;
#endif
    if(useFixedTimeStep) {
        spdlog::info("Use fixed timestep: {}", dt);
        TIMING(t, onUpdate(dt), "End onUpdate. Time cost: {} s");
        s_time += dt;
    } else {
        spdlog::info("Use adaptive sub-timestep. Current frame: {}", cFrame);

        REAL remainingTime = dt;
        while (remainingTime > EPSILON) {
            TIMING(t, int numSteps = numberOfSubTimeSteps(remainingTime), "Compute number of SubTimeSteps. Time costs: {} s");
            // NOTE: no need to do this now. CFL will consider effects of force in next step.
            // static bool firstIter = true;  // hack to avoid big time step in first update
            // if(firstIter) {
            //     numSteps = 4;
            //     firstIter = false;
            // }
            REAL actualTimeInterval = remainingTime / static_cast<REAL>(numSteps);

            spdlog::info("Number of remaining sub-timesteps: {}, dt:{}", 
                                numSteps, actualTimeInterval);
            spdlog::info("memory usage: peak : {} G, current: {} G", memutils::getPeakRSS() / 1024.0 / 1024.0 / 1024.0, memutils::getCurrentRSS() / 1024.0 / 1024.0 / 1024.0);

            TIMING(t, onUpdate(actualTimeInterval), "End onUpdate. Time cost: {} s");

            s_time += actualTimeInterval;
            remainingTime -= actualTimeInterval;
        }
    }
    endOfFrame();
    saveResult(++cFrame);
}
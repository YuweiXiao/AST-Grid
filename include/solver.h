#pragma once
#include "general.h"
#include "constants.h"
#include "timer.h"

namespace Omni {

class Solver {
public:
    Solver(const std::string& _name="solver")
        : name(_name) {
        Constants::init();
        spdlog::info(name);
    }

    virtual void initialize() = 0;

    void update(REAL dt);
    virtual int numberOfSubTimeSteps(REAL dt) const = 0;

    bool useFixedTimeStep = false;
    int getCurrentFrame() {return cFrame;}
    std::string getName() {return name;}

protected:
    virtual void initData() = 0;
    virtual void initObstacle() = 0;
    virtual void onUpdate(REAL dt) = 0;
    // save current frame's result
    virtual void saveResult(int index) = 0;
    virtual void endOfFrame() {};
    
    std::string name;
    REAL s_time = 0; // total simulation time
    int cFrame = 0;
    bool is_init=false;
    Timer timer;
};

} // end of namespace Omni
#pragma once
#include <mutex>

namespace Omni {

using namespace std;

class AtomicCounter {
public:
    void increase(int inc = 1) {
        m.lock();
        count += inc;
        m.unlock();
    }

    void decrease(int dec = 1) {
        m.lock();
        count -= dec;
        m.unlock();
    }

    volatile int count = 0;
private:
    mutex m;

};


} // Omni

#pragma once
#include "general.h"

namespace ast {
    
template<int d>
class BoxIterator {
    using T_Index = typename Eigen::Matrix<int, d, 1>;
public:
    BoxIterator(const T_Index& _idx_min, const T_Index& _idx_max)
        : idx_min(_idx_min), idx_max(_idx_max)
    {
        reset();
    }

    void next() {
        for(int i = 0; i < d; ++i) {
            if(idx(i) < idx_max(i) || i == d-1) {
                ++idx(i);
                return;
            } else idx(i) = idx_min(i);
        }
        // for(int i = d-1; i >= 0; --i) {
        //     if(idx(i) < idx_max(i) || i == 0) {
        //         ++idx(i);
        //         return;
        //     } else idx(i) = idx_min(i);
        // }
    }
    // bool valid() const { return idx(0) <= idx_max(0); }
    bool valid() const { return idx(d-1) <= idx_max(d-1); }
    const T_Index& index() const { return idx; }
    void reset() { idx = idx_min; }

private:
    T_Index idx_min, idx_max;
    T_Index idx;

};

} // namespace ast

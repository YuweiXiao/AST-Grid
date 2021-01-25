#pragma once
#include "general.h"

namespace Omni {

class Grid2 {
public:
    Grid2(const Size2& resolution, REAL dx, const Vector2f& domain_origin = Vector2f::Zero());
    virtual ~Grid2() {}

    const Size2& resolution() const {return res;}
    const REAL& gridSpacing() const {return spacing;}
    void checkCoord(bool valid, int xIdx, int yIdx, const std::string& message = "Grid2") const {
        if(!valid) {
            spdlog::error("[{}], Get index out of bound. x:{}, y:{}", message, xIdx, yIdx);
            throw std::runtime_error("Grid::Get index out of bound.");
        }
    }

private:
    Size2 res;
    REAL spacing;
    Vector2f domain_origin, domain_max;
};

class Grid3 {
public:
    Grid3(const Size3& _res, Real _spacing)
        : res(_res), spacing(_spacing)
    {}
    virtual ~Grid3() {}

    const Size3& resolution() const {return res;}
    Real gridSpacing() const {return spacing;}

private:
    Size3 res;
    Real spacing;
};

class GridIterator {
public:
    GridIterator(){};
    virtual void next() = 0;
    virtual bool isValid() = 0;
    virtual void reset() = 0;
    virtual int index() = 0;
    virtual ~GridIterator() {};
};

typedef std::shared_ptr<Grid3> Grid3Ptr;
typedef std::shared_ptr<Grid2> Grid2Ptr;

}
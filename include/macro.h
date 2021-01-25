#pragma once

#define USE_TASKFLOW
// #if __linux__
//     #define USE_OPENMP
// #endif

#ifdef __APPLE__
    // #define DEBUG 1
#endif
// #define DEBUG 1
#define USE_DOUBLE 1
#define OCTREE_LEVEL_RES_FACTOR 2
#define OCTREE_CORRECTION_BASED

#define SQRT2 1.4142135623730951
#define SQRT3 1.7320508075688772
#define PI 3.141592653589793
#define INV_SQRT2 1.0 / SQRT2
#define INV_SQRT3 1.0 / SQRT3
#define SOLVER_ACCURACY 1e-5

// #define NO_TILT
// #define TILT_DOMAIN 96

#define MASK(x, mask) ((x) & (mask))
#define SHIFT_MASK(x, offset, mask) (((x)>>(offset)) & (mask))
#define LEFT_SHIFT_MASK(x, offset, mask) (((x)<<(offset)) & (mask))
// #define PACK_LCOOR(level, x, y) (LEFT_SHIFT_MASK(level, 28, 0xf0000000) | LEFT_SHIFT_MASK(x, 14, 0x0fffc000) | MASK(y, 0x3fff))
#define PACK_LCOOR(level, x, y) (((level)<<28) | ((x) << 14) | (y))
#define PACK_COOR(x, y) (((x) << 14) | (y))
#define UNPACK_LEVEL(levelCoor) SHIFT_MASK(levelCoor, 28, 0xf)
#define UNPACK_COOR(levelCoor) MASK(levelCoor, 0x0fffffff)
#define UNPACK_X(levelCoor) SHIFT_MASK(levelCoor, 14, 0x3fff)
#define UNPACK_Y(levelCoor) SHIFT_MASK(levelCoor, 0, 0x3fff)

#define LCOOR_T int64_t

#define PACK_COOR3(x, y, z) (((((LCOOR_T)(x)) << 40) | (((LCOOR_T)(y))<<20)) | ((LCOOR_T)(z)))
#define PACK_LCOOR3(level, x, y, z) ((((LCOOR_T)(level))<<60) | PACK_COOR3(x, y, z))
#define UNPACK_LEVEL3(levelCoor) (int)SHIFT_MASK(levelCoor, 60, 0xf)
#define UNPACK_COOR3(levelCoor) MASK(levelCoor, (((LCOOR_T)(1)<<60)-1))
#define UNPACK_X3(levelCoor) (int)(SHIFT_MASK(levelCoor, 40, 0xfffff))
#define UNPACK_Y3(levelCoor) (int)(SHIFT_MASK(levelCoor, 20, 0xfffff))
#define UNPACK_Z3(levelCoor) (int)(SHIFT_MASK(levelCoor, 0, 0xfffff))

#define UNUSED_VARIABLE(a) ((void)a);
// #define USE_HIGH_BACKTRACE
// #define PACK_COOR3(x, y, z) ((((x) << 20) | ((y)<<10)) | (z))
// #define PACK_LCOOR3(level, x, y, z) (((level)<<30) | PACK_COOR3(x, y, z))
// #define UNPACK_LEVEL3(levelCoor) SHIFT_MASK(levelCoor, 30, 0x3)
// #define UNPACK_COOR3(levelCoor) MASK(levelCoor, 0x3fffffff)
// #define UNPACK_X3(levelCoor) SHIFT_MASK(levelCoor, 20, 0x3ff)
// #define UNPACK_Y3(levelCoor) SHIFT_MASK(levelCoor, 10, 0x3ff)
// #define UNPACK_Z3(levelCoor) SHIFT_MASK(levelCoor, 0, 0x3ff)
#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <assert.h>
#include <stdexcept>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/ostream.h>
#include "macro.h"
// using Eigen::Vector2f;
// using Eigen::Vector3f;
// using Eigen::Vector4f;
// using Eigen::VectorXf;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
typedef Eigen::Matrix<double, 8, 1> Vector8d;
typedef Eigen::Matrix<int, 8, 1> Vector8i;
typedef Eigen::Matrix<LCOOR_T, 8, 1> Vector8LCoor;
typedef Eigen::Matrix<double, 7, 1> Vector7d;
typedef Eigen::Matrix<double, 3, 7> Matrix37d;
typedef Eigen::Matrix<double, 8, 3> Matrix83d;
typedef Eigen::Matrix<double, 3, 4> Matrix34d;

using Eigen::Vector2i;
using Eigen::Vector3i;
using Eigen::Vector4i;

#define SPDLOG_TRACE_ON

#ifdef USE_DOUBLE
    using FLOAT = double;
    using REAL = double;
    using Real = double;
#else
    #define FLOAT float
    #define REAL float
    #define Real float
#endif

namespace Omni {

typedef unsigned long long ULL;
#define EPSILON std::numeric_limits<Real>::epsilon()
#define MAX_FLOAT std::numeric_limits<Real>::max()
typedef Eigen::Triplet<Real> FTriplet;
typedef Vector2i Size2;
typedef Vector3i Size3;
typedef Vector2d Vector2f;
typedef Vector3d Vector3f;
typedef Vector4d Vector4f;
typedef Vector8d Vector8f;
typedef VectorXd VectorXf;
typedef Matrix3d Matrix3f;
typedef Matrix4d Matrix4f;
typedef Vector7d Vector7f;
typedef Matrix37d Matrix37f;
typedef Matrix83d Matrix83f;

typedef std::function<void(int, int)> LoopFunc;
typedef std::function<void(int, int, int)> LoopFunc3;
typedef std::function<void(int, int, int, int)> LoopFunc4;

#ifdef DEBUG
#   define DEBUG_ONLY(a) a
#   define ASSERT(a, b) if(!(a)) throw std::runtime_error(b);
#else
#   define DEBUG_ONLY(a)
#   define ASSERT(a, b)
#endif

#ifdef WIN32
	#define u_int16_t uint16_t
#endif
}
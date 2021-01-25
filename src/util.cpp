#include "util.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h> 
#include <filesystem>
namespace fs = std::filesystem;

tf::Taskflow TaskManager::tf;
tf::Taskflow JobManager::tf;
tf::BasicTaskflow<tf::WorkStealingThreadpool> RealTaskManager::tf;
tf::BasicTaskflow<tf::WorkStealingThreadpool> ResultManager::tf(2);

namespace Omni {

std::string resToStr(int x, int y, int z) {
    std::string ret = std::to_string(x) + "_" + std::to_string(y);
    if(z > 0) {
        ret += "_" + std::to_string(z);
    }
    return ret;
}

// e.g input: "3", '0', 4, output : "0003"
std::string paddingStr(const std::string &str, char c, int target_length) {
    int num = target_length - str.size();
    if(num < 0)
        throw std::runtime_error("[padding Str] num is negative :" + std::to_string(num));
    std::string front = std::string(num, c);
    return front + str;
}

Vector8f buildVector8f(const Real& a, const Real& b, const Real& c, const Real& d,
        const Real& e, const Real& f, const Real& g, const Real& h) {
    Vector8f ret;
    ret << a, b, c, d, e, f, g, h;
    return ret;
}

bool isDirExist(const char* filename) {
    return fs::exists(std::string(filename));
}

bool isDirExist(const std::string filename) {
    return fs::exists(filename);
}

// filename : xxx/xxx/xxx.pdf
bool createDir(std::string filename) {
    auto p = filename.find_last_of("/");
    if(p != std::string::npos) {
        filename = filename.substr(0, p);
    } else {
        return true;
    }
    return fs::create_directories(filename);
}

// xxx/xxx/xxx.pdf
bool createDir(const char* filename) {
    return createDir(std::string(filename));
}

// p0 is out of boundary
// p1 is inside of grid
// Vector2f clampToBoundary(Grid2Ptr grid, const Vector2f& p0, const Vector2f& p1) {
//     // spdlog::info("{},{},{},{}", p0.transpose(), p1.transpose(), grid->gridDomainX().transpose(), grid->gridDomainY().transpose());
//     // assert(!grid->isInDomain(p0) && grid->isInDomain(p1));

    // Vector2f ret = p0, delta_p = p0 - p1;
    // // clamp x first
    // ret.x() = fmin(fmax(0.0, ret.x()), grid->gridDomain().x() - EPSILON);
    // ret.y() = p1.y() + delta_p.y() * (ret.x() - p1.x()) / (delta_p.x() + EPSILON);

    // // clamp y
    // ret.y() = fmin(fmax(0.0, ret.y()), grid->gridDomain().y() - EPSILON);
    // ret.x() = p1.x() + delta_p.x() * (ret.y() - p1.y()) / (delta_p.y() + EPSILON);

    // return ret;
    // Vector2f left = p0;
    // Vector2f right = p1;
    // int precision = 6;
    // for(int i = 0; i < precision; ++i) {
    //     Vector2f mid = (left + right) * 0.5;
    //     if(grid->isInDomain(mid)) {
    //         right = mid;
    //     } else {
    //         left = mid;
    //     }
    // }
    // spdlog::info("p0:{}, p1:{}, distance:{}", p0.transpose(), p1.transpose(), (p1-p0).norm());
    // return right;
// }

// Size2 clampIndexToInner(const Size2& idx, const Size2& xBound, const Size2& yBound) {
//     int x = std::max(xBound(0), idx.x());
//     x = std::min(xBound(1), x);
//     int y = std::max(yBound(0), idx.y());
//     y = std::min(yBound(1), y);
//     return Size2(x, y);
// }

// Real clampPosToInner(Real x, const Vector2f& bound) {
//     return fmax(bound(0), fmin(x, bound(1)));
// }

// Vector2f clampPosToInner(const Vector2f& idx, const Vector2f& xBound, const Vector2f& yBound) {
//     Real x = clampPosToInner(idx.x(), xBound);
//     Real y = clampPosToInner(idx.y(), yBound);
//     return Vector2f(x, y);
// }

Vector3f grid3ToOpenGL(Vector3f p)
{
	return Vector3f(p.x(), p.z(), -p.y());
}

bool isSymmetrical(Eigen::SparseMatrix<Real>& A)
{
	bool ret = true;
	// int n = A.rows();
	Eigen::SparseMatrix<Real> res = Eigen::SparseMatrix<Real>(A.transpose()) - A;
	for (unsigned int k = 0; k < res.outerSize(); k++) {
		for (Eigen::SparseMatrix<Real>::InnerIterator it(res, k); it; ++it)
		{
			if (std::fabs(it.value()) > EPSILON) {
				spdlog::info("row:{} col:{} value: {}", it.row(), it.col(), A.coeffRef(it.row(), it.col()));
				return false;
			}
		}
	}
	return ret;
}

// TODO try to use more colorful colormap
// ! https://github.com/libigl/libigl/blob/master/include/igl/colormap.cpp
// from https://stackoverflow.com/a/7811134
Vector3f getHeatColor(Real v,Real vmin,Real vmax) {
    // x->r, y->g, z->b
    Vector3f c(1.0,1.0,1.0); // white
    Real dv;

    v = v < vmin ? vmin : v;
    v = v > vmax ? vmax : v;
    dv = vmax - vmin;

    if (v < (vmin + 0.25 * dv)) {
        c.x() = 0;
        c.y() = 4 * (v - vmin) / dv;
    } else if (v < (vmin + 0.5 * dv)) {
        c.x() = 0;
        c.z() = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
    } else if (v < (vmin + 0.75 * dv)) {
        c.x() = 4 * (v - vmin - 0.5 * dv) / dv;
        c.z() = 0;
    } else {
        c.y() = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
        c.z() = 0;
    }

    return c;
}

REAL interpolateMLS(Vector3f pos, Vector3f p[4], REAL value[4]) {
    Vector4f v(value[0], value[1], value[2], value[3]);
    Vector4f w;
    REAL delta = (0.03 * (p[1]-p[2]).norm());
    delta *= delta;
    Matrix4d B; // B^T
    for(int i = 0; i < 4; ++i) {
        REAL l = (p[i] - pos).norm();
        if(l < EPSILON) return value[i];
        w[i] = 1.0 / (l * l + delta);
        w[i] *= w[i];
    }
    B << Vector4f(1.0, p[0].x(), p[0].y(), p[0].z()),
        Vector4f(1.0, p[1].x(), p[1].y(), p[1].z()),
        Vector4f(1.0, p[2].x(), p[2].y(), p[2].z()),
        Vector4f(1.0, p[3].x(), p[3].y(), p[3].z());

    Matrix4d W2(w.asDiagonal());   // W^2
    
    REAL ret = Vector4f(1.0, pos.x(), pos.y(), pos.z()).transpose() * 
        (B * W2 * B.transpose()).inverse() * B * W2 * v;
    // if(std::isnan(ret) || debug) {
    //     std::cout<<"------------------------"<<std::endl;
    //     std::cout<<B<<std::endl;
    //     std::cout<<W2<<std::endl;
    //     std::cout<<v.transpose()<<std::endl;
    //     std::cout<<(B * W2 * B.transpose()).inverse()<<std::endl;
    // }
    return ret;
}

Vector3f interpolateMLS(Vector3f pos, Vector3f p[4], Vector3f value[4]) {
    Eigen::Matrix<REAL, 4, 3> v;
    v << value[0].transpose(), value[1].transpose(), value[2].transpose(), value[3].transpose();
    Vector4f w;
    REAL delta = (0.03 * (p[1]-p[2]).norm());
    delta *= delta;
    Matrix4d B; // B^T
    for(int i = 0; i < 4; ++i) {
        REAL l = (p[i] - pos).norm();
        if(l < EPSILON)
            return value[i];
        w[i] = 1.0 / (l * l + delta);
        w[i] *= w[i];
    }
    B << Vector4f(1.0, p[0].x(), p[0].y(), p[0].z()),
        Vector4f(1.0, p[1].x(), p[1].y(), p[1].z()),
        Vector4f(1.0, p[2].x(), p[2].y(), p[2].z()),
        Vector4f(1.0, p[3].x(), p[3].y(), p[3].z());
    Matrix4d W2(w.asDiagonal());   // W^2
    
    // std::cout<<B<<std::endl;
    // std::cout<<W2<<std::endl;
    // std::cout<<(B * W2 * B.transpose()).inverse()<<std::endl;

    return Vector4f(1.0, pos.x(), pos.y(), pos.z()).transpose() * 
        (B * W2 * B.transpose()).inverse() * B * W2 * v;
}

} // end of namespace Omni
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

} // end of namespace Omni
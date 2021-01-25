#pragma once
#include <cmath>
#include <string>
#include <fstream>
#include <execinfo.h>
#include "general.h"
#include "grid.h"
#include <taskflow/taskflow.hpp>


class TaskManager {
public:
	static tf::Taskflow tf;

	template<typename Callable>
	static void dispatch(Callable func) {
		tf.emplace(func);
		tf.dispatch();
	}
};

class JobManager {
public:
	static tf::Taskflow tf;
	
	template<typename Callable>
	static void dispatch(Callable func) {
		tf.emplace(func);
		tf.dispatch();
	}
};

class RealTaskManager {
public:
	static tf::BasicTaskflow<tf::WorkStealingThreadpool> tf;
	
	template<typename Callable>
	static void dispatch(Callable func) {
		tf.emplace(func);
		tf.dispatch();
	}
};

class ResultManager {
public:
	static tf::BasicTaskflow<tf::WorkStealingThreadpool> tf;

	template<typename Callable>
	static void dispatch(Callable func) {
		tf.emplace(func);
		tf.dispatch();
	}
};

namespace Omni {

std::string resToStr(int x, int y, int z = -1);
// e.g input: "3", '0', 4, output : "0003"
std::string paddingStr(const std::string &str, char c, int target_length);

// Vector8f buildVector8f(const Real& a=0, const Real& b=0, const Real& c=0, const Real& d=0,
						// const Real& e=0, const Real& f=0, const Real& g=0, const Real& h=0);

#pragma region SYSTEM_UTIL

// filename: e.g. ./build/hello/xxxx.png
bool isDirExist(const char* filename);

// filename: e.g. ./build/hello/xxxx.png
bool createDir(const char* filename);

bool isDirExist(const std::string filename);

bool createDir(const std::string filename);

#pragma endregion SYSTEM_UTIL


// Siqi Wang
// Input: coordinate of p, (x,y) are both [0,1], [0.5,0.6]->0, but cannot be [0.5,0.5]!
// Output: the code of region (clockwise 0-7)
inline int getRegion2D(const Vector2f& p) {
	int region = (p(0) > 0.5 || (p(0) == 0.5 && p(1) > 0.5) ? 0 : 4) 
                + ((p(1)-0.5) * (p(0)-0.5) > 0 || p(0) == 0.5 ? 0 : 2) 
                + ((abs(p(0)-0.5) - abs(p(1)-0.5)) * (p(1)-0.5) * (p(0)-0.5) > 0 ? 1 : 0);

    DEBUG_ONLY(if(region < 0 || region > 7) {
        throw std::runtime_error("[getRegion2D] out of region:" 
                + std::to_string(region) + ", :" 
                + std::to_string(p.x()) + ", :" 
                + std::to_string(p.y()));
    });

    return region;
}


// Input: normalized x, y, z. \in [0, 1]
// Output: region in 3D cube. 
//		positive x, y, z : 0, 1, 2
//		negative x, y, z : 3, 4, 5
inline int getRegion3D(Real nx, Real ny, Real nz) {
	DEBUG_ONLY(if(nx < 0 || nx > 1 || ny < 0 || ny > 1 || nz < 0 || nz > 1) {
        throw std::runtime_error("[getRegion3D] invalid normalized input: x " + 
            std::to_string(nx) + ", y " + std::to_string(ny) + ", z " + std::to_string(nz));
    });

    nx -= 0.5;ny -= 0.5;nz -= 0.5;
    Real abs_nx = fabs(nx);
    Real abs_ny = fabs(ny);
    Real abs_nz = fabs(nz);
    int r = 0;
    if(abs_ny >= abs_nx && abs_ny >= abs_nz) {
        r = ny >= 0 ? 1 : 4;
    } else if(abs_nz >= ny && abs_nz >= abs_nx) {
        r = nz >= 0 ? 2 : 5;
    } else {
        r = nx >= 0 ? 0 : 3;
    }
    
    return r;
}

// Chen Siyu
// Input : 3 point p0, p1, p2 and their value v0, v1,v2
// Output: value in position p 
template<typename T>
T interpolate2D(const Vector2f& p, const Vector2f& p0, const Vector2f& p1, const Vector2f& p2, 
			const T& v0, const T& v1, const T& v2)
{
	// denominator
	Real denom = (p1(1) - p2(1)) * (p0(0) - p2(0)) +
		(p2(0) - p1(0)) * (p0(1) - p2(1));
	
	Real w0 = (p1(1) - p2(1)) * (p(0) - p2(0))
		+ (p2(0) - p1(0)) * (p(1) - p2(1));
	w0 = w0 / denom;

	Real w1 = (p2(1) - p0(1)) * (p(0) - p2(0)) +
		(p0(0) - p2(0)) * (p(1) - p2(1));
	w1 = w1 / denom;

	Real w2 = 1 - w0 - w1;

    return w0 * v0 + w1 * v1 + w2 * v2;
}

template<typename T>
T interpolateMAC(const Vector2f& p, const Vector2f& p0, const Vector2f& p1, const Vector2f& p2,
	const Vector2f& p3, const T& v0, const T& v1, const T& v2, const T& v3)
{
	// denominator
	Real denom = fabs((p2(1) - p0(1)) * (p2(0) - p0(0))) + EPSILON;
	
	// weights
	Real w0 = fabs((p2(0) - p(0)) * (p2(1) - p(1))) / denom;  
	Real w1 = fabs((p3(0) - p(0)) * (p3(1) - p(1))) / denom;
	Real w2 = fabs((p0(0) - p(0)) * (p0(1) - p(1))) / denom;
	Real w3 = fabs((p1(0) - p(0)) * (p1(1) - p(1))) / denom;

    return w0 * v0 + w1 * v1 + w2 * v2 + w3 * v3;
}

template<typename S, typename T>
inline S lerp(const S& v0, const S& v1, T f)
{
	// f = (x - x1) / (x2 - x1)
    return v0 * (1 - f) + v1 * f;
}

template<typename S, typename T>
inline S bilerp(const S& v00, const S& v01, const S& v11, const S& v10, T fx, T fy)
{
	// fx = (x - x00) / (x11 - x00)
    return lerp(lerp(v00, v10, fx), lerp(v01, v11, fx), fy);
}

template<typename S, typename T>
inline S trilerp(
    const S& f000, const S& f100, const S& f010, const S& f110, 
	const S& f001, const S& f101, const S& f011, const S& f111,
    T tx, T ty, T fz) {
    return lerp(
        bilerp(f000, f010, f110, f100, tx, ty),
        bilerp(f001, f011, f111, f101, tx, ty),
        fz);
}

template<typename S>
inline S trilerp(
    const S& f000, const S& f100, const S& f010, const S& f110, 
	const S& f001, const S& f101, const S& f011, const S& f111,
    const Vector3f& f) {
    return lerp(
        bilerp(f000, f010, f110, f100, f.x(), f.y()),
        bilerp(f001, f011, f111, f101, f.x(), f.y()),
        f.z());
}

template<typename T>
T interpolateByE(const Real& e0, const Real& e1, const Vector2f& spacing,
	const T& v0, const T& v1) {
	Real t0 = spacing.x() - e0;
	Real t1 = spacing.x() - e1;
	return (v0 * t1 + v1 * t0) / (t0 + t1 + EPSILON);
}

// inverse distance interpolation in pentagon
template<typename T>
T interpolatePentagon(const Vector2f& p, const Vector2f& p0, const Vector2f& p1, 
		const Vector2f& p2, const Vector2f& p3, const Vector2f& p4, 
		const T& v0, const T& v1, const T& v2, const T& v3, const T& v4) {
	Real dis0 = (p0 - p).norm();
	Real dis1 = (p1 - p).norm();
	Real dis2 = (p2 - p).norm();
	Real dis3 = (p3 - p).norm();
	Real dis4 = (p4 - p).norm();
	Real totalDis = dis0 + dis1 + dis2 + dis3 + dis4;

	dis0 /= dis0 / totalDis  + EPSILON;
	dis1 /= dis1 / totalDis + EPSILON;
	dis2 /= dis2 / totalDis + EPSILON;
	dis3 /= dis3 / totalDis + EPSILON;
	dis4 /= dis4 / totalDis + EPSILON;

	return v0 * (1.0/dis0) + v1 * (1.0/dis1) + v2 * (1.0/dis2) + 
			v3 * (1.0/dis3) + v4 * (1.0/dis4);
}

template<typename T>
T interpolateTilt(const Vector2f& p, const Vector2f& p0, const Vector2f& p1, const Vector2f& p2,
	const Vector2f& p3, const T& v0, const T& v1, const T& v2, const T& v3)
{
	Vector2f length = Vector2f((p1-p0).norm(), (p3-p0).norm());
	Vector2f dir = (p - p0);
	Vector2f unitY = (p3-p0).normalized();
	Vector2f unitX = (p1-p0).normalized();
	Vector2f normalized = Vector2f(dir.dot(unitX) / length.x(),
									dir.dot(unitY) / length.y()) ;
	return bilerp(v0, v1, v2, v3, normalized.y(), normalized.x());
}

template<typename Derived1, typename Derived2, typename Index1, typename Index2>
inline void getBaryCentric(const Eigen::MatrixBase<Derived1>& x,  
	const Eigen::MatrixBase<Index1>& high,
	Eigen::MatrixBase<Index2>& i, Eigen::MatrixBase<Derived2>& f)
{
	static_assert(std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value, "should have same scalar type");
	static_assert(std::is_same<typename Index1::Scalar, typename Index2::Scalar>::value, "should have same scalar type");
	i = (x.template cast<typename Index1::Scalar>())
		 .cwiseMax(Eigen::MatrixBase<Index2>::Zero())
		 .cwiseMin(high-Eigen::MatrixBase<Index2>::Ones());
	f = (x - i.template cast<typename Derived1::Scalar>())
		 .cwiseMax(Eigen::MatrixBase<Derived2>::Zero())
		 .cwiseMin(Eigen::MatrixBase<Derived2>::Ones());
}

template<typename T>
inline void getBarycentric(T x, int iLow, int iHigh,
    int* i, T* f) {
    T s = std::floor(x);
    *i = static_cast<int>(s);

    int offset = -iLow;
    iLow += offset;
    iHigh += offset;

    if (iLow == iHigh) {
        *i = iLow;
        *f = 0;
    } else if (*i < iLow) {
        *i = iLow;
        *f = 0;
    } else if (*i > iHigh - 1) {
        *i = iHigh - 1;
        *f = 1;
    } else {
        *f = static_cast<T>(x - s);
    }

    *i -= offset;
}

// http://mathworld.wolfram.com/Tetrahedron.html
// https://math.stackexchange.com/questions/1603651/volume-of-tetrahedron-using-cross-and-dot-product
template<typename T>
inline T interpolateTetrahedron(const Vector3f& p, 
		const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const Vector3f& p3,
		const T& v0, const T& v1, const T& v2, const T& v3/* , bool& valid */) {
	// Real volume = fabs((p1 - p0).dot((p2 - p0).cross(p3 - p0)));

	// Real w0 = fabs((p1 - p).dot((p2-p).cross(p3-p)));
	// Real w1 = fabs((p0 - p).dot((p2-p).cross(p3-p)));
	// Real w2 = fabs((p1 - p).dot((p0-p).cross(p3-p)));
	// Real w3 = fabs((p1 - p).dot((p2-p).cross(p0-p)));

	// DEBUG_ONLY(if(volume < EPSILON) {
	// 	spdlog::error("[interpolateTetrahedron] zero volume:{}. p:{}, p0:{}, p1:{}, p2:{}, p3:{}",
	// 		volume, p.transpose(), p0.transpose(), p1.transpose(), p2.transpose(), p3.transpose());
	// 	throw std::runtime_error("[interpolateTetrahedron] zero volume.");
	// });

	// T ret = (w0 * v0 + w1 * v1 + w2 * v2 + w3 * v3) / volume;
	Eigen::Matrix3d m;
	m << (p0-p3), (p1-p3), (p2-p3);
	Vector3f l = m.inverse() * (p - p3);
	return v0 * l[0] + v1 * l[1] + v2 * l[2] + v3 * (1 - l.sum());

	// return ret;
}

REAL interpolateMLS(Vector3f pos, Vector3f p[4], REAL v[4]);
Vector3f interpolateMLS(Vector3f pos, Vector3f p[4], Vector3f v[4]);


Vector3f grid3ToOpenGL(Vector3f p);

Vector3f getHeatColor(Real v,Real vmin,Real vmax);


template <typename T> int sign(T val) {return (val > T(0)) ? 1: -1;}
template <typename T> inline T squared(T e) {return e*e;}
template <typename T> inline T cubed(T e) {return e*e*e;}
template <typename T> inline bool isInsideSdf(T phi) {return phi <= 0;}


template<class ArrayT> bool All_Less(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]>=a1[i])return false;}return true;}
template<class ArrayT> bool All_Less_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]>a1[i])return false;}return true;}
template<class ArrayT> bool All_Greater(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]<=a1[i])return false;}return true;}
template<class ArrayT> bool All_Greater_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]<a1[i])return false;}return true;}
template<class ArrayT> bool Has_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++)if(a0[i]==a1[i])return true;return false;}
template<class ArrayT> bool Has_Less_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++)if(a0[i]<=a1[i])return true;return false;}
template<class ArrayT> bool Has_Greater_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++)if(a0[i]>=a1[i])return true;return false;}

// coord >= 0 && coord < res
inline bool Valid(const Size2& coord,const Size2& res) {return All_Greater_Equal<Size2>(coord, Size2::Zero()) && All_Less(coord, res);}
inline bool Valid(const Size3& coord,const Size3& res) {return All_Greater_Equal<Size3>(coord, Size3::Zero()) && All_Less(coord, res);}

//Given two signed distance values, determine what fraction of a connecting segment is "inside"
inline REAL fractionInside(REAL phi_left, REAL phi_right) {
	if(phi_left < 0 && phi_right < 0)
		return 1;
	if (phi_left < 0 && phi_right >= 0)
		return phi_left / (phi_left - phi_right);
	if(phi_left >= 0 && phi_right < 0)
		return phi_right / (phi_right - phi_left);
	else
		return 0;
}

inline int reverseTiltIdx(int x) {
	return (x&1?8:6) - x;
}

inline REAL theta(const REAL phi_1,const REAL phi_2) {
    return phi_1/(phi_1-phi_2);
}

// NOTE: calculate number of 1 bits in <<n>>
inline int hammingWeight(int n) {
    int sum = 0;
    while (n != 0) {
        sum++;
        n &= (n - 1);
    }
    return sum;
}

// TODO fix, don't compile in windows
inline void printCallStack() {
	// https://linux.die.net/man/3/backtrace_symbols
	void* buffer[200];
	int nptrs = backtrace(buffer, 200);
	printf("backtrace() returned %d addresses\n", nptrs);
	char **strings = backtrace_symbols(buffer, nptrs);
	for (int j = 0; j < nptrs; j++)
			printf("%s\n", strings[j]);
}

template<class T> inline T zero() { return (T)0; }
template<> inline Vector2f zero<Vector2f>() { return Vector2f::Zero(); }
template<> inline Vector3f zero<Vector3f>() { return Vector3f::Zero(); }
template<> inline Vector7f zero<Vector7f>() { return Vector7f::Zero(); }
template<> inline Vector8f zero<Vector8f>() { return Vector8f::Zero(); }
template<> inline Eigen::Matrix2d zero<Eigen::Matrix2d>() { return Eigen::Matrix2d::Zero(); }

}


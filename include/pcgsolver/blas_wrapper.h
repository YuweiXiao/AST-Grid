#pragma once

// Simplified BLAS wrapper (overloaded readable names)
// with stride==1 assumed, and for std::vector's.
// For the moment, no complex number support, and many routines have been dropped.

#include <vector>
#include <numeric>
#include "general.h"
#include "thread_pool.h"
#if defined(USE_OPENMP)
#include <omp.h>
#endif

namespace BLAS{

inline double dot(const std::vector<double> &x, const std::vector<double> &y) {
    double sum = 0;
#ifdef USE_TASKFLOW
	int coreNum = TaskManager::tf.num_workers();
	std::vector<double> sums(coreNum + 1);

	ThreadPool::parallelForTFWithCoreID(0, (int)x.size(), [&](int i, int coreID){
		sums[coreID] += x[i] * y[i];
	});

    sum = std::accumulate(sums.begin(), sums.end(), 0.0);
#elif defined(USE_OPENMP)
	#pragma omp parallel for reduction(+: sum)
	for (int i = 0; i < x.size(); i++) {
		sum += x[i] * y[i];
	}
#else
	for (int i = 0; i < x.size(); i++) {
		sum += x[i] * y[i];
	}
#endif
	return sum; 
}


// inf-norm (maximum absolute value)
inline double abs_max(const std::vector<double> &x) {
    double max_val = 0;
#ifdef USE_TASKFLOW
	// tf::Taskflow tf;
	TaskManager::tf.transform_reduce(x.begin(), x.end(), max_val, 
			[](double a, double b) { return std::max(a, b);	},
			[](double c) { return std::fabs(c); } 
			);

	TaskManager::tf.wait_for_all();
#elif defined(USE_OPENMP)
	#pragma omp parallel for reduction(max:max_val) 
	for (int i = 0; i < x.size(); i++) {
		if (std::fabs(x[i]) > max_val)
			max_val = std::fabs(x[i]);
	}	
#else
	for (int i = 0; i < x.size(); i++) {
		if (std::fabs(x[i]) > max_val)
			max_val = std::fabs(x[i]);
	}
#endif
	return max_val;
}


// (y=alpha*x+y)
inline void add_scaled(double alpha, const std::vector<double> &x, std::vector<double> &y) {
#ifdef USE_TASKFLOW
	ThreadPool::parallelForTF((size_t)0, x.size(), [&](size_t i){
		y[i] += alpha * x[i];
	});
#elif defined(USE_OPENMP)
	#pragma omp parallel for 
	for (int i = 0; i < x.size(); i++) {
		y[i] = alpha * x[i] + y[i];
	}
#else
	for (int i = 0; i < x.size(); i++) {
		y[i] = alpha * x[i] + y[i];
	}
#endif
}


// (x=alpha*x)
inline void scale(double alpha, std::vector<double> &x) {
#ifdef USE_TASKFLOW
    ThreadPool::parallelForTF((size_t)0, x.size(), [&](size_t i){
		x[i] *= alpha;
	});
#else
	for (int i = 0; i < x.size(); i++) {
		y[i] = alpha * x[i] + y[i];
	}
#endif
}

inline double sum(const std::vector<double>& data) {
	std::vector<double> sums(TaskManager::tf.num_workers(), 0);
	ThreadPool::parallelForPlainTF(0, (int)data.size(), [&](int start, int end, int coreID){
		double tmp = 0;
		for(int i = start; i < end; ++i)
			tmp += data[i];
		sums[coreID] = tmp;
	});
	return std::accumulate(sums.begin(), sums.end(), 0.0);
}

inline double dot(const Eigen::VectorXd &x, const Eigen::VectorXd &y) { return x.dot(y); }
inline double abs_max(const Eigen::VectorXd &x) { return x.cwiseAbs().maxCoeff(); }
inline void add_scaled(double alpha, const Eigen::VectorXd &x, Eigen::VectorXd &y) { y.noalias() += alpha * x; }
inline void scale(double alpha, Eigen::VectorXd &x) { x *= alpha; }

} // end of namespace BLAS

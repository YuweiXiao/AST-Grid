// ref https://ideone.com/Z7zldb
#pragma once
////////////////////////////////////////////////////////////////////////////////
#include <thread>
#include <vector>
#include <cmath>
#include "util.h"
#include "box_iterator.h"
////////////////////////////////////////////////////////////////////////////////

class ThreadPool {
public:
	
	template <typename I>
	static std::vector<I> initParallelIterator(I iter, int slice = 0) {
		int coreNum = TaskManager::tf.num_workers();
		int nodeNum = iter.totalSize();
		auto tmp = iter;
		slice = (slice == 0) ? ((nodeNum + coreNum - 1) / coreNum) : slice;

		std::vector<I> res;
		while(iter.isValid()) {
			for (int j = 0; j < slice && tmp.isValid(); j++, tmp.next()) {
				//spdlog::info("tmp idx{} level{} x{} y{}", tmp.index(), tmp.level, tmp.xIdx(), tmp.yIdx());
				;
			}
			auto q = iter;
			q.setUpperBound(slice);
			res.push_back(q);
			iter = tmp;
		}
		return res;
	}


	template <typename I, typename F>
	static void parallelIterateGrid(std::vector<I> iterators, const F& f) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		DEBUG_ONLY(if(TaskManager::tf.num_workers() < iterators.size()) { spdlog::info("Num of iterators is larger than num of process, Please consider use parallelIterateGridRealTask instead!");});
		for(I iter : iterators) {
			auto task = TaskManager::tf.silent_emplace([iter, f]() mutable {
				for (; iter.isValid(); iter.next()) {
					//spdlog::info("iter idx{} level{} x{} y{}", iter.index(), iter.level, iter.xIdx(), iter.yIdx());
					f(iter);
				}
			});
			source.precede(task);
			task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}

	// parallel API For iterators num larger than process number
	template <typename I, typename F>
	static void parallelIterateGridRealTask(std::vector<I> iterators, const F &f) {
		auto source = RealTaskManager::tf.placeholder();
		auto target = RealTaskManager::tf.placeholder();
		DEBUG_ONLY(if(RealTaskManager::tf.num_workers() >= iterators.size()) { spdlog::info("Num of iterators is smaller than num of process, Please consider use parallelIterateGrid instead!");});
		for(I iter : iterators) {
			auto task = RealTaskManager::tf.silent_emplace([iter, f]() mutable {
				for (; iter.isValid(); iter.next()) {
					//spdlog::info("iter idx{} level{} x{} y{}", iter.index(), iter.level, iter.xIdx(), iter.yIdx());
					f(iter);
				}
			});
			source.precede(task);
			task.precede(target);
		}
		RealTaskManager::tf.wait_for_all();
	}


	template <typename I, typename F>
	static void parallelIterateGridWithCoreIdx(std::vector<I> iterators, const F &f) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();

		if (iterators.size() != coreNum) {
			throw std::runtime_error("parallelIterateGridWithCoreIdx:: Num of iterators should be same with num of process");
		}
		for (int i = 0; i < coreNum; i++) {
			I iter = iterators[i];
			auto task = TaskManager::tf.silent_emplace([iter, &f, i]() mutable {
				for (; iter.isValid(); iter.next()) {
					//spdlog::info("iter idx{} level{} x{} y{}", iter.index(), iter.level, iter.xIdx(), iter.yIdx());
					f(iter, i);
				}
			});
			source.precede(task);
			task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}


	template<typename F>
	static void taskParallelFor(int start, int end, const F &func, int slice=0) {
		auto source = RealTaskManager::tf.placeholder();
		auto target = RealTaskManager::tf.placeholder();
		int coreNum = RealTaskManager::tf.num_workers();

		// Size of a slice for the range functions
		int n = end - start + 1;
		if(slice == 0) { slice = (int) std::round(n / static_cast<double> (coreNum)); }
		slice = std::max(slice, 1);

		int i1 = start;
		int i2 = std::min(start + slice, end);

		for (int i = 0; i1 < end; i++) {
			auto task = RealTaskManager::tf.silent_emplace([i1, i2, func]() mutable {
				for (int k = i1; k < i2; k++) {
						func(k);
				}
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, end);
		}
		if (i1 < end) {
			auto task = RealTaskManager::tf.silent_emplace([i1, end, func]() mutable {
				for (int k = i1; k < end; k++) {
						func(k);
				}
			});
			source.precede(task); task.precede(target);
		}
		RealTaskManager::tf.wait_for_all();
	}
	

	// I: OctreeGridIterator2 / OctreeOmniGridIterator2
	template <typename I, typename F>
	static void parallelIterateGrid(I iter, const F &f) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		auto tmp = iter;
		int coreNum = TaskManager::tf.num_workers();
		//spdlog::info("core num {}", coreNum);
		int nodeNum = iter.totalSize();
		//spdlog::info("node num {}", nodeNum);
		int size = (nodeNum + coreNum - 1) / coreNum;
		for (int i = 0; i < coreNum; i++) {
			auto task = TaskManager::tf.silent_emplace([iter, size, f]() mutable {
				for (int j = 0; j < size && iter.isValid(); j++, iter.next()) {
					//spdlog::info("iter idx{} level{} x{} y{}", iter.index(), iter.level, iter.xIdx(), iter.yIdx());
					f(iter);
				}
			});
			source.precede(task);
			task.precede(target);
			for (int j = 0; j < size && tmp.isValid(); j++, tmp.next()) {
				//spdlog::info("tmp idx{} level{} x{} y{}", tmp.index(), tmp.level, tmp.xIdx(), tmp.yIdx());
				;
			}
			iter = tmp;
		}
		TaskManager::tf.wait_for_all();
	}
	
	// Iter: std iterator
	template <typename I, typename F>
	static void parallelIterate(I iter, I end, int nodeNum, const F &f) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();
		int size = (nodeNum + coreNum - 1) / coreNum;
		int c = 0;
		for (int i = 0; i < coreNum; i++) {
			auto task = TaskManager::tf.silent_emplace([iter, end, size, f, c]() mutable {
				for (int j = 0; j < size && iter != end; j++, c++, iter++) {
					//spdlog::info("iter idx{} level{} x{} y{}", iter.index(), iter.level, iter.xIdx(), iter.yIdx());
					f(iter, c);
				}
			});
			source.precede(task);
			task.precede(target);
			for (int j = 0; j < size && iter != end; j++, c++, iter++) {
				//spdlog::info("tmp idx{} level{} x{} y{}", tmp.index(), tmp.level, tmp.xIdx(), tmp.yIdx());
				;
			}
		}
		TaskManager::tf.wait_for_all();
	}

	// I: OctreeGridIterator2 / OctreeOmniGridIterator2
	// if (slice == 0) slice = totalnum / corenum
	template <typename I, typename F>
	static void parallelIterateGrid(I iter, int slice, const F &f) {
		auto source = RealTaskManager::tf.placeholder();
		auto target = RealTaskManager::tf.placeholder();
		auto tmp = iter;
		int coreNum = RealTaskManager::tf.num_workers();
		//spdlog::info("core num {}", coreNum);
		int nodeNum = iter.totalSize();
		//spdlog::info("node num {}", nodeNum);
		int size = (slice > 0) ? slice : (nodeNum + coreNum - 1) / coreNum;
		int n = (nodeNum + size - 1) / size;
		for (int i = 0; i < n; i++) {
			auto task = RealTaskManager::tf.silent_emplace([iter, size, f]() mutable {
				for (int j = 0; j < size && iter.isValid(); j++, iter.next()) {
					f(iter);
				}
			});
			source.precede(task);
			task.precede(target);
			for (int j = 0; j < size && tmp.isValid(); j++, tmp.next()) {
				;
			}
			iter = tmp;
		}
		RealTaskManager::tf.wait_for_all();
	}

	template <typename I, typename F>
	static void parallelIterateGridWithCoreIdx(I iter, const F &f) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		I tmp = iter;
		int coreNum = TaskManager::tf.num_workers();
		//spdlog::info("core num {}", coreNum);
		int nodeNum = iter.totalSize();
		//spdlog::info("node num {}", nodeNum);
		int size = (nodeNum + coreNum - 1) / coreNum;
		for (int i = 0; i < coreNum; i++) {
			auto task = TaskManager::tf.silent_emplace([iter, size, &f, i]() mutable {
				for (int j = 0; j < size && iter.isValid(); j++, iter.next()) {
					//spdlog::info("iter idx{} level{} x{} y{}", iter.index(), iter.level, iter.xIdx(), iter.yIdx());
					f(iter, i);
				}
			});
			source.precede(task);
			task.precede(target);
			for (int j = 0; j < size && tmp.isValid(); j++, tmp.next()) {
				//spdlog::info("tmp idx{} level{} x{} y{}", tmp.index(), tmp.level, tmp.xIdx(), tmp.yIdx());
				;
			}
			iter = tmp;
		}
		TaskManager::tf.wait_for_all();
	}


	template<typename Index, typename F>
	static void parallelForTF(Index start, Index end, const F &func) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();

		// Size of a slice for the range functions
		Index n = end - start + 1;
		Index slice = (Index) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, Index(1));
		
		Index i1 = start;
		Index i2 = std::min(start + slice, end);

		for (int i = 0; i+1 < coreNum && i1 < end; i++) {
			auto task = TaskManager::tf.silent_emplace([i1, i2, &func]() mutable {
				for (Index k = i1; k < i2; k++) {
					func(k);
				}
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, end);
		}
		if (i1 < end) {
			auto task = TaskManager::tf.silent_emplace([i1, end, &func]() mutable {
				for (Index k = i1; k < end; k++) {
					func(k);
				}
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}


	template<typename Element, typename F>
	static void parallelForTF(const std::vector<Element>& arr, const F &func) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();

		// Size of a slice for the range functions
		auto end = arr.size();
		auto n = end + 1;
		auto slice = (size_t) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, size_t(1));
		
		size_t i1 = 0;
		auto i2 = std::min(slice, end);

		for (size_t i = 0; i+1 < coreNum && i1 < end; i++) {
			auto task = TaskManager::tf.silent_emplace([i1, i2, &arr, &func]() mutable {
				for (size_t k = i1; k < i2; k++) {
					func(arr[k]);
				}
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, end);
		}
		if (i1 < end) {
			auto task = TaskManager::tf.silent_emplace([i1, end, &arr, &func]() mutable {
				for (size_t k = i1; k < end; k++) {
					func(arr[k]);
				}
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}

	template<typename Index, typename F>
	static void parallelForPlainTF(Index start, Index end, const F &func) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();

		// Size of a slice for the range functions
		Index n = end - start + 1;
		Index slice = (Index) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, Index(1));
		
		Index i1 = start;
		Index i2 = std::min(start + slice, end);

		int i = 0;
		for (; i+1 < coreNum && i1 < end; i++) {
			auto task = TaskManager::tf.silent_emplace([i1, i2, i, &func]() mutable {
				func(i1, i2, i);
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, end);
		}
		if (i1 < end) {
			auto task = TaskManager::tf.silent_emplace([i1, end, i, &func]() mutable {
				func(i1, end, i);
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}

	template<int d, typename F>
	static void parallelForBoxIteratorTF(const Eigen::Matrix<int, d, 1>& min_coor,
		const Eigen::Matrix<int, d, 1>& max_coor, const F &func) {
		using IndexT = Eigen::Matrix<int, d, 1>;
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();

		int start = min_coor[d-1];
		int end = max_coor[d-1] + 1;
		// Size of a slice for the range functions
		int n = end - start + 1;
		int slice = (int) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, int(1));
		
		int i1 = start;
		int i2 = std::min(start + slice, end);

		for (int i = 0; i+1 < coreNum && i1 < end; i++) {
			auto task = TaskManager::tf.silent_emplace([i1, i2, min_coor, max_coor, &func]() mutable {
				IndexT first_coor = min_coor, last_coor = max_coor; 
				first_coor(d-1) = i1; last_coor(d-1) = i2-1;
				for(auto iter=ast::BoxIterator<d>(first_coor, last_coor); iter.valid(); iter.next())
					func(iter.index());
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, end);
		}
		if (i1 < end) {
			auto task = TaskManager::tf.silent_emplace([i1, min_coor, max_coor, &func]() mutable {
				IndexT first_coor = min_coor; first_coor(d-1) = i1;
				for(auto iter=ast::BoxIterator<d>(first_coor, max_coor); iter.valid(); iter.next())
					func(iter.index());
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}

	template<typename Index, typename F>
	static void parallelForTFWithCoreID(Index start, Index end, const F &func) {
		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();

		// Size of a slice for the range functions
		Index n = end - start + 1;
		Index slice = (Index) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, Index(1));
		
		Index i1 = start;
		Index i2 = std::min(start + slice, end);

		int i = 0;
		for (; i+1 < coreNum && i1 < end; i++) {
			auto task = TaskManager::tf.silent_emplace([i, i1, i2, &func]() mutable {
				for (Index k = i1; k < i2; k++) {
					func(k, i);
				}
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, end);
		}
		if (i1 < end) {
			auto task = TaskManager::tf.silent_emplace([i, i1, end, &func]() mutable {
				for (Index k = i1; k < end; k++) {
					func(k, i);
				}
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();
	}

	template<typename Index, typename T>
	static T reduceMax(T* start, Index size) {
		ASSERT(size > 0, "at least one element");
		if(size < 100) { T m = start[0]; for(Index i = 0; i < size; ++i) { m = std::max(m, start[i]); } return m; }

		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();
		
		std::vector<T> ele(coreNum);
		Index n = size - 0 + 1;
		Index slice = (Index) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, Index(1));
		
		Index i1 = 0;
		Index i2 = std::min(0 + slice, size);

		int i = 0;
		for (; i+1 < coreNum && i1 < size; i++) {
			auto task = TaskManager::tf.silent_emplace([i, i1, i2, &start, &ele]() mutable {
				T m = start[i1];
				for (Index k = i1; k < i2; k++) { m = std::max(m, start[k]); }
				ele[i] = m;
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, size);
		}
		if (i1 < size) {
			auto task = TaskManager::tf.silent_emplace([i, i1, size, &start, &ele]() mutable {
				T m = start[i1];
				for (Index k = i1; k < size; k++) { m = std::max(m, start[k]); }
				ele[i] = m;
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();

		T m = ele[0];
		for(auto& e: ele) {
			m = std::max(m, e);
		}
		return m;
	}

	template<typename T, typename Func>
	static auto reduceTransformMax(const std::vector<T>& value, const Func& f) {
		int size = value.size();
		ASSERT(size > 0, "at least one element");
		typedef decltype(f(std::declval<T>())) ReturnType;
		if(size < 100) { ReturnType m = f(value[0]); for(int i = 0; i < size; ++i) { m = std::max(m, f(value[i])); } return m; }

		auto source = TaskManager::tf.placeholder();
		auto target = TaskManager::tf.placeholder();
		int coreNum = TaskManager::tf.num_workers();
		
		std::vector<ReturnType> ele(coreNum);
		int n = size - 0 + 1;
		int slice = (int) std::round(n / static_cast<double> (coreNum));
		slice = std::max(slice, int(1));
		
		int i1 = 0;
		int i2 = std::min(0 + slice, size);

		int i = 0;
		for (; i+1 < coreNum && i1 < size; i++) {
			auto task = TaskManager::tf.silent_emplace([i, i1, i2, &value, &ele, &f]() mutable {
				ReturnType m = f(value[i1]);
				for (int k = i1; k < i2; k++) { m = std::max(m, f(value[k])); }
				ele[i] = m;
			});
			source.precede(task); task.precede(target);
			i1 = i2;
			i2 = std::min(i2+slice, size);
		}
		if (i1 < size) {
			auto task = TaskManager::tf.silent_emplace([i, i1, size, &value, &ele, &f]() mutable {
				ReturnType m = f(value[i1]);
				for (int k = i1; k < size; k++) { m = std::max(m, f(value[k])); }
				ele[i] = m;
			});
			source.precede(task); task.precede(target);
		}
		TaskManager::tf.wait_for_all();

		ReturnType m = ele[0];
		for(auto& e: ele) { m = std::max(m, e); }
		return m;
	}

	// template<typename Index, typename Callable>
	// static void parallelFor(Index start, Index end, Callable func) {
	// 	// Estimate number of threads in the pool
	// 	const static unsigned nb_threads_hint = std::thread::hardware_concurrency();
	// 	const static unsigned nb_threads = (nb_threads_hint == 0u ? 8u : nb_threads_hint);
	// 	// Size of a slice for the range functions
	// 	Index n = end - start + 1;
	// 	Index slice = (Index) std::round(n / static_cast<double> (nb_threads));
	// 	slice = std::max(slice, Index(1));

	// 	// [Helper] Inner loop
	// 	auto launchRange = [&func] (int k1, int k2) {
	// 		for (Index k = k1; k < k2; k++) {
	// 			func(k);
	// 		}
	// 	};

	// 	// Create pool and launch jobs
	// 	std::vector<std::thread> pool;
	// 	pool.reserve(nb_threads);
	// 	Index i1 = start;
	// 	Index i2 = std::min(start + slice, end);
	// 	for (unsigned i = 0; i + 1 < nb_threads && i1 < end; ++i) {
	// 		pool.emplace_back(launchRange, i1, i2);
	// 		i1 = i2;
	// 		i2 = std::min(i2 + slice, end);
	// 	}
	// 	if (i1 < end) {
	// 		pool.emplace_back(launchRange, i1, end);
	// 	}

	// 	// Wait for jobs to finish
	// 	for (std::thread &t : pool) {
	// 		if (t.joinable()) {
	// 			t.join();
	// 		}
	// 	}
	// }

	// Serial version for easy comparison
	// template<typename Index, typename Callable>
	// static void sequentialFor(Index start, Index end, Callable func) {
	// 	for (Index i = start; i < end; i++) {
	// 		func(i);
	// 	}
	// }

	template <typename Index, typename Function>
	static void parallelForTF(Index beginIndexX, Index endIndexX,
					Index beginIndexY, Index endIndexY,
					const Function& function) {
		parallelForTF(beginIndexY, endIndexY,
			[&](Index y) {
				for (Index i = beginIndexX; i < endIndexX; ++i) {
					function(i, y);
				}
		});
	}

	template <typename Index, typename Function>
	static void parallelForTF(Index beginIndexX, Index endIndexX,
					Index beginIndexY, Index endIndexY,
					Index beginIndexZ, Index endIndexZ,
					const Function& function) {
		parallelForTF(beginIndexZ, endIndexZ,
			[&](Index k) {
				for (Index j = beginIndexY; j < endIndexY; ++j) {
					for (Index i = beginIndexX; i < endIndexX; ++i) {
						function(i, j, k);
					}
			}
		});
	}

};

////////////////////////////////////////////////////////////////////////////////
// #include <iostream>
// #include <mutex>
////////////////////////////////////////////////////////////////////////////////

// int main() {
// 	std::mutex critical;
// 	ThreadPool::ParallelFor(0, 16, [&] (int i) {
// 		std::lock_guard<std::mutex> lock(critical);
// 		std::cout << i << std::endl;
// 	});
// 	return 0;
// }

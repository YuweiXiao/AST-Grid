#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <iostream>
#include <vector>
#include "util.h"
#include "macro.h"
#include "general.h"
#include "thread_pool.h"

//============================================================================
// Dynamic compressed sparse row matrix.
namespace PCG {
	template<class T>
	struct SparseMatrix
	{
		unsigned int n; // dimension
		std::vector<std::vector<unsigned int> > index; // for each row, a list of all column indices (sorted)
		std::vector<std::vector<T> > value; // values corresponding to index

		explicit SparseMatrix(unsigned int n_ = 0, unsigned int expected_nonzeros_per_row = 7)
			: n(n_), index(n_), value(n_)
		{
			for (unsigned int i = 0; i < n; ++i) {
				index[i].reserve(expected_nonzeros_per_row);
				value[i].reserve(expected_nonzeros_per_row);
			}
		}

		void clear(void)
		{
			n = 0;
			index.clear();
			value.clear();
		}

		void zero(void)
		{
			for (unsigned int i = 0; i < n; ++i) {
				index[i].resize(0);
				value[i].resize(0);
			}
		}

		void resize(int n_)
		{
			n = n_;
			index.resize(n);
			value.resize(n);
		}

		void scale(T s) {
			for(int i = 0; i < n; ++i) {
				for(int k = 0; k < value[i].size(); ++k) {
					value[i][k] *= s;
				}
			}
		}

		T operator()(unsigned int i, unsigned int j) const
		{
			for (unsigned int k = 0; k < index[i].size(); ++k) {
				if (index[i][k] == j) return value[i][k];
				else if (index[i][k] > j) return 0;
			}
			return 0;
		}

		void set_element(unsigned int i, unsigned int j, T new_value)
		{
			unsigned int k = 0;
			for (; k < index[i].size(); ++k) {
				if (index[i][k] == j) {
					value[i][k] = new_value;
					return;
				}
				else if (index[i][k] > j) {
					insert(index[i], k, j);
					insert(value[i], k, new_value);
					return;
				}
			}
			index[i].push_back(j);
			value[i].push_back(new_value);
		}

		void add_to_element(unsigned int i, unsigned int j, T increment_value)
		{
			unsigned int k = 0;
			for (; k < index[i].size(); ++k) {
				if (index[i][k] == j) {
					value[i][k] += increment_value;
					return;
				}
				else if (index[i][k] > j) {
					insert(index[i], k, j);
					insert(value[i], k, increment_value);
					return;
				}
			}
			index[i].push_back(j);
			value[i].push_back(increment_value);
		}

		// assumes indices is already sorted
		void add_sparse_row(unsigned int i, const std::vector<unsigned int> &indices, const std::vector<T> &values)
		{
			unsigned int j = 0, k = 0;
			while (j < indices.size() && k < index[i].size()) {
				if (index[i][k] < indices[j]) {
					++k;
				}
				else if (index[i][k] > indices[j]) {
					insert(index[i], k, indices[j]);
					insert(value[i], k, values[j]);
					++j;
				}
				else {
					value[i][k] += values[j];
					++j;
					++k;
				}
			}
			for (; j < indices.size(); ++j) {
				index[i].push_back(indices[j]);
				value[i].push_back(values[j]);
			}
		}

		// assumes matrix has symmetric structure - so the indices in row i tell us which columns to delete i from
		void symmetric_remove_row_and_column(unsigned int i)
		{
			for (unsigned int a = 0; a < index[i].size(); ++a) {
				unsigned int j = index[i][a]; // 
				for (unsigned int b = 0; b < index[j].size(); ++b) {
					if (index[j][b] == i) {
						erase(index[j], b);
						erase(value[j], b);
						break;
					}
				}
			}
			index[i].resize(0);
			value[i].resize(0);
		}

		void write_matlab(std::ostream &output, const char *variable_name)
		{
			output << variable_name << "=sparse([";
			for (unsigned int i = 0; i < n; ++i) {
				for (unsigned int j = 0; j < index[i].size(); ++j) {
					output << i + 1 << " ";
				}
			}
			output << "],...\n  [";
			for (unsigned int i = 0; i < n; ++i) {
				for (unsigned int j = 0; j < index[i].size(); ++j) {
					output << index[i][j] + 1 << " ";
				}
			}
			output << "],...\n  [";
			for (unsigned int i = 0; i < n; ++i) {
				for (unsigned int j = 0; j < value[i].size(); ++j) {
					output << value[i][j] << " ";
				}
			}
			output << "], " << n << ", " << n << ");" << std::endl;
		}


		void printReadable(std::ostream& out) const {
			for (unsigned int i = 0; i < n; i++) {
				for (unsigned int j = 0; j < index[i].size(); j++) {	
					for(unsigned int k = (j == 0 ? 0 : index[i][j-1]+1); k < index[i][j]; ++k) 
						out<<0<<' ';
					out<<value[i][j]<<' ';
				}
				for(unsigned int k = (index[i].size() > 0 ? index[i][index[i].size()-1]+1 : 0) ; k < n; ++k)
					out<<0<<' ';
				out<<std::endl;
			}
		}
	};

	typedef SparseMatrix<float> SparseMatrixf;
	typedef SparseMatrix<double> SparseMatrixd;

	template<class T>
	bool isSymmetrical(const PCG::SparseMatrix<T> &A) {
		bool ret = true;
		int n = A.n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < A.index[i].size(); j++) {	
				int k = A.index[i][j];
				if (std::fabs(A(i, k) - A(k, i) ) > 1e-6) {
					spdlog::info("row:{} col:{}, value(i, j): {}, value(j, i): {}", i, k, A(i, k), A(k, i));
					ret = false;
				}
			}
		}
		return ret;
	}

	// perform result=matrix*x
	template<class T>
	void multiply(const SparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
	{
		assert(matrix.n == x.size());
		result.resize(matrix.n);
#ifdef USE_OPENMP
#pragma omp parallel for
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T temp = 0;
			for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
				temp += matrix.value[i][j] * x[matrix.index[i][j]];
			}
			result[i] = temp;
		}
#else
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T temp = 0;
			for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
				temp += matrix.value[i][j] * x[matrix.index[i][j]];
			}
			result[i] = temp;
		}
#endif
	}
	
	// perform result=matrix*x
	template<class T>
	void multiply(const SparseMatrix<T> &matrix, const Eigen::Matrix<T, -1, 1> &x, Eigen::Matrix<T, -1, 1> &result)
	{
		ASSERT(matrix.n == x.rows(), "matrix dimension should be same");
		if(result.rows() != matrix.n) {
			result = Eigen::Matrix<T, -1, 1>(matrix.n);
		}
#ifdef USE_OPENMP
#pragma omp parallel for
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T temp = 0;
			for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
				temp += matrix.value[i][j] * x[matrix.index[i][j]];
			}
			result[i] = temp;
		}
#else
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T temp = 0;
			for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
				temp += matrix.value[i][j] * x[matrix.index[i][j]];
			}
			result[i] = temp;
		}
#endif
	}

	// perform result=result-matrix*x
	template<class T>
	void multiply_and_subtract(const SparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
	{
		assert(matrix.n == x.size());
		result.resize(matrix.n);
		for (unsigned int i = 0; i < matrix.n; ++i) {
			for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
				result[i] -= matrix.value[i][j] * x[matrix.index[i][j]];
			}
		}
	}

	//============================================================================
	// Fixed version of SparseMatrix. This is not a good structure for dynamically
	// modifying the matrix, but can be significantly faster for matrix-vector
	// multiplies due to better data locality.

	template<class T>
	struct FixedSparseMatrix
	{
		unsigned int n; // dimension
		std::vector<T> value; // nonzero values row by row
		std::vector<int> colindex; // corresponding column indices
		std::vector<int> rowstart; // where each row starts in value and colindex (and last entry is one past the end, the number of nonzeros)

		explicit FixedSparseMatrix(unsigned int n_ = 0)
			: n(n_), value(0), colindex(0), rowstart(n_ + 1)
		{}

		void clear(void)
		{
			n = 0;
			value.clear();
			colindex.clear();
			rowstart.clear();
		}

		void resize(int n_)
		{
			n = n_;
			rowstart.resize(n + 1);
		}

		void construct_from_matrix(const SparseMatrix<T> &matrix)
		{
			resize(matrix.n);
			rowstart[0] = 0;
			for (unsigned int i = 0; i < n; ++i) {
				rowstart[i + 1] = rowstart[i] + matrix.index[i].size();
			}
			value.resize(rowstart[n]);
			colindex.resize(rowstart[n]);
			ThreadPool::parallelForTF(0u, n, [&](unsigned int i){
				for(unsigned int k = 0; k < matrix.index[i].size(); ++k) {
					value[rowstart[i]+k] = matrix.value[i][k];
					colindex[rowstart[i]+k] = matrix.index[i][k];
				}
			});
			// unsigned int j = 0;
			// for (unsigned int i = 0; i < n; ++i) {
			// 	for (unsigned int k = 0; k < matrix.index[i].size(); ++k) {
			// 		value[j] = matrix.value[i][k];
			// 		colindex[j] = matrix.index[i][k];
			// 		++j;
			// 	}
			// }
		}

		void write_matlab(std::ostream &output, const char *variable_name)
		{
			output << variable_name << "=sparse([";
			for (unsigned int i = 0; i < n; ++i) {
				for (unsigned int j = rowstart[i]; j < rowstart[i + 1]; ++j) {
					output << i + 1 << " ";
				}
			}
			output << "],...\n  [";
			for (unsigned int i = 0; i < n; ++i) {
				for (unsigned int j = rowstart[i]; j < rowstart[i + 1]; ++j) {
					output << colindex[j] + 1 << " ";
				}
			}
			output << "],...\n  [";
			for (unsigned int i = 0; i < n; ++i) {
				for (unsigned int j = rowstart[i]; j < rowstart[i + 1]; ++j) {
					output << value[j] << " ";
				}
			}
			output << "], " << n << ", " << n << ");" << std::endl;
		}
	};

	typedef FixedSparseMatrix<float> FixedSparseMatrixf;
	typedef FixedSparseMatrix<double> FixedSparseMatrixd;

	// perform result=matrix*x
	template<class T>
	void multiply(const FixedSparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
	{
		assert(matrix.n == x.size());
		result.resize(matrix.n);
#ifdef USE_OPENMP
#pragma omp parallel for
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T t = 0;
			for (unsigned int j = matrix.rowstart[i]; j < matrix.rowstart[i + 1]; ++j) {
				t += matrix.value[j] * x[matrix.colindex[j]];
			}
			result[i] = t;
		}
#else
		ThreadPool::parallelForTF(0u, matrix.n, [&](unsigned int i){
			T t = 0;
			for (unsigned int j = matrix.rowstart[i]; j < matrix.rowstart[i + 1]; ++j) {
				t += matrix.value[j] * x[matrix.colindex[j]];
			}
			result[i] = t;
		});
#endif
	}
	
	// perform result=matrix*x
	template<class T>
	void multiply(const FixedSparseMatrix<T> &matrix, const Eigen::Matrix<T, -1, 1> &x,  Eigen::Matrix<T, -1, 1> &result)
	{
		assert(matrix.n == x.rows());
		if(matrix.n != result.rows()) {
			result = Eigen::Matrix<T, -1, 1>(matrix.n);
		}
#ifdef USE_OPENMP
#pragma omp parallel for
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T t = 0;
			for (unsigned int j = matrix.rowstart[i]; j < matrix.rowstart[i + 1]; ++j) {
				t += matrix.value[j] * x[matrix.colindex[j]];
			}
			result[i] = t;
		}
#elif defined(USE_TASKFLOW)
		ThreadPool::parallelForTF(0, (int)matrix.n, [&](int i) {
			T t = 0;
			for (unsigned int j = matrix.rowstart[i]; j < matrix.rowstart[i + 1]; ++j) {
				t += matrix.value[j] * x[matrix.colindex[j]];
			}
			result[i] = t;
		});
#else
		for (unsigned int i = 0; i < matrix.n; ++i) {
			T t = 0;
			for (unsigned int j = matrix.rowstart[i]; j < matrix.rowstart[i + 1]; ++j) {
				t += matrix.value[j] * x[matrix.colindex[j]];
			}
			result[i] = t;
		}
#endif
	}

	// perform result=result-matrix*x
	template<class T>
	void multiply_and_subtract(const FixedSparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
	{
		assert(matrix.n == x.size());
		result.resize(matrix.n);
		for (unsigned int i = 0; i < matrix.n; ++i) {
			for (unsigned int j = matrix.rowstart[i]; j < matrix.rowstart[i + 1]; ++j) {
				result[i] -= matrix.value[j] * x[matrix.colindex[j]];
			}
		}
	}

	// result = this*x - diagonal(this) * x
    template<class T>
    void multiply_non_diagonal(const FixedSparseMatrix<T> &matrix, const std::vector<T>& x, std::vector<T>& result) {
        ASSERT(matrix.n == x.size(), "matrix vector multiplication dimension mismatch");
        result.resize(matrix.n);
		ThreadPool::parallelForTF(0u, matrix.n, [&](uint i){
        // for(uint i = 0; i < matrix.n; ++i) {
            T tmp = 0;
            for(uint j = matrix.rowstart[i]; j < matrix.rowstart[i+1]; ++j) {
                if(matrix.colindex[j] != i) {
                    tmp += matrix.value[j] * x[matrix.colindex[j]];
                }
            }
            result[i] = tmp;
        });
    }


    // result = inv(diagonal(this)) * x
    template<class T>
    void multiply_inv_diagonal(const FixedSparseMatrix<T> &matrix, const std::vector<T>& x, std::vector<T>& result) {
        ASSERT(matrix.n == x.size(), "matrix vector multiplication dimension mismatch");
        result.resize(matrix.n);
		ThreadPool::parallelForTF(0u, matrix.n, [&](uint i){
        // for(uint i = 0; i < matrix.n; ++i) {
			uint j = matrix.rowstart[i];
            for(; j < matrix.rowstart[i+1]; ++j) {
                if(matrix.colindex[j] == i) {
                    result[i] = x[i] / matrix.value[j];
                    break;
                }
            }
			if(j == matrix.rowstart[i+1]+1) {
				throw std::runtime_error("diagonal entry does not exist");
			}
        });
    }
}

#endif

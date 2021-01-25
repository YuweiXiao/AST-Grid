#ifndef PCG_SOLVER_H
#define PCG_SOLVER_H

// Implements PCG with Modified Incomplete Cholesky (0) preconditioner.
// PCGSolver<T> is the main class for setting up and solving a linear system.
// Note that this only handles symmetric positive (semi-)definite matrices,
// with guarantees made only for M-matrices (where off-diagonal entries are all
// non-positive, and row sums are non-negative).

#include <cmath>
#include "sparse_matrix.h"
#include "blas_wrapper.h"
#include "timer.h"
#include "global_benchmark.h"

#define JACOBIAN_P 0

//============================================================================
// A simple compressed sparse column data structure (with separate diagonal)
// for lower triangular matrices
using namespace PCG;
template<class T>
struct SparseColumnLowerFactor
{
   unsigned int n;
   std::vector<T> invdiag; // reciprocals of diagonal elements
   std::vector<T> value; // values below the diagonal, listed column by column
   std::vector<unsigned int> rowindex; // a list of all row indices, for each column in turn
   std::vector<unsigned int> colstart; // where each column begins in rowindex (plus an extra entry at the end, of #nonzeros)
   std::vector<T> adiag; // just used in factorization: minimum "safe" diagonal entry allowed

   explicit SparseColumnLowerFactor(unsigned int n_=0)
      : n(n_), invdiag(n_), colstart(n_+1), adiag(n_)
   {}

   void clear(void)
   {
      n=0;
      invdiag.clear();
      value.clear();
      rowindex.clear();
      colstart.clear();
      adiag.clear();
   }

   void resize(unsigned int n_)
   {
      n=n_;
      invdiag.resize(n);
      colstart.resize(n+1);
      adiag.resize(n);
   }

   void write_matlab(std::ostream &output, const char *variable_name)
   {
      output<<variable_name<<"=sparse([";
      for(unsigned int i=0; i<n; ++i){
         output<<" "<<i+1;
         for(unsigned int j=colstart[i]; j<colstart[i+1]; ++j){
            output<<" "<<rowindex[j]+1;
         }
      }
      output<<"],...\n  [";
      for(unsigned int i=0; i<n; ++i){
         output<<" "<<i+1;
         for(unsigned int j=colstart[i]; j<colstart[i+1]; ++j){
            output<<" "<<i+1;
         }
      }
      output<<"],...\n  [";
      for(unsigned int i=0; i<n; ++i){
         output<<" "<<(invdiag[i]!=0 ? 1/invdiag[i] : 0);
         for(unsigned int j=colstart[i]; j<colstart[i+1]; ++j){
            output<<" "<<value[j];
         }
      }
      output<<"], "<<n<<", "<<n<<");"<<std::endl;
   }
};

//============================================================================
// Incomplete Cholesky factorization, level zero, with option for modified version.
// Set modification_parameter between zero (regular incomplete Cholesky) and
// one (fully modified version), with values close to one usually giving the best
// results. The min_diagonal_ratio parameter is used to detect and correct
// problems in factorization: if a pivot is this much less than the diagonal
// entry from the original matrix, the original matrix entry is used instead.

template<class T>
void factor_modified_incomplete_cholesky0(const PCG::SparseMatrix<T> &matrix, SparseColumnLowerFactor<T> &factor,
                                          T modification_parameter=0.97, T min_diagonal_ratio=0.25)
{
   // first copy lower triangle of matrix into factor (Note: assuming A is symmetric of course!)
   factor.resize(matrix.n);
   zero(factor.invdiag); // important: eliminate old values from previous solves!
   factor.value.resize(0);
   factor.rowindex.resize(0);
   zero(factor.adiag);
   for(unsigned int i=0; i<matrix.n; ++i){
      factor.colstart[i]=(unsigned int)factor.rowindex.size();
      for(unsigned int j=0; j<matrix.index[i].size(); ++j){
         if(matrix.index[i][j]>i){
            factor.rowindex.push_back(matrix.index[i][j]);
            factor.value.push_back(matrix.value[i][j]);
         }else if(matrix.index[i][j]==i){
            factor.invdiag[i]=factor.adiag[i]=matrix.value[i][j];
         }
      }
   }
   factor.colstart[matrix.n]=(unsigned int)factor.rowindex.size();
   // now do the incomplete factorization (figure out numerical values)

   // MATLAB code:
   // L=tril(A);
   // for k=1:size(L,2)
   //   L(k,k)=sqrt(L(k,k));
   //   L(k+1:end,k)=L(k+1:end,k)/L(k,k);
   //   for j=find(L(:,k))'
   //     if j>k
   //       fullupdate=L(:,k)*L(j,k);
   //       incompleteupdate=fullupdate.*(A(:,j)~=0);
   //       missing=sum(fullupdate-incompleteupdate);
   //       L(j:end,j)=L(j:end,j)-incompleteupdate(j:end);
   //       L(j,j)=L(j,j)-omega*missing;
   //     end
   //   end
   // end

   for(unsigned int k=0; k<matrix.n; ++k){
      if(factor.adiag[k]==0) continue; // null row/column
      // figure out the final L(k,k) entry
      if(factor.invdiag[k]<min_diagonal_ratio*factor.adiag[k])
         factor.invdiag[k]=1/sqrt(factor.adiag[k]); // drop to Gauss-Seidel here if the pivot looks dangerously small
      else
         factor.invdiag[k]=1/sqrt(factor.invdiag[k]);
      // finalize the k'th column L(:,k)
      for(unsigned int p=factor.colstart[k]; p<factor.colstart[k+1]; ++p){
         factor.value[p]*=factor.invdiag[k];
      }
      // incompletely eliminate L(:,k) from future columns, modifying diagonals
      for(unsigned int p=factor.colstart[k]; p<factor.colstart[k+1]; ++p){
         unsigned int j=factor.rowindex[p]; // work on column j
         T multiplier=factor.value[p];
         T missing=0;
         unsigned int a=factor.colstart[k];
         // first look for contributions to missing from dropped entries above the diagonal in column j
         unsigned int b=0;
         while(a<factor.colstart[k+1] && factor.rowindex[a]<j){
            // look for factor.rowindex[a] in matrix.index[j] starting at b
            while(b<matrix.index[j].size()){
               if(matrix.index[j][b]<factor.rowindex[a])
                  ++b;
               else if(matrix.index[j][b]==factor.rowindex[a])
                  break;
               else{
                  missing+=factor.value[a];
                  break;
               }
            }
            ++a;
         }
         // adjust the diagonal j,j entry
         if(a<factor.colstart[k+1] && factor.rowindex[a]==j){
            factor.invdiag[j]-=multiplier*factor.value[a];
         }
         ++a;
         // and now eliminate from the nonzero entries below the diagonal in column j (or add to missing if we can't)
         b=factor.colstart[j];
         while(a<factor.colstart[k+1] && b<factor.colstart[j+1]){
            if(factor.rowindex[b]<factor.rowindex[a])
               ++b;
            else if(factor.rowindex[b]==factor.rowindex[a]){
               factor.value[b]-=multiplier*factor.value[a];
               ++a;
               ++b;
            }else{
               missing+=factor.value[a];
               ++a;
            }
         }
         // and if there's anything left to do, add it to missing
         while(a<factor.colstart[k+1]){
            missing+=factor.value[a];
            ++a;
         }
         // and do the final diagonal adjustment from the missing entries
         factor.invdiag[j]-=modification_parameter*multiplier*missing;
      }
   }
}

//============================================================================
// Solution routines with lower triangular matrix.

// solve L*result=rhs
template<class T>
void solve_lower(const SparseColumnLowerFactor<T> &factor, const std::vector<T> &rhs, std::vector<T> &result)
{
   ASSERT(factor.n==rhs.size(), "size should be same");
   ASSERT(factor.n==result.size(), "size should be same");
   result=rhs;
   for(unsigned int i=0; i<factor.n; ++i){
      result[i]*=factor.invdiag[i];
      for(unsigned int j=factor.colstart[i]; j<factor.colstart[i+1]; ++j){
         result[factor.rowindex[j]]-=factor.value[j]*result[i];
      }
   }
}

template<class T>
void solve_lower(const SparseColumnLowerFactor<T> &factor, const Eigen::Matrix<T, -1, 1> &rhs, Eigen::Matrix<T, -1, 1> &result)
{
    ASSERT(factor.n==rhs.rows(), "rhs should same dimension"); 
    ASSERT(factor.n==result.rows(), "result should have same dimension");
    result = rhs;
    for(unsigned int i=0; i < factor.n; ++i){
        result[i] *=factor.invdiag[i];
        for(unsigned int j=factor.colstart[i]; j<factor.colstart[i+1]; ++j){
            result[factor.rowindex[j]] -= factor.value[j]*result[i];
        }
    }
}

// solve L^T*result=rhs
template<class T>
void solve_lower_transpose_in_place(const SparseColumnLowerFactor<T> &factor, std::vector<T> &x)
{
   assert(factor.n==x.size());
   assert(factor.n>0);
   unsigned int i=factor.n;
   do{
      --i;
      for(unsigned int j=factor.colstart[i]; j<factor.colstart[i+1]; ++j){
         x[i]-=factor.value[j]*x[factor.rowindex[j]];
      }
      x[i]*=factor.invdiag[i];
   }while(i!=0);
}

template<class T>
void solve_lower_transpose_in_place(const SparseColumnLowerFactor<T> &factor, Eigen::Matrix<T, -1, 1> &x)
{

   ASSERT(factor.n==x.rows() && factor.n>0, "x should have same dimension and non-negative");
   unsigned int i=factor.n;
   do{
      --i;
      for(unsigned int j=factor.colstart[i]; j<factor.colstart[i+1]; ++j){
         x[i]-=factor.value[j]*x[factor.rowindex[j]];
      }
      x[i]*=factor.invdiag[i];
   }while(i!=0);
}

//============================================================================
// Encapsulates the Conjugate Gradient algorithm with incomplete Cholesky
// factorization preconditioner.

template <class T>
struct PCGSolver
{
   PCGSolver(void)
   {
      set_solver_parameters(1e-5, 100, 0.97, 0.25);
   }

   void set_solver_parameters(T tolerance_factor_, int max_iterations_, T modified_incomplete_cholesky_parameter_=0.97, T min_diagonal_ratio_=0.25)
   {
      tolerance_factor=tolerance_factor_;
      if(tolerance_factor<1e-30) tolerance_factor=1e-30;
      max_iterations=max_iterations_;
      modified_incomplete_cholesky_parameter=modified_incomplete_cholesky_parameter_;
      min_diagonal_ratio=min_diagonal_ratio_;
   }

   void weightedJacobianSmooth(const FixedSparseMatrix<T>& matrix, std::vector<T>& u, const std::vector<T>& b) {
      const T omega = 2/3.0;
      const T oneMinusOmega = 1.0 - omega;
      // x_{k+1} = \omega D^(-1) [b-(L+U)x_k] + (1-\omega) x_{k+1}
      // multiply_non_diagonal(matrix, u, m);
      // ThreadPool::parallelForTF(0, (int)m.size(), [&](int i) {
      //    m[i] = b[i] - m[i];
      // });
      // multiply_inv_diagonal(matrix, m, m);
      // ThreadPool::parallelForTF(0, (int)m.size(), [&](int i) {
      //    u[i] = oneMinusOmega * u[i] + omega * m[i];
      // });
      ThreadPool::parallelForTF(0u, matrix.n, [&](uint i){
            T tmp = 0, diag = 0;
            for(uint j = matrix.rowstart[i]; j < matrix.rowstart[i+1]; ++j) {
               if(matrix.colindex[j] != i) {
                  tmp += matrix.value[j] * u[matrix.colindex[j]];
               } else {
                  diag = matrix.value[j];
               }
            }
            tmp = (b[i] - tmp) / diag;
            u[i] = oneMinusOmega * u[i] + omega * tmp;
      });
   }

   bool solve(const SparseMatrix<T> &matrix, const std::vector<T> &rhs, std::vector<T> &result, T &residual_out, int &iterations_out) 
   {
      unsigned int n=matrix.n;
      if(m.size()!=n)
	   { 
         m.resize(n); s.resize(n); z.resize(n); r.resize(n); 
      }
      zero(result);
      r=rhs;
      residual_out=BLAS::abs_max(r);
      if(residual_out == 0) {
         iterations_out=0;
         return true;
      }
      double tol=tolerance_factor;//*residual_out;
      {  BENCHMARK_SCOPED_TIMER_SECTION t("construct fixed matrix");
         fixed_matrix.construct_from_matrix(matrix);
      }
#if JACOBIAN_P
      apply_preconditioner_jacobian(r, z);
#else
      form_preconditioner(matrix);
      apply_preconditioner(r, z);   // v0 = C^(-1)*r
#endif
      double rho = BLAS::dot(z, r);   // rho = <v0, r>
      if(rho==0 || rho!=rho) {
         iterations_out=0;
         return false;
      }

      s=z;
      int iteration;
      {  BENCHMARK_SCOPED_TIMER_SECTION t("pcg iteration");
         for(iteration=0; iteration<max_iterations; ++iteration){
            multiply(fixed_matrix, s, z);    // A*v0
            double alpha=rho/BLAS::dot(s, z);   // t = rho / <v0, A*v0>
            ThreadPool::parallelForTF((size_t)0, s.size(), [&](size_t i){
               result[i] += alpha * s[i];
               r[i] -= alpha * z[i];
            });
            // BLAS::add_scaled(alpha, s, result); // x = x + t * v0
            // BLAS::add_scaled(-alpha, z, r);     // r = r - t * A * v0
            residual_out=BLAS::abs_max(r);
            if(residual_out<=tol) {
               iterations_out=iteration+1;
               return true; 
            }
#if JACOBIAN_P
            apply_preconditioner_jacobian(r, z);
#else
            apply_preconditioner(r, z);         // z = C^(-1) * r
#endif
            double rho_new=BLAS::dot(z, r);     // rho_new = <C^(-1) * r, r>
            double beta=rho_new/rho;            // s = rho_new / rho
            BLAS::add_scaled(beta, s, z);       // v = C^(-1) * r + s * v0
            s.swap(z);
            rho=rho_new;
         }
      }
      iterations_out=iteration;
      return false;
   }

   // eigen version of pcg solve
   bool solve(const SparseMatrix<T> &matrix, const VectorXd &rhs, VectorXd &result, T &residual_out, int &iterations_out, 
      std::function<void(int iter)> callback = nullptr) {
         Omni::Timer timer;
         // REAL a1 = 0, a2 = 0, a3 = 0;
         unsigned int n = matrix.n;
         if (em.rows() != n) {
           em.resize(n);es.resize(n);ez.resize(n);er.resize(n);
         }
         result.setZero();
         er = rhs;
         residual_out = er.cwiseAbs().maxCoeff();
         if (residual_out == 0) {
               iterations_out = 0;
               return true;
         }
         double tol = tolerance_factor * residual_out;

         spdlog::warn("\t\t[time profile] basic init time : {}", timer.durationInSeconds()); timer.reset();
         form_preconditioner(matrix);
         spdlog::warn("\t\t[time profile] form precondition time : {}", timer.durationInSeconds()); timer.reset();
         apply_preconditioner(er, ez);
         double rho = BLAS::dot(ez, er);
         if (rho == 0 || rho != rho) {
            iterations_out = 0;
            return false;
         }

         es = ez;
         TIMING(timer, fixed_matrix.construct_from_matrix(matrix), "\t\t[time profile] construct form matrix time : {}");
         int iteration;
         for (iteration = 0; iteration < max_iterations; ++iteration) {
            timer.reset();
            if(callback) callback(iteration);
            multiply(fixed_matrix, es, ez);
            // a1 += timer.durationInSeconds(); timer.reset();
            double alpha = rho / es.dot(ez);
            result += alpha * es;
            er += -alpha * ez;
            // BLAS::add_scaled(alpha, es, result);
            // BLAS::add_scaled(-alpha, ez, er);
            residual_out = er.cwiseAbs().maxCoeff();// BLAS::abs_max(er);
            // a2 += timer.durationInSeconds(); timer.reset();
            if (residual_out <= tol) {
                iterations_out = iteration + 1;
               //  spdlog::warn("\t\t[time profile] multiply:{} s, BLAS:{} s, Precondition:{} s", a1, a2, a3);
                return true;
            }
            apply_preconditioner(er, ez);
            // a3 += timer.durationInSeconds(); timer.reset();
            double rho_new = ez.dot(er);// BLAS::dot(ez, er);
            double beta = rho_new / rho;
            ez += beta * es;
            // BLAS::add_scaled(beta, es, ez);
            es.swap(ez); // s=beta*s+z
            // a2 += timer.durationInSeconds(); timer.reset();
            rho = rho_new;
         }
         if(callback) callback(iteration);
         iterations_out = iteration;
         return false;
   }

 protected:

    // internal structures
    SparseColumnLowerFactor<T> ic_factor; // modified incomplete cholesky factor
    std::vector<T> m, z, s, r; // temporary vectors for PCG
    VectorXd em, ez, es, er; // temporary vectors for PCG
    FixedSparseMatrix<T> fixed_matrix; // used within loop

    // parameters
    T tolerance_factor;
    int max_iterations;
    T modified_incomplete_cholesky_parameter;
    T min_diagonal_ratio;

    void form_preconditioner(const SparseMatrix<T>& matrix) {
         BENCHMARK_SCOPED_TIMER_SECTION t("form precondition");
         factor_modified_incomplete_cholesky0(matrix, ic_factor, modified_incomplete_cholesky_parameter, min_diagonal_ratio);
    }

    void apply_preconditioner(const std::vector<T> &x, std::vector<T> &result) {
        solve_lower(ic_factor, x, result);
        solve_lower_transpose_in_place(ic_factor,result);
    }

    void apply_preconditioner_jacobian(const std::vector<T> &b, std::vector<T> &u) {
         const T omega = 2/3.0;
         const T oneMinusOmega = 1.0 - omega;
         // x_{k+1} = \omega D^(-1) [b-(L+U)x_k] + (1-\omega) x_{k+1}
         // multiply_non_diagonal(matrix, u, m);
         // ThreadPool::parallelForTF(0, (int)m.size(), [&](int i) {
         //    m[i] = b[i] - m[i];
         // });
         // multiply_inv_diagonal(matrix, m, m);
         // ThreadPool::parallelForTF(0, (int)m.size(), [&](int i) {
         //    u[i] = oneMinusOmega * u[i] + omega * m[i];
         // });
         for(int k = 0; k < JACOBIAN_P; ++k) { 
            ThreadPool::parallelForTF(0u, fixed_matrix.n, [&](uint i){
               T tmp = 0, diag = 0;
               for(uint j = fixed_matrix.rowstart[i]; j < fixed_matrix.rowstart[i+1]; ++j) {
                  if(fixed_matrix.colindex[j] != i) {
                     tmp += fixed_matrix.value[j] * u[fixed_matrix.colindex[j]];
                  } else {
                     diag = fixed_matrix.value[j];
                  }
               }
               tmp = (b[i] - tmp) / diag;
               u[i] = oneMinusOmega * u[i] + omega * tmp;
            });
         }
    }

    void apply_preconditioner(const Eigen::Matrix<T, -1, 1> &x, Eigen::Matrix<T, -1, 1> &result) {
        solve_lower(ic_factor, x, result);
        solve_lower_transpose_in_place(ic_factor,result);
    }
};

#endif

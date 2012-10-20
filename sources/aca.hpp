// Copyright (c) 2012, Matthias Messner
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
//  * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.

#ifndef ACA_HPP
#define ACA_HPP

#include "sources/blas.h"


/*! \defgroup mapprox Low-Rank Approximations Schemes

  \par Truncated Singular Value Decomposition
  A truncated singular value decomposition (tSVD) gives the optimal low-rank
  representation of a given dense matrix with respect to a prescribed
  approximation error. However, the cost grows cubically with the size of the
  matrix to be approximated. That is why we deploy the ACA (see the following
  section) to reduce the cost. Since it lead to only an almost optimal
  low-rank we have the possiblity to use a subsequent QR decomposition
  followed by a truncated SVD (QRtSVD) to obtain the optimal low-rank very
  efficiently.

  \par Adaptive Cross Approximation
  The adaptive cross approximation (ACA) is an iterative (semi-heuristic)
  algorithm which computes the low-rank approximation of a dense matrix
  satisfying a prescribed accuracy. There exist two variants, (1) the fully
  pivoted adaptive cross approximation (fACA) and (2) the partially pivoted
  adaptive cross approximation (pACA). For more detailed information about
  this algorithms we refer to the publication "Approximation of boundary
  element matrices" in <em>Numerische Mathematik, 86(4):565–589, 2000</em> by
  M. Bebendorf. It is well known that the singular value decomposition (SVD)
  provides the optimal low-rank approximation of a given matrix relative to a
  prescribed accuracy. The ACA gives the \e almost optimal low-rank
  approximation. Moreover, the partially pivoted ACA has a complexity of
  \f$\mathcal{O}(2Nk)\f$ (low-rank \f$k\f$) as opposed to
  \f$\mathcal{O}(N^3)\f$ of the tSVD. By a subsequent application of the tSVD
  the (don't have a mathematical proof of that, though :) optimal low-rank can
  be obtained very efficiently.
 */





#include <cstdlib>

namespace dfmm {
  namespace math {

    double abs(const double val)               { return     fabs(val); }
    float  abs(const float  val)               { return     fabs(val); }
    double abs(const std::complex<double> val) { return std::abs(val); }
    float  abs(const std::complex<float>  val) { return std::abs(val); }

    void conjugate(const unsigned int m, double *const A) {}
    void conjugate(const unsigned int m, std::complex<double> *const A)
    { for (unsigned int i=0; i<m; ++i) A[i] = std::conj(A[i]); }
  }
}


namespace dfmm {
  /*! \ingroup mapprox
    
    \brief Fully pivoted ACA

    The fully pivoted adaptive cross approximation (fACA) compresses a
    far-field interaction as \f$K\sim UV^\top\f$. The fACA requires all
    entries to be computed beforehand, then the compression follows in
    \f$\mathcal{O}(kn_x + kn_y)\f$ operations based on the required accuracy
    \f$\varepsilon\f$. The matrix K will be destroyed as a result.

    \attention This implementation can treat only matrices which have no zero
    blocks! This is sufficient, though, for our needs in this library, i.e.,
    we use it for the low-rank approximation of the M2L kernels which are per
    definition dense matrices which have low-rank.

    @tparam T value type

    @param[in] K far-field to be approximated
    @param[in] nx number of rows
    @param[in] ny number of cols
    @param[in] eps prescribed accuracy
    @param[out] U matrix containing \a k column vectors
    @param[out] V matrix containing \a k row vectors
    @param[out] k final low-rank depends on prescribed accuracy \a eps
  */
  template <typename T>
  void fACA(T *const K,
            const unsigned int nx, const unsigned int ny,
            const double eps, T* &U, T* &V, unsigned int &k)
  {
    // control vectors (true if not used, false if used)
    bool *const r = new bool[nx];
    bool *const c = new bool[ny];
    for (unsigned int i=0; i<nx; ++i) r[i] = true;
    for (unsigned int j=0; j<ny; ++j) c[j] = true;

    // compute Frobenius norm of original Matrix K
    double norm2K = 0;
    for (unsigned int j=0; j<ny; ++j) {
      const T *const colK = K + j*nx;
      norm2K += blas::dot(nx, colK);
    }

    // initialize rank k and UV'
    k = 0;
    //  const int maxk = (nx + ny) / 2;
    const int maxk = nx < ny ? nx : ny;
    U = new T [nx * maxk];
    V = new T [ny * maxk];
    blas::setzero(nx*maxk, U);
    blas::setzero(ny*maxk, V);
    double norm2R;

    ////////////////////////////////////////////////
    // start fully pivoted ACA
    do {
    
      // find max(K) and argmax(K)
      double maxK = 0.;
      int pi=0, pj=0;
      for (unsigned int j=0; j<ny; ++j)
        if (c[j]) {
          const T *const colK = K + j*nx;
          for (unsigned int i=0; i<nx; ++i)
            if (r[i] && maxK < math::abs(colK[i])) {
              maxK = math::abs(colK[i]);
              pi = i; 
              pj = j;
            }
        }

      // copy pivot cross into U and V
      T *const colU = U + k*nx;
      T *const colV = V + k*ny;
      const T pivot = K[pj*nx + pi];
      for (unsigned int i=0; i<nx; ++i)
        if (r[i]) colU[i] = K[pj*nx + i];
      for (unsigned int j=0; j<ny; ++j)
        if (c[j]) colV[j] = K[j *nx + pi] / pivot;
    
      // dont use these cols and rows anymore
      c[pj] = false;
      r[pi] = false;
    
      // subtract k-th outer product from K
      for (unsigned int j=0; j<ny; ++j)
        if (c[j]) {
          T *const colK = K + j*nx;
          blas::axpy(nx, T(-1. * colV[j]), colU, colK);
        }

      // compute Frobenius norm of updated K
      norm2R = 0.;
      for (unsigned int j=0; j<ny; ++j)
        if (c[j]) {
          const T *const colK = K + j*nx;
          norm2R += blas::dot(nx, colK);
        }

      // increment rank k
      k++;

      // throw if rank gets high (k<min(nx,ny))
      if (k >= maxk)
        throw std::runtime_error("ACA did not converge with " + k);

    } while (norm2R > eps*eps * norm2K);
    ////////////////////////////////////////////////

    // only needed in the complex valued case   
    math::conjugate(ny * k, V);

    // free memory
    delete [] r;
    delete [] c;
  }









  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////





  /*! \ingroup mapprox

    \brief Partially pivoted ACA

    The partially pivoted adaptive cross approximation (pACA) compresses a
    far-field interaction as \f$K\sim UV^\top\f$. The pACA computes the matrix
    entries on the fly, as they are needed. The compression follows in
    \f$\mathcal{O}(kn_x + kn_y)\f$ operations based on the required accuracy
    \f$\varepsilon\f$. The matrix K will be destroyed as a result.

    \attention This implementation can treat only matrices which have no zero
    blocks! This is sufficient, though, for our needs in this library, i.e.,
    we use it for the low-rank approximation of the M2L kernels which are per
    definition dense matrices which have low-rank.

    @tparam ComputerType the functor type which allows to compute matrix entries
  
    @param[in] Computer the entry-computer functor
    @param[in] eps prescribed accuracy
    @param[in] nx number of rows
    @param[in] ny number of cols
    @param[out] U matrix containing \a k column vectors
    @param[out] V matrix containing \a k row vectors
    @param[out] k final low-rank depends on prescribed accuracy \a eps
  */
  template <typename T,
            typename ComputerType>
  void pACA(const ComputerType& Computer,
            const unsigned int nx, const unsigned int ny,
            const double eps, T* &U, T* &V, unsigned int &k)
  {
    // control vectors (true if not used, false if used)
    bool *const r = new bool[nx];
    bool *const c = new bool[ny];
    for (unsigned int i=0; i<nx; ++i) r[i] = true;
    for (unsigned int j=0; j<ny; ++j) c[j] = true;
  
    // initialize rank k and UV'
    k = 0;
    const double eps2 = eps * eps;
    const int maxk = nx < ny ? nx : ny;
    //  const int maxk = (nx + ny) / 2;
    U = new T [nx * maxk];
    V = new T [ny * maxk];
  
    // initialize norm
    double norm2S(0.);
    double norm2uv(0.);
  
    ////////////////////////////////////////////////
    // start partially pivoted ACA with initial guess for I and J
    unsigned int J = 0, I = 0;
  
    //////////////////////////////////////////////////////////////////
    // if more rows then columns do row based ACA ////////////////////
    //////////////////////////////////////////////////////////////////
    if (nx <= ny) 
      do {
        T *const colU = U + nx*k;
        T *const colV = V + ny*k;
        
        ////////////////////////////////////////////
        // compute row I and its residual
        Computer(I, I+1, 0, ny, colV);
        r[I] = false;
        for (unsigned int l=0; l<k; ++l) {
          T *const u = U + nx*l;
          T *const v = V + ny*l;
          blas::axpy(ny, T(-1. * u[I]), v, colV);
        }
        
        // find max of residual and argmax
        double maxval = 0.;
        for (unsigned int j=0; j<ny; ++j) {
          const double abs_val = math::abs(colV[j]);
          if (c[j] && maxval < abs_val) {
            maxval = abs_val;
            J = j;
          }
        }
        // find pivot and scale column of V
        const T pivot = 1. / colV[J];
        blas::scal(ny, pivot, colV);
        
        ////////////////////////////////////////////
        // compute col J and its residual
        Computer(0, nx, J, J+1, colU);
        c[J] = false;
        for (unsigned int l=0; l<k; ++l) {
          T *const u = U + nx*l;
          T *const v = V + ny*l;
          blas::axpy(nx, T(-1. * v[J]), u, colU);
        }
        
        // find max of residual and argmax
        maxval = 0.;
        for (unsigned int i=0; i<nx; ++i) {
          const double abs_val = math::abs(colU[i]);
          if (r[i] && maxval < abs_val) {
            maxval = abs_val;
            I = i;
          }
        }
        
        ////////////////////////////////////////////
        // increment Frobenius norm: |Sk|^2 += |uk|^2 |vk|^2 + 2 sumj ukuj vjvk
        double normuuvv(0.);
        for (unsigned int l=0; l<k; ++l)
          normuuvv += math::abs(blas::scpr(nx, colU, U+nx*l) *
                                blas::scpr(ny, V+ny*l, colV));
        norm2uv = blas::dot(nx, colU) * blas::dot(ny, colV);
        norm2S += norm2uv + 2*normuuvv;
        
        ////////////////////////////////////////////
        // increment low-rank
        k++;
        
        // throw if rank gets high (k<min(nx,ny))
        if (k >= maxk)
          throw std::runtime_error("ACA did not converge with " + k);
        
      } while (norm2uv > eps2 * norm2S);

    //////////////////////////////////////////////////////////////////
    // if more rows then columns do col based ACA ////////////////////
    //////////////////////////////////////////////////////////////////
    else 
      do {
        T *const colU = U + nx*k;
        T *const colV = V + ny*k;
        
        ////////////////////////////////////////////
        // compute col J and its residual
        Computer(0, nx, J, J+1, colU);
        c[J] = false;
        for (unsigned int l=0; l<k; ++l) {
          T *const u = U + nx*l;
          T *const v = V + ny*l;
          blas::axpy(nx, T(-1. * v[J]), u, colU);
        }
        
        // find max of residual and argmax
        double maxval = 0.;
        for (unsigned int i=0; i<nx; ++i) {
          const double abs_val = math::abs(colU[i]);
          if (r[i] && maxval < abs_val) {
            maxval = abs_val;
            I = i;
          }
        }
        // find pivot and scale column of V
        const T pivot = 1. / colU[I];
        blas::scal(nx, pivot, colU);
        
        ////////////////////////////////////////////
        // compute col J and its residual
        Computer(I, I+1, 0, ny, colV);
        r[I] = false;
        for (unsigned int l=0; l<k; ++l) {
          T *const u = U + nx*l;
          T *const v = V + ny*l;
          blas::axpy(ny, T(-1. * u[I]), v, colV);
        }
        
        // find max of residual and argmax
        maxval = 0.;
        for (unsigned int j=0; j<ny; ++j) {
          const double abs_val = math::abs(colV[j]);
          if (c[j] && maxval < abs_val) {
            maxval = abs_val;
            J = j;
          }
        }
        
        ////////////////////////////////////////////
        // increment Frobenius norm: |Sk|^2 += |uk|^2 |vk|^2 + 2 sumj ukuj vjvk
        double normuuvv(0.);
        for (unsigned int l=0; l<k; ++l)
          normuuvv += math::abs(blas::scpr(nx, colU, U+nx*l) *
                                blas::scpr(ny, V+ny*l, colV));
        norm2uv = blas::dot(nx, colU) * blas::dot(ny, colV);
        norm2S += norm2uv + 2*normuuvv;
        
        ////////////////////////////////////////////
        // increment low-rank
        k++;
        
        // throw if rank gets high (k<min(nx,ny))
        if (k >= maxk)
          throw std::runtime_error("ACA did not converge with " + k);
        
      } while (norm2uv > eps2 * norm2S);
    
    // only needed in the complex valued case   
    math::conjugate(ny * k, V);
    
    // free memory
    delete [] r;
    delete [] c;
  }



} // end namespace dfmm

#endif

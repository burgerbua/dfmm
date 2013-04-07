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

#ifndef BLAS_H
#define BLAS_H

#include <assert.h>
#include <stdlib.h>

#include <complex>
#define dcomp std :: complex< double >
#define scomp std :: complex< float  >

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#define SIGN(a) ((a) >= 0 ? 1.0 : -1.0)



const double D_ZERO =  0.0;
const double D_ONE  =  1.0;
const double D_MONE = -1.0;
const float S_ZERO  =  0.0;
const float S_ONE   =  1.0;
const float S_MONE  = -1.0;
const scomp C_ZERO  =  0.0;
const scomp C_ONE   =  1.0;
const scomp C_MONE  = -1.0;
const dcomp Z_ZERO  =  0.0;
const dcomp Z_ONE   =  1.0;
const dcomp Z_MONE  = -1.0;

const double D_PREC = 1e-16;

const unsigned N_ONE = 1;
const int N_MONE = -1;
const char JOB_STR[] = "NTOSVULCRMD";





namespace dfmm
{
  extern "C"
  {
    /******************************************************************/
    //double precision real
    /******************************************************************/
    void dhad_(const unsigned int*, const double*, const double*,
               const unsigned int*, double*, const unsigned int*);

    unsigned idamax_(const unsigned*, const double*, const unsigned*);
    void dcopy_(const unsigned*, const double*, const unsigned*,
                double*, const unsigned*);
    void daxpy_(const unsigned*, const double*, const double*,
                const unsigned*, double*, const unsigned*);
    void dscal_(const unsigned*, const double*, const double*,
                const unsigned*);
    double ddot_(const unsigned*, const double*, const unsigned*,
                 const double*, const unsigned*);
    double dnrm2_(const unsigned*, const double*, const unsigned*);

    void dgemm_(const char*, const char*, const unsigned*, const unsigned*,
                const unsigned*, const double*, double*, const unsigned*,
                double*, const unsigned*, const double*, double*,
                const unsigned*);
    void dgemv_(const char*, const unsigned*, const unsigned*, const double*,
                const double*, const unsigned*, const double*, const unsigned*,
                const double*, double*, const unsigned*);
    void dorgqr_(const unsigned*, const unsigned*, const unsigned*,
                 double*, const unsigned*, double*, double*,
                 const unsigned*, int*);
    void dormqr_(const char*, const char*, const unsigned*, const unsigned*,
                 const unsigned*, double*, const unsigned*, double*,
                 double*, const unsigned*, double*, const unsigned*,
                 int*);
    void dgeqrf_(const unsigned*, const unsigned*, double*, const unsigned*,
                 double*, double*, const unsigned*, int*);
    void dgesvd_(const char*, const char*, const unsigned*, const unsigned*,
                 double*, const unsigned*, double*, double*, const unsigned*,
                 double*, const unsigned*, double*, const unsigned*, int*);
    void dgetrf_(const unsigned*, const unsigned*, double*, const unsigned*,
                 unsigned*, int*);
    void dgetrs_(const char*, const unsigned*, const unsigned*, double*,
                 const unsigned*, const unsigned*, double*, const unsigned*,
                 int*);
    void dgetri_(const unsigned*, double*, const unsigned*, unsigned*, double*,
                 const unsigned*, int*);
    void dsptrf_(const char*, const unsigned*, const double*, const unsigned*,
                 int*);
    void dsptri_(const char*, const unsigned*, const double*, const unsigned*,
                 double*, int*);

    /******************************************************************/
    //single precision real
    /******************************************************************/
    unsigned isamax_(const unsigned*, const float*, const unsigned*);
    void scopy_(const unsigned*, const float*, const unsigned*,
                float*, const unsigned*);
    void saxpy_(const unsigned*, const float*, const float*,
                const unsigned*, float*, const unsigned*);
    void sscal_(const unsigned*, const float*, const float*,
                const unsigned*);
    float sdot_(const unsigned*, const float*, const unsigned*,
                const float*, const unsigned*);
    float snrm2_(const unsigned*, const float*, const unsigned*);

    void sgemm_(const char*, const char*, const unsigned*, const unsigned*,
                const unsigned*, const float*, float*, const unsigned*,
                float*, const unsigned*, const float*, float*,
                const unsigned*);
    void sgemv_(const char*, const unsigned*, const unsigned*, const float*,
                const float*, const unsigned*, const float*, const unsigned*,
                const float*, float*, const unsigned*);
    void sorgqr_(const unsigned*, const unsigned*, const unsigned*,
                 float*, const unsigned*, float*, float*,
                 const unsigned*, int*);
    void sormqr_(const char*, const char*, const unsigned*, const unsigned*,
                 const unsigned*, float*, const unsigned*, float*,
                 float*, const unsigned*, float*, const unsigned*,
                 int*);
    void sgeqrf_(const unsigned*, const unsigned*, float*, const unsigned*,
                 float*, float*, const unsigned*, int*);
    void sgesvd_(const char*, const char*, const unsigned*, const unsigned*,
                 float*, const unsigned*, float*, float*, const unsigned*,
                 float*, const unsigned*, float*, const unsigned*, int*);
    void sgetrf_(const unsigned*, const unsigned*, float*, const unsigned*,
                 unsigned*, int*);
    void sgetri_(const unsigned*, float*, const unsigned*, unsigned*, float*,
                 const unsigned*, int*);
    void ssptrf_(const char*, const unsigned*, const float*, const unsigned*,
                 int*);
    void ssptri_(const char*, const unsigned*, const float*, const unsigned*,
                 float*, int*);

    /******************************************************************/
    //single precision complex
    /******************************************************************/

    void ccopy_(const unsigned*, const scomp*, const unsigned*,
                scomp*, const unsigned*);
    void caxpy_(const unsigned*, const scomp*, const scomp*,
                const unsigned*, scomp*, const unsigned*);
    void cscal_(const unsigned*, const scomp*, const scomp*,
                const unsigned*);

#ifdef DFMM_USE_MKL
    void cdotc_(scomp*, const unsigned*, const scomp*, const unsigned*,
                 const scomp*, const unsigned*);
#else
    scomp cdotc_(const unsigned*, const scomp*, const unsigned*,
                 const scomp*, const unsigned*);
#endif

    float scnrm2_(const unsigned*, const scomp*, const unsigned*);

    void cgemm_(const char*, const char*, const unsigned*, const unsigned*,
                const unsigned*, const scomp*, scomp*, const unsigned*,
                scomp*, const unsigned*, const scomp*, scomp*,
                const unsigned*);
    void cgemv_(const char*, const unsigned*, const unsigned*, const scomp*,
                const scomp*, const unsigned*, const scomp*, const unsigned*,
                const scomp*, scomp*, const unsigned*);
    void cungqr_(const unsigned*, const unsigned*, const unsigned*,
                 scomp*, const unsigned*, scomp*, scomp*,
                 const unsigned*, int*);
    void cunmqr_(const char*, const char*, const unsigned*, const unsigned*,
                 const unsigned*, scomp*, const unsigned*, scomp*,
                 scomp*, const unsigned*, scomp*, const unsigned*,
                 int*);
    void cgeqrf_(const unsigned*, const unsigned*, scomp*, const unsigned*,
                 scomp*, scomp*, const unsigned*, int*);
    void cgesvd_(const char*, const char*, const unsigned*, const unsigned*,
                 scomp*, const unsigned*, float*, scomp*, const unsigned*,
                 scomp*, const unsigned*, scomp*, const unsigned*, float*,
                 int*);
    void cgetrf_(const unsigned*, const unsigned*, scomp*, const unsigned*,
                 unsigned*, int*);
    void cgetri_(const unsigned*, scomp*, const unsigned*, unsigned*, scomp*,
                 const unsigned*, int*);
    void chpmv_(const char*, const unsigned*, const scomp*,
                const scomp*, const scomp*, const unsigned*,
                const scomp*, scomp*, const unsigned*);
    void csptrf_(const char*, const unsigned*, const scomp*, const unsigned*,
                 int*);
    void csptri_(const char*, const unsigned*, const scomp*, const unsigned*,
                 scomp*, int*);

    /******************************************************************/
    //double precision complex
    /******************************************************************/
    void zhad_(const char*, const unsigned int*, const dcomp*, const dcomp*,
               const unsigned int*, dcomp*, const unsigned int*);

    void zcopy_(const unsigned*, const dcomp*, const unsigned*,
                dcomp*, const unsigned*);
    void zaxpy_(const unsigned*, const dcomp*, const dcomp*,
                const unsigned*, dcomp*, const unsigned*);
    void zscal_(const unsigned*, const dcomp*, const dcomp*,
                const unsigned*);

#if defined(DFMM_USE_MKL) || defined(DFMM_APPLE)
    void zdotc_(dcomp*, const unsigned*, const dcomp*, const unsigned*,
                 const dcomp*, const unsigned*);
#else
    dcomp zdotc_(const unsigned*, const dcomp*, const unsigned*,
                 const dcomp*, const unsigned*);
#endif

    double dznrm2_(const unsigned*, const dcomp*, const unsigned*);

    void zgemm_(const char*, const char*, const unsigned*, const unsigned*,
                const unsigned*, const dcomp*, dcomp*, const unsigned*,
                dcomp*, const unsigned*, const dcomp*, dcomp*,
                const unsigned*);
    void zgemv_(const char*, const unsigned*, const unsigned*, const dcomp*,
                const dcomp*, const unsigned*, const dcomp*, const unsigned*,
                const dcomp*, dcomp*, const unsigned*);
    void zungqr_(const unsigned*, const unsigned*, const unsigned*,
                 dcomp*, const unsigned*, dcomp*, dcomp*,
                 const unsigned*, int*);
    void zunmqr_(const char*, const char*, const unsigned*, const unsigned*,
                 const unsigned*, dcomp*, const unsigned*, dcomp*,
                 dcomp*, const unsigned*, dcomp*, const unsigned*,
                 int*);
    void zgeqrf_(const unsigned*, const unsigned*, dcomp*, const unsigned*,
                 dcomp*, dcomp*, const unsigned*, int*);
    void zgesvd_(const char*, const char*, const unsigned*, const unsigned*,
                 dcomp*, const unsigned*, double*, dcomp*, const unsigned*,
                 dcomp*, const unsigned*, dcomp*, const unsigned*, double*,
                 int*);
    void zgetrf_(const unsigned*, const unsigned*, dcomp*, const unsigned*,
                 unsigned*, int*);
    void zgetri_(const unsigned*, dcomp*, const unsigned*, unsigned*, dcomp*,
                 const unsigned*, int*);
    void zhpmv_(const char*, const unsigned*, const dcomp*,
                const dcomp*, const dcomp*, const unsigned*,
                const dcomp*, dcomp*, const unsigned*);
    void zsptrf_(const char*, const unsigned*, const dcomp*, const unsigned*,
                 int*);
    void zsptri_(const char*, const unsigned*, const dcomp*, const unsigned*,
                 dcomp*, int*);
  }

} // end namespace dfmm


//+++++++++++++++++++++
using namespace dfmm;
//+++++++++++++++++++++

namespace blas
{

  inline unsigned maxi(const unsigned n, double* const v)
  { return idamax_(&n, v, &N_ONE); }
  inline unsigned maxi(const unsigned n, float* const v)
  { return isamax_(&n, v, &N_ONE); }

  inline void load(const unsigned n, double e, double* const v)
  { for (unsigned i=0; i<n; ++i) v[i] = e; }
  inline void load(const unsigned n, float e, float* const v)
  { for (unsigned i=0; i<n; ++i) v[i] = e; }
  inline void load(const unsigned n, scomp e, scomp* const v)
  { for (unsigned i=0; i<n; ++i) v[i] = e; }
  inline void load(const unsigned n, dcomp e, dcomp* const v)
  { for (unsigned i=0; i<n; ++i) v[i] = e; }

  inline void setzero(const unsigned n, unsigned* const v)
  {
    for (unsigned i=0; i<n; ++i) v[i] = 0;
  }
  inline void setzero(const unsigned n, double* const v)
  {
    load(n, D_ZERO, v);
  }
  inline void setzero(const unsigned n, float* const v)
  {
    load(n, S_ZERO, v);
  }
  inline void setzero(const unsigned n, scomp* const v)
  {
    load(n, C_ZERO, v);
  }
  inline void setzero(const unsigned n, dcomp* const v)
  {
    load(n, Z_ZERO, v);
  }

  inline double nrm2(const unsigned n, double* const v)
  {
    return dnrm2_(&n, v, &N_ONE);
  }
  inline float nrm2(const unsigned n, float* const v)
  {
    return snrm2_(&n, v, &N_ONE);
  }
  inline float nrm2(const unsigned n, scomp* const v)
  {
    return scnrm2_(&n, v, &N_ONE);
  }
  inline double nrm2(const unsigned n, dcomp* const v)
  {
    return dznrm2_(&n, v, &N_ONE);
  }

  inline void copy(const unsigned n, double* orig, double* dest)
  {
    dcopy_(&n, orig, &N_ONE, dest, &N_ONE);
  }
  inline void copy(const unsigned n, float* orig, float* dest)
  {
    scopy_(&n, orig, &N_ONE, dest, &N_ONE);
  }
  inline void copy(const unsigned n, scomp* orig, scomp* dest)
  {
    ccopy_(&n, orig, &N_ONE, dest, &N_ONE);
  }
  inline void copy(const unsigned n, dcomp* orig, dcomp* dest)
  {
    zcopy_(&n, orig, &N_ONE, dest, &N_ONE);
  }


  inline void copy(const unsigned n, double* orig, const unsigned inco,
                   double* dest, const unsigned incd)
  {
    dcopy_(&n, orig, &inco, dest, &incd);
  }
  inline void copy(const unsigned n, float* orig, const unsigned inco,
                   float* dest, const unsigned incd)
  {
    scopy_(&n, orig, &inco, dest, &incd);
  }
  inline void copy(const unsigned n, scomp* orig, const unsigned inco,
                   scomp* dest, const unsigned incd)
  {
    ccopy_(&n, orig, &inco, dest, &incd);
  }
  inline void copy(const unsigned n, dcomp* orig, const unsigned inco,
                   dcomp* dest, const unsigned incd)
  {
    zcopy_(&n, orig, &inco, dest, &incd);
  }

  // Conv.Copy double2float
  inline void copy(const unsigned n, double* orig, float* dest)
  {
    for (unsigned i=0; i<n; i++) dest[i] = (float) orig[i];
  }

  // Conv.Copy float2double
  inline void copy(const unsigned n, float* orig, double* dest)
  {
    for (unsigned i=0; i<n; i++) dest[i] = (double) orig[i];
  }

  // Conv.Copy dcomp2scomp
  inline void copy(const unsigned n, dcomp* orig, scomp* dest)
  {
    for (unsigned i=0; i<n; i++)
      dest[i].real()
        = (float) orig[i].real(), dest[i].imag() = (float) orig[i].imag(); 
  }

  // Conv.Copy scomp2dcomp
  inline void copy(const unsigned n, scomp* orig, dcomp* dest)
  {
    for (unsigned i=0; i<n; i++)
      dest[i].real()
        = (double) orig[i].real(), dest[i].imag() = (double) orig[i].imag();
  }

  // Scalar product conj(x)*y
  inline double scpr(const unsigned n, const double* const v1,
                     const double* const v2)
  { return ddot_(&n, v1, &N_ONE, v2, &N_ONE); }
  inline float scpr(const unsigned n, const float* const v1,
                    const float* const v2)
  { return sdot_(&n, v1, &N_ONE, v2, &N_ONE); }
  inline scomp scpr(const unsigned n, const scomp* const v1,
                    const scomp* const v2)
  {
#ifdef DFMM_USE_MKL
    scomp result;
    cdotc_(&result, &n, v1, &N_ONE, v2, &N_ONE);
    return result;
#else
    return cdotc_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
  }
  inline dcomp scpr(const unsigned n, const dcomp* const v1,
                    const dcomp* const v2)
  {
#if defined(DFMM_USE_MKL) || defined(DFMM_APPLE)
    dcomp result;
    zdotc_(&result, &n, v1, &N_ONE, v2, &N_ONE);
    return result;
#else
    return zdotc_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
  }

  // dot product
  static double dot(const unsigned n, const double *const v)
  { return scpr(n, v, v); }
  static float dot(const unsigned n, const float *const v)
  { return scpr(n, v, v); }
  static double dot(const unsigned n, const scomp *const v)
  { return std::real(scpr(n, v, v)); }
  static double dot(const unsigned n, const dcomp *const v)
  {
    dcomp value(0.);
    for (unsigned int i=0; i<n; ++i)
      value += v[i] * std::conj(v[i]);
    return std::real(value);
    //return std::real(scpr(n, v, v));
  }
  


  inline void had(const unsigned int n, const double d, const double *const x,
                  double *const y)
  {
    dhad_(&n, &d, x, &N_ONE, y, &N_ONE);
  }
  inline void hadm(const unsigned int n, const dcomp d, const dcomp *const x,
                   dcomp *const y)
  {
    zhad_(JOB_STR+9,&n, &d, x, &N_ONE, y, &N_ONE);
  }
  inline void hadd(const unsigned int n, const dcomp d, const dcomp *const x,
                   dcomp *const y)
  {
    zhad_(JOB_STR+10,&n, &d, x, &N_ONE, y, &N_ONE);
  }



  inline void add(const unsigned n, double* const x, double* const y)
  {
    daxpy_(&n, &D_ONE, x, &N_ONE, y, &N_ONE);
  }
  inline void add(const unsigned n, float* const x, float* const y)
  {
    saxpy_(&n, &S_ONE, x, &N_ONE, y, &N_ONE);
  }
  inline void add(const unsigned n, scomp* const x, scomp* const y)
  {
    caxpy_(&n, &C_ONE, x, &N_ONE, y, &N_ONE);
  }
  inline void add(const unsigned n, dcomp* const x,dcomp* const y)
  {
    zaxpy_(&n, &Z_ONE, x, &N_ONE, y, &N_ONE);
  }

  inline void axpy(const unsigned n, const double d, const double* const x,
                   double* const y)
  {
    daxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
  }
  inline void axpy(const unsigned n, const float d, const float* const x,
                   float* const y)
  {
    saxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
  }
  inline void axpy(const unsigned n, const scomp d, const scomp* const x,
                   scomp* const y)
  {
    caxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
  }
  inline void axpy(const unsigned n, const dcomp d, const dcomp* const x,
                   dcomp* const y)
  {
    zaxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
  }

  // scale (unary increment)
  inline void scal(const unsigned n, const double d, double* const x)
  { dscal_(&n, &d, x, &N_ONE); }
  inline void scal(const unsigned n, const float d, float* const x)
  { sscal_(&n, &d, x, &N_ONE); }
  inline void scal(const unsigned n, const scomp d, scomp* const x)
  { cscal_(&n, &d, x, &N_ONE); }
  inline void scal(const unsigned n, const dcomp d, dcomp* const x)
  { zscal_(&n, &d, x, &N_ONE); }

  // scale (variable increment)
  inline void scal(const unsigned n, const double d, double* const x,
                   const unsigned incd)
  { dscal_(&n, &d, x, &incd); }
  inline void scal(const unsigned n, const float d, float* const x,
                   const unsigned incd)
  { sscal_(&n, &d, x, &incd); }
  inline void scal(const unsigned n, const dcomp d, dcomp* const x,
                   const unsigned incd)
  { zscal_(&n, &d, x, &incd); }
  inline void scal(const unsigned n, const scomp d, scomp* const x,
                   const unsigned incd)
  { cscal_(&n, &d, x, &incd); }

  
  template<class T> inline void normalize(unsigned n, T* x)
    {
      T s = 1.0/blas::nrm2(n, x);
      blas::scal(n, s, x);
    }
  
  template<class T> inline void mkOrth(unsigned n, const T* v1, T* v2)
    {
      T s = -blas::scpr(n, v1, v2);
      blas::axpy(n, s, v1, v2);
    }

  template<class T>
    inline void scal(const unsigned n, const T d, T* const x, T* const y)
    {
      for  (unsigned i=0; i<n; ++i) y[i] = d * x[i];
    }


  // y = d Ax
  inline void gemv(const unsigned m, const unsigned n, double d, double* A,
                   double *x, double *y)
  {
    dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE);
  }
  inline void gemv(const unsigned m, const unsigned n, float d, float* A,
                   float *x, float *y)
  {
    sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE);
  }
  inline void gemv(const unsigned m, const unsigned n, scomp d, scomp* A,
                   scomp *x, scomp *y)
  {
    cgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE);
  }
  inline void gemv(const unsigned m, const unsigned n, dcomp d, dcomp* A,
                   dcomp *x, dcomp *y)
  {
    zgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE);
  }


  // y += d Ax
  inline void gemva(const unsigned m, const unsigned n, double d, double* A,
                    double *x, double *y)
  {
    dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);
  }
  inline void gemva(const unsigned m, const unsigned n, float d, float* A,
                    float *x, float *y)
  {
    sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);
  }
  inline void gemva(const unsigned m, const unsigned n, scomp d, scomp* A,
                    scomp *x, scomp *y)
  {
    cgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);
  }
  inline void gemva(const unsigned m, const unsigned n, dcomp d, dcomp* A,
                    dcomp *x, dcomp *y)
  {
    zgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);
  }

  // y = d A^H x
  inline void gemhv(const unsigned m, const unsigned n, double d, double* A,
                    double *x, double *y)
  {
    dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE);
  }
  inline void gemhv(const unsigned m, const unsigned n, float d, float* A,
                    float *x, float *y)
  {
    sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE);
  }
  inline void gemhv(const unsigned m, const unsigned n, scomp d, scomp* A,
                    scomp *x, scomp *y)
  {
    cgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE);
  }
  inline void gemhv(const unsigned m, const unsigned n, dcomp d, dcomp* A,
                    dcomp *x, dcomp *y)
  {
    zgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE);
  }

  // y += d A^H x
  inline void gemhva(const unsigned m, const unsigned n, double d, double* A,
                     double *x, double *y)
  {
    dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);
  }
  inline void gemhva(const unsigned m, const unsigned n, float d, float* A,
                     float *x, float *y)
  {
    sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);
  }
  inline void gemhva(const unsigned m, const unsigned n, scomp d, scomp* A,
                     scomp *x, scomp *y)
  {
    cgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);
  }
  inline void gemhva(const unsigned m, const unsigned n, dcomp d, dcomp* A,
                     dcomp *x, dcomp *y)
  {
    zgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);
  }

  // y += d A^T x
  inline void gemtva(const unsigned m, const unsigned n, double d, double* A,
                     double *x, double *y)
  {
    dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);
  }
  inline void gemtva(const unsigned m, const unsigned n, float d, float* A,
                     float *x, float *y)
  {
    sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);
  }
  inline void gemtva(const unsigned m, const unsigned n, scomp d, scomp* A,
                     scomp *x, scomp *y)
  {
    cgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);
  }
  inline void gemtva(const unsigned m, const unsigned n, dcomp d, dcomp* A,
                     dcomp *x, dcomp *y)
  {
    zgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);
  }

  // C = d A B, A is m x p, B is p x n
  inline void gemm(unsigned m, unsigned p, unsigned n, double d,
                   double* A, unsigned ldA, double* B, unsigned ldB,
                   double* C, unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &D_ZERO, C, &ldC);
  }
  inline void gemm(unsigned m, unsigned p, unsigned n, float d,
                   float* A, unsigned ldA, float* B, unsigned ldB,
                   float* C, unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &S_ZERO, C, &ldC);
  }
  inline void gemm(unsigned m, unsigned p, unsigned n, scomp d,
                   scomp* A, unsigned ldA, scomp* B, unsigned ldB,
                   scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &C_ZERO, C, &ldC);
  }
  inline void gemm(unsigned m, unsigned p, unsigned n, dcomp d,
                   dcomp* A, unsigned ldA, dcomp* B, unsigned ldB,
                   dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &Z_ZERO, C, &ldC);
  }


  // C += d A B, A is m x p, B is p x n
  inline void gemma(unsigned m, unsigned p, unsigned n, double d,
                    double* A, unsigned ldA, double* B, unsigned ldB,
                    double* C, unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &D_ONE, C, &ldC);
  }
  inline void gemma(unsigned m, unsigned p, unsigned n, float d,
                    float* A, unsigned ldA, float* B, unsigned ldB,
                    float* C, unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &S_ONE, C, &ldC);
  }
  inline void gemma(unsigned m, unsigned p, unsigned n, scomp d,
                    scomp* A, unsigned ldA, scomp* B, unsigned ldB,
                    scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &C_ONE, C, &ldC);
  }
  inline void gemma(unsigned m, unsigned p, unsigned n, dcomp d,
                    dcomp* A, unsigned ldA, dcomp* B, unsigned ldB,
                    dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &Z_ONE, C, &ldC);
  }

  // C = d A^H B, A is m x p, B is m x n
  inline void gemhm(unsigned m, unsigned p, unsigned n, double d,
                    double* A, unsigned ldA, double *B, unsigned ldB,
                    double* C, unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &D_ZERO, C, &ldC);
  }
  inline void gemhm(unsigned m, unsigned p, unsigned n, float d,
                    float* A, unsigned ldA, float *B, unsigned ldB,
                    float* C, unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &S_ZERO, C, &ldC);
  }
  inline void gemhm(unsigned m, unsigned p, unsigned n, scomp d,
                    scomp* A, unsigned ldA, scomp *B, unsigned ldB,
                    scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &C_ZERO, C, &ldC);
  }
  inline void gemhm(unsigned m, unsigned p, unsigned n, dcomp d,
                    dcomp* A, unsigned ldA, dcomp *B, unsigned ldB,
                    dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &Z_ZERO, C, &ldC);
  }

  //////////////////////////////////////////////////////////
  // C = d A^H B, A is m x p, B is m x n
  inline void gemtm(unsigned m, unsigned p, unsigned n, double d,
                    double* A, unsigned ldA, double *B, unsigned ldB,
                    double* C, unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &D_ZERO, C, &ldC);
  }
  inline void gemtm(unsigned m, unsigned p, unsigned n, float d,
                    float* A, unsigned ldA, float *B, unsigned ldB,
                    float* C, unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &S_ZERO, C, &ldC);
  }
  inline void gemtm(unsigned m, unsigned p, unsigned n, scomp d,
                    scomp* A, unsigned ldA, scomp *B, unsigned ldB,
                    scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &C_ZERO, C, &ldC);
  }
  inline void gemtm(unsigned m, unsigned p, unsigned n, dcomp d,
                    dcomp* A, unsigned ldA, dcomp *B, unsigned ldB,
                    dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &Z_ZERO, C, &ldC);
  }

  // C += d A^H B, A is m x p, B is m x n
  inline void gemhma(unsigned m, unsigned p, unsigned n, double d,
                     double* A, unsigned ldA, double *B, unsigned ldB,
                     double* C, unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &D_ONE, C, &ldC);
  }
  inline void gemhma(unsigned m, unsigned p, unsigned n, float d,
                     float* A, unsigned ldA, float *B, unsigned ldB,
                     float* C, unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &S_ONE, C, &ldC);
  }
  inline void gemhma(unsigned m, unsigned p, unsigned n, scomp d,
                     scomp* A, unsigned ldA, scomp *B, unsigned ldB,
                     scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &C_ONE, C, &ldC);
  }
  inline void gemhma(unsigned m, unsigned p, unsigned n, dcomp d,
                     dcomp* A, unsigned ldA, dcomp *B, unsigned ldB,
                     dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &Z_ONE, C, &ldC);
  }

  // C += d A^T B, A is m x p, B is m x n
  inline void gemtma(unsigned m, unsigned p, unsigned n, double d,
                     double* A, unsigned ldA, double *B, unsigned ldB,
                     double* C, unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &D_ONE, C, &ldC);
  }
  inline void gemtma(unsigned m, unsigned p, unsigned n, float d,
                     float* A, unsigned ldA, float *B, unsigned ldB,
                     float* C, unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &S_ONE, C, &ldC);
  }
  inline void gemtma(unsigned m, unsigned p, unsigned n, scomp d,
                     scomp* A, unsigned ldA, scomp *B, unsigned ldB,
                     scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &C_ONE, C, &ldC);
  }
  inline void gemtma(unsigned m, unsigned p, unsigned n, dcomp d,
                     dcomp* A, unsigned ldA, dcomp *B, unsigned ldB,
                     dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
           &Z_ONE, C, &ldC);
  }

  // C = d A B^H, A is m x p, B is n x p
  inline void gemmh(unsigned m, unsigned p, unsigned n, double d,
                    double* A, unsigned ldA, double *B, unsigned ldB,
                    double* C, unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &D_ZERO, C, &ldC);
  }
  inline void gemmh(unsigned m, unsigned p, unsigned n, float d,
                    float* A, unsigned ldA, float *B, unsigned ldB,
                    float* C, unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &S_ZERO, C, &ldC);
  }
  inline void gemmh(unsigned m, unsigned p, unsigned n, scomp d,
                    scomp* A, unsigned ldA, scomp *B, unsigned ldB,
                    scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &C_ZERO, C, &ldC);
  }
  inline void gemmh(unsigned m, unsigned p, unsigned n, dcomp d,
                    dcomp* A, unsigned ldA, dcomp *B, unsigned ldB,
                    dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &Z_ZERO, C, &ldC);
  }

  // C += d A B^H, A is m x p, B is n x p
  inline void gemmha(unsigned m, unsigned p, unsigned n, double d,
                     double* A, unsigned ldA, double *B, unsigned ldB,
                     double* C, unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
           &D_ONE, C, &ldC);
  }
  inline void gemmha(unsigned m, unsigned p, unsigned n,
                     double* A, unsigned ldA, double *B, unsigned ldB,
                     double* C, unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &D_ONE, A, &ldA, B, &ldB,
           &D_ONE, C, &ldC);
  }
  inline void gemmha(unsigned m, unsigned p, unsigned n,
                     float* A, unsigned ldA, float *B, unsigned ldB,
                     float* C, unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &S_ONE, A, &ldA, B, &ldB,
           &S_ONE, C, &ldC);
  }
  inline void gemmha(unsigned m, unsigned p, unsigned n,
                     scomp* A, unsigned ldA, scomp *B, unsigned ldB,
                     scomp* C, unsigned ldC)
  {
    cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &C_ONE, A, &ldA, B, &ldB,
           &C_ONE, C, &ldC);
  }
  inline void gemmha(unsigned m, unsigned p, unsigned n,
                     dcomp* A, unsigned ldA, dcomp *B, unsigned ldB,
                     dcomp* C, unsigned ldC)
  {
    zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &Z_ONE, A, &ldA, B, &ldB,
           &Z_ONE, C, &ldC);
  }


  // Singular Value Decomposition
  inline int gesvdS(unsigned m, unsigned n, double* A, double* S,
                    double* U, unsigned ldU, double* VT, unsigned ldVT,
                    unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+3, JOB_STR+3, &m, &n, A, &m, S, U, &ldU, VT, &ldVT,
            wk, &nwk, &INF);
    return INF;
  }
  inline int gesvd(unsigned m, unsigned n, double* A, double* S,
                   double* VT, unsigned ldVT, unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
            wk, &nwk, &INF);
    return INF;
  }
  inline int gesvd(unsigned m, unsigned n, float* A, float* S,
                   float* VT, unsigned ldVT, unsigned nwk, float* wk)
  {
    int INF;
    sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
            wk, &nwk, &INF);
    return INF;
  }
  inline int gesvd(unsigned m, unsigned n, scomp* A, float* S,
                   scomp* VT, unsigned ldVT, unsigned nwk, scomp* wk)
  {
    int INF;
    float* rwk = new float[5*MIN(m,n)];
    cgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
            wk, &nwk, rwk, &INF);
    delete [] rwk;
    return INF;
  }
  inline int gesvd(unsigned m, unsigned n, dcomp* A, double* S,
                   dcomp* VT, unsigned ldVT, unsigned nwk, dcomp* wk)
  {
    int INF;
    double* rwk = new double[5*MIN(m,n)];
    zgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
            wk, &nwk, rwk, &INF);
    delete [] rwk;
    return INF;
  }

  // Singular Value Decomposition (SO)
  inline int gesvdSO(unsigned m, unsigned n, double* A, double* S,
                     double* U, unsigned ldU, unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU,
            wk, &nwk, &INF);
    return INF;
  }
  inline int gesvdSO(unsigned m, unsigned n, float* A, float* S,
                     float* U, unsigned ldU, unsigned nwk, float* wk)
  {
    int INF;
    sgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU,
            wk, &nwk, &INF);
    return INF;
  }
  inline int gesvdSO(unsigned m, unsigned n, scomp* A, float* S,
                     scomp* U, unsigned ldU, unsigned nwk, scomp* wk)
  {
    int INF;
    float* rwk = new float[5*MIN(m,n)];
    cgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU,
            wk, &nwk, rwk, &INF);
    delete [] rwk;
    return INF;
  }
  inline int gesvdSO(unsigned m, unsigned n, dcomp* A, double* S,
                     dcomp* U, unsigned ldU, unsigned nwk, dcomp* wk)
  {
    int INF;
    double* rwk = new double[5*MIN(m,n)];
    zgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU,
            wk, &nwk, rwk, &INF);
    delete [] rwk;
    return INF;
  }

  // compute singular values
  inline int svals(unsigned m, unsigned n, double* A, double* S,
                   unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n, wk, &nwk, &INF);
    return INF;
  }
  inline int svals(unsigned m, unsigned n, float* A, float* S,
                   unsigned nwk, float* wk)
  {
    int INF;
    sgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n, wk, &nwk, &INF);
    return INF;
  }
  inline int svals(unsigned m, unsigned n, scomp* A, float* S,
                   unsigned nwk, scomp* wk)
  {
    int INF;
    float *rwk = new float[5*MIN(m,n)];
    cgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n,
            wk, &nwk, rwk, &INF);
    delete [] rwk;
    return INF;
  }
  inline int svals(unsigned m, unsigned n, dcomp* A, double* S,
                   unsigned nwk, dcomp* wk)
  {
    int INF;
    double *rwk = new double[5*MIN(m,n)];
    zgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n,
            wk, &nwk, rwk, &INF);
    delete [] rwk;
    return INF;
  }

  // triangular factorisation
  inline int getrf(const unsigned n, double* A, unsigned* ipiv)
  {
    int INF;
    dgetrf_(&n, &n, A, &n, ipiv, &INF);
    return INF;
  }
  inline int getrf(const unsigned n, float* A, unsigned* ipiv)
  {
    int INF;
    sgetrf_(&n, &n, A, &n, ipiv, &INF);
    return INF;
  }
  inline int getrf(const unsigned n, scomp* A, unsigned* ipiv)
  {
    int INF;
    cgetrf_(&n, &n, A, &n, ipiv, &INF);
    return INF;
  }
  inline int getrf(const unsigned n, dcomp* A, unsigned* ipiv)
  {
    int INF;
    zgetrf_(&n, &n, A, &n, ipiv, &INF);
    return INF;
  }

  // QR factorisation
  inline int geqrf(const unsigned m, const unsigned n, double* A,
                   double* tau, unsigned nwk, double* wk)
  {
    int INF;
    dgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int geqrf(const unsigned m, const unsigned n, float* A,
                   float* tau, unsigned nwk, float* wk)
  {
    int INF;
    sgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }

  inline int geqrf(const unsigned m, const unsigned n, scomp* A,
                   scomp* tau, unsigned nwk, scomp* wk)
  {
    int INF;
    cgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }

  inline int geqrf(const unsigned m, const unsigned n, dcomp* A,
                   dcomp* tau, unsigned nwk, dcomp* wk)
  {
    int INF;
    zgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }

  // Multiply a general Matrix with the Q-Matrix (QR factorization)
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                   double* A, double* tau, double* C,
                   unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                   float* A, float* tau, float* C,
                   unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                   scomp* A, scomp* tau, scomp* C,
                   unsigned nwk, scomp* wk)
  {
    int INF;
    cunmqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                   dcomp* A, dcomp* tau, dcomp* C,
                   unsigned nwk, dcomp* wk)
  {
    int INF;
    zunmqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                    double* A, double* tau, double* C,
                    unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6,JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                    float* A, float* tau, float* C,
                    unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6,JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                    scomp* A, scomp* tau, scomp* C,
                    unsigned nwk, scomp* wk)
  {
    int INF;
    cunmqr_(JOB_STR+6,JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                    dcomp* A, dcomp* tau, dcomp* C,
                    unsigned nwk, dcomp* wk)
  {
    int INF;
    zunmqr_(JOB_STR+6,JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }

  // return Q-Matrix (QR factorization) in A
  inline int orgqr(const unsigned m, const unsigned n, double* A, double* tau,
                   unsigned nwk, double* wk)
  {
    int INF;
    dorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int orgqr(const unsigned m, const unsigned n, float* A, float* tau,
                   unsigned nwk, float* wk)
  {
    int INF;
    sorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int orgqr(const unsigned m, const unsigned n, scomp* A, scomp* tau,
                   unsigned nwk, scomp* wk)
  {
    int INF;
    cungqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int orgqr(const unsigned m, const unsigned n, dcomp* A, dcomp* tau,
                   unsigned nwk, dcomp* wk)
  {
    int INF;
    zungqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }

  // B=A^H, A in R^(m,n)
  inline void transpose(unsigned m, unsigned n, double* A, double* B)
  {
    for (unsigned i=0; i<m; ++i) dcopy_(&n, A+i, &m, B+i*n, &N_ONE);
  }
  inline void transpose(unsigned m, unsigned n, float* A, float* B)
  {
    for (unsigned i=0; i<m; ++i) scopy_(&n, A+i, &m, B+i*n, &N_ONE);
  }
  inline void transpose(unsigned m, unsigned n, scomp* A, scomp* B)
  {
    for (unsigned i=0; i<m; i++)
      for (unsigned j=0; j<n; j++) B[j+i*n] = conj(A[i+j*m]);
  }
  inline void transpose(unsigned m, unsigned n, dcomp* A, dcomp* B)
  {
    for (unsigned i=0; i<m; i++)
      for (unsigned j=0; j<n; j++) B[j+i*n] = conj(A[i+j*m]);
  }

  // B=A^H, A in R^(ldA,n) B in R^(n,m) transpose only first m rows of A 
  inline void transpose(unsigned m, unsigned n,
                        double* A, unsigned ldA, double* B)
  {
    for (unsigned i=0; i<m; ++i) dcopy_(&n, A+i, &ldA, B+i*n, &N_ONE);
  }
  inline void transpose(unsigned m, unsigned n,
                        float* A, unsigned ldA, float* B)
  {
    for (unsigned i=0; i<m; ++i) scopy_(&n, A+i, &ldA, B+i*n, &N_ONE);
  }
  inline void transpose(unsigned m, unsigned n,
                        scomp* A, unsigned ldA, scomp* B)
  {
    for (unsigned i=0; i<m; i++)
      for (unsigned j=0; j<n; j++) B[j+i*n] = conj(A[i+j*ldA]);
  }
  inline void transpose(unsigned m, unsigned n,
                        dcomp* A, unsigned ldA, dcomp* B)
  {
    for (unsigned i=0; i<m; i++)
      for (unsigned j=0; j<n; j++) B[j+i*n] = conj(A[i+j*ldA]);
  }

}


namespace lapack
{
  
  // general inversion
  inline void geinv(const unsigned n, double* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    double* const work = new double[lwork];
    assert(work!=NULL);

    int INFO = blas::getrf(n, A, ipiv);
    assert(INFO==0);

    dgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv(const unsigned n, float* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    float* const work = new float[lwork];
    assert(work!=NULL);

    int INFO = blas::getrf(n, A, ipiv);
    assert(INFO==0);

    sgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv(const unsigned n, scomp* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    scomp* const work = new scomp[lwork];
    assert(work!=NULL);

    int INFO = blas::getrf(n, A, ipiv);
    assert(INFO==0);

    cgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv(const unsigned n, dcomp* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    dcomp* const work = new dcomp[lwork];
    assert(work!=NULL);

    int INFO = blas::getrf(n, A, ipiv);
    assert(INFO==0);

    zgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }

  inline void geinv_sym(const unsigned n, double* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    double* const work = new double[lwork];
    assert(work!=NULL);

    int INFO;
    dsptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
    assert(INFO==0);

    dsptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv_sym(const unsigned n, float* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    float* const work = new float[lwork];
    assert(work!=NULL);

    int INFO;
    ssptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
    assert(INFO==0);

    ssptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv_sym(const unsigned n, scomp* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    scomp* const work = new scomp[lwork];
    assert(work!=NULL);

    int INFO;
    csptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
    assert(INFO==0);

    csptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv_sym(const unsigned n, dcomp* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    dcomp* const work = new dcomp[lwork];
    assert(work!=NULL);

    int INFO;
    zsptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
    assert(INFO==0);

    zsptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }

}

#endif

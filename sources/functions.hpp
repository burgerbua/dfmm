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

#ifndef functions_hpp
#define functions_hpp

#include <stdexcept>

#include "sources/blas.h"


// relative point-wise L2 error
template <typename T>
const double computeL2norm(unsigned int N, T *const u, T *const v)
{
  T *const w = new T [N];
  blas::copy(N, v, w);
  blas::axpy(N, -1., u, w);
  double      dot = blas::dot(N, u);
  double diff_dot = blas::dot(N, w);
  delete [] w;
  return sqrt(diff_dot / dot);
}


// relative point-wise inf error
template <typename T>
const double computeINFnorm(unsigned int N, T *const u, T *const v)
{
  double      max = 0.;
  double diff_max = 0.;
  for (unsigned int n=0; n<N; ++n) {
    if (     max<std::abs(u[n]))           max = std::abs(u[n]);
    if (diff_max<std::abs(u[n]-v[n])) diff_max = std::abs(u[n]-v[n]);
  }
  return diff_max / max;
}


// check if \a data has nan entries
template <typename T>
bool check_nan(const unsigned int n, const T *const data)
{
  for (unsigned int i=0; i<n; ++i)
    if (data[i]!=data[i]) {
      std::cout << i << " of " << n << " gives " << data[i] << std::endl;
      return false;
    }
  return true;
}


// determine rank singular-values cutoff based on euclidian error with eps
unsigned int getRank(const double *const singular_values,
                     const unsigned int size,
                     const double eps)
{
  const double nrm2 = blas::scpr(size, singular_values, singular_values);
  double nrm2k(0.);
  for (unsigned int k=size; k>0; --k) {
    nrm2k += singular_values[k-1] * singular_values[k-1];
    if (nrm2k > eps*eps * nrm2) return k;
  }
  throw std::runtime_error("rank cannot be larger than size");
  return 0;
}





#ifdef DFMM_COUNT_FLOPS
class FlopCounter
{
  std::string description;
  unsigned int nlevel;
  unsigned long long flops;
public:
  FlopCounter(const FlopCounter& other)
  {
    description = other.description;
    nlevel      = other.nlevel;
    flops = 0;
  }
  

  FlopCounter& operator= (const FlopCounter& other)
  {
    description = other.description;
    nlevel      = other.nlevel;
    flops = 0;
    return *this;
  }

  explicit FlopCounter(const std::string _description,
                       const unsigned int _nlevel)
    : description(_description), nlevel(_nlevel), flops(0)
  {}
  
  ~FlopCounter()
  {
    if (flops)
      std:: cout << " (" << description << "[" << nlevel
                 << "] flops = " << flops << ") " << std::endl;
  }
  
  void addFlops(const unsigned int _flops)
  { flops += _flops; }

};
#else
class FlopCounter
{
public:
  explicit FlopCounter(const std::string&,
                       const unsigned int) {}
  void addFlops(const unsigned int) {}
};
#endif




#endif

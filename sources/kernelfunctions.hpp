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

#ifndef KERNELFUNCTIONS_HPP
#define KERNELFUNCTIONS_HPP



#include <complex>
#include <math.h>

#include <boost/math/constants/constants.hpp>

const double PI = boost::math::constants::pi<double>();
const double PI4 = 4 * PI;


// extendable
enum KERNEL_FUNCTION_IDENTIFIER {LAPLACE3D,
                                 HELMHOLTZ3D,
                                 LEONARD_JONES_POTENTIAL};

// probably not extedable :)
enum KERNEL_FUNCTION_TYPE {HOMOGENEOUS, NON_HOMOGENEOUS};



template <KERNEL_FUNCTION_IDENTIFIER>
class KernelFunction;


/// Provides \f$K(x,y)\f$ for the Laplace kernel
template <>
struct KernelFunction<LAPLACE3D>
{
  typedef double value_type;
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;

  template <typename point_type>
  value_type operator()(const point_type& x, const point_type& y) const
  {
    return 1. / (PI4 * x.distance(y));
  }

  value_type getScaleFactor(const value_type root_cluster_extension,
                            const unsigned int nlevel) const
  {
    const value_type extension = root_cluster_extension / pow(2, nlevel);
    return getScaleFactor(extension);
  }

  value_type getScaleFactor(const value_type extension) const
  {
    return 2. / extension;
  }

};


/// Provides \f$K(x,y)\f$ for the Helmholtz kernel
template<>
struct KernelFunction<HELMHOLTZ3D>
{
  typedef std::complex<double> value_type;
  static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;

  explicit KernelFunction(const value_type _s) : s(_s) {}

  template <typename point_type>
  value_type operator()(const point_type& x, const point_type& y) const
  {
    const double r = x.distance(y);
    return exp(-s*r) / (PI4*r);
  }
  
  double getScaleFactor(const double, const unsigned int) const
  {
    // return 1 because non homogeneous kernel functions cannot be scaled!!!
    return 1.;
  }

  double getScaleFactor(const double) const
  {
    // return 1 because non homogeneous kernel functions cannot be scaled!!!
    return 1.;
  }

private:
  const value_type s; //!< laplace parameter
};









/*!  \class EntryComputer 

  This class is a functor which provides the interface to assemble a matrix
  based on the number of rows and cols and on the coordinates x and y and the
  type of the generating matrix-kernel function.

  \brief Evaluates a given \f$K(x,y)\f$

 */
template <typename kernelfunction_type,
          typename point_type>
class EntryComputer
{
  typedef typename kernelfunction_type::value_type value_type; 
  const kernelfunction_type& kernel;

  const unsigned int nx, ny;
  const point_type *const px, *const py;

  const double *const weights;

public:
  explicit EntryComputer(const unsigned int _nx, const point_type *const _px,
                         const unsigned int _ny, const point_type *const _py,
                         const kernelfunction_type& _kernel,
                         const double *const _weights = NULL)
    : kernel(_kernel),  nx(_nx), ny(_ny), px(_px), py(_py), weights(_weights) {}


  void operator()(const unsigned int xbeg, const unsigned int xend,
                  const unsigned int ybeg, const unsigned int yend,
                  value_type *const data) const
  {
    unsigned int idx = 0;
    if (weights) {
      for (unsigned int j=ybeg; j<yend; ++j)
        for (unsigned int i=xbeg; i<xend; ++i)
          data[idx++] = weights[i] * weights[j] * kernel(px[i], py[j]);
    } else {
      for (unsigned int j=ybeg; j<yend; ++j)
        for (unsigned int i=xbeg; i<xend; ++i)
          data[idx++] = kernel(px[i], py[j]);
    }

  }
};















#endif // FCHEBMATRIXKERNEL_HPP

// [--END--]

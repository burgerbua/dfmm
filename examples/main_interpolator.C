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

// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <utility>
#include <algorithm>
#include <omp.h> 


// qbdfmm
#include "sources/dimordertraits.hpp"

#include "sources/blas.h"
#include "sources/functions.hpp"

#include "sources/mapping.hpp"
#include "sources/chebyshev.hpp"
#include "sources/tensor.hpp"

#include "sources/interpolator.hpp"
#include "sources/aca.hpp"

#include "sources/kernelfunctions.hpp"




// give random double number r with max(abs(r)) = radius
double random(const double radius)
{
  const double val = ((double)rand()/(double)RAND_MAX - (double).5) * radius;
  return val;
}


// y = UV'x where U is of size M times k and V is of size N times k
template <typename T>
void applyUV(const unsigned int M, const unsigned int N, const unsigned int k,
             const T *const U, const T *const V,
             const T *const x, T *const y)
{
  T *const temp = new T [k];
  blas::gemhv(N, k, 1., const_cast<T*const>(V), const_cast<T*const>(x), temp);
  blas::gemv( M, k, 1., const_cast<T*const>(U), temp, y);
  delete [] temp;
}



int main( int argc, char * argv[ ] )
{
  // start timer
  const double t_0 = omp_get_wtime( );

  const unsigned int DIM   = 3;  
  const unsigned int ORDER = 7;

  std::cout << "\nOrder = " << ORDER << std::endl;

  const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

  typedef DimTraits<DIM>::point_type point_type;
  typedef Chebyshev< ORDER>          basis_type;
  typedef Tensor<DIM,ORDER>         tensor_type; 


  const double diam = 2.;
  const point_type cx(0.,      0., 0.);
  const point_type cy(2.*diam, 0., 0.);
  const unsigned int M = 10000;
  const unsigned int N = 10000;

  // local to global mapper
  map_loc_glob<DIM> map_l2g_x(cx, diam);
  map_loc_glob<DIM> map_l2g_y(cy, diam);

  // interpolation points and particles
  point_type X[nnodes], Y[nnodes];
  point_type *const x = new point_type [M];
  point_type *const y = new point_type [N];
  for (unsigned int d=0; d<DIM; ++d) {
    for (unsigned int n=0; n<nnodes; ++n) {
      X[n][d] = map_l2g_x(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
      Y[n][d] = map_l2g_y(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
    }
    for (unsigned int n=0; n<M; ++n) x[n][d] = map_l2g_x(d, random(2.));
    for (unsigned int n=0; n<N; ++n) y[n][d] = map_l2g_y(d, random(2.));
  }
  
  // target and random source vector
  double *const p = new double [M];
  blas::setzero(M, p);
  double *const w = new double [N];
  for (unsigned int n=0; n<N; ++n) w[n] = (double)rand()/(double)RAND_MAX;

  // define kernel function
  typedef KernelFunction<LAPLACE3D> kernel_type;
  kernel_type kernel;
  EntryComputer<kernel_type,point_type> computer(M, x, N, y, kernel);

  // direct computation
  double *const K = new double [M * N];
  computer(0, M, 0, N, K);
  blas::gemv(M, N, 1., K, w, p);
  
  // computer P2M and L2P operators
  IOcomputerLF<DIM,ORDER,double> iocx(cx, diam);
  IOcomputerLF<DIM,ORDER,double> iocy(cy, diam);
  double *const Sx = new double [M * nnodes];
  double *const Sy = new double [N * nnodes];
  const double t0 = omp_get_wtime();
  iocx(M, x, Sx);
  const double t1 = omp_get_wtime();
  iocy(N, y, Sy);
  const double t2 = omp_get_wtime();
  std::cout << "\nL2P computed in "    << t1-t0
            << "s\nP2M computed in " << t2-t1 << "s" << std::endl;
  
  // compute M2L operator
  EntryComputer<kernel_type,point_type> m2lc(nnodes, X, nnodes, Y, kernel);
  double *const R = new double [nnodes * nnodes];
  m2lc(0, nnodes, 0, nnodes, R);

  //
  double W[nnodes];
  blas::gemhv(N, nnodes, 1., Sy, w, W);
  double P[nnodes];
  blas::gemv(nnodes, nnodes, 1., R, W, P);
  double *const pp = new double [M];
  blas::gemv( M, nnodes, 1., Sx, P, pp);
  

  std::cout << "\nL2  error = "   << computeL2norm( M, p, pp)
            << "\nInf error = " << computeINFnorm(M, p, pp)
            << std::endl;


  delete [] Sx;
  delete [] Sy;
  delete [] R;
  delete [] pp;

  delete [] K;
  
  delete [] w;
  delete [] p;
  
  delete [] x;
  delete [] y;



  ////////////////////////////////////////////////////
  std::cout << "\n- Overall time " << omp_get_wtime( ) - t_0
            << " s" << std::endl;
  ////////////////////////////////////////////////////


  return 0;
}   

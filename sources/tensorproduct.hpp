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

#ifndef tensorproduct_hpp
#define tensorproduct_hpp


/*!
  \ingroup mtranslation

  \class TensorProductHandler 

  \brief Handles the tensor-product M2M and L2L operator application

  \note So far only the TensorProductHandler is only implemented for the three
  dimensional case

  @tparam DIM spatial dimension
  @tparam ORDER interpolation order

  \par Tensor-product

  We assume the uniform case in three dimensions, i.e., each cell has 8 child
  cells. We address them based on Morton ordering the local position of the
  child cells is given by their local position 000 until 111, if we translate
  this binary representation we obtain indices from 0 to 7.

  If we consider the one-dimensional case each cell has two child cells 0 and
  1 and \f$\ell\f$ interpolation points. We set up an interpolation matrix
  \f$\mathsf{S}_0, \mathsf{S}_1 \in \mathsf{R}^{\ell\times\ell}\f$ for both of
  them. The entries are computed as \f[(\mathsf{S}_i)_{nm} = S_\ell(x_n,\bar
  x_m) \quad \mathrm{with} \quad x_n = \frac{\bar x_n }{2} +
  \frac{(-1)^{i+1}}{2} \quad \mathrm{for} \quad i=1,2 \quad \mathrm{and} \quad
  n,m=1,\dots,\ell,\f] with the %Chebyshev roots, defined as \f$\bar x_m =
  \cos\left(\frac{(2m-1)\pi}{2\ell}\right)\f$. Finally, the M2M and L2L
  kernels consist of applying both matrices \f$ \mathsf{S}_0 \f$ and \f$
  \mathsf{S}_1 \f$ (for child 0 and child 1), respectively.

  Next, we consider the three dimensional case: each cell has 8 child cells
  and the number of interpolation points in each cell is \f$\ell^3\f$. The
  most straight forward (but also the most naive) way of implementing the M2M
  and L2L kernels would be to compute 8 interpolation matrices of size
  \f$\ell^3\times\ell^3\f$.  However, by exploiting the tensor-product
  approach we can construct a more efficient scheme. We consider the child
  cell of index 1 (in binary of Morton notation 001). The full interpolation
  matrix \f$\mathsf{S}_{001}\in\mathsf{R}^{\ell^3\times\ell^3}\f$ can be
  defined as a tensor product (also Kronecker product) and reads as
  \f[\mathsf{S}_{001} = \mathsf{S}_0 \times \mathsf{S}_0 \times
  \mathsf{S}_1\f] Again, we do never compute \f$\mathsf{S}_{001}\f$
  explicitly, but we exploit properties of the tensor product. In fact, it
  allows us to reduce the cost \f$\mathcal{O}(\ell^6)\f$ of applying the
  matrix \f$\mathsf{S}_{001}\f$ to only \f$\mathcal{O}(3\ell^4)\f$. We use the
  above equation to write the M2M kernel as \f$\mathsf{w}^\prime = \left(
  \mathsf{S}_0 \times \mathsf{S}_0 \times \mathsf{S}_1 \right)^\top
  \mathsf{w}\f$ and the L2L kernel as \f$\mathsf{f} = \left( \mathsf{S}_0
  \times \mathsf{S}_0 \times \mathsf{S}_1 \right)
  \mathsf{f}^\prime\f$. Recall, we have constructed also the interpolation
  points based on the tensor product approach as \f$\bar x_m = \bar x_i \times
  \bar x_j \times \bar x_k\f$ with \f$m = k\ell^2 + j\ell + i\f$. In the
  following, we use the indices \f$i, j, k\f$ and \f$i^\prime, j^\prime,
  k^\prime\f$ to address them (indices in \f$x,y,z\f$ direction and
  \f$(^\prime)\f$ denotes the upper level). According to these indices we can
  reassemble the vectors \f$\mathsf{w,w}^\prime\f$ and
  \f$\mathsf{f,f}^\prime\f$ as tensors of size \f$\ell \times \ell \times
  \ell\f$ with the \f$ijk\f$-th entry being given by the \f$m\f$-th entry of
  the vectors, respectively. Finally, by using index notation we can write the
  M2M and L2L kernels as \f[(\mathsf{w}^\prime)_{i^\prime j^\prime k^\prime} =
  (\mathsf{S}_0)_{kk^\prime} (\mathsf{S}_0)_{jj^\prime}
  (\mathsf{S}_1)_{ii^\prime} (\mathsf{w})_{ijk} \quad \mathrm{and} \quad
  (\mathsf{f})_{ijk} = (\mathsf{S}_0)_{kk^\prime} (\mathsf{S}_0)_{jj^\prime}
  (\mathsf{S}_1)_{ii^\prime} (\mathsf{f}^\prime)_{i^\prime j^\prime
  k^\prime},\f] which can be implemented as three consequent matrix-matrix
  products followed by a permutation of the resulting matrices. Each product
  has a cost of \f$\mathcal{O}(\ell^4)\f$. This operation has to be repeated
  for each child cell.

  There are symmetries in the arrangement of the child clusters with respect
  to their parent, the same holds true for the interpolation points. Hence, it
  is also possible to express the matrix \f$ \mathsf{S}_1 \f$ as a permutation
  of \f$ \mathsf{S}_0 \f$. If we follow this idea we could group the 8 times 3
  consecutive matrix-matrix products to 8 larger matrix-matrix products. In
  reality, though, the time required for the M2M and L2L is vanishingly small
  compared to all other kernels. That is why we did not further investigate in
  this optimization.

 */
template <int DIM,
          int ORDER>
class TensorProductHandler;


template <int ORDER>
class TensorProductHandler<3,ORDER>
{
  static const unsigned int DIM = 3;
  static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

  unsigned int perm[DIM][nnodes];

public:
  TensorProductHandler() {
    for (unsigned int i=0; i<ORDER; ++i)
      for (unsigned int j=0; j<ORDER; ++j) 
        for (unsigned int k=0; k<ORDER; ++k) {
          const unsigned int index = k*ORDER*ORDER + j*ORDER + i;
          perm[0][index] = k*ORDER*ORDER + j*ORDER + i;
          perm[1][index] = i*ORDER*ORDER + k*ORDER + j;
          perm[2][index] = j*ORDER*ORDER + i*ORDER + k;
        }
  }


  // 3 * N * (ORDER*ORDER*(2*ORDER-1)) + nnodes
  template <typename T>
  void m2m(const T *const S, const unsigned int N,
           const T *const y, T *const Y) const
  {
    T exp[N*nnodes], pexp[N*nnodes];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemhm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
                const_cast<T*>(S), ORDER,
                const_cast<T*>(y), ORDER, pexp, ORDER);
    for (unsigned int e=0; e<N; ++e)
      for (unsigned int n=0; n<nnodes; ++n)
        exp[e*nnodes + n] = pexp[e*nnodes + perm[1][n]];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemhm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
                const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
                exp, ORDER, pexp, ORDER);
    for (unsigned int e=0; e<N; ++e)
      for (unsigned int n=0; n<nnodes; ++n)
        exp[e*nnodes + perm[1][n]] = pexp[e*nnodes + perm[2][n]];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemhm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
                const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
                exp, ORDER, pexp, ORDER);
    // N * nnodes
    for (unsigned int e=0; e<N; ++e)
      for (unsigned int n=0; n<nnodes; ++n)
        Y[e*nnodes + perm[2][n]] += pexp[e*nnodes + n];
  }


  // 3 * N * (ORDER*ORDER*(2*ORDER-1)) + nnodes
  template <typename T>
  void l2l(const T *const S, const unsigned int N,
           const T *const X, T *const x) const
  {
    T exp[N*nnodes], pexp[N*nnodes];
    // ORDER * (2*ORDER-1) * N*ORDER
    blas::gemm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
               const_cast<T*>(S), ORDER,
               const_cast<T*>(X), ORDER, pexp, ORDER);
    for (unsigned int e=0; e<N; ++e)
      for (unsigned int n=0; n<nnodes; ++n)
        exp[e*nnodes + n] = pexp[e*nnodes + perm[1][n]];
    // ORDER * (2*ORDER-1) * N*ORDER
    blas::gemm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
               const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
               exp, ORDER, pexp, ORDER);
    for (unsigned int e=0; e<N; ++e)
      for (unsigned int n=0; n<nnodes; ++n)
        exp[e*nnodes + perm[1][n]] = pexp[e*nnodes + perm[2][n]];
    // ORDER * (2*ORDER-1) * N*ORDER
    blas::gemm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
               const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
               exp, ORDER, pexp, ORDER);
    // N * ORDER*ORDER*ORDER
    for (unsigned int e=0; e<N; ++e)
      for (unsigned int n=0; n<nnodes; ++n)
        x[e*nnodes + perm[2][n]] += pexp[e*nnodes + n];
  }

  
  // 3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes
  template <typename T>
  void m2m(const T *const S,  const T *const y, T *const Y) const
  {
    T exp[nnodes], pexp[nnodes];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemhm(ORDER, ORDER, ORDER*ORDER, T(1.),
                const_cast<T*>(S), ORDER,
                const_cast<T*>(y), ORDER, pexp, ORDER);
    for (unsigned int n=0; n<nnodes; ++n) exp[n] = pexp[perm[1][n]];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemhm(ORDER, ORDER, ORDER*ORDER, T(1.),
                const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
                exp, ORDER, pexp, ORDER);
    for (unsigned int n=0; n<nnodes; ++n) exp[perm[1][n]] = pexp[perm[2][n]];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemhm(ORDER, ORDER, ORDER*ORDER, T(1.),
                const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
                exp, ORDER, pexp, ORDER);
    // ORDER*ORDER*ORDER
    for (unsigned int n=0; n<nnodes; ++n) Y[perm[2][n]] += pexp[n];
  }


  // 3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes
  template <typename T>
  void l2l(const T *const S, const T *const X, T *const x) const
  {
    T exp[nnodes], pexp[nnodes];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemm(ORDER, ORDER, ORDER*ORDER, T(1.),
               const_cast<T*>(S), ORDER,
               const_cast<T*>(X), ORDER, pexp, ORDER);
    for (unsigned int n=0; n<nnodes; ++n) exp[n] = pexp[perm[1][n]];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemm(ORDER, ORDER, ORDER*ORDER, T(1.),
               const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
               exp, ORDER, pexp, ORDER);
    for (unsigned int n=0; n<nnodes; ++n) exp[perm[1][n]] = pexp[perm[2][n]];
    // ORDER * (2*ORDER-1) * ORDER
    blas::gemm(ORDER, ORDER, ORDER*ORDER, T(1.),
               const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
               exp, ORDER, pexp, ORDER);
    // ORDER*ORDER*ORDER
    for (unsigned int n=0; n<nnodes; ++n) x[perm[2][n]] += pexp[n];
  }

};



#endif

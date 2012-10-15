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

#ifndef tensor_hpp
#define tensor_hpp


#include <stdexcept>

#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>




// 1d
template <int ORDER>
static
const unsigned int tensor_idx_1d(const unsigned int d,
                                 const unsigned int idx)
{
  assert(idx<ORDER);
  assert(d<1);
  return idx;
}


// 2d
template <int ORDER>
static
const unsigned int tensor_idx_2d(const unsigned int d,
                                 const unsigned int idx)
{
  enum {nnodes = BasisTraits<ORDER,2>::nnodes};
  assert(idx<nnodes);
  switch (d) {
  case 0: return (idx%ORDER); break;
  case 1: return (idx/ORDER); break;
  default: 
    std::cerr << "Not in 2D: " << d << std::endl;
    exit(-1);
  }
}


// 3d
template <int ORDER>
static
const unsigned int tensor_idx_3d(const unsigned int d,
                                 const unsigned int idx)
{
  enum {nnodes = BasisTraits<ORDER,3>::nnodes};
  assert(idx<nnodes);
  switch (d) {
  case 0: return ( idx        % ORDER);  break;
  case 1: return ((idx/ORDER) % ORDER);  break;
  case 2: return ( idx/(ORDER * ORDER)); break;
  default: 
    std::cerr << "Not in 3D: " << d << std::endl;
    exit(-1);
  }
}






// tensor
template <int DIM, int ORDER>
class Tensor
{
  enum {nnodes = BasisTraits<ORDER,DIM>::nnodes};
  typedef boost::array<unsigned int[DIM],nnodes> node_id_array;
	typedef double weights_array[nnodes];

  static const node_id_array getBasisNodeIds()
  {
    node_id_array ids;
    for (unsigned int d=0; d<DIM; ++d)
      for (unsigned int n=0; n<nnodes; ++n)
        switch (DIM) {
        case 1: ids[n][d] = tensor_idx_1d<ORDER>(d,n); break;
        case 2: ids[n][d] = tensor_idx_2d<ORDER>(d,n); break;
        case 3: ids[n][d] = tensor_idx_3d<ORDER>(d,n); break;}
    return ids;
  }
  


public:
  const static node_id_array node_ids;



	/*!  Sets the roots of the Chebyshev quadrature weights defined as \f$w_i =
	  \frac{\pi}{\ell}\sqrt{1-\bar x_i^2}\f$ with the Chebyshev roots \f$\bar
	  x\f$.
	 
	  @param weights[out] the root of the weights \f$\sqrt{w_i}\f$
	 */
	static void setRootOfWeights(double weights[nnodes])
	{
		// weights in 1 dimension
		typedef Chebyshev<ORDER> Cheb;
		const double fracPiORDER = boost::math::constants::pi<double>() / ORDER;
		double weights1D[ORDER];
		for (unsigned int o=0; o<ORDER; ++o)
			weights1D[o] = fracPiORDER * sqrt(1.-Cheb::nodes[o]*Cheb::nodes[o]);
		
		// weights in DIM dimensions (tensor structure)
		for (unsigned int n=0; n<nnodes; ++n) {
			const unsigned int *const node_id = node_ids[n];
			double w = 1.;
			for (unsigned int d=0; d<DIM; ++d) w *= weights1D[d];
			weights[n] = sqrt(w);
		}
	}


};


template <int DIM, int ORDER>
const typename Tensor<DIM,ORDER>::node_id_array
Tensor<DIM,ORDER>::node_ids = Tensor<DIM,ORDER>::getBasisNodeIds();









#endif

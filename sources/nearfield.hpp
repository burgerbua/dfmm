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

#ifndef nearfield_hpp
#define nearfield_hpp


#include <boost/foreach.hpp>


//template <typename T>
//class NearField
//{
//  typedef std::vector<const T*> values_vector;
//  values_vector nf;
//
//public:
//  typedef T value_type;
//
//  explicit NearField()
//  { nf.clear(); }
//
//  void reserve(const unsigned int size)
//  { nf.reserve(size); }
//
//  void push_back(const T *const K)
//  { nf.push_back(K); }
//
//  const values_vector& get() const
//  { return nf; }
//
//  ~NearField()
//  {
//    for (typename values_vector::iterator it=nf.begin(); it!=nf.end(); ++it)
//      delete [] *it;
//    nf.clear();
//  }
//};
//
//
//
//template <typename cluster_type,
//					typename particle_type,
//          typename relations_type,
//          typename nearfield_type,
//          typename kernel_type>
//class NFsetter
//  : public std::unary_function<cluster_type, void>
//{
//  typedef typename nearfield_type::value_type value_type;
//
//  const kernel_type& kernel;
//  const std::vector<relations_type>& relations;
//  std::vector<nearfield_type>& nearfield;
//  const unsigned int nx, ny;
//	const unsigned int *const xidx, *const yidx;
//	const particle_type *const px, *const py; 
//  
//public:
//  explicit NFsetter(const kernel_type& _kernel,
//                    const std::vector<relations_type>& _relations,
//                    std::vector<nearfield_type>& _nearfield,
//                    const unsigned int _nx,
//										const unsigned int *const _xidx,
//										const particle_type *const _px,
//                    const unsigned int _ny,
//										const unsigned int *const _yidx,
//										const particle_type *const _py)
//    : kernel(_kernel), relations(_relations), nearfield(_nearfield),
//      nx(_nx), ny(_ny),	xidx(_xidx), yidx(_yidx), px(_px), py(_py) {}
//  
//  
//  void operator()(const cluster_type *const clx) const
//  {
//    const relations_type& rl = relations.at(clx->cidx);
//    nearfield_type& nf = nearfield.at(clx->lidx);
//    nf.reserve(rl.getNeighbors().size());
//    BOOST_FOREACH(const cluster_type *const cly, rl.getNeighbors()) {
//      value_type *const K = new value_type [clx->size * cly->size];
//      if (clx->cidx != cly->cidx)
//				for (unsigned int i=0; i<clx->size; ++i)
//					for (unsigned int j=0; j<cly->size; ++j)
//						K[j*clx->size+i] = kernel(px[xidx[clx->nbeg + i]].getPoint(),
//																			py[yidx[cly->nbeg + j]].getPoint());
////        kernel.compute(clx->size, xdofs+clx->nbeg,
////									 cly->size, ydofs+cly->nbeg, K);
//      assert(check_nan(clx->size * cly->size, K)); 
//      nf.push_back(K);
//    }
//  }
//  
//};
//
//      
//
//
//
//
//
//
//template <typename cluster_type,
//          typename relations_type,
//          typename nearfield_type>
//class NFapplyer
//  : public std::unary_function<cluster_type, void>
//{
//  typedef typename nearfield_type::value_type value_type;
//
//  const std::vector<relations_type>& relations;
//  const std::vector<nearfield_type>& nearfield;
//  value_type *const x;
//  const value_type *const y;
//
//public:
//  NFapplyer(const std::vector<relations_type>& _relations,
//            const std::vector<nearfield_type>& _nearfield,
//            value_type *const _x, const value_type *const _y)
//    : relations(_relations), nearfield(_nearfield), x(_x), y(_y)
//  {}
//
//  void operator()(cluster_type *const clx) const
//  {
//    const relations_type& rl = relations.at(clx->cidx);
//    const unsigned int size = rl.getNeighbors().size();
//    assert(size == nearfield.at(clx->lidx).get().size());
//    for (unsigned int i=0; i<size; ++i) {
//      const cluster_type *const cly = rl.getNeighbors().at(i); 
//      applyNearField(clx->size, cly->size,
//                     nearfield.at(clx->lidx).get().at(i),
//                     y + cly->nbeg, x + clx->nbeg);
//      assert(check_nan(clx->size, x+clx->nbeg)); 
//    }
//
//  }
//  
//  // near field
//  static void applyNearField(const unsigned int nx,
//                             const unsigned int ny,
//                             const value_type *const K,
//                             const value_type *const y,
//                             value_type *const x)
//  {
//    blas::gemva(nx, ny, 1., const_cast<value_type*>(K),
//                const_cast<value_type*>(y), x);
//  }
//};




template <typename cluster_type,
					typename particle_type,
          typename relations_type,
          typename kernel_type>
class NFadder
  : public std::unary_function<cluster_type, void>
{
  typedef typename kernel_type::value_type value_type;

  const kernel_type& kernel;
  const std::vector<relations_type>& relations;
	const unsigned int *const xidx, *const yidx;
	const particle_type *const px, *const py; 
  value_type *const x;
  const value_type *const y;

public:
	
  explicit NFadder(const kernel_type& _kernel,
									 const std::vector<relations_type>& _relations,
									 const unsigned int *const _xidx,
									 const particle_type *const _px,
									 const unsigned int *const _yidx,
									 const particle_type *const _py,
									 value_type *const _x,
									 const value_type *const _y)
		: kernel(_kernel), relations(_relations),
			xidx(_xidx), yidx(_yidx),
			px(_px), py(_py),
			x(_x), y(_y) {}

  void operator()(const cluster_type *const clx) const
  {
		const unsigned int xcidx = clx->getCidx();
		const unsigned int xsize = clx->getSize();
		const unsigned int xnbeg = clx->getNbeg();
    const relations_type& rl = relations.at(xcidx);
    BOOST_FOREACH(const cluster_type *const cly, rl.getNeighbors()) {
			const unsigned int ycidx = cly->getCidx();
			const unsigned int ysize = cly->getSize();
			const unsigned int ynbeg = cly->getNbeg();
      if (xcidx != ycidx)
				for (unsigned int i=0; i<xsize; ++i)
					for (unsigned int j=0; j<ysize; ++j)
						x[xnbeg+i] += kernel(px[xidx[xnbeg+i]].getPoint(),
																 py[yidx[ynbeg+j]].getPoint()) * y[ynbeg+j];
		}
	}
  
};







//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


#include <omp.h>



//template <typename cluster_type,
//          typename nearfield_cmp_type>
//void ComputeNearField(const std::vector<cluster_type*>& leafs,
//                      const nearfield_cmp_type& nfc)
//{
//  assert(leafs.front()->isleaf());
//  
//  std::cout << "\n- Evaluating near field operators of "
//            << leafs.size() << " leaf clusters " << std::flush;
//	const double t = omp_get_wtime();
//
//  // non parallel mode
//  std::for_each(leafs.begin(), leafs.end(), nfc);
//  
//	//  // parallel mode
//	//#pragma omp parallel for schedule(dynamic) // guided
//	//  for (typename std::vector<cluster_type*>::const_iterator
//	//         it=leafs.begin(); it < leafs.end(); ++it) nfc(*it);
//
//  std::cout << "took " << omp_get_wtime() - t << " s" << std::endl;
//}


template <typename cluster_type,
          typename nearfield_cmp_type>
void ApplyNearField(const std::vector<cluster_type*>& leafs,
										const nearfield_cmp_type& nfc)
{
  assert(leafs.front()->isLeaf());
  
  std::cout << "\n- Adding near field contributions of "
            << leafs.size() << " leaf clusters " << std::flush;
	const double t = omp_get_wtime();

#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
	for (int i=0; i<leafs.size(); ++i) nfc(leafs.at(i));
#else
	std::for_each(leafs.begin(), leafs.end(), nfc);
#endif

  std::cout << "took " << omp_get_wtime() - t << " s" << std::endl;
}



//template <typename cluster_type,
//          typename relations_type,
//          typename nearfield_type,
//          typename value_type>
//void AddNearFieldContribution(const std::vector<cluster_type*>& leafs,
//                              const std::vector<relations_type>& relations,
//                              const std::vector<nearfield_type>& nearfield,
//                              value_type *const y,
//                              value_type *const x)
//{
//  assert(leafs.front()->isleaf());
//  std::for_each(leafs.begin(), leafs.end(),
//                NFapplyer<cluster_type,relations_type,nearfield_type>
//                (relations, nearfield, x, y));
//}



#endif

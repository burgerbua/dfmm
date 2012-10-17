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



#endif

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

#ifndef l2lapplyer_tpl
#define l2lapplyer_tpl


#include <boost/foreach.hpp>

#include "sources/planewave.hpp"
#include "sources/tensorproduct.hpp"


/*! 
  \addtogroup mtranslation

  \par L2L: Local-To-Local Operators
  The L2L (and L2P) operators translate the local expansions from
  parant-to-child cluster levels. They differ between low-, respectively
  high-frequency regime.

*/




/*! \class L2Lapplyer
  Abstract base class for all local-to-local operators

  \brief Abstract base class for all L2L operators

  \ingroup mtranslation
 */
template <typename cluster_type,
          typename m2l_handler_type>
class L2Lapplyer
  : public std::unary_function<cluster_type, void>
{
  const m2l_handler_type& m2l;

protected:

  mutable FlopCounter fcount;


  L2Lapplyer(const m2l_handler_type& _m2l,
             const std::string& description)
    : m2l(_m2l), fcount(description, m2l.getLevel()->getNlevel())
  {}

  // does something only for Qu Ci Qb' convolution
  template <typename value_type>
  void expand(value_type *const X,
              const unsigned int c=0) const
  { m2l.applyL2Lpart(X, c); }

public:
  virtual void operator()(const cluster_type *const clx) const = 0;

  // l2l
  template <typename value_type>
  void applyL2L(const value_type *const S,
                const unsigned int nX,
                const value_type *const X,
                const unsigned int nx,
                value_type *const x) const
  {
    blas::gemva(nx, nX, 1., const_cast<value_type*>(S),
                const_cast<value_type*>(X), x);
    assert(check_nan(nx, x));

    // add flops
    fcount.addFlops(nx*(2*nX-1) + nx);
  }
};






/*! \class L2PapplyerLF
  Low-frequency (non-directional) local-to-particle operator (L2P)

  \brief L2P in the low-frequency regime 

  \ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
          typename expansionhandler_type,
          typename m2l_handler_type>
class L2PapplyerLF : L2Lapplyer<cluster_type,m2l_handler_type>
{
  enum {nnodes = clusterbasis_type::nnodes};
  typedef typename clusterbasis_type::value_type value_type;

  using L2Lapplyer<cluster_type,m2l_handler_type>::fcount;

  std::vector<clusterbasis_type>&     rbasis;
  std::vector<expansionhandler_type>& rexph;

  value_type *const x;

public:
  L2PapplyerLF(const m2l_handler_type& m2l,
                std::vector<clusterbasis_type>& _rbasis,
                std::vector<expansionhandler_type>& _rexph,
                value_type *const _x)
    : L2Lapplyer<cluster_type,m2l_handler_type>(m2l, "L2PapplyerLF"),
      rbasis(_rbasis), rexph(_rexph), x(_x) {}

  void operator()(const cluster_type *const clx) const
  {
    // asserts
    assert(clx->isLeaf());                            // leaf
    assert(rexph.at(clx->getCidx()).getExpansions().size()==1); // low freq

    // unique cluster index
    const unsigned int cidx = clx->getCidx();

    // does something only for Qu Ci Qb' convolution
    value_type *const X = rexph.at(cidx).getExpansions().at(0)->getvals();
    this->expand(X);

    this->applyL2L(rbasis.at(cidx).getS(), nnodes, X, clx->getSize(),
									 x+clx->getNbeg());
  }
};






/*! \class L2LapplyerLF
  Low-frequency (non-directional) local-to-local operator (L2L)

  \brief L2L in the low-frequency regime 

  \ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
          typename expansionhandler_type,
          typename m2l_handler_type>
class L2LapplyerLF : L2Lapplyer<cluster_type,m2l_handler_type>
{
  enum {DIM = cluster_type::dim,
        ORDER = clusterbasis_type::order,
        nboxes = cluster_type::nboxes,
        nnodes = clusterbasis_type::nnodes};
  typedef typename clusterbasis_type::value_type value_type;
  typedef typename expansionhandler_type::expansion_type expansion_type;

  using L2Lapplyer<cluster_type,m2l_handler_type>::fcount;

  TensorProductHandler<DIM,ORDER> phandler;

  std::vector<clusterbasis_type>&     rbasis;
  std::vector<expansionhandler_type>& rexph;

public:
  L2LapplyerLF(const m2l_handler_type& m2l,
               std::vector<clusterbasis_type>& _rbasis,
               std::vector<expansionhandler_type>& _rexph)
    : L2Lapplyer<cluster_type,m2l_handler_type>(m2l, "L2LapplyerLF"),
      rbasis(_rbasis), rexph(_rexph) {}

  void operator()(const cluster_type *const clx) const
  {
    // asserts
    assert(!clx->isLeaf());                           // no leaf
    assert(rexph.at(clx->getCidx()).getExpansions().size()==1); // low freq

    // unique cluster index
    const unsigned int cidx = clx->getCidx();

    // does something only for Qu Ci Qb' convolution
    expansion_type *const expansion = rexph.at(cidx).getExpansions().at(0);
    value_type *const X = expansion->getvals();
    this->expand(X);

    // apply l2l
    const value_type *const S = rbasis.at(cidx).getS();
    for (unsigned int b=0; b<nboxes; ++b) {
      value_type *const x = expansion->getcvals(b);
      if (x) {
        phandler.l2l(S + b * 3*ORDER*ORDER, X, x);
        fcount.addFlops(3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes);
      }
    }

  }
  
};






/*! \class L2PapplyerHF
  High-frequency (directional) local-to-particle operator (L2P)

  \brief L2P in the high-frequency regime 

  \ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
          typename expansionhandler_type,
          typename m2l_handler_type,
          typename particle_type>
class L2PapplyerHF
  : L2Lapplyer<cluster_type,m2l_handler_type>
{
  enum {DIM = cluster_type::dim,
        ORDER = clusterbasis_type::order,
        nnodes = clusterbasis_type::nnodes};
  typedef typename cluster_type::point_type              point_type;
  typedef typename clusterbasis_type::value_type         value_type;

  typedef value_type T;

  typedef typename expansionhandler_type::expansion_type expansion_type;
  typedef typename expansionhandler_type::expansion_map  expansion_map;

  using L2Lapplyer<cluster_type,m2l_handler_type>::fcount;

  const PlaneWave<         DIM,ORDER,T>&                pwave;
  const PlaneWaveParticles<DIM,ORDER,T,particle_type>& ppwave;

  std::vector<clusterbasis_type>& rbasis;
  std::vector<expansionhandler_type>& rexph;

  value_type *const x;

public:
  explicit
  L2PapplyerHF(const m2l_handler_type& m2l,
               std::vector<clusterbasis_type>& _rbasis,
               std::vector<expansionhandler_type>& _rexph,
               value_type *const _x,
               const PlaneWave<DIM,ORDER,T>& _pwave,
               const PlaneWaveParticles<DIM,ORDER,T,particle_type>& _ppwave)
    : L2Lapplyer<cluster_type,m2l_handler_type>(m2l, "L2PapplyerHF"),
      pwave(_pwave), ppwave(_ppwave), rbasis(_rbasis), rexph(_rexph), x(_x) {}

  void operator()(const cluster_type *const clx) const
  {
    // assert
    assert(clx->isLeaf()); // leaf

    // unique cluster index and center
    const unsigned int cidx  = clx->getCidx();
    const point_type& center = clx->getCenter();

    // 1. Expand, apply right plane wave and copy //////////
    expansion_map& rexp = rexph.at(cidx).getExpansions();
    const unsigned int nexpansions = rexp.size();
    value_type *const X_all = new value_type [nnodes * nexpansions];
    unsigned int counter = 0;
    BOOST_FOREACH(typename expansion_map::value_type& ep, rexp) {
      expansion_type *const expansion = ep.second;
      
      // does something only for Qu Ci Qb' convolution
      const unsigned int c = ep.first;
      value_type *const X = expansion->getvals();
      this->expand(X, c);
      
      // apply plane wave
      const value_type scal = pwave.plane_wave(center, expansion->getu());
      blas::hadd(nnodes, scal, pwave.getppw(c), X);
      fcount.addFlops(7 + 2*nnodes);
      blas::copy(nnodes, X, X_all + (counter++)*nnodes);
    }
    
    // 2. Apply l2l ////////////////////////////////////////
    const unsigned int size = clx->getSize();
    value_type *const x_all = new value_type [size * nexpansions];
    blas::gemm(size, nnodes, nexpansions, 1.,
               rbasis.at(cidx).getS(), size, X_all, nnodes, x_all, size);
    fcount.addFlops(size*(2*nnodes-1)*nexpansions);
    delete [] X_all;
    
    // 3. Apply plane wave and add contribution ////////////
    value_type *const cpw = new value_type [size];
    counter = 0;
    BOOST_FOREACH(typename expansion_map::value_type& ep, rexp) {
      expansion_type *const expansion = ep.second;
      ppwave.setcpw(clx, expansion->getu(), cpw);
      blas::hadm(size, 1., cpw, x_all+counter*size);
      blas::add(size, x_all+counter*size, x+clx->getNbeg());
      fcount.addFlops(10*size);
      counter++;
    }
    delete [] cpw;
    delete [] x_all;

  }

};





/*!
  \class L2LapplyerHF
  High-frequency (directional) local-to-local operator (L2L)

  \brief L2L in the high-frequency regime 

  \ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
          typename expansionhandler_type,
          typename m2l_handler_type>
class L2LapplyerHF
  : L2Lapplyer<cluster_type,m2l_handler_type>
{
  enum {DIM = cluster_type::dim,
        ORDER = clusterbasis_type::order,
        nboxes = cluster_type::nboxes,
        nnodes = clusterbasis_type::nnodes};

  typedef typename cluster_type::point_type              point_type;
  typedef typename clusterbasis_type::value_type         value_type;

  typedef value_type T;

  typedef typename expansionhandler_type::expansion_type expansion_type;
  typedef typename expansionhandler_type::expansion_map  expansion_map;

  using L2Lapplyer<cluster_type,m2l_handler_type>::fcount;

  TensorProductHandler<DIM,ORDER> phandler;

  const PlaneWave<DIM,ORDER,T>& pwave;

  std::vector<clusterbasis_type>& rbasis;
  std::vector<expansionhandler_type>& rexph;

public:
  L2LapplyerHF(const m2l_handler_type& m2l,
               std::vector<clusterbasis_type>& _rbasis,
               std::vector<expansionhandler_type>& _rexph,
               const PlaneWave<DIM,ORDER,T>& _pwave)
    : L2Lapplyer<cluster_type,m2l_handler_type>(m2l, "L2LapplyerHF"),
      pwave(_pwave), rbasis(_rbasis), rexph(_rexph)
  {}

  void operator()(const cluster_type *const clx) const
  {
    // assert    
    assert(!clx->isLeaf()); // no leaf

    // unique cluster index and center
    const unsigned int cidx = clx->getCidx();
    const point_type& center = clx->getCenter();

    // 1. Expand, apply right plane wave and copy //////////
    expansion_map& rexp = rexph.at(cidx).getExpansions();
    const unsigned int nexpansions = rexp.size();
    unsigned int counter[nboxes];
    value_type* X_all[  nboxes];
    for (unsigned int b=0; b<nboxes; ++b) {
      X_all[b] = new value_type [nnodes * nexpansions];
      counter[b] = 0;
    }
    BOOST_FOREACH(typename expansion_map::value_type& ep, rexp) {
      expansion_type *const expansion = ep.second;
      
      // does something only for Qu Ci Qb' convolution
      const unsigned int c = ep.first;
      value_type *const X = expansion->getvals();
      this->expand(X, c);
      
      // apply plane wave
      const point_type& u = expansion->getu();
      const value_type scal = pwave.plane_wave(center, u);
      blas::hadd(nnodes, scal, pwave.getppw(c), X);
      fcount.addFlops(2*nnodes + 7);

      // copy only needed ones
      for (unsigned int b=0; b<nboxes; ++b) {
        value_type *const x_ = expansion->getcvals(b);
        if (x_) blas::copy(nnodes, X, X_all[b] + (counter[b]++) * nnodes);
      }
    }

    // 2. Apply l2l ////////////////////////////////////////
    value_type* x_all[nboxes];
    value_type *const S = rbasis.at(cidx).getS();
    for (unsigned int b=0; b<nboxes; ++b) {
      x_all[b] = new value_type [nnodes * counter[b]];
      phandler.l2l(S + b * 3*ORDER*ORDER, counter[b], X_all[b], x_all[b]);
      fcount.addFlops(3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes);
      delete [] X_all[b];
      counter[b] = 0;
    }

    // 3. Apply plane wave and add contribution ////////////
    BOOST_FOREACH(typename expansion_map::value_type& ep, rexp) {
      expansion_type *const expansion = ep.second; 
      const unsigned int c = ep.first;
      for (unsigned int b=0; b<nboxes; ++b) {
        value_type *const x_ = expansion->getcvals(b);
        if (x_) {
          const value_type *const cpw = pwave.getcpw(c);
          const value_type scal = pwave.plane_wave(center, expansion->getu());
          blas::hadm(nnodes, scal, cpw+b*nnodes, x_all[b]+counter[b]*nnodes);
          blas::add(nnodes, x_all[b]+counter[b]*nnodes, x_);
          fcount.addFlops(7 + 3*nnodes);
          counter[b]++;
        }
      }
    }

    // free memory
    for(unsigned int b=0; b<nboxes; ++b)
      delete [] x_all[b];
  }

};








//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


/*! 
  \ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
          typename expansionhandler_type,
          typename m2l_handler_type,
          typename particle_type>
void ApplyL2Loperators
(const std::vector<std::vector<cluster_type*> >& rclvv,
 std::vector<clusterbasis_type>& rbasis,
 std::vector<expansionhandler_type>& rexph,
 std::vector<m2l_handler_type>& m2lv,
 typename clusterbasis_type::value_type *const x,
 const particle_type *const particles,
 const unsigned int  *const pindices,
 const double wavenum)
{
  enum {DIM = cluster_type::dim,
        ORDER = clusterbasis_type::order};
  typedef typename clusterbasis_type::value_type T;

  typedef L2PapplyerLF<cluster_type,clusterbasis_type,
                        expansionhandler_type,m2l_handler_type>
    l2l_lfl_type;
  typedef L2LapplyerLF< cluster_type,clusterbasis_type,
                        expansionhandler_type,m2l_handler_type>
    l2l_lf_type;
  typedef L2PapplyerHF<cluster_type,clusterbasis_type,
                       expansionhandler_type,m2l_handler_type,particle_type>
    l2l_hfl_type;
  typedef L2LapplyerHF< cluster_type,clusterbasis_type,
                        expansionhandler_type,m2l_handler_type>
    l2l_hf_type;
  
  for (unsigned int l=0; l<rclvv.size(); ++l) {
    const std::vector<cluster_type*>& rclv = rclvv.at(l);
    m2l_handler_type& m2l = m2lv.at(l);
    if (m2l.hasInteractions()) {
      
      const double t = omp_get_wtime();
      std::cout << "  at " << l << std::flush;

      switch (m2l.getLevel()->is_in_high_frequency_regime()) {
      case true: {

        const PlaneWave<DIM,ORDER,T>
          pwave(wavenum, m2l.getLevel());

        switch (m2l.getLevel()->isLeaf()) {
        case true: {

          const PlaneWaveParticles<DIM,ORDER,T,particle_type>
            ppwave(wavenum, particles, pindices);

          l2l_hfl_type l2l(m2l, rbasis, rexph, x, pwave, ppwave);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
          for (int i=0; i<rclv.size(); ++i) l2l(rclv.at(i));
#else
          std::for_each(rclv.begin(), rclv.end(), l2l);
#endif
        } break;
        case false: {
          l2l_hf_type l2l(m2l, rbasis, rexph, pwave);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
          for (int i=0; i<rclv.size(); ++i) l2l(rclv.at(i));
#else
          std::for_each(rclv.begin(), rclv.end(), l2l);
#endif
        } break; }
      } break;
      case false: {
        switch (m2l.getLevel()->isLeaf()) {
        case true: {
          l2l_lfl_type l2l(m2l, rbasis, rexph, x);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
          for (int i=0; i<rclv.size(); ++i) l2l(rclv.at(i));
#else
          std::for_each(rclv.begin(), rclv.end(), l2l);
#endif
        } break;
        case false: {
          l2l_lf_type l2l(m2l, rbasis, rexph);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
          for (int i=0; i<rclv.size(); ++i) l2l(rclv.at(i));
#else
          std::for_each(rclv.begin(), rclv.end(), l2l);
#endif
        } break; }
      } break;
      }
      
      std::cout << " took " << omp_get_wtime() - t << " s" << std::endl;

    }
  }

}






/*! 
  \ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
          typename expansionhandler_type,
          typename m2l_handler_type>
void ApplyL2Loperators
(const std::vector<std::vector<cluster_type*> >& rclvv,
 std::vector<clusterbasis_type>& rbasis,
 std::vector<expansionhandler_type>& rexph,
 std::vector<m2l_handler_type>& m2lv,
 typename clusterbasis_type::value_type *const x)
{
  enum {DIM = cluster_type::dim,
        ORDER = clusterbasis_type::order};
  typedef typename clusterbasis_type::value_type T;

  typedef L2PapplyerLF<cluster_type,clusterbasis_type,
                        expansionhandler_type,m2l_handler_type> l2l_lfl_type;
  typedef L2LapplyerLF< cluster_type,clusterbasis_type,
                        expansionhandler_type,m2l_handler_type> l2l_lf_type;
  
  for (unsigned int l=0; l<rclvv.size(); ++l) {
    const std::vector<cluster_type*>& rclv = rclvv.at(l);
    m2l_handler_type& m2l = m2lv.at(l);
    if (m2l.hasInteractions()) {
      
      const double t = omp_get_wtime();
      std::cout << "  at " << l << std::flush;

        switch (m2l.getLevel()->isLeaf()) {
        case true: {
          l2l_lfl_type l2l(m2l, rbasis, rexph, x);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
          for (int i=0; i<rclv.size(); ++i) l2l(rclv.at(i));
#else
          std::for_each(rclv.begin(), rclv.end(), l2l);
#endif
        } break;
        case false: {
          l2l_lf_type l2l(m2l, rbasis, rexph);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
          for (int i=0; i<rclv.size(); ++i) l2l(rclv.at(i));
#else
          std::for_each(rclv.begin(), rclv.end(), l2l);
#endif
        } break; }
      
      std::cout << " took " << omp_get_wtime() - t << " s" << std::endl;

    }
  }

}




#endif

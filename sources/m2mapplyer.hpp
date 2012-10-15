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

#ifndef m2mapplyer_tpl
#define m2mapplyer_tpl


#include <boost/foreach.hpp>

#include "sources/planewave.hpp"
#include "sources/permutations.hpp"


/*! 
	\addtogroup mtranslation

	\par M2M: Multipole-To-Multipole Operators
	The M2M (and P2M) operators translate the multipole expansions from
	child-to-parant cluster levels. They differ between low-, respectively
	high-frequency regime.

*/



/*! \class M2Mapplyer
	Abstract base class for all multipole-to-multipole operators

	\brief Abstract base class for all M2M operators

	\ingroup mtranslation
 */
template <typename cluster_type,
          typename m2l_handler_type>
class M2Mapplyer
  : public std::unary_function<cluster_type, void>
{
  const m2l_handler_type& m2l;


protected:
	mutable FlopCounter fcount;

  M2Mapplyer(const m2l_handler_type& _m2l,
						 const std::string& description)
		: m2l(_m2l), fcount(description, m2l.getLevel()->getNlevel())
	{}

  // does something only for Qu Ci Qb' convolution
  template <typename value_type>
  void compress(value_type *const Y,
                const unsigned int c=0) const
  { m2l.applyM2Mpart(Y, c); }

public:
  virtual void operator()(const cluster_type *const cly) const = 0;

  template <typename value_type>
  void applyM2M(const value_type *const S,
								const unsigned int ny,
								const value_type *const y,
								const unsigned int nY,
								value_type *const Y) const
  {
    blas::gemhv(ny, nY, 1., const_cast<value_type*>(S),
                const_cast<value_type*>(y), Y);

		// add flops
		fcount.addFlops(nY*(2*ny-1));
  }
};






/*! \class P2MapplyerLF
	Low-frequency (non-directional) particle-to-multipole operator (P2M)

	\brief P2M in the low-frequency regime 

	\ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
					typename expansionhandler_type,
          typename m2l_handler_type>
class P2MapplyerLF
  : M2Mapplyer<cluster_type, m2l_handler_type>
{
  enum {nnodes = clusterbasis_type::nnodes};
  typedef typename clusterbasis_type::value_type value_type;

	using M2Mapplyer<cluster_type,m2l_handler_type>::fcount;

  std::vector<clusterbasis_type>&     cbasis;
  std::vector<expansionhandler_type>& cexph;

  const value_type *const y;

public:
  explicit P2MapplyerLF(const m2l_handler_type& m2l,
                         std::vector<clusterbasis_type>&     _cbasis,
                         std::vector<expansionhandler_type>& _cexph,
                         const value_type *const _y)
    : M2Mapplyer<cluster_type,m2l_handler_type>(m2l, "P2MapplyerLF"),
      cbasis(_cbasis), cexph(_cexph), y(_y) {}

  void operator()(const cluster_type *const cly) const
  {
		// asserts
    assert(cly->isLeaf());                                    // leaf 
    assert(cexph.at(cly->getCidx()).getExpansions().size()==1); // low freq

		// unique cluster index
		const unsigned int cidx = cly->getCidx();

    // apply M2M
    value_type *const Y = cexph.at(cidx).getExpansions().at(0)->getvals();
    applyM2M(cbasis.at(cidx).getS(),cly->getSize(),y+cly->getNbeg(),nnodes,Y);
    assert(check_nan(nnodes, Y)); 

    // do something only for Qu Ci Qb' convolution
    compress(Y);
    assert(check_nan(nnodes, Y+nnodes)); 
  }
};






/*! \class M2MapplyerLF
	Low-frequency (non-directional) multipole-to-multipole operator (M2M)

	\brief M2M in the low-frequency regime 

	\ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
					typename expansionhandler_type,
          typename m2l_handler_type>
class M2MapplyerLF
  : M2Mapplyer<cluster_type,m2l_handler_type>
{
  enum {DIM = cluster_type::dim,
				ORDER = clusterbasis_type::order,
				nboxes = cluster_type::nboxes,
        nnodes = clusterbasis_type::nnodes};

  typedef typename clusterbasis_type::value_type value_type;

  typedef typename expansionhandler_type::expansion_type expansion_type;

	using M2Mapplyer<cluster_type,m2l_handler_type>::fcount;

	PermutationHandler<DIM,ORDER> phandler;

  std::vector<clusterbasis_type>& cbasis;
  std::vector<expansionhandler_type>& cexph;

public:
  M2MapplyerLF(const m2l_handler_type& m2l,
               std::vector<clusterbasis_type>& _cbasis,
							 std::vector<expansionhandler_type>& _cexph)
    : M2Mapplyer<cluster_type,m2l_handler_type>(m2l, "M2MapplyerLF"),
      cbasis(_cbasis), cexph(_cexph)
  {}

  void operator()(const cluster_type *const cly) const
  {
		// asserts
    assert(!cly->isLeaf());                                     // no leaf  
    assert(cexph.at(cly->getCidx()).getExpansions().size()==1); // low freq

		// unique cluster index
		const unsigned int cidx = cly->getCidx();

    // apply m2m
    const value_type *const S = cbasis.at(cidx).getS();
    expansion_type *const expansion = cexph.at(cidx).getExpansions().at(0);
    value_type *const Y = expansion->getvals();
    for (unsigned int b=0; b<nboxes; ++b) {
      value_type *const y = expansion->getcvals(b);
      if (y) {
				phandler.m2m(S + b * 3*ORDER*ORDER, y, Y);
				fcount.addFlops(3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes);
			}
    }
    assert(check_nan(nnodes, Y)); 

    // do something only for Qu Ci Qb' convolution
    compress(Y);
    assert(check_nan(nnodes, Y+nnodes)); 
  }

};








/*! \class P2MapplyerHF
	High-frequency (directional) particle-to-multipole operator (P2M)

	\brief P2M in the high-frequency regime 

	\ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
					typename expansionhandler_type,
          typename m2l_handler_type,
					typename particle_type>
class P2MapplyerHF
  : M2Mapplyer<cluster_type,m2l_handler_type>
{
  enum {DIM = cluster_type::dim,
				ORDER = clusterbasis_type::order,
				nnodes = clusterbasis_type::nnodes};

  typedef typename cluster_type::point_type              point_type;
  typedef typename clusterbasis_type::value_type         value_type;

	typedef value_type T;

  typedef typename expansionhandler_type::expansion_type expansion_type;
  typedef typename expansionhandler_type::expansion_map  expansion_map;

	using M2Mapplyer<cluster_type,m2l_handler_type>::fcount;

	const PlaneWave<         DIM,ORDER,T>&                pwave;
	const PlaneWaveParticles<DIM,ORDER,T,particle_type>& ppwave;

  std::vector<clusterbasis_type>&     cbasis;
  std::vector<expansionhandler_type>& cexph;

  const value_type *const y;

public:
  explicit
	P2MapplyerHF(const m2l_handler_type& m2l,
							 std::vector<clusterbasis_type>& _cbasis,
							 std::vector<expansionhandler_type>& _cexph,
							 const value_type *const _y,
							 const PlaneWave<         DIM,ORDER,T>&  _pwave,
							 const PlaneWaveParticles<DIM,ORDER,T,particle_type>& _ppwave)
    : M2Mapplyer<cluster_type,m2l_handler_type>(m2l, "P2MapplyerHF"),
			pwave(_pwave), ppwave(_ppwave), cbasis(_cbasis), cexph(_cexph), y(_y) {}

  void operator()(const cluster_type *const cly) const
  {
		// assert
    assert(cly->isLeaf()); // leaf

		// unique cluster index
		const unsigned int cidx = cly->getCidx();

		const unsigned int size = cly->getSize();


		expansion_map& cexp = cexph.at(cidx).getExpansions();
		value_type *const y_ = const_cast<value_type *const>(y+cly->getNbeg());

		// 1. copy and apply right plane wave ////////
		value_type *const cpw = new value_type [size];
		const unsigned int nexpansions = cexp.size();
    value_type *const y_all = new value_type [size * nexpansions];
		unsigned int cntr_y = 0;
    BOOST_FOREACH(typename expansion_map::value_type& pY, cexp) {
      expansion_type *const expansion = pY.second;
      blas::copy(size, y_,                      y_all + cntr_y*size);
			ppwave.setcpw(cly, expansion->getu(), cpw);
			blas::hadd(size, 1., cpw, y_all + cntr_y*size);
			fcount.addFlops(9*size); // computing and had-product of cpw
			cntr_y++;
		}
		delete [] cpw;

		// 2. apply m2m //////////////////////////////
		value_type *const Y_all = new value_type [nnodes * nexpansions];
		blas::gemhm(size, nnodes, nexpansions, 1.,
								cbasis.at(cidx).getS(), size, y_all, size, Y_all, nnodes);
		fcount.addFlops(size * (2*nnodes) * nexpansions);
    delete [] y_all; // free memory
		
		// 3. apply left plane wave and compress /////
		cntr_y = 0;
    BOOST_FOREACH(typename expansion_map::value_type& pY, cexp) {
			
      // apply plane wave
      expansion_type *const expansion = pY.second;
      const unsigned int c = pY.first;
      value_type *const Y = expansion->getvals();
			blas::copy(nnodes, Y_all + (cntr_y++)*nnodes, Y);
			const point_type& u = expansion->getu();
			const point_type& center = cly->getCenter();
			const value_type scal = pwave.plane_wave(center, u);
      blas::hadm(nnodes, scal, pwave.getppw(c), Y);
			fcount.addFlops(2*nnodes + 7); // had-product of ppw
      
      // do something only for Qu Ci Qb' convolution
      compress(Y, c);
		}
		delete [] Y_all; // free memory
   
  }
  
};




/*! \class M2MapplyerHF
	High-frequency (directional) multipole-to-multipole operator (M2M)
	
	\brief M2M in the high-frequency regime
	
	\ingroup mtranslation
 */
template <typename cluster_type,
          typename clusterbasis_type,
					typename expansionhandler_type,
          typename m2l_handler_type>
class M2MapplyerHF
  : M2Mapplyer<cluster_type,m2l_handler_type>
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

	using M2Mapplyer<cluster_type,m2l_handler_type>::fcount;

	PermutationHandler<DIM,ORDER> phandler;

	const PlaneWave<DIM,ORDER,T>& pwave;

  std::vector<clusterbasis_type>&     cbasis;
  std::vector<expansionhandler_type>& cexph;


public:
  M2MapplyerHF(const m2l_handler_type& m2l,
               std::vector<clusterbasis_type>& _cbasis,
							 std::vector<expansionhandler_type>& _cexph,
							 const PlaneWave<DIM,ORDER,T>& _pwave)
    : M2Mapplyer<cluster_type,m2l_handler_type>(m2l, "M2MapplyerHF"),
			pwave(_pwave), cbasis(_cbasis), cexph(_cexph)
  {}

  void operator()(const cluster_type *const cly) const
  {
		// assert
    assert(!cly->isLeaf()); // no leaf

		// unique cluster index and center
		const unsigned int cidx  = cly->getCidx();
		const point_type& center = cly->getCenter();

		expansion_map& cexp = cexph.at(cidx).getExpansions();

		// init box counter and child sources
		unsigned int cntr_y[nboxes];
		value_type* y_all[  nboxes];
		for (unsigned int b=0; b<nboxes; ++b) {
			cntr_y[b] = 0;
			y_all[ b] = new value_type [nnodes * cexp.size()];
		}

		// 1. gather child sources: copy and add right plane wave ////////
    BOOST_FOREACH(typename expansion_map::value_type& pY, cexp) {
      expansion_type *const expansion = pY.second;
			const unsigned int c = pY.first;
			const value_type scal = pwave.plane_wave(center, expansion->getu());
			fcount.addFlops(7);
      const value_type *const cpw = pwave.getcpw(c);
      for (unsigned int b=0; b<nboxes; ++b) {
        value_type *const y_ = expansion->getcvals(b);
        if (y_) {
					blas::copy(nnodes, y_,                   y_all[b] + cntr_y[b]*nnodes);
					blas::hadd(nnodes, scal, cpw + b*nnodes, y_all[b] + cntr_y[b]*nnodes);
					fcount.addFlops(2*nnodes);
					cntr_y[b]++;
        }
      }
    }

		// 2. anterpolate child sources //////////////////////////////////
		value_type *const S = cbasis.at(cidx).getS();
		value_type* Y_all[nboxes];
		for(unsigned int b=0; b<nboxes; ++b) {	
			Y_all[b] = new value_type [nnodes * cntr_y[b]];
			phandler.m2m(S + b * 3*ORDER*ORDER, cntr_y[b], y_all[b], Y_all[b]);
			fcount.addFlops(3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes);
			delete [] y_all[b]; // free memory
			cntr_y[b] = 0;      // set box counter to 0
		}

		// 3. apply left plane wave //////////////////////////////////////
    BOOST_FOREACH(typename expansion_map::value_type& pY, cexp) {
      expansion_type *const expansion = pY.second;

			// gather contributions from children
      value_type *const Y = expansion->getvals();
			for (unsigned int b=0; b<nboxes; ++b) {
        value_type *const y_ = expansion->getcvals(b);
        if (y_)	{
					blas::add(nnodes, Y_all[b] + (cntr_y[b]++) * nnodes, Y);
					fcount.addFlops(nnodes);
				}
      }

      // apply plane wave
      const value_type *const ppw = pwave.getppw(pY.first);
			const value_type scal = pwave.plane_wave(center, expansion->getu());
      blas::hadm(nnodes, scal, ppw, Y);
			fcount.addFlops(2*nnodes);

      // do something only for Qu Ci Qb' convolution
      const unsigned int c = pY.first;
      compress(Y, c);
    }

		// free memory
		for(unsigned int b=0; b<nboxes; ++b) delete [] Y_all[b];

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
					typename particle_type,
          typename clusterbasis_type,
					typename expansionhandler_type,
          typename m2l_handler_type>
void ApplyM2Moperators
(const std::vector<std::vector<cluster_type*> >& cclvv,
 std::vector<clusterbasis_type>& cbasis,
 std::vector<expansionhandler_type>& cexph,
 std::vector<m2l_handler_type>& m2lv,
 const typename clusterbasis_type::value_type *const y,
 const particle_type *const particles,
 const unsigned int  *const pindices,
 const double wavenum)
{
	enum {DIM = cluster_type::dim,
				ORDER = clusterbasis_type::order};
	typedef typename clusterbasis_type::value_type T;

  typedef P2MapplyerLF<cluster_type,clusterbasis_type,
											 expansionhandler_type,m2l_handler_type>
		m2m_lfl_type;
  typedef M2MapplyerLF<cluster_type,clusterbasis_type,
											 expansionhandler_type,m2l_handler_type>
		m2m_lf_type;
  typedef P2MapplyerHF<cluster_type,clusterbasis_type,
											 expansionhandler_type,m2l_handler_type,particle_type>
		m2m_hfl_type;
  typedef M2MapplyerHF<cluster_type,clusterbasis_type,
											 expansionhandler_type,m2l_handler_type>
		m2m_hf_type;

  int l = cclvv.size();
  while (--l>=0) {
    const std::vector<cluster_type*>& cclv = cclvv.at(l);
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
					
					m2m_hfl_type m2m(m2l, cbasis, cexph, y, pwave, ppwave);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
					for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
#else
          std::for_each(cclv.begin(), cclv.end(), m2m);
#endif
				} break;
				case false: {
					m2m_hf_type m2m(m2l, cbasis, cexph, pwave);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
					for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
#else
          std::for_each(cclv.begin(), cclv.end(), m2m);
#endif
				} break; }
			} break;
			case false: {
				switch (m2l.getLevel()->isLeaf()) {
				case true: {
					m2m_lfl_type m2m(m2l, cbasis, cexph, y);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
					for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
#else
          std::for_each(cclv.begin(), cclv.end(), m2m);
#endif
				} break;
				case false: {
					m2m_lf_type m2m(m2l, cbasis, cexph);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
					for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
#else
          std::for_each(cclv.begin(), cclv.end(), m2m);
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
void ApplyM2Moperators
(const std::vector<std::vector<cluster_type*> >& cclvv,
 std::vector<clusterbasis_type>& cbasis,
 std::vector<expansionhandler_type>& cexph,
 std::vector<m2l_handler_type>& m2lv,
 const typename clusterbasis_type::value_type *const y)
{
	enum {DIM = cluster_type::dim,
				ORDER = clusterbasis_type::order};
	typedef typename clusterbasis_type::value_type T;

  typedef P2MapplyerLF<cluster_type,clusterbasis_type,
											 expansionhandler_type,m2l_handler_type> m2m_lfl_type;
  typedef M2MapplyerLF<cluster_type,clusterbasis_type,
											 expansionhandler_type,m2l_handler_type> m2m_lf_type;

  int l = cclvv.size();
  while (--l>=0) {
    const std::vector<cluster_type*>& cclv = cclvv.at(l);
    m2l_handler_type& m2l = m2lv.at(l);
		if (m2l.hasInteractions()) {

			const double t = omp_get_wtime();
			std::cout << "  at " << l << std::flush;

				switch (m2l.getLevel()->isLeaf()) {
				case true: {
					m2m_lfl_type m2m(m2l, cbasis, cexph, y);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
					for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
#else
          std::for_each(cclv.begin(), cclv.end(), m2m);
#endif
				} break;
				case false: {
					m2m_lf_type m2m(m2l, cbasis, cexph);
#ifdef DFMM_USE_OMP
#pragma omp parallel for schedule(dynamic) // guided
					for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
#else
          std::for_each(cclv.begin(), cclv.end(), m2m);
#endif
				} break; }

			std::cout << " took " << omp_get_wtime() - t << " s" << std::endl;

		}

	}

}


#endif



//  int l = cclvv.size();
//  while (--l>=0) {
//    const std::vector<cluster_type*>& cclv = cclvv.at(l);
//    m2l_handler_type& m2l = m2lv.at(l);
//		if (m2l.hasInteractions()) {
//      switch (m2l.getLevel()->isLeaf()) {
//      case true:
//        switch (m2l.getLevel()->is_in_high_frequency_regime()) {
//        case true: {
//          m2m_hfl_type m2m(m2l, cbasis, cexph, y, wavenum);
////#pragma omp parallel for schedule(dynamic) // guided
////          for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
//          std::for_each(cclv.begin(), cclv.end(), m2m);
//        } break;
//        case false: {
//          m2m_lfl_type m2m(m2l, cbasis, cexph, y);
////#pragma omp parallel for schedule(dynamic) // guided
////          for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
//          std::for_each(cclv.begin(), cclv.end(), m2m);
//        } break;
//        } break;
//      case false: 
//        switch (m2l.getLevel()->is_in_high_frequency_regime()) {
//        case true: {
//          m2m_hf_type m2m(m2l, cbasis, cexph, wavenum);
////#pragma omp parallel for schedule(dynamic) // guided
////          for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
//          std::for_each(cclv.begin(), cclv.end(), m2m);
//        } break;
//        case false: {
//          m2m_lf_type m2m(m2l, cbasis, cexph);
////#pragma omp parallel for schedule(dynamic) // guided
////          for (int i=0; i<cclv.size(); ++i) m2m(cclv.at(i));
//          std::for_each(cclv.begin(), cclv.end(), m2m);
//        } break;
//        } break;
//
//      }
//    }
//  }

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

#ifndef m2lapplyer_tpl
#define m2lapplyer_tpl


#include <boost/foreach.hpp>




template <int> class Cluster;
template <int> class Relations;
template <int,int,typename> class ExpansionHandler;



/*! 
	\defgroup mm2l M2L: Multipole-To-Local Operators

	The M2L operator is one of the central components of the FMM (besides L2L
	and M2M).

*/


/*!
	\ingroup mm2l

	\class M2Lapplyer 
	Multipole-to-local operator (M2L)

	\brief Applies all M2L operators 
	
  The M2L applyer does not differ between low- and high-frequency regime as
  opposed to the L2L- and M2M-applyer.

 */
template <int DIM,
					int ORDER,
					typename T,
          typename m2lhandler_type>
class M2Lapplyer
  : public std::unary_function<Cluster<DIM>, void>
{
  enum {nnodes = BasisTraits<ORDER,DIM>::nnodes};
	typedef ExpansionHandler<DIM,ORDER,T> expansionhandler_type;
	typedef Relations<DIM>                       relations_type;

  typedef typename expansionhandler_type::value_type value_type;
	typedef typename expansionhandler_type::expansion_map expansion_map;
  typedef typename relations_type::cluster_pair_vector_map ilist_map;
  typedef typename ilist_map::value_type                   ilist_pair;
  typedef typename ilist_map::mapped_type                  ilist_pair_vec;
  
  std::vector<expansionhandler_type>& rexph;
  const std::vector<expansionhandler_type>& cexph;

  const std::vector<relations_type>& relations;

  const m2lhandler_type& m2l;

  
public:
  M2Lapplyer(std::vector<expansionhandler_type>& _rexph,
             const std::vector<expansionhandler_type>& _cexph,
             const std::vector<relations_type>& _relations,
             const m2lhandler_type& _m2l)
    : rexph(_rexph), cexph(_cexph), relations(_relations), m2l(_m2l)
  {}



  void operator()(const Cluster<DIM> *const clx)
  {
		expansion_map& rexp	= rexph.at(clx->getCidx()).getExpansions();
    const relations_type& rl = relations.at(clx->getCidx());
    BOOST_FOREACH(const ilist_pair& ilp, rl.getInteractions()) {
      const unsigned int c = ilp.first;

      value_type *const X = rexp.at(c)->getvals();

      BOOST_FOREACH(const typename ilist_pair_vec::value_type& ip,
                    ilp.second) {
        value_type *const Y
          = cexph.at((ip.second)->getCidx()).getExpansions().at(c)->getvals();
				m2l.apply(ip.first, Y, X, c);
      }
			
    }
    
  }


};





//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#include <omp.h>


/*!
	\ingroup mm2l
 */
template <int DIM,
					int ORDER,
					typename T,
					typename m2lhandler_type>
void ApplyM2Loperators
(const std::vector<std::vector<Cluster<DIM>*> >& rclvv,
 std::vector<ExpansionHandler<DIM,ORDER,T> >& rexph,
 const std::vector<ExpansionHandler<DIM,ORDER,T> >& cexph,
 const std::vector<Relations<DIM> >& relations,
 std::vector<m2lhandler_type>& m2lv)
{
  typedef M2Lapplyer<DIM,ORDER,T,m2lhandler_type> m2lapplyer_type;
  for (unsigned int l=0; l<rclvv.size(); ++l) {
    const std::vector<Cluster<DIM>*>& rclv = rclvv.at(l);
    m2lhandler_type& m2l = m2lv.at(l);

		if (m2l.hasInteractions()) {

			const double t = omp_get_wtime();
			std::cout << "  at " << l << std::flush;

			m2l.activate();

#ifdef DFMM_USE_OMP
			m2lapplyer_type apply(rexph, cexph, relations, m2l);
#pragma omp parallel for schedule(dynamic) // guided
			for (typename std::vector<Cluster<DIM>*>::const_iterator
						 it=rclv.begin(); it<rclv.end(); ++it) apply(*it);
#else
			std::for_each(rclv.begin(), rclv.end(),
										m2lapplyer_type(rexph, cexph, relations, m2l));
#endif
			
			m2l.deactivate();

			std::cout << " took " << omp_get_wtime() - t << " s" << std::endl;

		}
      
  }
}



#endif

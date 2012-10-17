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

#ifndef relationsetter_tpl
#define relationsetter_tpl


#include <boost/foreach.hpp>



/*!
  \ingroup mcluster

  \class RelationSetter
  
  \brief Abstract base class for relation-setter classes
*/
template <typename cluster_type,
          typename level_type,
          typename expansionhandler_type,
          typename relations_type,
          typename m2lhandler_type>
class RelationSetter
  : public std::unary_function<cluster_type, void>
{
  enum {dim    = cluster_type::dim,
        nboxes = cluster_type::nboxes,
        nnodes = expansionhandler_type::nnodes};

  typedef typename cluster_type::point_type point_type;

  typedef typename expansionhandler_type::value_type value_type;
  typedef typename expansionhandler_type::expansion_map      expansion_map;
  typedef typename expansionhandler_type::expansion_map_flag expansion_map_flag;

  std::vector<expansionhandler_type>& rexph;
  std::vector<expansionhandler_type>& cexph;

  const std::vector<level_type*>& levels;

  m2lhandler_type& m2l;
  

  value_type *const alloc(const cluster_type *const cl,
                          const unsigned int c,
                          const point_type& u,
                          std::vector<expansionhandler_type>& exph)
  {
    // allocate expansion
    point_type uu(0.);
    level_type *const level = levels.at(cl->getNlevel());
    if (level->is_in_high_frequency_regime()) {
      uu = level->getConeDirection(c);
      level->insert(c); // cone c exists, insert it
    }
    expansion_map_flag expansion = exph.at(cl->getCidx()).insert(c, uu);

    // if direction not yet inserted and cluster is no leaf and cluster is in
    // high frequency regime insert cluster into child clusters recursively
    if(expansion.second && !cl->isLeaf())
      for (unsigned int b=0; b<nboxes; ++b) {
        const cluster_type *const ccl = cl->getChild(b);
        assert(!expansion.first->second->getcvals(b));
        if (ccl) {
          unsigned int cc = 0;
          const level_type *const clevel = levels.at(ccl->getNlevel());
          if (clevel->is_in_high_frequency_regime())
            cc = clevel->getConeIdx(u);
          expansion.first->second->setcvals(b, alloc(ccl, cc, u, exph));
        }
      }
    
    return expansion.first->second->getvals();
  }
 


protected:
  std::vector<relations_type>& relations;


  explicit RelationSetter(std::vector<expansionhandler_type>& _rexph,
                          std::vector<expansionhandler_type>& _cexph,
                          std::vector<relations_type>& _relations,
                          const std::vector<level_type*>& _levels,
                          m2lhandler_type& _m2l)
    : rexph(_rexph),
      cexph(_cexph),
      levels(_levels),
      m2l(_m2l),
      relations(_relations)
  {}

  virtual void operator()(const cluster_type *const cl) = 0;

  void setRelation(const cluster_type *const clx, 
                   const cluster_type *const cly)
  {
    typedef boost::tuple<bool,unsigned int,point_type> admissibility_info;
    
    assert(clx && clx->getNlevel() == m2l.getLevel()->getNlevel());
    assert(cly && cly->getNlevel() == m2l.getLevel()->getNlevel());
    
    // relation info
    const point_type tvec(cly->getCenter() - clx->getCenter());
    const admissibility_info info = m2l.getLevel()->getRelation(tvec);
    const bool admissible = boost::get<0>(info);
    
    // get relations
    relations_type& rl = relations.at(clx->getCidx());

    // if admissible - interaction
    if (admissible) {
      const unsigned int c = boost::get<1>(info);
      // insert unique transfer vector
      const unsigned int t_id = m2l.addTransferVector(clx, cly, c);
      // insert interaction
      rl.push_back_interaction(cly, t_id, c);
      // recursively insert directional expansions
      alloc(clx, c, boost::get<2>(info), rexph);
      alloc(cly, c, boost::get<2>(info), cexph);
    }
    // else insert neighbor
    else rl.push_back_neighbor(cly);
  }

};






/*!
  \ingroup mcluster

  \class RelationSetterCousins
  
  \brief Checks relation between cousins (parents neighbors)
*/
template <typename cluster_type,
          typename level_type,
          typename expansionhandler_type,
          typename relations_type,
          typename m2lhandler_type>
class RelationSetterCousins
  : RelationSetter<cluster_type,level_type,expansionhandler_type,
                   relations_type,m2lhandler_type>
{
  using RelationSetter<cluster_type,level_type,expansionhandler_type,
                       relations_type,m2lhandler_type>::relations;

public:
  RelationSetterCousins(std::vector<expansionhandler_type>& rexph,
                        std::vector<expansionhandler_type>& cexph,
                        std::vector<relations_type>& relations,
                        const std::vector<level_type*>& levels,
                        m2lhandler_type& m2l)
    : RelationSetter<cluster_type,level_type,expansionhandler_type,
                     relations_type,m2lhandler_type>
      (rexph, cexph, relations, levels, m2l)
  {}

  void operator()(const cluster_type *const clx)
  {
    // all neighbors of the parent cluster are either neighbors of
    // interactions
    BOOST_FOREACH(const cluster_type *const clY,
                  relations.at(clx->getParent()->getCidx()).getNeighbors()) {
      for (unsigned int b=0; b<cluster_type::nboxes; ++b) {
        const cluster_type *const cly = clY->getChild(b);
        if (cly) setRelation(clx, cly);
      }
    }
  }
};



/*!
  \ingroup mcluster

  \class RelationSetterPeers
  
  \brief Checks relation between peers (clusters of same level)

  \warning Not tested for a long time
*/
template <typename cluster_type,
          typename level_type,
          typename expansionhandler_type,
          typename relations_type,
          typename m2lhandler_type>
class RelationSetterPeers
  : RelationSetter<cluster_type,level_type,expansionhandler_type,
                   relations_type,m2lhandler_type>
{
  std::vector<cluster_type*>& pr;
  
public:
  explicit RelationSetterPeers(std::vector<expansionhandler_type>& rexph,
                               std::vector<expansionhandler_type>& cexph,
                               std::vector<relations_type>& relations,
                               const std::vector<level_type*>& levels,
                               m2lhandler_type& m2l,
                               const std::vector<cluster_type*>& _pr)
    : RelationSetter<cluster_type,level_type,expansionhandler_type,
                     relations_type,m2lhandler_type>
      (rexph, cexph, relations, levels, m2l), pr(_pr)
  {}

  void operator()(const cluster_type *const clx)
  {
    // all siblings are either neighbors or interactions
    BOOST_FOREACH(const cluster_type *const cly, pr)
      setRelation(clx, cly);
  }
};








//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*!
  \ingroup mcluster
*/
template <typename cluster_type,
          typename level_type,
          typename expansionhandler_type,
          typename relations_type,
          typename m2lhandler_type>
void SetupClusterRelations
(const std::vector<std::vector<cluster_type*> >& rclvv,
 std::vector<expansionhandler_type>& rexph,
 const std::vector<std::vector<cluster_type*> >& cclvv,
 std::vector<expansionhandler_type>& cexph,
 std::vector<relations_type>& relations,
 const std::vector<level_type*>& levels,
 std::vector<m2lhandler_type>& m2lv)
{
  // setup multilevel relations
  typedef RelationSetterCousins<cluster_type,level_type,expansionhandler_type,
                                relations_type,m2lhandler_type>
    cousin_relation_setter;

  assert(rclvv.front().size() == 1);
  assert(cclvv.front().size() == 1);
  const cluster_type *const rcl_root = *(rclvv.front().begin());
  const cluster_type *const ccl_root = *(cclvv.front().begin());
  relations.at(rcl_root->getCidx()).push_back_neighbor(ccl_root);
  
  for (unsigned int l=1; l<rclvv.size(); ++l)
    std::for_each(rclvv.at(l).begin(), rclvv.at(l).end(),
                  cousin_relation_setter(rexph, cexph,
                                         relations, levels, m2lv.at(l)));
}







#endif

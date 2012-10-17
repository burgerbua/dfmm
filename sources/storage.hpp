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

#ifndef storage_hpp
#define storage_hpp

#include <map>
#include <vector>


#include "sources/level.hpp"


/*! \class Storage
  The class storage is thought as a kind of singleton; it initializes and
  stores all levels of the octree.
 */
template <int DIM>
class Storage
  : boost::noncopyable
{
  enum {dim    = DIM,
        nboxes = DimTraits<dim>::nboxes};

  typedef typename DimTraits<dim>::point_type point_type;

  typedef Level<dim>               level_type;
  typedef std::vector<level_type*> level_vector;

  level_vector levels;

  double     root_cluster_diam;
  point_type root_cluster_center;




 public:
  explicit Storage(const std::pair<point_type,point_type> bbx,
                   const unsigned int lmax);
  

  ~Storage()
  {
    for (typename level_vector::iterator it=levels.begin();
         it!=levels.end(); ++it) delete *it;
    levels.clear();
  }

  const unsigned int getNumLevels() const
  {
    assert(!levels.empty());
    return levels.size();
  }

  const level_vector& getLevels() const
  {
    return levels;
  }

  const level_type *const getLevel(const unsigned int l) const
  {
    assert(!levels.empty());
    return levels.at(l);
  }

  level_type *const getLevel(const unsigned int l)
  {
    assert(!levels.empty());
    return levels.at(l);
  }

  
  //! \name Set either levels for a dynamic or static analysis
  /* @{ */
  void initLevels()
  {
    for (unsigned int l=0; l<levels.size(); ++l) {
      level_type *const level = levels.at(l);
      level->setLowFrequencyRegime();
    }
  }

  void initLevels(const double wavenum)
  {
    const double delimiter = FREQUENCY_REGIME_DELIMITER;
    for (unsigned int l=0; l<levels.size(); ++l) {
      level_type *const level = levels.at(l);
      const double diam = level->getDiam();
      if (diam*wavenum > delimiter) level->setHighFrequencyRegime(wavenum);
      else                          level->setLowFrequencyRegime();
    }
  }
  /* @} */
  
  

  const double getRootClusterDiam()
  {
    return root_cluster_diam;
  }
  
  const point_type getRootClusterCenter()
  {
    return root_cluster_center;
  }




  template <typename m2l_handler_type>
  void initialize(std::vector<m2l_handler_type>& operators) const
  {
    for (typename level_vector::const_iterator it=levels.begin();
         it!=levels.end(); ++it)
      operators.push_back(m2l_handler_type(*it));
  }
  



  template <typename m2l_handler_type>
  void writeInfo(const std::vector<m2l_handler_type>& handler) const
  {
    assert(handler.size() == levels.size());

    // root cluster info
    std::cout << "\nBounding box of diam " << root_cluster_diam
              << " is centered at " << root_cluster_center << std::endl;

    // level info
    std::cout << "- Level info:" << std::endl;
    for (unsigned int l=0; l<levels.size(); ++l) {
      levels.at(l)->writeInfo();
      handler.at(l).writeInfo();
    }
  }
  
};





#include "sources/storage.tpl"


#endif

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
  
  \brief Manages levels of quad/octree and initializes M2L operators
  
  The class storage is thought as a kind of singleton; it initializes and
  stores all levels of the quad/octree and initializes all M2L operators.
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

  /*! The constructor sets root cluster center and width and creates all
      quad/octree levels
      
      @param[in] bbx bounding box of root cluster
      @param[in] lmax maximum number of levels
  */
  explicit Storage(const std::pair<point_type,point_type> bbx,
                   const unsigned int lmax);
  
  
  ~Storage()
  {
    for (typename level_vector::iterator it=levels.begin();
         it!=levels.end(); ++it) delete *it;
    levels.clear();
  }

  
  /*! @return number of quad\octree levels */
  const unsigned int getNumLevels() const
  {
    assert(!levels.empty());
    return levels.size();
  }


  /*! @return reference to vector containing all levels */
  const level_vector& getLevels() const
  {
    return levels;
  }

  
  /*! @param[in] l level index
    @return const pointer to level l
  */
  const level_type *const getLevel(const unsigned int l) const
  {
    assert(!levels.empty());
    return levels.at(l);
  }

  /*! @param[in] l level index
    @return pointer to level l
  */
  level_type *const getLevel(const unsigned int l)
  {
    assert(!levels.empty());
    return levels.at(l);
  }
  

  /*! Set levels for static analysis */
  void initLevels()
  {
    for (unsigned int l=0; l<levels.size(); ++l) {
      level_type *const level = levels.at(l);
      level->setLowFrequencyRegime();
    }
  }


  /*!  Set levels for dynamic analysis

    \note Based on experience and numerical studies we have found that the
    threshold between low- and high-frequency regime is given by \f$wk >
    0.8\f$. The value for \f$0.8\f$ can be chosen differently if desired (it
    is not recommended, though).
      
    @param[in] wavenum wavenumber
  */
  void initLevels(const double wavenum)
  {
    for (unsigned int l=0; l<levels.size(); ++l) {
      level_type *const level = levels.at(l);
      const double diam = level->getDiam();
      if (diam*wavenum > 0.8) level->setHighFrequencyRegime(wavenum);
      else                    level->setLowFrequencyRegime();
    }
  }
  
  
  /*! @return root cluster width */
  const double getRootClusterDiam()
  {
    return root_cluster_diam;
  }
  

  /*! @return root cluster center */
  const point_type getRootClusterCenter()
  {
    return root_cluster_center;
  }


  /*! Instantiate an M2L operator per level
    
    @param[out] operators vector containing M2L operators
  */
  template <typename m2l_handler_type>
  void initialize(std::vector<m2l_handler_type>& operators) const
  {
    operators.clear();
    for (typename level_vector::const_iterator it=levels.begin();
         it!=levels.end(); ++it)
      operators.push_back(m2l_handler_type(*it));
  }
  

  /*! Write out information of all M2L operators and all levels of the
      quad/octree
      
      @param[in] operators vector containing M2L operators
  */
  template <typename m2l_handler_type>
  void writeInfo(const std::vector<m2l_handler_type>& operators) const
  {
    assert(operators.size() == levels.size());

    // root cluster info
    std::cout << "\nBounding box of diam " << root_cluster_diam
              << " is centered at " << root_cluster_center << std::endl;

    // level info
    std::cout << "- Level info:" << std::endl;
    for (unsigned int l=0; l<levels.size(); ++l) {
      levels.at(l)->writeInfo();
      operators.at(l).writeInfo();
    }
  }
  
};





#include "sources/storage.tpl"


#endif

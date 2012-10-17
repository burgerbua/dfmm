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

#ifndef level_hpp
#define level_hpp


#include "sources/frequencyregime.hpp"








/**
 * \ingroup mlevel
 *
 * \class Level
 *
 * This class contains the information which is constant along a level in the
 * quad/oct-tree. Mainly it handles the low- and high-frequency regime status
 * of a given level. It also stores the diameter of all clusters in a given
 * level (we deal only with uniform trees). There exists only one instance of
 * a Level container per cluster level.
 *
 */

template <int DIM>
class Level
{
  // typedefs
  typedef typename DimTraits<DIM>::point_type point_type;
  typedef typename DimTraits<DIM>::point_int_type point_int_type;
  typedef typename Frequency_Regime<DIM>::admissibility_info admissibility_info;

  
  // private member variables  
  Frequency_Regime<DIM>* regime; //<! frequency regime (high or low)
  const unsigned int nlevel;     //!< level index
  const bool is_leaf;            //!< is leaf level or not
  const double diam;             //!< cluster diameter
  
  
public:
  explicit Level(const unsigned int _nlevel,
                 const double _diam,
                 const bool _is_leaf)
    : regime(new Frequency_Regime<DIM>()),
      nlevel(_nlevel),
      is_leaf(_is_leaf),
      diam(_diam)
  {}
  
  
  
  ~Level()
  {
    if (regime) delete regime;
    regime = NULL;
  }
  
  
  
  /**
   * This functions convert a transfer vector @p dvec with floating point
   * coordinates into a transfer vector @p ivec with integral coordinates and
   * vice versa.
   */
  /* @{ */
  const point_int_type convert(const point_type& dvec) const
  {
    point_int_type ivec;
    for (unsigned int d=0; d<DIM; ++d)
      ivec[d] = std::floor(dvec[d]/diam + 0.5);
    return ivec;
  }
  /* @} */
  
  
  
  
  /*! @return level index */
  const unsigned int getNlevel() const
  { return nlevel; }
  
  /*! @return if cluster is leaf or not */
  const bool isLeaf() const
  { return is_leaf; }
  
  /*! @return diameter of all clusters at current level */
  const double getDiam() const
  { return diam;  }
  
  
  //! \name Control functions
  /* @{ */
  const bool is_frequency_regime_set() const
  { return regime->is_frequency_regime_set(); }
  const bool is_in_high_frequency_regime() const
  { return regime->is_in_high_frequency_regime(); }
  /* @} */
  
  
  
  //! \name Set regimes of level and all clusters
  /* @{ */
  void setHighFrequencyRegime(const double wavenum)
  {
    assert(!regime->is_frequency_regime_set());
    delete regime;
    regime = new High_Frequency_Regime<DIM>(diam, wavenum);
  } 
  
  void setLowFrequencyRegime()
  {
    assert(!regime->is_frequency_regime_set());
    delete regime;
    regime = new Low_Frequency_Regime<DIM>(diam);
  }
  /* @} */
  


  //! @return unique id for a given relative position vector
  const unsigned int getTransferId(const point_int_type& relpos) const
  {
    const unsigned int size = 2 << (nlevel-1);
    unsigned int tid = 0;
    unsigned int scale = 1;
    for (unsigned int d=0; d<DIM; ++d) {
      tid += (relpos[d] + size) * scale;
      scale *= size;
    }
    return tid;
  }

  

  
  //! set relation to other cluster
  const admissibility_info getRelation(const point_type& transfer) const
  {
    return regime->getRelation(transfer);
  }
  
  



  
  //! \name High frequency regime only
  /* @{ */
  const unsigned int getConeIdx(const point_type& direction) const
  {
    return regime->getConeIdx(direction);
  }
  const unsigned int getNCones() const
  {
    return regime->getNCones();
  }
  const point_type getConeDirection(const unsigned int c) const
  {
    return regime->getConeDirection(c);
  }
  void insert(const unsigned int c)
  {
    regime->insert(c);
  }
  const std::map<unsigned int,point_type>& getExistDirs() const
  {
    return regime->getExistDirs();
  }
  /* @} */
  





  //! Write infos
  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    out << "  " << nlevel
        << (is_in_high_frequency_regime() ? " - HF " : " - LF ")
        << "[d,mdist]:[" << diam << ",";
    regime->writeInfo(out);
    return out;
  }


};






#endif

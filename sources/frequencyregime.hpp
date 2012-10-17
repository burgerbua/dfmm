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

#ifndef frequencyregime_hpp
#define frequencyregime_hpp

#include <stdexcept>
#include <iomanip>
#include <set>

#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/static_assert.hpp>

#include "./dimordertraits.hpp"


/*!
  \defgroup mfreqregime Frequency Regimes
 
  There exist two frequency regimes: the low- and the high-frequency
  regime. In our case the low-frequency regime correspond to where we all
  kernels (P2M, M2M, L2L and L2P) are non-directional. Hence, if the kernel
  function is a non-oscillatory function (eg. the Laplace kernel) all levels
  are in the low-frequency (non-directional) regime. On the other hand in the
  high-frequency regime all kernels are directional and is needed when
  oscillatory kernels are treated (eg. Helmholtz kernel). For the treshold
  between low- and high-frequency regime we refer to the paper mentioned on
  the <a href="index.html">main page</a>.
 */


/*!
  \ingroup mfreqregime
 
  \class Frequency_Regime
 
  \brief Abstract base class for low- and high frequency regime
 
  This class is an (not completely true) abstract base class, in the sense
  that it throws exeptions when here declared functions are not defined in
  either the high or the low frequency regime. We used this workaround since
  any tree level \b must either be in the low- or in the high-frequency
  regime.

  @tparam DIM spatial dimension
 */
template <int DIM>
class Frequency_Regime
{
protected:
  typedef typename DimTraits<DIM>::point_type point_type;


public:
  typedef boost::tuple<bool,unsigned int,point_type> admissibility_info;


  virtual ~Frequency_Regime() {}


  // set relation to other cluster
  virtual const admissibility_info getRelation(const point_type& t) const
  { throw std::runtime_error("Frequency regime not yet set!"); }

  /**
   * \name Control functions
   */
  /* @{ */
  virtual const bool is_frequency_regime_set() const
  { return false; }
  virtual const bool is_in_high_frequency_regime() const
  { throw std::runtime_error("is_high_freq: frequency regime not yet set!"); }
  /* @} */
  

  //! \name High frequency regime only
  /* @{ */
  virtual const unsigned int getConeIdx(const point_type& direction) const
  { throw std::runtime_error("getConeIdx: not high frequency!"); }
  virtual const unsigned int getNCones() const
  { throw std::runtime_error("getNCones: not high frequency!"); }
  virtual const point_type getConeDirection(const unsigned int c) const
  { throw std::runtime_error("getConeDirection: not high frequency!"); }
  virtual void insert(const unsigned int c)
  { throw std::runtime_error("insert: not high frequency!"); }
  virtual const std::map<unsigned int,point_type>& getExistDirs() const
  { throw std::runtime_error("getExistDirs: not high frequency!"); }
  /* @} */


  virtual std::ostream& writeInfo(std::ostream &out) const
  {
    out << "Frequency regime not yet set." << std::endl;
    return out;
  }

};


















/*!
  \ingroup mfreqregime
 
  \class Low_Frequency_Regime
 
  \brief Low-frequency functionalities for Level
 
  @tparam DIM spatial dimension
 */
template <int DIM>
class Low_Frequency_Regime
  : public Frequency_Regime<DIM>,
    private boost::noncopyable
{
  typedef typename Frequency_Regime<DIM>::point_type point_type;
  typedef typename Frequency_Regime<DIM>::admissibility_info admissibility_info;

  double mindist; //!< minimal admissible distance between two clusters



public:
  Low_Frequency_Regime(const double extension)
    : Frequency_Regime<DIM>()
  {
#if defined   LF_MIN_DISTANCE_1
    mindist = extension * 1.1 * sqrt(3);
#elif defined LF_MIN_DISTANCE_2
    mindist = 2. * extension;
#elif defined LF_MIN_DISTANCE_3
    mindist = 2. * extension * sqrt(2);
#else
    throw std::runtime_error("Frequency regime delimiter not defined!");
#endif
  }

  
  
  //! \name Control functions
  /* @{ */
  const bool is_frequency_regime_set() const
  {
    return true;
  }
  
  const bool is_in_high_frequency_regime() const
  {
    return false;
  }
  /* @} */


  const admissibility_info getRelation(const point_type& transfer) const
  {
    const double dist = transfer.norm();
    const unsigned int c = 0;
    const point_type u(0.);
    if (dist>mindist) return boost::make_tuple(true,  c, u);
    else              return boost::make_tuple(false, c, u);
  }



  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    out << mindist << "]," << std::setw(20) << "\t[ndirs,nilists,ni]:[" << 1;
    return out;
  }

};















#include "sources/conefactory.hpp"



/*!
  \ingroup mfreqregime
 
  \class High_Frequency_Regime
 
  \brief High-frequency functionalities for Level

  @tparam DIM spatial dimension
 
 */
template <int DIM>
class High_Frequency_Regime
  : public Frequency_Regime<DIM>,
    private boost::noncopyable
{
  // typedef
  typedef ConeFactory<DIM> cone_factory_type;
  typedef typename Frequency_Regime<DIM>::point_type point_type;
  typedef typename Frequency_Regime<DIM>::admissibility_info admissibility_info;

  

  /**
   * The ConeFactory contains all cone informations for a given cluster level
   * and a give wavenumber.
   */
  const cone_factory_type cf;

  double mindist; //!< minimal admissible distance between two clusters
 
  std::map<unsigned int,point_type> existDirs;
  

public:
  High_Frequency_Regime(const double extension,
                        const double wavenum)
    : Frequency_Regime<DIM>(),
      cf(extension, wavenum)
  {
#if defined   HF_MIN_DISTANCE_1
    const double mdlf = extension * 1.1 * sqrt(3);
#elif defined   HF_MIN_DISTANCE_2
    const double mdlf = 2. * extension;
#elif defined   HF_MIN_DISTANCE_3
    const double mdlf = 2. * extension * sqrt(2);
#else
    throw std::runtime_error("Frequency regime delimiter not defined!");
#endif
    const double mdhf = extension*extension * wavenum;
    mindist = ((mdhf>mdlf) ? mdhf : mdlf);
  }
  


  //! \name Control functions
  /* @{ */
  const bool is_frequency_regime_set() const
  {
    return true;
  }
  
  const bool is_in_high_frequency_regime() const
  {
    return true;
  }
  /* @} */



  //! \name Informations about cones
  /* @{ */
  const unsigned int getNCones() const
  {
    return cf.getNCones();
  }
  const point_type getConeDirection(const unsigned int c) const
  {
    assert(c==cf.getConeIdx(cf.getConeDirection(c)));
    return cf.getConeDirection(c);
  }
  const unsigned int getConeIdx(const point_type& direction) const
  {
    return cf.getConeIdx(direction);
  }
  /* @} */



  void insert(const unsigned int c)
  { 
    existDirs.insert(std::make_pair(c, cf.getConeDirection(c)));
  }
  
  const std::map<unsigned int,point_type>& getExistDirs() const
  {
    return existDirs;
  }


  const admissibility_info getRelation(const point_type& transfer) const
  {
    const double dist = transfer.norm();
    const unsigned int c = cf.getConeIdx(transfer);
    const point_type u = cf.getConeDirection(c); 
    if (dist>mindist) return boost::make_tuple(true,  c, u);
    else              return boost::make_tuple(false, c, u);
  }

  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    out << mindist << "]," << std::setw(30) << "\t[ndirs,nilists,ni]:["
        << cf.getNCones();
    return out;
  }

};




#endif

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

#ifndef mapping_hpp
#define mapping_hpp


#include <boost/utility.hpp> 



/**
 * Base class for cluster mapper
 */
template <int DIM>
class Mapping
  : boost::noncopyable
{
  typedef typename DimTraits<DIM>::point_type point_type;

private:
  void setBoundingBox(const point_type& center,
                      const double diam)
  {
    for (unsigned int d=0; d<DIM; ++d) {
      a[d] = center[d] - diam / 2;
      b[d] = center[d] + diam / 2;
    }
  }


protected:
  point_type a;
  point_type b;
 
  explicit Mapping(const point_type& center,
                   const double diam)
  { setBoundingBox(center, diam); }

  explicit Mapping()
  { for (unsigned int d=0; d<DIM; ++d) {a[d]=-1.; b[d]=1.;} }


  virtual const double operator()(const int d, const double x) const = 0;

public:
  template <typename POINT>
  bool is_in_cluster(const POINT& point) const
  {
    for (unsigned int d=0; d<DIM; ++d)
      if (a[d]-point[d] > 10. * std::numeric_limits<double>::epsilon() ||
          point[d]-b[d] > 10. * std::numeric_limits<double>::epsilon()) {
        std::cout << a[d]-point[d] << "\t" << point[d]-b[d]
                  << "\t" << std::numeric_limits<double>::epsilon()
                  << std::endl;
        return false;
      }
    return true;
  }
};


/**
 * Global to local cluster mapper
 */
template <int DIM>
class map_glob_loc
  : public Mapping<DIM>
{
  typedef typename DimTraits<DIM>::point_type point_type;

  using Mapping<DIM>::a;
  using Mapping<DIM>::b;

public:
  explicit map_glob_loc()
    : Mapping<DIM>() {}
  explicit map_glob_loc(const point_type& center,
                        const double diam)
    : Mapping<DIM>(center, diam) {}
  
  const double operator()(const int d, const double x) const
  { return (2*x - b[d] - a[d]) / (b[d] - a[d]); }

  const double getJacobian(const int d) const
  { return (2. / (b[d] - a[d])); }
};


/**
 * Local to global cluster mapper
 */
template <int DIM>
class map_loc_glob
  : public Mapping<DIM>
{
  typedef typename DimTraits<DIM>::point_type point_type;

  using Mapping<DIM>::a;
  using Mapping<DIM>::b;

public:
  explicit map_loc_glob()
    : Mapping<DIM>() {}
  explicit map_loc_glob(const point_type& center,
                        const double diam)
    : Mapping<DIM>(center, diam) {}
  
  const double operator()(const int d, const double x) const
  { return (a[d] + b[d]) / 2. + (b[d] - a[d]) * x / 2.; }
};



#endif

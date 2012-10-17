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

#ifndef expansion_hpp
#define expansion_hpp





/*!
  \ingroup mcluster

  \class Expansion
  
  \brief Stores multipole and local expansions
*/
template <int DIM, int ORDER, typename T>
class Expansion
{

  static const unsigned int FAC = 2;

  // typedefs and compile time constants
  enum {nboxes = DimTraits<DIM>::nboxes,
        nnodes = BasisTraits<ORDER,DIM>::nnodes};
  typedef T value_type;
  typedef typename DimTraits<DIM>::point_type point_type;

  // members
  const point_type u;
  value_type *const vals;
  boost::array<value_type*,nboxes> cvals;

  
public:
  explicit Expansion(const point_type& _u)
    : u(_u), vals(new value_type [FAC*nnodes])
  {
    this->zero();
    for (unsigned int b=0; b<nboxes; ++b) cvals.at(b) = NULL;
  }
  
  ~Expansion()
  { delete [] vals; }

  void zero()
  { blas::setzero(FAC*nnodes, vals); }

  void setcvals(const unsigned int b, value_type *const cval)
  { cvals.at(b) = cval; }

  const point_type& getu() const
  { return u; }
  const value_type *const getvals() const
  { return vals; }
  value_type *const getvals()
  { return vals; }
  const value_type *const getcvals(const unsigned int b) const
  { return cvals.at(b); }
  value_type *const getcvals(const unsigned int b)
  { return cvals.at(b); }

};







/*!
  \ingroup mcluster

  \class ExpansionHandler
  
  \brief Handles multipole and local expansions of a cluster for all directions
*/
template <int DIM, int ORDER, typename T>
class ExpansionHandler
{
  typedef typename DimTraits<DIM>::point_type point_type;

public:
  enum {dim    = DIM,
        order  = ORDER,
        nnodes = BasisTraits<ORDER,DIM>::nnodes};

  typedef T value_type;

  typedef Expansion<DIM,ORDER,T>                      expansion_type;
  typedef std::map<unsigned int,expansion_type*const> expansion_map;
  typedef typename expansion_map::iterator            expansion_map_iterator;
  typedef std::pair<expansion_map_iterator,bool>      expansion_map_flag;


private:
  expansion_map expansions;


public:


  explicit ExpansionHandler()
  {
    expansions.clear();
  }

  
  ~ExpansionHandler()
  {
    BOOST_FOREACH(typename expansion_map::value_type expansion_pair,
                  expansions) delete expansion_pair.second;
    expansions.clear();
  }

  
  std::pair<typename expansion_map::iterator,bool>
  insert(const unsigned int c, const point_type& u)
  {
    typename expansion_map::iterator exp = expansions.find(c);
    if (exp!=expansions.end()) return std::make_pair(exp, false);
    else return expansions.insert(std::make_pair(c, new expansion_type(u)));
  }
  
  
  expansion_map& getExpansions()
  { return expansions; }
  const expansion_map& getExpansions() const
  { return expansions; }


  void clearExpansions()
  {
    BOOST_FOREACH(typename expansion_map::value_type& ep, expansions)
      ep.second.zero();
  }

};

#endif

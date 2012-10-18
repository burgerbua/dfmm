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

#ifndef transferhandler_hpp
#define transferhandler_hpp





#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>




/*! \struct point_int_typeHash
  Unary function which computes the hash value of a point whose coordinates
  are integral values.
 */
template <int DIM>
struct point_int_typeHash
  : std::unary_function<typename DimTraits<DIM>::point_int_type, std::size_t>
{
  typedef typename DimTraits<DIM>::point_int_type point_int_type;

  std::size_t operator()(const point_int_type& p) const
  {
    std::size_t seed = 0;
    for (unsigned int d=0; d<DIM; ++d)
      boost::hash_combine(seed, p[d]);
    return seed;
  }
};


/*! \struct point_int_typeCmpr 
  Binary function which compares two points where the coordinates are integral
  numbers
 */
template <int DIM>
struct point_int_typeCmpr
  : std::binary_function<typename DimTraits<DIM>::point_int_type,
                         typename DimTraits<DIM>::point_int_type, bool>
{
  typedef typename DimTraits<DIM>::point_int_type point_int_type;
  
  bool operator()(const point_int_type& v1,
                  const point_int_type& v2) const
  {
    return (v1 == v2);
  }
};






/*! /class ExistHandler
 This class provides various maps to store unique (existing) transfer-vectors
  (points whose coordinates have integral values).  
 */
template <int DIM>
class ExistHandler
{
  typedef point_int_typeCmpr<DIM> Cmpr;
  typedef point_int_typeHash<DIM> Hash;
  
  typedef typename DimTraits<DIM>::point_int_type point_int_type;
  

protected:
  ExistHandler()
    : existing(false) 
  {
    exist_orig.clear();
    exist_perm.clear();
    exist_dirs_orig.clear();
  }

  typedef boost::unordered_map<point_int_type,
                               const unsigned int,
                               Hash,Cmpr>                  exist_map;
  typedef std::pair<typename exist_map::iterator,bool>     exist_map_key;
  typedef boost::unordered_map<unsigned int,exist_map>     exist_map_map;
  typedef std::pair<typename exist_map_map::iterator,bool> exist_map_map_key;

  bool existing; /*!< true if interaction list is not empty */
  exist_map exist_orig; /*!< stores pairs of original transfer-vector and its
                           incremental id */
  exist_map exist_perm;  /*!< stores pairs of permuted transfer-vector and its
                           incremental id */
  exist_map_map exist_dirs_orig; /*!< needed for dFMM only, stores pairs cone
                                    indices and the corresponding directional
                                    interactions lists*/
};









//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*! \addtogroup mm2l
  
  \par Transfer-handler
  The transfer-handler are responsible for collecting and storing unique (we
  assume the kernel functions are translation invariant) transfer vectors
  which are then used to compute by \ref mm2lcomputer to compute the M2L
  operators. In all transfer handler we use counter which provide unique
  indices for each transfer vector and, hence, allow us to use vectors instead
  of maps (faster access).
 */


/*! \ingroup mm2l

  \class TransferHandler

  \brief Abstract base class for all transfer handler
 
  @tparam DIM spatial dimension
 */
template <int DIM>
class TransferHandler
{
protected:
  virtual const bool hasInteractions() const = 0;
  virtual std::ostream& writeInfo(std::ostream&) const = 0;
  virtual const unsigned int addTransferVector(const Cluster<DIM> *const,
                                               const Cluster<DIM> *const,
                                               const unsigned int) = 0;
};














/*! \ingroup mm2l

  \class TransferHandlerSA

  \brief Handles unique transfer vectors for all cones
  
  This transfer handler is needed for the \f$K_i \sim UC_iB_i^*\f$
  optimization of the M2L operators. It requires the unified treatment of an
  entire set of transfer vectors, hence, here they are collected in a
  cone-wise manner.

  @tparam DIM spatial dimension
 */
template <int DIM>
class TransferHandlerSA
  : public TransferHandler<DIM>,
    public ExistHandler<DIM>
{
private:
  typedef typename DimTraits<DIM>::point_type     point_type;
  typedef typename DimTraits<DIM>::point_int_type point_int_type;

  typedef std::vector<point_int_type>            trans_vec;
  typedef std::map<unsigned int,trans_vec>       trans_vec_map;
  typedef typename trans_vec_map::value_type     trans_vec_map_pair;
  typedef typename trans_vec_map::iterator       trans_vec_map_iterator;
  typedef std::pair<trans_vec_map_iterator,bool> trans_vec_map_key;

  // used members
  using ExistHandler<DIM>::exist_dirs_orig;
  using ExistHandler<DIM>::existing;


protected:

  // own members  
  trans_vec_map ffield; //!< contains union of all possible interactions


public:

  /// Destructor
  virtual ~TransferHandlerSA() {}

  
  const unsigned int addTransferVector(const Cluster<DIM> *const clx,
                                       const Cluster<DIM> *const cly,
                                       const unsigned int c)
  {
    // interactions are existing
    existing = true;

    // typedefs
    typedef typename ExistHandler<DIM>::exist_map     exist_map;
    typedef typename ExistHandler<DIM>::exist_map_key exist_map_key;

    // insert transfer list c
    typename ExistHandler<DIM>::exist_map_map_key em
      = exist_dirs_orig.insert(std::make_pair(c,exist_map()));
    trans_vec_map_key tv = ffield.insert(trans_vec_map_pair(c,trans_vec()));
    assert(em.second == tv.second);

    // insert transfer vector t
    exist_map& e_map = em.first->second;
    const point_int_type t = cly->getPosition() - clx->getPosition();
    exist_map_key e = e_map.insert(std::make_pair(t, e_map.size()));

    // unique transfer vector id
    const unsigned int t_id = e.first->second;

    // push back transfer vector t if not yet existing
    if (e.second) {
      trans_vec& t_vec = tv.first->second;
      assert(t_vec.size() == t_id);
      t_vec.push_back(t);
    }

    // return transfer vector id
    return t_id;
  }



  const bool hasInteractions() const
  { return existing; }



  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    unsigned int counter = 0;
    for (typename trans_vec_map::const_iterator it=ffield.begin();
         it!=ffield.end(); ++it) counter += it->second.size();
    
    out << "," << ffield.size() << "," << counter << "]" << std::endl;
    return out;
  }
  
};









/*! \ingroup mm2l

  \class TransferHandlerNoSym

  \brief Handles unique transfer vectors is symmetries are take into account
 
  This transfer handler is used if symmetries in the M2L operators are \b not
  considered.
  
  @tparam DIM spatial dimension
 */
template <int DIM>
class TransferHandlerNoSym
  : public TransferHandler<DIM>,
    public ExistHandler<DIM>
{
private:
  using ExistHandler<DIM>::exist_orig;
  using ExistHandler<DIM>::existing;

protected:
  typedef typename DimTraits<DIM>::point_int_type point_int_type;
  typedef std::vector<point_int_type> transfer_vec;

  transfer_vec ffield;

public:
  const unsigned int addTransferVector(const Cluster<DIM> *const clx,
                                       const Cluster<DIM> *const cly,
                                       const unsigned int /* not needed */)
  {
    // interactions are existing
    existing = true;

    // level number needed to compute the unique id
    assert(clx->getNlevel() == cly->getNlevel());
    const unsigned int nlevel = clx->getNlevel();

    // insert transfer vector t
    typedef typename ExistHandler<DIM>::exist_map exist_map;
    typedef std::pair<typename exist_map::iterator,bool> exist_map_key;
    const point_int_type t_orig = cly->getPosition() - clx->getPosition();
    exist_map_key ek = exist_orig.insert(std::make_pair(t_orig,
                                                        exist_orig.size()));

    // unique transfer vector id
    const unsigned int id_orig = ek.first->second;

    // push back transfer vector t if not yet existing
    if (ek.second) {
      assert(ffield.size() == id_orig);
      ffield.push_back(t_orig);
    }

    // return transfer vector id
    return id_orig;
  }



  const bool hasInteractions() const
  { return existing; }



  void clearFfield()
  { 
    ffield.clear();
    exist_orig.clear();
  }


  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    out << ",1," << ffield.size() << "]" << std::endl;
    return out;
  }

};



















#include "sources/symmetries.hpp"


/*! \ingroup mm2l

  \class TransferHandlerSym

  \brief Handles unique transfer vectors if symmetries are considered
 
  This transfer handler is used if symmetries in the M2L operators are
  considered.
  
  @tparam ORDER interpolation order
  @tparam DIM spatial dimension
 */
template <int ORDER, int DIM>
class TransferHandlerSym
  : public TransferHandler<DIM>,
    public ExistHandler<DIM>
{
private:
  using ExistHandler<DIM>::exist_orig;
  using ExistHandler<DIM>::exist_perm;
  using ExistHandler<DIM>::existing;

protected:
  typedef typename DimTraits<DIM>::point_int_type point_int_type;
  typedef std::vector<point_int_type> transfer_vec;
  typedef std::pair<unsigned int,const unsigned int*> perm_pair;
  typedef std::vector<perm_pair> perm_vec;


  Symmetries<ORDER,DIM> symmetries; /*!< Provides the permutations and
                                      surjective mapping from \f$T\f$ to
                                      \f$T_{sym}\f$ */
  transfer_vec ffield; /*!< Stores transfer-vectors whose M2L operators need
                          to be computed. */
  perm_vec permut; /*!< Stores original original key of transfer-vector and
                      corresponding permutation vector */



public:

  /// Add unique transfer vectors
  const unsigned int addTransferVector(const Cluster<DIM> *const clx,
                                       const Cluster<DIM> *const cly,
                                       const unsigned int /* not needed */)
  {
    // interactions are existing
    existing = true;

    // typedefs
    typedef typename ExistHandler<DIM>::exist_map exist_map;
    typedef std::pair<typename exist_map::iterator,bool> exist_map_key;

    // level number needed to compute the unique id
    assert(clx->getNlevel() == cly->getNlevel());

    // original transfer vector and id
    const point_int_type t_orig = cly->getPosition() - clx->getPosition();
    exist_map_key ek_orig
      = exist_orig.insert(std::make_pair(t_orig, exist_orig.size()));
    const unsigned int id_orig = ek_orig.first->second;

    // permuted transfer vector, id and permutation array
    point_int_type t_perm;
    const unsigned int* pvec = symmetries.getPermutation(t_orig, t_perm);
    exist_map_key ek_perm
      = exist_perm.insert(std::make_pair(t_perm, exist_perm.size()));
    const unsigned int id_perm = ek_perm.first->second;

    // push back pair(id_perm, pvec) if not yet existing
    if (ek_orig.second) {
      assert(permut.size() == id_orig);
      permut.push_back(std::make_pair(id_perm, pvec));
    }
    
    // push back t_perm if not yet existing
    if (ek_perm.second) {
      assert(ffield.size() == id_perm);
      ffield.push_back(t_perm);
    }

    // return original transfer vector id
    return id_orig;
  }


  void clearFfield()
  {
    ffield.clear();
    exist_orig.clear();
    exist_perm.clear();
  }


  const bool hasInteractions() const
  { return existing; }



  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    out << ",1," << ffield.size() << "(sym)]" << std::endl;
    return out;
  }

};

#endif

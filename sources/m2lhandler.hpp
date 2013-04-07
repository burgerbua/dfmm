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

#ifndef m2lhandler_hpp
#define m2lhandler_hpp


#include "boost/tuple/tuple.hpp"


#include "sources/transferhandler.hpp"
#include "sources/m2lcomputer.hpp"
#include "sources/blas.h"




/*! 
  \addtogroup mm2l

  \par M2L-handler
  The M2L-handler is responsible for identifying, computing and storing all
  unique M2L operators.

*/



template <int> class map_loc_glob;
template <int> class Chebyshev;
template <int,int> class Tensor;
template <typename,typename> class EntryComputer;




/*! 
  \ingroup mm2l
  
  \brief Abstract base class for all m2l-handler classes
 */
template <int DIM, typename T>
class M2LHandler
{
protected:
  const Level<DIM>* level;

  mutable FlopCounter fcount;

  /// constructor
  M2LHandler(const Level<DIM> *const _level,
             const std::string& description)
    : level(_level), fcount(description, level->getNlevel())
  {}

  /// destructor
  virtual ~M2LHandler() {}

  /// Compressed M2L operation
  virtual void apply(const unsigned int,
                     const T *const, T *const, const unsigned int) const = 0;


public:
  /// Constructor
  const Level<DIM> *const getLevel() const { return level; }

  /// Apply part of **M2L** which can be shifted to **L2L**
  virtual void applyL2Lpart(T *const, const unsigned int) const {}

  /// Apply part of **M2L** which can be shifted to **M2M**
  virtual void applyM2Mpart(T *const, const unsigned int) const {}

  /// Needed for blocked m2l kernel application
  virtual void   activate() const {}
  virtual void deactivate() const {}
};










//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// 2 BIG SVD's PER CONE (SA) ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////




/*! 
  \ingroup mm2l

  \brief \f$K_i \sim UC_iB^*\f$ approximation per interaction list (ie. per
  cone/direction)
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerSArcmp
  : public M2LHandler<DIM,T>,
    public TransferHandlerSA<DIM>,
    public M2LComputerSA<DIM,ORDER,T>
{
  static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
  typedef typename M2LComputerSA<DIM,ORDER,T>::m2l_map::mapped_type data_tuple;

  // used members
  using M2LHandler<DIM,T>::level;
  using TransferHandlerSA<DIM>::ffield;
  using M2LComputerSA<DIM,ORDER,T>::operators;


public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;


  M2LHandlerSArcmp(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerSArcmp"),
      M2LComputerSA<DIM,ORDER,T>(M2LHandler<DIM,T>::fcount)
  {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the SARCMP variant " << std::endl; }


  /// Compute **M2L**
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double eps)
  {
    M2LComputerSA<DIM,ORDER,T>::compute(ffield, level, kernel, eps);
    M2LComputerSA<DIM,ORDER,T>::recompress(eps);
  }



  //! \name Qu C Qb'
  /* @{ */
  /// Expand potentials X -> x
  void applyL2Lpart(T *const X, const unsigned int c=0) const
  { M2LComputerSA<DIM,ORDER,T>::applyU(X, c); }

  /// Compressed M2L operation
  void apply(const unsigned int t_id, const T *const Y, T *const X,
             const unsigned int c=0) const
  { M2LComputerSA<DIM,ORDER,T>::applyC(t_id, Y, X, c); }

  /// Compress densities Y -> Y+nnodes
  void applyM2Mpart(T *const Y, const unsigned int c=0) const
  { M2LComputerSA<DIM,ORDER,T>::applyB(Y, c); }
  /* @} */

};









//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// VARIANTS OF MULTIPLE APPROXIMATIONS ///////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// FULL M2L OPERTOR PER INTERACTION //////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*! 
  \ingroup mm2l

  \brief Full set of \f$K\f$
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerNA
  : public M2LHandler<DIM,T>,
    public TransferHandlerNoSym<DIM>,
    public M2LComputerNA<DIM,ORDER,T>
{
  // used member
  using M2LHandler<DIM,T>::level;
  using M2LHandler<DIM,T>::fcount;
  using TransferHandlerNoSym<DIM>::ffield;
  using M2LComputerNA<DIM,ORDER,T>::operators;


public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;

  /// Constructor
  M2LHandlerNA(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerNA") {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the NA variant " << std::endl; }


  /// Call **M2L** computer
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double /* not needed */)
  {
    M2LComputerNA<DIM,ORDER,T>::compute(ffield, level, kernel);
    TransferHandlerNoSym<DIM>::clearFfield();
  }


  /// Compute X += K Y
  void apply(const unsigned int t_id,
             const T *const Y,
             T *const X,
             const unsigned int /* not needed */) const
  {
    static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
    const T *const K = operators.at(t_id);
    this->applyK(nnodes, nnodes, const_cast<T *const>(K),
								 const_cast<T *const>(Y), X);

    // count flops
    fcount.addFlops(nnodes*(2*nnodes-1) + nnodes);
  }

};









//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// LOW RANK REPRESENTATION PER CONE (IA) /////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*! 
  \ingroup mm2l

  \brief Full set of \f$K \sim UV^*\f$ approximation
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerIA
  : public M2LHandler<DIM,T>,
    public TransferHandlerNoSym<DIM>,
    public M2LComputerIA<DIM,ORDER,T>
{
  // used member
  using M2LHandler<DIM,T>::level;
  using M2LHandler<DIM,T>::fcount;
  using TransferHandlerNoSym<DIM>::ffield;
  using M2LComputerIA<DIM,ORDER,T>::operators;


public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;


  /// Constructor
  M2LHandlerIA(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerIA") {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the IA variant " << std::endl; }


  /// Compute **M2L**
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double eps)
  {
    M2LComputerIA<DIM,ORDER,T>::compute(ffield, level, kernel, eps);
    TransferHandlerNoSym<DIM>::clearFfield();
  }


  /// Compute X += UV' Y  
  void apply(const unsigned int t_id,
             const T *const Y,
             T *const X,
             const unsigned int /* not needed */) const
  {
    static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
    typedef typename M2LComputerIA<DIM,ORDER,T>::m2l_vec::value_type data_pair;

    const data_pair pair = operators.at(t_id);
    const unsigned int k = pair.first;
    const T *const U = pair.second;
    const T *const V = U + k*nnodes;
    applyUV(k, nnodes, U, nnodes, V, Y, X);

    // count flops
    fcount.addFlops(nnodes*(2*k-1) + k*(2*nnodes-1) + nnodes);
  }

};














//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// FULL M2L OPERTOR PER INTERACTION WITH SYMMETRIES //////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


/*! 
  \ingroup mm2l

  \brief Reduced set of \f$K\f$ due to the consideration of symmetries
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerNAsym
  : public M2LHandler<DIM,T>,
    public TransferHandlerSym<ORDER,DIM>,
    public M2LComputerNA<DIM,ORDER,T>
{
  // used member
  using M2LHandler<DIM,T>::level;
  using M2LHandler<DIM,T>::fcount;
  using TransferHandlerSym<ORDER,DIM>::ffield;
  using TransferHandlerSym<ORDER,DIM>::permut;
  using M2LComputerNA<DIM,ORDER,T>::operators;


public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;


  /// Constructor
  M2LHandlerNAsym(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerNAsym")
  {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the NASYM variant " << std::endl; }


  /// Compute **M2L**
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double /* not needed */)
  {
    M2LComputerNA<DIM,ORDER,T>::compute(ffield, level, kernel);
    TransferHandlerSym<ORDER,DIM>::clearFfield();
  }


  /// Compute X += K Y (with permutations)
  void apply(const unsigned int id_orig,
             const T *const Y,
             T *const X,
             const unsigned int /* not needed */) const
  {
    static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
    typedef typename TransferHandlerSym<ORDER,DIM>::perm_pair perm_pair;
    const perm_pair& pp = permut.at(id_orig);
    
    const T *const K = operators.at(pp.first);
    const unsigned int *const pvec = pp.second;

    T Xp[nnodes] = {T(0.)};
    T Yp[nnodes];

    // permute
    for (unsigned int n=0; n<nnodes; ++n) Yp[pvec[n]] = Y[n];
    // apply K
    this->applyK(nnodes, nnodes, K, Yp, Xp);
    // permute back
    for (unsigned int n=0; n<nnodes; ++n) X[n] += Xp[pvec[n]];

    // count flops
    fcount.addFlops(nnodes*(2*nnodes-1) + nnodes);
  }

};







//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// LOW RANK REPRESENTATION PER CONE (IA) WITH SYMMETRIES /////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*! 
  \ingroup mm2l

  \brief Reduced set of \f$K \sim UV^*\f$ approximations due to the
  consideration of symmetries
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerIAsym
  : public M2LHandler<DIM,T>,
    public TransferHandlerSym<ORDER,DIM>,
    public M2LComputerIA<DIM,ORDER,T>
{
  // used member
  using M2LHandler<DIM,T>::level;
  using M2LHandler<DIM,T>::fcount;
  using TransferHandlerSym<ORDER,DIM>::ffield;
  using TransferHandlerSym<ORDER,DIM>::permut;
  using M2LComputerIA<DIM,ORDER,T>::operators;


public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;


  //! Constructor
  M2LHandlerIAsym(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerIAsym")
  {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the IASYM variant " << std::endl; }


  /// Compute **M2L**
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double eps)
  {
    M2LComputerIA<DIM,ORDER,T>::compute(ffield, level, kernel, eps);
    TransferHandlerSym<ORDER,DIM>::clearFfield();
  }


  /// Compute X += UV' Y (with permutations)
  void apply(const unsigned int id_orig,
             const T *const Y,
             T *const X,
             const unsigned int /* not needed */) const
  {
    static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
    typedef typename TransferHandlerSym<ORDER,DIM>::perm_pair perm_pair;
    typedef typename M2LComputerIA<DIM,ORDER,T>::m2l_vec::value_type data_pair;

    const perm_pair& pp = permut.at(id_orig);
    
    const data_pair pair = operators.at(pp.first);
    const unsigned int k = pair.first;
    const T *const U = pair.second;
    const T *const V = U + k*nnodes;

    const unsigned int *const pvec = pp.second;

    T Xp[nnodes] = {T(0.)};
    T Yp[nnodes];

    // permute
    for (unsigned int n=0; n<nnodes; ++n) Yp[pvec[n]] = Y[n];
    // apply UV'
    applyUV(k, nnodes, U, nnodes, V, Yp, Xp);
    // permute back
    for (unsigned int n=0; n<nnodes; ++n) X[n] += Xp[pvec[n]];

    // count flops
    fcount.addFlops(nnodes*(2*k-1) + k*(2*nnodes-1) + nnodes);
  }

};










//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// WITH SYMMETRIES AND BLOCKING //////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



/*!  \class BlockingHandler
  
  \ingroup mm2l
  
  \brief Handles the blocking scheme of expansions if symmetries are used

  The \a BlockingHandler handles the blocking of expansions if the symmetries
  in the m2l kernels are considered. By its means single \c gemm - functions
  can be used instead of multiple \c gemv - functions (blas3 instead of blas2
  functions).
 */
template <int DIM,
          int ORDER,
          typename T>
class BlockingHandler
{
  static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

  static const unsigned int NEXP = 128;
  mutable T* Xp;

  /*!  /typedef EXP_TUPLE_ACCESSOR
    
    This enumerator is used to access the \a exp_tuple. \a CEXP stands for
    expansion counter (how many expansion vectors are inserted yet); \a PVEC
    stands for permutation vectors; \a MEXP stands for \b permuted multipole
    expansions; \a LEXP stands for \a non-permuted local expansions.
   */
  enum EXP_TUPLE_ACCESSOR {CEXP=0, PVEC=1, MEXP=2, LEXP=3};
  typedef boost::tuple<
    unsigned int,         // (CEXP) expansion counter
    const unsigned int**, // (PVEC) array of permutation vectors
    T*,                   // (MEXP) matrix of permuted multipole expansions
    T**                   // (LEXP) matrix of (non permuted) local expansions
    > exp_tuple;

  /*!
    \var blocked
    
    The size of this vector corresponds to the reduced number of unique m2l
    kernels due to the symmetries. Each entry contains a tuple of type
    exp_typle. It is \a mutable because blocked expansions need to be
    changed. Eventually a better workaround will be implemented.
   */
  mutable std::vector<exp_tuple> blocked;



  /// Differs for K and UV', hence, pure virtual
  virtual void applyBlockedM2L(const unsigned int, const unsigned int,
                               const T *const, T *const) const = 0;



  /// apply blocked m2l operator
  void applyBlockedM2L(exp_tuple& et, const unsigned int id_perm) const
  {
    // get expansion counter and permuted multipole expansions
    const unsigned int cexp = boost::get<CEXP>(et);
    T *const Yp = boost::get<MEXP>(et);
    applyBlockedM2L(id_perm, cexp, Yp, Xp);

    // get arrays of permutation vectors and local expansions
    const unsigned int *const *const pvecs = boost::get<PVEC>(et);
    T *const *const Xs                     = boost::get<LEXP>(et);

    // permute back
    for (unsigned int p=0; p<cexp; ++p) {
      const unsigned int *const pvec = pvecs[p];
      T *const X = Xs[p];
      const T *const Xp_ = Xp + p*nnodes;
      for (unsigned int n=0; n<nnodes; ++n) X[n] += Xp_[pvec[n]];
    }

    // reset multipole expansion counter to 0
    boost::get<CEXP>(et) = 0;
  }



protected:
  
  /*! Constructor
    
    \attention The BlockingHandler is not yet thread-safe!
    
    \todo The BlockingHandler needs to be made tread-safe; an idea is to
    instantiate one per thread.
   */
  BlockingHandler()
  {
#ifdef DFMM_USE_OMP
    if (omp_get_max_threads() > 1)
      throw std::runtime_error("BlockingHandler is not yet thread-safe! Switch off OpenMP or use only one thread!"); 
#endif
  }

  /*!
    \brief Activate blocked expansions

    @param[in] nm2l number of m2l kernels
  */
  void activate(const unsigned int nm2l) const
  {
    assert(blocked.empty());

    // allocate array for permuted local expansions
    Xp = new T [nnodes * NEXP];

    // allocate blocked expansions for all existing m2l operators
    blocked.clear();
    for (unsigned int p=0; p<nm2l; ++p)
      blocked.push_back(boost::make_tuple(0,                             // CEXP
                                          new const unsigned int* [NEXP],// PVEC
                                          new T  [nnodes * NEXP],        // MEXP
                                          new T* [NEXP]));               // LEXP
  }

  
  /// Deactivate blocked expansions
  void deactivate() const
  {
    // finish applying blocked m2l operators
    const unsigned int nblocked = blocked.size();
    for (unsigned int id_perm=0; id_perm<nblocked; ++id_perm)
      applyBlockedM2L(blocked.at(id_perm), id_perm);

    // free allocated memory
    if (Xp) delete [] Xp;
    for (unsigned int p=0; p<nblocked; ++p) {
      exp_tuple& et = blocked.at(p);
      delete [] boost::get<PVEC>(et);
      delete [] boost::get<MEXP>(et);
      delete [] boost::get<LEXP>(et);
    }
    blocked.clear();
  }


  /*! 
    \brief Compute X += K Y (with permutations)

    @param[in] id_perm id of m2l kernel
    @param[in] pvec permutation vector
    @param[in] Y multipole expansion
    @param[in,out] X local expansion
   */
  void apply(const unsigned int id_perm,
             const unsigned int *const pvec,
             const T *const Y,
             T *const X) const
  {
    // get blocked expansion tuple
    exp_tuple& et = blocked.at(id_perm);

    // get multipole expansion counter
    const unsigned int cexp = boost::get<CEXP>(et);

    // store (non permuted) local expansion and permutation vector
    boost::get<LEXP>(et)[cexp] = X;
    boost::get<PVEC>(et)[cexp] = pvec;

    // permute multipole expansion
    T*const Yp = boost::get<MEXP>(et) + nnodes * cexp;
    for (unsigned int n=0; n<nnodes; ++n) Yp[pvec[n]] = Y[n];

    // apply blocked m2l operator if block max-block-size is achieved
    if (++(boost::get<CEXP>(et)) == NEXP) applyBlockedM2L(et, id_perm);
  }

};








/*! 
  \ingroup mm2l

  \brief Blocked reduced set of \f$K\f$ due to the consideration of symmetries
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerNAblk
  : public M2LHandler<DIM,T>,
    public TransferHandlerSym<ORDER,DIM>,
    public M2LComputerNA<DIM,ORDER,T>,
    public BlockingHandler<DIM,ORDER,T>
{
  static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

  // used member
  using M2LHandler<DIM,T>::level;
  using M2LHandler<DIM,T>::fcount;
  using TransferHandlerSym<ORDER,DIM>::ffield;
  using TransferHandlerSym<ORDER,DIM>::permut;
  using M2LComputerNA<DIM,ORDER,T>::operators;


  /// apply blocked m2l operator (K)
  void applyBlockedM2L(const unsigned int id_perm,
                       const unsigned int cexp,
                       const T *const Yp, T *const Xp) const
  {
    blas::gemm(nnodes, nnodes, cexp, T(1.),
               const_cast<T*>(operators.at(id_perm)), nnodes,
               const_cast<T*>(Yp), nnodes, Xp, nnodes);

    // count flops
    fcount.addFlops((nnodes*(2*nnodes-1)*cexp) + cexp*nnodes);
  }

public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;


  /// Constructor
  M2LHandlerNAblk(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerNAblk") {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the NABLK variant " << std::endl; }


  /// Activate blocked expansions
  void activate() const
  { BlockingHandler<DIM,ORDER,T>::activate(operators.size()); }


  /// Deactivate blocked expansions
  void deactivate() const
  { BlockingHandler<DIM,ORDER,T>::deactivate(); }


  /// Compute m2l kernels
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double /* not needed */)
  { 
    M2LComputerNA<DIM,ORDER,T>::compute(ffield, level, kernel);
    TransferHandlerSym<ORDER,DIM>::clearFfield();
  }
  
  
  /// Compute X += K Y (with permutations)
  void apply(const unsigned int id_orig,
             const T *const Y,
             T *const X,
             const unsigned int /* not needed */) const
  {
    typedef typename TransferHandlerSym<ORDER,DIM>::perm_pair perm_pair;
    const perm_pair& pp = permut.at(id_orig);
    BlockingHandler<DIM,ORDER,T>::apply(pp.first, pp.second, Y, X);
  }

};










/*! 
  \ingroup mm2l

  \brief Blocked reduced set of \f$K \sim UV^*\f$ approximations due to the
  consideration of symmetries
 */
template <int DIM,
          int ORDER,
          typename T>
class M2LHandlerIAblk
  : public M2LHandler<DIM,T>,
    public TransferHandlerSym<ORDER,DIM>,
    public M2LComputerIA<DIM,ORDER,T>,
    public BlockingHandler<DIM,ORDER,T>
{
  static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

  // used member
  using M2LHandler<DIM,T>::level;
  using M2LHandler<DIM,T>::fcount;
  using TransferHandlerSym<ORDER,DIM>::ffield;
  using TransferHandlerSym<ORDER,DIM>::permut;
  using M2LComputerIA<DIM,ORDER,T>::operators;


  /// apply blocked m2l operator (UV)
  void applyBlockedM2L(const unsigned int id_perm,
                       const unsigned int cexp,
                       const T *const Yp, T *const Xp) const
  {
    typedef typename M2LComputerIA<DIM,ORDER,T>::m2l_vec::value_type data_pair;
    const data_pair pair = operators.at(id_perm);
    const unsigned int k = pair.first;
    const T *const U = pair.second;
    const T *const V = U + k*nnodes;
    
    T *const YX = new T [k * cexp];
    blas::gemhm(nnodes, k, cexp, T(1.), const_cast<T*>(V), nnodes,
                const_cast<T*>(Yp), nnodes, YX, k);
    blas::gemm( nnodes, k, cexp, T(1.), const_cast<T*>(U), nnodes,
                YX, k, Xp, nnodes);
    delete [] YX;

    // count flops
    fcount.addFlops( ((nnodes*(2*k-1) + k*(2*nnodes-1)) + nnodes) * cexp);
  }





public:
  enum {dim = DIM, order = ORDER};
  typedef T value_type;


  //! Constructor
  M2LHandlerIAblk(const Level<DIM> *const level)
    : M2LHandler<DIM,T>(level, "M2LHandlerIAblk")
  {}


  /// Write out info
  static void getInfo() 
  { std::cout << "\n- Using the IABLK variant " << std::endl; }


  /// Activate blocked expansions
  void activate() const
  { BlockingHandler<DIM,ORDER,T>::activate(operators.size()); }


  /// Deactivate blocked expansions
  void deactivate() const
  { BlockingHandler<DIM,ORDER,T>::deactivate(); }


  /// Compute **M2L**
  template <typename kernel_type>
  void computeM2L(const kernel_type& kernel, const double eps)
  {
    M2LComputerIA<DIM,ORDER,T>::compute(ffield, level, kernel, eps);
    TransferHandlerSym<ORDER,DIM>::clearFfield();
  }


  /// Compute X += UV' Y (with permutations)
  void apply(const unsigned int id_orig,
             const T *const Y,
             T *const X,
             const unsigned int /* not needed */) const
  {
    typedef typename TransferHandlerSym<ORDER,DIM>::perm_pair perm_pair;
    const perm_pair& pp = permut.at(id_orig);
    BlockingHandler<DIM,ORDER,T>::apply(pp.first, pp.second, Y, X);
  }

};



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////




#include <boost/foreach.hpp>
#include <omp.h>




template <typename kernel_type,
          typename m2l_handler_type>
void ComputeM2Loperators(std::vector<m2l_handler_type>& m2lv,
                         const kernel_type& kernel,
                         const double eps)
{
  std::cout << "\n- Evaluating M2L operators " << std::endl;
  const double t0 = omp_get_wtime();
  for (unsigned int l=0; l<m2lv.size(); ++l) {
    if (m2lv.at(l).hasInteractions()) {
      const double t = omp_get_wtime();
      std::cout << "  at " << l << " with ";
      m2lv.at(l).computeM2L(kernel, eps);
      std::cout << "took " << omp_get_wtime() - t << " s" << std::endl;
    }
  }
  std::cout << "  - overall " << omp_get_wtime() - t0 << " s" << std::endl;
}


#endif

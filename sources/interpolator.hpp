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

#ifndef interpolator_hpp
#define interpolator_hpp


#include <vector>

#include <boost/shared_array.hpp>
#include <boost/foreach.hpp>


template <int> class map_loc_glob;
template <int> class map_glob_loc;
template <int> struct Chebyshev;
template <int,int> class Tensor;



/*! \defgroup mtranslation M2M, L2L: Translation Operators
 
  The translation operators are

  \arg P2M: Particle-To-Multipole Operator
  \arg M2M: Multipole-To-Multipole Operator
  \arg L2L: Local-To-Local Operator
  \arg L2P: Local-To-Particle Operator

  and they translate the \a multipole and \local expansions from
  child-to-parent (upward pass), respectively, from parent-to-child (downward
  pass) level.
 */






/*! \ingroup mtranslation

  \class IOcomputerLF 
  
  \brief Computes interpolation operators for leaf clusters

  \todo Find a better name and clean up or even better refactor

  @tparam DIM spatial dimension
  @tparam ORDER interpolation order
  @tparam T value type
 */
template <int DIM, int ORDER, typename T>
class IOcomputerLF :
  boost::noncopyable
{
  // compile time constants and types
  enum {nnodes = BasisTraits<ORDER,DIM>::nnodes};

  typedef typename DimTraits<DIM>::point_type point_type;

  typedef Chebyshev< ORDER>  basis_type;
  typedef Tensor<DIM,ORDER> tensor_type;
  typedef boost::array<double[ORDER],ORDER> order_order_array;
  typedef boost::array<double[DIM],  ORDER>   order_dim_array;

  // members
  const map_glob_loc<DIM> map;

  const order_order_array T_of_roots;
  const static order_order_array eval_T_of_roots()
  {
    order_order_array values;
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        values[o][j] = basis_type::T(o, basis_type::nodes[j]);
    return values;
  }
  
  
  
public:

  // ctor for a cluster of given center and diam
  explicit IOcomputerLF(const point_type& center, const double diam)
    : map(center, diam), T_of_roots(eval_T_of_roots()) {}


  // ctor for unit cluster [-1,1]^d
  explicit IOcomputerLF()
    : map(), T_of_roots(eval_T_of_roots()) {}


  // public overloaded operator ()
  void operator()(const point_type& x, T S[nnodes]) const
  {
    // check if x is in current cluster
    assert(map.is_in_cluster(x));
    
    // evaluate T of x
    order_dim_array T_of_x;
    for (unsigned int d=0; d<DIM; ++d) {
      const double xi = map(d, x[d]);
      for (unsigned int o=1; o<ORDER; ++o)
        T_of_x[o][d] = basis_type::T(o, xi);
    }
    
    // assemble S for x and all interpolation nodes
    for (unsigned int n=0; n<nnodes; ++n) {
      S[n] = 1.;
      for (unsigned int d=0; d<DIM; ++d) {
        const unsigned int j = tensor_type::node_ids[n][d];
        double S_d = 1. / ORDER;
        for (unsigned int o=1; o<ORDER; ++o)
          S_d += 2. / ORDER * T_of_x[o][d] * T_of_roots[o][j];
        S[n] *= S_d;
      }
    }
    
  }

  
  /*!
    Assembles the P2M/L2P opertor for particles
    @param[in] nparticles number of particles
    @param[in] pindices reordered particle indices
    @param[in] particles array of particles
    @param[out] S P2M/L2P matrix to assemble
   */
  template <typename value_type,
            typename particle_type>
  void operator()(const unsigned int nparticles,
                  const unsigned int *const pindices,
                  const particle_type *const particles,
                  value_type *const S) const
  {
    value_type S_[nnodes];
    for (unsigned int i=0; i<nparticles; ++i) {
      operator()(particles[pindices[i]].getPoint(), S_);
      assert(check_nan(nnodes, S_));
      for (unsigned int j=0; j<nnodes; ++j)
        S[j*nparticles + i] = S_[j];
    }
  }


  /*!
    Assembles the P2M/L2P opertor for particles
    @param[in] npoints number of points
    @param[in] points array of points
    @param[out] S P2M/L2P matrix to assemble
   */
  template <typename value_type>
  void operator()(const unsigned int npoints,
                  const point_type *const points,
                  value_type *const S) const
  {
    value_type S_[nnodes];
    for (unsigned int i=0; i<npoints; ++i) {
      operator()(points[i], S_);
      assert(check_nan(nnodes, S_));
      for (unsigned int j=0; j<nnodes; ++j)
        S[j*npoints + i] = S_[j];
    }
  }

};





/*! \ingroup mtranslation

  \class IOsetterNL
  
  \brief Computes interpolation operators not for leaf clusters

  Interpolation operator setter of (N)on (L)eaf clusters. Since the grid is
  regular, all non leaf clusters of the same level have identical
  interpolation operators. Hence, memory for the interpolation operator S is
  only allocated here but not freed.

  \todo Find a better name and clean up or even better refactor

  @tparam DIM spatial dimension
  @tparam ORDER interpolation order
  @tparam T value type
 */
template <int DIM,
          int ORDER,
          typename T>
class IOsetterNL
{
  enum {nboxes = DimTraits<DIM>::nboxes,
        nnodes = BasisTraits<ORDER,DIM>::nnodes};
  typedef typename DimTraits<DIM>::point_type point_type;

  typedef Chebyshev< ORDER>  basis_type;
  typedef Tensor<DIM,ORDER> tensor_type; 

  boost::shared_array<T> S;
  
public:
  explicit IOsetterNL() : S(new T [nboxes*nnodes * nnodes])
  {
    // setup interpolation nodes of child clusters
    point_type x [nboxes*nnodes];
    for (unsigned int b=0; b<nboxes; ++b) {
      point_type center;
      for (unsigned int d=0; d<DIM; ++d)
        center[d] = 0.5 * DimTraits<DIM>::child_pos[b][d];
      const map_loc_glob<DIM> map(center, 1.0);
      
      for (unsigned int n=0; n<nnodes; ++n)
        for (unsigned int d=0; d<DIM; ++d)
          x[b*nnodes + n][d]
            = map(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
    }

    // compute S
    IOcomputerLF<DIM,ORDER,T> cmp;
    for (unsigned int b=0; b<nboxes; ++b)
      cmp(nnodes, x + b*nnodes, S.get() + b*nnodes*nnodes);

  }

  const boost::shared_array<T>& getS() const
  { return S; }

};




template <int ORDER,
          typename T>
class IOsetterNL<3,ORDER,T>
{
  enum {DIM = 3,
        nboxes = DimTraits<DIM>::nboxes,
        nnodes = BasisTraits<ORDER,DIM>::nnodes};
  //typedef typename DimTraits<DIM>::point_type point_type;
  typedef boost::array<double[ORDER],ORDER> order_order_array;
  
  typedef Chebyshev< ORDER>  basis_type;
  typedef Tensor<DIM,ORDER> tensor_type; 

  // own members
  boost::shared_array<T> S;
  const order_order_array T_of_roots;


  static
  const order_order_array eval_T_of_roots()
  {
    order_order_array values;
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        values[o][j] = basis_type::T(o, basis_type::nodes[j]);
    return values;
  }


  void assemble(const double *const x, T *const S) const
  {
    // values of chebyshev polynomials of source particle: T_o(x_i)
    double T_of_x[ORDER];
    // loop: local points (mapped in [-1,1])
    for (unsigned int m=0; m<ORDER; ++m) {
      // evaluate chebyshev polynomials at local points
      for (unsigned int o=1; o<ORDER; ++o) T_of_x[o] = basis_type::T(o, x[m]);
      // assemble S
      for (unsigned int n=0; n<ORDER; ++n) {
        S[n*ORDER + m] = 1. / ORDER;
        for (unsigned int o=1; o<ORDER; ++o)
          S[n*ORDER + m] += 2. / ORDER * T_of_x[o] * T_of_roots[o][n];
      }
    }
  }
  

  
public:
  explicit IOsetterNL()
    : S(new T [nboxes * 3*ORDER*ORDER]),
      T_of_roots(eval_T_of_roots())
  {
    // loop over all child boxed
    for (unsigned int b=0; b<nboxes; ++b) {
      typename DimTraits<DIM>::point_type center;
      for (unsigned int d=0; d<DIM; ++d)
        center[d] = 0.5 * DimTraits<DIM>::child_pos[b][d];
      const map_loc_glob<DIM> map(center, 1.0);
      
      // assemble S for each child box
      T *const Sb = S.get() + b * 3*ORDER*ORDER;
      for (unsigned int d=0; d<DIM; ++d) {
        double x[ORDER];
        for (unsigned int o=0; o<ORDER; ++o)
          x[o] = map(d, basis_type::nodes[o]);
        assemble(x, Sb + d * ORDER*ORDER);
      }
    }
  }

  const boost::shared_array<T>& getS() const
  { return S; }

};















//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



/*!
  \ingroup mtranslation
*/
template <typename cluster_type,
          typename particle_type,
          typename clusterbasis_type,
          typename m2l_handler_type>
void setIO(const std::vector<std::vector<cluster_type*> >& clvv,
           std::vector<clusterbasis_type>& basis,
           std::vector<m2l_handler_type>& m2lv,
           const unsigned int *const pindices,
           const particle_type *const particles)
{
  enum {dim    = clusterbasis_type::dim,
        order  = clusterbasis_type::order,
        nnodes = clusterbasis_type::nnodes};
  typedef typename clusterbasis_type::value_type value_type;
  typedef IOcomputerLF<dim,order,value_type> ioc_type;
  
  for (unsigned int l=0; l<clvv.size(); ++l) {
    const std::vector<cluster_type*>& clv = clvv.at(l);
    m2l_handler_type& m2l = m2lv.at(l);
    const double ext = m2l.getLevel()->getDiam();
    if (m2l.hasInteractions()) {
      switch (m2l.getLevel()->isLeaf()) {
      case true: {
        BOOST_FOREACH(const cluster_type *const cl, clv) {
          ioc_type ioc(cl->getCenter(), ext);
          boost::shared_array<value_type>
            S(new value_type [cl->getSize()*nnodes]);
          ioc(cl->getSize(), pindices+cl->getNbeg(), particles, S.get());
          basis.at(cl->getCidx()).assign(S);
        }
      } break;
      case false: {
        IOsetterNL<dim,order,value_type> ios;
        BOOST_FOREACH(const cluster_type *const cl, clv)
          basis.at(cl->getCidx()).assign(ios.getS());
      } break;
      }
    }
  }
}



#endif

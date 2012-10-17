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

#ifndef cluster_hpp
#define cluster_hpp



#include <boost/noncopyable.hpp>
#include <boost/timer.hpp>
#include <boost/array.hpp>

#include "sources/dimordertraits.hpp"


/*! 
  \defgroup mcluster Clusters and Related Things
 */


/*! \ingroup mcluster

  \class Cluster
 
  \brief Partitions of the bounding box including the domain of interest
 
  Clusters are the main ingredient of the tree partition of the bounding box
  which includes the computational domain. The root-cluster corresponds to
  this bounding box. Each cluster knows only information which is not shared
  with other clusters, all other informations are stored in Level.
 
  @tparam DIM spatial dimension
 */
template <int DIM>
class Cluster
  : private boost::noncopyable
{
public:
  enum {dim = DIM, nboxes = DimTraits<DIM>::nboxes};

  typedef typename DimTraits<DIM>::point_type     point_type;
  typedef typename DimTraits<DIM>::point_int_type point_int_type;


private:
  typedef Level<DIM> level_type;
  typedef boost::array<Cluster*,nboxes> cluster_array;


  Cluster *const parent;       //!< parent cluster
  const unsigned int cidx;     //!< overall cluster index in tree
  const unsigned int bidx;     //!< child box index
  const unsigned int lidx;     //!< overall cluster index in level
  const unsigned int nbeg;     //!< starting index in particle index array
  const unsigned int nend;     //!< ending index in particle index array
  const unsigned int size;     //!< number of particles in cluster
  const point_type center;     //!< center of cluster
  const point_int_type pos;    //!< relative level coordinate of cluster
  const unsigned int nlevel;   //!< level of cluster
  const bool is_leaf;          // !< is cluster leaf or not


  /**
   * List of pointers to child clusters: The size of the list depends only on
   * nboxes, which is a constant depending only on the spatial dimension (see
   * \ref mtraits module). Note that if a pointer points to @p NULL there
   * exists no child cluster with the corresponding box index.
   */
  cluster_array child_list;



  template <typename particle_type>
  static void subdivide(const std::vector<level_type*>& levels,
                        std::vector<std::vector<Cluster*> >& clvv,
                        Cluster *const cl,
                        const particle_type *const particles,
                        unsigned int *const pindices,
                        const unsigned int lmax,
                        unsigned int& cluster_counter);


public:


  /**
   * This constructor is used within a tree (not allowed for the root
   * cluster). Then @p parent cluster and box index @p boxidx are relevant.
   *
   * @param[in] parent parent cluster
   * @param[in] cidx unique cluster index (overall)
   * @param[in] bidx child index relative to parent cluster
   * @param[in] lidx unique cluster index (level)
   * @param[in] nbeg index indicating beginning in particle array
   * @param[in] nend index indicating ending in particle array
   * @param[in] center center of cluster
   * @param[in] position position of cluster at current level
   * @param[in] nlevel in tree
   * @param[in] is_leaf at leaf-level or not
   */
  explicit Cluster(Cluster* parent,
                   const unsigned int cidx,
                   const unsigned int bidx,
                   const unsigned int lidx,
                   const unsigned int nbeg,
                   const unsigned int nend,
                   const point_type center,
                   const point_int_type position,
                   const unsigned int nlevel,
                   const bool is_leaf);
  


  /**
   * This constructor is used for the root cluster. Only local dof iterators,
   * the root center and root diam are provided.
   *
   * @param[in] N number of particles
   * @param[in] root_center center of bounding box (corresponds to root cluster)
   */
  explicit Cluster(const unsigned int N,
                   const point_type root_center);
  
  


  // Destructor
  virtual ~Cluster()
  {
    for (unsigned int b=0; b<nboxes; ++b)
      if (child_list[b]) { delete child_list[b]; child_list[b] = NULL; }
  }
  



  // call recursive subdivision and fill level based cluster list
  template <typename particle_type>
  const unsigned int subdivide(const std::vector<level_type*>& levels,
                               std::vector<std::vector<Cluster*> >& clvv,
                               const particle_type *const particles,
                               unsigned int *const pindices,
                               const unsigned int lmax)
  {
    assert(cidx==0);
    assert(clvv.empty());

    for (typename std::vector<level_type*>::const_iterator it=levels.begin();
         it!=levels.end(); ++it) clvv.push_back(std::vector<Cluster*>());
    
    clvv.front().push_back(this);

    unsigned int cluster_counter = 1;
    subdivide(levels, clvv, this, particles, pindices, lmax, cluster_counter);
    return cluster_counter;
  }


  /*! @return level index */
  const unsigned int getNlevel() const
  { return nlevel; }

  /*! @return unique cluster index */
  const unsigned int getCidx() const
  { return cidx; }

  /*! @return child index */
  const unsigned int getBidx() const
  { return bidx; }

  /*! @return unique cluster index at its level */
  const unsigned int getLidx() const
  { return lidx; }

  /*! @return starting index */
  const unsigned int getNbeg() const
  { return nbeg; }

  /*! @return ending index */
  const unsigned int getNend() const
  { return nend; }

  /*! @return cluster size */
  const unsigned int getSize() const
  { return size; }

  /*! @return cluster center */
  const point_type& getCenter() const
  { return center; }

  /*! @return cluster position */
  const point_int_type& getPosition() const
  { return pos; }

  /*! @return if cluster is leaf or not */
  const bool isLeaf() const
  { return is_leaf; }

  /*! @return parent cluster */
  const Cluster *const getParent() const
  { return parent; }
  /*! @return parent cluster */
  Cluster *const getParent()
  { return parent; }

  /*! @return child cluster of \a bidx =  b */
  const Cluster *const getChild(const unsigned int b) const
  { return child_list.at(b); }
  /*! @return child cluster of \a bidx =  b */
  Cluster *const getChild(const unsigned int b)
  { return child_list.at(b); }

  
//  static
//  const unsigned int getTransferId(const unsigned int nlevel,
//                                   const point_int_type& relpos)
//  {
//    const unsigned int size = 2 << (nlevel-1);
//    unsigned int tid = 0;
//    unsigned int scale = 1;
//    for (unsigned int d=0; d<DIM; ++d) {
//      tid += (relpos[d] + size) * scale;
//      scale *= size;
//    }
//    return tid;
//  }


//  const unsigned int getTransferId(const Cluster *const other) const
//  {
//    assert(nlevel == other->nlevel);
//    const unsigned int size = 2 << (nlevel-1);
//    const point_int_type relpos = other->getPosition() - pos;
//    
//    unsigned int tid = 0;
//    unsigned int scale = 1;
//    for (unsigned int d=0; d<DIM; ++d) {
//      tid += (relpos[d] + size) * scale;
//      scale *= size;
//    }
//    return tid;
//  }


  /*! Writes out info of the cluster */
  std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    unsigned int nc = 0;
    for (unsigned int b=0; b<nboxes; ++b) if (child_list.at(b)) nc++;

    out << "Cluster " << cidx
        << " at level " << nlevel
        << " of bidx " << bidx << " centered at " << center
        << " has " << size << " dofs starting at " << nbeg
        << " and children " << nc << std::endl;
    return out;
  }

};






/*! \class AreParticlesInCluster

  This class (functor) checks if a provided cluster containes the particles
  given by the pindices array.
 */
template <typename particle_type>
class AreParticlesInCluster
{
  typedef typename particle_type::point_type point_type;

  const particle_type *const particles;
  const unsigned int  *const pindices;
  
public:
  /*! constructor
   @param[in] _particles array of particles
   @param[in] _pindices reordered array of indices to access particles
  */
  AreParticlesInCluster(const particle_type *const _particles,
                        const unsigned int  *const _pindices)
    : particles(_particles), pindices(_pindices)
  {}
  
  /*! This functions checks whether the particles given by \a
      pindices[cluster->nbeg] until \a pindices[cluster->nend] are in the
      given cluster at the given level */
  template <typename cluster_type,
            typename level_type>
  bool operator()(const cluster_type *const cluster,
                  const level_type& level) const
  {
    const point_type& center = cluster->getCenter();
    const double diam   = level->getDiam();
    for (unsigned int i=0; i<cluster->getSize(); ++i) {
      const point_type& point
        = particles[pindices[cluster->getNbeg()+i]].getPoint();
      for (unsigned int d=0; d<point_type::dim; ++d)
        if (fabs(center[d]-point[d])-(diam/2.) >
            std::numeric_limits<double>::epsilon()) {
          cluster->writeInfo();
          std::cout << " diam/2=" << diam/2.
                    << " |x-c|=" << fabs(center[d]-point[d])
                    << " diff=" << (diam/2. -fabs(center[d]-point[d]))
                    << " idx=" << i << " dim=" << d << "\t point "
                    << point << std::endl;
          return false;
        }
    }
    return true;
  }
};











// make member definition visible to this header
#include "sources/cluster.tpl"



#endif












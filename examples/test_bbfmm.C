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

// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <utility>
#include <algorithm>
#include <omp.h>

#define LF_MIN_DISTANCE_1
#define HF_MIN_DISTANCE_2


// dfmm
#include "sources/dimordertraits.hpp"
#include "sources/dof.hpp"
#include "sources/plot.hpp"
#include "sources/storage.hpp"
#include "sources/cluster.hpp"

#include "sources/m2lhandler.hpp"
#include "sources/kernelfunctions.hpp"
#include "sources/clusterbasis.hpp"
#include "sources/relations.hpp"
#include "sources/relationsetter.hpp"
#include "sources/expansion.hpp"

#include "sources/psdraw.hpp"

#include "sources/mapping.hpp"
#include "sources/chebyshev.hpp"
#include "sources/tensor.hpp"

#include "sources/interpolator.hpp"

#include "sources/m2mapplyer.hpp"
#include "sources/m2lapplyer.hpp"
#include "sources/l2lapplyer.hpp"
#include "sources/nearfield.hpp"

#include "sources/functions.hpp"





template <typename m2l_type,
          typename kernel_type>
class bbFMM
{
  enum {dim   = m2l_type::dim,
        order = m2l_type::order};
  typedef typename m2l_type::value_type value_type;
  typedef std::vector<m2l_type> m2l_handler_type_vector;
  typedef ClusterBasis<dim,order,value_type> clusterbasis_type;
  typedef ExpansionHandler<dim,order,value_type> expansion_type;
  typedef Relations<dim> relations_type;
  typedef std::vector<std::vector<Cluster<dim>*> > octree_cluster_vector;
  typedef Dof<dim> particle_type;

  const kernel_type& kernel;
  const octree_cluster_vector& cluster_tree;
  const particle_type *const particles;
  const unsigned int *const pindices;

  std::vector<m2l_type> all_m2l_handler;
  std::vector<clusterbasis_type> rbasis, cbasis;
  std::vector<expansion_type> rexph, cexph;
  std::vector<relations_type> relations;


public:

  bbFMM(const unsigned int ncl,
       const kernel_type& _kernel,
       const octree_cluster_vector& _cluster_tree,
       const particle_type *const _particles,
       const unsigned int *const _pindices)
    : kernel(_kernel),
      cluster_tree(_cluster_tree),
      particles(_particles),
      pindices(_pindices),
      rbasis(ncl), cbasis(ncl),
      rexph(ncl), cexph(ncl),
      relations(ncl)
  {}


  void precomputeM2LOperators(const Storage<dim>& storage, const double eps)
  {
    // init m2l handler 
    storage.initialize(all_m2l_handler);

    std::cout << "\n/////////////////////////////////////////////" << std::endl;
    m2l_type::getInfo();

    // setup multilevel relations
    SetupClusterRelations(cluster_tree, rexph,
                          cluster_tree, cexph, 
                          relations, storage.getLevels(), all_m2l_handler);

    // write out info
    storage.writeInfo(all_m2l_handler);

    // compute m2l operators
    ComputeM2Loperators(all_m2l_handler, kernel, eps);

    // l2l/m2m operators
    std::cout << "\n- Evaluating interpolation operators " << std::flush;
    const double t_io = omp_get_wtime();
    setIO(cluster_tree, rbasis, all_m2l_handler, pindices, particles);
    setIO(cluster_tree, cbasis, all_m2l_handler, pindices, particles);
    std::cout << " took " << omp_get_wtime() - t_io << " s" << std::endl;
  }
  
  
  void multiply(const value_type *const y, value_type *const x)
  {
    // m2m
    std::cout << "\n- Applying M2M operators\n" << std::flush;
    const double t_m2m = omp_get_wtime();
    ApplyM2Moperators(cluster_tree, cbasis, cexph, all_m2l_handler, y);
    std::cout << "  - overall " << omp_get_wtime() - t_m2m << " s" << std::endl;

    // m2l
    std::cout << "\n- Applying M2L operators\n" << std::flush;
    const double t_m2l = omp_get_wtime();
    ApplyM2Loperators(cluster_tree, rexph, cexph, relations, all_m2l_handler);
    std::cout << "  - overall " << omp_get_wtime() - t_m2l << " s" << std::endl;

    // l2l
    std::cout << "\n- Applying L2L operators\n" << std::flush;
    const double t_l2l = omp_get_wtime();
    ApplyL2Loperators(cluster_tree, rbasis, rexph, all_m2l_handler, x);
    std::cout << "  - overall " << omp_get_wtime() - t_l2l << " s" << std::endl;

    // apply near field on the fly
    typedef NFadder<Cluster<dim>,particle_type,relations_type,kernel_type>
      nfc_type;
    ApplyNearField(cluster_tree.back(), nfc_type(kernel, relations,
                                                 pindices, particles,
                                                 pindices, particles,
                                                 x, y));
  }

};







int main( int argc, char * argv[ ] )
{
  // start timer
  const double t_0 = omp_get_wtime( );

  const unsigned int dim   = 3;  
  const unsigned int order = DFMM_ORDER;
  const double eps         = DFMM_EPSILON;

  typedef DimTraits<dim>::point_type point_type;
  typedef Dof<dim>                particle_type;

  // required command line arguments
#ifdef DFMM_USE_OMP
  if (argc != 5) {
#else
  if (argc != 4) {
#endif
    std::cerr << "Wrong number of command line arguments" << std::endl;
    exit(-1);
  }
  const std::string filename(    argv[1]);
  const unsigned long N   = atoi(argv[2]);
  const unsigned int lmax = atoi(argv[3]);

#ifdef DFMM_USE_OMP
  const unsigned int nthreads = atoi(argv[4]);
  omp_set_num_threads(nthreads);
  std::cout << "Using " << omp_get_max_threads() << " threads." << std::endl;
#endif

  // read particles
  particle_type *const particles = new particle_type [N];
  unsigned int  *const pindices  = new unsigned int  [N];
  ReadParticles(filename, N, particles);

  // fill pindices
  for (unsigned int p=0; p<N; ++p)
    pindices[p] = particles[p].getId();

  // get bounding box of root cluster
  const std::pair<point_type,point_type> bbx = GetBoundingBox(particles, N);

  // storage: setup and stores level infos, stores kernel wrapper, etc.
  typedef Storage<dim> storage_type;
  storage_type storage(bbx, lmax);


  // setup root cluster
  typedef Cluster<dim> cluster_type;
  cluster_type* cl = new cluster_type(N, storage.getRootClusterCenter());
  
  // subdivide root cluster, then generated level based cluster lists
  std::vector<std::vector<cluster_type*> > all_cluster_vectors;
  const unsigned int ncl = cl->subdivide(storage.getLevels(),
                                         all_cluster_vectors,
                                         particles, pindices, lmax);
  std::cout << "\n- Number of clusters " << ncl << std::endl;

  // check if particles are in correct cluster
  AreParticlesInCluster<particle_type> areThey(particles,pindices);
  bool all_particles_are_in_cluster = true;
  BOOST_FOREACH(const cluster_type *const lcl, all_cluster_vectors.back())
    if (!areThey(lcl, storage.getLevels().at(lcl->getNlevel()))) {
      std::cout << "Particles are not in cluster they belong to." << std::endl;
      all_particles_are_in_cluster = false;
    }
  assert(all_particles_are_in_cluster);



  // LAPLACE
  typedef KernelFunction<LAPLACE3D> kernel_type;
  typedef kernel_type::value_type value_type;
  kernel_type kernel;
  storage.initLevels();


  // densities and potentials
  value_type *const y = new value_type [N];
  value_type *const x = new value_type [N];

  for (unsigned int n=0; n<N; ++n)
    y[n] = (double)rand()/(double)RAND_MAX;

  for (unsigned int n=0; n<N; ++n) x[n] = 0.;


  // exact clusterX
  std::cout << "\n- Exact evaluation of cluster X " << std::flush;
  const double t_ex = omp_get_wtime();
  cluster_type *const clx = all_cluster_vectors.back().front();
  assert(clx->getNbeg() == 0);
  const unsigned int nx = clx->getSize();
  const unsigned int ny = N;
  value_type *const x0 = new value_type [nx];
  blas::setzero(nx, x0);
  for (unsigned int i=0; i<nx; ++i)
    for (unsigned int j=nx; j<ny; ++j)
      x0[i] += kernel(particles[pindices[i]].getPoint(),
                      particles[pindices[j]].getPoint()) * y[j];
  std::cout << "took " << omp_get_wtime() - t_ex << " s" << std::endl;


  { // SArcmp ////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerSArcmp<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }

  { // NA ////////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerNA<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }

  { // NAsym /////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerNAsym<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }

  { // NAblk /////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerNAblk<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }

  { // IA ////////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerIA<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }

  { // IAsym /////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerIAsym<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }

  { // IAblk /////////////////////////////////////////////////////////
    for (unsigned int n=0; n<N; ++n) x[n] = 0.;
    typedef M2LHandlerIAblk<dim,order,value_type> m2l_handler_type;
    bbFMM<m2l_handler_type, kernel_type>
      fmm(ncl, kernel, all_cluster_vectors, particles, pindices);
    fmm.precomputeM2LOperators(storage, eps);
    fmm.multiply(y, x);
    std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
    std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;
    if (computeL2norm( nx, x0, x) > DFMM_EPSILON*10) return 1;
    if (computeINFnorm(nx, x0, x) > DFMM_EPSILON*10) return 1;
  }



  std::cout << "\n- INTERPOLATION_ORDER = " << order
            << "\n- EPSILON = " << eps << std::endl;



  // clear memory
  delete [] x;
  delete [] y;
  delete [] x0;

  delete cl;
  delete [] particles;                         
  delete [] pindices;
  


  ////////////////////////////////////////////////////
  std::cout << "\n- Overall time " << omp_get_wtime( ) - t_0
            << " s" << std::endl;
  ////////////////////////////////////////////////////


  return 0;
}   

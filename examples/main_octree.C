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

#define FREQUENCY_REGIME_DELIMITER .8
#define LF_MIN_DISTANCE_1
#define HF_MIN_DISTANCE_2
#define SCALE_CONE_APERTURE 2.


// qbdfmm
#include "sources/dimordertraits.hpp"
#include "sources/dof.hpp"
#include "sources/plot.hpp"
#include "sources/storage.hpp"
#include "sources/cluster.hpp"





// reads input files (examples can be found in cmk/grid/...)
template<typename particle_type>
void ReadParticles(const std::string filename,
                   const unsigned int N,
                   particle_type *const particles)
{
  double coords[particle_type::point_type::dim];
  std::ifstream is( filename.c_str() );
  unsigned int id = 0;
  for (unsigned int i=0; i<N; ++i) {
    for (unsigned int d=0; d<particle_type::dim; ++d) is >> coords[d];
    particles[i]
      = particle_type(typename particle_type::point_type(coords),id++);
  }
  assert(id == N);
  is.close();
}







template <typename particle_type>
const std::pair<typename particle_type::point_type,
                typename particle_type::point_type>
GetBoundingBox(const particle_type *const particles, const unsigned int N)
{
  const unsigned int dim = particle_type::point_type::dim;
  typename particle_type::point_type max, min;
  for (unsigned int d=0; d<dim; ++d) {
    max[d] = particles[0].getPoint()[d];
    min[d] = particles[0].getPoint()[d];
    for (unsigned long p=1; p<N; ++p) {
      const double value = particles[p].getPoint()[d];
      if (max[d]<value) max[d] = value;
      if (min[d]>value) min[d] = value;
    }
  }
  return std::make_pair(min,max);
}





int main( int argc, char * argv[ ] )
{
  // start timer
  const double time = omp_get_wtime( );

  const unsigned int dim = 3;  

  typedef DimTraits<dim>::point_type point_type;
  typedef Dof<dim>                particle_type;

  // required command line arguments
  if (argc != 4) {
    std::cerr << "Wrong number of command line arguments" << std::endl;
    exit(-1);
  }
  const std::string filename(    argv[1]);
  const unsigned long N   = atoi(argv[2]);
  const unsigned int lmax = atoi(argv[3]);


  //dof** dofs = new dof* [ndofs];
  particle_type *const particles = new particle_type [N];
  unsigned int  *const pindices  = new unsigned int  [N];
  ReadParticles(filename, N, particles);

  // fill pindices
  for (unsigned int p=0; p<N; ++p)
    pindices[p] = particles[p].getId();

  //// plot particles
  //Plot<particle_type> plotter(particles, N);
  //plotter.plot();


  // get bounding box of root cluster
  const std::pair<point_type,point_type> bbx = GetBoundingBox(particles, N);
  std::cout << "bounding box = [" << bbx.first << "," << bbx.second << "]" 
            << std::endl;

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
  //const cluster_type *const cluster = all_cluster_vectors.back().front();
  BOOST_FOREACH(const cluster_type *const lcl, all_cluster_vectors.back()) {
    //std::cout <<         lcl->getCoord() << " - "
    //          << cluster    ->getCoord() << " = "
    //          << lcl->getCoord() - cluster->getCoord() << std::endl;
    if (!areThey(lcl, storage.getLevels().at(lcl->getNlevel()))) {
      std::cout << "Particles are not in cluster they belong to." << std::endl;
      all_particles_are_in_cluster = false;
    }
  }
  if (all_particles_are_in_cluster)
    std::cout << "All particles are in the clusters they belong to."
              << std::endl;

  // free memory
  delete cl;
  delete [] particles;                         
  delete [] pindices;


  ////////////////////////////////////////////////////
  std::cout << "\n- Overall time " << omp_get_wtime( ) - time
            << " s" << std::endl;
  ////////////////////////////////////////////////////


  return 0;
}   

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

#ifndef cluster_tpl
#define cluster_tpl

#include <boost/foreach.hpp>




// ctor: root cluster
template <int DIM>
Cluster<DIM>::Cluster(const unsigned int root_ndofs,
                      const point_type   root_center)
  : parent(NULL),
    cidx(0),
    bidx(-1),
    lidx(0),
    nbeg(0),
    nend(root_ndofs),
    size(nend-nbeg),
    center(root_center),
		pos(0),
		nlevel(0),
		is_leaf(false)
{
	for (unsigned int c=0; c<nboxes; ++c) child_list[c] = NULL;
}




// ctor
template <int DIM>
Cluster<DIM>::Cluster(Cluster*             _parent,
                      const unsigned int   _cidx,
                      const unsigned int   _bidx,
                      const unsigned int   _lidx,
                      const unsigned int   _nbeg,
                      const unsigned int   _nend,
                      const point_type     _center,
											const point_int_type _pos,
											const unsigned int   _nlevel,
											const bool           _is_leaf)
  : parent(_parent),
    cidx(  _cidx),
		bidx(  _bidx),
    lidx(  _lidx),
    nbeg(  _nbeg),
    nend(  _nend),
    size(nend-nbeg),
		center(_center),
		pos( _pos),
		nlevel(_nlevel),
		is_leaf(_is_leaf)
{
	for (unsigned int c=0; c<nboxes; ++c) child_list[c] = NULL;
}






// child box index
template <int DIM>
static
const unsigned int
getBoxIdx(const typename DimTraits<DIM>::point_type& pos)
{
  for (unsigned int b=0; b<DimTraits<DIM>::nboxes; ++b) {
    bool is_it = true;
    for (unsigned int d=0; d<DIM; ++d)
      if (DimTraits<DIM>::child_pos[b][d] * pos[d] < 0.) is_it = false;
    if (is_it) return b;
  }
  throw std::runtime_error("X is in no box!");
  return std::numeric_limits<unsigned int>::signaling_NaN();
}


// child center
template <int DIM>
static
const typename DimTraits<DIM>::point_type
updateCenter(const typename DimTraits<DIM>::point_type& pos,
						 const double diam,
						 const unsigned int b)
{
  typename DimTraits<DIM>::point_type cpos;
  for (unsigned int d=0; d<DIM; ++d)
    cpos[d] = pos[d] + DimTraits<DIM>::child_pos[b][d] * diam;
  return cpos;
}


// child position
template <int DIM>
static
const typename DimTraits<DIM>::point_int_type
updatePosition(const typename DimTraits<DIM>::point_int_type& pos,
						const unsigned int b)
{
  typename DimTraits<DIM>::point_int_type cpos;
  for (unsigned int d=0; d<DIM; ++d)
    cpos[d]
			= 2 * pos[d] + (DimTraits<DIM>::child_pos[b][d]==1 ? 1 : 0);
  return cpos;
}









// recursive subdivision
template <int DIM>
template <typename particle_type>
void Cluster<DIM>::subdivide
(const std::vector<level_type*>& levels,
 std::vector<std::vector<Cluster<DIM>* > >& clvv,
 Cluster<DIM> *const cl,
 const particle_type *const particles,
 unsigned int *const pindices,
 const unsigned int lmax,
 unsigned int& cluster_counter)
{
  const point_type& center   = cl->getCenter();
	const point_int_type& pos  = cl->getPosition();
  const unsigned int nbeg    = cl->getNbeg();
  const unsigned int nend    = cl->getNend();
  const unsigned int size    = cl->getSize();
  const unsigned int nlevel  = cl->getNlevel();
  const unsigned int cnlevel = nlevel + 1;
  const double diam          = levels.at(nlevel)->getDiam();

  // initialization of child boxes 
  unsigned int* boxes[nboxes];
  unsigned int ndofs_per_box[nboxes];
  for (unsigned int b=0; b<nboxes; ++b) {
    boxes[        b] = new unsigned [size];
    ndofs_per_box[b] = 0;
  }

  // assign particles to respective child box
  for (unsigned int p=nbeg; p<nend; ++p) {
		const particle_type& particle = particles[pindices[p]];
    const unsigned int b = getBoxIdx<DIM>(particle.getPoint() - center);
    boxes[b][ndofs_per_box[b]++] = particle.getId();
  }

	//////////////////////////////////////////////////////////
  // create and subdivide child boxes
  unsigned int dof_counter = 0;
  unsigned int box_counter = 0;
  BOOST_FOREACH(Cluster* &child, cl->child_list) {
    
    // create child only if it contains particles
    assert(child==NULL);
    if (ndofs_per_box[box_counter] > 0) {
      
      // update particle info
      for (unsigned int p=0; p<ndofs_per_box[box_counter]; ++p) {
        const unsigned int idx = nbeg + dof_counter++;
        pindices[idx] = boxes[box_counter][p];
      }
      
      // get child center
      const point_type ccenter
        = updateCenter<DIM>(center, diam/4, box_counter);
      
			// get child position
			const point_int_type cpos
				= updatePosition<DIM>(pos, box_counter);

			const bool is_leaf = (cnlevel<lmax) ? false : true;

      // create child
      child = new Cluster(cl,
                          cluster_counter++,
                          box_counter,
                          clvv.at(cnlevel).size(),
                          nbeg + dof_counter - ndofs_per_box[box_counter],
                          nbeg + dof_counter,
                          ccenter,
													cpos,
                          cnlevel,
													is_leaf);

			// push child into vector of its level (stores clusters sequentially)
      clvv.at(cnlevel).push_back(child);
			
      // subdivide child box if not yet leaf
      if (!is_leaf)
        subdivide(levels, clvv, child, particles, pindices,
									lmax, cluster_counter);
    }
    
    delete [] boxes[box_counter++];
  } ////////////////////////////////////////////////////////
  
  // control
  if (dof_counter != size) {
    std::cerr << "number of particles not coinciding" << std::endl;
    exit(-1);
  }

}




#endif

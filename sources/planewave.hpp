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

#ifndef planewave_hpp
#define planewave_hpp


#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>


template <int> class map_loc_glob;
template <int> struct Chebyshev;
template <int,int> class Tensor;






template <typename T>
class PlaneWaveBase
{
  const T k; //!< laplace parameter

protected:
  explicit PlaneWaveBase(const double wavenum) : k(0., wavenum) {}

public:
	template <typename point_type>
  T plane_wave(const point_type& point, const point_type& u) const
	{	return exp(k * u.dotProduct(point)); } // 7 flops
};








template <int DIM, int ORDER,	typename T>
class PlaneWave
	: boost::noncopyable, public PlaneWaveBase<T>
{
  enum {nboxes = DimTraits<DIM>::nboxes,
        nnodes = BasisTraits<ORDER,DIM>::nnodes};

  typedef Chebyshev< ORDER>  basis_type;
  typedef Tensor<DIM,ORDER> tensor_type; 
	typedef typename DimTraits<DIM>::point_type point_type;


  const T k; //!< laplace parameter
	
	std::vector<T*> pwaves;

public:
	template <typename level_type>
	explicit PlaneWave(const double wavenum, const level_type *const level)
		: PlaneWaveBase<T>(wavenum)
	{
		assert(level->is_in_high_frequency_regime());

		const unsigned int ncones = level->getNCones();
		pwaves.resize(ncones, NULL);
		for (unsigned int c=0; c<ncones; ++c) pwaves.at(c) = NULL;

		typedef std::map<unsigned int,point_type> direction_map;
		const direction_map& existDirs = level->getExistDirs();
		BOOST_FOREACH(const typename direction_map::value_type& dir, existDirs) {

			const unsigned int c = dir.first;
			const point_type&  u = dir.second;

			pwaves.at(c) = new T [(1+nboxes) * nnodes];

			T *const ppw = pwaves.at(c);
			T *const cpw = ppw + nnodes;

			const double diam = level->getDiam();

			// parent cluster plane wave
			
			boost::array<point_type,nnodes> x; //!< parent cluster roots
			map_loc_glob<DIM> map_l2g(point_type(0.), diam);
			for (unsigned int n=0; n<nnodes; ++n) {
				for (unsigned int d=0; d<DIM; ++d)
					x[n][d] = map_l2g(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
				ppw[n] = plane_wave(x[n], u);
			}
			
			// child cluster plane wave
			for (unsigned int b=0; b<nboxes; ++b) {
				// child center
				point_type cc;
				for (unsigned int d=0; d<DIM; ++d)
					cc[d] = DimTraits<DIM>::child_pos[b][d] * (diam / 4.);
				// child roots
				for (unsigned int n=0; n<nnodes; ++n) {
					const point_type cx = (x[n] / 2.) + cc;
					cpw[b*nnodes + n] = plane_wave(cx, u);
				}
			}
			
		}
		
	}


	~PlaneWave()
	{	BOOST_FOREACH(T *const pw, pwaves)	if (pw) delete [] pw;	}

	const T *const getppw(const unsigned int c) const
	{ return pwaves.at(c); }

	const T *const getcpw(const unsigned int c) const
	{ return pwaves.at(c) + nnodes; }

};








template <int DIM, int ORDER,	typename T,
					typename particle_type>
class PlaneWaveParticles
	: boost::noncopyable, public PlaneWaveBase<T>
{
	enum {nnodes = BasisTraits<ORDER,DIM>::nnodes};
	typedef typename DimTraits<DIM>::point_type point_type;

	const particle_type *const particles;
	const unsigned int  *const pindices;


public:
	explicit PlaneWaveParticles(const double wavenum,
															const particle_type *const _particles,
															const unsigned int  *const _pindices)
		: PlaneWaveBase<T>(wavenum),
			particles(_particles),
			pindices( _pindices)
	{}


	void setcpw(const Cluster<DIM> *const cl,	const point_type& u,
							T *const cpw) const
	{
		assert(cl->isLeaf());
		
		const unsigned int nbeg = cl->getNbeg();
		for (unsigned int n=0; n<cl->getSize(); ++n) {
			const point_type& point = particles[pindices[nbeg+n]].getPoint();
			cpw[n] = plane_wave(point, u);
		}
	}


};




#endif

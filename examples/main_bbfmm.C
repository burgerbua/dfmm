// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <utility>
#include <algorithm>
#include <omp.h>


// dfmm
#include "sources/dimordertraits.hpp"

#include "sources/blas.h"
#include "sources/functions.hpp"

#include "sources/mapping.hpp"
#include "sources/chebyshev.hpp"
#include "sources/tensor.hpp"

#include "sources/aca.hpp"

#include "sources/kernelfunctions.hpp"



// y = UV'x where U is of size M times k and V is of size N times k
template <typename T>
void applyUV(const unsigned int M, const unsigned int N, const unsigned int k,
						 const T *const U, const T *const V,
						 const T *const x, T *const y)
{
	T *const temp = new T [k];
	blas::gemhv(N, k, 1., const_cast<T*const>(V), const_cast<T*const>(x), temp);
	blas::gemv( M, k, 1., const_cast<T*const>(U), temp, y);
	delete [] temp;
}





template <typename kernel_type,
					typename point_type>
void testACA(const kernel_type& kernel,
						 const double eps,
						 const unsigned int nx, const point_type *const x,
						 const unsigned int ny, const point_type *const y,
						 const typename kernel_type::value_type *const w)
{
	// typedefs
	typedef typename kernel_type::value_type value_type;

	// entry computer
	EntryComputer<kernel_type,point_type>
		computer(nx, x, ny, y, kernel);

	// direct computation
	value_type *const K = new value_type [nx * ny];
	computer(0, nx, 0, ny, K);
	value_type *const p = new value_type [nx];
	blas::gemv(nx, ny, value_type(1.), K, const_cast<value_type*>(w), p);

	// fully pivoted ACA
	value_type *const R = new value_type [nx * ny];
	blas::copy(nx * ny, K, R);
	const double tf0 = omp_get_wtime();
	value_type *Uf, *Vf;
	unsigned int rankf;
	dfmm::fACA(R, nx, ny, eps, Uf, Vf, rankf);
	value_type *const pf = new value_type [nx];
	applyUV(nx, ny, rankf, Uf, Vf, w, pf);
	const double tf1 = omp_get_wtime();

	// partially pivoted ACA
	const double tp0 = omp_get_wtime();
	value_type *Up, *Vp;
	unsigned int rankp;
	dfmm::pACA(computer, nx, ny, eps, Up, Vp, rankp);
	value_type *const pp = new value_type [nx];
	applyUV(nx, ny, rankp, Up, Vp, w, pp);
	const double tp1 = omp_get_wtime();
	
	std::cout << "- fACA (k="<< rankf
						<<  ", t=" << tf1-tf0 << "s): \tL2 error = "
						<< computeL2norm(nx, p, pf)
						<< "\t Inf error = " << computeINFnorm(nx, p, pf)
						<< std::endl;
	std::cout << "- pACA (k=" << rankp
						<<  ", t=" << tp1-tp0 << "s): \tL2 error = "
						<< computeL2norm(nx, p, pp)
						<< "\t Inf error = " << computeINFnorm(nx, p, pp)
						<< std::endl;


	delete [] K;
	delete [] R;
	delete [] Uf;
	delete [] Vf;
	delete [] Up;
	delete [] Vp;

	delete [] p;
	delete [] pf;
	delete [] pp;
}






int main( int argc, char * argv[ ] )
{
	// start timer
	const double t_0 = omp_get_wtime( );

  const unsigned int DIM   = 3;  
	const unsigned int ORDER = 4;
	const double eps = 1e-4;

	std::cout << "\n= (order,eps) = (" << ORDER << "," << eps << ")"
						<< std::endl;

	const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

	typedef DimTraits<DIM>::point_type point_type;
	typedef Chebyshev< ORDER>          basis_type;
  typedef Tensor<DIM,ORDER>         tensor_type; 


	const double diam = 2.;
	const point_type cx(0., 0., 0.);


	// local to global mapper
	map_loc_glob<DIM> map_l2g_x(cx, diam);

	// target nodes (x)
	point_type x[nnodes];
	for (unsigned int n=0; n<nnodes; ++n)
		for (unsigned int d=0; d<DIM; ++d)
			x[n][d] = map_l2g_x(d, basis_type::nodes[tensor_type::node_ids[n][d]]);

	// source nodes (y)
	unsigned int cells = 0;
	point_type *const y = new point_type [nnodes * 316];
	for (int i=-3; i<=3; ++i)
		for (int j=-3; j<=3; ++j)
			for (int k=-3; k<=3; ++k)
				if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
					// center of Y cluster
					const point_type cy((double)i*diam,	(double)j*diam,	(double)k*diam);
					const map_loc_glob<DIM> map_l2g_y(cy, diam);
					
					for (unsigned int n=0; n<nnodes; ++n)
						for (unsigned int d=0; d<DIM; ++d)
							y[cells*nnodes + n][d]
								= map_l2g_y(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
					cells++;
				}
	
	std::cout << "Number of interactions = " << cells << std::endl;
	

	{
		////////////////////////////////////////////////////////////////////
		std::cout << "\nLaplace kernel -> real valued ACA" << std::endl;
		////////////////////////////////////////////////////////////////////

		typedef KernelFunction<LAPLACE3D> kernel_type;
		typedef kernel_type::value_type    value_type;
		kernel_type kernel;

		// random source vector (for y)
		value_type *const wy = new value_type [nnodes*cells];
		for (unsigned int n=0; n<nnodes*cells; ++n)
			wy[n] = (double)rand()/(double)RAND_MAX;
		testACA(kernel, eps, nnodes, x, nnodes*cells, y, wy);
		delete [] wy;
		
		// random source vector (for x)
		value_type *const wx = new value_type [nnodes];
		for (unsigned int n=0; n<nnodes; ++n)
			wx[n] = (double)rand()/(double)RAND_MAX;
		testACA(kernel, eps, nnodes*cells, y, nnodes, x, wx);
		delete [] wx;
	}
	
		
	{
		////////////////////////////////////////////////////////////////////
		std::cout << "\nHelmholtz kernel -> complex valued ACA" << std::endl;
		////////////////////////////////////////////////////////////////////
		
		// required command line arguments
		if (argc != 2) {
			std::cerr << "Wrong number of command line arguments" << std::endl;
			exit(-1);
		}
		const double wavenum = atof(argv[1]);
		typedef KernelFunction<HELMHOLTZ3D> kernel_type;
		typedef kernel_type::value_type      value_type;
		kernel_type kernel(value_type(0., wavenum));

		// random source vector (for y)
		value_type *const wy = new value_type [nnodes*cells];
		for (unsigned int n=0; n<nnodes*cells; ++n)
			wy[n] = value_type((double)rand()/(double)RAND_MAX,
												 (double)rand()/(double)RAND_MAX);
		testACA(kernel, eps, nnodes, x, nnodes*cells, y, wy);
		delete [] wy;

		// random source vector (for x)
		value_type *const wx = new value_type [nnodes];
		for (unsigned int n=0; n<nnodes; ++n)
			wx[n] = value_type((double)rand()/(double)RAND_MAX,
												 (double)rand()/(double)RAND_MAX);
		testACA(kernel, eps, nnodes*cells, y, nnodes, x, wx);
		delete [] wx;
	}


	delete [] y;

	////////////////////////////////////////////////////
	std::cout << "\n- Overall time " << omp_get_wtime( ) - t_0
            << " s" << std::endl;
	////////////////////////////////////////////////////


	return 0;
} 	
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

	particle_type *const particles = new particle_type [N];
	unsigned int  *const pindices  = new unsigned int  [N];
  ReadParticles(filename, N, particles);

	// fill pindices
	for (unsigned int p=0; p<N; ++p)
		pindices[p] = particles[p].getId();

  // plot particles
  //Plot<particle_type> plotter(particles, N);
  //plotter.plot();

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



  //typedef M2LHandlerSArcmp<dim,order,value_type> m2l_handler_type;
  //typedef M2LHandlerNA<dim,order,value_type> m2l_handler_type;
	//typedef M2LHandlerNAsym<dim,order,value_type> m2l_handler_type;
	//typedef M2LHandlerNAblk<dim,order,value_type> m2l_handler_type;
  //typedef M2LHandlerIA<dim,order,value_type> m2l_handler_type;
  //typedef M2LHandlerIAsym<dim,order,value_type> m2l_handler_type;
	typedef M2LHandlerIAblk<dim,order,value_type> m2l_handler_type;


  // init m2l handler 
  typedef std::vector<m2l_handler_type> m2l_handler_type_vector;
  m2l_handler_type_vector all_m2l_handler;
  storage.initialize(all_m2l_handler);
	m2l_handler_type::getInfo();




  // cluster bases
  typedef ClusterBasis<dim,order,value_type> clusterbasis_type;
  std::vector<clusterbasis_type> rbasis(ncl), cbasis(ncl);

  // cluster expansions
  typedef ExpansionHandler<dim,order,value_type> expansion_type;
  std::vector<expansion_type> rexph(ncl), cexph(ncl);
  
  // cluster relations
  typedef Relations<dim> relations_type;
  std::vector<relations_type> relations(ncl);

  // setup multilevel relations
  SetupClusterRelations(all_cluster_vectors, rexph,
                        all_cluster_vectors, cexph, 
                        relations, storage.getLevels(), all_m2l_handler);


  // write out info
  storage.writeInfo(all_m2l_handler);


	////////////////////////////////////////////////////////////////////
  // draw cross section through cluster oct-tree at "loc" in direction "dir"
	point_type loc;
	for (unsigned int i=0; i<dim; ++i)
		loc[i] = (bbx.first[i]+bbx.second[i]) / 2;
	std::cout << "\n- Center of particles at " << loc << std::endl;
	std::cout << " - a = " << bbx.first << "\tb = " << bbx.second << std::endl;
	
	//	std::cout << "\n- Origin of PS pictures at " << loc << std::endl;
	
	//const unsigned int dir0 = 0;
	//drawCrossSection(filename, loc, dir0, 
	//								 storage.getRootClusterCenter(),
	//								 storage.getRootClusterExtension(),
	//								 all_cluster_vectors, storage.getLevels(), all_m2l_handler);
	//
	//const unsigned int dir1 = 1;
	//drawCrossSection(filename, loc, dir1, 
	//								 storage.getRootClusterCenter(),
	//								 storage.getRootClusterExtension(),
	//								 all_cluster_vectors, storage.getLevels(), all_m2l_handler);
	//
	//const unsigned int dir2 = 2;
	//drawCrossSection(filename, loc, dir2, 
	//								 storage.getRootClusterCenter(),
	//								 storage.getRootClusterExtension(),
	//								 all_cluster_vectors, storage.getLevels(), all_m2l_handler);
	//////////////////////////////////////////////////////////////////////



  // compute m2l operators
  ComputeM2Loperators(all_m2l_handler, kernel, eps);



  // l2l/m2m operators
  std::cout << "\n- Evaluating interpolation operators " << std::flush;
  const double t_io = omp_get_wtime();
  setIO(all_cluster_vectors, rbasis, all_m2l_handler,	pindices, particles);
  setIO(all_cluster_vectors, cbasis, all_m2l_handler,	pindices, particles);
  std::cout << " took " << omp_get_wtime() - t_io << " s" << std::endl;



  // densities and potentials
  value_type *const y = new value_type [N];
  for (unsigned int n=0; n<N; ++n)
		y[n] = (double)rand()/(double)RAND_MAX;
  value_type *const x = new value_type [N];
  for (unsigned int n=0; n<N; ++n) x[n] = 0.;





  // m2m
  std::cout << "\n- Applying M2M operators\n" << std::flush;
  const double t_m2m = omp_get_wtime();
  ApplyM2Moperators(all_cluster_vectors, cbasis, cexph, all_m2l_handler, y);
  std::cout << "  - overall " << omp_get_wtime() - t_m2m << " s" << std::endl;

  // m2l
  std::cout << "\n- Applying M2L operators\n" << std::flush;
  const double t_m2l = omp_get_wtime();
  ApplyM2Loperators(all_cluster_vectors, rexph, cexph, relations,
                    all_m2l_handler);
  std::cout << "  - overall " << omp_get_wtime() - t_m2l << " s" << std::endl;

  // l2l
  std::cout << "\n- Applying L2L operators\n" << std::flush;
  const double t_l2l = omp_get_wtime();
  ApplyL2Loperators(all_cluster_vectors, rbasis, rexph, all_m2l_handler, x);
  std::cout << "  - overall " << omp_get_wtime() - t_l2l << " s" << std::endl;

  // apply near field on the fly
  typedef NFadder<cluster_type,particle_type,relations_type,kernel_type>
		nfc_type;
  ApplyNearField(all_cluster_vectors.back(),
								 nfc_type(kernel, relations,
													pindices, particles,
													pindices, particles,
													x, y));



//	std::cout << "\n- Size of value_type is " << sizeof(value_type)
//						<< " byte." << std::endl;
//
//	// direct calculation
//	const long size = ndofs*ndofs;
//	const long size_in_bytes = size * sizeof(value_type);
//	const double mem_dns = (long)size_in_bytes / (1024.*1024.);
//  std::cout << "\n- Direct calculation of " << ndofs << " times " << ndofs
//						<< " particles requires " 
//						<< (long)size_in_bytes / (1024.*1024.)
//						<< " Mb."<< std::endl;
//	const double t_direct = omp_get_wtime();
//  value_type *const x0 = new value_type [ndofs];
//	blas::setzero(ndofs, x0);
//	std::cout << " - Memory allocation " << std::flush;
//	const double t_alloc = omp_get_wtime();
//  value_type *const K = new value_type [ndofs * ndofs];
//	const double t_calc = omp_get_wtime();
//  std::cout << "took " << t_calc - t_alloc << " s" << std::endl;
//	std::cout << " - Matrix evaluation " << std::flush;
//  storage.getFundSol().compute(ndofs, dofs, ndofs, dofs, K);
//	const double t_matvec = omp_get_wtime();
//  std::cout << "took " << t_matvec - t_calc << " s" << std::endl;
//  // matrix vector product
//	std::cout << " - Matrix vector multiplication " << std::flush;
//  blas::gemva(ndofs, ndofs, 1., K, y, x0);
//  std::cout << "took " << omp_get_wtime() - t_matvec << " s" << std::endl;
//  std::cout << "  took " << omp_get_wtime() - t_direct << " s" << std::endl;


	
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


//	value_type sum_x  = 0;
//	value_type sum_x0 = 0;
//  for (unsigned int n=0; n<nx; ++n) {
//		sum_x  += x[ n];
//		sum_x0 += x0[n];
//	}
//	std::cout << "\nSize of cluster X = " << nx
//						<< "\tsumx = " << sum_x
//						<< "\tsumx0 = " << sum_x0 << std::endl;
//  for (unsigned int n=0; n<(10<nx?10:nx); ++n)
//    std::cout << x[n] << "\t" << x0[n] << "\t" << x[n]-x0[n] << std::endl;
  
  std::cout << "\n- L2  error " << computeL2norm( nx, x0, x) << std::endl;
  std::cout <<   "- Inf error " << computeINFnorm(nx, x0, x) << std::endl;

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

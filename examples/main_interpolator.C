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


// qbdfmm
#include "sources/dimordertraits.hpp"

#include "sources/blas.h"
#include "sources/functions.hpp"

#include "sources/mapping.hpp"
#include "sources/chebyshev.hpp"
#include "sources/tensor.hpp"

#include "sources/interpolator.hpp"
#include "sources/aca.hpp"

#include "sources/kernelfunctions.hpp"




// give random double number r with max(abs(r)) = radius
double random(const double radius)
{
	const double val = ((double)rand()/(double)RAND_MAX - (double).5) * radius;
	return val;
}


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



int main( int argc, char * argv[ ] )
{
	// start timer
	const double t_0 = omp_get_wtime( );

  const unsigned int DIM   = 3;  
	const unsigned int ORDER = 7;

	std::cout << "\nOrder = " << ORDER << std::endl;

	const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

	typedef DimTraits<DIM>::point_type point_type;
	typedef Chebyshev< ORDER>          basis_type;
  typedef Tensor<DIM,ORDER>         tensor_type; 


	const double diam = 2.;
	const point_type cx(0.,      0., 0.);
	const point_type cy(2.*diam, 0., 0.);
	const unsigned int M = 10000;
	const unsigned int N = 10000;

	// local to global mapper
	map_loc_glob<DIM> map_l2g_x(cx, diam);
	map_loc_glob<DIM> map_l2g_y(cy, diam);

	// interpolation points and particles
	point_type X[nnodes], Y[nnodes];
	point_type *const x = new point_type [M];
	point_type *const y = new point_type [N];
	for (unsigned int d=0; d<DIM; ++d) {
		for (unsigned int n=0; n<nnodes; ++n) {
			X[n][d] = map_l2g_x(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
			Y[n][d] = map_l2g_y(d, basis_type::nodes[tensor_type::node_ids[n][d]]);
		}
		for (unsigned int n=0; n<M; ++n) x[n][d] = map_l2g_x(d, random(2.));
		for (unsigned int n=0; n<N; ++n) y[n][d] = map_l2g_y(d, random(2.));
	}
	
	// target and random source vector
	double *const p = new double [M];
	blas::setzero(M, p);
	double *const w = new double [N];
	for (unsigned int n=0; n<N; ++n) w[n] = (double)rand()/(double)RAND_MAX;

	// define kernel function
	typedef KernelFunction<LAPLACE3D> kernel_type;
	kernel_type kernel;
	EntryComputer<kernel_type,point_type>	computer(M, x, N, y, kernel);

	// direct computation
	double *const K = new double [M * N];
	computer(0, M, 0, N, K);
	blas::gemv(M, N, 1., K, w, p);
	
  // computer P2M and L2P operators
	IOcomputerLF<DIM,ORDER,double> iocx(cx, diam);
	IOcomputerLF<DIM,ORDER,double> iocy(cy, diam);
	double *const Sx = new double [M * nnodes];
	double *const Sy = new double [N * nnodes];
	const double t0 = omp_get_wtime();
	iocx(M, x, Sx);
	const double t1 = omp_get_wtime();
	iocy(N, y, Sy);
	const double t2 = omp_get_wtime();
	std::cout << "\nL2P computed in "    << t1-t0
						<< "s\nP2M computed in " << t2-t1 << "s" << std::endl;
	
	// compute M2L operator
	EntryComputer<kernel_type,point_type>	m2lc(nnodes, X, nnodes, Y, kernel);
	double *const R = new double [nnodes * nnodes];
	m2lc(0, nnodes, 0, nnodes, R);

	//
	double W[nnodes];
	blas::gemhv(N, nnodes, 1., Sy, w, W);
	double P[nnodes];
	blas::gemv(nnodes, nnodes, 1., R, W, P);
	double *const pp = new double [M];
	blas::gemv( M, nnodes, 1., Sx, P, pp);
	

	std::cout << "\nL2  error = "   << computeL2norm( M, p, pp)
						<< "\nInf error = " << computeINFnorm(M, p, pp)
						<< std::endl;


	delete [] Sx;
	delete [] Sy;
	delete [] R;
	delete [] pp;

	delete [] K;
	
	delete [] w;
	delete [] p;
	
	delete [] x;
	delete [] y;



	////////////////////////////////////////////////////
	std::cout << "\n- Overall time " << omp_get_wtime( ) - t_0
            << " s" << std::endl;
	////////////////////////////////////////////////////


	return 0;
} 	

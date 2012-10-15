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
#include <iostream>
#include <omp.h>

#include "sources/dimordertraits.hpp"
#include "sources/symmetries.hpp"

void write(const unsigned int nnodes,
					 const unsigned int *const P)
{
	for (unsigned int i=0; i<nnodes; ++i) {
		for (unsigned int j=0; j<nnodes; ++j)
			std::cout << P[j*nnodes + i] << " ";
		std::cout << std::endl;
	}
}

void p2P(const unsigned int nnodes,
				 const unsigned int *const p,
				 unsigned int *const P)
{
	for (unsigned int i=0; i<nnodes; ++i)
		for (unsigned int j=0; j<nnodes; ++j)
			if (p[i] == j) P[j*nnodes + i] = 1;
			else           P[j*nnodes + i] = 0;
}

void P2p(const unsigned int nnodes,
				 const unsigned int *const P,
				 unsigned int *const p)
{
	for (unsigned int i=0; i<nnodes; ++i)
		for (unsigned int j=0; j<nnodes; ++j)
			if (P[j*nnodes + i] == 1) p[i] = j;
}

// P = P1 * P2
void mult(const unsigned int nnodes,
					const unsigned int *const P1,
					const unsigned int *const P2,
					unsigned int *const P)
{
	for (unsigned int i=0; i<nnodes; ++i)
		for (unsigned int j=0; j<nnodes; ++j) {
			P[j*nnodes + i] = 0;
			for (unsigned int n=0; n<nnodes; ++n)
				P[j*nnodes + i] += P1[n*nnodes + i] * P2[j*nnodes + n];
		}
}



int main( int argc, char * argv[ ] )
{
  const unsigned int dim = 3;  
	const unsigned int order = 3;

	const unsigned int nnodes = BasisTraits<order,dim>::nnodes;
	
	typedef DimTraits<dim>::point_int_type point_int_type;

	Symmetries<order,dim> sym;

	point_int_type ip, op, oop;
	std::cout << "X - coordinate = "; std::cin >> ip[0];
	std::cout << "Y - coordinate = "; std::cin >> ip[1];
	std::cout << "Z - coordinate = "; std::cin >> ip[2];
	std::cout << "t = " << ip << std::endl;

	const unsigned int* perm;

	perm = sym.getPermutation(ip, op);
	std::cout << "p(t) = " << op << std::endl;
	for (unsigned int i=0; i<nnodes; ++i)	std::cout << perm[i] << " ";
	std::cout << std::endl;
	
	unsigned int *const P = new unsigned int [nnodes * nnodes];
	p2P(nnodes, perm, P);
	//write(nnodes, P);

	perm = sym.getOctant(ip, oop);
	//std::cout << "p(t) = " << oop << std::endl;
	//for (unsigned int i=0; i<nnodes; ++i)	std::cout << perm[i] << " ";
	//std::cout << std::endl;

	unsigned int *const P1 = new unsigned int [nnodes * nnodes];
	p2P(nnodes, perm, P1);
	//write(nnodes, P1);

	perm = sym.getPermutation(oop, op);
	std::cout << "p(t) = " << op << std::endl;
	//for (unsigned int i=0; i<nnodes; ++i)	std::cout << perm[i] << " ";
	//std::cout << std::endl;

	unsigned int *const P2 = new unsigned int [nnodes * nnodes];
	p2P(nnodes, perm, P2);
	//write(nnodes, P2);
	mult(nnodes, P1, P2, P);
	//write(nnodes, P);

	unsigned int *const p = new unsigned int [nnodes];
	P2p(nnodes, P, p);
	for (unsigned int i=0; i<nnodes; ++i)	std::cout << p[i] << " ";
	std::cout << std::endl;
	

	return 0;
} 	
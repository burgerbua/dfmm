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

#ifndef functions_hpp
#define functions_hpp

#include <stdexcept>

#include "sources/blas.h"


// relative point-wise L2 error
template <typename T>
const double computeL2norm(unsigned int N, T *const u, T *const v)
{
	T *const w = new T [N];
	blas::copy(N, v, w);
	blas::axpy(N, -1., u, w);
	double      dot = blas::dot(N, u);
	double diff_dot = blas::dot(N, w);
	delete [] w;
	return sqrt(diff_dot / dot);
}


// relative point-wise inf error
template <typename T>
const double computeINFnorm(unsigned int N, T *const u, T *const v)
{
	double      max = 0.;
	double diff_max = 0.;
	for (unsigned int n=0; n<N; ++n) {
		if (     max<std::abs(u[n]))           max = std::abs(u[n]);
		if (diff_max<std::abs(u[n]-v[n])) diff_max = std::abs(u[n]-v[n]);
	}
	return diff_max / max;
}


// check if \a data has nan entries
template <typename T>
bool check_nan(const unsigned int n, const T *const data)
{
  for (unsigned int i=0; i<n; ++i)
		if (data[i]!=data[i]) {
			std::cout << i << " of " << n << " gives " << data[i] << std::endl;
			return false;
		}
  return true;
}


// determine rank singular-values cutoff based on euclidian error with eps
unsigned int getRank(const double *const singular_values,
										 const unsigned int size,
										 const double eps)
{
	const double nrm2 = blas::scpr(size, singular_values, singular_values);
	double nrm2k(0.);
	for (unsigned int k=size; k>0; --k) {
		nrm2k += singular_values[k-1] * singular_values[k-1];
		if (nrm2k > eps*eps * nrm2)	return k;
	}
	throw std::runtime_error("rank cannot be larger than size");
	return 0;
}





#ifdef DFMM_COUNT_FLOPS
class FlopCounter
{
	std::string description;
	unsigned int nlevel;
	unsigned long long flops;
public:
	FlopCounter(const FlopCounter& other)
	{
		description = other.description;
		nlevel      = other.nlevel;
		flops = 0;
	}
	

	FlopCounter& operator= (const FlopCounter& other)
	{
		description = other.description;
		nlevel      = other.nlevel;
		flops = 0;
		return *this;
	}

	explicit FlopCounter(const std::string _description,
											 const unsigned int _nlevel)
		: description(_description), nlevel(_nlevel), flops(0)
	{}
	
	~FlopCounter()
	{
		if (flops)
			std:: cout << " (" << description << "[" << nlevel
								 << "] flops = " << flops << ") " << std::endl;
	}
	
	void addFlops(const unsigned int _flops)
	{	flops += _flops; }

};
#else
class FlopCounter
{
public:
	explicit FlopCounter(const std::string&,
											 const unsigned int) {}
	void addFlops(const unsigned int) {}
};
#endif













//#include "sources/tensor.hpp"
//
//#include <boost/unordered_set.hpp>
//#include <boost/unordered_map.hpp>
//
//#include <numeric>


//// forward declaration
//class dof;
//template <int> class Dof;



//template <typename T>
//void write(const unsigned int n, const T *const data)
//{
//  for (unsigned int i=0; i<n; ++i)
//    std::cout << data[i] << std::endl;
//}




//template <typename T>
//struct add_square : public std::binary_function<T,T,T>
//{
//  T operator()(T result, T value) const
//  { return result + value*value; }
//};



//template <typename T>
//const unsigned int getRank(const T *const singular_values,
//                           const unsigned int n,
//                           const double eps)
//{
//  const double sum
//    = std::accumulate(singular_values, singular_values+n,
//                      0., add_square<double>());
//  
//  unsigned int k = 0;
//  double sum_k = sum;
//  while (sqrt(sum_k) > eps*sqrt(sum))
//    sum_k -= singular_values[k]*singular_values[k++];
//  return k;
//}




///**
// * Generate a sequential cluster list
// */
//template <typename CL>
//void generateSequentialClusterList(CL *const cl,
//                                   std::vector<CL*>& clList)
//{
//  clList.push_back(cl);
//  for (unsigned int i=0; i<CL::nboxes; ++i) {
//    CL *const child = cl->getChild(i);
//    if (child) generateSequentialClusterList(child, clList);
//  }
//}





///**
// * Generate a sequential leaf cluster list
// */
//template <typename CL>
//void generateSequentialLeafClusterList(CL *const cl,
//                                       std::vector<CL*>& clList)
//{
//  if (cl->isleaf()) clList.push_back(cl);
//  for (typename CL::child_iterator it=cl->beginChild();
//       it!=cl->endChild(); ++it)
//    if (*it) generateSequentialLeafClusterList(static_cast<CL*>(*it), clList);
//}


///**
// * Generate a level based sequential cluster list 
// */
//template <typename CL>
//void generateLevelBasedClusterLists(CL *const cl,
//                                    std::vector<std::vector<CL*> >& clist)
//{
//  clist.at(cl->getLevel()->getNLevel()).push_back(cl);
//  for (typename CL::child_iterator it=cl->beginChild();
//       it!=cl->endChild(); ++it)
//    if (*it) generateLevelBasedClusterLists(static_cast<CL*>(*it), clist);
//}
















//template <typename CL,
//          typename KERNEL>
//void evaluate(const KERNEL& kernel,
//							CL* clX, CL* clY,
//							typename KERNEL::value_type *const f,
//							typename KERNEL::value_type *const &p)
//{
//	typedef typename KERNEL::value_type value_type;
//	const unsigned int Nx = clX->getN();
//	const unsigned int Ny = clY->getN();
//	value_type* K = new value_type [Nx*Ny];
//	kernel.evaluate_nondirectional(Nx, (const dof**)clX->getParticles(),
//																 Ny, (const dof**)clY->getParticles(), K);
//	blas::gemv(Nx, Ny, 1., K, f, p);
//	delete [] K;
//}






//template <typename CL,
//          typename KERNEL>
//void approximate(KERNEL& kernel,
//								 typename DimTraits<CL::dim>::point_type direction,
//								 CL* clX, CL* clY,
//								 typename KERNEL::value_type *const f,
//								 typename KERNEL::value_type *const &p)
//{
//	typedef typename KERNEL::base_type   base_type;
//	typedef typename KERNEL::value_type value_type;
//	enum {nnodes = BasisTraits<CL::order,CL::dim>::nnodes};
//	const unsigned int Nx = clX->getN();
//	const unsigned int Ny = clY->getN();
//
//	kernel.upward_directional(Ny, (const dof**)clY->getParticles(),
//														direction, f);
//
//	value_type F[nnodes];
//	blas::gemhv(Ny, nnodes, 1., clY->getIAoperator(), f, F);
//
//	unsigned int k;
//	value_type *U, *V;
//	approximate<value_type>(&base_type::evaluate_directional,	&kernel,
//													nnodes, clX->getBasisParticles(),
//													nnodes, clY->getBasisParticles(),
//													k, U, V);
//	value_type P[nnodes];
//	blas::setzero(nnodes, P);
//	for (unsigned int l=0; l<k; ++l) {
//		const value_type e = blas::scpr(nnodes, V+l*nnodes, F);
//		blas::axpy(nnodes, e, U+l*nnodes, P);
//	}
//	delete [] U;
//	delete [] V;
//
//	blas::gemv(Nx, nnodes, 1., clX->getIAoperator(), P, p);
//
//	kernel.downward_directional(Nx, (const dof**)clX->getParticles(),
//															direction, p);
//}








//template <typename T>
//bool checkIAoperator(const T *const f, const unsigned int nf,
//										 const T *const F, const unsigned int nF)
//{
//	T sum_f = 0., sum_F = 0.;
//	for (unsigned int n=0; n<nf; ++n) sum_f += f[n];
//	for (unsigned int n=0; n<nF; ++n) sum_F += F[n];
//	return (std::abs(sum_f-sum_F)	<= 1e-10);
//}











//template <typename superelement_type>
//class SuperelementExtractor
//{
//  typedef typename superelement_type::gdof_type gdof_type;
//  typedef typename superelement_type::ldof_type ldof_type;
//
//
//public:
//  void operator()(const unsigned int ndofs, dof** dofs,
//                  std::vector<const superelement_type*>& elements) const
//  {
//    boost::unordered_set<unsigned int> eids;
//    eids.clear(); 
//    elements.clear();
//    elements.reserve(ndofs);
//    
//    for (unsigned int i=0; i<ndofs; ++i) {
//      const gdof_type* gd = static_cast<const ldof_type*>(dofs[i])->getGDof();
//      for (unsigned int e=0; e<gd->getNumberSupportElements(); ++e) {
//        const superelement_type* el = gd->getSuperElement(e);
//        if ((eids.insert(el->getID())).second) elements.push_back(el);
//      }
//    }
//
//  }
//
//};



//template <typename cluster_type,
//          typename superelement_type>
//class ExtensionEvaluator
//  : public std::unary_function<cluster_type, void>
//{
//  enum {dim       = cluster_type::dim,
//        num_nodes = superelement_type::element_type::num_nodes};
//  typedef typename DimTraits<dim>::point_type point_type;
//  
//
//  const SuperelementExtractor<superelement_type> extractor;
//  dof** dofs;
//
//  double& extension;
//  
//public:
//  explicit ExtensionEvaluator(dof** _dofs, double& _extension)
//    : extractor(), dofs(_dofs), extension(_extension)
//  { extension = 0.; }
//  
//  void operator()(cluster_type *const cl)
//  {
//    assert(cl->isleaf());
//
//    // extract superelements
//    std::vector<const superelement_type*> elements;
//    extractor(cl->size, dofs+cl->nbeg, elements);
//
//    BOOST_FOREACH(const superelement_type *const el, elements) {
//      for (unsigned int k=0; k<num_nodes; ++k) {
//        const point_type& point = *(el->getNode(k));
//        for (unsigned int d=0; d<dim; ++d) {
//          const double ext = 2. * std::abs(point[d] - cl->center[d]);
//          if (ext > extension) extension = ext;
//        }
//      }
//    }
//
//  }
//
//};




//template <int> class Cluster;
//template <int> class map_loc_glob;
//
//template <int DIM, int ORDER>
//void getInterpolationNodes(const Cluster<DIM> *const cl,
//                           const dof* X[BasisTraits<ORDER,DIM>::nnodes])
//{
//  // compile time info
//  typedef typename DimTraits<DIM>::point_type point_type;
//  typedef Chebyshev< ORDER>  basis_type;
//  typedef Tensor<DIM,ORDER> tensor_type; 
//  enum {nnodes = BasisTraits<ORDER,DIM>::nnodes};
//
//  // local to global mapper
//  map_loc_glob<DIM> map_l2g(cl->getCenter(),
//                            cl->getLevel()->getExtensionBoundingBox());
//  
//  // coordinates
//  boost::array<point_type, nnodes> x;
//  for (unsigned int n=0; n<nnodes; ++n)
//    for (unsigned int d=0; d<DIM; ++d)
//      x[n][d] = map_l2g(d,basis_type::nodes[tensor_type::node_ids[n][d]]);
//  
//  // X (field points)
//  for (unsigned int n=0; n<nnodes; ++n)
//    X[n] = new Dof<DIM>(x[n]);
//}






#endif

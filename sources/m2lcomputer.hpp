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

#ifndef m2lcomputer_hpp
#define m2lcomputer_hpp


#include "boost/tuple/tuple.hpp"


#include "sources/blas.h"
#include "sources/functions.hpp"
#include "sources/aca.hpp"


/*! \def DFMM_USE_CHEB_WEIGHTS
	
	If \a DFMM_USE_CHEB_WEIGHTS is defined Chebyshev weights are used
 */
#define DFMM_USE_CHEB_WEIGHTS false


/*! \ingroup mm2l
	
	\name Low-rank approximation schemes
	
	\attention {DEFINE ONLY ONE SCHEME}
 */
/* @{ */
#define DFMM_USE_TSVD false //<! \def DFMM_USE_TSVD Use truncated SVD
#define DFMM_USE_FACA false //<! \def DFMM_USE_FACA Use fully pivoted ACA
#define DFMM_USE_PACA true  //<! \def DFMM_USE_PACA Use partially pivoted ACA

#if DFMM_USE_PACA || DFMM_USE_FACA
#define DFMM_USE_QRTSVD true //<! \def DFMM_USE_QRTSVD Use QR & tSVD after ACA
#endif
/* @} */


// only one low-rank approximation scheme is permitted
#if (DFMM_USE_TSVD + DFMM_USE_FACA + DFMM_USE_PACA) != 1
NOT ALLOWED
#endif






/*! \addtogroup mm2l

	\par M2L-computer
	The M2L-computers are responsible for computing and assembling all M2L
	operators (either dense of by means of approximation schemes \ref mapprox)
	
 */




namespace dfmm {
	
	//! compute x += UV' y	
	template <typename T>
  static void applyUV(const unsigned int k,
                      const unsigned int nx,
                      const T *const U,
                      const unsigned int ny,
                      const T *const V,
                      const T *const y,
                      T *const x)
  {
    for (unsigned int l=0; l<k; ++l) {
      const T e = blas::scpr(ny, V + l*ny, y);
      blas::axpy(nx, e, U + l*nx, x);
    }
  }




	/*! \ingroup mapprox

		\brief Truncated SVD

		The truncated singular value decomposition (tSVD) finds the optimal
		low-rank approximation of a matrix \f$K \sim UV^*\f$ such that
		\f$|K-UV^*|_{L_2} \le \varepsilon |K|_{L_2}\f$ with the prescribed
		accuracy \f$\varepsilon\f$ is satisfied.

		@tparam T value type (either real or complex single or double precision)

		@param[in] nx number of rows of the matrix
		@param[in] ny number of cols of the matrix
		@param[in] K matrix (destroyed on output)
		@param[in] eps prescribed accuracy
		@param[out] U matrix of cols (\a nx times \a k)
		@param[out] V matrix of rows (\a ny times \a k)
		@param[out] k rank of low-rank approximation
	*/
	template <typename T>
	void tSVD(const unsigned int nx,
						const unsigned int ny,
						const T *const K,
						const double eps,
						T * &U,
						T * &V,
						unsigned int &k)
	{
		unsigned int INFO;
		const unsigned int nmax = nx > ny ? nx : ny;
		const unsigned int nmin = nx < ny ? nx : ny;
		const unsigned int n1 = 3*nmin + nmax;
		const unsigned int n2 = 5*nmin;
		const unsigned int LWORK = 2 * (n1>n2?n1:n2);
		T *const WORK = new T [LWORK];
	
		U = new T [nx*nx];
		V = new T [ny*ny];
		double *const S  = new double [nmin];
		T *const VT = new T [ny*ny];
		blas::copy(nx*ny, const_cast<T *const>(K), U);
		INFO = blas::gesvd(nx, ny, U, S, VT, ny, LWORK, WORK);
		if (INFO!=0) throw std::runtime_error("SVD did not converge with " + INFO);
		k = getRank(S, nmin, eps);
	
		// U <- U * S
		for (unsigned int r=0; r<k; ++r)	blas::scal(nx, S[r], U + r*nx);
		// V <- hermitian transpose(VT)
		blas::transpose(ny, ny, VT, V);

		delete [] WORK;
		delete [] S;
		delete [] VT;
	}



	/*! \ingroup mapprox
		
		\brief ACA recompression via subsequent QR-decomposition and truncated SVD

		This function is used if a (non optimal) low-rank representation
		satisfying \f$|K-UV^*|_{L_2} \le \varepsilon |K|_{L_2}\f$ of rank \a k is
		existing (\a U and \a V are non-unitary; \a U, \a V and \a k on the
		input). This function computes QR - decompositions and with a subsequent
		SVD we obtain the optimal low-rank representation (\a U, \a V and \a k on
		the output).

		@tparam T value type (either real or complex single or double precision)

		@param[in] nx number of rows of the matrix
		@param[in] ny number of cols of the matrix
		@param[in] eps prescribed accuracy
		@param[in,out] U matrix of cols (\a nx times \a k)
		@param[in,out] V matrix of rows (\a ny times \a k)
		@param[in,out] k of low-rank approximation
	*/
	template <typename T>
	void QRtSVD(const unsigned int nx,
							const unsigned int ny,
							const double eps,
							T *const U,
							T *const V,
							unsigned int &k)
	{
		unsigned int INFO;
		const unsigned int nmax = nx > ny ? nx : ny;
		const unsigned int nmin = nx < ny ? nx : ny;
		const unsigned int n1 = 3*nmin + nmax;
		const unsigned int n2 = 5*nmin;
		const unsigned int LWORK = 2 * (n1>n2?n1:n2);
		T *const WORK = new T [LWORK];

		// some copying
		T *const QU = new T [nx*k]; blas::copy(nx*k, U, QU);
		T *const QV = new T [ny*k]; blas::copy(ny*k, V, QV);

		// QR decomposition
		T *const phi = new T [k*k];
		{
			// QR of U and V
			T* tauU = new T [k];
			INFO = blas::geqrf(nx, k, QU, tauU, LWORK, WORK);
			assert(INFO==0);
			T* tauV = new T [k];
			INFO = blas::geqrf(ny, k, QV, tauV, LWORK, WORK);
			assert(INFO==0);
			// phi = Ru Rv'
			T* RU = new T [2 * k*k];
			T* RV = RU + k*k;
			blas::setzero(2 * k*k, RU);
			for (unsigned int l=0; l<k; ++l) {
				blas::copy(l+1, QU + l*nx, RU + l*k);
				blas::copy(l+1, QV + l*ny, RV + l*k);
			}
			blas::gemmh(k, k, k, T(1.), RU, k, RV, k, phi, k);
			delete [] RU;
			// get Qu and Qv
			INFO = blas::orgqr(nx, k, QU, tauU, LWORK, WORK);
			assert(INFO==0);
			INFO = blas::orgqr(ny, k, QV, tauV, LWORK, WORK);
			assert(INFO==0);
			delete [] tauU;
			delete [] tauV;
		}
		
		// store low-rank as low-rank obtained via ACA
		const unsigned int aca_k = k;
		
		// further compression
		{
			// SVD
			T *const VT  = new T [k*k];
			double *const S = new double [k];
			INFO = blas::gesvd(aca_k, aca_k, phi, S, VT, aca_k, LWORK, WORK);
			if (INFO!=0)
				throw std::runtime_error("SVD did not converge with " + INFO);
			k = getRank(S, aca_k, eps);

			// (U Sigma)
			for (unsigned int r=0; r<k; ++r)
				blas::scal(aca_k, S[r], phi + r*aca_k);
			delete [] S;
			
			// Qu (U Sigma) 
			blas::gemm(nx, aca_k, k, T(1.), QU, nx, phi, aca_k, U, nx);
			delete [] phi;
			
			// V <- hermitian transpose(VT)
			T *const VV = new T [aca_k * aca_k];
			blas::transpose(aca_k, aca_k, VT, VV);
			delete [] VT;

			// Qv V
			blas::gemm(ny, aca_k, k, T(1.), QV, ny, VV, aca_k, V, ny);
			delete [] VV;
		}

		// free memory
		delete [] QU;
		delete [] QV;
		delete [] WORK;
		
	}







	/*! \ingroup mapprox
 
		\brief Compresses \f$K_i \sim Q_u * C_i * Q_b^*\f$
	
		@tparam T value type
		
		@param[in] nnodes number of interpolation points
		@param[in] ninteractions number of interactions
		@param[in] eps approximation accuracy
		@param[out] Qu matrix of size \a nnodes times \a k
		@param[in,out] C input: all K matrices (as row); output: all C matrices
		@param[out] Qb matrix of size \a nnodes time \a k
		@param[out] k low-rank and each \a C is of size \a k times \a k
	*/
	template <typename T>
	void UCiB(const unsigned int nnodes,
						const unsigned int ninteractions,
						const double eps,
						T* &Qu, T* &C, T* &Qb,
						unsigned int& k)
	{
		// init SVD
		const unsigned int LWORK = 2 * (3*nnodes + ninteractions*nnodes);
		T *const WORK = new T [LWORK];

		// K_col ///////////////////////////////////////////////////////////
		T *const K_col = new T [ninteractions * nnodes*nnodes]; 
		for (unsigned int i=0; i<ninteractions; ++i)
			for (unsigned int j=0; j<nnodes; ++j)
				blas::copy(nnodes,
									 C     + i*nnodes*nnodes + j*nnodes,
									 K_col + j*ninteractions*nnodes + i*nnodes);
		// singular value decomposition
		T *const Q = new T [nnodes*nnodes];
		double *const S = new double [nnodes];
		const unsigned int info_col
			= blas::gesvd(ninteractions*nnodes, nnodes, K_col, S, Q, nnodes,
										LWORK, WORK);
		if (info_col!=0)
			throw std::runtime_error("SVD did not converge with " + info_col);
		delete [] K_col;
		const unsigned int k_col = getRank(S, nnodes, eps);


//		//////////////////////////////////////////////////////////
//		std::ofstream sing_vals;
//		sing_vals.open("e2_16_singular_values.txt",
//									 std::ios::out | std::ios::app);
//		for (unsigned int i=0; i<nnodes; ++i) sing_vals << S[i] << " ";
//		sing_vals << std::endl;
//		sing_vals.close();
//		//////////////////////////////////////////////////////////


		// Q' -> Qb (transpose only first k_col rows of Q', discard the rest)
		Qb = new T [nnodes*k_col];
		blas::transpose(k_col, nnodes, Q, nnodes, Qb);

		// K_row //////////////////////////////////////////////////////////////
		T *const K_row = C;

		const unsigned int info_row
			= blas::gesvdSO(nnodes, ninteractions*nnodes, K_row, S, Q, nnodes,
											LWORK, WORK);
		if (info_col!=0)
			throw std::runtime_error("SVD did not converge with " + info_row);
		const unsigned int k_row = getRank(S, nnodes, eps);
		delete [] WORK;

		// Q -> Qu
		Qu = Q;

		// V' -> V (transpose only first k_row rows of V'(K_row), discard the rest)
		T *const V = new T [nnodes*ninteractions * k_row];
		blas::transpose(k_row, nnodes*ninteractions, K_row, nnodes, V);

		// rank k(eps) /////////////////////////////////////////////////////
		k = k_row < k_col ? k_row : k_col;

		// C_row ///////////////////////////////////////////////////////////
		C = new T [ninteractions * k*k];
		for (unsigned int i=0; i<k; ++i) {
			blas::scal(nnodes*ninteractions, S[i], V + i*nnodes*ninteractions);
			for (unsigned int m=0; m<ninteractions; ++m)
				for (unsigned int j=0; j<k; ++j)
					C[m*k*k + j*k + i]
						= blas::scpr(nnodes,
												 V + i*nnodes*ninteractions + m*nnodes,
												 Qb + j*nnodes);
		}

		delete [] V;
		delete [] S;
		delete [] K_row;
	}









	/**
	 * The M2L operator matrix for the ith transfer operator is either K_i = U
	 * V'_i or K_i = A_i B'. Hence the operator can be compressed: A QR
	 * decomposition of all U = Qu Ru, V_i = Qv_i Rv, A_i = Qa_i Ra, B = Qb Rb is
	 * performed. Next phi = Ru Rv' (the same as psi = Ra Rb'). The matrices Qu,
	 * Qv_i, Qa_i and Qb are unitary, thus the ith M2L operator matrix can be
	 * rewritten as K_i = Qu C_i Qb, with C_i = phi Qv'_i Qb (the same as C_i =
	 * Qu' Qa_i psi). Hence two 'nnodes*k' matrices Qu and Qb and ninteraction
	 * 'k*k' transfer matrices have to be computed.
	 *
	 * Qu is stored in U
	 *
	 * Qb is stored in B
	 *
	 * C_i is stored in V_i 
	 */
	template <typename T>
	void UCiB(const unsigned int nnodes,
						const unsigned int ninteractions,
						T *const U, T *const V, T *const A, T *const B,
						const unsigned int k)
	{
		// compute compressed M2L operators
		unsigned int INFO;
		const unsigned int lwork = k*k;
		T* work = new T [lwork];
	
	
		// QR of U and V_i
		T* tauU = new T [k];
		INFO = blas::geqrf(nnodes, k, U, tauU, lwork, work);
		assert(INFO==0);
		T* tauV = new T [k];
		INFO = blas::geqrf(nnodes*ninteractions, k, V, tauV, lwork, work);
		assert(INFO==0);
		// phi = Ru Rv'
		T* rU = new T [2 * k*k];
		T* rV = rU + k*k;
		blas::setzero(2* k*k, rU);
		for (unsigned int l=0; l<k; ++l) {
			blas::copy(l+1, U + l*nnodes,               rU + l*k);
			blas::copy(l+1, V + l*nnodes*ninteractions, rV + l*k);
		}
		T* phi = new T [k*k];
		blas::gemmh(k, k, k, T(1.), rU, k, rV, k, phi, k);
		delete [] rU;
		// get Qu and Qv
		INFO = blas::orgqr(nnodes,               k, U, tauU, lwork, work);
		assert(INFO==0);
		INFO = blas::orgqr(nnodes*ninteractions, k, V, tauV, lwork, work);
		assert(INFO==0);
		delete [] tauU;
		delete [] tauV;
	
	
		// QR of A_i and B
		T* tauA = new T [k];
		INFO = blas::geqrf(nnodes*ninteractions, k, A, tauA, lwork, work);
		assert(INFO==0);
		T* tauB = new T [k];
		INFO = blas::geqrf(nnodes, k, B, tauB, lwork, work);
		assert(INFO==0);
		// get Qa and Qb
		INFO = blas::orgqr(nnodes*ninteractions, k, A, tauA, lwork, work);
		assert(INFO==0);
		INFO = blas::orgqr(nnodes,               k, B, tauB, lwork, work);
		assert(INFO==0);
		delete [] tauA;
		delete [] tauB;
	
	
		// compressed operator Cuv stored in V
		T* aux = new T [ninteractions * k*k];
		for (unsigned int i=0; i<ninteractions; ++i)
			for (unsigned int n=0; n<k; ++n)
				for (unsigned int m=0; m<k; ++m)
					aux[i*k*k + n*k + m] = 
						blas::scpr(nnodes, V + m*ninteractions*nnodes + i*nnodes,
											 B + n*nnodes);
		for (unsigned int i=0; i<ninteractions; ++i)
			blas::gemm(k, k, k, T(1.),
								 phi, k, aux + i*k*k, k, V + i*k*k, k);
	
		delete [] aux;
		delete [] work;
		delete [] phi;
	}









} // namespace dfmm









//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



template <int> class map_loc_glob;
template <int> class Chebyshev;
template <int,int> class Tensor;
template <typename,typename> class EntryComputer;






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// COMPUTE Kt ~ UCB' /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



/*! \ingroup mm2l

	\class M2LComputerSA

	\brief Computes M2L operator as \f$K_i \sim U C_i B^*\f$

	The M2L operator, also a matrix \f$ K_i \in C^{\ell^d \times \ell^d} \f$ is
	approximated by using the adaptive cross approximation (ACA) and further
	compressed by using a QR decomposition to \f$K_i \approx Q_U C_i Q_B^*\f$,
	with \f$Q_U \in C^{\ell^d\times k}\f$, \f$C_i \in C^{k\times k}\f$ and
	\f$Q_B \in C^{\ell^d\times k}\f$.

	@tparam DIM spatial dimension
	@tparam ORDER interpolation order
	@tparam T value type
*/
template <int DIM,
					int ORDER,
					typename T>
class M2LComputerSA
{
	static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
	typedef typename DimTraits<DIM>::point_type     point_type;
	typedef typename DimTraits<DIM>::point_int_type point_int_type;


	enum {LRNK=0, ISZE=1, RCMP=2, DATA=3};
	typedef boost::tuple<unsigned int,unsigned int,unsigned int*,T*> data_tuple;


	/*!
		Store compressed M2L operators
	*/
	void setUCB(const unsigned int ninteractions,
							const unsigned int k,
              T *const Qu,
              T *const C,
              T *const Qb,
              const unsigned int c=0)
	{
    assert(k>0);
		static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
    T *const data = new T [k*k*ninteractions + 2 * nnodes*k];
    blas::copy(k*nnodes,            Qu, data);                // Qu
    blas::copy(k*nnodes,            Qb, data + nnodes*k);     // Qb
		blas::copy(k*k * ninteractions, C,  data + nnodes*k * 2); // C

		unsigned int *const rcmp = new unsigned int [ninteractions];
		for (unsigned int i=0; i<ninteractions; ++i) rcmp[i] = 0;

    operators.insert(std::make_pair(c, data_tuple(k,ninteractions,rcmp,data)));
  }


	template <typename kernel_type>
	void usetSVD(const point_type *const X,
							 const std::vector<point_int_type>& interactions,
							 const double diam,
							 const kernel_type& kernel,
							 const double eps,
							 T* &Qu, T* &C, T* &Qb, unsigned int &k)
	{
		// compute all K_i
		const unsigned int ninteractions = interactions.size();
		C = new T [nnodes*nnodes * ninteractions];
		point_type Y[nnodes];
		for (unsigned int t=0; t<ninteractions; ++t) {
			for (unsigned int n=0; n<nnodes; ++n)
				for (unsigned int d=0; d<DIM; ++d)
					Y[n][d] = X[n][d] + interactions.at(t)[d] * diam;
			EntryComputer<kernel_type,point_type> ec(nnodes, X, nnodes, Y, kernel);
			ec(0, nnodes, 0, nnodes, C + t*nnodes*nnodes);
			assert(check_nan(nnodes*nnodes, C + t*nnodes*nnodes)); 
		}
		
		// K_i = Qu * C * Qb'
		dfmm::UCiB(nnodes, ninteractions, eps, Qu, C, Qb, k);
	}
	

	template <typename kernel_type>
	void usefACA(const point_type *const X,
							 const std::vector<point_int_type>& interactions,
							 const double diam,
							 const kernel_type& kernel,
							 const double eps,
							 T* &Qu, T* &C, T* &Qb, unsigned int &k)
	{
		// compute full block row Kr and full block column Kc
		const unsigned int ninteractions = interactions.size();
		T *const Kr = new T [nnodes*nnodes * ninteractions];
		T *const Kc = new T [ninteractions * nnodes*nnodes];
		point_type Y[nnodes];
		for (unsigned int t=0; t<ninteractions; ++t) {
			for (unsigned int n=0; n<nnodes; ++n)
				for (unsigned int d=0; d<DIM; ++d)
					Y[n][d] = X[n][d] + interactions.at(t)[d] * diam;
			EntryComputer<kernel_type,point_type> ec(nnodes, X, nnodes, Y, kernel);
			// row
			ec(0, nnodes, 0, nnodes, Kr + t*nnodes*nnodes);
			assert(check_nan(nnodes*nnodes, Kr + t*nnodes*nnodes)); 
			// col
			for (unsigned int n=0; n<nnodes; ++n)
				blas::copy(nnodes,
									 Kr + t*nnodes*nnodes + n*nnodes,
									 Kc + n*nnodes*ninteractions + t*nnodes);
		}

		unsigned int kr, kc;
		T *Ur, *Vr, *Uc, *Vc; 

		dfmm::fACA(Kr, nnodes, nnodes*ninteractions, eps, Ur, Vr, kr);
		dfmm::fACA(Kc, nnodes*ninteractions, nnodes, eps, Uc, Vc, kc);

		delete [] Kr;
		delete [] Kc;
		
#if DFMM_USE_QRTSVD
		dfmm::QRtSVD(nnodes, nnodes*ninteractions, eps, Ur, Vr, kr);
		dfmm::QRtSVD(nnodes*ninteractions, nnodes, eps, Uc, Vc, kc);
#endif

		// set the rank to min(kr, kc);
		k = kr < kc ? kr : kc;

		dfmm::UCiB(nnodes, ninteractions, Ur, Vr, Uc, Vc, k);

		Qu = Ur;
		C  = Vr;
		delete [] Uc;
		Qb = Vc;
	}



	template <typename kernel_type>
	void usepACA(const point_type *const X,
							 const std::vector<point_int_type>& interactions,
							 const double diam,
							 const kernel_type& kernel,
							 const double eps,
							 T* &Qu, T* &C, T* &Qb, unsigned int &k)
	{
		// compute full block row Kr and full block column Kc
		const unsigned int ninteractions = interactions.size();
		point_type *const Yr = new point_type [nnodes*ninteractions];
		point_type *const Yc = new point_type [ninteractions*nnodes];
		for (unsigned int t=0; t<ninteractions; ++t)
			for (unsigned int n=0; n<nnodes; ++n)
				for (unsigned int d=0; d<DIM; ++d) {
					Yr[t*nnodes + n][d] = X[n][d] + interactions.at(t)[d] * diam;
					Yc[t*nnodes + n][d] = X[n][d] - interactions.at(t)[d] * diam;
				}

		unsigned int kr, kc;
		T *Ur, *Vr, *Uc, *Vc; 

		EntryComputer<kernel_type,point_type>
			ecr(nnodes, X, nnodes*ninteractions, Yr, kernel);
		dfmm::pACA(ecr, nnodes, nnodes*ninteractions, eps, Ur, Vr, kr);

		EntryComputer<kernel_type,point_type>
			ecc(nnodes*ninteractions, Yc, nnodes, X, kernel);
		dfmm::pACA(ecc, nnodes*ninteractions, nnodes, eps, Uc, Vc, kc);

#if DFMM_USE_QRTSVD
		dfmm::QRtSVD(nnodes, nnodes*ninteractions, eps, Ur, Vr, kr);
		dfmm::QRtSVD(nnodes*ninteractions, nnodes, eps, Uc, Vc, kc);
#endif

		delete [] Yr;
		delete [] Yc;

		// set the rank to min(kr, kc);
		k = kr < kc ? kr : kc;

		dfmm::UCiB(nnodes, ninteractions, Ur, Vr, Uc, Vc, k);

		Qu = Ur;
		C  = Vr;
		delete [] Uc;
		Qb = Vc;
	}



	mutable FlopCounter fcount; 


protected:
	explicit M2LComputerSA(FlopCounter& _fcount)
		: fcount(_fcount)
	{}

	M2LComputerSA& operator= (const M2LComputerSA& other)
	{
		fcount = other.fcount;
		return *this;
	}

  typedef std::map<unsigned int,data_tuple> m2l_map;


  m2l_map operators;


	~M2LComputerSA()
  {
		BOOST_FOREACH(typename m2l_map::value_type& dp, operators) {
			delete [] boost::get<RCMP>(dp.second);
			delete [] boost::get<DATA>(dp.second);
		}
  }


	/// Compute **M2L**
	template <typename transfer_map,
						typename kernel_type>
	void compute(const transfer_map& ffield,
							 const Level<DIM> *const level,
							 const kernel_type& kernel,
							 const double eps)
	{
		// don't do anything if no cone has interactions
		if (ffield.empty()) return;

		// typedefs
		typedef Chebyshev< ORDER>  basis_type;
		typedef Tensor<DIM,ORDER> tensor_type; 

    // target points (X)
		const double diam = level->getDiam();
    map_loc_glob<DIM> map_l2g(point_type(0.), diam);
		point_type X[nnodes];
    for (unsigned int n=0; n<nnodes; ++n)
      for (unsigned int d=0; d<DIM; ++d)
        X[n][d] = map_l2g(d,basis_type::nodes[tensor_type::node_ids[n][d]]);

		// loop over all cones
    unsigned int kmax =  0;
    unsigned int kmin = -1;
    unsigned int ksum =  0;
    BOOST_FOREACH(const typename transfer_map::value_type& tp, ffield) {
			
			const std::vector<point_int_type>& interactions = tp.second;
			const unsigned int ninteractions = interactions.size();
			
			// memory for compressed M2L kernels
			unsigned int k;
			T *Qu, *C, *Qb;
			
#if DFMM_USE_TSVD
			usetSVD(X, interactions, diam, kernel, eps, Qu, C, Qb, k);
#elif DFMM_USE_FACA
			usefACA(X, interactions, diam, kernel, eps, Qu, C, Qb, k);
#elif DFMM_USE_PACA
			usepACA(X, interactions, diam, kernel, eps, Qu, C, Qb, k);
#endif			

			// rank logging
			ksum += k;
			if (k>kmax) kmax = k;
			if (k<kmin) kmin = k;
			
			// set UCB
			const unsigned int c = tp.first;
			this->setUCB(ninteractions, k, Qu, C, Qb, c);
			delete [] Qu;
			delete [] C;
			delete [] Qb;
			
		}

		// to be verbose
		const double kavg = (double)ksum / (double)ffield.size();
		std::cout << "[kavg,kmax,kmin]:[" << kavg << "," << kmax << ","
							<< kmin << "] ";

  }





	//! \name Qu C Qb'
	/* @{ */
  /// Expand potentials X -> x
  void applyU(T *const X,	const unsigned int c=0) const
  {
    const data_tuple dt = operators.at(c);
    const unsigned int k = boost::get<LRNK>(dt);
    const T *const Qu    = boost::get<DATA>(dt);
    blas::gemva(nnodes, k, 1., const_cast<T *const>(Qu), X+nnodes, X);

		// count flops
		fcount.addFlops(k*(2*nnodes-1) + k);
  }

  /// Compressed M2L operation
  void applyC(const unsigned int t_id,
						 const T *const Y,
						 T *const X,
						 const unsigned int c=0) const
  {
    const data_tuple dt = operators.at(c);
    const unsigned int k    = boost::get<LRNK>(dt);
    T *const C              = boost::get<DATA>(dt) + 2*nnodes*k + t_id*k*k;
		const unsigned int rcmp = boost::get<RCMP>(dt)[t_id];

    assert(check_nan(k, X+nnodes)); 
    assert(check_nan(k, Y+nnodes)); 

		// not recompressed
		if (!rcmp) {
			blas::gemva(k, k, 1., C, const_cast<T*>(Y+nnodes), X+nnodes);
			fcount.addFlops(k*(2*k-1) + k);
		}
		// recompressed
		else {
			dfmm::applyUV(rcmp, k, C, k, C+rcmp*k, Y+nnodes, X+nnodes);
			fcount.addFlops((k*(2*rcmp-1) + rcmp*(2*k-1)) + k);
		}

    assert(check_nan(k, X+nnodes)); 
    assert(check_nan(k, Y+nnodes)); 
  }

  /// Compress densities Y -> Y+nnodes
  void applyB(T *const Y,	const unsigned int c=0) const
  {
    const data_tuple dt = operators.at(c);
    const unsigned int k = boost::get<LRNK>(dt);
    const T *const Qb    = boost::get<DATA>(dt) + nnodes*k;
    blas::gemhv(nnodes, k, 1., const_cast<T *const>(Qb), Y, Y+nnodes);

		// count flops
		fcount.addFlops(nnodes*(2*k-1));
  }
	/* @} */



	/*!
		\todo The fully pivoted ACA doesn't work here, need to find out why!
	*/
	void recompress(const double eps) 
	{
		const double t = omp_get_wtime();

		unsigned int enentries = 0; // effective number of entries
		unsigned int onentries = 0; //   overall number of entries

		BOOST_FOREACH(typename m2l_map::value_type& dp, operators) {
			data_tuple& dt = dp.second;
			
			const unsigned int k             = boost::get<LRNK>(dt);
			const unsigned int ninteractions = boost::get<ISZE>(dt);
			unsigned int *const rcmp         = boost::get<RCMP>(dt);
			T *const C                       = boost::get<DATA>(dt) + nnodes*k * 2;
			
			for (unsigned int i=0; i<ninteractions; ++i) {

				if (rcmp[i]) throw std::runtime_error("Recompress only once!");
				
				T *const Ci = C + i * k*k;
				T *U, *V;
				unsigned int kr;

				///////////////////////////////////////////////////////////////
				// fACA does not work here - need to find out why! ////////////
				///////////////////////////////////////////////////////////////
				dfmm::tSVD(k, k, Ci, eps, U, V, kr);
				//dfmm::fACA(Ci, k, k, eps, U, V, kr);
				//dfmm::QRtSVD(  k, k, eps, U, V, kr);

				const unsigned int maxrank = ((k-1) / 2) + 1;

				if (kr < maxrank) {
					blas::copy(k*kr, U, Ci);
					blas::copy(k*kr, V, Ci + k*kr);
					rcmp[i] = kr;
					enentries += 2 * k*kr;
				} else enentries += k*k;

				onentries += k*k;

				delete [] U;
				delete [] V;
			}
		}

		std::cout << "\t(rcmp in " << omp_get_wtime() - t << "s: "
							<< enentries << "/" << onentries << " = "
							<< 100. * (double)enentries / (double)onentries << "%)\t"
							<< std::flush;

	}

	
};







//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// COMPUTE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////




/*! \ingroup mm2l

	\class M2LComputerNA

	\brief Computes M2L operator as \f$K_i\f$

	This class computes, stores and applies dense M2L operators.

	@tparam DIM spatial dimension
	@tparam ORDER interpolation order
	@tparam T value type
*/
template <int DIM,
					int ORDER,
					typename T>
class M2LComputerNA
{
  typedef std::vector<const T*> m2l_vec;


protected:
  m2l_vec operators; //!< m2l operators

	/// Destructor
	~M2LComputerNA()
  { BOOST_FOREACH(const T *const K, operators) delete [] K; }


	/// Compute **M2L**
	template <typename transfer_vec,
						typename kernel_type>
	void compute(const transfer_vec& ffield,
							 const Level<DIM> *const level,
							 const kernel_type& kernel)
	{
		// don't do anything if no cone has interactions
		if (ffield.empty()) return;

		// typedefs
		static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
		typedef typename DimTraits<DIM>::point_int_type point_int_type;
		typedef typename DimTraits<DIM>::point_type     point_type;
		typedef Chebyshev< ORDER>  basis_type;
		typedef Tensor<DIM,ORDER> tensor_type; 

    // target points (X)
		const double diam = level->getDiam();
    map_loc_glob<DIM> map_l2g(point_type(0.), diam);
		point_type X[nnodes];
    for (unsigned int n=0; n<nnodes; ++n)
      for (unsigned int d=0; d<DIM; ++d)
        X[n][d] = map_l2g(d,basis_type::nodes[tensor_type::node_ids[n][d]]);

		// loop over all interactions and compute K
		assert(operators.empty());
		operators.clear();
    BOOST_FOREACH(const point_int_type& t, ffield) {
			point_type Y[nnodes];
			for (unsigned int n=0; n<nnodes; ++n)
				for (unsigned int d=0; d<DIM; ++d)
					Y[n][d] = X[n][d] + t[d] * diam;
			T *const K = new T [nnodes * nnodes];
			EntryComputer<kernel_type,point_type> ec(nnodes, X, nnodes, Y, kernel);
			ec(0, nnodes, 0, nnodes, K);
			assert(check_nan(nnodes*nnodes, K)); 
			operators.push_back(K);
		}
		
		// to be verbose
		std::cout << "[dense] ";

  }


	//! compute x += K y
  static void applyK(const unsigned int nx,
                     const unsigned int ny,
                     const T *const K,
                     T *const y,
                     T *const x)
  { blas::gemva(nx, ny, 1., const_cast<T *const>(K), y, x); }

	
};








//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// COMPUTE K ~ UV' ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////










/*! \ingroup mm2l

	\class M2LComputerIA

	\brief Computes M2L operator as \f$K_i \sim U_i V_i^*\f$

	This class computes, stores and applies low-rank M2L operators.

	@tparam DIM spatial dimension
	@tparam ORDER interpolation order
	@tparam T value type
*/
template <int DIM,
					int ORDER,
					typename T>
class M2LComputerIA
{
	typedef std::pair<unsigned int,const T*> data_pair;


	//! Store UV
	void pushback(const unsigned int k,	const T *const U,	const T *const V)
	{
		// rank must be larger then 0
		assert(k>0);

		// allocate memory for compressed m2l operator
		static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
		T *const data = new T [2 * nnodes*k];
		
    // copy U
  	blas::copy(k*nnodes, const_cast<T *const>(U), data);
		// copy V
		blas::copy(k*nnodes, const_cast<T *const>(V), data + k*nnodes);

		// assign low-rank and compressed m2l operator
		operators.push_back(std::make_pair(k, data));
  }


protected:
  typedef std::vector<data_pair> m2l_vec;
  m2l_vec operators; //!< m2l operators

	/// Destructor
	~M2LComputerIA()
  { BOOST_FOREACH(data_pair& dp, operators)	delete [] dp.second; }


	/// Compute **M2L**
	template <typename transfer_vec,
						typename kernel_type>
	void compute(const transfer_vec& ffield,
							 const Level<DIM> *const level,
							 const kernel_type& kernel,
							 const double eps)
	{
		// don't do anything if no cone has interactions
		if (ffield.empty()) return;

		// typedefs
		static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;
		typedef typename DimTraits<DIM>::point_type     point_type;
		typedef typename DimTraits<DIM>::point_int_type point_int_type;
		typedef Chebyshev< ORDER>  basis_type;
		typedef Tensor<DIM,ORDER> tensor_type; 

    // target points (X)
		const double diam = level->getDiam();
    map_loc_glob<DIM> map_l2g(point_type(0.), diam);
		point_type X[nnodes];
    for (unsigned int n=0; n<nnodes; ++n)
      for (unsigned int d=0; d<DIM; ++d)
        X[n][d] = map_l2g(d,basis_type::nodes[tensor_type::node_ids[n][d]]);

		// loop over all interactions
    unsigned int kmax =  0;
    unsigned int kmin = -1;
    unsigned int ksum =  0;
    BOOST_FOREACH(const point_int_type& t, ffield) {
			
			// setup interpolation point Y in source cluster
			point_type Y[nnodes];
			for (unsigned int n=0; n<nnodes; ++n)
				for (unsigned int d=0; d<DIM; ++d)
					Y[n][d] = X[n][d] + t[d] * diam;

			// setup entry computer
#if DFMM_USE_CHEB_WEIGHTS
			double weights[nnodes];
			tensor_type::setRootOfWeights(weights);
			EntryComputer<kernel_type,point_type> ec(nnodes, X,
																							 nnodes, Y, kernel, weights);
#else
			EntryComputer<kernel_type,point_type> ec(nnodes, X, nnodes, Y, kernel);
#endif

			// declare variables for low-rank representation K ~ UV'
			unsigned int k;
			T *U, *V;

#if DFMM_USE_TSVD
			// Truncated singular value decomposition
			T *const K = new T [nnodes*nnodes];
			ec(0, nnodes, 0, nnodes, K);
			assert(check_nan(nnodes*nnodes, K)); 
			dfmm::tSVD(nnodes, nnodes, K, eps, U, V, k);
			delete [] K;
#elif DFMM_USE_FACA
			// Fully pivoted adaptive cross approximation
			T *const K = new T [nnodes*nnodes];
			ec(0, nnodes, 0, nnodes, K);
			assert(check_nan(nnodes*nnodes, K)); 
			dfmm::fACA(K, nnodes, nnodes, eps, U, V, k);
			delete [] K;
#endif			


#if DFMM_USE_PACA
			// Partially pivoted adaptive cross approximation
			dfmm::pACA(ec, nnodes, nnodes, eps, U, V, k);
#endif

#if DFMM_USE_QRTSVD
			// QR & truncated SVD after (any form of) ACA
			dfmm::QRtSVD(   nnodes, nnodes, eps, U, V, k);
#endif
		
#if DFMM_USE_CHEB_WEIGHTS
			for (unsigned int n=0; n<nnodes; ++n) {
				blas::scal(k, 1./weights[n], U+n, nnodes); // un-weight U (cols of K)
				blas::scal(k, 1./weights[n], V+n, nnodes); // un-weight V (rows of K)
			}
#endif

			// store UV
			this->pushback(k, U, V);
			delete [] U;
			delete [] V;
			
			// rank logging
			ksum += k;
			if (k>kmax) kmax = k;
			if (k<kmin) kmin = k;
		}

		// to be verbose
		const double kavg = (double)ksum / (double)ffield.size();
		std::cout << "[kavg,kmax,kmin]:[" << kavg << "," << kmax << ","
							<< kmin << "] ";
  }

	
	//! compute x += UV' y	
  static void applyUV(const unsigned int k,
                      const unsigned int nx,
                      const T *const U,
                      const unsigned int ny,
                      const T *const V,
                      const T *const y,
                      T *const x)
  {
    for (unsigned int l=0; l<k; ++l) {
      const T e = blas::scpr(ny, V + l*ny, y);
      blas::axpy(nx, e, U + l*nx, x);
    }
  }

};


#endif

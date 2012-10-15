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

#ifndef permutations_hpp
#define permutations_hpp


template <int DIM,
					int ORDER>
class PermutationHandler;


/*!
	\ingroup mtranslation

	\class PermutationHandler 

	\brief Handles permutations and tensor-product M2M and L2L operator
	application

	@tparam DIM spatial dimension
	@tparam ORDER interpolation order

 */
template <int ORDER>
class PermutationHandler<3,ORDER>
{
	static const unsigned int DIM = 3;
	static const unsigned int nnodes = BasisTraits<ORDER,DIM>::nnodes;

	unsigned int perm[DIM][nnodes];

public:
	PermutationHandler() {
		for (unsigned int i=0; i<ORDER; ++i)
			for (unsigned int j=0; j<ORDER; ++j) 
				for (unsigned int k=0; k<ORDER; ++k) {
					const unsigned int index = k*ORDER*ORDER + j*ORDER + i;
					perm[0][index] = k*ORDER*ORDER + j*ORDER + i;
					perm[1][index] = i*ORDER*ORDER + k*ORDER + j;
					perm[2][index] = j*ORDER*ORDER + i*ORDER + k;
				}
	}


	// 3 * N * (ORDER*ORDER*(2*ORDER-1)) + nnodes
	template <typename T>
	void m2m(const T *const S, const unsigned int N,
					 const T *const y, T *const Y) const
	{
		T exp[N*nnodes], pexp[N*nnodes];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemhm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
								const_cast<T*>(S), ORDER,
								const_cast<T*>(y), ORDER, pexp, ORDER);
		for (unsigned int e=0; e<N; ++e)
			for (unsigned int n=0; n<nnodes; ++n)
				exp[e*nnodes + n] = pexp[e*nnodes + perm[1][n]];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemhm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
								const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
								exp, ORDER, pexp, ORDER);
		for (unsigned int e=0; e<N; ++e)
			for (unsigned int n=0; n<nnodes; ++n)
				exp[e*nnodes + perm[1][n]] = pexp[e*nnodes + perm[2][n]];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemhm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
								const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
								exp, ORDER, pexp, ORDER);
		// N * nnodes
		for (unsigned int e=0; e<N; ++e)
			for (unsigned int n=0; n<nnodes; ++n)
				Y[e*nnodes + perm[2][n]] += pexp[e*nnodes + n];
	}


	// 3 * N * (ORDER*ORDER*(2*ORDER-1)) + nnodes
	template <typename T>
	void l2l(const T *const S, const unsigned int N,
					 const T *const X, T *const x) const
	{
		T exp[N*nnodes], pexp[N*nnodes];
		// ORDER * (2*ORDER-1) * N*ORDER
		blas::gemm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
							 const_cast<T*>(S), ORDER,
							 const_cast<T*>(X), ORDER, pexp, ORDER);
		for (unsigned int e=0; e<N; ++e)
			for (unsigned int n=0; n<nnodes; ++n)
				exp[e*nnodes + n] = pexp[e*nnodes + perm[1][n]];
		// ORDER * (2*ORDER-1) * N*ORDER
		blas::gemm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
							 const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
							 exp, ORDER, pexp, ORDER);
		for (unsigned int e=0; e<N; ++e)
			for (unsigned int n=0; n<nnodes; ++n)
				exp[e*nnodes + perm[1][n]] = pexp[e*nnodes + perm[2][n]];
		// ORDER * (2*ORDER-1) * N*ORDER
		blas::gemm(ORDER, ORDER, N * ORDER*ORDER, T(1.),
							 const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
							 exp, ORDER, pexp, ORDER);
		// N * ORDER*ORDER*ORDER
		for (unsigned int e=0; e<N; ++e)
			for (unsigned int n=0; n<nnodes; ++n)
				x[e*nnodes + perm[2][n]] += pexp[e*nnodes + n];
	}

	
	// 3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes
	template <typename T>
	void m2m(const T *const S,	const T *const y,	T *const Y) const
	{
		T exp[nnodes], pexp[nnodes];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemhm(ORDER, ORDER, ORDER*ORDER, T(1.),
								const_cast<T*>(S), ORDER,
								const_cast<T*>(y), ORDER, pexp, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	exp[n] = pexp[perm[1][n]];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemhm(ORDER, ORDER, ORDER*ORDER, T(1.),
								const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
								exp, ORDER, pexp, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	exp[perm[1][n]] = pexp[perm[2][n]];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemhm(ORDER, ORDER, ORDER*ORDER, T(1.),
								const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
								exp, ORDER, pexp, ORDER);
		// ORDER*ORDER*ORDER
		for (unsigned int n=0; n<nnodes; ++n)	Y[perm[2][n]] += pexp[n];
	}


	// 3 * (ORDER*ORDER*(2*ORDER-1)) + nnodes
	template <typename T>
	void l2l(const T *const S, const T *const X, T *const x) const
	{
		T exp[nnodes], pexp[nnodes];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemm(ORDER, ORDER, ORDER*ORDER, T(1.),
							 const_cast<T*>(S), ORDER,
							 const_cast<T*>(X), ORDER, pexp, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	exp[n] = pexp[perm[1][n]];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemm(ORDER, ORDER, ORDER*ORDER, T(1.),
							 const_cast<T*>(S) + 2 * ORDER*ORDER, ORDER,
							 exp, ORDER, pexp, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	exp[perm[1][n]] = pexp[perm[2][n]];
		// ORDER * (2*ORDER-1) * ORDER
		blas::gemm(ORDER, ORDER, ORDER*ORDER, T(1.),
							 const_cast<T*>(S) + 1 * ORDER*ORDER, ORDER,
							 exp, ORDER, pexp, ORDER);
		// ORDER*ORDER*ORDER
		for (unsigned int n=0; n<nnodes; ++n)	x[perm[2][n]] += pexp[n];
	}

};



#endif

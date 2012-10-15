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

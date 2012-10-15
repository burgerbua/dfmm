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

#ifndef dof_hpp
#define dof_hpp


#include <boost/noncopyable.hpp>

#include "sources/dimordertraits.hpp"

/*! \class Dof
	Stands for degree of freedom (represents particles) and provides their
	geometric location and identification

	\brief degree of freedom (particles), stores location and its id
 */
template <int DIM>
class Dof
{
public:
	static const unsigned int dim = DIM;
  typedef typename DimTraits<DIM>::point_type point_type;

	explicit Dof() : x(point_type()), id(-1) {}
	explicit Dof(const point_type _x)	: x(_x), id(-1)	{}
	explicit Dof(const point_type _x, const unsigned int _id) : x(_x), id(_id) {}
  
	double getCoordinate(unsigned int i) const
	{	return x[i]; }

	const point_type& getPoint() const
	{	return x;	}

	unsigned int getId() const
	{ assert(id != -1); return id; }

	std::ostream& writeInfo(std::ostream &out = std::cout) const
  {
    out << x << std::endl;
    return out;
  }


private:
	point_type x;
	unsigned int id;
};


#endif

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

#ifndef conefactory_hpp
#define conefactory_hpp

#include <boost/math/constants/constants.hpp>

#include "./dimordertraits.hpp"

template <int> class ConeFactory;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// 2D
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <>
class ConeFactory<2>
{
	enum {dim = 2};
	typedef DimTraits<dim>::point_type point_type;

public:
	ConeFactory(const double _extension,
							const double _wavenum)
		: extension(_extension),
			wavenum(  _wavenum)
	{
		// set angle and number of cones
		ncones     = 1;
		cone_angle = 2. * boost::math::constants::pi<double>();
		while(cone_angle > 1./(wavenum*extension)) {
			cone_angle /= 2.;
			ncones     *= 2.;
		}
		// assertions
		assert(cone_angle<=1./(wavenum*extension));
		assert(ncones>0);
	}


	const unsigned int getNCones() const
	{
		return ncones;
	}


	const point_type getConeDirection(const unsigned int c) const
	{
		const double angle = (c + .5) * cone_angle;
		return point_type(cos(angle), sin(angle));
	}
	
	
	const unsigned int getConeIdx(point_type vecXY) const
	{
		// if within cone (the weird .99999... comes due to the fact that (acos(.)
		// <= .) does somehow not recognise if it is equal?!?)
		for (unsigned int c=0; c<ncones; ++c) {
			const double dotCONExXY
        = getConeDirection(c).dotProduct(vecXY /= vecXY.norm());
			if (acos(dotCONExXY*.9999999)*.99 <= cone_angle/2.)	return c;
		}
		// not possible
		throw std::runtime_error("Each direction is in a cone!");
		return -1;
	}


	const double getConeAngle() const
	{
		return cone_angle;
	}



private:
	const double extension;
	const double wavenum;

	unsigned int ncones;
	double cone_angle;
};






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// 1D
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <>
class ConeFactory<3>
{
	enum {dim = 3};
	typedef DimTraits<dim>::point_type point_type;

public:
	ConeFactory(const double _extension,
							const double _wavenum)
		: extension(_extension),
			wavenum(  _wavenum)
	{
		const double scale_cone_aperture = SCALE_CONE_APERTURE;

		// set angle and number of cones
		ncones     = 1;
		cone_angle = 2. * boost::math::constants::pi<double>();
//		while(cone_angle > 1./(wavenum*extension)) {
		while(cone_angle > scale_cone_aperture/(wavenum*extension) || ncones<4) {
			cone_angle /= 2.;
			ncones     *= 2.;
		}

		// number of cones per direction per pyramid
		ncones /= 4;

		// assertions
		assert(ncones>=1);
	}


	const unsigned int getNCones() const
	{
		return ncones*ncones * 6;	// number of cones is 6 times ncones^2
	}


	const point_type getConeDirection(const unsigned int c) const
	{
		assert(c<getNCones());
		const unsigned int pyramid = c / (ncones*ncones);
		const unsigned int coneidx = c % (ncones*ncones);
		const unsigned int x1 = coneidx / ncones;
		const unsigned int x2 = coneidx % ncones;
		const double angle1 = cone_angle * (-(double)ncones/2 + (double)x1 + .5);
		const double angle2 = cone_angle * (-(double)ncones/2 + (double)x2 + .5);

    point_type u;
		switch (pyramid) {
		case 0:  u[0]= 1.;         u[1]=tan(angle1); u[2]=tan(angle2); break;
    case 1:  u[0]=-1.;         u[1]=tan(angle1); u[2]=tan(angle2); break;
		case 2:  u[0]=tan(angle2); u[1]= 1.;         u[2]=tan(angle1); break;
		case 3:  u[0]=tan(angle2); u[1]=-1.;         u[2]=tan(angle1); break;
		case 4:  u[0]=tan(angle1); u[1]=tan(angle2); u[2]= 1.;         break;
		case 5:  u[0]=tan(angle1); u[1]=tan(angle2); u[2]=-1.;         break;
		default: throw std::runtime_error("getConeDirection: No valid pyramid!");
		};
    return (u /= u.norm());

	}


	const unsigned int getConeIdx(const point_type& u) const
	{
		unsigned int pyramid;
		if      (fabs(u[2])<= u[0] && fabs(u[1])<= u[0]) pyramid=0;
		else if (fabs(u[2])<=-u[0] && fabs(u[1])<=-u[0]) pyramid=1;
		else if (fabs(u[0])<= u[1] && fabs(u[2])<= u[1]) pyramid=2;
		else if (fabs(u[0])<=-u[1] && fabs(u[2])<=-u[1]) pyramid=3;
		else if (fabs(u[1])<= u[2] && fabs(u[0])<= u[2]) pyramid=4;
		else if (fabs(u[1])<=-u[2] && fabs(u[0])<=-u[2]) pyramid=5;
		else throw std::runtime_error("getConeIdx (pyramid): No valid pyramid!");

		if (ncones==1) return pyramid;

		double angle1, angle2;
		switch (pyramid) {
		case 0:
			angle1 =  atan(u[1]/u[0]);
			angle2 =  atan(u[2]/u[0]);
			break;
		case 1:
			angle1 = -atan(u[1]/u[0]);
			angle2 = -atan(u[2]/u[0]);
			break;
		case 2:
			angle1 =  atan(u[2]/u[1]);
			angle2 =  atan(u[0]/u[1]);
			break;
		case 3:
			angle1 = -atan(u[2]/u[1]);
			angle2 = -atan(u[0]/u[1]);
			break;
		case 4:
			angle1 =  atan(u[0]/u[2]);
			angle2 =  atan(u[1]/u[2]);
			break;
		case 5:
			angle1 = -atan(u[0]/u[2]);
			angle2 = -atan(u[1]/u[2]);
			break;
		default:
			throw std::runtime_error("getConeIdx(angle1,angle2): No valid pyramid!");
		};
		
		const unsigned int x1 = floor(angle1*.99/cone_angle) + ncones/2;
		const unsigned int x2 = floor(angle2*.99/cone_angle) + ncones/2;
		
		assert(x1 < ncones);
		assert(x2 < ncones);

		return pyramid * ncones*ncones + x1*ncones + x2;
	}


	const double getConeAngle() const
	{
		return cone_angle;
	}



private:
	const double extension;
	const double wavenum;

	unsigned int ncones;
	double cone_angle;
};



#endif

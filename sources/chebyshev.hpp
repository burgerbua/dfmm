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

#ifndef chebyshev_hpp
#define chebyshev_hpp

#include <cmath>




/*! \class Chebyshev
 
  \brief Provides the %Chebyshev polynomials of 1st and 2nd kind and roots

  This class provides the trigonometric definition of the Chebyshev poynomials
  of 1st and 2nd kind in [-1,1].
  
  @tparam ORDER interpolation order
 */
template <int ORDER>
struct Chebyshev
{
  enum {order = ORDER};

  /*! Chebyshev roots of given ORDER in [-1,1] */
  const static double nodes[];

  /**
   * Chebyshev polynomials of first kind:
   * \f$ T_n(x) = \cos(n \arccos(x)) \f$
   */
  const static double T(const unsigned int n, double x)
  {
    assert(std::fabs(x)-1.<10.*std::numeric_limits<double>::epsilon());
    if (std::fabs(x)>1.) {
      x = (x >  1. ?  1. : x);
      x = (x < -1. ? -1. : x);
    }

    return (cos(n * acos(x)));
  }


  /**
   * Chebyshev polynomials of second kind:
   * \f$ U_{n-1}(x) = \frac{\sin(n arccos(x))}{\sqrt{1-x^2}} \f$
   * \f$ \frac{\mathrm{d} T_n(x)}{\mathrm{d}x} = n U_{n-1}(x) \f$
   */
  const static double U(const unsigned int n, double x)
  {
    assert(std::fabs(x)-1.<10.*std::numeric_limits<double>::epsilon());
    if (std::fabs(x)>1.) {
      x = (x >  1. ?  1. : x);
      x = (x < -1. ? -1. : x);
    }
  
    return (n * (sin(n * acos(x))) / sqrt(1 - x*x));
//    return ((sin(n * acos(x))) / sin(acos(x)));
  }
};





// order 2
template<> const double Chebyshev<2>::nodes[] = {0.707106781186548,
                                                 -0.707106781186547};

/// order 3
template<> const double Chebyshev<3>::nodes[] = {8.66025403784439e-01,
                                                 0.0,
                                                 -8.66025403784438e-01};

// order 4
template<> const double Chebyshev<4>::nodes[] = {0.923879532511287,
                                                 0.382683432365090,
                                                 -0.382683432365090,
                                                 -0.923879532511287};

// order 5
template<> const double Chebyshev<5>::nodes[] = {9.51056516295154e-01,
                                                 5.87785252292473e-01,
                                                 0.0,
                                                 -5.87785252292473e-01,
                                                 -9.51056516295154e-01};

// order 6
template<> const double Chebyshev<6>::nodes[] = {0.965925826289068,
                                                 0.707106781186548,
                                                 0.258819045102521,
                                                 -0.258819045102521,
                                                 -0.707106781186547,
                                                 -0.965925826289068};

// order 7
template<> const double Chebyshev<7>::nodes[] = {9.74927912181824e-01,
                                                 7.81831482468030e-01,
                                                 4.33883739117558e-01,
                                                 0.0,
                                                 -4.33883739117558e-01,
                                                 -7.81831482468030e-01,
                                                 -9.74927912181824e-01};

// order 8
template<> const double Chebyshev<8>::nodes[] = {0.980785280403230,
                                                 0.831469612302545,
                                                 0.555570233019602,
                                                 0.195090322016128,
                                                 -0.195090322016128,
                                                 -0.555570233019602,
                                                 -0.831469612302545,
                                                 -0.980785280403230};

// order 9
template<> const double Chebyshev<9>::nodes[] = {9.84807753012208e-01,
                                                 8.66025403784439e-01,
                                                 6.42787609686539e-01,
                                                 3.42020143325669e-01,
                                                 0.0,
                                                 -3.42020143325669e-01,
                                                 -6.42787609686539e-01,
                                                 -8.66025403784438e-01,
                                                 -9.84807753012208e-01,};

// order 10
template<> const double Chebyshev<10>::nodes[] = {0.987688340595138,
                                                  0.891006524188368,
                                                  0.707106781186548,
                                                  0.453990499739547,
                                                  0.156434465040231,
                                                  -0.156434465040231,
                                                  -0.453990499739547,
                                                  -0.707106781186547,
                                                  -0.891006524188368,
                                                  -0.987688340595138};

// order 11
template<> const double Chebyshev<11>::nodes[] = {9.89821441880933e-01,
                                                  9.09631995354518e-01,
                                                  7.55749574354258e-01,
                                                  5.40640817455598e-01,
                                                  2.81732556841430e-01,
                                                  2.83276944882399e-16,
                                                  -2.81732556841430e-01,
                                                  -5.40640817455597e-01,
                                                  -7.55749574354258e-01,
                                                  -9.09631995354518e-01,
                                                  -9.89821441880933e-01};

// order 12
template<> const double Chebyshev<12>::nodes[] = {0.991444861373810,
                                                  0.923879532511287,
                                                  0.793353340291235,
                                                  0.608761429008721,
                                                  0.382683432365090,
                                                  0.130526192220052,
                                                  -0.130526192220052,
                                                  -0.382683432365090,
                                                  -0.608761429008721,
                                                  -0.793353340291235,
                                                  -0.923879532511287,
                                                  -0.991444861373810};


// order 13
template<> const double Chebyshev<13>::nodes[] = {9.92708874098054e-01,
                                                  9.35016242685415e-01,
                                                  8.22983865893656e-01,
                                                  6.63122658240795e-01,
                                                  4.64723172043769e-01,
                                                  2.39315664287558e-01,
                                                  -1.60812264967664e-16,
                                                  -2.39315664287557e-01,
                                                  -4.64723172043769e-01,
                                                  -6.63122658240795e-01,
                                                  -8.22983865893656e-01,
                                                  -9.35016242685415e-01,
                                                  -9.92708874098054e-01};






#endif

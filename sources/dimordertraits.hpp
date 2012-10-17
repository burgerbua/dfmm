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

#ifndef dimordertraits_hpp
#define dimordertraits_hpp


#include "./point.hpp"



/*! 
  \defgroup mtraits Traits classes: compile-time constants and types
 */




/*! \class DimTraits
  Compile time constants and data types depending on the spatial dimension

  \ingroup mtraits
  
  @tparam DIM spatial dimension
 */
template <int DIM> struct DimTraits;

// 1D
template <> struct DimTraits<1>
{ enum {nboxes = 2};
  typedef Point<double,1> point_type; 
  typedef Point<int,   1> point_int_type; 
  static const int child_pos[][1]; //!< relative child-cluster positions
};

// 2D
template <> struct DimTraits<2>
{ enum {nboxes = 4};
  typedef Point<double,2> point_type; 
  typedef Point<int,   2> point_int_type; 
  static const int child_pos[][2]; //!< relative child-cluster positions
};

// 3D
template <> struct DimTraits<3>
{ enum {nboxes = 8};
  typedef Point<double,3> point_type; 
  typedef Point<int,   3> point_int_type; 
  static const int child_pos[][3]; //!< relative child-cluster positions
};



// 1D
const int DimTraits<1>::child_pos[][1] = { { 1},
                                           {-1} };
// 2D
const int DimTraits<2>::child_pos[][2] = { { 1,  1},
                                           {-1,  1},
                                           {-1, -1},
                                           { 1, -1} };
// 3D
const int DimTraits<3>::child_pos[][3] = { { 1,  1,  1},
                                           {-1,  1,  1},
                                           {-1, -1,  1},
                                           { 1, -1,  1},
                                           { 1,  1, -1},
                                           {-1,  1, -1},
                                           {-1, -1, -1},
                                           { 1, -1, -1} };



/*! \class BasisTraits
  Compile time contants depending on the interpolation basis and spatial
  dimension

  \ingroup mtraits
  
  @tparam ORDER interpolation order
  @tparam DIM spatial dimension
 */
/* @{ */
template <int ORDER, int DIM> struct BasisTraits;

// 1D
template <int ORDER> struct BasisTraits<ORDER, 1>
{
  enum {nnodes = ORDER};
};

// 2D
template <int ORDER> struct BasisTraits<ORDER, 2>
{
  enum {nnodes = ORDER*ORDER};
};

// 3D
template <int ORDER> struct BasisTraits<ORDER, 3>
{
  enum {nnodes = ORDER*ORDER*ORDER};
};
/* @} */
#endif



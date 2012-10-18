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

#ifndef point_hpp
#define point_hpp



// system includes
#include <ostream>
#include <cmath>
#include <limits>
#include <cassert>




/*! /class Point
  Geometric point whose dimension and type of coordinates is given by template
  parameters T and DIM.

  \brief %Points can be used as position vectors, vectors, etc.

  @tparam T value type of coordinates
  @tparam DIM spatial dimension
 */
template<typename T, int DIM>
class Point
{
  template <typename U>
  friend std::ostream& operator << (std::ostream& output,
                                    const Point<U,DIM>& p);

public:

  typedef T value_type;                  
  static const unsigned int dim = DIM;

  explicit Point()
  { for (unsigned int d=0; d<dim; ++d) coords[d] = value_type(0.); }

  explicit Point(const value_type _coord)
  { for(unsigned int d=0; d<dim; ++d) coords[d] = _coord; }

  explicit Point(const value_type _coords[dim])
  { for(unsigned int d=0; d<dim; ++d) coords[d] = _coords[d]; }

  Point(const Point& other)
  { for(unsigned int d=0; d<dim; ++d) coords[d] = other.coords[d]; }


  explicit Point(const double x, const double y);
  explicit Point(const double x, const double y, const double z);

  explicit Point(const int x, const int y);
  explicit Point(const int x, const int y, const int z);




  const value_type& operator [] (const unsigned int d) const
  { assert(d<dim); return coords[d]; }
  value_type& operator [] (const unsigned int d)
  { assert(d<dim); return coords[d]; }



  Point& operator = (const Point& other)
  {
    for(unsigned int d=0; d<dim; ++d)
      coords[d] = other.coords[d];
    return *this;
  }


  bool operator == (const Point& other) const;




  /**
   * add one point to the other
   */
  Point& operator += (const Point& other)
  {
    for(unsigned int d=0; d<dim; ++d) coords[d] += other.coords[d];
    return *this;
  }



  /**
   * subtract one point from the other
   */
  Point& operator -= (const Point& other)
  {
    for(unsigned int d=0; d<dim; ++d) coords[d] -= other.coords[d];
    return *this;
  }



  /**
   * scalar multiplication
   * NOTE: only post-multiplications of type "point *= scalar" are possible
   */
  Point& operator *= (const value_type scal)
  {
    for(unsigned int d=0; d<dim; ++d) coords[d] *= scal;
    return *this;
  }



  /**
   * scalar division
   * NOTE: only post-divisions of type "point /= scalar" are possible
   */
  Point& operator /= (const value_type scal)
  {
    for(unsigned int d=0; d<dim; ++d) coords[d] /= scal;
    return *this;
  }



  /**
   * + operator: add two point
   */
  Point operator + (const Point& other) const
  {
    value_type res[dim];
    for(unsigned int d=0; d<dim; ++d) res[d] = coords[d] + other.coords[d];
    return Point(res);
  }



  /**
   * - operator: subtract two points
   */
  Point operator - (const Point& other) const
  {
    value_type res[dim];
    for(unsigned int d=0; d<dim; ++d) res[d] = coords[d] - other.coords[d];
    return Point(res);
  }



  /**
   * Unary minus operator. Negate all entries of a point.
   */
  Point operator - () const
  {
    value_type res[dim];
    for(unsigned int d=0; d<dim; ++d) res[d] = -coords[d];
    return Point(res);
  }



  /**
   * scalar multiplication
   * NOTE: only post-multiplications of type "point * scalar" are possible
   */
  Point operator * (const value_type scal) const
  {
    value_type res[dim];
    for(unsigned int d=0; d<dim; ++d) res[d] = coords[d] * scal;
    return Point(res);
  }



  /**
   * scalar division
   * NOTE: only post-divisions of type "point / scalar" are possible
   */
  Point operator / (const value_type scal) const
  {
    value_type res[dim];
    for(unsigned int d=0; d<dim; ++d) res[d] = coords[d] / scal;
    return Point(res);
  }







  /**
   * dot (inner) product of two vectors.
   * \f$ a = <\mathbf{u},\mathbf{v}> = {\sum\limits_i u_i v_i} \f$
   * @param[in] v vector to be multiplied
   * @return value of norm
   */
  value_type dotProduct(const Point& other) const;

  /**
   * \f$ L_2\f$ -norm of a vector.
   * \f$ |\mathbf{v}|_2 = \sqrt{\sum\limits_i v_i^2} \f$
   * @return value of norm
   */
  double norm() const;

  /**
   * Returns the Euclidian distance of
   * <tt>this</tt> to <tt>p</tt>, i.e. the <tt>L_2</tt> norm
   * of the difference vector between the two geometric points.
   */
  double distance(const Point& other) const;



  


private:

  value_type coords[dim]; //!< coords, stored in a simple C-array

};




// double point
template<> inline
Point<double,2>::Point(const double x, const double y)
{ coords[0] = x; coords[1] = y; }
template<> inline
Point<double,3>::Point(const double x, const double y, const double z)
{ coords[0] = x; coords[1] = y; coords[2] = z; }

// int point
template<> inline
Point<int,2>::Point(const int x, const int y)
{ coords[0] = x; coords[1] = y; }
template<> inline
Point<int,3>::Point(const int x, const int y, const int z)
{ coords[0] = x; coords[1] = y; coords[2] = z; }

// euclidian distance
template<> inline
double Point<double,2>::distance(const Point& other) const
{
  const double diff0 = coords[0] - other.coords[0];
  const double diff1 = coords[1] - other.coords[1];
  return sqrt(diff0*diff0 + diff1*diff1);
}
template<> inline
double Point<double,3>::distance(const Point& other) const
{
  const double diff0 = coords[0] - other.coords[0];
  const double diff1 = coords[1] - other.coords[1];
  const double diff2 = coords[2] - other.coords[2];
  return sqrt(diff0*diff0 + diff1*diff1 + diff2*diff2);
}

// norm
template<> inline
double Point<double,2>::norm() const
{
  const double n0 = fabs(coords[0]);
  const double n1 = fabs(coords[1]);
  return sqrt(n0*n0 + n1*n1);
}
template<> inline
double Point<double,3>::norm() const
{
  const double n0 = fabs(coords[0]);
  const double n1 = fabs(coords[1]);
  const double n2 = fabs(coords[2]);
  return sqrt(n0*n0 + n1*n1 + n2*n2);
}

// dot product
template <> inline
double Point<double,2>::dotProduct(const Point& other) const
{
  const double d0 = coords[0] * other.coords[0];
  const double d1 = coords[1] * other.coords[1];
  return d0 + d1;
}
template <> inline
double Point<double,3>::dotProduct(const Point& other) const
{
  const double d0 = coords[0] * other.coords[0];
  const double d1 = coords[1] * other.coords[1];
  const double d2 = coords[2] * other.coords[2];
  return d0 + d1 + d2;
}

// equal operator
template<> inline
bool Point<double,2>::operator == (const Point& other) const
{
  return (fabs(coords[0]-other.coords[0]) < 
          std::numeric_limits<double>::epsilon() &&
          fabs(coords[1]-other.coords[1]) <
          std::numeric_limits<double>::epsilon());
}
template<> inline
bool Point<double,3>::operator == (const Point& other) const
{
  return (fabs(coords[0]-other.coords[0]) <
          std::numeric_limits<double>::epsilon() &&
          fabs(coords[1]-other.coords[1]) <
          std::numeric_limits<double>::epsilon() &&
          fabs(coords[2]-other.coords[2]) <
          std::numeric_limits<double>::epsilon());
}

template<> inline
bool Point<int,2>::operator == (const Point& other) const
{
  return (coords[0] == other.coords[0] &&
          coords[1] == other.coords[1]);
}
template<> inline
bool Point<int,3>::operator == (const Point& other) const
{
  return (coords[0] == other.coords[0] &&
          coords[1] == other.coords[1] &&
          coords[2] == other.coords[2]);
}


// friend function
template <typename T>
std::ostream& operator << (std::ostream& output, const Point<T,2>& p)
{
  output << "(" <<  p.coords[0] << "," << p.coords[1] << ")";
  return output;
}
template <typename T>
std::ostream& operator << (std::ostream& output, const Point<T,3>& p)
{
  output << "(" <<  p.coords[0] << "," << p.coords[1] << "," << p.coords[2]
         << ")";
  return output;
}




#endif //include guard

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

// system includes
#include <iostream>
#include <omp.h>



// qbdfmm
#include "sources/dimordertraits.hpp"
#include "sources/point.hpp"









int main( int argc, char * argv[ ] )
{
  // start timer
  const double t_0 = omp_get_wtime( );

  const unsigned int dim = 3;  
  
  typedef DimTraits<dim>::point_type     point_type;
  typedef DimTraits<dim>::point_int_type point_int_type;

  point_type p1(2.,1.,2.);
  point_type p2(2.,1.,2.);

  std::cout << "p1 = " << p1 << std::endl;
  std::cout << "p2 = " << p2 << std::endl;

  std::cout << "norm(p1) = " << p1.norm() << std::endl;
  std::cout << "p1.distance(p2) = " << p1.distance(p2) << std::endl;

  std::cout << "p1==p2 = " << (p1==p2) << std::endl;

  ////////////////////////////////////////////////////
  std::cout << "\n- Overall time " << omp_get_wtime( ) - t_0
            << " s" << std::endl;
  ////////////////////////////////////////////////////


  return 0;
}   

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

#ifndef plot_hpp
#define plot_hpp

// system
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdexcept>


/*!
	\class Plot
	
	The class Plot uses Gnuplot to plot particles. I'm sure, a smarter way of
	visualizing particles exists.
 */
template<typename particle_type>
class Plot
{
	static const unsigned int dim = particle_type::dim;
public:
	explicit Plot(const particle_type *const _X0,
								const unsigned int _N0)
		: X0(_X0), N0(_N0)
	{
    printf ("Checking if processor is available...");
    if (system(NULL)) puts("Ok");
    else exit (1);
  }
	
	~Plot() { }
	
	void plot() const;
	void plot(const particle_type *const X1, const unsigned int N1) const;
	void plot(const particle_type *const X1, const unsigned int N1,
						const particle_type *const X2, const unsigned int N2) const;


private:
	const particle_type *const X0;
	const unsigned int N0;
	
	void writeDataFile(const std::string file,
										 const particle_type *const X,
										 const unsigned int N) const;
};






template<typename particle_type>
void Plot<particle_type>::writeDataFile(const std::string file,
																				const particle_type *const X,
																				const unsigned int N) const
{
	std::ofstream dstr;
	dstr.open( file.c_str( ) );
	for (unsigned int i=0; i<N; ++i) {
		for (unsigned int d=0; d<dim; ++d)
			dstr << X[i].getPoint()[d] << "\t";
		dstr << std::endl;
	}
	dstr.close();
}






template<typename particle_type>
void Plot<particle_type>::plot() const
{
	// data files
	std::string data0("_data0");
	writeDataFile(data0, X0, N0);
	
	// control file
	std::string ctrl("_ctrl");
	std::ofstream cstr;
	cstr.open( ctrl.c_str() ); 	
	cstr << "set key off " << std::endl;
	cstr << "set size ratio 1" << std::endl;
	cstr << "set grid" << std::endl;

	switch (dim) {
	case 2: cstr << "plot \'";  break;
	case 3: cstr << "splot \'"; break;
	default: throw std::runtime_error("No valid dimension!"); }
	cstr << data0 << "\' w d" << std::endl;

	cstr << "pause -1 'Press any key to continue '"<< std::endl;
	cstr.close();
	
	// call plot
	const std::string gplot = "gnuplot ";
	system( (gplot + ctrl).c_str() ); 
	
	// remove files
	const std :: string rm = "rm ";
	system( (rm + data0 ).c_str() );
	system( (rm + ctrl  ).c_str() );
}






template<typename particle_type>
void Plot<particle_type>::plot(const particle_type *const X1,
															 const unsigned int N1) const
{
	// data files
	std::string data0("_data0");
	writeDataFile(data0, X0, N0);
	std::string data1("_data1");
	writeDataFile(data1, X1, N1);
	
	// control file
	std::string ctrl("_ctrl");
	std::ofstream cstr;
	cstr.open( ctrl.c_str() ); 	
	cstr << "set key off " << std::endl;
	cstr << "set size ratio 1" << std::endl;
	cstr << "set grid" << std::endl;

	switch (dim) {
	case 2: cstr << "plot \'";  break;
	case 3: cstr << "splot \'"; break;
	default: throw std::runtime_error("No valid dimension!"); }
	cstr << data0 << "\' w d, \'" << data1 << "\' w d" << std::endl;

	//	cstr << "plot \'" << data0 << "\', \'" << data1 << "\'" << std::endl;

	cstr << "pause -1 'Press any key to continue '"<< std::endl;
	cstr.close();
	
	// call plot
	const std::string gplot = "gnuplot ";
	system( (gplot + ctrl).c_str() ); 
	
	// remove files
	const std :: string rm = "rm ";
	system( (rm + data0 ).c_str() );
	system( (rm + data1 ).c_str() );
	system( (rm + ctrl  ).c_str() );
}






template<typename particle_type>
void Plot<particle_type>::plot(const particle_type *const X1,
															 const unsigned int N1,
															 const particle_type *const X2,
															 const unsigned int N2) const
{
	// data files
	std::string data0("_data0");
	writeDataFile(data0, X0, N0);
	std::string data1("_data1");
	writeDataFile(data1, X1, N1);
	std::string data2("_data2");
	writeDataFile(data2, X2, N2);
	
	// control file
	std::string ctrl("_ctrl");
	std::ofstream cstr;
	cstr.open( ctrl.c_str() ); 	
	cstr << "set key off " << std::endl;
	cstr << "set size ratio 1" << std::endl;
	cstr << "set grid" << std::endl;

	switch (dim) {
	case 2: cstr << "plot \'";  break;
	case 3: cstr << "splot \'"; break;
	default: throw std::runtime_error("No valid dimension!"); }
	cstr << data0 << "\' w d, \'" << data1 << "\' w d, \'" << data2 << "\' w d"
			 << std::endl;

	//	cstr << "plot \'"
	//			 << data0 << "\', \'"
	//			 << data1 << "\', \'"
	//			 << data2 << "\'" << std::endl;

	cstr << "pause -1 'Press any key to continue '" << std::endl;
	cstr << "quit" << std::endl;
	cstr.close();
	
	// call plot
	const std::string gplot = "gnuplot ";
	system( (gplot + ctrl).c_str() ); 
	
	// remove files
	const std :: string rm = "rm ";
	system( (rm + data0 ).c_str() );
	system( (rm + data1 ).c_str() );
	system( (rm + data2 ).c_str() );
	system( (rm + ctrl  ).c_str() );
}










#endif



//template<>
//void Plot<3>::plot() const
//{
//	// data files
//	std::string data0("_data0");
//	writeDataFile(data0, X0, N0);
//	
//	// control file
//	std::string ctrl("_ctrl");
//	std::ofstream cstr;
//	cstr.open( ctrl.c_str() ); 	
//	cstr << "set key off " << std::endl;
//	cstr << "set size ratio 1" << std::endl;
//	cstr << "set grid" << std::endl;
//	cstr << "splot \'" << data0 << "\' w d" << std::endl;
//	cstr << "pause -1 'Press any key to continue '"<< std::endl;
//	cstr.close();
//	
//	// call plot
//	const std::string gplot = "gnuplot ";
//	system( (gplot + ctrl).c_str() ); 
//	
//	// remove files
//	const std :: string rm = "rm ";
//	system( (rm + data0 ).c_str() );
//	system( (rm + ctrl  ).c_str() );
//}

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

#ifndef psdraw_hpp
#define psdraw_hpp

#include <fstream>


/*! \class DrawNearField
	This class is not directly needed within the FMM algorithm; it provides the
	functionality of drawing the near-field (dense blocks) in the system matrix
	to a postscript file.

	\brief Auxillary class draws near-field of system matrix to postscript file

	@tparam cluster_type type of cluster
	@tparam relations_type type of cluster-relation
 */
template <typename cluster_type,
          typename relations_type>
class DrawNearField
  : public std::unary_function<cluster_type,void>
{
  static const double origin = 45;
  static const double maxpts = 500;
	std::ofstream& os;
  const unsigned int N;
  const std::vector<relations_type>& relations;
  const double scale;
  const double x0;

public:
  explicit DrawNearField(std::ofstream& _os,
                         const unsigned int _N,
                         const std::vector<relations_type>& _relations)
    : os(_os), N(_N), relations(_relations),
      scale(maxpts/static_cast<double>(N)),
      x0(origin/scale)
  {
    assert(os.is_open());
    os << "%!PS-Adobe-2.0-2.0 EPSF-2.0" << std::endl;
    os << "%%BoundingBox: 0 0 " << maxpts << " " << maxpts << std::endl;
    os << scale << " " << scale << " scale" << std::endl; 
    os << "newpath" << std::endl;
    os << x0 << " " << x0+N << " moveto" << std::endl;
    os << x0 << " " << x0 << " lineto" << std::endl;
    os << x0+N << " " << x0 << " lineto" << std::endl;
    os << x0+N << " " << x0+N << " lineto" << std::endl;
    os << "closepath stroke" << std::endl;
  }

  ~DrawNearField()
  {
    os << "showpage" << std::endl;
  }

  void operator()(cluster_type *const clx) const
  {
    const unsigned int a = x0 + N - clx->nbeg;
    const unsigned int m = clx->size;

    BOOST_FOREACH(const cluster_type *const cly,
                  relations.at(clx->cidx).getNeighbors()) {
      const unsigned int b = x0 + cly->nbeg;
      const unsigned int n = cly->size;

      os << b << " " << a << " moveto" << std::endl;
      os << "0 -" << m << " rlineto" << std::endl;
      os << n << " 0 rlineto" << std::endl;
      os << "0 " << m << " rlineto" << std::endl;
      os << "closepath gsave 0.8 setgray fill grestore 0.5 setlinewidth 0 setgray stroke" << std::endl;
    }
  }
};




/*! \class DrawCrossSection
	This class is not directly needed within the FMM algorithm; it provides the
	functionality of drawing the cross section through a octree into a
	postscript file.

	\brief Auxillary class draws cross section of octree to postscript file

	@tparam cluster_type type of cluster
	@tparam level_type type of level
	@tparam m2l_handler_type type of m2l-handler
 */
template <typename cluster_type,
					typename level_type,
					typename m2l_handler_type>
class DrawCrossSection
  : public std::unary_function<cluster_type,void>
{
  enum {dim = cluster_type::dim};
  typedef typename cluster_type::point_type point_type;

  static const double origin = 40;
  static const double maxpts = 595;
	std::ofstream& os;
  const point_type rcenter;
  const double rext;
  const point_type loc;
  const unsigned int dir0;
	const std::vector<level_type*>& levels;
	const std::vector<m2l_handler_type>& m2lv;
  const double scale;
  const double lw;
  const double rc;

public:
  explicit DrawCrossSection(std::ofstream& _os,
                            const point_type _rcenter,
                            const double _rext,
                            const point_type _loc,
                            const unsigned int _dir0,
														const std::vector<level_type*>& _levels,
														const std::vector<m2l_handler_type>& _m2lv)
    : os(_os),
      rcenter(_rcenter), rext(_rext), loc(_loc), dir0(_dir0), levels(_levels),
			m2lv(_m2lv), scale((maxpts-2*origin)/rext),
//      lw(2/scale),
      lw(.5/scale),
      rc(maxpts/(2*scale))
  {
    assert(os.is_open());
    os << "%!PS-Adobe-2.0-2.0 EPSF-2.0" << std::endl;
    os << "%%BoundingBox: 0 0 " << maxpts << " " << maxpts << std::endl;
    os << scale << " " << scale << " scale" << std::endl; 
    os << lw << " setlinewidth" << std::endl;
    os << "newpath" << std::endl;
    os << rc-rext/2 << " " << rc-rext/2 << " moveto" << std::endl;
    os << rc+rext/2 << " " << rc-rext/2 << " lineto" << std::endl;
    os << rc+rext/2 << " " << rc+rext/2 << " lineto" << std::endl;
    os << rc-rext/2 << " " << rc+rext/2 << " lineto" << std::endl;
    os << "closepath stroke" << std::endl;
  }

//  void operator()(cluster_type *const cl) const
//  {
//    const unsigned int dir1 = (dir0+2) % dim;
//    const unsigned int dir2 = (dir0+1) % dim;
//    const point_type& center = cl->getCenter();
//    const double ext = cl->getLevel()->getExtension();
//
//    const double lw_ = lw / (cl->getLevel()->getNLevel()+1);
//
//    if ( ((center[dir0]-loc[dir0]) <= ext/2) &&
//         ((center[dir0]-loc[dir0]) > -ext/2) ) {
//
//      const double a = center[dir1] - rcenter[dir1] - ext/2;
//      const double b = center[dir2] - rcenter[dir2] - ext/2;
//      os << rc+a     << " " << rc+b     << " moveto" << std::endl;
//      os << rc+a+ext << " " << rc+b     << " lineto" << std::endl;
//      os << rc+a+ext << " " << rc+b+ext << " lineto" << std::endl;
//      os << rc+a     << " " << rc+b+ext << " lineto" << std::endl;
//      os << "closepath " << lw_ << " setlinewidth stroke" << std::endl;
//
//      if (cl->isleaf()) {
//        os << "/Times-Roman findfont" << std::endl;
//        os << ext*0.8 << " scalefont" << std::endl;
//        os << "setfont" << std::endl;
//        os << rc+a+ext*.05 << " " << rc+b+ext*.1 << " moveto" << std::endl;
//        os << "(" << cl->getsize() << ") show" << std::endl;
//      }
//
//    }
//  }

  void operator()(cluster_type *const cl) const
  {
    const unsigned int dir1 = (dir0+2) % dim;
    const unsigned int dir2 = (dir0+1) % dim;
    const point_type& center = cl->getCenter();
		const level_type *const level = levels.at(cl->getNumLevel());
    const double ext = level->getExtension();

		if ( ((center[dir0]-loc[dir0]) <= ext/2) &&
         ((center[dir0]-loc[dir0]) > -ext/2) ) {

      const double a = center[dir1] - rcenter[dir1] - ext/2;
      const double b = center[dir2] - rcenter[dir2] - ext/2;
      os << rc+a     << " " << rc+b     << " moveto" << std::endl;
      os << rc+a+ext << " " << rc+b     << " lineto" << std::endl;
      os << rc+a+ext << " " << rc+b+ext << " lineto" << std::endl;
      os << rc+a     << " " << rc+b+ext << " lineto" << std::endl;
			
			if (!cl->isLeaf() && m2lv.at(cl->getNumLevel()).size()>0)
				os << "closepath gsave 0.9 setgray fill grestore 0 setgray stroke"
					 << std::endl;
			else if (cl->isLeaf())
				os << "closepath gsave 0.6 setgray fill grestore 0 setgray stroke"
					 << std::endl;
			else
				os << "closepath stroke" << std::endl;

//      if (cl->isleaf()) {
//        os << "/Times-Roman findfont" << std::endl;
//        os << ext*0.8 << " scalefont" << std::endl;
//        os << "setfont" << std::endl;
//        os << rc+a+ext*.05 << " " << rc+b+ext*.1 << " moveto" << std::endl;
//        os << "(" << cl->getsize() << ") show" << std::endl;
//      }

    }
  }

};


#include <boost/lexical_cast.hpp>


template <typename cluster_type,
					typename level_type,
					typename m2l_handler_type>
void drawCrossSection(const std::string& filename,
											const typename cluster_type::point_type location,
											const unsigned int direction,
											const typename cluster_type::point_type center,
											const double extension,
											const std::vector<std::vector<cluster_type*> >& clvv,
											const std::vector<level_type*>& levels,
											const std::vector<m2l_handler_type>& m2lv)
{
	typedef DrawCrossSection<cluster_type,level_type,m2l_handler_type>
		drawer_type;

  std::ofstream psout
		((filename+"_d"+boost::lexical_cast<std::string>(direction)+".ps").c_str());
  drawer_type draw(psout, center, extension, location, direction, levels, m2lv);
  BOOST_FOREACH(const std::vector<cluster_type*>& clv, clvv)
    std::for_each(clv.begin(), clv.end(), draw);
  psout.close();
}




#endif

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

#ifndef storage_tpl
#define storage_tpl



template <int DIM>
Storage<DIM>::Storage(const std::pair<point_type,point_type> bbx,
                      const unsigned int lmax)
{
  // set root cluster //////////////////
  root_cluster_diam = 0.;
  for (unsigned int d=0; d<dim; ++d)
    if (root_cluster_diam < bbx.second[d]-bbx.first[d])
      root_cluster_diam   = bbx.second[d]-bbx.first[d];
  root_cluster_center = (bbx.first+bbx.second) / 2.;
  
  // set levels /////////////
  double diam = root_cluster_diam;
  for (unsigned int l=0; l<=lmax; ++l) {
    const bool is_leaf = (l == lmax ? true : false);
    levels.push_back(new level_type(l, diam, is_leaf));
    diam /= 2.;
  }
}




#endif

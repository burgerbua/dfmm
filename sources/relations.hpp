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

#ifndef relations_hpp
#define relations_hpp



template <int> class Cluster;


template <int DIM>
class Relations
{
  template <int,int,typename,typename> friend class M2Lapplyer;

  typedef Cluster<DIM>                                cluster_type;
	typedef std::vector<const cluster_type*>            cluster_vector;
  typedef std::pair<unsigned int,const cluster_type*> cluster_pair;
  typedef std::vector<cluster_pair>                   cluster_pair_vector;
  typedef std::map<unsigned int,cluster_pair_vector>  cluster_pair_vector_map;

	cluster_vector             neighbor_list;
  cluster_pair_vector_map interaction_lists;

public:
  explicit Relations()
  {
    neighbor_list.clear();
    interaction_lists.clear();
  }

  void push_back_neighbor(const cluster_type* cl)
  { neighbor_list.push_back(cl); }

  void push_back_interaction(const cluster_type* cl,
                             const unsigned int t_id,
                             const unsigned int c)
  {
    interaction_lists.insert(std::make_pair(c, cluster_pair_vector()));
    interaction_lists[c].push_back(cluster_pair(t_id, cl));
  }

  const cluster_vector& getNeighbors() const
  { return neighbor_list; }

  const cluster_pair_vector_map& getInteractions() const
  { return interaction_lists; }
};


#endif


#include "base/logger.hpp"
#include "base/change_type.hpp"
#include "base/fasta_reader.hpp"
#include "base/residue_seq.hpp"
#include "base/proteoform_factory.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph.hpp"

#include <boost/timer.hpp>

namespace prot {

ProteoGraph::ProteoGraph(FastaSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_ptr_vec,
                         MassGraphPtr graph_ptr, bool is_nme,
                         double convert_ratio, int max_mod_num,
                         int max_ptm_sum_mass) {
  db_proteo_ptr_ = ProteoformFactory::geneDbProteoformPtr(fasta_seq_ptr, fix_mod_ptr_vec);
  graph_ptr_ = graph_ptr;
  is_nme_ = is_nme;
  node_num_ = num_vertices(*graph_ptr.get());
  LOG_DEBUG("node num " << node_num_);
  pair_num_ = node_num_ * (node_num_ + 1) /2;
  compSeqMasses(convert_ratio);
  compDistances(convert_ratio, max_mod_num, max_ptm_sum_mass);
}

int ProteoGraph::getVecIndex(int v1, int v2) {
  int index =  (node_num_ + node_num_ - v1 + 1) * v1 /2  + (v2 - v1);
  return index;
}

int ProteoGraph::getSeqMass(int v1, int v2) {
  int index = getVecIndex(v1, v2);
  return seq_masses_[index];
}

void ProteoGraph::compSeqMasses(double convert_ratio) {
  ResSeqPtr res_seq_ptr = db_proteo_ptr_->getResSeqPtr();
  seq_masses_= std::vector<int>(pair_num_, 0);
  for (int i = 0; i < node_num_; i ++) {
    int mass = 0;
    for (int j = i+1; j < node_num_; j++) {
      int cur_mass 
          = std::round(res_seq_ptr->getResiduePtr(j-1)->getMass() * convert_ratio) ;
      mass += cur_mass;
      int index = getVecIndex(i, j);
      seq_masses_[index] = mass;
    }
  }
}

void ProteoGraph::compDistances(double convert_ratio, int max_mod_num, int max_ptm_sum_mass) {
  boost::timer t;

  MassGraph *g_p = graph_ptr_.get();
  //get mass without ptms

  std::vector<std::vector<std::set<int>>> dist_vecs; 
  for (int i = 0; i < pair_num_; i++) {; 
    std::set<int> empty_set;
    std::vector<std::set<int>> one_pair_vec;
    for (int j = 0; j < max_mod_num + 1; j ++) {
      one_pair_vec.push_back(empty_set);
    }
    dist_vecs.push_back(one_pair_vec);
  }
  // initialize pair (i, i)
  for (int i = 0; i < node_num_; i++) {
    int index = getVecIndex(i, i);
    dist_vecs[index][0].insert(0);
  }
  for (int i = 0; i < node_num_ - 1; i++) {
    for (int j = i+1; j < node_num_; j++) {
      Vertex v2 = vertex(j, *g_p);
      Vertex pre_v2 = vertex(j-1, *g_p);
      int index = getVecIndex(i, j);
      int pre_index = getVecIndex(i, j-1);
      boost::graph_traits<MassGraph>::out_edge_iterator ei, ei_end;
      boost::tie(ei, ei_end) = out_edges(pre_v2, *g_p);
      for( ; ei != ei_end; ++ei) {
        if(target(*ei, *g_p) == v2 ) {
          MassGraph::edge_descriptor e = *ei;
          int d =(*g_p)[e].int_mass_;
          int change = (*g_p)[e].change_type_;
          for (int k = 0; k < max_mod_num + 1; k++) {
            if (k == max_mod_num && 
                (change == ChangeType::PROTEIN_VARIABLE->getId() || change == ChangeType::VARIABLE->getId())) {
              continue;
            }
            for (std::set<int>::iterator it=dist_vecs[pre_index][k].begin(); 
                 it!=dist_vecs[pre_index][k].end(); it++) {
              int new_d = d + *it;
              if (std::abs(new_d - seq_masses_[index]) <= max_ptm_sum_mass) {
                if (change == ChangeType::PROTEIN_VARIABLE->getId() || change == ChangeType::VARIABLE->getId()) {
                  dist_vecs[index][k+1].insert(new_d);
                } else {
                  dist_vecs[index][k].insert(new_d);
                }
              }
            }
          }
        }
      }
    }
  }

  DistVec tmp;
  for (int k = 0; k < max_mod_num + 1; k++) {
    dist_vec_.push_back(tmp);
  }

  for (int k = 0; k < max_mod_num + 1; k++) {
    addToDistVec(graph_ptr_, dist_vecs, node_num_, k, dist_vec_[k]);
  }
  
  std::cout << "compDistances time: " << t.elapsed() << std::endl;
}

}


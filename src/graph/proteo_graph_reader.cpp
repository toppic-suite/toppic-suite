#include "base/logger.hpp"
#include "graph/proteo_graph_reader.hpp"

namespace prot {

ProteoGraphReader::ProteoGraphReader(const std::string &db_file_name,
                                     const ModPtrVec &fix_mod_ptr_vec, 
                                     const ProtModPtrVec &prot_mod_ptr_vec,
                                     const ModPtrVec &var_mod_ptr_vec,
                                     double convert_ratio, int max_mod_num,
                                     int max_ptm_sum_mass, int proteo_graph_gap) {
  fix_mod_ptr_vec_ = fix_mod_ptr_vec;
  convert_ratio_ = convert_ratio;
  max_mod_num_ = max_mod_num;
  max_ptm_sum_mass_ = max_ptm_sum_mass;
  reader_ptr_ = FastaReaderPtr(new FastaReader(db_file_name));
  proteo_anno_ptr_ = ProteoAnnoPtr(
      new ProteoAnno(fix_mod_ptr_vec, prot_mod_ptr_vec, var_mod_ptr_vec));
  proteo_graph_gap_ = proteo_graph_gap;
}

MassGraphPtr ProteoGraphReader::getMassGraphPtr() {
  MassGraphPtr graph_ptr = MassGraphPtr(new MassGraph());
  int seq_len  = proteo_anno_ptr_->getLen();
  for (int i = 0; i < seq_len + 1; i++) {
    VertexInfo v(i);
    add_vertex(v, *graph_ptr.get());
  }

  for (int i = 0; i < seq_len; i++) {
    Vertex v1, v2;
    v1 = vertex(i, *graph_ptr.get());
    v2 = vertex(i+1, *graph_ptr.get());
    ResiduePtrVec res_ptr_vec = proteo_anno_ptr_->getResiduePtrVec(i);
    std::vector<int> change_vec = proteo_anno_ptr_->getChangeVec(i);
    for (size_t j=0; j < res_ptr_vec.size(); j++) { 
      EdgeInfo edge_info(res_ptr_vec[j], change_vec[j], convert_ratio_);
      add_edge(v1, v2, edge_info , *graph_ptr.get());
    }
  }
  return graph_ptr;
}

ProteoGraphPtr ProteoGraphReader::getNextProteoGraphPtr() {
  FastaSeqPtr seq_ptr = reader_ptr_->getNextSeq();
  if (seq_ptr.get() == nullptr) {
    return ProteoGraphPtr(nullptr);
  }
  LOG_DEBUG("name " << seq_ptr->getName() << " seq " << seq_ptr->getRawSeq());
  proteo_anno_ptr_->anno(seq_ptr->getRawSeq());
  MassGraphPtr graph_ptr = getMassGraphPtr(); 

  return ProteoGraphPtr(new ProteoGraph(seq_ptr, fix_mod_ptr_vec_, graph_ptr, 
                                        proteo_anno_ptr_->isNme(),
                                        convert_ratio_, max_mod_num_,
                                        max_ptm_sum_mass_, proteo_graph_gap_));
}

ProteoGraphPtr ProteoGraphReader::getProteoGraphPtrBySeq(FastaSeqPtr seq_ptr) {
  proteo_anno_ptr_->anno(seq_ptr->getRawSeq());
  MassGraphPtr graph_ptr = getMassGraphPtr(); 
  return ProteoGraphPtr(new ProteoGraph(seq_ptr, fix_mod_ptr_vec_, graph_ptr, 
                                        proteo_anno_ptr_->isNme(),
                                        convert_ratio_, max_mod_num_,
                                        max_ptm_sum_mass_, proteo_graph_gap_));	
}

}


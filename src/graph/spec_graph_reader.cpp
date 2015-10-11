#include "graph/spec_graph_reader.hpp"

namespace prot {

SpecGraphReader::SpecGraphReader(const std::string &sp_file_name,
                                 int group_sp_num,
                                 double convert_ratio,
                                 SpParaPtr sp_para_ptr) {
  ms_reader_ptr_ 
      = MsAlignReaderPtr(new MsAlignReader(sp_file_name, group_sp_num));
  group_sp_num_ = group_sp_num;
  convert_ratio_ = convert_ratio;
  sp_para_ptr_ = sp_para_ptr;
}

MassGraphPtr SpecGraphReader::getMassGraphPtr(const PrmPeakPtrVec &peak_vec) {

  LOG_DEBUG("start mass graph");
  MassGraphPtr graph_ptr = MassGraphPtr(new MassGraph());
  
  //*graph_ptr.get()[boost::graph_bundle].name_ = spec_set_ptr->getDeconvMsPtr()
  //    ->getHeaderPtr()->getScansString();
  //LOG_DEBUG("spec_graph name " << spec_graph[boost::graph_bundle].name_);
  //PrmMsPtrVec ms_six_vec = spec_set_ptr->getMsSixPtrVec();
  //PrmPeakPtrVec peak_vec = getPrmPeakPtrs(ms_six_vec, sp_para_ptr_->getPeakTolerancePtr());

  // add mass 0/start nod
  VertexInfo v(0);
  add_vertex(v, *graph_ptr.get());

  for (size_t i = 1; i < peak_vec.size(); i++) {
    // add a new node for the prm
    VertexInfo cur_v(i);
    add_vertex(cur_v, *graph_ptr.get());

    Vertex v1, v2;
    v1 = vertex(i-1, *graph_ptr.get());
    v2 = vertex(i, *graph_ptr.get());

    double dist = peak_vec[i]->getMonoMass() - peak_vec[i-1]->getMonoMass();

    EdgeInfo edge_info(dist, convert_ratio_);
    add_edge(v1, v2, edge_info , *graph_ptr.get());
  }

  return graph_ptr;
}

SpecGraphPtr SpecGraphReader::getNextSpecGraphPtr() {
  SpectrumSetPtr spec_set_ptr = ms_reader_ptr_->getNextSpectrumSet(sp_para_ptr_);
  LOG_DEBUG("get spec set ");
  if (spec_set_ptr  == nullptr) {
    return SpecGraphPtr(nullptr);
  }
  if (spec_set_ptr->isValid()) {
    PrmMsPtrVec ms_six_vec = spec_set_ptr->getMsSixPtrVec();
    PrmPeakPtrVec peak_vec = getPrmPeakPtrs(ms_six_vec, sp_para_ptr_->getPeakTolerancePtr());
    MassGraphPtr graph_ptr = getMassGraphPtr(peak_vec); 
    LOG_DEBUG("graph complete");
    return SpecGraphPtr(new SpecGraph(spec_set_ptr, peak_vec, graph_ptr, convert_ratio_));
  }
  else {
    return SpecGraphPtr(new SpecGraph(spec_set_ptr));
  }
}

}


#include "graph/spec_graph_reader.hpp"

namespace prot {

SpecGraphReader::SpecGraphReader(const std::string &sp_file_name,
                                 int group_sp_num,
                                 double convert_ratio,
                                 SpParaPtr sp_para_ptr) {
  ms_reader_ptr_ 
      = MsAlignReaderPtr(new MsAlignReader(sp_file_name, group_sp_num,
                                           sp_para_ptr->getActivationPtr()));
  group_sp_num_ = group_sp_num;
  convert_ratio_ = convert_ratio;
  sp_para_ptr_ = sp_para_ptr;
}

MassGraphPtr SpecGraphReader::getMassGraphPtr(const PrmPeakPtrVec &peak_vec) {

  LOG_DEBUG("start mass graph");
  MassGraphPtr graph_ptr = MassGraphPtr(new MassGraph());

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

SpecGraphPtrVec SpecGraphReader::getNextSpecGraphPtrVec(int error) {
  SpectrumSetPtr spec_set_ptr = ms_reader_ptr_->getNextSpectrumSet(sp_para_ptr_);
  LOG_DEBUG("get spec set ");
  SpecGraphPtrVec graph_ptr_vec;
  if (spec_set_ptr  == nullptr) {
    return graph_ptr_vec;
  }
  std::vector<double> prec_errors;
  prec_errors.push_back(0);
  for (int i = 1; i <= error; i++) {
    prec_errors.push_back(- i * MassConstant::getIsotopeMass());
    prec_errors.push_back(i * MassConstant::getIsotopeMass());
  }

  DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
  //LOG_DEBUG("deconv ms size " << deconv_ms_ptr_vec.size());
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  LOG_DEBUG("prec_mono_mass  " << prec_mono_mass);
  if (spec_set_ptr->isValid()) {
    LOG_DEBUG("valid");
    for (size_t i = 0; i < prec_errors.size(); i++) {
      SpectrumSetPtr adjusted_spec_set_ptr(
          new SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr_, prec_mono_mass + prec_errors[i]));
      PrmMsPtrVec ms_six_vec = adjusted_spec_set_ptr->getMsSixPtrVec();
      PrmPeakPtrVec peak_vec = PrmMs::getPrmPeakPtrs(ms_six_vec, sp_para_ptr_->getPeakTolerancePtr());
      MassGraphPtr graph_ptr = getMassGraphPtr(peak_vec); 
      LOG_DEBUG("graph complete");
      SpecGraphPtr spec_graph_ptr = SpecGraphPtr(new SpecGraph(adjusted_spec_set_ptr, peak_vec, graph_ptr, convert_ratio_));
      graph_ptr_vec.push_back(spec_graph_ptr);
    }
  }
  else {
    LOG_DEBUG("no valid");
    SpecGraphPtr spec_graph_ptr = SpecGraphPtr(new SpecGraph(spec_set_ptr));
    graph_ptr_vec.push_back(spec_graph_ptr);
  }
  LOG_DEBUG("set geneneted");
  return graph_ptr_vec;
}

}


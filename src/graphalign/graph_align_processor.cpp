#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph_reader.hpp"
#include "graph/spec_graph_reader.hpp"
#include "graphalign/graph_align.hpp"
#include "graphalign/graph_align_processor.hpp"
#include "threadpool.hpp"

#define NUM_THREAD 4

namespace prot {

GraphAlignProcessor::GraphAlignProcessor(GraphAlignMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
}

void GraphAlignProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;
  LOG_DEBUG("start reading " << var_mod_file_name);
  ModPtrVec var_mod_ptr_vec = ModUtil::readModTxt(var_mod_file_name)[2];

  LOG_DEBUG("end reading " << var_mod_file_name);

  ProteoGraphReader reader(db_file_name, 
                           prsm_para_ptr->getFixModPtrVec(),
                           prsm_para_ptr->getProtModPtrVec(),
                           var_mod_ptr_vec,
                           mng_ptr_->convert_ratio_,
                           mng_ptr_->max_known_mods_,
                           mng_ptr_->getIntMaxPtmSumMass(),
                           mng_ptr_->proteo_graph_gap_);
  LOG_DEBUG("init reader complete");
  ProteoGraphPtrVec proteo_ptrs;
  ProteoGraphPtr proteo_ptr;
  int count = 0;
  while ((proteo_ptr = reader.getNextProteoGraphPtr()) != nullptr) {
    count++;
    proteo_ptrs.push_back(proteo_ptr);
  }

  LOG_DEBUG("Prot graph number " << count);
  std::string output_file_name = FileUtil::basename(sp_file_name)+"."+mng_ptr_->output_file_ext_;
  PrsmXmlWriter prsm_writer(output_file_name);

  LOG_DEBUG("start init spec reader");
  SpecGraphReader spec_reader(sp_file_name, 
                              prsm_para_ptr->getGroupSpecNum(),
                              mng_ptr_->convert_ratio_,
                              sp_para_ptr);
  LOG_DEBUG("init spec reader complete");

  SpecGraphPtrVec spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr_->prec_error_);
  LOG_DEBUG("spec ptr reading complete");
  LOG_DEBUG("spec_ptr_vec " << spec_ptr_vec.size());
  int spectrum_num = MsAlignUtil::getSpNum (prsm_para_ptr->getSpectrumFileName());
  int sp_count = 0;
  ThreadPool pool(NUM_THREAD);
  while (spec_ptr_vec.size() != 0) {
    sp_count++;
    //std::cout << std::flush << "Mass graph is processing " << sp_count << " of " << spectrum_num << " spectra.\r";
    //LOG_DEBUG("spectrum id " << spec_ptr_vec[0]->getSpectrumSetPtr()->getSpecId());
    spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr_->prec_error_);
    pool.Enqueue([sp_count](){
                 std::cout << "Processed: " << sp_count<< std::endl; 
                 });
    spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr_->prec_error_);
  }
  prsm_writer.close();
}

}


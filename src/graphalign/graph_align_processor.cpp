#include <iostream>
#include <fstream>

#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "base/threadpool.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph_reader.hpp"
#include "graph/spec_graph_reader.hpp"
#include "graphalign/graph_align.hpp"
#include "graphalign/graph_align_processor.hpp"

namespace prot {


GraphAlignProcessor::GraphAlignProcessor(GraphAlignMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
}

std::function<void()> geneTask(GraphAlignMngPtr mng_ptr,
                               ProteoGraphPtr proteo_ptr,
                               SpecGraphPtr spec_ptr,
                               ThreadPoolPtr pool_ptr) {
  return [mng_ptr, proteo_ptr, spec_ptr, pool_ptr]() {
    GraphAlignPtr graph_align 
        = std::make_shared<GraphAlign>(mng_ptr, proteo_ptr, spec_ptr);
    graph_align->process();
    boost::thread::id thread_id = boost::this_thread::get_id();
    PrsmXmlWriterPtr writer_ptr = pool_ptr->getWriter(thread_id);
    for (int shift = 0; shift <= mng_ptr->n_unknown_shift_; shift++) {
      PrsmPtr prsm_ptr = graph_align->geneResult(shift);
      if (prsm_ptr != nullptr) {
        writer_ptr->write(prsm_ptr);
      }
    }
    graph_align = nullptr;
  };
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

  int spectrum_num = MsAlignUtil::getSpNum (prsm_para_ptr->getSpectrumFileName());
  std::string base_file_name = FileUtil::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  //ProteoGraphPtrVec proteo_ptrs;
  ProteoGraphPtr proteo_ptr;
  int proteo_count = 0;
  while ((proteo_ptr = reader.getNextProteoGraphPtr()) != nullptr) {
    std::string output_file_name = base_file_name + "_" + std::to_string(proteo_count);
    int sp_count = 0;
    ThreadPoolPtr pool_ptr(new ThreadPool(mng_ptr_->thread_num_, output_file_name));
    SpecGraphReader spec_reader(sp_file_name, 
                                prsm_para_ptr->getGroupSpecNum(),
                                mng_ptr_->convert_ratio_,
                                sp_para_ptr);
    LOG_DEBUG("init spec reader complete");
    SpecGraphPtrVec spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr_->prec_error_);
    LOG_DEBUG("spec ptr reading complete");
    LOG_DEBUG("spec_ptr_vec " << spec_ptr_vec.size());
    while (spec_ptr_vec.size() != 0) {
      sp_count++;
      LOG_DEBUG("spectrum id " << spec_ptr_vec[0]->getSpectrumSetPtr()->getSpecId());
      std::cout << std::flush << "Mass graph is processing protein " << proteo_count 
          << " spectrum " << sp_count << " of " << spectrum_num << " spectra.\r";
      for (size_t i = 0; i < spec_ptr_vec.size(); i++) {
        if (spec_ptr_vec[i]->getSpectrumSetPtr()->isValid()) {
          while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
            boost::this_thread::sleep( boost::posix_time::milliseconds(100) );
          }
          pool_ptr->Enqueue(geneTask(mng_ptr_, proteo_ptr,  spec_ptr_vec[i], pool_ptr));
        }
      }
      spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr_->prec_error_);
    }

    pool_ptr->ShutDown();
    std::cout << std::endl;

    // combine result files
    std::vector<std::string> input_exts;
    std::string cur_output_ext = mng_ptr_->output_file_ext_ + "_" + std::to_string(proteo_count);
    for (int i = 0; i < mng_ptr_->thread_num_; i++) {
      std::string fname = cur_output_ext + "_" + std::to_string(i); 
      //std::cout << "file name " << fname << std::endl;
      input_exts.push_back(fname);
    }
    int top_num = 1;
    PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts,  cur_output_ext, top_num));
    combine_ptr->process();
    combine_ptr = nullptr;
    proteo_count++;
  }

  std::vector<std::string> input_exts;
  for (int i = 0; i < proteo_count; i++) {
    std::string fname = mng_ptr_->output_file_ext_ + "_" + std::to_string(i); 
    input_exts.push_back(fname);
  }
  int top_num = 1;
  PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, mng_ptr_->output_file_ext_, top_num));
  combine_ptr->process();
  combine_ptr = nullptr;
}

}


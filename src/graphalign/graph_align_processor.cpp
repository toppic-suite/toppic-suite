#include <iostream>
#include <fstream>

#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph_reader.hpp"
#include "graph/spec_graph_reader.hpp"
#include "graphalign/graph_align.hpp"
#include "graphalign/graph_align_processor.hpp"
#include "threadpool.hpp"

#define NUM_THREAD 4

namespace prot {

boost::mutex mtx;

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

  FileUtil::createFolder("mass_graph");
  std::string output_file_name = "mass_graph" + FileUtil::getFileSeparator() + FileUtil::basename(sp_file_name);

  MsAlignReader ms_reader(sp_file_name, prsm_para_ptr->getGroupSpecNum(), sp_para_ptr->getActivationPtr());
  std::vector<std::string> ms_spec = ms_reader.readOneSpectrum();
  int sp_count = 0;
  while(ms_spec.size() > 0) {
    sp_count++;
    std::ofstream spec_file;
    spec_file.open(output_file_name + "_" + std::to_string(sp_count) + ".msalign");
    for (size_t i = 0; i < ms_spec.size(); i++) {
      spec_file << ms_spec[i] << std::endl;
    }
    spec_file.close();
    ms_spec = ms_reader.readOneSpectrum();
  }

  int spectrum_num = MsAlignUtil::getSpNum (prsm_para_ptr->getSpectrumFileName());

  ThreadPool pool(NUM_THREAD);
  GraphAlignMngPtr mng_ptr = mng_ptr_;

  for (sp_count = 0; sp_count <= spectrum_num; sp_count++) {
    pool.Enqueue([&proteo_ptrs, &mng_ptr, &prsm_para_ptr, &sp_para_ptr, output_file_name, sp_count, spectrum_num](){
      PrsmXmlWriter prsm_writer(output_file_name + "_" + std::to_string(sp_count) + "." + mng_ptr->output_file_ext_);
      mtx.lock();
      std::cout << std::flush << "Mass graph is processing " << sp_count << " of " << spectrum_num << " spectra.\r";
      mtx.unlock();
      SpecGraphReader spec_reader(output_file_name + "_" + std::to_string(sp_count) + ".msalign",
                                  prsm_para_ptr->getGroupSpecNum(),
                                  mng_ptr->convert_ratio_,
                                  sp_para_ptr);
      SpecGraphPtrVec spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr->prec_error_);
      for (size_t spec = 0; spec < spec_ptr_vec.size(); spec++) {
        if (spec_ptr_vec[0]->getSpectrumSetPtr()->isValid()) {
          for (size_t i = 0; i < proteo_ptrs.size(); i++) {
            GraphAlignPtr graph_align 
              = std::make_shared<GraphAlign>(mng_ptr, proteo_ptrs[i], spec_ptr_vec[0]);
            graph_align->process();
            PrsmPtr prsm_ptr = graph_align->geneResult(0);
            if (prsm_ptr != nullptr) {
              prsm_writer.write(prsm_ptr);
            } 
          }
        }
      }
      prsm_writer.close();
    });
  }

  pool.ShutDown();
  std::cout << std::endl;
  FastaIndexReaderPtr seq_reader(new FastaIndexReader(db_file_name));
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
  PrsmPtrVec prsm_vec;
  for (int i = 0; i <= spectrum_num; i++) {
    std::string fname = output_file_name + "_" + std::to_string(i) + "." + mng_ptr_->output_file_ext_; 
    PrsmReader prsm_reader(fname);
    PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
    while (prsm_ptr != nullptr) {
      prsm_vec.push_back(prsm_ptr);
      prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
    }
  }
  PrsmXmlWriter prsm_writer(FileUtil::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_);
  for (size_t i = 0; i < prsm_vec.size(); i++) {
    prsm_writer.write(prsm_vec[i]);
  }
  prsm_writer.close();
}

}


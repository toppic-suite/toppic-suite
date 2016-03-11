#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph_reader.hpp"
#include "graph/spec_graph_reader.hpp"
#include "graphalign/graph_align.hpp"
#include "graphalign/graph_align_processor.hpp"

namespace prot {

GraphAlignProcessor::GraphAlignProcessor(
    GraphAlignMngPtr mng_ptr) {
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

  FastaReader fasta_reader(db_file_name);
  std::vector<FastaSeqPtr> fasta_vec;
  FastaSeqPtr fasta_ptr;
  LOG_DEBUG("init reader complete");
  int count = 0;
  while ((fasta_ptr = fasta_reader.getNextSeq()) != nullptr) {
    count++;
    fasta_vec.push_back(fasta_ptr);
  }

  LOG_DEBUG("Prot graph number " << count);
  std::string output_file_name = FileUtil::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;
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
  while (spec_ptr_vec.size() != 0) {
    sp_count++;
    LOG_DEBUG("spectrum id " << spec_ptr_vec[0]->getSpectrumSetPtr()->getSpecId());
    std::cout << std::flush << "Mass graph is processing " << sp_count 
        << " of " << spectrum_num << " spectra.\r";
    for (size_t i = 0; i < fasta_vec.size(); i++) {
      ProteoGraphPtr proteo_ptr = reader.getProteoGraphPtrBySeq(fasta_vec[i]);

      for (size_t spec = 0; spec < spec_ptr_vec.size(); spec++) {
        if (spec_ptr_vec[spec]->getSpectrumSetPtr()->isValid() && proteo_ptr != nullptr) {
          GraphAlignPtr graph_align 
              = std::make_shared<GraphAlign>(mng_ptr_, proteo_ptr , spec_ptr_vec[spec]);
          graph_align->process();
          LOG_DEBUG("align process complete");
          for (int shift = 0; shift <= mng_ptr_->n_unknown_shift_; shift++) {
            PrsmPtr prsm_ptr = graph_align->geneResult(shift);
            if (prsm_ptr != nullptr) {
              prsm_writer.write(prsm_ptr);
            }
          }
          graph_align = nullptr;
        }
      }
    }
    spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(mng_ptr_->prec_error_);
  }
  prsm_writer.close();
}

}


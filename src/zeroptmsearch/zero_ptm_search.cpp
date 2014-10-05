#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

namespace prot {

void zeroPtmSearch(SpectrumSetPtr spec_set_ptr, 
                   SemiAlignTypePtr type_ptr,
                   ProteoformPtrVec &proteoform_ptr_vec, 
                   ZeroPtmMngPtr mng_ptr, 
                   PrsmPtrVec &prsms) {
  ExtendMsPtr ms_three = spec_set_ptr->getMsThreePtr();

  ZpFastMatchPtrVec fast_matches 
      = zeroPtmFastFilter(type_ptr, ms_three, proteoform_ptr_vec, 
                          mng_ptr->zero_ptm_filter_result_num_);

  //LOG_DEBUG("fast_match ended size " << fast_matches.size());
  DeconvMsPtr deconv_ms = spec_set_ptr->getDeconvMsPtr();
  ZpSlowMatchPtrVec slow_matches 
      = zeroPtmSlowFilter(deconv_ms, fast_matches, mng_ptr); 

  //LOG_DEBUG("slow_match ended size " << slow_matches.size());
  for (size_t i = 0; i < slow_matches.size(); i++) {
      prsms.push_back(slow_matches[i]->geneResult());
  }
  //LOG_DEBUG("prsm generation ended size " << prsms.size());

  std::sort(prsms.begin(), prsms.end(), prsmMatchFragmentDown);
  if (prsms.size() > 0) {
    prsms.erase(prsms.begin() + 1, prsms.end());
  }
}

void zeroPtmSearchProcess(ZeroPtmMngPtr mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr->getSearchDbFileName(), 
                              prsm_para_ptr->getFixModResiduePtrVec());

  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, prsm_para_ptr->getAllowProtModPtrVec());

  int spectra_num = countSpNum (prsm_para_ptr->getSpectrumFileName());
  //LOG_DEBUG("spectra_number  " << spectra_num);

  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName());
  std::string output_file_name = basename(prsm_para_ptr->getSpectrumFileName())
                                          + "." + mng_ptr->output_file_ext_;

  std::string log_file_name = prsm_para_ptr->getLogFileName();

  std::ofstream logfile;

  if (log_file_name.length() != 0){
    logfile.open(log_file_name, std::ios::out | std::ios::app);
  }

  PrsmWriter comp_writer(output_file_name + "_" 
                         + SemiAlignTypeFactory::getCompletePtr()->getName());
  PrsmWriter pref_writer(output_file_name + "_"
                         + SemiAlignTypeFactory::getPrefixPtr()->getName());
  PrsmWriter suff_writer(output_file_name + "_"
                         + SemiAlignTypeFactory::getSuffixPtr()->getName());
  PrsmWriter internal_writer(output_file_name + "_"
      + SemiAlignTypeFactory::getInternalPtr()->getName());
  PrsmWriter all_writer(output_file_name);

  //LOG_DEBUG("start reading");
  int n = 0;
  DeconvMsPtr ms_ptr = reader.getNextMs();
  LOG_DEBUG("init ms_ptr");
  double delta = 0;
  while (ms_ptr.get() != nullptr) {
    n++;
    SpectrumSetPtr spec_set_ptr 
        = getSpectrumSet(ms_ptr, delta, prsm_para_ptr->getSpParaPtr());
    if (spec_set_ptr.get() != nullptr) {
      PrsmPtrVec comp_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getCompletePtr(), 
                    prot_mod_forms, mng_ptr, comp_prsms);
      comp_writer.writeVector(comp_prsms);
      all_writer.writeVector(comp_prsms);
      PrsmPtrVec pref_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getPrefixPtr(), 
                    prot_mod_forms, mng_ptr, pref_prsms);
      pref_writer.writeVector(pref_prsms);
      all_writer.writeVector(pref_prsms);
      PrsmPtrVec suff_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getSuffixPtr(), 
                    raw_forms, mng_ptr, suff_prsms);
      suff_writer.writeVector(suff_prsms);
      all_writer.writeVector(suff_prsms);
      PrsmPtrVec internal_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getInternalPtr(), 
                    raw_forms, mng_ptr, internal_prsms);
      internal_writer.writeVector(internal_prsms);
      all_writer.writeVector(internal_prsms);
      if (log_file_name.length() != 0){
        if (logfile.is_open()) {
          logfile << (double) n / spectra_num * 0.053 << std::endl;
        }
      }
      std::cout << std::flush << "Zero ptm searching is processing " << n << " of " << spectra_num << " spectra.\r";
    }
    ms_ptr = reader.getNextMs();
    //LOG_DEBUG("spectrum " << n);
  }

  reader.close();
  logfile.close();

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
  all_writer.close();

  std::cout << std::endl << "Zero ptm searching finished." << std::endl;
}

} // end namespace

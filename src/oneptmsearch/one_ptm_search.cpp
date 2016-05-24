#include <chrono>

#include "base/logger.hpp"
#include "base/web_logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/db_block.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "oneptmsearch/one_ptm_slow_match.hpp"
#include "oneptmsearch/one_ptm_search.hpp"

namespace prot {

void onePtmSearchOneSpec(SpectrumSetPtr spec_set_ptr, 
                         SimplePrsmPtrVec &simple_prsm_ptr_vec,
                         //CompShiftLowMemPtr comp_shift_ptr,
                         FastaIndexReaderPtr reader_ptr,
                         PtmSearchMngPtr mng_ptr, 
                         AlignTypePtr type_ptr,
                         PrsmPtrVec &prsms) {

  ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ProtModPtrVec prot_mod_ptr_vec = mng_ptr->prsm_para_ptr_->getProtModPtrVec();
  ProteoformPtrVec proteoform_ptr_vec;
  SimplePrsmPtrVec prsm_vec;
  for (size_t i = 0; i < simple_prsm_ptr_vec.size(); i++) {
    std::string seq_name = simple_prsm_ptr_vec[i]->getSeqName();
    std::string seq_desc = simple_prsm_ptr_vec[i]->getSeqDesc();
    ProteoformPtr proteo_ptr = ProteoformFactory::readFastaToProteoformPtr(
        reader_ptr, seq_name, seq_desc, fix_mod_list);
    if (type_ptr == AlignType::COMPLETE || type_ptr == AlignType::PREFIX) {
      ProteoformPtrVec mod_form_ptr_vec = ProteoformFactory::geneProtModProteoform(
          proteo_ptr, prot_mod_ptr_vec);
      proteoform_ptr_vec.insert(proteoform_ptr_vec.end(), mod_form_ptr_vec.begin(), 
                                mod_form_ptr_vec.end());
      prsm_vec.insert(prsm_vec.end(), mod_form_ptr_vec.size(), simple_prsm_ptr_vec[i]);
    }
    else {
      proteoform_ptr_vec.push_back(proteo_ptr);
      prsm_vec.push_back(simple_prsm_ptr_vec[i]);
    }
  }

  for (size_t i = 0; i < proteoform_ptr_vec.size(); i++) { 
    //auto start = std::chrono::high_resolution_clock::now();
    /*
       PtmSlowMatch slow_match(proteoform_ptr_vec[i],
       spec_set_ptr, type_ptr,
       comp_shift_ptr, mng_ptr);
       */
    OnePtmSlowMatch slow_match(proteoform_ptr_vec[i], spec_set_ptr,
                               prsm_vec[i], type_ptr, mng_ptr);
    //auto step_1 = std::chrono::high_resolution_clock::now();
    //LOG_DEBUG("Init time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(step_1-start).count());
    PrsmPtr tmp = slow_match.compute(type_ptr, 1);
    //auto step_2 = std::chrono::high_resolution_clock::now();
    //LOG_DEBUG("Alignment time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(step_2-step_1).count());

    if (tmp != nullptr)
      prsms.push_back(tmp);
  }
  //LOG_DEBUG("prsm generation ended size " << prsms.size());
  std::sort(prsms.begin(), prsms.end(), Prsm::cmpMatchFragmentDecMatchPeakDec);
  if (prsms.size() > 0) {
    prsms.erase(prsms.begin() + 1, prsms.end());
  }
}

void OnePtmSearch::process(PtmSearchMngPtr mng_ptr){
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = FileUtil::basename(sp_file_name)+"."+mng_ptr->input_file_ext_;
  SimplePrsmReader comp_prsm_reader(input_file_name + "_" + AlignType::COMPLETE->getName());
  SimplePrsmReader pref_prsm_reader(input_file_name + "_" + AlignType::PREFIX->getName());
  SimplePrsmReader suff_prsm_reader(input_file_name + "_" + AlignType::SUFFIX->getName());
  SimplePrsmReader internal_prsm_reader(input_file_name + "_" + AlignType::INTERNAL->getName());
  SimplePrsmPtr comp_prsm_ptr = comp_prsm_reader.readOnePrsm();
  SimplePrsmPtr pref_prsm_ptr = pref_prsm_reader.readOnePrsm();
  SimplePrsmPtr suff_prsm_ptr = suff_prsm_reader.readOnePrsm();
  SimplePrsmPtr internal_prsm_ptr = internal_prsm_reader.readOnePrsm();

  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_;  
  PrsmXmlWriter comp_writer(output_file_name + "_" + AlignType::COMPLETE->getName());
  PrsmXmlWriter pref_writer(output_file_name + "_" + AlignType::PREFIX->getName());
  PrsmXmlWriter suff_writer(output_file_name + "_" + AlignType::SUFFIX->getName());
  PrsmXmlWriter internal_writer(output_file_name + "_" + AlignType::INTERNAL->getName());
  PrsmXmlWriter all_writer(output_file_name);

  //CompShiftLowMemPtr comp_shift_ptr = CompShiftLowMemPtr(new CompShiftLowMem());

  //init variables
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  FastaIndexReaderPtr reader_ptr(new FastaIndexReader(db_file_name));
  int spectrum_num = MsAlignUtil::getSpNum (sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          sp_para_ptr->getActivationPtr());
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;
  //LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpecId();
      // complete
      SimplePrsmPtrVec comp_selected_prsm_ptrs;
      while (comp_prsm_ptr != nullptr && comp_prsm_ptr->getSpectrumId() == spec_id) {
        comp_selected_prsm_ptrs.push_back(comp_prsm_ptr);
        comp_prsm_ptr = comp_prsm_reader.readOnePrsm();
      }
      if (comp_selected_prsm_ptrs.size() > 0) {
        //LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms;
        onePtmSearchOneSpec(spec_set_ptr, comp_selected_prsm_ptrs, 
                            reader_ptr, mng_ptr, AlignType::COMPLETE, prsms);
        comp_writer.writeVector(prsms);
        all_writer.writeVector(prsms);
      }

      // prefix
      SimplePrsmPtrVec pref_selected_prsm_ptrs;
      while (pref_prsm_ptr != nullptr && pref_prsm_ptr->getSpectrumId() == spec_id) {
        pref_selected_prsm_ptrs.push_back(pref_prsm_ptr);
        pref_prsm_ptr = pref_prsm_reader.readOnePrsm();
      }
      if (pref_selected_prsm_ptrs.size() > 0) {
        //LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms;
        onePtmSearchOneSpec(spec_set_ptr, pref_selected_prsm_ptrs, 
                            reader_ptr, mng_ptr, AlignType::PREFIX, prsms);
        pref_writer.writeVector(prsms);
        all_writer.writeVector(prsms);
      }

      // suffix
      SimplePrsmPtrVec suff_selected_prsm_ptrs;
      while (suff_prsm_ptr != nullptr && suff_prsm_ptr->getSpectrumId() == spec_id) {
        suff_selected_prsm_ptrs.push_back(suff_prsm_ptr);
        suff_prsm_ptr = suff_prsm_reader.readOnePrsm();
      }
      if (suff_selected_prsm_ptrs.size() > 0) {
        //LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms;
        onePtmSearchOneSpec(spec_set_ptr, suff_selected_prsm_ptrs, 
                            reader_ptr, mng_ptr, AlignType::SUFFIX, prsms);
        suff_writer.writeVector(prsms);
        all_writer.writeVector(prsms);
      }

      // internal
      SimplePrsmPtrVec internal_selected_prsm_ptrs;
      while (internal_prsm_ptr != nullptr && internal_prsm_ptr->getSpectrumId() == spec_id) {
        internal_selected_prsm_ptrs.push_back(internal_prsm_ptr);
        internal_prsm_ptr = internal_prsm_reader.readOnePrsm();
      }
      if (internal_selected_prsm_ptrs.size() > 0) {
        //LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms;
        onePtmSearchOneSpec(spec_set_ptr, internal_selected_prsm_ptrs,  
                            reader_ptr, mng_ptr, AlignType::INTERNAL, prsms);
        internal_writer.writeVector(prsms);
        all_writer.writeVector(prsms);
      }
    }
    std::cout << std::flush <<  "One PTM search is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
    WebLog::percentLog(cnt, spectrum_num, WebLog::OnePtmSearchTime());
  }
  sp_reader.close();
  comp_prsm_reader.close();
  pref_prsm_reader.close();
  suff_prsm_reader.close();
  internal_prsm_reader.close();
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
  all_writer.close();
  std::cout << std::endl;
}

}

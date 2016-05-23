#include "base/prot_mod.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "base/web_logger.hpp"
#include "base/threadpool.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_search_processor.hpp"
#include "ptmsearch/ptm_search_slow_filter.hpp"

namespace prot {

template <int N>
PrsmXmlWriterSet<N>::PrsmXmlWriterSet(const std::string & output_file_name){
  all_writer_ptr_ = PrsmXmlWriterPtr(new PrsmXmlWriter (output_file_name));
  for (int s = 2; s <= N; s++) {
    std::string file_name = output_file_name+"_"+ StringUtil::StringUtil::convertToString(s)
        +"_"+ AlignType::COMPLETE->getName();
    PrsmXmlWriterPtr complete_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
    complete_writer_ptrs_.push_back(complete_writer_ptr);
    file_name = output_file_name+"_"+ StringUtil::convertToString(s)
        +"_"+ AlignType::PREFIX->getName();
    PrsmXmlWriterPtr prefix_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
    prefix_writer_ptrs_.push_back(prefix_writer_ptr);
    file_name = output_file_name+"_"+ StringUtil::convertToString(s)
        +"_"+ AlignType::SUFFIX->getName();
    PrsmXmlWriterPtr suffix_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
    suffix_writer_ptrs_.push_back(suffix_writer_ptr);
    file_name = output_file_name+"_"+ StringUtil::convertToString(s)
        +"_"+ AlignType::INTERNAL->getName();
    PrsmXmlWriterPtr internal_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
    internal_writer_ptrs_.push_back(internal_writer_ptr);
  }
}

template<int N>
void PrsmXmlWriterSet<N>::close() {
  all_writer_ptr_->close();
  for (int s = 2; s <= N; s++) {
    complete_writer_ptrs_[s-2]->close();
    prefix_writer_ptrs_[s-2]->close();
    suffix_writer_ptrs_[s-2]->close();
    internal_writer_ptrs_[s-2]->close();
  }
}

void seleTopPrsms(const PrsmPtrVec &all_prsm_ptrs, 
                  PrsmPtrVec &sele_prsm_ptrs, int n_report) {
  int match_size = all_prsm_ptrs.size();
  if(all_prsm_ptrs.size()!=0){
    for(int r=0;r< n_report;r++){
      if(r >= match_size){
        break;
      }
      if(all_prsm_ptrs[r]->getMatchFragNum() > 0){
        sele_prsm_ptrs.push_back(all_prsm_ptrs[r]);
      }
    }
  }
  std::sort(sele_prsm_ptrs.begin(), sele_prsm_ptrs.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);
}

template <int N>
std::function<void()> geneTask(SpectrumSetPtr spectrum_set_ptr, 
                               SimplePrsmPtrVec ori_simple_prsm_ptrs,
                               PtmSearchMngPtr mng_ptr,
                               CompShiftLowMem comp_shift,
                               std::shared_ptr<ThreadPool<PrsmXmlWriterSet<N>>> pool_ptr) {

  return [spectrum_set_ptr, ori_simple_prsm_ptrs, mng_ptr, comp_shift, pool_ptr]() {
    SimplePrsmPtrVec simple_prsm_ptrs = 
        SimplePrsmUtil::getUniqueMatches(ori_simple_prsm_ptrs);
    PtmSearchSlowFilterPtr slow_filter_ptr = 
        std::make_shared<PtmSearchSlowFilter>(spectrum_set_ptr, simple_prsm_ptrs,
                                              comp_shift, mng_ptr);
    boost::thread::id thread_id = boost::this_thread::get_id();
    std::shared_ptr<PrsmXmlWriterSet<N>> writer_ptr = pool_ptr->getWriter(thread_id); 
    for (int s = 2; s <= N; s++) {
      PrsmPtrVec complete_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, AlignType::COMPLETE);
      std::sort(complete_prsm_ptrs.begin(), complete_prsm_ptrs.end(), 
                Prsm::cmpMatchFragDecStartPosInc);
      PrsmPtrVec sele_complete_prsm_ptrs;
      seleTopPrsms(complete_prsm_ptrs, sele_complete_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->complete_writer_ptrs_[s-2]->writeVector(sele_complete_prsm_ptrs);
      writer_ptr->all_writer_ptr_->writeVector(sele_complete_prsm_ptrs);

      PrsmPtrVec prefix_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, AlignType::PREFIX);
      std::sort(prefix_prsm_ptrs.begin(), prefix_prsm_ptrs.end(), 
                Prsm::cmpMatchFragDecStartPosInc);
      PrsmPtrVec sele_prefix_prsm_ptrs;
      seleTopPrsms(prefix_prsm_ptrs, sele_prefix_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->prefix_writer_ptrs_[s-2]->writeVector(sele_prefix_prsm_ptrs);
      writer_ptr->all_writer_ptr_->writeVector(sele_prefix_prsm_ptrs);

      PrsmPtrVec suffix_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, AlignType::SUFFIX);
      std::sort(suffix_prsm_ptrs.begin(), suffix_prsm_ptrs.end(), Prsm::cmpMatchFragmentDec);
      PrsmPtrVec sele_suffix_prsm_ptrs;
      seleTopPrsms(suffix_prsm_ptrs, sele_suffix_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->suffix_writer_ptrs_[s-2]->writeVector(sele_suffix_prsm_ptrs);
      writer_ptr->all_writer_ptr_->writeVector(sele_suffix_prsm_ptrs);

      PrsmPtrVec internal_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, AlignType::INTERNAL);
      std::sort(internal_prsm_ptrs.begin(), internal_prsm_ptrs.end(), Prsm::cmpMatchFragmentDec);
      PrsmPtrVec sele_internal_prsm_ptrs;
      seleTopPrsms(internal_prsm_ptrs, sele_internal_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->internal_writer_ptrs_[s-2]->writeVector(sele_internal_prsm_ptrs);
      writer_ptr->all_writer_ptr_->writeVector(sele_internal_prsm_ptrs);
    }
  };
}

PtmSearchProcessor::PtmSearchProcessor(PtmSearchMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  comp_shift_ = CompShiftLowMem();
}


// process ptm search
void PtmSearchProcessor::process(){
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = FileUtil::basename(sp_file_name)+"."+mng_ptr_->input_file_ext_;
  SimplePrsmReader simple_prsm_reader(input_file_name);
  SimplePrsmPtr prsm_ptr = simple_prsm_reader.readOnePrsm();

  //init variables
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  FastaIndexReaderPtr reader_ptr(new FastaIndexReader(db_file_name));
  int spectrum_num = MsAlignUtil::getSpNum (sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  std::string output_file_name = FileUtil::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          sp_para_ptr->getActivationPtr());

  const int n_unknonw_shift = 2;

  std::shared_ptr<ThreadPool<PrsmXmlWriterSet<n_unknonw_shift>>> pool_ptr 
      = std::make_shared<ThreadPool<PrsmXmlWriterSet<n_unknonw_shift>>>(mng_ptr_->thread_num_ , output_file_name);

  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;
  //LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpecId();
      SimplePrsmPtrVec selected_prsm_ptrs;
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        selected_prsm_ptrs.push_back(prsm_ptr);
        prsm_ptr = simple_prsm_reader.readOnePrsm();
      }
      if (selected_prsm_ptrs.size() > 0) {
        while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
          boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(spec_set_ptr, selected_prsm_ptrs,
                                   mng_ptr_, comp_shift_, pool_ptr));
      }
    }
    std::cout << std::flush <<  "PTM search is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";

    WebLog::percentLog(cnt, spectrum_num, WebLog::TwoPtmSearchTime());
  }
  pool_ptr->ShutDown();
  LOG_DEBUG("Search completed");
  sp_reader.close();
  simple_prsm_reader.close();

  std::cout << std::endl;
  PrsmXmlWriterPtr all_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name);
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    PrsmReaderPtr all_reader_ptr = std::make_shared<PrsmReader>(output_file_name + "_" + std::to_string(i)); 
    PrsmPtr p = all_reader_ptr->readOnePrsm(reader_ptr, fix_mod_ptr_vec);
    while(p != nullptr) {
      all_writer_ptr->write(p); 
      p = all_reader_ptr->readOnePrsm(reader_ptr, fix_mod_ptr_vec); 
    }
    all_reader_ptr->close();
  }
  all_writer_ptr->close();
}

} /* namespace prot */

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/db_block.hpp"
#include "base/proteoform_factory.hpp"
#include "base/file_util.hpp"
#include "base/web_logger.hpp"
#include "base/threadpool.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "diagfilter/mass_diag_filter.hpp"
#include "diagfilter/diag_filter_processor.hpp"

namespace prot {

typedef std::shared_ptr<ThreadPool<SimplePrsmXmlWriter>> SimplePrsmThreadPoolPtr;

std::function<void()> geneTask(MassDiagFilterPtr filter_ptr,
                               const PrmMsPtrVec & ms_ptr_vec,
                               SimplePrsmThreadPoolPtr  pool_ptr) {
  return [filter_ptr, ms_ptr_vec, pool_ptr]() {
    SimplePrsmPtrVec match_ptrs = filter_ptr->getBestMatch(ms_ptr_vec);
    boost::thread::id thread_id = boost::this_thread::get_id();
    std::shared_ptr<SimplePrsmXmlWriter> writer_ptr = pool_ptr->getWriter(thread_id);
    writer_ptr->write(match_ptrs);
  };
}

DiagFilterProcessor::DiagFilterProcessor(DiagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
}

void DiagFilterProcessor::process(){
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  for(size_t i=0; i< db_block_ptr_vec.size(); i++){
    std::cout << "Diagonal filtering block " << (i+1) << " out of " 
        << db_block_ptr_vec.size() << " started." << std::endl; 
    processBlock(db_block_ptr_vec[i], db_block_ptr_vec.size());
    std::cout << "Diagonal filtering block " << (i +1) 
        << " finished. " << std::endl;
  }

  std::cout << "Diagonal filtering: combining blocks started." << std::endl; 

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int block_num = db_block_ptr_vec.size();

  SimplePrsmStrCombinePtr combine_ptr(new SimplePrsmStrCombine(sp_file_name, mng_ptr_->output_file_ext_,
                                                               block_num, mng_ptr_->output_file_ext_, 
                                                               mng_ptr_->filter_result_num_));
  combine_ptr->process();

  std::cout << "Diagonal filtering: combining blocks finished." << std::endl; 
}

void DiagFilterProcessor::processBlock(DbBlockPtr block_ptr, int total_block_num) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName() 
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms 
      = ProteoformFactory::readFastaToProteoformPtrVec(db_block_file_name, 
                                                       prsm_para_ptr->getFixModPtrVec());
  MassDiagFilterPtr filter_ptr(new MassDiagFilter(raw_forms, mng_ptr_));

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       sp_para_ptr->getActivationPtr());

  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+ std::to_string(block_ptr->getBlockIdx());

  SimplePrsmThreadPoolPtr pool_ptr
      = std::make_shared<ThreadPool<SimplePrsmXmlWriter>>(mng_ptr_->thread_num_, output_file_name);

  SpectrumSetPtr spec_set_ptr;
  int spectrum_num = MsAlignUtil::getSpNum (prsm_para_ptr->getSpectrumFileName());
  int cnt = 0;
  while((spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)) != nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      PrmMsPtrVec ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
      while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
      }
      pool_ptr->Enqueue(geneTask(filter_ptr, ms_ptr_vec, pool_ptr));
    }
    WebLog::percentLog(cnt, spectrum_num, block_ptr->getBlockIdx(), total_block_num, 
                       WebLog::DiagFilterTime());
    std::cout << std::flush << "Diagonal filtering is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;
  reader.close();
  std::vector<std::string> input_exts;
  std::string cur_output_ext = mng_ptr_->output_file_ext_+"_"+ std::to_string(block_ptr->getBlockIdx());
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    std::string fname = cur_output_ext + "_" + std::to_string(i); 
    input_exts.push_back(fname);
  }

  SimplePrsmStrCombinePtr combine_ptr
      = std::make_shared<SimplePrsmStrCombine>(mng_ptr_->prsm_para_ptr_->getSpectrumFileName(),
                                               input_exts, cur_output_ext,  INT_MAX);
  combine_ptr->process();
  combine_ptr = nullptr; 
}

} /* namespace prot */

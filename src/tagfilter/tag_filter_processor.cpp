
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/db_block.hpp"
#include "base/file_util.hpp"
#include "base/proteoform_factory.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "spec/msalign_util.hpp"
#include "tag_filter_processor.hpp"
#include "tag_filter.hpp"

namespace prot{

TagFilterProcessor::TagFilterProcessor(TagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
}

void TagFilterProcessor::process(){
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  for(size_t i=0; i< db_block_ptr_vec.size(); i++){
    processBlock(db_block_ptr_vec[i], db_block_ptr_vec.size());
  }
}

void TagFilterProcessor::processBlock(DbBlockPtr block_ptr, 
                                      int total_block_num) {

  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName() 
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms 
      = ProteoformFactory::readFastaToProteoformPtrVec(db_block_file_name, 
                                                       prsm_para_ptr->getFixModPtrVec());

  TagFilterPtr filter_ptr = std::make_shared<TagFilter>(raw_forms, mng_ptr_);

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       sp_para_ptr->getActivationPtr());

  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+ std::to_string(block_ptr->getBlockIdx());

  std::ofstream writer;
  writer.open(output_file_name);

  std::vector<double> prec_errors;
  prec_errors.push_back(0);
  for (int i = 1; i <= mng_ptr_->prec_error_; i++) {
    prec_errors.push_back(- i * MassConstant::getIsotopeMass());
    prec_errors.push_back(i * MassConstant::getIsotopeMass());
  }

  SpectrumSetPtr spec_set_ptr;

  while((spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)) != nullptr){
    DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
    double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
    std::cout << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId() << std::endl;
    if (spec_set_ptr->isValid()) {
      for (size_t i = 0; i < prec_errors.size(); i++) {
        SpectrumSetPtr adjusted_spec_set_ptr(
            new SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass + prec_errors[i]));
        PrmMsPtrVec ms_six_vec = adjusted_spec_set_ptr->getMsSixPtrVec();
        //std::vector<std::string> seq_vec = 
        filter_ptr->getBestMatch(ms_six_vec);
        /*for (size_t j = 0; j < seq_vec.size(); j++) {*/
        //writer << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId() << "\t" 
        //<< seq_vec[j] << "\t"
        //<< prec_errors[i] << std::endl;
        /*}*/
      }
    }
  }
  reader.close();
  writer.close();
}

}

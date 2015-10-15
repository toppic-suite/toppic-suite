#include <map>

#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_species.hpp"
#include "prsmview/xml_generator.hpp"

namespace prot {
XmlGenerator::XmlGenerator(PrsmParaPtr prsm_para_ptr, 
                           const std::string &exec_dir, 
                           const std::string &input_file_ext) {
  input_file_ext_ = input_file_ext;
  mng_ptr_ = PrsmViewMngPtr(new PrsmViewMng(prsm_para_ptr, exec_dir));
  anno_view_ptr_ = AnnoViewPtr(new AnnoView());
}

void XmlGenerator::outputPrsms(const PrsmPtrVec &prsm_ptrs){
  for(unsigned int i=0;i<prsm_ptrs.size();i++){
    std::string file_name = mng_ptr_->xml_path_+ FILE_SEPARATOR + 
        "prsms" + FILE_SEPARATOR + "prsm"+convertToString(prsm_ptrs[i]->getId())+".xml";
    XmlWriter writer(file_name,"");
    writer.write(genePrsmView(writer.getDoc(),prsm_ptrs[i], mng_ptr_));
    writer.close();

    LOG_DEBUG("writer prsm completed");

    std::vector<std::string> file_info;
    file_info.push_back(file_name);
    file_info.push_back(mng_ptr_->executive_dir_ + FILE_SEPARATOR + "toppic_resources" + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "prsm.xsl");
    file_info.push_back(mng_ptr_->html_path_+ FILE_SEPARATOR + "prsms" + FILE_SEPARATOR 
                        + "prsm"+ convertToString(prsm_ptrs[i]->getId())+".html");
    anno_view_ptr_->file_list_.push_back(file_info);

  }

}
void XmlGenerator::outputAllPrsms(const PrsmPtrVec &prsm_ptrs){
  std::string file_name = mng_ptr_->xml_path_+ FILE_SEPARATOR + "prsms.xml";
  XmlWriter writer(file_name,"prsm_list");
  for(unsigned int i=0;i<prsm_ptrs.size();i++){
    writer.write(genePrsmView(writer.getDoc(),prsm_ptrs[i], 
                              mng_ptr_));
    writer.close();
  }
}

void XmlGenerator::outputProteoforms(const PrsmPtrVec &prsm_ptrs){

  std::vector<int> species_ids = getSpeciesIds(prsm_ptrs);
  for(unsigned int i=0;i<species_ids.size();i++){
    PrsmPtrVec select_prsm_ptrs = selectSpeciesPrsms(prsm_ptrs,species_ids[i]);
    if(select_prsm_ptrs.size()>0){
      std::string file_name = mng_ptr_->xml_path_+ FILE_SEPARATOR + "proteoforms" 
          + FILE_SEPARATOR + "proteoform"+convertToString(species_ids[i])+".xml";
      XmlWriter writer(file_name,"");
      std::sort(select_prsm_ptrs.begin(),select_prsm_ptrs.end(),prsmEValueUp);
      writer.write(proteoformToXml(writer.getDoc(),select_prsm_ptrs, mng_ptr_));
      writer.close();

      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->executive_dir_ + FILE_SEPARATOR + "toppic_resources" + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "proteoform.xsl");
      file_info.push_back(mng_ptr_->html_path_+ FILE_SEPARATOR+ "proteoforms" + FILE_SEPARATOR 
                          + "proteoform"+convertToString(species_ids[i])+".html");
      anno_view_ptr_->file_list_.push_back(file_info);
    }
  }
}

void XmlGenerator::outputProteins(const PrsmPtrVec &prsm_ptrs){
  //LOG_DEBUG("prsm number " << prsm_ptrs.size());
  ProteoformReader reader(mng_ptr_->prsm_para_ptr_->getSearchDbFileName());
  ResiduePtrVec residue_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModResiduePtrVec();
  ProteoformPtr proteo_ptr = reader.getNextProteoformPtr(residue_ptr_vec);

  while (proteo_ptr != nullptr) { 
    std::vector<int> species = getSpeciesIds(prsm_ptrs,proteo_ptr->getDbResSeqPtr()->getId());
    //LOG_DEBUG("species size " << species.size());
    if(species.size()>0){
      std::string file_name = mng_ptr_->xml_path_ + FILE_SEPARATOR +"proteins" 
          +FILE_SEPARATOR+ "protein"+convertToString(proteo_ptr->getDbResSeqPtr()->getId())+".xml";
      XmlWriter writer(file_name,"");
      writer.write(proteinToXml(writer.getDoc(),prsm_ptrs,proteo_ptr,species, mng_ptr_));
      writer.close();
      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->executive_dir_ + FILE_SEPARATOR + "toppic_resources" + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "protein.xsl");
      file_info.push_back(mng_ptr_->html_path_+ FILE_SEPARATOR + "proteins" + FILE_SEPARATOR 
                          + "protein"+convertToString(proteo_ptr->getDbResSeqPtr()->getId())+".html");
      anno_view_ptr_->file_list_.push_back(file_info);
    }
    proteo_ptr = reader.getNextProteoformPtr(residue_ptr_vec);
  }
}

void XmlGenerator::outputAllProteins(const PrsmPtrVec &prsm_ptrs){

  std::string file_name = mng_ptr_->xml_path_+ FILE_SEPARATOR +"proteins.xml";
  XmlWriter writer(file_name,"protein_list");
  writer.write(allProteinToXml(writer.getDoc(),prsm_ptrs, mng_ptr_));
  writer.close();
  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_ptr_->executive_dir_ + FILE_SEPARATOR + "toppic_resources" + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR +"proteins.xsl");
  file_info.push_back(mng_ptr_->html_path_+ FILE_SEPARATOR + "proteins.html");
  anno_view_ptr_->file_list_.push_back(file_info);
}

void XmlGenerator::outputFileList(){
  std::string file_name = mng_ptr_->xml_path_+ FILE_SEPARATOR + "files.xml";
  XmlWriter writer(file_name,"");
  writer.write(anno_view_ptr_->geneFileList(writer.getDoc()));
  writer.close();
}
void XmlGenerator::process(){
  LOG_DEBUG("process start");
  prot::createFolder(mng_ptr_->xml_path_ + FILE_SEPARATOR + "proteoforms");
  prot::createFolder(mng_ptr_->xml_path_ + FILE_SEPARATOR + "prsms");
  prot::createFolder(mng_ptr_->xml_path_ + FILE_SEPARATOR + "proteins");
  LOG_DEBUG("fold created");

  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string spectrum_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = basename(spectrum_file_name) + "." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  ResiduePtrVec residue_ptr_vec = prsm_para_ptr->getFixModResiduePtrVec();

  PrsmPtrVec prsm_ptrs = readAllPrsms(input_file_name, db_file_name,
                                      residue_ptr_vec);
  LOG_DEBUG("prsm loaded");

  addSpectrumPtrsToPrsms(prsm_ptrs, prsm_para_ptr);
  LOG_DEBUG("spectrum added");

  std::sort(prsm_ptrs.begin(), prsm_ptrs.end(), prsmEValueUp); 

  outputPrsms(prsm_ptrs);
  outputAllPrsms(prsm_ptrs);
  outputProteoforms(prsm_ptrs);
  outputProteins(prsm_ptrs);
  outputAllProteins(prsm_ptrs);
  outputFileList();
}

} /* namespace prot */

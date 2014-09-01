#include <map>

#include "base/file_util.hpp"
#include "prsm/prsm_species.hpp"
#include "prsmview/xml_generator.hpp"

namespace prot {
XmlGenerator::XmlGenerator(PrsmParaPtr prsm_para_ptr, 
                           std::string exec_dir, 
                           std::string input_file) {
  input_file_ = input_file;
  mng_ = PrsmViewMngPtr(new PrsmViewMng(prsm_para_ptr, exec_dir));
  anno_view_ = AnnoViewPtr(new AnnoView());
}

void XmlGenerator::outputPrsms(PrsmPtrVec prsms){
  for(unsigned int i=0;i<prsms.size();i++){
    std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + 
        "prsms" + FILE_SEPARATOR + "prsm"+convertToString(prsms[i]->getId())+".xml";
    XmlWriter writer(file_name,"");
    writer.write(genePrsmView(writer.getDoc(),prsms[i], mng_));
    writer.close();

    std::vector<std::string> file_info;
    file_info.push_back(file_name);
    file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "prsm.xsl");
    file_info.push_back(mng_->html_path_+ FILE_SEPARATOR + "prsms" + FILE_SEPARATOR 
                        + "prsm"+ convertToString(prsms[i]->getId())+".html");
    anno_view_->file_list_.push_back(file_info);

  }

}
void XmlGenerator::outputAllPrsms(PrsmPtrVec prsms){
  std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + "prsms.xml";
  XmlWriter writer(file_name,"prsm_list");
  for(unsigned int i=0;i<prsms.size();i++){
    writer.write(genePrsmView(writer.getDoc(),prsms[i], 
                              mng_));
    writer.close();
  }
}

void XmlGenerator::outputProteins(PrsmPtrVec prsms){

  for(unsigned int i=0;i<raw_forms_.size();i++){
    std::vector<int> species = getSpeciesIds(prsms,raw_forms_[i]->getDbResSeqPtr()->getId());
    if(species.size()>0){
      std::string file_name = mng_->xml_path_ + FILE_SEPARATOR +"proteins" 
          +FILE_SEPARATOR+ "protein"+convertToString(raw_forms_[i]->getDbResSeqPtr()->getId())+".xml";
      XmlWriter writer(file_name,"");
      writer.write(proteinToXml(writer.getDoc(),prsms,raw_forms_[i],species, 
                              mng_));
      writer.close();
      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "protein.xsl");
      file_info.push_back(mng_->html_path_+ FILE_SEPARATOR + "proteins" + FILE_SEPARATOR 
                          + "protein"+convertToString(raw_forms_[i]->getDbResSeqPtr()->getId())+".html");
      anno_view_->file_list_.push_back(file_info);
    }
  }
}
void XmlGenerator::outputAllProteins(PrsmPtrVec prsms){

  std::string file_name = mng_->xml_path_+ FILE_SEPARATOR +"proteins.xml";
  XmlWriter writer(file_name,"protein_list");
  writer.write(allProteinToXml(writer.getDoc(),prsms,raw_forms_, 
                               mng_));
  writer.close();
  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR +"proteins.xsl");
  file_info.push_back(mng_->html_path_+ FILE_SEPARATOR + "proteins.html");
  anno_view_->file_list_.push_back(file_info);
}

void XmlGenerator::outputSpecies(PrsmPtrVec prsms){

  std::vector<int> species = getSpeciesIds(prsms);
  for(unsigned int i=0;i<species.size();i++){
    PrsmPtrVec select_prsms = selectSpeciesPrsms(prsms,species[i]);
    if(select_prsms.size()>0){
      std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + "species" 
          + FILE_SEPARATOR + "species"+convertToString(species[i])+".xml";
      XmlWriter writer(file_name,"");
      std::sort(select_prsms.begin(),select_prsms.end(),prsmEValueUp);
      writer.write(speciesToXml(writer.getDoc(),select_prsms, 
                                mng_));
      writer.close();

      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "species.xsl");
      file_info.push_back(mng_->html_path_+ FILE_SEPARATOR+ "species" + FILE_SEPARATOR 
                          + "species"+convertToString(species[i])+".html");
      anno_view_->file_list_.push_back(file_info);
    }
  }
}
void XmlGenerator::outputFileList(){
  std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + "files.xml";
  XmlWriter writer(file_name,"");
  writer.write(anno_view_->geneFileList(writer.getDoc()));
  writer.close();
}
void XmlGenerator::process(){
  prot::createFolder(mng_->xml_path_ + FILE_SEPARATOR + "species");
  prot::createFolder(mng_->xml_path_ + FILE_SEPARATOR + "prsms");
  prot::createFolder(mng_->xml_path_ + FILE_SEPARATOR + "proteins");

  PrsmParaPtr prsm_para_ptr = mng_->prsm_para_ptr_;
  std::string spectrum_file_name = prsm_para_ptr->getSpectrumFileName();
  raw_forms_ = readFastaToProteoform(prsm_para_ptr->getSearchDbFileName(),
                                     prsm_para_ptr->getFixModResiduePtrVec());
  std::string input_name = basename(spectrum_file_name)+"."+input_file_;
  PrsmPtrVec prsms = readPrsm(basename(spectrum_file_name)+"."+input_file_,raw_forms_);

  addSpectrumPtrsToPrsms(prsms, prsm_para_ptr);
  setSpeciesId(prsms,prsm_para_ptr->getSpParaPtr()->getPeakTolerancePtr()->getPpo());
  outputPrsms(prsms);
  outputAllPrsms(prsms);
  outputSpecies(prsms);
  outputProteins(prsms);
  outputAllProteins(prsms);
  outputFileList();

}

} /* namespace prot */

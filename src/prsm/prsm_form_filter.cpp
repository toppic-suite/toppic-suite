#include "base/file_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_form_filter.hpp"


namespace prot {

PrsmFormFilter::PrsmFormFilter(const std::string &db_file_name,
                               const std::string &spec_file_name,
                               const std::string &input_file_ext,
                               const std::string &output_file_ext,
                               const std::string &output_file_ext_2):
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext),
    output_file_ext_2_(output_file_ext_2) {
    }

void PrsmFormFilter::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;

  //PrsmStrPtrVec prsms = PrsmReader::readAllPrsmStrs(input_file_name);
  ModPtrVec fix_mod_list;
  PrsmPtrVec prsms = PrsmReader::readAllPrsms(input_file_name, db_file_name_, fix_mod_list );

  sort(prsms.begin(),prsms.end(), Prsm::cmpEValueInc);

  PrsmPtrVec selected_prsms;

  PrsmPtrVec selected_forms;

  for(size_t i=0; i<prsms.size(); i++){
    //std::cout << "prsm " << i << std::endl;
    bool found = false; 
    for (size_t j = 0; j < selected_forms.size(); j++) {
      if (selected_forms[j]->getProteoformPtr()->getSpeciesId() == 
          prsms[i]->getProteoformPtr()->getSpeciesId()) {
        found = true;
        break;
      }
    }
    if (found) {
      selected_prsms.push_back(prsms[i]);
    }
    else {
      bool keep = true;
      std::string form = prsms[i]->getProteoformPtr()->getProteinMatchSeq();
      for (size_t j = 0; j < selected_forms.size(); j++) {
        if (selected_forms[j]->getPrecFeatureId() == prsms[i]->getPrecFeatureId()) {
          //std::cout << "scan " << prsms[i]->getSpectrumScan() << " removed by scan " << selected_forms[j]->getSpectrumScan() << std::endl;
          keep = false;
          break;
        }
        if (selected_forms[j]->getProteoformPtr()->getProteinMatchSeq() == form) {
          //std::cout << "scan " << prsms[i]->getSpectrumScan() << " removed by scan " << selected_forms[j]->getSpectrumScan() << std::endl;
          keep = false;
          break;
        }
      }
      if (keep) {
        selected_forms.push_back(prsms[i]);
        selected_prsms.push_back(prsms[i]);
      }
    }
  }

  //output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  sort(selected_prsms.begin(), selected_prsms.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);
  writer.writeVector(selected_prsms);
  writer.close();

  output_file_name = base_name + "." + output_file_ext_2_;
  PrsmXmlWriter writer2(output_file_name);
  sort(selected_forms.begin(), selected_forms.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);
  writer2.writeVector(selected_forms);
  writer2.close();
}

} /* namespace prot */

/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <base/logger.hpp>

#include "base/base_data.hpp"
#include "string_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

BaseData::BaseData(std::string config_file_name) {

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    LOG_DEBUG("config_file_name: " << config_file_name);
    XmlDOMDocument* doc = new XmlDOMDocument(parser, config_file_name.c_str());
    LOG_DEBUG("doc " << doc);
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      LOG_DEBUG("root " << root);
      std::string acid_file_name = getChildValue(root, "acid_list_file_name",
                                                 0);
      LOG_DEBUG("acid file name: " << acid_file_name);
      AcidFactory::initFactory(acid_file_name);
      LOG_DEBUG("acid initialized ");

      std::string ptm_file_name = getChildValue(root, "ptm_list_file_name", 0);
      LOG_DEBUG("ptm file name: " << ptm_file_name);
      PtmFactory::initFactory(ptm_file_name);
      LOG_DEBUG("ptm initialized");

      std::string residue_file_name = getChildValue(root,
                                                    "residue_list_file_name",
                                                    0);
      LOG_DEBUG("residue file name: " << residue_file_name);
      ResidueFactory::initFactory(residue_file_name);
      LOG_DEBUG("residue initialized");

      std::string trunc_file_name = getChildValue(root, "trunc_list_file_name",
                                                  0);
      LOG_DEBUG("trunc file name: " << trunc_file_name);
      TruncFactory::initFactory(trunc_file_name);
      LOG_DEBUG("trunc initialized ");

      std::string prot_mod_file_name = getChildValue(root,
                                                     "prot_mod_list_file_name",
                                                     0);
      LOG_DEBUG("prot mod file name: " << prot_mod_file_name);
      ProtModFactory::initFactory(prot_mod_file_name);
      LOG_DEBUG("prot mod initialized ");

      std::string ion_type_file_name = getChildValue(root,
                                                     "ion_type_list_file_name",
                                                     0);
      LOG_DEBUG("ion type file name: " << ion_type_file_name);
      IonTypeFactory::initFactory(ion_type_file_name);
      LOG_DEBUG("ion type initialized ");

      std::string neutral_loss_file_name = getChildValue(
          root, "neutral_loss_list_file_name", 0);
      LOG_DEBUG("neutral loss file name: " << neutral_loss_file_name);
      NeutralLossFactory::initFactory(neutral_loss_file_name);
      LOG_DEBUG("neutral loss initialized ");

      std::string activation_file_name = getChildValue(
          root, "activation_list_file_name", 0);
      LOG_DEBUG("activation file name: " << activation_file_name);
      ActivationFactory::initFactory(activation_file_name);
      LOG_DEBUG("activation initialized ");

      std::string sp_type_file_name = getChildValue(
          root, "support_peak_type_file_name", 0);
      LOG_DEBUG("support_peak_type_file_name: " << sp_type_file_name);
      SPTypeFactory::initFactory(sp_type_file_name);
      LOG_DEBUG("support peak type initialized ");

      std::string fix_mod_residue_file_name = getChildValue(
          root, "fix_mod_residue_file_name", 0);
      LOG_DEBUG("fix mod residue file name: " << fix_mod_residue_file_name);
      fix_mod_residue_list_ = ResidueFactory::getResiduePtrVecInstance(
          fix_mod_residue_file_name);
      LOG_DEBUG("fix mod residue initialized ");

      LOG_DEBUG("allow prot mods initialization ");
      xercesc::DOMElement* parent = getChildElement(root, "allow_prot_mod_list",
                                                    0);
      int prot_mod_num = getChildCount(parent, "prot_mod");
      for (int i = 0; i < prot_mod_num; i++) {
        std::string mod_name = getChildValue(parent, "prot_mod", i);
        ProtModPtr ptr = ProtModFactory::getBaseProtModPtrByName(mod_name);
        allow_prot_mod_list_.push_back(ptr);
      }
      LOG_DEBUG("allow prot mods initialized ");

      std::string activation_type = getChildValue(root, "activation_type", 0);
      LOG_DEBUG("acitivation type: " << activation_type);
      activation_ptr_ = ActivationFactory::getBaseActivationPtrByName(
          activation_type);
    }
    delete doc;
  }
}

BaseData::BaseData(std::map<std::string, std::string> arguments) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    LOG_DEBUG("config_file_name: " << arguments["configuration"]);
    XmlDOMDocument* doc = new XmlDOMDocument(parser, arguments["configuration"].c_str());
    LOG_DEBUG("doc " << doc);
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      LOG_DEBUG("root " << root);
      std::string acid_file_name = getChildValue(root, "acid_list_file_name",
                                                 0);
      LOG_DEBUG("acid file name: " << acid_file_name);
      AcidFactory::initFactory(acid_file_name);
      LOG_DEBUG("acid initialized ");

      std::string ptm_file_name = getChildValue(root, "ptm_list_file_name", 0);
      LOG_DEBUG("ptm file name: " << ptm_file_name);
      PtmFactory::initFactory(ptm_file_name);
      LOG_DEBUG("ptm initialized");

      std::string residue_file_name = getChildValue(root,
                                                    "residue_list_file_name",
                                                    0);
      LOG_DEBUG("residue file name: " << residue_file_name);
      ResidueFactory::initFactory(residue_file_name);
      LOG_DEBUG("residue initialized");

      std::string trunc_file_name = getChildValue(root, "trunc_list_file_name",
                                                  0);
      LOG_DEBUG("trunc file name: " << trunc_file_name);
      TruncFactory::initFactory(trunc_file_name);
      LOG_DEBUG("trunc initialized ");

      std::string prot_mod_file_name = getChildValue(root,
                                                     "prot_mod_list_file_name",
                                                     0);
      LOG_DEBUG("prot mod file name: " << prot_mod_file_name);
      ProtModFactory::initFactory(prot_mod_file_name);
      LOG_DEBUG("prot mod initialized ");

      std::string ion_type_file_name = getChildValue(root,
                                                     "ion_type_list_file_name",
                                                     0);
      LOG_DEBUG("ion type file name: " << ion_type_file_name);
      IonTypeFactory::initFactory(ion_type_file_name);
      LOG_DEBUG("ion type initialized ");

      std::string neutral_loss_file_name = getChildValue(
          root, "neutral_loss_list_file_name", 0);
      LOG_DEBUG("neutral loss file name: " << neutral_loss_file_name);
      NeutralLossFactory::initFactory(neutral_loss_file_name);
      LOG_DEBUG("neutral loss initialized ");

      std::string activation_file_name = getChildValue(
          root, "activation_list_file_name", 0);
      LOG_DEBUG("activation file name: " << activation_file_name);
      ActivationFactory::initFactory(activation_file_name);
      LOG_DEBUG("activation initialized ");

      std::string sp_type_file_name = getChildValue(
          root, "support_peak_type_file_name", 0);
      LOG_DEBUG("support_peak_type_file_name: " << sp_type_file_name);
      SPTypeFactory::initFactory(sp_type_file_name);
      LOG_DEBUG("support peak type initialized ");

      std::string fix_mod_residue_file_name = getChildValue(
          root, "fix_mod_residue_file_name", 0);
      LOG_DEBUG("fix mod residue file name: " << fix_mod_residue_file_name);
      fix_mod_residue_list_ = ResidueFactory::getResiduePtrVecInstance(
          fix_mod_residue_file_name);
      LOG_DEBUG("fix mod residue initialized ");

      LOG_DEBUG("allow prot mods initialization ");
      std::string allowed_ptm_args = arguments["n-terminal_variable_ptm"];
      if (allowed_ptm_args.length() > 0) {
        char spliter = '=';
        std::vector<std::string> allowed_ptms = prot::split(allowed_ptm_args,
                                                      spliter);
        for (unsigned int i = 0; i < allowed_ptms.size(); i++) {
          std::string mod_name = allowed_ptms[i];
          ProtModPtr ptr = ProtModFactory::getBaseProtModPtrByName(mod_name);
          allow_prot_mod_list_.push_back(ptr);
        }
      }
      else {
        xercesc::DOMElement* parent = getChildElement(root,
                                                      "allow_prot_mod_list", 0);
        int prot_mod_num = getChildCount(parent, "prot_mod");
        for (int i = 0; i < prot_mod_num; i++) {
          std::string mod_name = getChildValue(parent, "prot_mod", i);
          ProtModPtr ptr = ProtModFactory::getBaseProtModPtrByName(mod_name);
          allow_prot_mod_list_.push_back(ptr);
        }
      }
      LOG_DEBUG("allow prot mods initialized ");

      std::string activation_type = arguments["activation"];
      LOG_DEBUG("acitivation type: " << activation_type);
      activation_ptr_ = ActivationFactory::getBaseActivationPtrByName(
          activation_type);
    }
    delete doc;
  }
}

}


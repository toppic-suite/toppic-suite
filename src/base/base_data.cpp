/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */

#include <map>

#include "base/base_data.hpp"
#include "base/residue.hpp"
#include "base/prot_mod.hpp"
#include "base/ion_type.hpp"
#include "base/neutral_loss.hpp"
#include "base/activation.hpp"
#include "base/support_peak_type.hpp"
#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/file_util.hpp"

namespace prot {

void initBaseData(const std::string &exe_dir) {
  std::string conf_dir = exe_dir + FILE_SEPARATOR + CONFIG_DIR;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    std::string config_file_name = conf_dir + FILE_SEPARATOR + CONFIG_FILE_NAME;
    LOG_DEBUG("config_file_name: " << config_file_name);
    XmlDOMDocument doc(parser, config_file_name.c_str());
    xercesc::DOMElement* root = doc.getDocumentElement();
    LOG_DEBUG("root " << root);
    std::string acid_file_name = getChildValue(root, "acid_list_file_name",
                                               0);
    acid_file_name = conf_dir + FILE_SEPARATOR + acid_file_name;
    LOG_DEBUG("acid file name: " << acid_file_name);
    AcidFactory::initFactory(acid_file_name);
    LOG_DEBUG("acid initialized ");

    std::string ptm_file_name = getChildValue(root, "ptm_list_file_name", 0);
    ptm_file_name = conf_dir + FILE_SEPARATOR + ptm_file_name;
    LOG_DEBUG("ptm file name: " << ptm_file_name);
    PtmFactory::initFactory(ptm_file_name);
    LOG_DEBUG("ptm initialized");

    std::string residue_file_name = getChildValue(root,
                                                  "residue_list_file_name",
                                                  0);
    residue_file_name = conf_dir + FILE_SEPARATOR + residue_file_name;
    LOG_DEBUG("residue file name: " << residue_file_name);
    ResidueFactory::initFactory(residue_file_name);
    LOG_DEBUG("residue initialized");

    std::string trunc_file_name = getChildValue(root, "trunc_list_file_name",
                                                0);
    trunc_file_name = conf_dir + FILE_SEPARATOR + trunc_file_name;
    LOG_DEBUG("trunc file name: " << trunc_file_name);
    TruncFactory::initFactory(trunc_file_name);
    LOG_DEBUG("trunc initialized ");

    std::string prot_mod_file_name = getChildValue(root,
                                                   "prot_mod_list_file_name",
                                                   0);
    prot_mod_file_name = conf_dir + FILE_SEPARATOR + prot_mod_file_name;
    LOG_DEBUG("prot mod file name: " << prot_mod_file_name);
    ProtModFactory::initFactory(prot_mod_file_name);
    LOG_DEBUG("prot mod initialized ");

    std::string ion_type_file_name = getChildValue(root,
                                                   "ion_type_list_file_name",
                                                   0);
    ion_type_file_name = conf_dir + FILE_SEPARATOR + ion_type_file_name;
    LOG_DEBUG("ion type file name: " << ion_type_file_name);
    IonTypeFactory::initFactory(ion_type_file_name);
    LOG_DEBUG("ion type initialized ");

    std::string neutral_loss_file_name = getChildValue(
        root, "neutral_loss_list_file_name", 0);
    neutral_loss_file_name = conf_dir + FILE_SEPARATOR + neutral_loss_file_name;
    LOG_DEBUG("neutral loss file name: " << neutral_loss_file_name);
    NeutralLossFactory::initFactory(neutral_loss_file_name);
    LOG_DEBUG("neutral loss initialized ");

    std::string activation_file_name = getChildValue(
        root, "activation_list_file_name", 0);
    activation_file_name = conf_dir + FILE_SEPARATOR + activation_file_name;
    LOG_DEBUG("activation file name: " << activation_file_name);
    ActivationFactory::initFactory(activation_file_name);
    LOG_DEBUG("activation initialized ");

    std::string sp_type_file_name = getChildValue(
        root, "support_peak_type_file_name", 0);
    sp_type_file_name = conf_dir + FILE_SEPARATOR + sp_type_file_name;
    LOG_DEBUG("support_peak_type_file_name: " << sp_type_file_name);
    SPTypeFactory::initFactory(sp_type_file_name);
    LOG_DEBUG("support peak type initialized ");

    std::string build_in_fix_mod_residue_list_file_name = getChildValue(
        root, "build_in_fix_mod_residue_list_file_name", 0);
    build_in_fix_mod_residue_list_file_name = conf_dir + FILE_SEPARATOR 
        + build_in_fix_mod_residue_list_file_name;
    LOG_DEBUG("build in fix mod residue list file name: " << build_in_fix_mod_residue_list_file_name);
    FixResidueFactory::initFactory(build_in_fix_mod_residue_list_file_name);
    LOG_DEBUG("fix mod residue initialized ");
  }
}

}


// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <map>

#include "base/base_data.hpp"
#include "base/logger.hpp"
#include "base/xml_dom_document.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/file_util.hpp"

#include "base/acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/trunc_base.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/ion_type_base.hpp"
#include "base/neutral_loss_base.hpp"
#include "base/activation_base.hpp"
#include "base/support_peak_type_base.hpp"

namespace prot {

void BaseData::init(const std::string &exe_dir) {
  std::string separator = FileUtil::getFileSeparator();
  std::string base_data_dir = exe_dir + separator + BaseData::getBaseDataDir();
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    std::string config_file_name = base_data_dir + separator 
        + BaseData::getBaseDataConfigFileName();
    LOG_DEBUG("config_file_name: " << config_file_name);
    XmlDOMDocument doc(parser, config_file_name.c_str());
    xercesc::DOMElement* root = doc.getDocumentElement();
    LOG_DEBUG("root " << root);
    std::string acid_file_name = XmlDomUtil::getChildValue(root, "acid_list_file_name", 0);
    acid_file_name = base_data_dir + separator + acid_file_name;
    LOG_DEBUG("acid file name: " << acid_file_name);
    AcidBase::initBase(acid_file_name);
    LOG_DEBUG("acid initialized ");

    std::string ptm_file_name = XmlDomUtil::getChildValue(root, "ptm_list_file_name", 0);
    ptm_file_name = base_data_dir + separator + ptm_file_name;
    LOG_DEBUG("ptm file name: " << ptm_file_name);
    PtmBase::initBase(ptm_file_name);

    std::string residue_file_name = XmlDomUtil::getChildValue(root, "residue_list_file_name", 0);
    residue_file_name = base_data_dir + separator + residue_file_name;
    LOG_DEBUG("residue file name: " << residue_file_name);
    ResidueBase::initBase(residue_file_name);
    LOG_DEBUG("residue initialized");

    std::string trunc_file_name 
        = XmlDomUtil::getChildValue(root, "trunc_list_file_name", 0);
    trunc_file_name = base_data_dir + separator + trunc_file_name;
    LOG_DEBUG("trunc file name: " << trunc_file_name);
    TruncBase::initBase(trunc_file_name);
    LOG_DEBUG("trunc initialized ");

    std::string mod_file_name 
        = XmlDomUtil::getChildValue(root, "mod_list_file_name", 0);
    mod_file_name = base_data_dir + separator + mod_file_name;
    LOG_DEBUG("mod file name: " << mod_file_name);
    ModBase::initBase(mod_file_name);
    LOG_DEBUG("mod initialized ");

    std::string prot_mod_file_name
        = XmlDomUtil::getChildValue(root, "prot_mod_list_file_name", 0);
    prot_mod_file_name = base_data_dir + separator + prot_mod_file_name;
    LOG_DEBUG("prot mod file name: " << prot_mod_file_name);
    ProtModBase::initBase(prot_mod_file_name);
    LOG_DEBUG("prot mod initialized ");

    std::string ion_type_file_name 
        = XmlDomUtil::getChildValue(root, "ion_type_list_file_name", 0);
    ion_type_file_name = base_data_dir + separator + ion_type_file_name;
    LOG_DEBUG("ion type file name: " << ion_type_file_name);
    IonTypeBase::initBase(ion_type_file_name);
    LOG_DEBUG("ion type initialized ");

    std::string neutral_loss_file_name 
        = XmlDomUtil::getChildValue(root, "neutral_loss_list_file_name", 0);
    neutral_loss_file_name = base_data_dir + separator + neutral_loss_file_name;
    LOG_DEBUG("neutral loss file name: " << neutral_loss_file_name);
    NeutralLossBase::initBase(neutral_loss_file_name);
    LOG_DEBUG("neutral loss initialized ");

    std::string activation_file_name 
        = XmlDomUtil::getChildValue(root, "activation_list_file_name", 0);
    activation_file_name = base_data_dir + separator + activation_file_name;
    LOG_DEBUG("activation file name: " << activation_file_name);
    ActivationBase::initBase(activation_file_name);
    LOG_DEBUG("activation initialized ");

    std::string sp_type_file_name 
        = XmlDomUtil::getChildValue(root, "support_peak_type_file_name", 0);
    sp_type_file_name = base_data_dir + separator + sp_type_file_name;
    LOG_DEBUG("support_peak_type_file_name: " << sp_type_file_name);
    SPTypeBase::initBase(sp_type_file_name);
    LOG_DEBUG("support peak type initialized ");
  }
}

}


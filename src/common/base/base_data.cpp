//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "common/util/logger.hpp"

#include "common/base/amino_acid_base.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/trunc_base.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/ion_type_base.hpp"
#include "common/base/neutral_loss_base.hpp"
#include "common/base/activation_base.hpp"
#include "common/base/support_peak_type_base.hpp"
#include "common/base/base_data.hpp"

namespace toppic {

namespace base_data {

bool base_data_init_ = false;

void init() {
  // base data only need to be init once
  if (base_data_init_) { 
    return; 
  }

  AminoAcidBase::initBase();
  LOG_DEBUG("acid initialized ");

  PtmBase::initBase();
  LOG_DEBUG("ptm initialized");

  ResidueBase::initBase();
  LOG_DEBUG("residue initialized");

  TruncBase::initBase();
  LOG_DEBUG("trunc initialized ");

  ModBase::initBase();
  LOG_DEBUG("mod initialized ");

  ProtModBase::initBase();
  LOG_DEBUG("prot mod initialized ");

  IonTypeBase::initBase();
  LOG_DEBUG("ion type initialized ");

  NeutralLossBase::initBase();
  LOG_DEBUG("neutral loss initialized ");

  ActivationBase::initBase();
  LOG_DEBUG("activation initialized ");

  SPTypeBase::initBase();
  LOG_DEBUG("support peak type initialized ");

  base_data_init_ = true;
}

} // namespace base_data

} // namespace toppic


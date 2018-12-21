//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_SPEC_EXTEND_MS_HPP_
#define PROT_SPEC_EXTEND_MS_HPP_

#include "spec/extend_peak.hpp"
#include "spec/ms.hpp"

namespace toppic {

typedef std::shared_ptr<Ms<ExtendPeakPtr> > ExtendMsPtr;

typedef std::vector<ExtendMsPtr> ExtendMsPtrVec;

namespace extend_ms {

std::vector<double> getExtendMassVec(ExtendMsPtr extend_ms_ptr);

std::vector<std::pair<int, int> > getExtendIntMassErrorList(const ExtendMsPtrVec &ext_ms_ptr_vec,
                                                            bool pref, double scale);

std::vector<std::pair<double, double> > getExtendMassToleranceList(ExtendMsPtr extend_ms_ptr);

}  // namespace extend_ms

}  // namespace toppic

#endif

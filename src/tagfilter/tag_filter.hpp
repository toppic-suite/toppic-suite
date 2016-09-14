// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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



#ifndef PROT_TAG_FILTER
#define PROT_TAG_FILTER

#include "base/proteoform.hpp"
#include "base/residue_util.hpp"
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "tag_filter_mng.hpp"
#include "tag.hpp"

namespace prot {

typedef std::map<std::string, std::vector<std::vector<double> > > SeqTag;

class TagFilter {
 public:
  TagFilter(const ProteoformPtrVec &proteo_ptrs,
            TagFilterMngPtr mng_ptr);

  std::vector<std::string> getBestMatch(const PrmMsPtrVec &ms_ptr_vec);

 private:
  TagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  std::vector<SeqTag> seq_tag_vec_;
  std::vector<std::string> seq_name_vec_;
  std::vector<std::pair<std::string, double> > residue_mass_list_;
  std::vector<SpecTag> geneSpecTag(const PrmPeakPtrVec & prm_peaks); 

  SeqTag geneSeqTag(ProteoformPtr proteoform);
};

typedef std::shared_ptr<TagFilter> TagFilterPtr;
} /* namespace prot */

#endif

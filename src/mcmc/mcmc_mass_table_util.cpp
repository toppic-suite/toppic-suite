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

#include <utility>
#include <map>
#include <string>
#include <algorithm>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "base/residue_util.hpp"

#include "mcmc/mcmc_mass_table_util.hpp"

namespace toppic {
namespace mass_table_util {

std::map<int, std::vector<std::string> > readMassTable(MCMCMngPtr mng_ptr) {
  std::map<int, std::vector<std::string> > mass_table;
  std::string mass_table_file
      = mng_ptr->prsm_para_ptr_->getResourceDir() + "/base_data/mass_table.txt";

  std::ifstream infile(mass_table_file);
  std::string line;
  while (std::getline(infile, line)) {
    std::vector<std::string> strs;
    boost::split(strs, line, boost::is_any_of("\t"));
    int m = std::stoi(strs[0]);
    for (size_t i = 1; i < strs.size(); i++) {
      mass_table[m].push_back(strs[i]);
    }
  }
  infile.close();

  return mass_table;
}

std::map<int, std::vector<std::string> > geneMassTable(MCMCMngPtr mng_ptr) {
  if (mng_ptr->prsm_para_ptr_->getFixModPtrVec().empty()) {
    return readMassTable(mng_ptr);
  } else {
    return geneMassTableFixMod(mng_ptr);
  }
}

std::map<int, std::vector<std::string> > geneMassTableFixMod(MCMCMngPtr mng_ptr) {
  double convert_ratio = mng_ptr->convert_ratio_;

  ResiduePtrVec res_vec
      = residue_util::convertStrToResiduePtrVec("GASPVTCLNDQKEMHFURYWO",
                                                mng_ptr->prsm_para_ptr_->getFixModPtrVec());

  std::map<int, std::vector<std::string> > aa_pair_mass_map_;

  for (size_t i = 0; i < res_vec.size(); i++) {
    for (size_t j = 0; j < res_vec.size(); j++) {
      for (size_t k = 0; k < res_vec.size(); k++) {
        std::vector<std::string> aa_vec_tmp;
        aa_vec_tmp.push_back(res_vec[i]->getAminoAcidPtr()->getOneLetter());
        aa_vec_tmp.push_back(res_vec[j]->getAminoAcidPtr()->getOneLetter());
        aa_vec_tmp.push_back(res_vec[k]->getAminoAcidPtr()->getOneLetter());
        std::sort(aa_vec_tmp.begin(), aa_vec_tmp.end());
        int m = std::round((res_vec[i]->getMass() + res_vec[j]->getMass() + res_vec[k]->getMass()) * convert_ratio);
        std::string str_tmp;
        for (std::string s : aa_vec_tmp) {
          str_tmp += s;
        }
        aa_pair_mass_map_[m].push_back(str_tmp);
      }
    }
  }

  std::vector<std::pair<int, std::vector<std::string> > > aa_pair_mass_vec;

  std::vector<std::pair<int, std::vector<std::string> > > aa_pair_mass_vec2;

  for (auto it = aa_pair_mass_map_.begin(); it != aa_pair_mass_map_.end(); it++) {
    std::vector<std::string> seq_vec = it->second;
    std::sort(seq_vec.begin(), seq_vec.end());
    seq_vec.erase(std::unique(seq_vec.begin(), seq_vec.end()), seq_vec.end());
    aa_pair_mass_vec.push_back(std::make_pair(it->first, seq_vec));
  }

  for (size_t i = 0; i < aa_pair_mass_vec.size(); i++) {
    int m = aa_pair_mass_vec[i].first;
    std::vector<std::string> seq_vec = aa_pair_mass_vec[i].second;
    for (int k = i - 1; k >= 0; k--) {
      if (m - aa_pair_mass_vec[k].first > 13) {
        break;
      }
      seq_vec.insert(seq_vec.end(), aa_pair_mass_vec[k].second.begin(), aa_pair_mass_vec[k].second.end());
    }

    for (size_t k = i + 1; k < aa_pair_mass_vec.size(); k++) {
      if (aa_pair_mass_vec[k].first - m > 13) {
        break;
      }
      seq_vec.insert(seq_vec.end(), aa_pair_mass_vec[k].second.begin(), aa_pair_mass_vec[k].second.end());
    }

    aa_pair_mass_vec2.push_back(std::make_pair(m, seq_vec));
  }

  std::map<int, std::vector<std::string> > mass_table;

  for (size_t i = 0; i < aa_pair_mass_vec2.size(); i++) {
    if (aa_pair_mass_vec2[i].second.size() == 1) continue;
    mass_table[aa_pair_mass_vec2[i].first] = aa_pair_mass_vec2[i].second;
  }

  return mass_table;
}

}  // namespace mass_table_util
}  // namespace toppic

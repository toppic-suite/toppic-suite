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

#include <map>
#include <iostream>
#include <string>
#include <algorithm>

#include "base/base_data.hpp"
#include "base/file_util.hpp"
#include "base/string_util.hpp"
#include "base/fasta_util.hpp"

#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_stat.hpp"

#include "spec/msalign_util.hpp"

#include "console/topquant_argument.hpp"

using namespace prot;

void writeCurSeq(std::vector<PrsmPtr> & prsm_vec,
                 const std::vector<std::string> & spec_file_lst,
                 std::ofstream & quant_output) {
  std::sort(prsm_vec.begin(), prsm_vec.end(), [] (PrsmPtr a, PrsmPtr b) { 
              return a->getPrecFeatureInte() > b->getPrecFeatureInte();
            });

  std::vector<std::vector<PrsmPtr> > clusters;

  for (size_t k = 0; k < prsm_vec.size(); k++) {
    bool is_found = false;
    PrsmPtr cur_prsm = prsm_vec[k];

    for (size_t i = 0; i < clusters.size(); i++) {
      PrsmPtr ref_prsm = clusters[i][0]; 
      if (std::abs(cur_prsm->getAdjustedPrecMass() - ref_prsm->getAdjustedPrecMass()) < 1.2) {
        clusters[i].push_back(prsm_vec[k]);
        is_found = true;
        break; 
      }
    }

    if (!is_found) {
      std::vector<PrsmPtr> new_cluster;
      new_cluster.push_back(cur_prsm); 
      clusters.push_back(new_cluster);
    }
  }

  for (size_t k = 0; k < clusters.size(); k++) {
    quant_output << clusters[k][0]->getProteoformPtr()->getSeqName() << "\t"
        << clusters[k][0]->getProteoformPtr()->getProteinMatchSeq() << "\t"
        << clusters[k][0]->getAdjustedPrecMass();
    std::map<std::string, double> inten_map;
    for (size_t i = 0; i < clusters[k].size(); i++) {
      inten_map[clusters[k][i]->getFileName()] = clusters[k][i]->getPrecFeatureInte();
    }

    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      if (inten_map[spec_file_lst[i]] > 0 ) {
        quant_output << "\t" << inten_map[spec_file_lst[i]];
      } else {
        quant_output << "\t" << "-" ;
      }
    }

    quant_output << std::endl;
  }
}

int main(int argc, char* argv[]) {
  Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);

  if (!success) {
    return EXIT_FAILURE;
  }

  std::map<std::string, std::string> arguments = argu_processor.getArguments();

  std::string exe_dir = file_util::getExecutiveDir(argv[0]);

  arguments["executiveDir"] = exe_dir;

  std::vector<std::string> spec_file_lst = argu_processor.getSpecFileList();

  std::string resource_dir = arguments["resourceDir"];

  base_data::init(resource_dir);

  std::string db_file_name = arguments["databaseFileName"];
  std::string ori_db_file_name = arguments["oriDatabaseFileName"];

  int db_block_size = std::stoi(arguments["databaseBlockSize"]);

  fasta_util::dbPreprocess(ori_db_file_name, db_file_name, false, db_block_size);

  PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);

  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  std::vector<PrsmPtr> prsm_vec;

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    std::string input_file_name = prot::file_util::basename(spec_file_lst[k]) + ".FORM_RESULT";
    PrsmReader prsm_reader(input_file_name);
    PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
    while (prsm_ptr != nullptr) {
      prsm_vec.push_back(prsm_ptr); 
      prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
    }
    prsm_reader.close();
  }

  std::sort(prsm_vec.begin(), prsm_vec.end(), [] (PrsmPtr a, PrsmPtr b) { 
              if (a->getProteoformPtr()->getSeqName() < b->getProteoformPtr()->getSeqName()) {
                return true; 
              } else if (a->getProteoformPtr()->getSeqName() > b->getProteoformPtr()->getSeqName()) {
                return false;
              }

              if (a->getAdjustedPrecMass() < b->getAdjustedPrecMass()) {
                return true;
              } else if (a->getAdjustedPrecMass() > b->getAdjustedPrecMass()) {
                return false;
              }

              return a->getFileName() < b->getFileName();
            });

  std::ofstream quant_output;
  quant_output.open(arguments["combinedOutputName"] + ".QUANT_OUTPUT_TABLE");

  quant_output << "Protein name" << "\t"
      << "Proteoform" << "\t"
      << "Proteoform mass";

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    boost::filesystem::path p(spec_file_lst[k]);
    spec_file_lst[k] = p.stem().string();
    spec_file_lst[k] += p.extension().string();
    quant_output << "\t" << spec_file_lst[k];
  }

  quant_output << std::endl;

  std::string cur_seq_name = prsm_vec[0]->getProteoformPtr()->getSeqName();

  std::vector<PrsmPtr> cur_prsm_vec;

  for (size_t k = 0; k < prsm_vec.size(); k++) {
    if (prsm_vec[k]->getProteoformPtr()->getSeqName() == cur_seq_name) {
      cur_prsm_vec.push_back(prsm_vec[k]);
    } else {
      writeCurSeq(cur_prsm_vec, spec_file_lst, quant_output);
      cur_prsm_vec.clear();
      cur_seq_name = prsm_vec[k]->getProteoformPtr()->getSeqName();
      cur_prsm_vec.push_back(prsm_vec[k]);
    } 
  }

  writeCurSeq(cur_prsm_vec, spec_file_lst, quant_output);

  quant_output.close();

  return EXIT_SUCCESS;
}

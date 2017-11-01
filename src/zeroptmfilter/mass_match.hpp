//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef ZERO_PTM_FILTER_MASS_MATCH_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"

namespace prot {

class MassMatch {
 public:
  MassMatch(std::vector<std::vector<int>> &mass_2d, 
            std::vector<std::vector<double>> &real_shift_2d,
            std::vector<std::vector<int>> &pos_2d,
            double max_proteoform_mass, double scale);  

  ~MassMatch();

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  std::vector<short> &scores);

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  int start, double shift, std::vector<short> &scores);

  void compMatchScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                       const std::pair<int,int> &prec_minus_water_mass_error, 
                       std::vector<short> &scores);

  int getRowNum() {return row_num_;}

  static int getPrecursorMatchScore() {return 10000;}

  std::vector<int>& getProteoRowBegins() {return proteo_row_begins_;}
  std::vector<int>& getProteoRowEnds() {return proteo_row_ends_;}
  std::vector<double>& getTruncShifts() {return trunc_shifts_;}

 private:
  double scale_;
  int proteo_num_;
  int col_num_;
  int row_num_;

  // the first row of each proteoform  
  std::vector<int> proteo_row_begins_;
  std::vector<int> proteo_row_ends_;
  // the proteoform id of each row
  std::vector<int> row_proteo_ids_;
  std::vector<double> trunc_shifts_;

  std::vector<int> col_index_begins_;
  std::vector<int> col_index_ends_;
  std::vector<int> col_indexes_;

  void initProteoformBeginEnds(std::vector<std::vector<int>> &mass_2d,
                               std::vector<std::vector<double>> &shift_2d);

  void initIndexes(std::vector<std::vector<int>> &mass_2d,
                   std::vector<std::vector<double>> &real_shift_2d,
                   std::vector<std::vector<int>> &pos_2d);     

  void compColumnMatchNums(std::vector<std::vector<int>> &mass_2d,
                           std::vector<std::vector<int>> &shift_2d,     
                           std::vector<std::vector<int>> &pos_2d,     
                           std::vector<int> &col_match_nums);

  void fillColumnIndex(std::vector<std::vector<int>> &mass_2d,
                       std::vector<std::vector<int>> &shift_2d,     
                       std::vector<std::vector<int>> &pos_2d,     
                       std::vector<int> &col_index_pnts);
};

typedef std::shared_ptr<MassMatch> MassMatchPtr;

} /* namespace prot */

#endif /* COMP_SHIFT_HPP_ */

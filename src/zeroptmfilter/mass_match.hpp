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

#ifndef ZERO_PTM_FILTER_MASS_MATCH_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"

namespace prot {

class MassMatch {
 public:
  MassMatch(const ProteoformPtrVec &proteo_ptrs, int scale, double
            max_proteoform_mass);  

  MassMatch(const ProteoformPtrVec &proteo_ptrs, 
            const std::vector<std::vector<int>> &shift_2d,
            int scale, double max_proteoform_mass);

  ~MassMatch();

  static FilterProteinPtrVec geneResults(std::vector<std::pair<int,int>>
                                         &single_type_results, double threshold,
                                         int single_type_num);

  FilterProteinPtrVec findTopProteins(std::vector<short> &scores, int num);

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  std::vector<short> &scores);

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  int start, std::vector<short> &scores);

  int getRowNum() {return row_num_;}

 private:
  int scale_;
  int proteo_num_;
  int col_num_;
  int row_num_;

  // the first row of each proteoform  
  std::vector<int> proteo_row_begins_;
  std::vector<int> proteo_row_ends_;
  // the proteoform id of each row
  std::vector<int> row_proteo_ids_;
  //std::vector<double> n_trunc_shifts_;
  //std::vector<double> c_trunc_shifts_;

  std::vector<int> col_index_begins_;
  std::vector<int> col_index_ends_;
  std::vector<int> col_indexes_;

  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs,
                               const std::vector<int> &index_nums);

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, 
                             std::vector<int> &col_match_nums);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs);

  void updateColumnMatchNums(ProteoformPtr proteo_ptr,
                             const std::vector<int> &int_shifts,
                             std::vector<int> &col_match_nums);

  void initIndexes(const ProteoformPtrVec &proteo_ptrs, 
                   const std::vector<std::vector<int>> &shift_2d);
};

typedef std::shared_ptr<MassMatch> MassMatchPtr;

} /* namespace prot */

#endif /* COMP_SHIFT_HPP_ */

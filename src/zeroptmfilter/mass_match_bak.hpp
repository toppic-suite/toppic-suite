#ifndef ZERO_PTM_FILTER_MASS_MATCH_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"

namespace prot {

class MassMatch {
 public:
  MassMatch(const ProteoformPtrVec &proteo_ptrs, int scale, 
            double max_proteoform_mass, bool rev);  

  MassMatch(const ProteoformPtrVec &proteo_ptrs, 
            const std::vector<std::vector<int>> &shift_2d,
            int scale, double max_proteoform_mass, bool rev);

  ~MassMatch();

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  std::vector<short> &scores);

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  int start, std::vector<short> &scores);


  void compMatchScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                       const std::pair<int,int> &prec_minus_water_mass_error, 
                       std::vector<short> &scores);

  int getRowNum() {return row_num_;}

  static int getPrecursorMatchScore() {return 10000;}

  std::vector<int>& getProteoRowBegins() {return proteo_row_begins_;}
  std::vector<int>& getProteoRowEnds() {return proteo_row_ends_;}

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
  std::vector<double> trunc_shifts_;

  std::vector<int> col_index_begins_;
  std::vector<int> col_index_ends_;
  std::vector<int> col_indexes_;

  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs,
                               const std::vector<int> &index_nums);

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, 
                             std::vector<int> &col_match_nums, 
                             bool rev);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs, bool rev);

  void updateColumnMatchNums(ProteoformPtr proteo_ptr,
                             const std::vector<int> &int_shifts,
                             std::vector<int> &col_match_nums, bool rev);

  void initIndexes(const ProteoformPtrVec &proteo_ptrs, 
                   const std::vector<std::vector<int>> &shift_2d,
                   bool rev);
};

typedef std::shared_ptr<MassMatch> MassMatchPtr;

} /* namespace prot */

#endif /* COMP_SHIFT_HPP_ */

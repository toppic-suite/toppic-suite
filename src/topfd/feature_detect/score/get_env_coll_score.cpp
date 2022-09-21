//
// Created by abbash on 9/21/22.
//

#include "get_env_coll_score.hpp"

void generateTensors(std::vector<fdeep::tensors> &tensorsL, const std::vector<double> &matrix) {
  fdeep::tensor_shape tensor_shape(8);
  fdeep::tensor t(tensor_shape, 0.0f);
  for (int y = 0; y < 8; ++y)
    t.set(0, 0, 0, 0, y, matrix[y]);

  std::vector<fdeep::tensor> tensors;
  tensors.push_back(t);
  tensorsL.push_back(tensors);
}

double toppic::env_coll_score::get_env_coll_score(fdeep::model& model, std::vector<double>& data) {
  std::vector<fdeep::tensors> tensorsL;
//    write_out_files::write_env_cnn_matrix(envcnn_data_matrix);
  generateTensors(tensorsL, data);
  if (!tensorsL.empty()) {
    std::vector<fdeep::tensors> pred_scores = model.predict_multi(tensorsL, false);
    return pred_scores[0][0].get(0, 0, 0, 0, 0);
  }
  return 0.0;

}
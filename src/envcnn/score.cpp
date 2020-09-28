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

#include <cmath>
#include <src/envcnn/fdeep/tensor5.hpp>
#include "generate_matrix.hpp"
#include "ms/spec/baseline_util.hpp"
#include "score.hpp"
#include "common/util/file_util.hpp"

namespace toppic {

    MatchEnvPtrVec MatchEnvFilterCNN::filter_using_cnn(MatchEnvPtrVec &ori_envs, PeakPtrVec &peak_list, fdeep::model model) {
        std::vector<fdeep::tensor5s> tensorsL = GenerateMatrix::getTensor(ori_envs, peak_list);
        if (!tensorsL.empty()){
            std::vector<fdeep::tensor5s> pred_scores = model.predict_multi(tensorsL, true);
            for (size_t i = 0; i < ori_envs.size(); i++) {
                ori_envs[i]->setScore(pred_scores[i][0].get(0, 0, 0, 0, 0));
            }
        }
        return ori_envs;
    }

    fdeep::model MatchEnvFilterCNN::loadModel(std::string path) {
        std::string model_file = path + file_util::getFileSeparator() + "my_model.json";
        const auto model = fdeep::load_model(model_file);
        return model;
    }
}


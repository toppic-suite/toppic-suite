// Copyright 2016, Tobias Hermann.
// https://github.com/Dobiasd/frugally-deep
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "activation_layer.hpp"

#include <string>

namespace fdeep { namespace internal
{

class elu_layer : public activation_layer
{
public:
    explicit elu_layer(const std::string& name, float_type alpha)
        : activation_layer(name), alpha_(alpha)
    {
    }
protected:
    float_type alpha_;
    static float_type activation_function(float_type alpha, float_type x)
    {
        return x >= 0 ? x : alpha * (std::exp(x) - 1);
    }
    tensor5 transform_input(const tensor5& in_vol) const override
    {
        return transform_tensor5(
            fplus::bind_1st_of_2(activation_function, alpha_),
            in_vol);
    }
};

} } // namespace fdeep, namespace internal
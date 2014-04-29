/*
 * transformer.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: xunlikun
 */

#ifndef TRANSFORMER_HPP_
#define TRANSFORMER_HPP_

#include <memory>
#include <xalanc/Include/PlatformDefinitions.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xalanc/XalanTransformer/XalanTransformer.hpp>


class Transformer {
 public:
  Transformer();
  virtual ~Transformer();
  void translate();
};

typedef std::shared_ptr<Transformer> TransformerPtr;

#endif /* TRANSFORMER_HPP_ */

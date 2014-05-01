/*
 * transformer.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: xunlikun
 */

#ifndef TRANSFORMER_HPP_
#define TRANSFORMER_HPP_

#include <memory>
#include <map>
#include <string>
#include <xalanc/Include/PlatformDefinitions.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xalanc/XalanTransformer/XalanTransformer.hpp>

namespace prot {

void translate(std::map<std::string,std::string> arguments);

}

#endif /* TRANSFORMER_HPP_ */

#ifndef PROT_TRANSFORMER_HPP_
#define PROT_TRANSFORMER_HPP_

#include <memory>
#include <map>
#include <string>
#include <xercesc/util/PlatformUtils.hpp>
#include <xalanc/Include/PlatformDefinitions.hpp>
#include <xalanc/XalanTransformer/XalanTransformer.hpp>

namespace prot {

void translate(std::map<std::string,std::string> &arguments);

}

#endif /* PROT_TRANSFORMER_HPP_ */

//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"

#include "ms/feature/single_charge_feature.hpp"

namespace toppic {

  SingleChargeFeature::SingleChargeFeature() :
      charge_(-1),
      time_begin_(-1),
      time_end_(-1),
      scan_begin_(-1),
      scan_end_(-1),
      spec_id_begin_(-1),
      spec_id_end_(-1),
      intensity_(-1),
      env_num_(0),
      mass_(-1) {
  }

  SingleChargeFeature::SingleChargeFeature(int charge,
                                           double time_begin, double time_end,
                                           int scan_begin, int scan_end,
                                           double intensity, int env_num,
                                           int spec_id_begin, int spec_id_end,
                                           double mass, std::vector<double> xic_inte,
                                           std::vector<double> envelopeMass, std::vector<double> aggregateEnvelopeInte) :
      charge_(charge),
      time_begin_(time_begin),
      time_end_(time_end),
      scan_begin_(scan_begin),
      scan_end_(scan_end),
      spec_id_begin_(spec_id_begin),
      spec_id_end_(spec_id_end),
      intensity_(intensity),
      env_num_(env_num),
      mass_(mass) {
    for (auto inte: xic_inte)
      xic_inte_.push_back(inte);

    for (auto inte: envelopeMass)
      envelopeMass_.push_back(inte);

    for (auto inte: aggregateEnvelopeInte)
      aggregateEnvelopeInte_.push_back(inte);
  }

  SingleChargeFeature::SingleChargeFeature(int charge,
                                           double time_begin, double time_end,
                                           int scan_begin, int scan_end,
                                           double intensity, int env_num) :
      charge_(charge),
      time_begin_(time_begin),
      time_end_(time_end),
      scan_begin_(scan_begin),
      scan_end_(scan_end),
      intensity_(intensity),
      env_num_(env_num) {
  }

  SingleChargeFeature::SingleChargeFeature(std::shared_ptr<SingleChargeFeature> f) {
    charge_ = f->getCharge();
    time_begin_ = f->getTimeBegin();
    time_end_ = f->getTimeEnd();
    scan_begin_ = f->getScanBegin();
    scan_end_ = f->getScanEnd();
    spec_id_begin_ = f->getSpecIDBegin();
    spec_id_end_ = f->getSpecIDEnd();
    intensity_ = f->getIntensity();
    mass_ = f->getMass();
    env_num_ = f->getEnvNum();
    for (auto inte: f->getXicInte())
      xic_inte_.push_back(inte);

    for (auto inte: f->getAggregateEnvelopeInte())
      aggregateEnvelopeInte_.push_back(inte);

    for (auto inte: f->getEnvelopeMass())
      envelopeMass_.push_back(inte);
  }

  SingleChargeFeature::SingleChargeFeature(XmlDOMElement *element) {
    charge_ = xml_dom_util::getIntChildValue(element, "charge", 0);
    time_begin_ = xml_dom_util::getDoubleChildValue(element, "time_begin", 0);
    time_end_ = xml_dom_util::getDoubleChildValue(element, "time_end", 0);
    scan_begin_ = xml_dom_util::getIntChildValue(element, "scan_begin", 0);
    scan_end_ = xml_dom_util::getIntChildValue(element, "scan_end", 0);
    intensity_ = xml_dom_util::getDoubleChildValue(element, "intensity", 0);
    env_num_ = xml_dom_util::getIntChildValue(element, "envelope_num", 0);
  }

  void SingleChargeFeature::appendToXml(XmlDOMDocument *xml_doc, XmlDOMElement *parent) {
    std::string element_name = SingleChargeFeature::getXmlElementName();
    XmlDOMElement *element = xml_doc->createElement(element_name.c_str());
    std::string str = str_util::toString(charge_);
    xml_doc->addElement(element, "charge", str.c_str());
    str = str_util::toString(time_begin_);
    xml_doc->addElement(element, "time_begin", str.c_str());
    str = str_util::toString(time_end_);
    xml_doc->addElement(element, "time_end", str.c_str());
    str = str_util::toString(scan_begin_);
    xml_doc->addElement(element, "scan_begin", str.c_str());
    str = str_util::toString(scan_end_);
    xml_doc->addElement(element, "scan_end", str.c_str());
    str = str_util::toString(intensity_);
    xml_doc->addElement(element, "intensity", str.c_str());
    str = str_util::toString(env_num_);
    xml_doc->addElement(element, "envelope_num", str.c_str());
    parent->appendChild(element);
  }

}

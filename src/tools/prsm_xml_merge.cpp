//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>

#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"

int main(int argc, const char *argv[]) {
  try {
    namespace po = boost::program_options;
    po::options_description desc("Options");
    int N = 100000;
    std::string output_filename = "";
    std::vector<std::string> prsm_xml_list;
    desc.add_options()
        ("help,h", "Help message")
        ("num,N", po::value<int>(&N), "N")
        ("prsm-xml,p", po::value<std::vector<std::string> >()->multitoken()->required(), "msalign file list")
        ("output,o", po::value<std::string>(&output_filename)->required(), "output msalign file name");
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
    if (vm.count("help")) {
      std::cout << desc << "\n";
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    if (vm.count("prsm-xml")) {
      prsm_xml_list = vm["prsm-xml"].as<std::vector<std::string> >();
      if (prsm_xml_list.size() == 1) {
        return EXIT_SUCCESS;
      }
    }

    prot::PrsmXmlWriterPtr prsm_writer = std::make_shared<prot::PrsmXmlWriter>(output_filename);

    std::cout << "Merging ";
    for (size_t i = 0; i < prsm_xml_list.size(); i++) {
      std::cout << prsm_xml_list[i] << " ";
      prot::PrsmReaderPtr prsm_reader = std::make_shared<prot::PrsmReader>(prsm_xml_list[i]);
      prot::PrsmStrPtr prsm = prsm_reader->readOnePrsmStr();
      while (prsm != nullptr) {
        prsm->setSpectrumId(i * N + prsm->getSpectrumId());
        if (prsm->getPrecFeatureId() > 0) {
          prsm->setPrecFeatureId(i * N + prsm->getPrecFeatureId());
        }
        if (prsm->getPrecursorId() > 0) {
          prsm->setPrecursorId(i * N + prsm->getPrecursorId());
        }
        prsm_writer->write(prsm);
        prsm = prsm_reader->readOnePrsmStr(); 
      }
    }
    std::cout << std::endl;
    prsm_writer->close();

  } catch (const std::exception &ex) {
    std::cerr << ex.what() << std::endl;
  }
  return EXIT_SUCCESS;
}


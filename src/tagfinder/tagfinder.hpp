
#ifndef __TAGFINDER_HPP__
#define __TAGFINDER_HPP__

#include "base/residue_util.hpp"
#include "prsm/prsm_para.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_util.hpp"

#include "peak_node.hpp"

namespace prot {

const int TAG_SIZE = 4;

class Tagfinder{

 public:
  Tagfinder(PrsmParaPtr prsm_para_ptr, int gap);
  void process();

 private:
  std::vector<PeakNodePtr> generateGapEdges(const PrmPeakPtrVec & prm_peaks, int gap);
  std::vector<std::string> getTags(std::vector<std::vector<PeakNodePtr> > componentsFromGraph);
  std::vector<std::pair<std::string, double>> mass_list_;
  PrsmParaPtr prsm_para_ptr_;
  int gap_;

};

typedef std::shared_ptr<Tagfinder> TagfinderPtr;

std::vector<std::vector<PeakNodePtr>> getComponentsFromGraph(std::vector<PeakNodePtr> peaks);

std::vector<PeakNodePtr> findBestTag(PeakNodePtr peak, std::vector<PeakNodePtr> best, int len, std::vector<PeakNodePtr> & prefix);

std::vector<PeakNodePtr> findBestTag(std::vector<PeakNodePtr> peaks);

}

#endif

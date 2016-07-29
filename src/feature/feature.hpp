#ifndef PROT_FEATURE_FEATURE_HPP_
#define PROT_FEATURE_FEATURE_HPP_

#include <memory>
#include <vector>

namespace prot {

class Feature {
 public:
  Feature(int id, double mono_mass, double inte,
          int scan_begin, int scan_end): 
      id_(id),
      mono_mass_(mono_mass),
      intensity_(inte),
      scan_begin_(scan_begin),
      scan_end_(scan_end) {
      }

  int getId() {return id_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  int getScanBegin() {return scan_begin_;}

  int getScanEnd() {return scan_end_;}


 private:
  int id_;
  double mono_mass_;
  double intensity_;
  int scan_begin_;
  int scan_end_;
};

typedef std::shared_ptr<Feature> FeaturePtr;
typedef std::vector<FeaturePtr> FeaturePtrVec;

}
#endif

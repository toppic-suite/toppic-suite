#include <map>
#include "base/logger.hpp"
#include "spec/ms.hpp"
#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

PrSM::PrSM(ProteoformPtr proteoform_ptr, DeconvMsPtr deconv_ms_ptr, 
           double adjusted_prec_mass, double calibration, 
           SpParaPtr sp_para_ptr) {
  proteoform_ptr_ = proteoform_ptr;
  deconv_ms_ptr_ = deconv_ms_ptr;
  MsHeaderPtr header_ptr = deconv_ms_ptr->getHeaderPtr();
  spectrum_id_ = header_ptr->getId();
  spectrum_scan_ = header_ptr->getScansString();
  precursor_id_ = header_ptr->getPrecId();
  ori_prec_mass_ = header_ptr->getPrecMonoMass();
  adjusted_prec_mass_ = adjusted_prec_mass;
  calibration_ = calibration;
  //sp_para_ptr_=sp_para_ptr;
  init(sp_para_ptr);
}

void PrSM::init(SpParaPtr sp_para_ptr) {
  double delta = adjusted_prec_mass_ - ori_prec_mass_;
  refine_ms_three_ = getMsThree(deconv_ms_ptr_, delta, sp_para_ptr);
  refine_ms_three_-> recalibrate(calibration_);
  initScores(sp_para_ptr);
}

void PrSM::initScores(SpParaPtr sp_para_ptr) {
  // refined one 
  PeakIonPairPtrVec pairs = getPeakIonPairs (proteoform_ptr_, refine_ms_three_, 
                                             sp_para_ptr->getMinMass());
  //LOG_DEBUG("peak ion pair size " << pairs.size());
  match_peak_num_ = 0;
  match_fragment_num_ = 0;
  TheoPeakPtr prev_ion(nullptr);;
  for (unsigned int i = 0; i < pairs.size(); i++) {
//    match_peak_num_ += pairs[i]->getRealPeakPtr()->getScore();
    if (pairs[i]->getTheoPeakPtr() != prev_ion) {
      prev_ion = pairs[i]->getTheoPeakPtr();
      match_fragment_num_ += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  std::sort(pairs.begin(),pairs.end(),peakIonPairUp);
  DeconvPeakPtr prev_deconv_peak(nullptr);
  for (unsigned int i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr() != prev_deconv_peak) {
      prev_deconv_peak = pairs[i]->getRealPeakPtr()->getBasePeakPtr();
      match_peak_num_ += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
}

xercesc::DOMElement* PrSM::toXmlElement(XmlDOMDocument* xml_doc){
	xercesc::DOMElement* element = xml_doc->createElement("prsm");
  //LOG_DEBUG("Element created");
	std::string str = convertToString(prsm_id_);
	xml_doc->addElement(element, "prsm_id", str.c_str());
	str = convertToString(spectrum_id_);
	xml_doc->addElement(element, "spectrum_id", str.c_str());
	xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
	str = convertToString(precursor_id_);
	xml_doc->addElement(element, "precursor_id", str.c_str());
	str = convertToString(ori_prec_mass_);
	xml_doc->addElement(element, "ori_prec_mass", str.c_str());
	str = convertToString(adjusted_prec_mass_);
	xml_doc->addElement(element, "adjusted_prec_mass", str.c_str());
	str = convertToString(calibration_);
	xml_doc->addElement(element, "calibration", str.c_str());
	str = convertToString(fdr_);
	xml_doc->addElement(element, "fdr", str.c_str());
	str = convertToString(match_peak_num_);
	xml_doc->addElement(element, "match_peak_num", str.c_str());
	str = convertToString(match_fragment_num_);
	xml_doc->addElement(element, "match_fragment_num", str.c_str());
	proteoform_ptr_->appendXml(xml_doc,element);
	if(prob_ptr_!=nullptr){
	  prob_ptr_->appendXml(xml_doc,element);
	}
	if(deconv_ms_ptr_!=nullptr){
	  deconv_ms_ptr_->appendXml(xml_doc,element);
	}
  return element;
}

void PrSM::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = toXmlElement(xml_doc);
  parent->appendChild(element);
}

PrSM::PrSM(xercesc::DOMElement* element,ProteoformPtrVec proteoforms){
  prsm_id_=getIntChildValue(element, "prsm_id", 0);
  spectrum_id_=getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_=getChildValue(element, "spectrum_scan", 0);
  precursor_id_=getIntChildValue(element, "precursor_id", 0);
  ori_prec_mass_=getDoubleChildValue(element, "ori_prec_mass", 0);
  adjusted_prec_mass_=getDoubleChildValue(element, "adjusted_prec_mass", 0);
  calibration_=getDoubleChildValue(element, "calibration", 0);
  fdr_=getDoubleChildValue(element, "fdr", 0);
  match_peak_num_=getDoubleChildValue(element, "match_peak_num", 0);
  match_fragment_num_=getDoubleChildValue(element, "match_fragment_num", 0);

  xercesc::DOMElement* proteoform_element
      = getChildElement(element,"proteoform",0);
  proteoform_ptr_ 
      = ProteoformPtr(new Proteoform(proteoform_element,proteoforms));

  int prob_count = getChildCount(element,"extreme_value");
  if(prob_count!=0){
    xercesc::DOMElement* prob_element 
        = getChildElement(element,"extreme_value",0);
    prob_ptr_ = ExtremeValuePtr(new ExtremeValue(prob_element));
  }

  xercesc::DOMElement* deconv_ms_element = getChildElement(element,"ms",0);
  xercesc::DOMElement* header_element
      = getChildElement(deconv_ms_element,"ms_header",0);
  MsHeaderPtr header_ptr  = MsHeaderPtr (new MsHeader(header_element));
//  xercesc::DOMElement* peak_element
//      = getChildElement(deconv_ms_element,"peaks",0);
  DeconvPeakPtrVec peaks;
//  int peak_num = getChildCount(peak_element,"deconv_peak");
//  for(int i=0;i<peak_num;i++){
//    xercesc::DOMElement* cur_ms_element
//        = getChildElement(deconv_ms_element,"deconv_peak",i);
//    peaks.push_back(DeconvPeakPtr(new DeconvPeak(cur_ms_element)));
//  }
  deconv_ms_ptr_ = DeconvMsPtr(new Ms<DeconvPeakPtr>(header_ptr,peaks));
  refine_ms_three_ = ExtendMsPtr(new Ms<ExtendPeakPtr>(header_ptr));
}

PrSMPtrVec readPrsm(std::string file_name,ProteoformPtrVec proteoforms){
  PrSMPtrVec results;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if(parser){
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name.c_str());
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      int simple_prsm_num = getChildCount(root, "prsm");
      for (int i = 0; i < simple_prsm_num; i++) {
        xercesc::DOMElement* prsm_element = getChildElement(root, "prsm", i);
        results.push_back(PrSMPtr(new PrSM(prsm_element,proteoforms)));

      }
    }
    delete doc;
  }
  return results;
}

bool isMatch(PrSMPtr prsm_ptr, MsHeaderPtr header_ptr) {
  int id = header_ptr->getId();
  std::string scan = header_ptr->getScansString();
  int prec_id = header_ptr->getPrecId();
  double prec_mass = header_ptr->getPrecMonoMass();
  if (id == prsm_ptr->getSpectrumId() && prec_id == prsm_ptr->getPrecurorId()) {
    LOG_DEBUG("scan " << scan << " prsm scan " << prsm_ptr->getSpectrumScan()
              << " prec mass " << prec_mass << " prsm mass " << prsm_ptr->getOriPrecMass());
    /*
    if (scan != prsm_ptr->getSpectrumScan()
        || prec_mass != prsm_ptr->getOriPrecMass()) {
      LOG_ERROR("Error in PrSM.");
    }
    */
    return true;
  } else {
    return false;
  }
}

void filterPrsms(PrSMPtrVec &prsms, MsHeaderPtr header_ptr, 
                 PrSMPtrVec &sele_prsms) {
  for (unsigned int i = 0; i < prsms.size(); i++) {
    if (isMatch(prsms[i], header_ptr)) {
      sele_prsms.push_back(prsms[i]);
    }
  }
}

}

/*
    public PrSM(Element element) throws Exception {
        prsmId = Integer.parseInt(element.getChildText("prsm_id"));
        spectrumId = Integer.parseInt(element.getChildText("spectrum_id"));
        spectrumScan = element.getChildText("spectrum_scan");
        precursorId = Integer.parseInt(element.getChildText("precursor_id"));
        OriPrecMass = Double.parseDouble(element
                .getChildText("original_precursor_mass"));
        adjustedPrecMass = Double.parseDouble(element
                .getChildText("adjusted_precursor_mass"));
        calibration = Double.parseDouble(element.getChildText("calibration"));
        Element probElement = element.getChild("probability");
        prob = new ExtremeValueProb(probElement);
        fdr = Double.parseDouble(element.getChildText("fdr"));
        Element protElement = element.getChild("annotated_protein");
        annoProtein = new AnnoProtein(protElement);
        Element spParaElement = element.getChild("sp_para");
        spPara = new SpPara(spParaElement);
        
    }

    public void process(Ms<DeconvPeak> ms, BpSpec seqs[])
            throws Exception {
        annoProtein.process(seqs);
        this.deconvMs = ms;
        if (!ms.getHeader().getScansString().equals(spectrumScan)
                || ms.getHeader().getPrecMonoMass() != OriPrecMass) {
            logger.error("Incorrect spectrum.");
            System.exit(1);
        }
        init();
    }

    public PeakIonPair[] getMatchedPairs() throws Exception {
        return annoProtein.getMatchPeak(refineMsThree, spPara.getMinMass());
    }


    public void outputMatchFragmentIon(PrintWriter writer) throws Exception {
        PeakIonPair pairs[] = annoProtein.getMatchPeak(refineMsThree,
                spPara.getMinMass());
        TheoPeak prevIon = null;
        int bp[] = new int[annoProtein.getSeq().getResSeq().getLen() - 1];
        for (int j = 0; j < pairs.length; j++) {
            nMatchPeak += pairs[j].getRealPeak().getScore();
            if (pairs[j].getTheoPeak() != prevIon) {
                prevIon = pairs[j].getTheoPeak();
                int pos = prevIon.getIon().getPos();
                if (prevIon.getIon().getIonType().isNTerm()) {
                    if (bp[pos] == 0) {
                        bp[pos] = 1;
                    } else if (bp[pos] == 2) {
                        bp[pos] = 3;
                    }
                } else {
                    if (bp[pos] == 0) {
                        bp[pos] = 2;
                    } else if (bp[pos] == 1) {
                        bp[pos] = 3;
                    }
                }
            }

        }
        for (int i = 0; i < bp.length; i++) {
            writer.print(bp[i]);
        }
    }

*/

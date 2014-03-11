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
  sp_para_ptr_=sp_para_ptr;
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
  PeakIonPairPtrVec pairs;
  getPeakIonPairs (proteoform_ptr_, refine_ms_three_, 
                   sp_para_ptr->getMinMass(), pairs);
  //LOG_DEBUG("peak ion pair size " << pairs.size());
  match_peak_num_ = 0;
  match_fragment_num_ = 0;
  TheoPeakPtr prev_ion(nullptr);;
  for (unsigned int i = 0; i < pairs.size(); i++) {
    match_peak_num_ += pairs[i]->getRealPeakPtr()->getScore();
    if (pairs[i]->getTheoPeakPtr() != prev_ion) {
      prev_ion = pairs[i]->getTheoPeakPtr();
      match_fragment_num_ += pairs[i]->getRealPeakPtr()->getScore();
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
//	proteoform_ptr_->appendXml(xml_doc,element);
//	if(prob_ptr_!=nullptr){
//	  prob_ptr_->appendXml(xml_doc,element);
//  }
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

//  xercesc::DOMElement* sp_para_element = getChildElement(element,"sp_para",0);
//  sp_para_ptr_ = SpParaPtr(new SpPara(sp_para_element));
  sp_para_ptr_ = SpParaPtr(new SpPara(10,50,nullptr,nullptr,nullptr));

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

xercesc::DOMElement* genePrSMView(XmlDOMDocument* xml_doc,PrSMPtr prsm){
  xercesc::DOMElement* element = xml_doc->createElement("prsm");
  std::string str = convertToString(prsm->getId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
  if(prsm->getProbPtr()->getPValue() != -std::numeric_limits<double>::max()){
    str=convertToString(prsm->getProbPtr()->getPValue());
    xml_doc->addElement(element, "p_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "p_value", "-1");
  }
  if(prsm->getProbPtr()->getEValue() != -std::numeric_limits<double>::max()){
    str=convertToString(prsm->getProbPtr()->getEValue());
    xml_doc->addElement(element, "e_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "e_value", "-1");
  }
  str=convertToString(prsm->getFdr());
  xml_doc->addElement(element, "fdr", str.c_str());
  str=convertToString(prsm->getMatchFragNum());
  xml_doc->addElement(element, "matched_fragment_number", str.c_str());
  str=convertToString(prsm->getMatchPeakNum());
  xml_doc->addElement(element, "matched_peak_number", str.c_str());
  str=convertToString(prsm->getAdjustedPrecMass());
  xml_doc->addElement(element, "adjusted_precursor_mass", str.c_str());
  str=convertToString(prsm->getCalibration());
  xml_doc->addElement(element, "calibration", str.c_str());

  //get ion_pair
  PeakIonPairPtrVec pairs;
  getPeakIonPairs (prsm->getProteoformPtr(), prsm->getRefineMs(),
                   prsm->getMinMass(), pairs);
  //peaks to view
  xercesc::DOMElement* ms_element = xml_doc->createElement("ms");
  prsm->getDeconvMsPtr()->getHeaderPtr()->appendXml(xml_doc,ms_element);//attention
  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
  ms_element->appendChild(peaks);
  for(unsigned int i=0;i<prsm->getDeconvMsPtr()->size();i++){
    xercesc::DOMElement* peak = xml_doc->createElement("peak");
    peaks->appendChild(peak);
    DeconvPeakPtr dp = prsm->getDeconvMsPtr()->getPeakPtr(i);
    str=convertToString(dp->getId());
    xml_doc->addElement(peak, "id", str.c_str());
    double mass = dp->getPosition();
    int charge = dp->getCharge();
    str=convertToString(mass);
    xml_doc->addElement(peak, "monoisotopic_mass", str.c_str());
    str=convertToString(mass/charge+MassConstant::getProtonMass());
    xml_doc->addElement(peak, "monoisotopic_mz", str.c_str());
    str=convertToString(dp->getIntensity());
    xml_doc->addElement(peak, "intensity", str.c_str());
    str=convertToString(charge);
    xml_doc->addElement(peak, "charge", str.c_str());
    PeakIonPairPtrVec selected_pairs;
    getMatchedPairs(pairs,dp->getId(),selected_pairs);
    if(selected_pairs.size()>0){
      xercesc::DOMElement* mi_element = xml_doc->createElement("matched_ions");
      peak->appendChild(mi_element);
      for(unsigned int j=0;j< selected_pairs.size();j++){
        selected_pairs[j]->appendIonToXml(xml_doc,mi_element);
      }
    }
  }
  element->appendChild(ms_element);

  //proteoform to view
  xercesc::DOMElement* prot_element = geneProteinView(xml_doc,
                                                      prsm->getProteoformPtr(),
                                                      prsm->getRefineMs(),
                                                      prsm->getMinMass());
  element->appendChild(prot_element);

  return element;
}

xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,
                                     ProteoformPtr proteoform_ptr,
                                     ExtendMsPtr refine_ms_three,
                                     double min_mass){
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  std::string str=convertToString(proteoform_ptr->getDbResSeqPtr()->getId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=convertToString(proteoform_ptr->getSpeciesId());
  xml_doc->addElement(prot_element, "species_id", str.c_str());
  str=proteoform_ptr->getDbResSeqPtr()->getName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  double mass = proteoform_ptr->getMass();
  str=convertToString(mass);
  xml_doc->addElement(prot_element, "protein_mass", str.c_str());
  str=convertToString(proteoform_ptr->getStartPos());
  xml_doc->addElement(prot_element, "first_residue_position", str.c_str());
  str=convertToString(proteoform_ptr->getEndPos());
  xml_doc->addElement(prot_element, "last_residue_position", str.c_str());
  str=convertToString(proteoform_ptr->getProtModPtr()->getPtmPtr()->isAcetylation());
  xml_doc->addElement(prot_element, "protein_mass", str.c_str());
  int know_shift_number = proteoform_ptr->getChangePtrVec().size()-proteoform_ptr->getUnexpectedChangeNum();
  str=convertToString(know_shift_number);
  xml_doc->addElement(prot_element, "know_shift_number", str.c_str());
  for(unsigned int i=0;i<proteoform_ptr->getChangePtrVec().size();i++){
    proteoform_ptr->getChangePtrVec()[i]->appendXml(xml_doc,prot_element);//attention
  }
  xercesc::DOMElement* annotation_element = xml_doc->createElement("annotation");
  prot_element->appendChild(annotation_element);
  CleavagePtrVec cleavages = getProteoCleavage(proteoform_ptr,refine_ms_three,min_mass);
  for(int i=0;i<proteoform_ptr->getResSeqPtr()->getLen();i++){
    cleavages[i]->appendXml(xml_doc,annotation_element);
    proteoform_ptr->getResSeqPtr()->getResiduePtr(i)->appendXml(xml_doc,annotation_element);
  }
  cleavages[cleavages.size()-1]->appendXml(xml_doc,annotation_element);
  return prot_element;
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

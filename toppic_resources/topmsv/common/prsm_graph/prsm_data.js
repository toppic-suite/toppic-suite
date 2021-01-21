class PrsmData {
  residues; 
  formFirstPos;
  formLastPos;
  breakPoints;

  sequence;
  fixedPtms;
  variablePtms;
  massShifts;
  annotations;
  proteoform;

  rowNum;
  displayFirstPos;
  displayLastPos;
  // if there is a skipping line at the beginning
  showStartSkipped; 
  startSkippedInfo = "";
  // if there is a skipping line at the end
  showEndSkipped;
  endSkippedInfo = "";

  initData = function(prsm, para) {
    //console.log(prsm);
    this.extractData(prsm, para);
    this.updatePara(para);
    this.addColor();
  }

  extractData = function(prsm, para) {
    this.residues = prsm.annotated_protein.annotation.residue; 
    this.formFirstPos = parseInt(prsm.annotated_protein.annotation.first_residue_position);
    this.formLastPos = parseInt(prsm.annotated_protein.annotation.last_residue_position);
    this.breakPoints = json2BreakPoints(prsm);
    [this.fixedPtms, this.protVarPtms, this.variablePtms] = json2Ptms(prsm);
    //console.log(this.fixedPtms, this.protVarPtms, this.variablePtms)
    this.massShifts = json2MassShifts(prsm);
    this.annotations = getAnnotations(this.protVarPtms, prsm);
    this.sequence = this.getAminoAcidSequence();
    this.proteoform = new Proteoform(this.sequence, this.formFirstPos, this.fixedPtms, 
      this.protVarPtms, this.variablePtms, this.massShifts);
  }

  setDataFromUserInput = function(residues, formFirstPos, formLastPos, breakPoints, proteoformObj){
    //because Inspect page does not have prsm object, it cannot use the function above
    //but it can use proteoform object instead, which contains similar information as prsm 
    //which is generated based on user input
    this.residues = residues;
    this.formFirstPos = formFirstPos;
    this.formLastPos = formLastPos;
    this.fixedPtms = proteoformObj.fixedPtms;
    this.protVariablePtms = proteoformObj.protVarPtms;
    this.variablePtms = proteoformObj.variablePtms;
    this.massShifts = proteoformObj.unexpectedMassShifts;
    this.sequence = proteoformObj.sequence;
    this.proteoform = proteoformObj;
    this.breakPoints = breakPoints;

    function generateAnnotation(protVarPtm, variablePtm, massShift){
      let anno = [];
      //annotation in prsm object contains only variable and unknown shifts
      for (let i = 0; i < protVarPtm.length; i++){
        let temp = {"annoText":protVarPtm[i].getAnnotation(), "leftPos":0, "rightPos":0};
        temp.leftPos = protVarPtm[i].getLeftPos();
        temp.rightPos = temp.leftPos + 1;
        anno.push(temp);
      }
      for (let i = 0; i < variablePtm.length; i++){
        let temp = {"annoText":variablePtm[i].getAnnotation(), "leftPos":0, "rightPos":0};
        temp.leftPos = variablePtm[i].getLeftPos();
        temp.rightPos = temp.leftPos + 1;
        anno.push(temp);
      }
      for (let i = 0; i < massShift.length; i++){
        let temp = {"annoText":massShift[i].getAnnotation(), "leftPos":massShift[i].getLeftPos(), "rightPos":massShift[i].getRightPos()};
        anno.push(temp);
      }
      return anno;
    }
    this.annotations = generateAnnotation(this.protVariablePtms, this.variablePtms, this.massShifts);
  }
  
  updatePara = function(para) {
    let len = this.residues.length; 
    //console.log(this.formFirstPos, this.formLastPos, len, para.rowLength);
    // Include 5 amino acids before and after the form
    this.displayFirstPos = Math.floor((this.formFirstPos - 5) / para.rowLength) * para.rowLength;
    if (this.displayFirstPos < 0) {
      this.displayFirstPos = 0;
    }
    this.displayLastPos = Math.ceil((this.formLastPos + 6) / para.rowLength) * para.rowLength - 1;
    //console.log("display last pos ", this.displayLastPos);
    if (this.displayLastPos > (len -1)) {
      this.displayLastPos = len -1;
    }
    this.rowNum = Math.ceil((this.displayLastPos - this.displayFirstPos + 1)/para.rowLength);
    
    // skipping line
    this.showStartSkipped = false;
    this.showEndSkipped = false;
    if (para.showSkippedLines) {
      if (this.displayFirstPos !== 0) {
        this.showStartSkipped = true;
		    this.startSkippedInfo = "... "+ this.displayFirstPos 
          + " amino acid residues are skipped at the N-terminus ... ";
      }
      if (this.displayLastPos !== len - 1) {
        this.showEndSkipped = true;
		    this.endSkippedInfo =  "... "+ (len - 1 - this.displayLastPos) 
          + " amino acid residues are skipped at the C-terminus ... ";
      }
    }

    //console.log(this.displayFirstPos, this.displayLastPos, 
    //  this.rowNum, this.showStartSkipped, this.showEndSkipped);
  }

  addColor = function() {
    for (let i = 0; i < this.residues.length; i++) {
      let residue = this.residues[i];
      let pos = residue.position;
      if (pos < this.formFirstPos || pos > this.formLastPos) {
        residue.color = "grey";
      }
      else {
        residue.color = "black";
      }
    }
    for (let i = 0; i < this.fixedPtms.length; i++) {
      let ptm = this.fixedPtms[i];
      let pos = ptm.getLeftPos();
      this.residues[pos].color = "red";
    }
  }

  getAminoAcidSequence = function() {
    let sequence = "";
    for (let i = this.formFirstPos; i <= this.formLastPos; i++) {
      sequence = sequence + this.residues[i].acid;
    }
    //console.log("sequence", sequence);
    return sequence;
  }
}

/**
 * Get the cleavage positions from the prsm data
 * @param {object} prsm - json obeject with complete prsm data 
 */
function json2BreakPoints(prsm) {
  let breakPoints = [] ;
  let dataBps = prsm.annotated_protein.annotation.cleavage; 
  for (let i = 0; i < dataBps.length; i++) {
    let dataBp = dataBps[i];
    if (dataBp.exist_n_ion == 0 && dataBp.exist_c_ion == 0) {
      continue;
    }
    let bp = {}; 
    bp.position = parseInt(dataBp.position) ;
    bp.existNIon = (dataBp.exist_n_ion == 1);
    bp.existCIon = (dataBp.exist_c_ion == 1);
    bp.anno = "";
    bp.masses = [];
    if(dataBp.matched_peaks != null) {
      let dataMasses = [];
      if(dataBp.matched_peaks.matched_peak.length > 1) {
        dataMasses = dataBp.matched_peaks.matched_peak; 
      }
      else {
        dataMasses.push(dataBp.matched_peaks.matched_peak);
      }
      for (let j = 0; j < dataMasses.length; j++) {
        let dataMass = dataMasses[j];
        let mass = {};
        // Ion type
        mass.ionType = dataMass.ion_type ;
        // Ion Display position
        mass.ionDispPos = parseInt(dataMass.ion_display_position) ;
        // Ion Charge
        mass.charge = parseInt(dataMass.peak_charge);
        // ion_position
        // mass.ionPos = parseInt(dataMass.ion_position);
        bp.masses.push(mass);
        if (bp.anno != "") {
          bp.anno = bp.anno + " ";
        }
        bp.anno = bp.anno + mass.ionType + mass.ionDispPos + " " + mass.charge + "+";
      }
    }
    breakPoints.push(bp);
  }
  return breakPoints;
}

function getJsonList(item) {
  let valueList = [];
  if (Array.isArray(item)) {
    valueList = item;
  }
  else {
    valueList.push(item);
  }
  return valueList;
}

function getAnnotations(protVarPtms, prsm) {
  let annos = [];
  for (let i = 0; i < protVarPtms.length; i++) {
    let ptm = protVarPtms[i]; 
    let pos = ptm.getLeftPos(); 
    let anno = {};
    anno.annoText = ptm.getAnnotation();
    anno.leftPos = pos;
    anno.rightPos = pos + 1;
    annos.push(anno);
  }
 
  if(prsm.annotated_protein.annotation.hasOwnProperty('mass_shift')) {
    let dataMassShifts = getJsonList(prsm.annotated_protein.annotation.mass_shift); 
    for (let i = 0; i < dataMassShifts.length; i++) {
      let dataShift = dataMassShifts[i];
      //console.log(dataShift);
      if(dataShift.right_position != "0") {
        let anno = {};
        anno.annoText = dataShift.anno; 
        anno.leftPos = dataShift.left_position; 
        anno.rightPos = dataShift.right_position;
        annos.push(anno);
      }
    }
  }
  annos.sort(function(x,y){
    return x.leftPos - y.leftPos;
  });
  return annos;
}
/**
 * Get occurence of fixed ptm positions
 * @param {object} prsm - json obeject with complete prsm data 
 */
function json2Ptms(prsm){
  let fixedPtmList = [];
  let protVarPtmList = [];
  let varPtmList = [];
  if(!prsm.annotated_protein.annotation.hasOwnProperty("ptm") ) {
    return [fixedPtmList, protVarPtmList, varPtmList];
  }
  let dataPtmList = getJsonList(prsm.annotated_protein.annotation.ptm); 

  for (let i = 0; i < dataPtmList.length; i++) {
    let dataPtm = dataPtmList[i];
    if(dataPtm.ptm_type == "Fixed" || dataPtm.ptm_type == "Protein variable" 
      || dataPtm.ptm_type == "Variable") {
      if (dataPtm.hasOwnProperty("occurence")) {
        let occList = getJsonList(dataPtm.occurence);
        //console.log(occList);
        for (let j = 0; j < occList.length; j++) {
          let occurence = occList[j];

          let ptm = new Ptm(occurence.anno, dataPtm.ptm.mono_mass, dataPtm.ptm.abbreviation);
          let massShift = new MassShift(occurence.left_pos, occurence.right_pos, ptm.getShift(), dataPtm.ptm_type, ptm.getName(), ptm);
          
          if (dataPtm.ptm_type == "Fixed") {
            fixedPtmList.push(massShift);
          }
          else if (dataPtm.ptm_type == "Protein variable") {
            protVarPtmList.push(massShift);
          }
          else {
            varPtmList.push(massShift);
          }
        }
      }
    }
    //console.log(dataPtm);
    
  }
  //console.log("fixed", fixedPtmList);
  //console.log("prot var", protVarPtmList);
  //console.log("var", varPtmList);
  return [fixedPtmList, protVarPtmList, varPtmList];
}

/**
 * Get left and right positions of background color and mass shift value
 * @param {object} prsm - json obeject with complete prsm data 
 */
function json2MassShifts(prsm) {
	let massShifts = [];
	if(prsm.annotated_protein.annotation.hasOwnProperty('mass_shift')) {
    let dataMassShifts = getJsonList(prsm.annotated_protein.annotation.mass_shift); 
    for (let i = 0; i < dataMassShifts.length; i++) {
      let dataShift = dataMassShifts[i];
			if(dataShift.shift_type == "unexpected" && dataShift.right_position != "0") {
        let massShift = new MassShift(dataShift.left_position, dataShift.right_position, dataShift.shift, dataShift.shift_type, dataShift.shift);
				massShifts.push(massShift) ;
			}
      else if (dataShift.right_position == 0) {
        console.log("Mass shift right position is 0!", dataShift);
      }
    }
	}
  return massShifts;
  /*
  // add protein N-terminal modifications			
  if(prsm.annotated_protein.annotation.hasOwnProperty('ptm')) {
    let ptms = getJsonList(prsm.annotated_protein.annotation.ptm); 
    for (let i = 0; i < ptms.length; i++) {
      let ptm = ptms[i];
      if(ptm.ptm_type != "Fixed" && ptm.hasOwnProperty("occurence")) {
        let occList = getJsonList(ptm.occurence);
        for (let j = 0; j < occList.length; j++) {
          let shift = {};
          shift.anno = ptm.ptm.abbreviation;
          shift.leftPos = occList[j].left_pos;
          shift.rightPos = occList[j].right_pos;
          massShifts.push(shift); 
        }
      }
    }
  }
  let noDupMassShift = [];
  let duplicate = false;
  //remove duplicate mass shifts
  for (let a = 0; a < massShifts.length; a++){
    let massShiftA = massShifts[a];
    for (let b = 0; b < noDupMassShift.length; b++){
      let massShiftB = noDupMassShift[b];
      if (massShiftA.anno == massShiftB.anno){
        if (massShiftA.leftPos == massShiftB.leftPos){
          if (massShiftA.rightPos == massShiftB.rightPos){
            duplicate = true;
          }
        }
      }
    }
    if (!duplicate){
      noDupMassShift.push(massShiftA);
      duplicate = false;
    }
  }
	return noDupMassShift ;
  */
}

class PrsmData {
  residues; 
  formFirstPos;
  formLastPos;
  breakPoints;

  sequence;
  fixedPtms;
  massShifts;
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
    this.extractData(prsm, para);
    this.updatePara(para);
    this.addColor();
  }

  extractData = function(prsm, para) {
    this.residues = prsm.annotated_protein.annotation.residue; 
    this.formFirstPos = parseInt(prsm.annotated_protein.annotation.first_residue_position);
    this.formLastPos = parseInt(prsm.annotated_protein.annotation.last_residue_position);
    this.breakPoints = json2BreakPoints(prsm);
    this.fixedPtms = json2FixedPtms(prsm);
    this.massShifts = json2MassShifts(prsm);
    this.sequence = this.getAminoAcidSequence();
    this.proteoform = new Proteoform(this.sequence, this.fixedPtms, this.massShifts);
  }

  setData = (residues, formFirstPos, formLastPos, fixedPtms, massShifts, sequence, breakPoints) => {
    this.residues = residues;
    this.formFirstPos = formFirstPos;
    this.formLastPos = formLastPos;
    this.fixedPtms = fixedPtms;
    this.massShifts = massShifts;
    this.sequence = sequence;
    this.proteoform = new Proteoform(this.sequence, this.fixedPtms, this.massShifts);
    this.breakPoints = breakPoints;
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
      for (let j = 0; j < ptm.posList.length; j++) {
        let pos = ptm.posList[j].pos;
        this.residues[pos].color = "red";
      }
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

/**
 * Get occurence of fixed ptm positions
 * @param {object} prsm - json obeject with complete prsm data 
 */
function json2FixedPtms(prsm){
  let ptmList = [] ;
  if(!prsm.annotated_protein.annotation.hasOwnProperty("ptm") ) {
    return ptmList;
  }
  let dataPtmList = getJsonList(prsm.annotated_protein.annotation.ptm); 
  for (let i = 0; i < dataPtmList.length; i++) {
    let dataPtm = dataPtmList[i];
    let ptm = {};
    ptm.name = dataPtm.ptm.abbreviation;
    ptm.posList = [];
    if(dataPtm.ptm_type == "Fixed" && dataPtm.hasOwnProperty("occurence")) {
      let occList = getJsonList(dataPtm.occurence);
      for (let j = 0; j < occList.length; j++) {
        let occurence = occList[j];
        let pos = {};
        pos.pos = occurence.left_pos;
        pos.acid = occurence.anno;
        ptm.posList.push(pos);
      }
    }
    ptmList.push(ptm);
  }
  //console.log(ptmList);
  return ptmList;
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
			if(dataShift.right_position != "0") {
        let shift = {};
        shift.anno = dataShift.anno; 
        shift.leftPos = dataShift.left_position; 
        shift.rightPos = dataShift.right_position;
				massShifts.push(shift) ;
			}
      else {
        console.log("Mass shift right position is 0!", mass_shift);
      }
    }
	}
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
}


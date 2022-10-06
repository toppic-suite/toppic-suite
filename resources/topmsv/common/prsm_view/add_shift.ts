//enable click-to-add a mass shift in a prsm sequence in the inspect page. 
class AddShift {
  static appliedPtm: number[] = [];//variable PTM applied to the sequence so far
  static unknownMassShift: {"pos": number, "mass": number}[] = [];
  static clickedLetter: string = "";
  static clickedPos: number = -1;

  constructor() {
    this.addTabEventListener();
  };

  updateInspectPage(sequence: string, unknownMassShiftList: MassShift[], protVarPtmsList: MassShift[], variablePtmsList: MassShift[]): void {
    let errorType: string | number | string[] | undefined = $(".error_dropdown").val();
    if (errorType) {
      let errorVal: string | number | string[] | undefined = $("#errorval").val();
      if (typeof(errorVal) == "string") {
        errorVal = parseFloat(errorVal.trim());

        setDataToSequence(sequence, unknownMassShiftList, protVarPtmsList, variablePtmsList);

        let executionObj: SeqOfExecution = new SeqOfExecution();
        if (typeof(errorType) == "string") {
          executionObj.sequenceOfExecution(errorType, errorVal, "");
        }
      }
    }
  }

  addUnknownShift(e: JQuery.ClickEvent<HTMLElement, undefined, HTMLElement, HTMLElement>): void {
    /*let ptmDiv: HTMLElement | null = document.getElementById("tooltip-pop");
    if (ptmDiv) {
      ptmDiv.style.display = "none";
    }*/
    $("#tooltip-pop").modal('hide');
    let val: string | number | string[] | undefined = $("#mass-value").val();

    if (typeof(val) != "string") {
      alert("Please enter a valid number for mass shift!");
      return;
    }
    let mass: number = parseFloat(val);

    if (isNaN(mass)) {
      alert("Please enter a valid number for mass shift!");
      return;
    }

    let unknownMassShiftList: MassShift[] = [];
    let protVarPtmsList: MassShift[] = [];
    let variablePtmsList: MassShift[] = [];

    let sequence: string = "";
    sequence = getSequenceFromUI();

    let parseResult: [string, MassShift[], MassShift[], MassShift[]] | null = parseSequenceMassShift(sequence);

    if (!parseResult) {
      return;
    }
    sequence = parseResult[0];
    unknownMassShiftList = parseResult[1];
    protVarPtmsList = parseResult[2];
    variablePtmsList = parseResult[3];
    
    //remove previous mass shift if exists on a same residue
    for (let i: number = 0; i < unknownMassShiftList.length; i++) {
      if (unknownMassShiftList[i].getLeftPos() == AddShift.clickedPos) {
        unknownMassShiftList.splice(i, 1);
        break;
      }
    }
    for (let j: number = 0; j < protVarPtmsList.length; j++) {
      if (protVarPtmsList[j].getLeftPos() == AddShift.clickedPos) {
        for (let a = 0; a < AddShift.appliedPtm.length; a++) {
          if (commonPtmList[AddShift.appliedPtm[a]].abbr  == protVarPtmsList[j].getAnnotation()) {
            AddShift.appliedPtm.splice(a, 1);
            break;
          }
        }
        protVarPtmsList.splice(j, 1);
        break;
      }
    }
    for (let k: number = 0; k < variablePtmsList.length; k++) {
      if (variablePtmsList[k].getLeftPos() == AddShift.clickedPos) {
        for (let a = 0; a < AddShift.appliedPtm.length; a++) {
          if (commonPtmList[AddShift.appliedPtm[a]].abbr == variablePtmsList[k].getAnnotation()) {
            AddShift.appliedPtm.splice(a, 1);
            break;
          }
        }
        variablePtmsList.splice(k, 1);
        break;
      }
    }
    //if None was selected, don't add it to ptm 
    if (mass != 0.0) {
      //add new variable ptm
      unknownMassShiftList.push(new MassShift(AddShift.clickedPos, AddShift.clickedPos+1, mass, "Unknown", mass.toString()));
    }
    this.updateInspectPage(sequence, unknownMassShiftList, protVarPtmsList, variablePtmsList);
  }

  addTabEventListener(): void {
    //switch between variable ptm and unknown shift tabs
    $("#var-ptm-link").on("click", function() {
      $("#var-ptm-tab").show();
      $("#unknown-mod-tab").hide();
    })
    $("#unknown-mod-link").on("click", function() {
      $("#var-ptm-tab").hide();
      $("#unknown-mod-tab").show();
    })
    $("#add-unknown-mod-btn").on("click", (e) => {
      this.addUnknownShift(e);
    })
    //event listener for search bar
  }

  applyShift(ptmIdx: number, letter: string, pos: number): void {
    /*let ptmDiv: HTMLElement | null = document.getElementById("tooltip-pop");
    if (ptmDiv) {
      ptmDiv.style.display = "none";
    }*/
    $("#tooltip-pop").modal('hide');
    let mass: number = commonPtmList[ptmIdx].mass;
    let unknownMassShiftList: MassShift[] = [];
    let protVarPtmsList: MassShift[] = [];
    let variablePtmsList: MassShift[] = [];


    let sequence: string = "";
    sequence = getSequenceFromUI();

    let parseResult: [string, MassShift[], MassShift[], MassShift[]] | null = parseSequenceMassShift(sequence);

    if (!parseResult) {
      return;
    }
    
    sequence = parseResult[0];
    unknownMassShiftList = parseResult[1];
    protVarPtmsList = parseResult[2];
    variablePtmsList = parseResult[3];
    
    //remove previous mass shift if exists on a same residue
    for (let i: number = 0; i < unknownMassShiftList.length; i++) {
      if (unknownMassShiftList[i].getLeftPos() == pos) {
        unknownMassShiftList.splice(i, 1);
        break;
      }
    }
    for (let j: number = 0; j < protVarPtmsList.length; j++) {
      if (protVarPtmsList[j].getLeftPos() == pos) {
        for (let a = 0; a < AddShift.appliedPtm.length; a++) {
          if (commonPtmList[AddShift.appliedPtm[a]].abbr  == protVarPtmsList[j].getAnnotation()) {
            AddShift.appliedPtm.splice(a, 1);
            break;
          }
        }
        protVarPtmsList.splice(j, 1);
        break;
      }
    }
    for (let k: number = 0; k < variablePtmsList.length; k++) {
      if (variablePtmsList[k].getLeftPos() == pos) {
        for (let a = 0; a < AddShift.appliedPtm.length; a++) {
          if (commonPtmList[AddShift.appliedPtm[a]].abbr == variablePtmsList[k].getAnnotation()) {
            AddShift.appliedPtm.splice(a, 1);
            break;
          }
        }
        variablePtmsList.splice(k, 1);
        break;
      }
    }
    //if None was selected, don't add it to ptm 
    if (mass != 0) {
      //add new variable ptm
      variablePtmsList.push(new MassShift(pos, pos+1, mass, "Variable", commonPtmList[ptmIdx].abbr, new Mod(letter, mass, commonPtmList[ptmIdx].name)));
      
      if (!AddShift.appliedPtm.includes(ptmIdx)) {
        AddShift.appliedPtm.push(ptmIdx);
      }
    }
    this.updateInspectPage(sequence, unknownMassShiftList, protVarPtmsList, variablePtmsList);
  }
  static searchBar(): void {
    let input: HTMLInputElement = <HTMLInputElement>document.getElementById("search-input");
    let filter: string = input.value.toUpperCase();
    let ul: HTMLUListElement = <HTMLUListElement>document.getElementById("ptm-list");
    let li: HTMLCollectionOf<HTMLLIElement> = ul.getElementsByTagName("li");
    for (let i: number = 0; i < li.length; i++) {
      let ptm: ChildNode | null = li[i].firstChild;
      if (!ptm) {
        break;
      }
      let ptmName: string = <string>ptm.nodeValue;
      if (ptmName.toUpperCase().indexOf(filter) > -1) {
        li[i].style.display = "";
      }
      else {
        li[i].style.display = "none";
      }
    }
  }

  handleOnClick(letter: string, pos: number): void{
    AddShift.clickedLetter = letter;
    AddShift.clickedPos = pos;

    let ptmDiv: HTMLElement | null = document.getElementById("tooltip-pop");
    if (ptmDiv != null) {
      $("#tooltip-pop").modal('show');
      ptmDiv.style.opacity = "1";

      $('#ptm-list').empty();
      $('#applied-ptm-list').empty();

      commonPtmList.forEach((ptm, idx) => {
        let entry: JQuery<HTMLLIElement> = $("<li></li>");
        entry.text(ptm.abbr + " (" + ptm.name + "); " + ptm.mass);
        entry.attr("id", "ptm" + idx);
        entry.attr("class", "text-center");
        entry.css("cursor", "pointer");
        entry.on("click", () => {this.applyShift(idx, letter, pos)});
        $("#ptm-list").append(entry);
      })

      AddShift.appliedPtm.forEach((ptmIdx, idx) => {
        let entry: JQuery<HTMLLIElement> = $("<li></li>");
        entry.text(commonPtmList[ptmIdx].abbr + " (" + commonPtmList[ptmIdx].name + "); " + commonPtmList[ptmIdx].mass);
        entry.attr("id", "applied-ptm" + commonPtmList[ptmIdx].abbr);
        entry.attr("class", "text-center");
        entry.css("cursor", "pointer");
        entry.on("click", () => {this.applyShift(ptmIdx, letter, pos)});
        $("#applied-ptm-list").append(entry);
      })
    }
  }
}
/**
 * Sets all the fixed ptms to UI under Fixed Ptm column
 */
function setFixedPtmListToUI(commonFixedPtmList: Mod[]): void{
    let dropdown: HTMLElement = domElements.dropDownMenuLink;
    if (!dropdown){
        console.error("ERROR: dropdown element is empty");
        return;
    }
    commonFixedPtmList.forEach((fixedPtm) => {
        let value: string = fixedPtm.getResidue()+" : "+ fixedPtm.getShift().toString();
        let option: HTMLOptionElement = document.createElement("option");
        let dropdown = domElements.dropDownMenuLink;
        
        option.setAttribute("value",value);
        option.innerHTML = value;
        dropdown!.appendChild(option);//already checked for null
    })
    jqueryElements.addFixedPtmRow.click(() => {
        let fixedPtm: string = (<HTMLInputElement>dropdown!).value;//already checked for null
        if(fixedPtm !== "fixed_ptm")
        {
            addNewFixedPtmRow(fixedPtm);
        }
    })
}

/**
 * @function addNewFixedPtmRow
 * @description 
 * On Click of '+' symbol under fixed ptms in the HTML, creates new block to add Acid and mass shift
 * @param {String} fixedPtm - Contains Acid name and mass shift seperated by ':'
 */
function addNewFixedPtmRow(fixedPtm: string){
    let acid: string = '';
    let mass: string = '';
    let ifExist: boolean = false;
    if(fixedPtm !== "other")
    {
        let splitVal: string[] = fixedPtm.split(":");
        acid = splitVal[0].trim();
        mass = splitVal[1].trim();
        let existingPtmList: Mod[] = getFixedPtmCheckList();
        existingPtmList.forEach((element) => {
            //  && element.mass === parseFloat(mass)
            if (element.getResidue() === acid) {
                if (element.getShift().toString() == mass) {
                    ifExist = true;
                    return;
                }
                let replace: boolean = confirm("Fixed ptm is already applied for acid " + acid + ". Do you want to replace it?")         
                if (replace){
                   $(".fixedptms").each((i, div) => {
                       console.log(div)
                    //let acid: string = div.val().toUpperCase();
                    if (acid == element.getResidue()){
                        div.remove();
                    }
                });
                }else{
                    ifExist = true;
                }    
            }
        });
    }
    if (ifExist) return;
    let fixedPtmListDiv: HTMLElement = <HTMLElement>domElements.fixedPtmList;
    let fixedptmsdiv = document.createElement("div");
    fixedptmsdiv.setAttribute("class","fixedptms");
    
    //Creating div with input fixed acid and fixed mass 
    let inputAcid: HTMLInputElement = document.createElement("input");
    inputAcid.setAttribute("type","text");
    inputAcid.setAttribute("class","form-control fixedptmacid");
    // inputAcid.setAttribute("id","fixedptmacid");
    inputAcid.setAttribute("value",acid);
    
    let span: HTMLSpanElement = document.createElement("span");
    span.innerHTML = "&nbsp;:&nbsp;";
    
    let inputMass: HTMLInputElement = document.createElement("input");
    inputMass.setAttribute("type","text");
    inputMass.setAttribute("class","form-control fixedptmmass");
    // inputMass.setAttribute("id","fixedptmmass");
    inputMass.setAttribute("value",mass);
    
    let span2: HTMLSpanElement = document.createElement("span");
    span2.innerHTML = "&nbsp;";
    
    let addButton: HTMLButtonElement = document.createElement("button");
    addButton.setAttribute("type","button");
    addButton.setAttribute("class","form-control btn btn-default btn-sm addnewrow");
    
    let iAddFrame: HTMLElement = document.createElement("i");
    iAddFrame.setAttribute("class","fa fa-plus");
    
    addButton.appendChild(iAddFrame);
    
    let span3: HTMLSpanElement = document.createElement("span");
    span3.innerHTML = "&nbsp;";
    
    let removeButton: HTMLButtonElement = document.createElement("button");
    removeButton.setAttribute("type","button");
    removeButton.setAttribute("class","form-control btn btn-default btn-sm removerow");

    let iRemoveFrame: HTMLElement = document.createElement("i");
    iRemoveFrame.setAttribute("class","fa fa-times");
    removeButton.appendChild(iRemoveFrame);
    
    fixedptmsdiv.appendChild(inputAcid);
    fixedptmsdiv.appendChild(span);
    fixedptmsdiv.appendChild(inputMass);
    fixedptmsdiv.appendChild(span2);
    fixedptmsdiv.appendChild(removeButton);
    fixedPtmListDiv.appendChild(fixedptmsdiv);	
    
    $('.removerow').click((e) => {
        let parents: JQuery<HTMLElement> = $(e.target).parents();
        for (let i: number = 0; i < parents.length; i++) {
            if (parents[i].className == "fixedptms") {
                let residue: HTMLInputElement = <HTMLInputElement>(parents[i].getElementsByClassName("fixedptmacid")[0]);
                let mass: HTMLInputElement = <HTMLInputElement>(parents[i].getElementsByClassName("fixedptmmass")[0]);
                for (let j = 0; j < USER_FIXED_PTM_LIST.length; j++) {
                    if (parseFloat(mass.value) == USER_FIXED_PTM_LIST[j].getShift()) {
                        if (residue.value == USER_FIXED_PTM_LIST[j].getResidue()) {
                            USER_FIXED_PTM_LIST.splice(j, 1);
                            break;
                        }
                    }
                }
                parents[i].remove();
            }
        }
        //temp code
        let errorVal: number;
        let errorType: any = jqueryElements.errorDropdown.val();
        if (typeof(errorType) == "string") {
            if(errorType === "masserror") {
                let val = jqueryElements.errorValue.val();
                if (typeof(val) == "string") {
                    errorVal = parseFloat(val.trim())
                }
            }
            else {
                let val = jqueryElements.errorValue.val();
                if (typeof(val) == "string") {
                    errorVal = parseFloat(val.trim());
                }
            }
        }
        // here
        // reload seqOfExecution to refresh result
        // let executionObj = new SeqOfExecution();
        // executionObj.sequenceOfExecution(errorType,errorVal,acid);
    })
}

/**
 * Set all the fixed masses on html
 */
function setFixedMasses(fixedPtmList: MassShift[]){
    if(fixedPtmList.length !=0)
    {
        // let commonfixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];
        for(let i=0;i<fixedPtmList.length;i++)
        {
            for(let j=0; j<COMMON_FIXED_PTM_LIST.length;j++)
            {
                if(fixedPtmList[i].getAnnotation().toUpperCase() ===  COMMON_FIXED_PTM_LIST[j].getName().toUpperCase()){
                    let existingPtmList: Mod[] = getFixedPtmCheckList();
                    let isNewPtm: boolean = true;
                    if (existingPtmList.length > 0){//prevent same ptm added again
                        existingPtmList.forEach(ptm => {
                            if (ptm.getShift().toFixed(4) == fixedPtmList[i].getShift().toString()){
                                isNewPtm = false;
                            }  
                        })
                    }
                    if (isNewPtm){
                        let fixedptm: string = COMMON_FIXED_PTM_LIST[j].getResidue() + ":" + COMMON_FIXED_PTM_LIST[j].getShift().toString();
                        addNewFixedPtmRow(fixedptm);
                        break;    
                    }
                }
            }
        }
    } 
}

/**
 * @return {Array} FixedPtmList - return all the selected fixed ptms with acid and mass
 */
function getFixedPtmCheckList(): Mod[]
{
    let result: Mod[] = [];
    $(".fixedptms").each((i, div) => {
        let acid: string | null = null;
        let mass: string | null = null;
        for (let i: number = 0; i < div.children.length; i++) {
            if (div.children[i].className == "form-control fixedptmacid") {
                let acidElement: HTMLInputElement = <HTMLInputElement> div.children[i];
                acid = acidElement.value;
            }
            else if (div.children[i].className == "form-control fixedptmmass") {
                let massElement: HTMLInputElement = <HTMLInputElement> div.children[i];
                mass = massElement.value;
            }
        }
        if (acid && mass) {
            let tempfixedptm: Mod = new Mod(acid.toUpperCase(), parseFloat(mass), "");
            result.push(tempfixedptm);
        }
    });
    return result;
}

/**
 * sets fixed PTM based on user input
 */
function addCustomPtm(): void {
    let acid: HTMLInputElement = <HTMLInputElement>document.getElementById("residue-name");
    let mass: HTMLInputElement = <HTMLInputElement>document.getElementById("residue-mass");
    let name: HTMLInputElement = <HTMLInputElement>document.getElementById("ptm-name");

    if (!acid || !mass) {
        alert("Please enter both acid and mass for this PTM!");
    }
    else {
        if (/\d/.test(acid.value)) {
            alert("Invalid amino acid!");
        }
        else if(isNaN(parseFloat(mass.value))) {
            alert("Invalid mass value!");
        }
        else {
            let fixedPtm = (acid.value).toUpperCase() + ":" + mass.value;
            USER_FIXED_PTM_LIST.push(new Mod((acid.value).toUpperCase(), parseFloat(mass.value), name.value))
            addNewFixedPtmRow(fixedPtm);
            $('.modal').modal('hide');
        }
    }
}
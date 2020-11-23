/**
 * Sets all the fixed ptms to UI under Fixed Ptm column
 */
function setFixedPtmListToUI(commonFixedPtmList){    
    // commonFixedPtmList EX. [{acid: XX, mass: XX}]
    commonFixedPtmList.forEach((fixedPtm) => {
        let value = fixedPtm.acid+" : "+fixedPtm.mass;
        let option = document.createElement("option");
        option.setAttribute("value",value);
        option.innerHTML = value;
        domElements.dropDownMenuLink.appendChild(option);
    })
    jqueryElements.addFixedPtmRow.click(() => {
        let fixedPtm = domElements.dropDownMenuLink.value;
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
function addNewFixedPtmRow(fixedPtm){
    // console.log("fixedptm : ", fixedptm);
    let acid = '';
    let mass = '';
    let ifExist = false;
    if(fixedPtm !== "other")
    {
        let splitVal = fixedPtm.split(":");
        acid = splitVal[0].trim();
        mass = splitVal[1].trim();
        let existingPtmList = getFixedPtmCheckList();
        existingPtmList.forEach((element) => {
            //  && element.mass === parseFloat(mass)
            if (element.acid === acid) {
                let replace = confirm("Fixed ptm is already applied for acid " + acid + ". Do you want to replace it?")         
                if (replace){
                   //let acid = element.acid;
                   // $(this).parent().remove();
                   $(".fixedPtms").each(function(){
                    let acid = $( this ).find('.fixedptmacid').val().toUpperCase();
                    if (acid == element.acid){
                        $( this ).remove();
                    }
                });
                }else{
                    ifExist = true;
                }    
            }
        });
    }
    if (ifExist) return;
    let fixedPtmListDiv = domElements.fixedPtmList;
    let fixedptmsdiv = document.createElement("div");
    fixedptmsdiv.setAttribute("class","fixedptms");
    
    //Creating div with input fixed acid and fixed mass 
    let inputAcid = document.createElement("input");
    inputAcid.setAttribute("type","text");
    inputAcid.setAttribute("class","form-control fixedptmacid");
    // inputAcid.setAttribute("id","fixedptmacid");
    inputAcid.setAttribute("value",acid);
    
    let span = document.createElement("span");
    span.innerHTML = "&nbsp;:&nbsp;";
    
    let inputMass = document.createElement("input");
    inputMass.setAttribute("type","text");
    inputMass.setAttribute("class","form-control fixedptmmass");
    // inputMass.setAttribute("id","fixedptmmass");
    inputMass.setAttribute("value",mass);
    
    let span2 = document.createElement("span");
    span2.innerHTML = "&nbsp;";
    
    let addButton = document.createElement("button");
    addButton.setAttribute("type","button");
    addButton.setAttribute("class","form-control btn btn-default btn-sm addnewrow");
    
    let iAddFrame = document.createElement("i");
    iAddFrame.setAttribute("class","fa fa-plus");
    addButton.appendChild(iAddFrame);
    
    let span3 = document.createElement("span");
    span3.innerHTML = "&nbsp;";
    
    let removeButton = document.createElement("button");
    removeButton.setAttribute("type","button");
    removeButton.setAttribute("class","form-control btn btn-default btn-sm removerow");
    
    let iRemoveFrame = document.createElement("i");
    iRemoveFrame.setAttribute("class","fa fa-times");
    removeButton.appendChild(iRemoveFrame);
    
    fixedptmsdiv.appendChild(inputAcid);
    fixedptmsdiv.appendChild(span);
    fixedptmsdiv.appendChild(inputMass);
    fixedptmsdiv.appendChild(span2);
    fixedptmsdiv.appendChild(removeButton);
    fixedPtmListDiv.appendChild(fixedptmsdiv);	
    
    $('.removerow').click(function(){
        let acid = $(this).parent().find(".fixedptmacid").val();
        $(this).parent().remove();
        //temp code
        let errorVal ;
        let errorType = jqueryElements.errorDropdown.val();
        if(errorType === "masserror") {
            errorVal = parseFloat(jqueryElements.errorValue.val().trim());
        }
        else {
            errorVal = parseFloat(jqueryElements.errorValue.val().trim());
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
function setFixedMasses(fixedPtmList){
    if(fixedPtmList.length !=0)
    {
        // let commonfixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];
        for(let i=0;i<fixedPtmList.length;i++)
        {
            for(let j=0; j<COMMON_FIXED_PTM_LIST.length;j++)
            {
                if(fixedPtmList[i].name.toUpperCase() ===  COMMON_FIXED_PTM_LIST[j].name.toUpperCase())
                {
                    let fixedptm = COMMON_FIXED_PTM_LIST[j].acid + ":" + COMMON_FIXED_PTM_LIST[j].mass;
                    addNewFixedPtmRow(fixedptm);
                    break;
                }
            }
        }
    } 
}

/**
 * @return {Array} FixedPtmList - return all the selected fixed ptms with acid and mass
 */
function getFixedPtmCheckList()
{
    let result = [];
    $(".fixedptms").each(function(){
        let acid = $( this ).find('.fixedptmacid').val().toUpperCase();
        let mass = parseFloat($( this ).find('.fixedptmmass').val());
        if(acid.length !== 0  && !isNaN(mass))
        {
            let tempfixedptm = {acid:acid,mass:mass}
            result.push(tempfixedptm);
        }
    });
    return result;
}
/**
 * Sets all the fixed ptms to UI under Fixed Ptm column
 */
function setFixedPtmListToUI(commonFixedPtmList){
    let dropDownMenuLink = domElements.dropDownMenuLink;
    
    commonFixedPtmList.forEach((fixedPtm) => {
        let value = fixedPtm.acid+" : "+fixedPtm.mass;
        let option = document.createElement("option");
        option.setAttribute("value",value);
        option.innerHTML = value;
        dropDownMenuLink.appendChild(option);
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
    if(fixedPtm !== "other")
    {
        let splitVal = fixedPtm.split(":");
        acid = splitVal[0].trim();
        mass = splitVal[1].trim();
    }
    
    let fixedPtmListDiv = domElements.fixedPtmList;
    let fixedptmsdiv = document.createElement("div");
    fixedptmsdiv.setAttribute("class","fixedptms");
    
    //Creating div with input fixed acid and fixed mass 
    let inputAcid = document.createElement("input");
    inputAcid.setAttribute("type","text");
    inputAcid.setAttribute("class","form-control");
    inputAcid.setAttribute("id","fixedptmacid");
    inputAcid.setAttribute("value",acid);
    
    let span = document.createElement("span");
    span.innerHTML = "&nbsp;:&nbsp;";
    
    let inputMass = document.createElement("input");
    inputMass.setAttribute("type","text");
    inputMass.setAttribute("class","form-control");
    inputMass.setAttribute("id","fixedptmmass");
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
    
    jqueryElements.removeFixedPtmRow.click(() => {
        let acid = $(this).parent().find("#fixedptmacid").val();
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
        let executionObj = new SeqOfExecution();
        executionObj.sequenceOfExecution(errorType,errorVal,acid);
    })
}

/**
 * Set all the fixed masses on html
 */
function setFixedMasses(fixedPtmList){
    if(fixedPtmList.length !=0)
    {
        let commonfixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];
        for(let i=0;i<fixedPtmList.length;i++)
        {
            for(let j=0; j<commonfixedPtmList.length;j++)
            {
                if(fixedPtmList[i].name.toUpperCase() ===  commonfixedPtmList[j].name.toUpperCase())
                {
                    let fixedptm = commonfixedPtmList[j].acid + ":" + commonfixedPtmList[j].mass;
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
getFixedPtmCheckList()
{
    let FixedPtmList = [];
    let divs = $( ".fixedptms").get();
    $( ".fixedptms" ).each(function( index ) {
        let acid = $( this ).find('#fixedptmacid').val().toUpperCase();
        let mass = parseFloat($( this ).find('#fixedptmmass').val());
        if(acid.length !=0  && mass.length != 0 && !isNaN(mass))
        {
            let tempfixedptm = {acid:acid,mass:mass}
            FixedPtmList.push(tempfixedptm);
        }
    });
    return FixedPtmList;
}
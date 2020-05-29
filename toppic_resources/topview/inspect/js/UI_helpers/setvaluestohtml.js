/**
 * Class to set all the data retrieved from local storage on to html
 */
class SetValuesToHtml{
    constructor(peakAndIntensityList,massAndIntensityList,sequence,l_fixedPtmList,unknownMassShiftList,precursorMass){
        this.peakAndIntensityList = peakAndIntensityList;
        this.massAndIntensityList = massAndIntensityList;
        this.sequence = sequence;
        this.fixedPtmList = l_fixedPtmList;
        this.unknownMassShiftList = unknownMassShiftList ;
        this.precursorMass = precursorMass ;
    }
    /**
     * Set Peak List on to html
     */
    setDataToPeakAndIntensity(){
        let peakdata = $("#peakdata");
        peakdata.val(this.peakAndIntensityList);
    }
    /**
     * Set Mass List on to html
     */
    setDataToMassAndIntensity(){
        let massdata = $("#data");
        massdata.val(this.massAndIntensityList);
    }
    /**
     * Set Sequence on to html
     */
    setDataToSequence(){
        let MassShiftsObj = new MassShifts();
        let modSequence = MassShiftsObj.formSequence(this.sequence,this.unknownMassShiftList);
        $("#sequencedata").val(modSequence);
    }
    /**
     * Set all the fixed masses on html
     */
    setFixedMasses(){
        if(this.fixedPtmList.length !=0)
        {
            let commonFixedPtmsObj = new commonFixedPtms();
            let comonfixedPtmList = commonFixedPtmsObj.fixedPtmList; 
            for(let i=0;i<this.fixedPtmList.length;i++)
            {
                for(let j=0; j<comonfixedPtmList.length;j++)
                {
                    if(this.fixedPtmList[i].name.toUpperCase() ==  comonfixedPtmList[j].name.toUpperCase())
                    {
                        let fixedptm = comonfixedPtmList[j].acid + ":" + comonfixedPtmList[j].mass;
                        commonFixedPtms.addNewFixedPtmRow(fixedptm);
                        break;
                    }
                }
            }
        } 
    }
    /**
     * Set Precursor mass on html 
     */
    setPrecursorMass()
    {
        document.getElementById("precursormass").innerHTML = this.precursorMass;
    }
    
}
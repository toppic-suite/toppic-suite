"use strict";
//let backGroundColorList_g = [];
class SeqOfExecution {
    /*constructor(){
      this.onClickMassShift = {};
      this.massShiftList = [];
    }*/
    /**
    * Function executes all the functionalities one by one and displays all the
    * needed contant on to HTML.
    * @param {string} errorType - This gives which type of error needed to be
    * considered when matched peaks are to be considered.
    * @param {float} errorVal - This gives the user entered threshhold value to
    * be considered when calculating matched peaks.
    * @param {char} removeAcid - This gives the acid up on which the fixed ptm
    * mass has to be removed when "X" is clicked at fixed ptms.
    */
    sequenceOfExecution(errorType: string, errorVal: number, removeAcid: string) {
        /**
        * unbind all the actions previously binded else each action will be
        * binded multiple times.
        */
        let fixedMassShiftList: MassShift[] = []; //contains fixedPTM
        let protVarPtmsList: MassShift[] = []; //contains protein variable PTM
        let variablePtmsList: MassShift[] = []; //contains non-protein variable PTM
        let completeShiftList: MassShift[] = []; //contains all 3 kinds of mass shifts
        let unknownMassShiftList: MassShift[] = []; //contains unknown mass shifts
        let modifiablePeakData: Peak[] = []; //will change value if shared peak
        let massErrorthVal: number;
        let matchedPeakList: MatchedUnMatchedPeak[] = [];
        let ppmErrorthVal: number;
        let spectrumGraphObj: SpectrumView = {} as SpectrumView;
        let monoMassGraphObj: SpectrumView = {} as SpectrumView;
        let precursorMass: number | null = getPrecursorMass();
        let ms2GraphList: SpectrumView[] = [];
        let prsmObj: Prsm = {} as Prsm;
        if (!precursorMass) {
            console.error("Error: precursor mass is invalid");
            return;
        }
        /* show submit button for precursor mass and add event handler*/
        jqueryElements.precursorMassSubmit.show();
        setPrecursorMassEventHandler();
        /* Hide everything when page launched before data is computed*/
        $("#" + Constants.SEQSVGID).hide();
        $("#" + Constants.SVGDOWNLOADID).hide();
        $("#" + Constants.SPECTRUMGRAPHID).hide();
        $("#" + Constants.GRAPHDOWNLOAD).hide();
        $("#" + Constants.DIVTABLECONTAINER).hide();
        $("#" + Constants.PEAKCOUNTID).hide();
        $("#" + Constants.MONOMASSGRAPHID).hide();
        /**
        * Get the parsed sequence after removing mass shift list from
        * the entered sequence.
        * Returns mass list embedded in [] in sequence of user entered sequence.
        */
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

        fixedMassShiftList = parseCheckedFixedPtm(sequence);
        /* If user removed fixed ptm mass, remove the mass from the list*/
        if (removeAcid !== "") {
            removeAcid = removeAcid.toUpperCase();
            for (let i: number = 0; i < fixedMassShiftList.length; i++) {
                let position: number = fixedMassShiftList[i].getLeftPos();
                if (sequence[position] === removeAcid) {
                    fixedMassShiftList.splice(i, 1);
                }
            }
        }
        /*create proteoform object containing mass shift information and sequence */
        let proteoformObj: Proteoform = new Proteoform("", "", "", sequence, "0", 0, sequence.length - 1, -1, unknownMassShiftList, fixedMassShiftList, protVarPtmsList, variablePtmsList);
        /**
        * Check if mass shift is entered by clicking on the acid. If entered
        * consider that mass shift and and append to the current mass shift list
        */
        /*if(!$.isEmptyObject(this.onClickMassShift)) {
          let tempPosition = this.onClickMassShift.position;
          let tempMass = this.onClickMassShift.mass;
          massShiftObj.appendtoMassShiftList(tempPosition,tempMass);
        }*/
        /**
        * Form sequence with mass shift embedded in []
        */
        // let seqToUI = massShiftObj.formSequence();
        /**
        * Write back to UI
        */
        // writeSeqToTextBox(seqToUI);
        /* Get all the Mass List data entered by the user.*/
        let monoMassList: Peak[] = getMassListFromUI();
        /* Get all the peak list data entered by the user */
        modifiablePeakData = getPeakListFromUI();
        let monoMassListLen: number = monoMassList.length;
        let seq: string = proteoformObj.getSeq();
        let matchedUnMatchedPeaks: MatchedUnMatchedPeak[] = [];
        let envelopeList: Envelope[];
        /* Setting masserror threshold value and ppm error threshhold value*/
        if (errorType === Constants.MASSERROR)
            massErrorthVal = errorVal;
        else
            ppmErrorthVal = errorVal;
        /*Get all the n, c terminus ions selected.*/
        let n_TerminusList: Ion[] = getNterminusCheckedList();
        let c_TerminusList: Ion[] = getCterminusCheckedList();
        /* create spectrum object*/
        let spectrum: Spectrum = new Spectrum("", "", 2, getPeakListFromUI(), [], [], n_TerminusList, c_TerminusList, precursorMass);
        /* Get all the matched peaks for all the n terminus fragmented ions selected.*/
        let calcMatchedPeaks: CalcMatchedPeaks = new CalcMatchedPeaks();
        spectrum.getNTerminalIon().forEach(function (ion) {
            let prefixMassList: TheoMass[] = proteoformObj.getNMasses(ion.getName());
            /* Get matched peak list*/
            let matchedPeaks: MatchedUnMatchedPeak[] = calcMatchedPeaks.getMatchedPeakList(prefixMassList, monoMassList, seq, massErrorthVal, ppmErrorthVal, ion);
            /* copy the matched peaks to a new list for each ion selected*/
            let temp_matchedPeaks: MatchedUnMatchedPeak[] = matchedPeaks.map(x => (Object.assign({}, x)));
            matchedPeakList = matchedPeakList.concat(temp_matchedPeaks);
        });
        /* Get all the matched peaks for all the c terminus fragmented ions selected.*/
        spectrum.getCTerminalIon().forEach(function (ion: Ion) {
            /* Get claculated suffix mass list*/
            let suffixMassList: TheoMass[] = proteoformObj.getCMasses(ion.getName());
            /* Get matched peak list*/
            let matchedPeaks: MatchedUnMatchedPeak[] = calcMatchedPeaks.getMatchedPeakList(suffixMassList, monoMassList, seq, massErrorthVal, ppmErrorthVal, ion);
            /*copy the matched peaks to a new list for each ion selected*/
            let temp_matchedPeaks: MatchedUnMatchedPeak[] = matchedPeaks.map(x => (Object.assign({}, x)));
            matchedPeakList = matchedPeakList.concat(temp_matchedPeaks);
        });
        /* Get combined list of both matched and unmatched peaks to write to table*/
        matchedUnMatchedPeaks = calcMatchedPeaks.getMatchedAndUnMatchedList(monoMassList, matchedPeakList);
        //add matchedUnmatchedPeaks as decovPeaks in spectrum object
        let decovPeaksList: Peak[] = [];
        let matchedPeakPairList: MatchedPeakEnvelopePair[] = [];
        matchedUnMatchedPeaks.forEach((peak: MatchedUnMatchedPeak) => {
            let mz: number = parseFloat((peak.mass / peak.charge + 1.007276466879).toFixed(4));
            let peakObj: Peak = new Peak(peak.peakId, peak.mass, mz, peak.intensity, peak.mass, peak.charge);
            decovPeaksList.push(peakObj);
            if (peak.matchedInd == "Y") {
                let ionObj: Ion = new Ion(peak.ion, peak.ion.slice(0, 1), "", -1, peak.massError, peak.PPMerror);
                matchedPeakPairList.push(new MatchedPeakEnvelopePair(peak.thMass, peakObj, ionObj));
            }
        });
        spectrum.setDeconvPeaks(decovPeaksList);
        /*Do the below function when Sequence entered is not empty*/
        if (seq.length !== 0) {
            /*Draw SVG of Sequence*/
            let breakPoints: BreakPoints[] = formBreakPoints(matchedPeakList);
            prsmObj = new Prsm("", proteoformObj, null, [spectrum], breakPoints, matchedPeakPairList);
            let prsmViewObj = new PrsmView(Constants.SEQSVGID, prsmObj, null, true);
            prsmViewObj.getPara().setRowLength(40);
            prsmViewObj.getPara().setLetterWidth(25);
            prsmViewObj.redraw();
            $("#" + Constants.SEQSVGID).show();
            $("#" + Constants.SVGDOWNLOADID).show();
            $("#" + Constants.GRAPHDOWNLOAD).show();
            /*Get total mass and wite to HTML*/
            let totalMass: number = proteoformObj.getMass();
            setTotalSeqMass(totalMass);
            //Set Mass Difference, precursorMass is a global variable form spectrum.html
            setMassDifference(precursorMass, totalMass);
            /*draw for the prsm download modal */
            let savePrsmObj: SavePrsm = new SavePrsm(prsmViewObj);
            savePrsmObj.main();
            /**
            * Do the below function when mono mass list entered is not empty
            */
            if (monoMassListLen !== 0) {
                /*jqueryElements.monoMassTableContainer.show();
                /**
                 * Display All-peaks/matched/non-matched buttons on click of submit
                 */
                // jqueryElements.peakCount.show();
                /**
                 * Create tabe to display the input mass list data and calculated data
                 */
                /*createMonoMassTable();
                /**
                 * 	Add data to the table
                 */
                /*addMassDataToTable(matchedUnMatchedPeaks, spectrumGraphObj);
                /**
                 * function to show the count of matched peaks, un matched peaks and All peaks
                 */
                /*jqueryElements.peakCount.show();
                /**
                 * Bootstrap syntax to keep the table box to 400px height
                 * and setting properties to the table.
                 */
                /*this.setBootStarpTableProperties();
                /*showPeakCounts(monoMassList, matchedPeakList);*/
            }
        }
        /**
         * calculate envelope distribution and draw spectrum graph
         */
        if (spectrum.getPeaks().length !== 0) {
            let calcMatchedPeaks: CalcMatchedPeaks = new CalcMatchedPeaks();
            //distributionList = matchedPeaksObj.getDistribution(peakDataList,sequence,matchedUnMatchedPeaks);
            envelopeList = calcMatchedPeaks.getDistribution(modifiablePeakData, matchedUnMatchedPeaks);
            spectrum.setEnvs(envelopeList);
            /**
             * Display the graph formed
             */
            $("#" + Constants.SPECTRUMGRAPHID).show();
            //$("#"+Constants.MONOMASSGRAPHID).show();
            /**
             * Call generateCorrespondingGraph which calls addSpectrum function in invokeSpectrum file to draw graph
             */
            let spectrumDataPeaks: SpectrumFunction = new SpectrumFunction();
            let spectrumDataEnvs: SpectrumFunction = new SpectrumFunction();
            spectrumDataPeaks.assignLevelPeaks(spectrum.getPeaks());
            spectrumDataEnvs.assignLevelEnvs(spectrum.getEnvs());
            let ionList: MatchedIon[] = getIonsSpectrumGraph(matchedPeakList, spectrum.getEnvs());
            spectrumGraphObj = new SpectrumView(Constants.SPECTRUMGRAPHID, spectrum.getPeaks());
            spectrumGraphObj.addRawSpectrumAnno(spectrum.getEnvs(), ionList);
            // console.log("envPeakList:", spectrumGraphObj.envPeakList);
            spectrumGraphObj.redraw();
            //ms2GraphList.push(spectrumGraphObj);
        }
        /**
         * Local function to set the actions on click of download button in HTML
         */
        //this.download();
        let completeListofMasswithMatchedInd: MatchedUnMatchedObj[] = [];
        let nIonType: string = "B";
        let cIonType: string = "Y";
        /**
         * Code to form the second table with all the prefix masses with matched
         * masses for each ion fragment selected.
         */
        spectrum.getNTerminalIon().forEach(function (ion: Ion) {
            let calcMatchedPeaks: CalcMatchedPeaks = new CalcMatchedPeaks();
            let prefixMassList: TheoMass[] = new Array();
            let matchedAndUnMatchedList: MatchedUnMatchedPeakSimplified[] = new Array();
            let matchedAndUnMatchedListObj: MatchedUnMatchedObj = {} as MatchedUnMatchedObj;
            //let massShift = parseFloat(ion.mass);
            if (ion.getName().indexOf("A") > -1 || ion.getName().indexOf("B") > -1 || ion.getName().indexOf("C") > -1) {
                nIonType = ion.getName();
            }
            /**
             * Get calculated prefix mass
             */
            prefixMassList = proteoformObj.getNMasses(nIonType);
            prefixMassList.shift();
            prefixMassList.pop();
            /**
             * Get Matched peaks
             */
            matchedAndUnMatchedList = calcMatchedPeaks.getMatchedAndUnmatchedPrefixAndSuffixMassList(prefixMassList, monoMassList, massErrorthVal, ppmErrorthVal, "prefix");
            matchedAndUnMatchedListObj = { ionFragment: ion.getName(), massList: matchedAndUnMatchedList };
            /**
             * Complete list of all the peaks for each ion fragment
             */
            completeListofMasswithMatchedInd.push(matchedAndUnMatchedListObj);
        });
        spectrum.getCTerminalIon().forEach(function (ion: Ion) {
            let calcMatchedPeaks: CalcMatchedPeaks = new CalcMatchedPeaks();
            let suffixMassList: TheoMass[] = new Array();
            let matchedAndUnMatchedList: MatchedUnMatchedPeakSimplified[] = new Array();
            let matchedAndUnMatchedListObj: {ionFragment: string, massList: MatchedUnMatchedPeakSimplified[]} = {} as {ionFragment: string, massList: MatchedUnMatchedPeak[]};
            //let massShift = parseFloat(ion.mass);
            if (ion.getName().indexOf("X") > -1 || ion.getName().indexOf("Y") > -1 || ion.getName().indexOf("Z") > -1 || ion.getName().indexOf("Z_DOT") > -1) {
                cIonType = ion.getName();
            }
            /**
             * Get calculated prefix mass
             */
            suffixMassList = proteoformObj.getCMasses(cIonType);
            suffixMassList.shift();
            suffixMassList.pop();
            // console.log("monoMassList:",monoMassList);
            /**
             * Get Matched peaks
             */
            matchedAndUnMatchedList = calcMatchedPeaks.getMatchedAndUnmatchedPrefixAndSuffixMassList(suffixMassList, monoMassList, massErrorthVal, ppmErrorthVal, "suffix");
            matchedAndUnMatchedListObj = { ionFragment: ion.getName(), massList: matchedAndUnMatchedList };
            /**
             * Complete list of all the peaks for each ion fragment
             */
            completeListofMasswithMatchedInd.push(matchedAndUnMatchedListObj);
        });
        // console.log("completeListofMasswithMatchedInd:", completeListofMasswithMatchedInd);
        // console.log("monomasslist:", monoMassList);
        if (completeListofMasswithMatchedInd.length !== 0) {
            $("#" + Constants.H_FRAGMENTEDTABLE).show();
        }
        $("#monoMasstitle").show();
        let ions: MatchedIon[] | null = getIonsMassGraph(matchedPeakList);
        let monoMassPeakList = [];
        for (let i = 0; i < monoMassList.length; i++) {
            let monoMass = monoMassList[i].getMonoMass();
            if (monoMass) {
                let peak = new Peak(monoMassList[i].getId(), monoMass, monoMass, monoMassList[i].getIntensity());
                monoMassPeakList.push(peak);
            }
        }
        let spectrumDataMonoPeaks = new SpectrumFunction();
        spectrumDataMonoPeaks.assignLevelPeaks(monoMassPeakList);
        monoMassGraphObj = new SpectrumView(Constants.MONOMASSGRAPHID, monoMassPeakList);
        // monoMassGraphObj.para.errorThreshold = 0.06;
        monoMassGraphObj.addMonoMassSpectrumAnno(ions, proteoformObj, nIonType, cIonType);
        monoMassGraphObj.getPara().setMonoMassGraph(true);
        monoMassGraphObj.redraw();
        /**
         * add download for mono mass and spectrum graph
         */
        let saveSpectrumObj = new SaveSpectrum([spectrumGraphObj], [monoMassGraphObj]);
        saveSpectrumObj.main();
        /*draw peaks table*/
        let dataTable = new DataTable(prsmObj, false, [spectrumGraphObj]);
        dataTable.setSpecSvgId("ms2_svg_div_graph_");
        dataTable.setMonoMassSvgId("ms2_svg_div_mono_graph_");
        dataTable.drawTable();

        /* create a nav bar and a tab for ms2 graph and mono mass graph */
        clearMs2NavElement(Constants.GRAPHTABNAV);
        createMs2NavElementInspect(0, Constants.GRAPHTABDIV, Constants.GRAPHTABNAV, "");
        addCheckboxTabInspect(Constants.GRAPHTABNAV);
        addEventNavBar(monoMassGraphObj);
        /*add event handlers for spectrum graph buttons*/
        $("#ms2_graph_show_btn").click(function () {
            if ($.trim($(this).text()) === 'Show Spectrum') {
                $("#ms2_graph_show_btn").text('Hide Spectrum');
                let helpBtn: HTMLElement | null = document.getElementById("ms2_graph_help_btn");
                let saveBtn: HTMLElement | null = document.getElementById("ms2_graph_save_btn");
                let svgDiv: HTMLElement | null = document.getElementById("ms2_svg_div");
                if (!helpBtn || !saveBtn || !svgDiv) {
                    console.error("ERROR: invalid button ID or SVG div ID");
                    return;
                }
                helpBtn.style.display = "block";
                saveBtn.style.display = "block";
                svgDiv.style.display = "block";
            }
            else {
                $("#ms2_graph_show_btn").text('Show Spectrum');
                let helpBtn: HTMLElement | null = document.getElementById("ms2_graph_help_btn");
                let saveBtn: HTMLElement | null = document.getElementById("ms2_graph_save_btn");
                let svgDiv: HTMLElement | null = document.getElementById("ms2_svg_div");
                if (!helpBtn || !saveBtn || !svgDiv) {
                    console.error("ERROR: invalid button ID or SVG div ID");
                    return;
                }
                helpBtn.style.display = "none";
                saveBtn.style.display = "none";
                svgDiv.style.display = "none";
            }
        });
        // MS2 graph help button 
        $("#ms2_graph_help_btn").click(function () {
            // @ts-ignore
            $("#ms2_graph_help_popup_window").draggable({
                appendTo: "body"
            });
        });
        /**
         * Disply the table of masses for all the fragmented ions
         */
        createTableForSelectedFragmentIons(sequence, completeListofMasswithMatchedInd, monoMassGraphObj);
        this.setBootStarpropertiesforFragmentIons();
    }
    /**
     * Sets the properties of bootstrap table
     */
    /*setBootStarpTableProperties()
    {
        $("#tableContainer").DataTable({
            "scrollY": Constants.TABLEHEIGHT,
            "scrollCollapse": true,
            "paging":         false,
            "bSortClasses": false,
            "searching": false,
            "bInfo" : false,
            "columns":[
                { "type": "num" },
                { "type": "num" },
                null,
                { "type": "num" },
                { "type": "num" },
                { "type": "num" },
                null,
                { "type": "num" },
                { "type": "num" },
                { "type": "num" }
            ]
        });
    }*/
    /**
     * Sets the properties of bootstrap table
     */
    setBootStarpropertiesforFragmentIons() {
        //to be correctly sorted, column type should be num for each ion column
        //this code will work regardless of number of ions selected
        let columnCnt: number = 0;
        let columnTypes: ({"type": string} | null)[] = [];
        $("#selectedIonTableContainer .th-sm").each(function () {
            columnCnt++;
        });
        for (let i: number = 0; i < columnCnt; i++) {
            let type: {"type": string} | null = null;
            if (i != 1) {
                type = { "type": "num" };
            }
            columnTypes.push(type);
        }
        //@ts-ignore
        $("#selectedIonTableContainer").DataTable({
            "scrollY": Constants.TABLEHEIGHT,
            "scrollCollapse": true,
            "paging": false,
            "bSortClasses": false,
            "searching": false,
            "bInfo": false,
            "columns": columnTypes
        });
    }
}

/*	Get data from global variable spectrum_data and utilities to manupulate----
 * 	the data-------------------------------------------------------------------*/
PeakData = function() {
    this.peak_list = [];
    this.envelope_list = [];
    this.circleColor_list = ["red","orange","blue","green"];
    this.maxPeakIntensity;
    this.maxMz;
    
    this.getPeakData = function(json_data){
        let peakList = [];
        let i = json_data.peaks.length ;
        /*	Pass the data into peakmass and peakintensity ----------------------------*/
        while(i--)
        {
            mz = parseFloat(json_data.peaks[i].mz);
            intensity = parseFloat(json_data.peaks[i].intensity);
            peak = {mz:mz, intensity:intensity};
            peakList[i] = peak;
        }
        peakList.sort(function(x,y){
            return d3.ascending(x.mz, y.mz);
        })
        return peakList;
    }
    this.getEnvelopeData = function(json_data){
        let envelopList = [];
        json_data.envelopes.sort(function(x,y){
            return d3.ascending(x.env_peaks[0].mz, y.env_peaks[0].mz);
        })
        
        let i = json_data.envelopes.length ;
        let colorListsize = this.circleColor_list.length;
        
        while(i--)
        {
            let env_peaks = [];
            let mono_mass = parseFloat(json_data.envelopes[i].mono_mass);
            let charge = parseFloat(json_data.envelopes[i].charge);
            let color = this.circleColor_list[i%colorListsize];
            j = json_data.envelopes[i].env_peaks.length ;
            while(j--){
                let mz = parseFloat(json_data.envelopes[i].env_peaks[j].mz);
                let intensity = parseFloat(json_data.envelopes[i].env_peaks[j].intensity);
                let env_peak = {mz:mz,intensity:intensity}
                env_peaks[j] = env_peak ;
            }
            let envelope = {mono_mass:mono_mass,charge:charge,env_peaks:env_peaks,color:color};
            envelopList[i] = envelope;
        }
        return envelopList;
    }
    this.getIonData = function(prsm_data,specId,json_data){
        let envelopes =  json_data.envelopes;
        let ionData = [];
        let intensity ;
        prsm_data.prsm.ms.peaks.peak.forEach(function(element){
            let ion = "";
            if(element.hasOwnProperty('matched_ions_num'))
            {   
                ion = element.matched_ions.matched_ion.ion_type + element.matched_ions.matched_ion.ion_display_position;
            }
            if(element.spec_id == specId)
            {
                for(let i=0;i<envelopes.length;i++)
                {
                    if(parseFloat(element.monoisotopic_mass).toFixed(3) == envelopes[i].mono_mass.toFixed(3))
                    {
                        
                        intensity = envelopes[i].env_peaks.sort(function(x,y){
                                        return d3.descending(x.intensity, y.intensity);
                                    })[0].intensity; //envelopes[i].env_peaks[0].intensity;
                        //Multiplying with 1.000484 to make the ions come to center of the max peak
                        ionDataTemp = {"mz":(parseFloat(element.monoisotopic_mz)*1.000484),"intensity":parseFloat(intensity),"ion":ion};
                        ionData.push(ionDataTemp);
                        break;
                    }
                }
            }
        });
        return ionData;
    }

}
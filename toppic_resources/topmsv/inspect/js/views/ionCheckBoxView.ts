/**
 * map the ion name from the prsm.js files to the ion checkbox id in the spectrum.html
 * @param {Char} ionType - Contains ion types from check box EX. "Y,B"
 */
function setIonCheckbox(ionType: string | null): void {
	const ionNames: {[key:string]:string} = {
		"B":"b",
		"A":"a",
		"Y":"y",
		"C": "c",
		"Z_DOT":"z_",
		"X":"x"		
	}

	let ionsAll: HTMLCollectionOf<HTMLInputElement> = domElements.customControlInput;

	if(!ionType) {
		for (let i = 0; i<ionsAll.length;i++){
			if ("b" === ionsAll[i].id || "y" === ionsAll[i].id){
				ionsAll[i].checked = true;
			}
		}
		return;
	};

	let ionTypeSplit = ionType.split(',');

	ionTypeSplit.forEach(function(ionTy: string){
		for (let i = 0; i<ionsAll.length;i++){
			if (ionNames[ionTy] === ionsAll[i].id){
				ionsAll[i].checked = true;
			}
		}
	})
}

const getIonType = (): string | null => {
	return window.localStorage.getItem('ionType');
}
/**
 * receives file name to parse and replace "" and [] to make it same as the result from original JSON.parse function
 * @param {String} dataName - contains parsePeakMass in a sting format as we stringify while storing in local storage variable
 */
function parsePeakMass(dataName){
	let data = window.localStorage.getItem(dataName);
	
	console.log(data == null)
	
	if (data != "" && data != null){
		data = data.replace("[", "")
		data = data.replace("]", "")
		data = data.replace(/,/g, "\n");
		data = data.replace(/"/g, "");
	}
	return data;
}
/**
 * receives file name to parse and replace "" and [] to make it same as the result from original JSON.parse function
 * @param {String} dataName - contains sequence in a sting format as we stringify while storing in local storage variable
 */
function parseSeq(dataName){
	let data = window.localStorage.getItem(dataName);
	
	if (data != "" && data != null){
		data = data.replace("[", "")
		data = data.replace("]", "")
		data = data.replace(/"/g, "");
	}
	return data;
}
/**
 * Function returns a json formatted after getting fixed ptm data from local storage
 * @param {String} dataName - contains Fixed PTMs in a sting format as we stringify while storing in local storage variable
 */
function parsePTM(dataName){
	//ptm is an object inside array, so to preserve its structure, using JSON parse
	//it is relatively small to other lists so the performance should not deteriorate much
	let data = JSON.parse(window.localStorage.getItem(dataName));
	return data;
}
/**
 * Function returns a json formatted data after getting unknow mass list data from local storage
 * @param {String} dataName - contains unkwon mass lists in a sting format as we stringify while storing in local storage variable
 */
function parseUnknowmassList(dataName){
	let data = JSON.parse(window.localStorage.getItem(dataName));
	return data;
}
/**
 * Function returns a json formatted data after getting precursor mass data from local storage
 * @param {String} dataName - contains precursor mass in a sting format as we stringify while storing in local storage variable
 */
function parsePrecursorMass(dataName){
	let data = parseFloat(JSON.parse(window.localStorage.getItem(dataName)));
	return data;
}

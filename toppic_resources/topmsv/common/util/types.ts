/**Custom object types used throughout the typescript codes
 * replace with classes in common folder whenever possible  
 */
type AminoDist = {
  "mass":number, 
  "intensity":number
}
type Atom = {
  "atom":string;
  "count":number;
}
type InteRt = {
  "rt": number,
  "inteSum": number,
  "intePercentage": number,
  "scanNum": string
}
type BreakPoints = {
  "anno": string,
  "existCIon": boolean,
  "existNIon": boolean,
  "masses": {"charge": number, "ionDispPos": number, "ionType": string}[],
  "position": string
}
type Annotation = {
  "annoText": string, 
  "leftPos": number, 
  "rightPos": number,
  "type": ModType
}
type Residue = {
  "position": number,
  "acid": string,
  "color": string
}
type MatchedIon = {
  "env"?:Envelope, 
  "error": number, 
  "intensity": number, 
  "mz": number, 
  "text": string, 
  "pos"?: string
}
type Padding = {
  "left":number, 
  "right":number, 
  "head":number, 
  "bottom": number
}
type MatchedUnMatchedObj = {
  "ionFragment": string,
  "massList": MatchedUnMatchedPeakSimplified[]
}
type MatchedUnMatchedPeakSimplified = {
  "position": number,
  "charge": number,
  "mass": number,
  "matchedInd": string
}
type MatchedUnMatchedPeak = {
  "peakId": string,
  "ion": string,
  "ionPos": string,
  "position": number,
  "massError": number,
  "charge": number,
  "mass": number,
  "PPMerror": number,
  "thMass": number,
  "intensity": number,
  "matchedInd": string
}
enum ModType {
  Fixed,
  Unexpected,
  ProteinVariable,
  Variable
};
  
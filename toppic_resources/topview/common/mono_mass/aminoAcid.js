
function getAcidMass(acid){
	
	let aminoAcid = AminoAcid();
	aminoAcid.aa = []
	aminoAcid.aa["A"] = aminoAcid.A ;
	aminoAcid.aa["R"] = aminoAcid.R;
	aminoAcid.aa["N"] = aminoAcid.N;
	aminoAcid.aa["D"] = aminoAcid.D;
	aminoAcid.aa["C"] = aminoAcid.C;
	aminoAcid.aa["E"] = aminoAcid.E;
	aminoAcid.aa["Q"] = aminoAcid.Q;
	aminoAcid.aa["G"] = aminoAcid.G;
	aminoAcid.aa["H"] = aminoAcid.H;
	aminoAcid.aa["I"] = aminoAcid.I;
	aminoAcid.aa["L"] = aminoAcid.L;
	aminoAcid.aa["K"] = aminoAcid.K;
	aminoAcid.aa["M"] = aminoAcid.M;
	aminoAcid.aa["F"] = aminoAcid.F;
	aminoAcid.aa["P"] = aminoAcid.P;
	aminoAcid.aa["S"] = aminoAcid.S;
	aminoAcid.aa["T"] = aminoAcid.T;
	aminoAcid.aa["W"] = aminoAcid.W;
	aminoAcid.aa["Y"] = aminoAcid.Y;
	aminoAcid.aa["V"] = aminoAcid.V;
	
	return aminoAcid.aa[acid];
}
AminoAcid = function(){
	this.A = {code:"A",	shortName:"Ala",	name:"Alanine",			mono:71.037113805,	avgMass:71.0779}
	this.R = {code:"R",	shortName:"Arg",	name:"Arginine",		mono:156.101111050,	avgMass:156.18568}
	this.N = {code:"N",	shortName:"Asn",	name:"Asparagine",		mono:114.042927470,	avgMass:114.10264}
	this.D = {code:"D",	shortName:"Asp",	name:"Aspartic Acid",	mono:115.026943065,	avgMass:115.0874}
	this.C = {code:"C",	shortName:"Cys",	name:"Cysteine",		mono:103.009184505,	avgMass:103.1429}
	this.E = {code:"E",	shortName:"Glu",	name:"Glutamine",		mono:129.042593135,	avgMass:129.11398}
	this.Q = {code:"Q", shortName:"Gln",	name:"Glutamic Acid", 	mono:128.058577540, avgMass:128.12922}
	this.G = {code:"G", shortName:"Gly",	name:"Glycine",			mono:57.021463735, 	avgMass:57.05132}
	this.H = {code:"H",	shortName:"His",	name:"Histidine",		mono:137.058911875,	avgMass:137.13928}
	this.I = {code:"I",	shortName:"Ile",	name:"Isoleucine",		mono:113.084064015,	avgMass:113.15764};
	this.L = {code:"L",	shortName:"Leu",	name:"Leucine",			mono:113.084064015,	avgMass:113.15764};
	this.K = {code:"K",	shortName:"Lys",	name:"Lysine",			mono:128.094963050, avgMass:128.17228};
	this.M = {code:"M", shortName:"Met",	name:"Methionine",		mono:131.040484645,	avgMass:131.19606};
	this.F = {code:"F", shortName:"Phe",	name:"Phenylalanine", 	mono:147.068413945, avgMass:147.17386};
	this.P = {code:"P", shortName:"Pro",	name:"Proline",			mono:97.052763875,	avgMass:97.11518};
	this.S = {code:"S", shortName:"Ser",	name:"Serine",			mono:87.032028435,	avgMass:87.0773};
	this.T = {code:"T", shortName:"Thr",	name:"Threonine",		mono:101.047678505,	avgMass:101.10388};
	this.W = {code:"W", shortName:"Trp",	name:"Tryptophan",		mono:186.079312980,	avgMass:186.2099};
	this.Y = {code:"Y", shortName:"Tyr",	name:"Tyrosine",		mono:163.063328575,	avgMass:163.063328575};
	this.V = {code:"V", shortName:"Val",	name:"Valine",			mono:99.068413945,	avgMass:99.13106};
	return this ;
}

	

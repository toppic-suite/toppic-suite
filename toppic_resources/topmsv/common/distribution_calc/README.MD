# Distribution calculation library

This library contains two classes, CalculateEmassAndDistribution and MolecularFormulae. CalculateEmassAndDistribution requires aminoaciddistribution.js. MolecularFormulae requires isotopes.js.

## Usage

Use CalculateEmassAndDistribution class:
```javascript
let calEmassAndDistObj = new CalculateEmassAndDistribution();
let distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,charge,pref_suffInd, massShiftList);
```

Use MolecularFormulae class:
```javascript
let molecularFormObj = new MolecularFormulae();
let peaks = molecularFormObj.emass(mass,charge, modifiedPeakList);
```
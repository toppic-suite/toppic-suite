"use strict";
const ION_MASS_SHIFT_LIST = {
    "A": -27.9949,
    "A-H2O": -46.0149,
    "A-NH3": -45.02542,
    "B": 0,
    "B-H2O": -18.02,
    "B-NH3": -17.03052,
    "C": 17.0266,
    "C-H2O": -0.9934,
    "C-NH3": -0.00392,
    "X": 43.99,
    "X-H2O": 25.97,
    "X-NH3": 26.95948,
    "Y": 18.0106,
    "Y-H2O": 0,
    "Y-NH3": 0.98008,
    "Z": 0.984,
    "Z-H2O": -17.026,
    "Z-NH3": -16.04652,
    "Z_DOT": 1.9919,
    "Z_DOT-H2O": -16.018664,
    "Z_DOT-NH3": -15.03862
};
function getIonMassShift(ionType) {
    if (ION_MASS_SHIFT_LIST.hasOwnProperty(ionType.toUpperCase())) {
        return ION_MASS_SHIFT_LIST[ionType.toUpperCase()];
    }
    else {
        return undefined;
    }
}
function getWaterMass() {
    return 18.01056468362;
}

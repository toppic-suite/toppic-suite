"use strict";
/**
 * @function getIsotopicMassOfAtom
 * @description
 * Code to get absolute mass and intensity values of atoms
 */
const getIsotopicMassOfAtom = function (atom) {
    if (atom == "C")
        return getCarbonIsotope();
    else if (atom == "H")
        return getHydrogenIsotope();
    else if (atom == "N")
        return getNitrogenIsotope();
    else if (atom == "O")
        return getOxygenIsotope();
    else if (atom == "S")
        return getSulfurIsotope();
    else
        return [];
};
const getCarbonIsotope = function () {
    let l_acarbonIsotopes = [];
    let C12 = { mass: 12.000000, intensity: 98.93 };
    let C13 = { mass: 13.003355, intensity: 1.07 };
    l_acarbonIsotopes.push(C12);
    l_acarbonIsotopes.push(C13);
    return l_acarbonIsotopes;
};
const getHydrogenIsotope = function () {
    let l_ahydrogenIsotopes = [];
    let H1 = { mass: 1.007825, intensity: 99.98 };
    let H2 = { mass: 2.014102, intensity: 0.02 };
    l_ahydrogenIsotopes.push(H1);
    l_ahydrogenIsotopes.push(H2);
    return l_ahydrogenIsotopes;
};
const getOxygenIsotope = function () {
    let l_aoxygenIsotopes = [];
    let O16 = { mass: 15.994915, intensity: 99.757 };
    let O17 = { mass: 16.999131, intensity: 0.038 };
    let O18 = { mass: 17.999159, intensity: 0.205 };
    l_aoxygenIsotopes.push(O16);
    l_aoxygenIsotopes.push(O17);
    l_aoxygenIsotopes.push(O18);
    return l_aoxygenIsotopes;
};
const getNitrogenIsotope = function () {
    let l_anitrogenIsotopes = [];
    let N14 = { mass: 14.003074, intensity: 99.63 };
    let N15 = { mass: 15.000109, intensity: 0.37 };
    l_anitrogenIsotopes.push(N14);
    l_anitrogenIsotopes.push(N15);
    return l_anitrogenIsotopes;
};
const getSulfurIsotope = function () {
    let l_asulfurIsotopes = [];
    let S32 = { mass: 31.972072, intensity: 94.99 };
    let S33 = { mass: 32.971459, intensity: 0.75 };
    let S34 = { mass: 33.967868, intensity: 4.25 };
    let S36 = { mass: 35.967079, intensity: 0.01 };
    l_asulfurIsotopes.push(S32);
    l_asulfurIsotopes.push(S33);
    l_asulfurIsotopes.push(S34);
    l_asulfurIsotopes.push(S36);
    return l_asulfurIsotopes;
};
/**
 * Code to calculate distributions by taking relative
 * values of mass and intensity of atoms.
 */
const getIsotopicMassRef = function (atom) {
    if (atom == "C")
        return getCarbonRef();
    else if (atom == "H")
        return getHydrogenRef();
    else if (atom == "N")
        return getNitrogenRef();
    else if (atom == "O")
        return getOxygenRef();
    else if (atom == "S")
        return getSulfurRef();
    else
        return [];
};
const getCarbonRef = function () {
    let l_acarbonIsotopes = [];
    let C12 = { mass: 0, intensity: 100 };
    let C13 = { mass: 1, intensity: 1.0815728 };
    l_acarbonIsotopes.push(C12);
    l_acarbonIsotopes.push(C13);
    return l_acarbonIsotopes;
};
const getHydrogenRef = function () {
    let l_ahydrogenIsotopes = [];
    let H1 = { mass: 0, intensity: 100 };
    let H2 = { mass: 1, intensity: 0.0115013 };
    l_ahydrogenIsotopes.push(H1);
    l_ahydrogenIsotopes.push(H2);
    return l_ahydrogenIsotopes;
};
const getOxygenRef = function () {
    let l_aoxygenIsotopes = [];
    let O16 = { mass: 0, intensity: 100 };
    let O17 = { mass: 1, intensity: 0.0380926 };
    let O18 = { mass: 2, intensity: 0.2054994 };
    l_aoxygenIsotopes.push(O16);
    l_aoxygenIsotopes.push(O17);
    l_aoxygenIsotopes.push(O18);
    return l_aoxygenIsotopes;
};
const getNitrogenRef = function () {
    let l_anitrogenIsotopes = [];
    let N14 = { mass: 0, intensity: 100 };
    let N15 = { mass: 1, intensity: 0.3653296 };
    l_anitrogenIsotopes.push(N14);
    l_anitrogenIsotopes.push(N15);
    return l_anitrogenIsotopes;
};
const getSulfurRef = function () {
    let l_asulfurIsotopes = [];
    let S32 = { mass: 0, intensity: 100 };
    let S33 = { mass: 1, intensity: 0.7895568 };
    let S34 = { mass: 2, intensity: 4.4741552 };
    let S36 = { mass: 4, intensity: 0.0105274 };
    l_asulfurIsotopes.push(S32);
    l_asulfurIsotopes.push(S33);
    l_asulfurIsotopes.push(S34);
    l_asulfurIsotopes.push(S36);
    return l_asulfurIsotopes;
};

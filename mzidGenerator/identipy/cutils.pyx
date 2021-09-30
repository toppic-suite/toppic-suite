cimport cython
from cpython.sequence cimport PySequence_GetSlice
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem

from pyteomics import cmass
from math import factorial

cimport pyteomics.cmass as cmass

from pyteomics import electrochem as ec
import numpy as np
cimport numpy as np

np.import_array()


cdef:
    dict std_aa_mass = cmass.std_aa_mass
    dict std_ion_comp = cmass.std_ion_comp
    dict nist_mass = cmass._nist_mass

cdef dict ion_shift_dict

ion_shift_dict = {
    'a': 46.00547930326002,
    'b': 18.010564683699954,
    'c': 0.984015582689949,
    'x': -25.979264555419945,
    'y': 0.0,
    'z': 17.026549101010005,
}

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
def RNHS_fast(set spectrum_fastset, dict spectrum_idict , dict theoretical_set, int min_matched):
    cdef int matched_approx_b, matched_approx_y, matched_approx
    cdef set matched_b, matched_y
    matched_b = spectrum_fastset.intersection(theoretical_set['b'])
    matched_y = spectrum_fastset.intersection(theoretical_set['y'])
    matched_approx_b = len(matched_b)
    matched_approx_y = len(matched_y)
    matched_approx = matched_approx_b + matched_approx_y
    if matched_approx >= min_matched:
        isum = 0
        for fr in matched_b:
            isum += spectrum_idict[fr]
        for fr in matched_y:
            isum += spectrum_idict[fr]
        return matched_approx, factorial(matched_approx_b) * factorial(matched_approx_y) * isum
    else:
        return 0, 0


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef float calc_ions_from_neutral_mass(str peptide, float nm, str ion_type, int charge, dict aa_mass, float cterm_mass, float nterm_mass):
    cdef float nmi
    if ion_type in 'abc':
        nmi = nm - aa_mass[peptide[-1]] - ion_shift_dict[ion_type] - (cterm_mass - 17.002735)
    else:
        nmi = nm - aa_mass[peptide[0]] - ion_shift_dict[ion_type] - (nterm_mass - 1.007825)
    return (nmi + 1.0072764667700085 * charge) / charge 

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef list get_n_ions(str peptide, float maxmass, int pl, int charge, dict k_aa_mass):
    cdef int i
    cdef list tmp
    tmp = [maxmass, ]
    for i in xrange(1, pl):
        tmp.append(tmp[-1] - k_aa_mass[peptide[-i-1]]/charge)
    return tmp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef list get_c_ions(str peptide, float maxmass, int pl, int charge, dict k_aa_mass):
    cdef int i
    cdef list tmp
    tmp = [maxmass, ]
    for i in xrange(pl-2, -1, -1):
        tmp.append(tmp[-1] - k_aa_mass[peptide[-(i+2)]]/charge)
    return tmp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef tuple ctheor_spectrum(str peptide, double acc_frag, double nterm_mass, double cterm_mass, tuple types,
                           int maxcharge, bint reshape, dict kwargs):
    cdef int pl, charge, i, n, i_type, n_types
    cdef bint nterminal
    cdef str ion_type, maxpart, part
    cdef float maxmass, part_mass, nm
    cdef dict peaks, theoretical_set
    cdef dict aa_mass, ion_comp, mass_data
    cdef set theoretical_set_item
    cdef list ions_scaled, marr
    cdef object marr_storage

    aa_mass = kwargs.get("aa_mass")
    if aa_mass is None:
        aa_mass = std_aa_mass
    ion_comp = kwargs.get("ion_comp")
    if ion_comp is None:
        ion_comp = std_ion_comp
    mass_data = kwargs.get("mass_data")
    if mass_data is None:
        mass_data = nist_mass
    nm = kwargs.get("nm")
    if nm is None:
        nm = cmass.fast_mass(peptide, **kwargs) + (nterm_mass - 1.007825) + (cterm_mass - 17.002735)
    peaks = {}
    theoretical_set = dict()

    pl = len(peptide) - 1
    n_types = len(types)
    for charge in range(1, maxcharge + 1):
        for i_type in range(n_types):
            ion_type = <str>PyTuple_GetItem(types, i_type)
            nterminal = ion_type[0] in 'abc'
            if nterminal:
                maxmass = calc_ions_from_neutral_mass(peptide, nm, ion_type=ion_type, charge=charge,
                                aa_mass=kwargs['aa_mass'], cterm_mass=cterm_mass, nterm_mass=nterm_mass)
                marr = get_n_ions(peptide, maxmass, pl, charge, kwargs['aa_mass'])
            else:
                maxmass = calc_ions_from_neutral_mass(peptide, nm, ion_type=ion_type, charge=charge,
                                aa_mass=kwargs['aa_mass'], cterm_mass=cterm_mass, nterm_mass=nterm_mass)
                marr = get_c_ions(peptide, maxmass, pl, charge, kwargs['aa_mass'])

            ions_scaled = [<int>(x / acc_frag) for x in marr]
            if ion_type in theoretical_set:
                theoretical_set_item = <set>PyDict_GetItem(theoretical_set, ion_type)
                theoretical_set_item.update(ions_scaled)
            else:
                theoretical_set[ion_type] = ions_scaled

            if reshape:
                marr_storage = np.array(marr)
                n = marr_storage.size
                marr_storage = marr_storage.reshape((n, 1))
                peaks[ion_type, charge] = marr_storage
            else:
                peaks[ion_type, charge] = marr
    return peaks, theoretical_set


def theor_spectrum(peptide, acc_frag, nterm_mass, cterm_mass, types=('b', 'y'), maxcharge=None, reshape=False, **kwargs):
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    return ctheor_spectrum(peptide, acc_frag, nterm_mass, cterm_mass, tuple(types), maxcharge, reshape, kwargs)

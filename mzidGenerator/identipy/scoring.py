from .utils import neutral_masses, theor_spectrum, get_aa_mass
from scipy.spatial import cKDTree
import numpy as np
from math import factorial
from copy import copy


def simple_score(spectrum, peptide, settings):
    acc = settings.getfloat('search', 'product accuracy')
    int_array = spectrum['intensity array']
    int_array = int_array / int_array.max() * 100
    charge = max(c for _, c in neutral_masses(spectrum, settings))
    theor = theor_spectrum(peptide, maxcharge=charge, aa_mass=get_aa_mass(settings),
        nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage'), cterm_mass=settings.getfloat('modifications', 'protein cterm cleavage'))
    fragments = np.concatenate(theor.values())
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
            (spectrum['m/z array'].size, 1)))
    dist, ind = spectrum['__KDTree'].query(fragments.reshape((fragments.size, 1)),
            distance_upper_bound=acc)
    mask = dist != np.inf
    if mask.size < settings.getint('scoring', 'minimum matched'):
        return -1
    return spectrum['intensity array'][ind[mask]].sum()


def get_fragment_mass_tol(spectrum, peptide, settings):
    """A function for obtaining optimal fragment mass tolerance, dynamic range"""
    acc = settings.getfloat('search', 'product accuracy')
    spectrum = copy(spectrum)
    idx = np.nonzero(spectrum['m/z array'] >= 150)
    spectrum['intensity array'] = spectrum['intensity array'][idx]
    spectrum['m/z array'] = spectrum['m/z array'][idx]
    int_array = spectrum['intensity array']
    int_array = int_array / int_array.max() * 100
    charge = 1#max(1, max(c for _, c in neutral_masses(spectrum, settings)) - 1)
    theor, _ = theor_spectrum(peptide, maxcharge=charge, aa_mass=get_aa_mass(settings), reshape=True, acc_frag=acc,
        nterm_mass=settings.getfloat('modifications', 'protein nterm cleavage'), cterm_mass=settings.getfloat('modifications', 'protein cterm cleavage'))
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
            (spectrum['m/z array'].size, 1)))

    dist_total, int_array_total = np.array([]), np.array([])
    for fragments in theor.values():
        n = fragments.size
        dist, ind = spectrum['__KDTree'].query(fragments.reshape((n, 1)),
            distance_upper_bound=acc)
        mask = (dist != np.inf)
        int_array_total = np.append(int_array_total, int_array[ind[mask]])
        dist_total = np.append(dist_total, dist[mask] / spectrum['m/z array'][ind[mask]] * 1e6)

    new_params = {}
    if dist_total.size:
        new_params['fmt'] = dist_total#2 * np.median(dist_total)
    else:
        new_params['fmt'] = []
    return new_params

def morpheusscore_fast(spectrum_fastset, spectrum_idict, theoretical_set, min_matched):
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
        return matched_approx, matched_approx + isum
        # return matched_approx, factorial(matched_approx_b) * (100 * matched_approx_b) + factorial(matched_approx_y) * (100 * matched_approx_y)
        # return matched_approx, factorial(matched_approx) * (100 * matched_approx)
    else:
        return 0, 0

def morpheusscore(spectrum, theoretical, acc, acc_ppm=False, position=False):
    if 'norm' not in spectrum:
        spectrum['norm'] = spectrum['Isum']#spectrum['intensity array'].sum()#spectrum['intensity array'].max() / 100.
    mz_array = spectrum['m/z array']
    score = 0
    match = {}
    match2 = {}
    total_matched = 0
    sumI = 0
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(mz_array.reshape((mz_array.size, 1)))

    dist_all = []
    for ion, fragments in theoretical.iteritems():
        dist, ind = spectrum['__KDTree'].query(fragments, distance_upper_bound=acc)
        mask1 = (dist != np.inf)
        if acc_ppm:
            mask2 = (dist[mask1] / spectrum['m/z array'][ind[mask1]] * 1e6 <= acc_ppm)
        else:
            mask2 = np.ones_like(dist[mask1], dtype=bool)
        nmatched = mask2.sum()
        if nmatched:
            total_matched += nmatched
            sumi = spectrum['intensity array'][ind[mask1][mask2]].sum()
            sumI += sumi
            score += sumi / spectrum['norm']
            dist_all.extend(dist[mask1][mask2])
        match[ion] = mask2
        match2[ion] = mask1
    if not total_matched:
        return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
    if position:
        yions = match2[('y', 1)]
        bions = match2[('b', 1)]
        plen = len(yions) + 1
        if position == 1:
            if not bions[0]:
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        elif position == plen:
            if not yions[0]:
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        else:
            if not (yions[plen - position] and yions[plen - position - 1]) or (bions[position - 1] and bions[position - 2]):
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}

    score += total_matched
    sumI = np.log10(sumI)

    return {'score': score, 'match': match, 'sumI': sumI, 'dist': dist_all, 'total_matched': total_matched}

def hyperscore_fast(spectrum_fastset, spectrum_idict, theoretical_set, min_matched):
    matched_b = spectrum_fastset.intersection(theoretical_set['b'])
    matched_y = spectrum_fastset.intersection(theoretical_set['y'])
    matched_approx_b = len(matched_b)
    matched_approx_y = len(matched_y)
    #matched_approx_b = len(spectrum_fastset.intersection(theoretical_set['b']))
    #matched_approx_y = len(spectrum_fastset.intersection(theoretical_set['y']))
    matched_approx = matched_approx_b + matched_approx_y
    if matched_approx >= min_matched:
        isum = 0
        for fr in matched_b:
            isum += spectrum_idict[fr]
        for fr in matched_y:
            isum += spectrum_idict[fr]
        # return matched_approx, factorial(matched_approx_b) * factorial(matched_approx_y)
        return matched_approx, factorial(matched_approx_b) * 100 * isum * (matched_approx_b + matched_approx_y) * factorial(matched_approx_y)
        # return matched_approx, factorial(matched_approx) * (100 * matched_approx)
    else:
        return 0, 0

def hyperscore(spectrum, theoretical, acc, acc_ppm=False, position=False):
    if 'norm' not in spectrum:
        spectrum['norm'] = spectrum['intensity array'].max() / 100.
    mz_array = spectrum['m/z array']
    score = 0
    mult = []
    match = {}
    match2 = {}
    total_matched = 0
    sumI = 0
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(mz_array.reshape((mz_array.size, 1)))

    dist_all = []
    for ion, fragments in theoretical.iteritems():
        dist, ind = spectrum['__KDTree'].query(fragments, distance_upper_bound=acc)
        mask1 = (dist != np.inf)
        if acc_ppm:
            mask2 = (dist[mask1] / spectrum['m/z array'][ind[mask1]] * 1e6 <= acc_ppm)
        else:
            mask2 = np.ones_like(dist[mask1], dtype=bool)
        nmatched = mask2.sum()
        if nmatched:
            total_matched += nmatched
            mult.append(factorial(nmatched))
            sumi = spectrum['intensity array'][ind[mask1][mask2]].sum()
            sumI += sumi
            score += sumi / spectrum['norm']
            dist_all.extend(dist[mask1][mask2])
        match[ion] = mask2
        match2[ion] = mask1
    if not total_matched:
        return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
    if position:
        yions = match2[('y', 1)]
        bions = match2[('b', 1)]
        plen = len(yions) + 1
        if position == 1:
            if not bions[0]:
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        elif position == plen:
            if not yions[0]:
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        else:
            if not (yions[plen - position] and yions[plen - position - 1]) or (bions[position - 1] and bions[position - 2]):
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}

    for m in mult:
        score *= m
    sumI = np.log10(sumI)

    return {'score': score, 'match': match, 'sumI': sumI, 'dist': dist_all, 'total_matched': total_matched}

def RNHS_fast(spectrum_fastset, spectrum_idict, theoretical_set, min_matched):
    matched_b = spectrum_fastset.intersection(theoretical_set['b'])
    matched_y = spectrum_fastset.intersection(theoretical_set['y'])
    matched_approx_b = len(matched_b)
    matched_approx_y = len(matched_y)
    #matched_approx_b = len(spectrum_fastset.intersection(theoretical_set['b']))
    #matched_approx_y = len(spectrum_fastset.intersection(theoretical_set['y']))
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

def RNHS(spectrum, theoretical, acc, acc_ppm=False, position=False):
    if 'norm' not in spectrum:
        spectrum['norm'] = spectrum['Isum']#spectrum['intensity array'].sum()#spectrum['intensity array'].max() / 100.
    mz_array = spectrum['m/z array']
    score = 0
    mult = []
    match = {}
    match2 = {}
    total_matched = 0
    sumI = 0
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(mz_array.reshape((mz_array.size, 1)))

    dist_all = []
    for ion, fragments in theoretical.iteritems():
        dist, ind = spectrum['__KDTree'].query(fragments, distance_upper_bound=acc)
        mask1 = (dist != np.inf)
        if acc_ppm:
            nacc = np.where(dist[mask1] / mz_array[ind[mask1]] * 1e6 > acc_ppm)[0]
            mask2 = mask1.copy()
            mask2[nacc] = False
        else:
            mask2 = mask1
        nmatched = mask2.sum()
        if nmatched:
            total_matched += nmatched
            mult.append(factorial(nmatched))
            sumi = spectrum['intensity array'][ind[mask2]].sum()
            sumI += sumi
            score += sumi / spectrum['norm']
            dist_all.extend(dist[mask2])
        match[ion] = mask2
        match2[ion] = mask2
    if not total_matched:
        return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
    if position:
        yions = match2[('y', 1)]
        bions = match2[('b', 1)]
        plen = len(yions)
        if position > plen + 1:
#           print 'Something wrong with aachange position'
            return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        if position == 1:
            if not bions[0]:
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        elif position == plen + 1:
            if not yions[0]:
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}
        else:
            if not (yions[plen - position + 1] and yions[plen - position]):
                return {'score': 0, 'match': None, 'sumI': 0, 'dist': [], 'total_matched': 0}

    for m in mult:
        score *= m
    sumI = np.log10(sumI)

    return {'score': score, 'match': match, 'sumI': sumI, 'dist': dist_all, 'total_matched': total_matched}

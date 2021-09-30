from scipy.stats import percentileofscore, scoreatpercentile
from scipy.optimize import curve_fit
from pyteomics import achrom, auxiliary as aux, parser, mass
from collections import Counter
from .main import *
from .scoring import get_fragment_mass_tol
import logging
logger = logging.getLogger(__name__)
import numpy as np
from .utils import get_info, get_aa_mass, get_enzyme, calculate_RT
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

def FDbinSize(X):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    X:  1D Data set
    Returns:
    h:  F-D bin size
    """
    X = np.sort(X)
    upperQuartile = scoreatpercentile(X, 75)
    lowerQuartile = scoreatpercentile(X, 25)
    IQR = upperQuartile - lowerQuartile
    h = 2. * IQR / len(X) ** (1. / 3.)
    return h


def get_subset(results, settings, fdr=0.01):
    """Filter results to given FDR using top 1 candidates"""
    subset = aux.filter(results, key=lambda x: x['e-values'][0],
            is_decoy = lambda x: x['candidates'][0][2] == 'd',
            fdr=fdr)
    return subset

def optimization(fname, settings):
    settings = settings.copy()
    settings.set('misc', 'first stage', '')
    efc = settings.get('scoring', 'e-values for candidates')
    settings.set('scoring', 'e-values for candidates', 1)
    left = settings.getfloat('search', 'precursor accuracy left')
    right = settings.getfloat('search', 'precursor accuracy right')
    wide = settings.getboolean('optimization', 'increase precursor mass tolerance')
    if settings.get('search', 'precursor accuracy unit') != 'ppm':
        left *= 1000
        right *= 1000
    if left < 100 and wide:
        settings.set('search', 'precursor accuracy left', 100)
    if right < 100 and wide:
        settings.set('search', 'precursor accuracy right', 100)
    settings.set('search', 'precursor accuracy unit', 'ppm')
    results = process_file(fname, settings, initial_run=False)
    filtered = get_subset(results, settings, fdr=0.01)
    logger.info('%s PSMs with 1%% FDR.', len(filtered))
    if len(filtered) < 50:
        if len(filtered) < 10:
            logger.warning('OPTIMIZATION ABORTED')
            return settings
        else:
            functions = [precursor_mass_optimization, fragment_mass_optimization,
                    missed_cleavages_optimization]
    else:
        functions = [
                rt_filtering,
                precursor_mass_optimization,
                fragment_mass_optimization,
#               missed_cleavages_optimization
                ]
    for func in functions:
        settings = func(filtered, settings)
    settings.set('scoring', 'e-values for candidates', efc)
    return settings


def charge_optimization(results, settings):
    settings = settings.copy()
    chargestates = np.array([get_info(res['spectrum'], res, settings)[1] for res in results])
    mincharge = chargestates.min()
    maxcharge = chargestates.max()
    
    for ch in range(mincharge, maxcharge+1):
        if float(chargestates[chargestates < ch].size) / chargestates.size < 0.01:
            mincharge = ch
    for ch in range(maxcharge, mincharge-1, -1):
        if float(chargestates[chargestates > ch].size) / chargestates.size < 0.01:
            maxcharge = ch
    logger.info('NEW charges = %s:%s', mincharge, maxcharge)
    settings.set('search', 'maximum charge', maxcharge)
    settings.set('search', 'minimum charge', mincharge)
    return settings

def calibrate_mass(bwidth, mass_left, mass_right, true_md):
    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]

    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
    mass_shift, mass_sigma = popt[1], np.abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]

def noisygaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

def precursor_mass_optimization(results, settings):
    settings_nopime = settings.copy()
    settings_nopime.set('search', 'precursor isotope mass error', '0')
    settings_nopime.set('search', 'shifts', '0')
#   results = get_output(results, settings_nopime)

    settings = settings.copy()
    mass_left = settings.getfloat('search', 'precursor accuracy left')
    mass_right = settings.getfloat('search', 'precursor accuracy right')
    massdif = np.array([res['candidates'][0][4]['mzdiff']['ppm'] for res in results])
    massdif = massdif[(massdif > -mass_left) & (massdif < mass_right)]
    if settings.get('search', 'precursor accuracy unit') != 'ppm':
        mass_left = mass_left * 1e6 / 400
        mass_right = mass_right * 1e6 / 400
    logger.info('mass_left, mass_right: %s, %s', mass_left, mass_right)
    try:
        mass_shift, mass_sigma, covvalue = calibrate_mass(0.1, mass_left, mass_right, massdif)
        if np.isinf(covvalue):
            mass_shift, mass_sigma, covvalue = calibrate_mass(0.01, mass_left, mass_right, massdif)
        logger.info('%s, %s -> %s +- 8 * %s; %s', mass_left, mass_right, mass_shift, mass_sigma, covvalue)
        best_par_mt_l = mass_shift - 8 * mass_sigma
        best_par_mt_r = mass_shift + 8 * mass_sigma
        logger.info('SMART MASS TOLERANCE = %s:%s', best_par_mt_l, best_par_mt_r)
    except RuntimeError:
        error = True
    else:
        error = False
    if not error and np.isinf(covvalue):
        error = True
        logger.warning('Double error when fitting precursor errors: %s', massdif)
    if error or (percentileofscore(massdif, best_par_mt_r) - percentileofscore(massdif, best_par_mt_l) < 95):
        best_par_mt_l = scoreatpercentile(massdif, 0.1)
        best_par_mt_r = scoreatpercentile(massdif, 99.9)
        logger.warning('Percentage sanity check FAILED. Falling back on percentage boundaries')
    else:
        best_par_mt_l = max(best_par_mt_l, scoreatpercentile(massdif, 0.1))
        best_par_mt_r = min(best_par_mt_r, scoreatpercentile(massdif, 99.9))

    logger.info('NEW PARENT MASS TOLERANCE = %s:%s', best_par_mt_l, best_par_mt_r)
    settings.set('search', 'precursor accuracy left', -best_par_mt_l)
    settings.set('search', 'precursor accuracy right', best_par_mt_r)
    settings.set('search', 'precursor accuracy unit', 'ppm')
    return settings

def missed_cleavages_optimization(results, settings):
    settings = settings.copy()
    missedcleavages = np.array([parser.num_sites(str(res['candidates'][0][1]), get_enzyme(str(settings.get('search', 'enzyme'))))
        for res in results])
    best_missedcleavages = missedcleavages.max()
    for mc in range(best_missedcleavages, -1, -1):
        if float(missedcleavages[missedcleavages > mc].size) / missedcleavages.size < 0.002:
            best_missedcleavages = mc
    logger.info('NEW miscleavages = %s', best_missedcleavages)
    settings.set('search', 'number of missed cleavages', best_missedcleavages)
    return settings

def fragment_mass_optimization(results, settings):
    settings = settings.copy()
    fragmassdif = []
    for res in results:
        fragmassdif.extend(get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings)['fmt'])
    fragmassdif = np.array(fragmassdif)

    best_frag_mt = scoreatpercentile(fragmassdif, 68) * 4    

    logger.info('NEW FRAGMENT MASS TOLERANCE ppm = %s', best_frag_mt)
    settings.set('search', 'product accuracy ppm', best_frag_mt)
    settings.set('search', 'product accuracy unit', 'ppm')
    return settings


def rt_filtering(results, settings):
    settings = settings.copy()
    if settings.has_option('misc', 'legend'):
        legend = settings.get('misc', 'legend')
    else:
        legend = None
    RTexp, seqs = zip(*[(utils.get_RT(res['spectrum']), res['candidates'][0][1]) for res in results])
    if legend is not None:
        stdl = set(parser.std_labels)
        newseqs = []
        for s in seqs:
            if parser.fast_valid(s):
                newseqs.append(list(s))
            else:
                seq = []
                c, n = False, False
                for c in s:
                    if c in stdl:
                        seq.append(c)
                    else:
                        mod, res, term = legend[c]
                        if res == '-':
                            if term == '[':
                                seq.append(mod+'-')
                                n = True
                            else:
                                seq.append('-'+mod)
                                c = True
                        else:
                            seq.append(mod+res)
                    if not n: seq.append(parser.std_nterm)
                    if not c: seq.append(parser.std_cterm)
                newseqs.append(seq)
        seqs = newseqs
    RTexp = [float(x) for x in RTexp]
    if np.allclose(RTexp, 0):
        logger.warning('RT is missing. Skipping RT optimization.')
        return settings
    RC_def = achrom.RCs_gilar_rp
    xdict = {}
    for key, val in RC_def['aa'].items():
        xdict[key] = [val, None]
    RC_dict = utils.get_RCs_vary_lcp(seqs, RTexp)
    RC_dict_new = dict()
    for key, val in RC_dict['aa'].items():
        xdict.setdefault(key, [val, None])[1] = val
    a, b, _, _ = aux.linear_regression([x[0] for x in xdict.values() if x[1] != None], [x[1] for x in xdict.values() if x[1] != None])
    for key, x in xdict.items():
        if x[1] == None:
            x[1] = x[0] * a + b
        RC_dict_new[key] = x[1]
    if legend is not None:
        for k, v in legend.items():
            if len(k) == 1: continue
            if k[-1] in '[]':
                if k[-2] == '-':
                    kk = ('-' + k[1:-1]) if k[-1] == ']' else (k[:-1])
                else:
                    kk = k[:-1]
            elif len(k) > 1:
                kk = k
            logger.debug('%s -> %s', k, kk)
            if kk in RC_dict_new:
                RC_dict_new[v] = RC_dict_new[kk]
            else:
                if kk[-1].isupper():
                    kkk = kk[-1]
                elif kk[-1] == '-':
                    kkk = parser.std_nterm
                elif kk[0] == '-':
                    kkk = parser.std_cterm
                RC_dict_new[v] = RC_dict_new.get(kkk, 0)
                logger.info('No RC for %s, using %s or 0: %s', kk, kkk, RC_dict_new[v])


    RC_dict['aa'] = RC_dict_new

    logger.debug('RC dict: %s', RC_dict)
    rtexp = np.array([np.mean(x) for x in RTexp])
    rttheor = np.array([calculate_RT(pep, RC_dict, raise_no_mod=False)
        for pep in seqs])
    deltaRT = rtexp - rttheor
    logger.debug('Linear regression: %s', aux.linear_regression(rtexp, rttheor))
    best_RT_l = scoreatpercentile(deltaRT, 0.05)
    best_RT_r = scoreatpercentile(deltaRT, 99.95)

    def condition(spectrum, cand, _, stored_value=False):
        if not stored_value:
            stored_value = calculate_RT(cand, RC_dict)
        rtd = spectrum['RT'] - stored_value
        return best_RT_l <= rtd <= best_RT_r, stored_value
    settings.set('scoring', 'condition', condition)
    return settings

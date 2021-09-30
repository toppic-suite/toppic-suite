import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, fasta, mzml, parser, mgf
import pandas as pd
import sys
from itertools import combinations, islice
from collections import defaultdict, Counter
import numpy as np
from multiprocessing import Queue, Process, cpu_count
import string
from copy import copy
try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser
import tempfile
import os
import logging
import itertools as it
try:
    from lxml import etree
except ImportError:
    etree = None
from time import strftime
from os import path
logger = logging.getLogger(__name__)

try:
    from pyteomics import cmass
except ImportError:
    logger.warning('pyteomics.cythonize not found. It is highly recommended for good performance.')
    cmass = mass
try:
    import pyximport; pyximport.install()
    from . import cparser
except:
    from . import customparser as cparser

default_tags = {
'tmt10plex': {
    'tmt_126': 126.1277261,
    'tmt_127N': 127.1247610,
    'tmt_128C': 128.1344357,
    'tmt_129N': 129.1314706,
    'tmt_130C': 130.1411453,
    'tmt_131': 131.1381802,
    'tmt_127C': 127.1310809,
    'tmt_128N': 128.1281158,
    'tmt_129C': 129.1377905,
    'tmt_130N': 130.1348254
},
'tmt11plex': {
    'tmt_126': 126.1277261,
    'tmt_127N': 127.1247610,
    'tmt_128C': 128.1344357,
    'tmt_129N': 129.1314706,
    'tmt_130C': 130.1411453,
    'tmt_131': 131.1381802,
    'tmt_127C': 127.1310809,
    'tmt_128N': 128.1281158,
    'tmt_129C': 129.1377905,
    'tmt_130N': 130.1348254,
    'tmt_131C': 130.144999
},
'tmt6plex':{
    'tmt_126': 126.1277261,
    'tmt_127N': 127.1247610,
    'tmt_128C': 128.1344357,
    'tmt_129N': 129.1314706,
    'tmt_130C': 130.1411453,
    'tmt_131': 131.1381802,
}
}
def get_peptides(prot_seq, enzyme, mc, minlen, maxlen, semitryptic=False):
    peptides = cparser._cleave(prot_seq, enzyme, mc)
    for pep, startposition in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen:
            if not semitryptic:
                yield pep, startposition, plen
            else:
                for i in range(plen-minlen+1):
                    yield pep[i:], startposition + i, plen - i
                for i in range(1, plen-minlen+1, 1):
                    yield pep[:-i], startposition, plen - i

seen_target = set()
seen_decoy = set()

def prot_peptides(prot_seq, enzyme, mc, minlen, maxlen, is_decoy, dont_use_seen_peptides=False, snp=False, desc=False, position=False, semitryptic=False):

    dont_use_fast_valid = parser.fast_valid(prot_seq)
    methionine_check = prot_seq[0] == 'M'
    if snp == 2:
        if desc:
            try:
                tmp = desc.split(' ')[0].split('|')
                pos = int(tmp[1]) - 1
                aach = tmp[2]
            except:
                desc = False
    # peptides = cparser._cleave(prot_seq, enzyme, mc)
    # for pep, startposition in peptides:
    #     plen = len(pep)
    for pep, startposition, plen in get_peptides(prot_seq, enzyme, mc, minlen, maxlen, semitryptic):
        loopcnt = 0
        if pep not in seen_target and pep not in seen_decoy and (dont_use_fast_valid or parser.fast_valid(pep)):
            loopcnt = 1
            if methionine_check and startposition == 0:
                if minlen <= plen - 2:
                    loopcnt = 3
                elif minlen <= plen - 1:
                    loopcnt = 2
        while loopcnt:
            f = pep[loopcnt-1:]
            if dont_use_seen_peptides:
                if snp == 1:
                    for ff, seq_new in custom_snp(f, startposition):
                        if not seq_new:
                            yield ff if not position else (ff, startposition)
                        else:
                            yield ff if not position else (ff, startposition)
                else:
                    yield f if not position else (f, startposition)
            else:
                if f not in seen_target and f not in seen_decoy:
                    if is_decoy:
                        seen_decoy.add(f)
                    else:
                        seen_target.add(f)
                    if snp == 1:
                        for ff, seq_new in custom_snp(f, startposition):
                            if not seq_new:
                                yield ff if not position else (ff, startposition)
                            if seq_new not in seen_decoy and seq_new not in seen_target:
                                yield ff if not position else (ff, startposition)
                    elif snp == 2:
                        if desc and startposition <= pos <= startposition + plen:
                            if len(aach) == 3 and aach[0] in parser.std_amino_acids and aach[2] in parser.std_amino_acids:
                                pos_diff = pos - startposition
                                f = f[:pos_diff] + 'snp%sto%sat%ssnp' % (aach.split('>')[0], aach.split('>')[-1], pos) + f[pos_diff+1:]
                                yield f if not position else (f, startposition)
                        else:
                            yield f if not position else (f, startposition)
                    else:
                        yield f if not position else (f, startposition)
            loopcnt -= 1
			
def build_pept_prot(settings, results):
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    isdecoy = is_decoy_function(settings)

    snp = settings.getint('search', 'snp')
    pept_prot = {}
    prots = {}
    peptides = set()
    pept_neighbors = {}
    pept_ntts = {}
    enzyme = settings.get('search', 'enzyme')
    semitryptic = settings.getint('search', 'semitryptic')
    for x in results:
        peptides.update(re.sub(r'[^A-Z]', '', normalize_mods(x['candidates'][i][1], settings)) for i in range(
            1 or len(x['candidates'])))
    seen_target.clear()
    seen_decoy.clear()
    enzyme_rule = get_enzyme(enzyme)
    for desc, prot in prot_gen(settings):
        dbinfo = desc.split(' ')[0]
        prots[dbinfo] = desc
        if semitryptic:
            cl_positions = set(z for z in it.chain([x.end() for x in re.finditer(enzyme_rule, prot)],
                   [0, 1, len(prot)]))
        for pep, startposition in prot_peptides(prot, enzyme_rule, mc, minlen, maxlen, isdecoy(desc), dont_use_seen_peptides=True, snp=snp, desc=desc, position=True, semitryptic=semitryptic):
            if snp:
                if 'snp' not in pep:
                    seqm = pep
                else:
                    tmp = pep.split('snp')
                    seqm = tmp[0] + tmp[1].split('at')[0].split('to')[-1] + tmp[2]
            else:
                seqm = pep
            if seqm in peptides:
                if not semitryptic:
                    pept_prot.setdefault(seqm, []).append(dbinfo)
                    pept_neighbors[seqm] = (prot[startposition-1] if startposition != 0 else 'N/A',
                        prot[startposition+len(seqm)] if startposition + len(seqm) < len(prot) else 'N/A',
                                            startposition, min(startposition + len(seqm), len(prot)))
                    pept_ntts[seqm] = 2
                else:
                    ntt = (startposition in cl_positions) + ((startposition + len(seqm)) in cl_positions)
                    if seqm in pept_ntts:
                        best_ntt = pept_ntts[seqm]
                        if best_ntt <= ntt:
                            if best_ntt < ntt:
                                del pept_prot[seqm]
                                del pept_ntts[seqm]
                            pept_prot.setdefault(seqm, []).append(dbinfo)
                            pept_neighbors[seqm] = (prot[startposition-1] if startposition != 0 else 'N/A',
                                prot[startposition+len(seqm)] if startposition + len(seqm) < len(prot) else 'N/A',
                                                    startposition, min(startposition + len(seqm), len(prot)))
                            pept_ntts[seqm] = ntt
                    else:
                        pept_prot.setdefault(seqm, []).append(dbinfo)
                        pept_neighbors[seqm] = (prot[startposition-1] if startposition != 0 else 'N/A',
                            prot[startposition+len(seqm)] if startposition + len(seqm) < len(prot) else 'N/A')
                        pept_ntts[seqm] = ntt

    return pept_prot, prots, pept_neighbors, pept_ntts

def write_output(inputfile, settings, results):
    formats = {'pepxml': write_pepxml, 'csv': write_csv, 'tsv': write_csv, 'pickle': write_pickle}
    of = settings.get('output', 'format')
    writer = formats[re.sub(r'[^a-z]', '', of.lower())]

    if settings.has_option('output', 'path'):
        outd = settings.get('output', 'path')
        if not os.path.isdir(outd):
            logger.info('Creating %s ...', outd)
            os.makedirs(outd)
    else:
        outpath = os.path.dirname(inputfile)
        settings.set('output', 'path', outpath)

    return writer(inputfile, settings, results)


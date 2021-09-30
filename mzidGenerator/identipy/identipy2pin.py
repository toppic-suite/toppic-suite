import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta
import numpy as np
from os import path
from collections import Counter
import argparse

def calc_RT(seq, RC):
    try:
        return achrom.calculate_RT(seq, RC)
    except:
        return 0
    
def is_decoy(proteins, decoy_prefix):
    return all(z.startswith(decoy_prefix) for z in proteins)

def parse_mods(df_raw):
    mods_counter = Counter()
    sequence, mods = df_raw['peptide'], df_raw['modifications']
    if isinstance(mods, list):
        for mod in mods:
            mod_mass, aa_ind = mod.split('@')
            mod_mass = float(mod_mass)
            aa_ind = int(aa_ind)
            if aa_ind == 0:
                aa = 'N_term'
                mod_mass = round(mod_mass - 1.007825, 3)
            elif aa_ind == len(sequence) + 1:
                aa = 'C_term'
                mod_mass = round(mod_mass - 17.002735, 3)
            else:
                aa = sequence[aa_ind-1]
                mod_mass = round(mod_mass - mass.std_aa_mass[aa], 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] += 1
    return mods_counter

def add_mod_info(df_raw, mod):
    sequence, mods_counter = df_raw['peptide'], df_raw['mods_counter']
    mod_aa = mod.split(' at ')[1]
    if 'term' not in mod_aa and mod_aa not in sequence:
        return -1
    else:
        return mods_counter.get(mod, 0)

def prepare_mods(df):
    all_mods = set()
    for cnt in df['mods_counter'].values:
        for k in cnt.keys():
            all_mods.add(k)
    for mod in all_mods:
        df[mod] = df.apply(add_mod_info, axis=1, mod=mod)
    return df

def getlabel(decoy):
    return -1 if decoy else 1

def prepare_dataframe(infile_path, decoy_prefix='DECOY_', use_rt=1):
    df1 = pepxml.DataFrame(infile_path, read_schema=False)
    df1['length'] = df1['peptide'].apply(len)
    try:
        df1['y-b_ions'] = df1['matched_y1_ions'] - df1['matched_b1_ions']
    except:
        pass
    df1 = df1[df1['length'] >= 6]
    df1['spectrum'] = df1['spectrum'].apply(lambda x: x.split(' RTINSECONDS')[0])
    df1['massdiff_int'] = df1['massdiff'].apply(lambda x: int(round(x, 0)))
    df1['massdiff_ppm'] = 1e6 * df1['massdiff'] / df1['calc_neutral_pep_mass']
    df1['decoy'] = df1['protein'].apply(is_decoy, decoy_prefix=decoy_prefix)
    df1['mods_counter'] = df1.apply(parse_mods, axis=1)
    df1 = prepare_mods(df1)
    
    if use_rt:
        try:
            df1['RT exp'] = df1['retention_time_sec'] / 60
            df1 = df1.drop(['retention_time_sec', ], axis=1)
            df1_f = aux.filter(df1, fdr=0.01, key='expect', is_decoy='decoy', correction=1)
            print('Default target-decoy filtering, 1%% PSM FDR: Number of target PSMs = %d' \
                    % (df1_f[~df1_f['decoy']].shape[0]))
            print('Calibrating retention model...')
            retention_coefficients = achrom.get_RCs_vary_lcp(df1_f['peptide'].values, \
                                                            df1_f['RT exp'].values)
            df1_f['RT pred'] = df1_f['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
            df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
            _, _, r_value, std_value = aux.linear_regression(df1_f['RT pred'], df1_f['RT exp'])
            print('R^2 = %f , std = %f' % (r_value**2, std_value))
            df1['RT diff'] = df1['RT pred'] - df1['RT exp']
        except:
            pass

    df1['Label'] = df1['decoy'].apply(getlabel)
    df1['SpecId'] = df1['index'] + 1
    df1['ScanNr'] = df1['index'] + 1
    try:
        prev_aa = df1['peptide_prev_aa'][0]
        next_aa = df1['peptide_next_aa'][0]
        df1['Peptide'] = df1['peptide'].apply(lambda x: prev_aa + '.' + x + '.' + next_aa)
    except:
        df1['Peptide'] = df1['peptide'].apply(lambda x: 'K.' + x + '.K')
    df1['Proteins'] = df1['protein']
    
    return df1

def get_features(dataframe):
    feature_columns = dataframe.columns
    columns_to_remove = []
    for feature in feature_columns:
        if feature not in ['expect', 'hyperscore', 'calc_neutral_pep_mass', 'bscore', 'yscore', \
                            'massdiff', 'massdiff_ppm', 'RT pred', 'RT diff', \
                            'sumI', 'RT exp', 'precursor_neutral_mass', 'massdiff_int', \
                            'num_missed_cleavages', 'tot_num_ions', 'num_matched_ions', 'length',\
                            'SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins',
                            'matched_y1_ions', 'matched_b1_ions', 'y-b_ions', 'fragmentMT']:
            if not feature.startswith('mass shift'):
                columns_to_remove.append(feature)
    feature_columns = feature_columns.drop(columns_to_remove)
    return feature_columns

def run():
    parser = argparse.ArgumentParser(
        description='Convert Identipy pep.xml to pin for Percolator',
        epilog='''

    Example usage
    -------------
    $ identipy2pin input.pep.xml
    -------------

    Also can be used for MSFragger and X!Tandem pep.xml files.
    ''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file',     help='input .pep.xml file')
    parser.add_argument('-out',     help='path to output .pin file. By default put pin file near the pep.xml', default='')
    parser.add_argument('-prefix',  help='decoy prefix', default='DECOY_')
    parser.add_argument('-rt',      help='add RT prediction to features. 1 or 0', default=1, type=int)


    args = vars(parser.parse_args())
    infile = args['file']
    prefix = args['prefix']
    use_rt = args['rt']
    out = args['out']
    if out:
        outfile = out
    else:
        outfile = infile.replace('.pep.xml', '.pin')
    df1 = prepare_dataframe(infile, decoy_prefix=prefix, use_rt=use_rt)
    df00 = df1[get_features(df1)]
    df00_col = list(df00.columns.values)
    df00_col.remove('SpecId')
    df00_col.remove('Label')
    df00_col.remove('ScanNr')
    df00_col.remove('Peptide')
    df00_col.remove('Proteins')

    df00_col.insert(0, 'ScanNr')
    df00_col.insert(0, 'Label')
    df00_col.insert(0, 'SpecId')
    df00_col.append('Peptide')
    df00_col.append('Proteins')

    dft = df00.reindex(columns=df00_col)
    dft['Proteins'] = dft['Proteins'].apply(lambda x: 'proteinsplittmp'.join(x))
    dft.to_csv(path_or_buf=outfile, index=False, sep='\t')

    lines = list(open(outfile, 'r').readlines())
    outf = open(outfile, 'w')
    for l in lines:
        tmp = l.split('\t')
        outf.write('\t'.join(tmp[:-1]) + '\t' + '\t'.join(tmp[-1].split('proteinsplittmp')))
    outf.close()
    
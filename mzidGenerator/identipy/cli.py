import argparse
import string
import logging.config
import os
import subprocess

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'formatters': {
        'verbose': {
            'format': '%(levelname)7s %(asctime)s %(module)s %(process)d %(thread)d %(message)s',
        },
        'simple': {
            'format': '%(levelname)7s: %(asctime)s %(message)s',
            'datefmt': '[%H:%M:%S]',
        },
    },
    'handlers': {
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple',
        },
    },
    'loggers': {
        'identipy': {
            'handlers': ['console'],
            'level': 'INFO',
        }
    }
}

logging.config.dictConfig(LOGGING)
import logging
logger = logging.getLogger(__name__)
from . import main, utils


def get_label(modmass, labels):
    abt = string.ascii_lowercase
    abt_l = len(abt) - 1
    if modmass in labels:
        return labels[modmass], labels, 0
    else:
        labels[modmass] = abt[labels['i']] + abt[labels['j']] + abt[labels['k']]
        labels['k'] += 1
        if labels['k'] > abt_l:
            labels['k'] = 0
            labels['j'] += 1
        if labels['j'] > abt_l:
            labels['j'] = 0
            labels['i'] += 1
        return labels[modmass], labels, 1

def _update(settings, section, name, value):
    if value is not None:
        settings.set(section, name, value)

def run():
    parser = argparse.ArgumentParser(
        description='Search proteins using LC-MS/MS spectra',
        epilog='''

    Example usage
    -------------
    $ identipy input.mgf -db human.fasta
    -------------
    ''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file',     help='input .mzML or .mgf file with MS/MS spectra')
    parser.add_argument('-db',      help='path to protein fasta file')
    parser.add_argument('-cfg',     help='path to file with parameters')
    parser.add_argument('-out',     help='output path')
    parser.add_argument('-of',      help='output format')
    parser.add_argument('-sep',     help='output column separator')
    parser.add_argument('-at',      help='Use auto-tuning of search parameters', action='store_true')
    parser.add_argument('-nopwide', help='Do not increase initial precursor mass accuracy for auto-tuning', action='store_true')
    parser.add_argument('-punit',   help='precursor mass tolerance unit. Can be ppm or Da', type=str)
    parser.add_argument('-ptol',    help='precursor mass tolerance', type=float)
    parser.add_argument('-lptol',   help='*left precursor mass tolerance', type=float)
    parser.add_argument('-rptol',   help='*right precursor mass tolerance', type=float)
    parser.add_argument('-funit',   help='fragment mass tolerance unit. Can be ppm or Da')
    parser.add_argument('-ftol',    help='fragment mass tolerance', type=float)
    parser.add_argument('-fminmz',  help='fragment min m/z', type=float)
    parser.add_argument('-lmin',    help='min length of peptides', type=int)
    parser.add_argument('-lmax',    help='max length of peptides', type=int)
    parser.add_argument('-massmin', help='min mass of peptides', type=int)
    parser.add_argument('-massmax', help='max mass of peptides', type=int)
    parser.add_argument('-e',       help='cleavage rule in quotes!. X!Tandem style for cleavage rules')
    parser.add_argument('-mc',      help='number of missed cleavages', type=int)
    parser.add_argument('-semi',    help='include semitryptic peptides', type=int)
    parser.add_argument('-cmin',    help='min precursor charge', type=int)
    parser.add_argument('-cmax',    help='max precursor charge', type=int)
    parser.add_argument('-cumin',   help='min unknown precursor charge', type=int)
    parser.add_argument('-cumax',   help='max unknown precursor charge', type=int)
    parser.add_argument('-ime',     help='precursor isotope mass error. The parent ion\
     mass tolerance is expanded by opening up multiple tolerance windows centered\
      on the given number of 13C isotope peaks for a peptide.', type=int)
    parser.add_argument('-shifts',  help='shifts. example: 0,16.000,23.000,12')
    parser.add_argument('-snp',     help='1 means make SNP changes for ALL peptides', type=int)
    parser.add_argument('-mm',      help='number of minimum matched ions', type=int)
    parser.add_argument('-ad',      help='add decoy', action='store_true')
    parser.add_argument('-prefix',  help='decoy prefix')
    parser.add_argument('-infix',   help='decoy infix')
    parser.add_argument('-method',  help='reverse or shuffle')
    parser.add_argument('-deis',    help='use MS/MS deisotoping. yes or no')
    parser.add_argument('-deistol', help='deisotope mass accuracy', type=float)
    parser.add_argument('-score',   help='used score. Can be RNHS, hyperscore or morpheusscore')
    parser.add_argument('-minp',    help='minumum peaks in MS/MS spectra', type=int)
    parser.add_argument('-maxp',    help='maximum peaks in MS/MS spectra', type=int)
    parser.add_argument('-dyn',     help='dynamic range', type=float)
    parser.add_argument('-mfc',     help='maximum fragment charge', type=int)
    parser.add_argument('-nproc',   help='number of processes. 0 means auto', type=int)
    parser.add_argument('-maxmods', help='maximum variable mods per sequence', type=int)
    parser.add_argument('-ncleave', help='protein nterm cleavage', type=float)
    parser.add_argument('-ccleave', help='protein cterm cleavage', type=float)
    parser.add_argument('-fmods',   help='fixed modifications. in mass1@aminoacid1,mass2@aminoacid2 format')
    parser.add_argument('-vmods',   help='variable modifications. in mass1@aminoacid1,mass2@aminoacid2 format')
    parser.add_argument('-tags',    help='Add quantitation tags to the pepXML output. Can be tmt10plex, tmt6plex, tmt11plex or custom format label1:mass1,label2:mass2...')
    parser.add_argument('-debug',  help='Print debugging messages', action='store_true')
    parser.add_argument('-dino', help='path to Dinosaur. Used for chimeric spectrum processing and MS1 Intensity calculation', default=False)

    args = vars(parser.parse_args())
    if args['debug']:
        logging.getLogger('identipy').setLevel(logging.DEBUG)

    if args['cfg']:
        settings = main.settings(args['cfg'])
    else:
        settings = main.settings()

    labels = {'i': 0, 'j': 0, 'k': 0}

    fmods_array = []
    vmods_array = []
    if args['fmods']:
        for mod in args['fmods'].split(','):
            modmass, modaa = mod.split('@')
            lbl, labels, flag = get_label(modmass, labels)
            if modaa == '[':
                ntermlabel, modaa, ctermlabel = '-', '', ''
            elif modaa == ']':
                ntermlabel, modaa, ctermlabel = '', '', '-'
            else:
                ntermlabel, ctermlabel = '', ''
            fmods_array.append(ctermlabel + lbl + modaa + ntermlabel)
            if flag:
                settings.set('modifications', lbl, modmass)
    if fmods_array or args['fmods'] is not None:
        settings.set('modifications', 'fixed', ','.join(fmods_array))

    if args['vmods']:
        for mod in args['vmods'].split(','):
            modmass, modaa = mod.split('@')
            lbl, labels, flag = get_label(modmass, labels)
            if modaa == '[':
                ntermlabel, modaa, ctermlabel = '-', '', ''
            elif modaa == ']':
                ntermlabel, modaa, ctermlabel = '', '', '-'
            else:
                ntermlabel, ctermlabel = '', ''
            vmods_array.append(ctermlabel + lbl + modaa + ntermlabel)
            if flag:
                settings.set('modifications', lbl, modmass)
    if vmods_array or args['vmods'] is not None:
        settings.set('modifications', 'variable', ','.join(vmods_array))


    _update(settings, 'input',  'database', args['db'])
    _update(settings, 'search', 'precursor accuracy unit', args['punit'])
    _update(settings, 'search', 'precursor accuracy left', (args['ptol'] if not args['lptol'] else args['lptol']))
    _update(settings, 'search', 'precursor accuracy right', (args['ptol'] if not args['rptol'] else args['rptol']))
    _update(settings, 'search', 'product accuracy unit', args['funit'])
    _update(settings, 'search', 'product accuracy', args['ftol'])
    _update(settings, 'search', 'product minimum m/z', args['fminmz'])
    _update(settings, 'search', 'peptide maximum length', args['lmax'])
    _update(settings, 'search', 'peptide minimum length', args['lmin'])
    _update(settings, 'search', 'peptide maximum mass', args['massmax'])
    _update(settings, 'search', 'peptide minimum mass', args['massmin'])
    _update(settings, 'search', 'enzyme', args['e'])
    _update(settings, 'search', 'number of missed cleavages', args['mc'])
    _update(settings, 'search', 'semitryptic', args['semi'])
    _update(settings, 'search', 'maximum charge', args['cmax'])
    _update(settings, 'search', 'minimum charge', args['cmin'])
    _update(settings, 'search', 'maximum unknown charge', args['cumax'])
    _update(settings, 'search', 'minimum unknown charge', args['cumin'])
    _update(settings, 'search', 'precursor isotope mass error', args['ime'])
    _update(settings, 'search', 'shifts', args['shifts'])
    _update(settings, 'search', 'snp', args['snp'])
    _update(settings, 'output', 'minimum matched', args['mm'])
    if args['ad']:
        _update(settings, 'input', 'add decoy', 'yes')
    _update(settings, 'input',  'decoy prefix', args['prefix'])
    _update(settings, 'input',  'decoy infix', args['infix'])
    _update(settings, 'input',  'decoy method', args['method'])
    _update(settings, 'input',  'deisotope', args['deis'])
    _update(settings, 'input',  'deisotoping mass tolerance', args['deistol'])
    if args['score']:
        _update(settings, 'scoring', 'score', 'identipy.scoring.' + args['score'])
    _update(settings, 'scoring', 'minimum peaks', args['minp'])
    _update(settings, 'scoring', 'maximum peaks', args['maxp'])
    _update(settings, 'scoring', 'dynamic range', args['dyn'])
    _update(settings, 'scoring', 'maximum fragment charge', args['mfc'])
    _update(settings, 'performance', 'processes', args['nproc'])
    _update(settings, 'modifications', 'maximum variable mods', args['maxmods'])
    _update(settings, 'modifications', 'protein nterm cleavage', args['ncleave'])
    _update(settings, 'modifications', 'protein cterm cleavage', args['ccleave'])
    _update(settings, 'output', 'path', args['out'])
    _update(settings, 'output', 'format', args['of'])
    _update(settings, 'output', 'separator', args['sep'])
    _update(settings, 'output', 'tags', args['tags'])
    if args['at']:
        ao_setting = 'identipy.extras.optimization'
        if args['nopwide']:
            _update(settings, 'optimization', 'increase precursor mass tolerance', 'no')
    else:
        ao_setting = None
    _update(settings, 'misc', 'first stage', ao_setting)

    inputfile = args['file']

    dino_path = args['dino']
    if dino_path:
        if os.path.splitext(inputfile)[1].lower() != '.mzml':
            logger.info('Only mzml supported for Dinosaur!\n')
        else:
            try:
                advpath = '--advParams=' + os.path.join(os.path.dirname(os.path.realpath(__file__)), 'adv.txt')
                logger.info('Start Dinosaur...\n')
                subprocess.call(['java', '-Djava.awt.headless=true', '-jar', os.path.realpath(dino_path), advpath, '--concurrency=12', inputfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                path_to_features = os.path.splitext(inputfile)[0] + os.extsep + 'features' + os.extsep + 'tsv'
                logger.info('Start demultiplexing...\n')
                path_to_mgf = utils.demix_chimeric(path_to_features, inputfile, 0.65)
                logger.info('Demultiplexing was finished...\n')
                utils.write_output(path_to_mgf, settings, main.process_file(path_to_mgf, settings))
                return
            except Exception as e:
                logger.error(e)

    utils.write_output(inputfile, settings, main.process_file(inputfile, settings))


if __name__ == '__main__':
    run()
